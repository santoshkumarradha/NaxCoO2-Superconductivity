#!/usr/bin/env bash
# GPU pipeline for the "cage-ensemble" set L/M.
#   Phase 1: set L -- one Born-Oppenheimer MD job (jobs_ensemble/md_cage),
#            water free / Na+CoO2 frozen, 290 K svr thermostat, run on ALL
#            4 GPUs at once via mpirun -np 4 pw.x -nk 2 (falls back to -nk 1
#            if pools misbehave).
#   Phase 2: parse the MD trajectory and write the set-M frozen-snapshot Na
#            scans (make_ensemble_scans.py; pod-side, needs the finished/
#            far-enough-along MD output).
#   Phase 3: run the 70 set-M SCF jobs through the flock 4-GPU pool, exactly
#            like run_mobile.sh.
# On completion of all jobs, writes the marker line "ENSEMBLE COMPLETE" to
# ensemble.log (this file's own stdout, redirected by the caller/bootstrap).
#
# Robustness: PMIX_MCA_gds=hash + /tmp/pmix* cleanup (stale PMIx sessions
# wedge multi-rank mpirun launches across pod restarts), the PATH/
# LD_LIBRARY_PATH exports used by every other run_*.sh in this repo, per-job
# 'JOB DONE' checks (idempotent -- safe to re-run), and phase 1 verifies the
# MD produced >= 3000 steps before proceeding (else logs "MD FAILED" and
# stops -- phases 2/3 never run on a dead or barely-started MD).
# Kill switch: create a file named STOP in this directory; the phase-3 pool
# workers exit after their current job (phase 1's single MD job is NOT
# interrupted by STOP -- kill it manually if needed).
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PATH=/usr/local/qe/bin:/usr/local/openmpi/bin:/usr/local/ucx/bin:/usr/local/nvidia/bin:/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib:/usr/local/cuda/lib64:/usr/local/fftw/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/lib:$LD_LIBRARY_PATH
export PMIX_MCA_gds=hash
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-14}"
PW_BIN="${PW_BIN:-pw.x}"
MIN_MD_STEPS=3000

rm -rf /tmp/pmix* 2>/dev/null || true

command -v "$PW_BIN" >/dev/null || { echo "ERROR: $PW_BIN not in PATH"; exit 1; }

NGPU=$(nvidia-smi -L 2>/dev/null | wc -l)
[ "$NGPU" -ge 1 ] || { echo "WARNING: no GPU detected, using 1 serial slot"; NGPU=1; }
echo "== $NGPU GPU(s) detected =="

[ -d pseudo ] && ls pseudo/*.UPF >/dev/null 2>&1 || bash get_pseudos.sh
[ -s pseudo/H.pbe-kjpaw_psl.1.0.0.UPF ] || bash get_pseudos.sh
[ -f manifest_ensemble.json ] || python3 generate_inputs.py --stage ensemble

step_done() {  # $1 = output file
  [ -f "$1" ] && tail -n 20 "$1" | grep -q "JOB DONE"
}

count_md_steps() {  # $1 = pw.out; cheap proxy, matches make_ensemble_scans.py
  [ -f "$1" ] && grep -c "ATOMIC_POSITIONS" "$1" || echo 0
}

# ------------------------------------------------------------- phase 1: MD --
echo "== phase 1: set L Born-Oppenheimer MD (md_cage, all $NGPU GPUs) =="
MD_DIR="jobs_ensemble/md_cage"
if step_done "$MD_DIR/pw.out"; then
  echo "  md_cage: already JOB DONE, skip"
else
  rm -rf /tmp/pmix* 2>/dev/null || true
  echo "  md_cage: mpirun -np $NGPU pw.x -nk 2 start $(date +%H:%M:%S)"
  ( cd "$MD_DIR" \
    && CUDA_VISIBLE_DEVICES=$(seq -s, 0 $((NGPU - 1))) \
       mpirun -np "$NGPU" "$PW_BIN" -nk 2 -input pw.in > pw.out 2>&1 )
  if ! step_done "$MD_DIR/pw.out"; then
    echo "  md_cage: -nk 2 did not finish cleanly, retrying with -nk 1"
    rm -rf /tmp/pmix* 2>/dev/null || true
    ( cd "$MD_DIR" \
      && CUDA_VISIBLE_DEVICES=$(seq -s, 0 $((NGPU - 1))) \
         mpirun -np "$NGPU" "$PW_BIN" -nk 1 -input pw.in > pw.out 2>&1 )
  fi
fi

MD_STEPS=$(count_md_steps "$MD_DIR/pw.out")
echo "  md_cage: $MD_STEPS ATOMIC_POSITIONS blocks in pw.out (need >= $MIN_MD_STEPS)"
if [ "$MD_STEPS" -lt "$MIN_MD_STEPS" ]; then
  echo "MD FAILED"
  exit 1
fi
echo "  md_cage: MD OK ($MD_STEPS steps) $(date +%H:%M:%S)"

# --------------------------------------------------- phase 2: snapshot gen --
echo "== phase 2: extract thermal snapshots + write set-M SCF inputs =="
python3 make_ensemble_scans.py
[ $? -eq 0 ] || { echo "ERROR: make_ensemble_scans.py failed"; exit 1; }

# -------------------------------------------- phase 3: set-M SCF job pool ---
job_done() {  # $1 = job dir; set-M jobs are single-step (pw.in -> pw.out)
  step_done "$1/pw.out"
}

run_job() {  # $1 = gpu id, $2 = job name
  local gpu="$1" name="$2" dir="jobs_ensemble/$2"
  if job_done "$dir"; then echo "[gpu$gpu] $name: already done, skip"; return 0; fi
  echo "[gpu$gpu] $name: pw.x scf start $(date +%H:%M:%S)"
  ( cd "$dir" && CUDA_VISIBLE_DEVICES="$gpu" "$PW_BIN" -nk 1 -input pw.in > pw.out 2>&1 )
  step_done "$dir/pw.out" \
    || { echo "[gpu$gpu] $name: pw.x FAILED (see $dir/pw.out)"; return 1; }
  echo "[gpu$gpu] $name: pw.x JOB DONE $(date +%H:%M:%S)"
}

run_pool() {  # $1 = file with one job name per line
  local queue="$1"
  [ -s "$queue" ] || { echo "  (queue empty)"; return 0; }
  local lock="$queue.lock"; : > "$lock"
  worker() {
    local gpu="$1" name
    while :; do
      [ -f STOP ] && { echo "[gpu$gpu] STOP file found, worker exiting"; break; }
      name=$( { flock 9; head -n 1 "$queue"; sed -i '1d' "$queue"; } 9>>"$lock" )
      [ -n "$name" ] || break
      run_job "$gpu" "$name"
    done
  }
  local g
  for g in $(seq 0 $((NGPU - 1))); do worker "$g" & done
  wait
  rm -f "$lock"
}

list_m_jobs() {  # set-M job names, manifest order
  python3 - <<'EOF'
import json
for j in json.load(open("manifest_ensemble.json"))["jobs"]:
    if j.get("set") == "M":
        print(j["name"])
EOF
}

echo "== phase 3: set M -- 70 frozen-snapshot Na SCF scans =="
list_m_jobs > .queue_ensemble
run_pool .queue_ensemble
[ -f STOP ] && { echo "stopped by STOP file"; exit 1; }

echo "== all ensemble jobs finished =="
echo "Summary:"
n_ok=0; n_bad=0
for d in jobs_ensemble/*/; do
  if job_done "$d"; then n_ok=$((n_ok+1)); else n_bad=$((n_bad+1)); echo "  incomplete: $d"; fi
done
echo "  $n_ok done, $n_bad incomplete"
echo "  analyze: theory/results_extra_analysis/analyze_ensemble.py fits each"
echo "  snapshot's E(delta) to E0 + alpha*delta^2 + beta*delta^4 and compares"
echo "  the ensemble mean curve against the set-G rigid-cage and set-J"
echo "  adiabatic results."
echo "ENSEMBLE COMPLETE"
