#!/usr/bin/env bash
# GPU pipeline for the "cage-ensemble" set L/M (4-walker design).
#   Phase 1: set L -- FOUR independent Born-Oppenheimer MD walkers
#            (jobs_ensemble/md_w1..md_w4), water free / Na+CoO2 frozen,
#            290 K svr thermostat, nstep=1200 each, run CONCURRENTLY, one
#            walker per GPU (CUDA_VISIBLE_DEVICES=$i, mpirun -np 1).
#            Gate: each walker must reach >= 1000 steps; walkers below the
#            gate are dropped (survivors-only). If NO walker survives, log
#            "MD FAILED" and stop.
#   Phase 2: parse the surviving trajectories and write the set-M
#            frozen-snapshot Na scans (make_ensemble_scans.py; 3 snapshots
#            per surviving walker at ~steps 600/900/1200, first 400
#            discarded) -- up to 12 snapshot units x 7 deltas = 84 SCF.
#   Phase 3: snapshot-unit warm-restart pool: each unit (w<K>s<J>, 7 SCFs)
#            is one work item on one GPU. delta=0.00 runs first; its ./out
#            (charge density + wavefunctions) is copied into each remaining
#            job of the unit, which carry startingwfc/startingpot='file',
#            in |delta|-increasing order. Full set-G production settings
#            (K_POINTS 6 6 2, conv_thr 1e-7). out/ dirs are deleted as the
#            unit completes (only pw.out is harvested; keeps the disk flat).
# On completion of all jobs, writes the marker line "ENSEMBLE COMPLETE" to
# stdout (redirected to run_ensemble.log by the bootstrap).
#
# Robustness: PMIX_MCA_gds=hash + /tmp/pmix* cleanup (stale PMIx sessions
# wedge mpirun launches across pod restarts), the PATH/LD_LIBRARY_PATH
# exports used by every other run_*.sh in this repo, per-job 'JOB DONE'
# checks (idempotent -- safe to re-run). Kill switch: create a file named
# STOP in this directory; the phase-3 pool workers exit after their current
# unit (phase 1's MD walkers are NOT interrupted by STOP -- kill them
# manually if needed).
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PATH=/usr/local/qe/bin:/usr/local/openmpi/bin:/usr/local/ucx/bin:/usr/local/nvidia/bin:/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib:/usr/local/cuda/lib64:/usr/local/fftw/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/lib:${LD_LIBRARY_PATH:-}
export PMIX_MCA_gds=hash
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-14}"
PW_BIN="${PW_BIN:-pw.x}"
MIN_WALKER_STEPS=1000

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

# ---------------------------------------------------- phase 1: MD walkers ---
echo "== phase 1: set L -- 4 BOMD walkers, one per GPU =="
i=0
for d in jobs_ensemble/md_w*/; do
  name=$(basename "$d")
  gpu=$((i % NGPU))
  i=$((i + 1))
  if step_done "$d/pw.out"; then
    echo "  $name: already JOB DONE, skip"
    continue
  fi
  steps_now=$(count_md_steps "$d/pw.out")
  if [ "$steps_now" -ge "$MIN_WALKER_STEPS" ]; then
    echo "  $name: already has $steps_now steps (>= $MIN_WALKER_STEPS), skip"
    continue
  fi
  echo "  $name: start on gpu$gpu $(date +%H:%M:%S)"
  ( cd "$d" && CUDA_VISIBLE_DEVICES="$gpu" \
      mpirun -np 1 "$PW_BIN" -nk 1 -input pw.in > pw.out 2>&1 ) &
done
wait

SURVIVORS=0
for d in jobs_ensemble/md_w*/; do
  name=$(basename "$d")
  steps=$(count_md_steps "$d/pw.out")
  if [ "$steps" -ge "$MIN_WALKER_STEPS" ]; then
    echo "  $name: $steps steps -- SURVIVOR"
    SURVIVORS=$((SURVIVORS + 1))
  else
    echo "  $name: $steps steps (< $MIN_WALKER_STEPS) -- dropped"
  fi
done
if [ "$SURVIVORS" -lt 1 ]; then
  echo "MD FAILED"
  exit 1
fi
echo "  phase 1 OK: $SURVIVORS surviving walker(s) $(date +%H:%M:%S)"

# --------------------------------------------------- phase 2: snapshot gen --
echo "== phase 2: extract thermal snapshots + write set-M SCF inputs =="
python3 make_ensemble_scans.py || { echo "ERROR: make_ensemble_scans.py failed"; exit 1; }

# ---------------------------------- phase 3: warm-restart snapshot units ----
job_done() {  # $1 = job dir
  step_done "$1/pw.out"
}

run_unit() {  # $1 = gpu id, $2 = unit name (w<K>s<J>); 7 SCFs, warm-restart
  local gpu="$1" unit="$2"
  local d0="jobs_ensemble/${unit}_d0.00"
  # delta order: 0.00 first (produces the warm-start density), then
  # |delta|-increasing
  local tag dir
  for tag in 0.00 +0.15 -0.15 +0.30 -0.30 +0.45 -0.45; do
    dir="jobs_ensemble/${unit}_d${tag}"
    [ -d "$dir" ] || continue
    if job_done "$dir"; then echo "[gpu$gpu] ${unit}_d${tag}: already done, skip"; continue; fi
    if [ "$tag" != "0.00" ] && [ -d "$d0/out" ]; then
      rm -rf "$dir/out"
      cp -r "$d0/out" "$dir/out"
    fi
    echo "[gpu$gpu] ${unit}_d${tag}: pw.x scf start $(date +%H:%M:%S)"
    ( cd "$dir" && CUDA_VISIBLE_DEVICES="$gpu" "$PW_BIN" -nk 1 -input pw.in > pw.out 2>&1 )
    if job_done "$dir"; then
      echo "[gpu$gpu] ${unit}_d${tag}: JOB DONE $(date +%H:%M:%S)"
    else
      echo "[gpu$gpu] ${unit}_d${tag}: pw.x FAILED (see $dir/pw.out)"
    fi
    # free disk: the warm-start copies are large and only pw.out is harvested
    [ "$tag" != "0.00" ] && rm -rf "$dir/out"
  done
  rm -rf "$d0/out"
}

run_pool() {  # $1 = file with one unit name per line
  local queue="$1"
  [ -s "$queue" ] || { echo "  (queue empty)"; return 0; }
  local lock="$queue.lock"; : > "$lock"
  worker() {
    local gpu="$1" unit
    while :; do
      [ -f STOP ] && { echo "[gpu$gpu] STOP file found, worker exiting"; break; }
      unit=$( { flock 9; head -n 1 "$queue"; sed -i '1d' "$queue"; } 9>>"$lock" )
      [ -n "$unit" ] || break
      run_unit "$gpu" "$unit"
    done
  }
  local g
  for g in $(seq 0 $((NGPU - 1))); do worker "$g" & done
  wait
  rm -f "$lock"
}

list_units() {  # unique set-M snapshot units, manifest order
  python3 - <<'EOF'
import json
seen = []
for j in json.load(open("manifest_ensemble.json"))["jobs"]:
    u = j.get("unit")
    if j.get("set") == "M" and u and u not in seen:
        seen.append(u)
print("\n".join(seen))
EOF
}

echo "== phase 3: set M -- frozen-snapshot Na scans (warm-restart units) =="
list_units > .queue_ensemble
echo "  $(wc -l < .queue_ensemble) snapshot units queued"
run_pool .queue_ensemble
[ -f STOP ] && { echo "stopped by STOP file"; exit 1; }

echo "== all ensemble jobs finished =="
echo "Summary:"
n_ok=0; n_bad=0
for d in jobs_ensemble/*/; do
  case "$(basename "$d")" in md_w*)
    steps=$(count_md_steps "$d/pw.out")
    [ "$steps" -ge "$MIN_WALKER_STEPS" ] || { n_bad=$((n_bad+1)); echo "  incomplete walker: $d ($steps steps)"; continue; }
    n_ok=$((n_ok+1)); continue ;;
  esac
  if job_done "$d"; then n_ok=$((n_ok+1)); else n_bad=$((n_bad+1)); echo "  incomplete: $d"; fi
done
echo "  $n_ok done, $n_bad incomplete"
echo "  analyze: theory/results_extra_analysis/analyze_ensemble.py fits each"
echo "  snapshot's E(delta) to E0 + alpha*delta^2 + beta*delta^4 and compares"
echo "  the ensemble mean curve against the set-G rigid-cage and set-J"
echo "  adiabatic results."
echo "ENSEMBLE COMPLETE"
