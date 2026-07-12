#!/usr/bin/env bash
# GPU-parallel launcher for the z_O sensitivity job set K in jobs_zscan/.
# z_O = 0.90/1.02 A, each with delta = 0.00/0.15/0.30/0.50/0.75 A:
# 10 vacuum sqrt3 x sqrt3 Na_1/3CoO2 SCF jobs at fixed c = 9.9 A.
# One pw.x SCF job per GPU at a time, pulled from a shared queue so all
# 10 jobs interleave across the available GPUs.  Idempotent: any job whose
# pw.out already ends in "JOB DONE" is skipped, so re-running this script
# resumes where it left off.  Kill switch: create a file named STOP in this
# directory; workers exit after their current job.
#
# There is NO PP-tools bootstrap here -- set K is pw.x-only (one SCF per job).
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PATH=/usr/local/qe/bin:/usr/local/openmpi/bin:/usr/local/ucx/bin:/usr/local/nvidia/bin:/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib:/usr/local/cuda/lib64:/usr/local/fftw/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/lib:$LD_LIBRARY_PATH
export PMIX_MCA_gds=hash
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-14}"
PW_BIN="${PW_BIN:-pw.x}"

command -v "$PW_BIN" >/dev/null || { echo "ERROR: $PW_BIN not in PATH"; exit 1; }

NGPU=$(nvidia-smi -L 2>/dev/null | wc -l)
[ "$NGPU" -ge 1 ] || { echo "WARNING: no GPU detected, using 1 serial slot"; NGPU=1; }
echo "== $NGPU GPU(s) detected =="

[ -d pseudo ] && ls pseudo/*.UPF >/dev/null 2>&1 || bash get_pseudos.sh
[ -f manifest_zscan.json ] || python3 generate_inputs.py --stage zscan

step_done() {  # $1 = output file
  [ -f "$1" ] && tail -n 20 "$1" | grep -q "JOB DONE"
}

job_done() {  # $1 = job dir; set-K jobs are single-step (pw.in -> pw.out)
  step_done "$1/pw.out"
}

run_job() {  # $1 = gpu id, $2 = job name
  local gpu="$1" name="$2" dir="jobs_zscan/$2"
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

list_jobs() {  # all set-K job names, manifest order
  python3 - <<'EOF'
import json
for j in json.load(open("manifest_zscan.json"))["jobs"]:
    print(j["name"])
EOF
}

echo "== set K: z_O sensitivity SCF jobs =="
list_jobs > .queue_zscan
run_pool .queue_zscan
[ -f STOP ] && { echo "stopped by STOP file"; exit 1; }

echo "== all zscan-set jobs finished =="
echo "Summary:"
n_ok=0; n_bad=0
for d in jobs_zscan/*/; do
  if job_done "$d"; then n_ok=$((n_ok+1)); else n_bad=$((n_bad+1)); echo "  incomplete: $d"; fi
done
echo "  $n_ok done, $n_bad incomplete"
echo "  analyze: fit alpha(z_O) at fixed c across each delta scan, then compare"
echo "  d(alpha)/d(z_O) with the set D single-point z_O sensitivity;"
echo "  see the 'analyze' note in manifest_zscan.json."
echo "ZSCAN COMPLETE"
