#!/usr/bin/env bash
# GPU-parallel launcher for the NaxCoO2 pw.x job set.
#   Phase 1: all SCF jobs (sets A + B), one pw.x per GPU concurrently.
#   Phase 2: generate NSCF jobs at the E(delta) minimum per (element,cell,c).
#   Phase 3: NSCF jobs + dos.x.
# Idempotent: any job whose pw.out already ends in "JOB DONE" is skipped, so
# re-running this script resumes where it left off.  Kill switch: create a
# file named STOP in this directory; workers exit after their current job.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-14}"
PW_BIN="${PW_BIN:-pw.x}"
DOS_BIN="${DOS_BIN:-dos.x}"

command -v "$PW_BIN" >/dev/null || { echo "ERROR: $PW_BIN not in PATH"; exit 1; }

NGPU=$(nvidia-smi -L 2>/dev/null | wc -l)
[ "$NGPU" -ge 1 ] || { echo "WARNING: no GPU detected, using 1 serial slot"; NGPU=1; }
echo "== $NGPU GPU(s) detected =="

[ -d pseudo ] && ls pseudo/*.UPF >/dev/null 2>&1 || bash get_pseudos.sh
[ -f manifest.json ] || python3 generate_inputs.py --stage scf

job_done() {  # $1 = job dir
  [ -f "$1/pw.out" ] && tail -n 20 "$1/pw.out" | grep -q "JOB DONE"
}

run_job() {  # $1 = gpu id, $2 = job name
  local gpu="$1" name="$2" dir="jobs/$2"
  if job_done "$dir"; then echo "[gpu$gpu] $name: already done, skip"; return 0; fi
  echo "[gpu$gpu] $name: start $(date +%H:%M:%S)"
  ( cd "$dir" && CUDA_VISIBLE_DEVICES="$gpu" "$PW_BIN" -nk 1 -input pw.in > pw.out 2>&1 )
  if job_done "$dir"; then
    echo "[gpu$gpu] $name: JOB DONE $(date +%H:%M:%S)"
    if [ -f "$dir/dos.in" ]; then
      ( cd "$dir" && "$DOS_BIN" < dos.in > dos.out 2>&1 ) \
        && echo "[gpu$gpu] $name: dos.x done" \
        || echo "[gpu$gpu] $name: dos.x FAILED"
    fi
  else
    echo "[gpu$gpu] $name: FAILED (see $dir/pw.out)"
  fi
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

list_jobs() {  # $1 = type (scf|nscf) -> job names, undone first-checked later
  python3 - "$1" <<'EOF'
import json, sys
jobs = json.load(open("manifest.json"))["jobs"]
for j in jobs:
    if j["type"] == sys.argv[1]:
        print(j["name"])
EOF
}

echo "== Phase 1: SCF jobs =="
list_jobs scf > .queue_scf
run_pool .queue_scf
[ -f STOP ] && { echo "stopped by STOP file"; exit 1; }

echo "== Phase 2: generating NSCF inputs at E(delta) minima =="
python3 generate_inputs.py --stage nscf

echo "== Phase 3: NSCF + DOS jobs =="
list_jobs nscf > .queue_nscf
run_pool .queue_nscf
[ -f STOP ] && { echo "stopped by STOP file"; exit 1; }

echo "== all phases finished =="
echo "Summary:"
n_ok=0; n_bad=0
for d in jobs/*/; do
  if job_done "$d"; then n_ok=$((n_ok+1)); else n_bad=$((n_bad+1)); echo "  incomplete: $d"; fi
done
echo "  $n_ok done, $n_bad incomplete"
echo "Next: bash bader_setup.sh   (Bader charges), then rsync results back."
