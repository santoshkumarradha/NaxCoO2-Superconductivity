#!/usr/bin/env bash
# GPU-parallel launcher for the referee-response job set I in jobs_bands/.
#   I1: Na fatband chains (reviewer figure + e_g' pocket check)   5 chains
#       per geometry: pw.x scf -> pw.x 'bands' (Gamma-M-K-Gamma, 60 pts)
#                     -> bands.x (spin up + dw) -> projwfc.x (all atoms)
#   I2: hydrate relaxed-water check (BFGS, Na z + H2O free)       2 relax
#   I3: hydrate vdW grimme-d3 spot-check                          3 SCF
# One job (= full chain) per GPU at a time, pulled from a shared queue so all
# 10 jobs interleave across the available GPUs.  Idempotent PER STEP: any
# step whose output file already ends in "JOB DONE" is skipped, so re-running
# this script resumes mid-chain.  Kill switch: create a file named STOP in
# this directory; workers exit after their current job.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

export PATH=/usr/local/qe/bin:/usr/local/openmpi/bin:/usr/local/ucx/bin:/usr/local/nvidia/bin:/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib:/usr/local/cuda/lib64:/usr/local/fftw/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/lib:$LD_LIBRARY_PATH
export PMIX_MCA_gds=hash
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-14}"
PW_BIN="${PW_BIN:-pw.x}"
BANDS_BIN="${BANDS_BIN:-bands.x}"
PROJWFC_BIN="${PROJWFC_BIN:-projwfc.x}"

command -v "$PW_BIN" >/dev/null || { echo "ERROR: $PW_BIN not in PATH"; exit 1; }

NGPU=$(nvidia-smi -L 2>/dev/null | wc -l)
[ "$NGPU" -ge 1 ] || { echo "WARNING: no GPU detected, using 1 serial slot"; NGPU=1; }
echo "== $NGPU GPU(s) detected =="

[ -d pseudo ] && ls pseudo/*.UPF >/dev/null 2>&1 || bash get_pseudos.sh
# set I needs the H pseudo (hydrate cells); fetch if pseudo/ predates it
[ -s pseudo/H.pbe-kjpaw_psl.1.0.0.UPF ] || bash get_pseudos.sh
[ -f manifest_bands.json ] || python3 generate_inputs.py --stage bands

# bands.x / projwfc.x are needed DURING the I1 chains, so ensure the CPU
# PP-tools build (same fallback as run_extra.sh) BEFORE the pool starts.
# NOTE: `make pp` builds bands.x and projwfc.x alongside dos.x/pp.x, but the
# run_extra.sh copy line only installed dos.x/pp.x -- the copy here is
# extended to bands.x/projwfc.x and reuses an existing /opt build if present.
ensure_pp_tools() {
  command -v "$BANDS_BIN" >/dev/null 2>&1 \
    && command -v "$PROJWFC_BIN" >/dev/null 2>&1 && return 0
  if [ -x /opt/q-e-qe-7.3.1/bin/bands.x ] && [ -x /opt/q-e-qe-7.3.1/bin/projwfc.x ]; then
    echo "== reusing existing PP-tools build in /opt/q-e-qe-7.3.1 =="
    cp /opt/q-e-qe-7.3.1/bin/bands.x /opt/q-e-qe-7.3.1/bin/projwfc.x /usr/local/bin/ && return 0
  fi
  echo "== compiling QE 7.3.1 PP tools (CPU) =="
  apt-get install -y -qq gfortran gcc make libfftw3-dev >/dev/null 2>&1
  ( cd /opt \
    && { [ -d q-e-qe-7.3.1 ] || { curl -fsSL -o qe.tar.gz https://gitlab.com/QEF/q-e/-/archive/qe-7.3.1/q-e-qe-7.3.1.tar.gz \
           && tar xzf qe.tar.gz; }; } \
    && cd q-e-qe-7.3.1 \
    && { [ -f make.inc ] || FC=gfortran F90=gfortran CC=gcc ./configure --disable-parallel > configure.log 2>&1; } \
    && make -j"$(nproc)" pp > make_pp.log 2>&1 \
    && cp bin/bands.x bin/projwfc.x bin/dos.x bin/pp.x /usr/local/bin/ ) \
    || { echo "PP TOOLS BUILD FAILED"; return 1; }
}
ensure_pp_tools || echo "WARNING: bands.x/projwfc.x unavailable; chains will stop after pw.x"

step_done() {  # $1 = output file
  [ -f "$1" ] && tail -n 20 "$1" | grep -q "JOB DONE"
}

job_done() {  # $1 = job dir; done = LAST step of that dir's chain is done
  local last="$1/pw.out"
  [ -f "$1/pw_bands.in" ] && last="$1/pw_bands.out"
  [ -f "$1/bands_up.in" ] && last="$1/bands_up.out"
  [ -f "$1/bands_dw.in" ] && last="$1/bands_dw.out"
  [ -f "$1/projwfc.in" ] && last="$1/projwfc.out"
  step_done "$last"
}

run_job() {  # $1 = gpu id, $2 = job name; runs the full chain for the dir
  local gpu="$1" name="$2" dir="jobs_bands/$2"
  if job_done "$dir"; then echo "[gpu$gpu] $name: already done, skip"; return 0; fi
  # step 1: pw.x scf (or relax for set I2)
  if step_done "$dir/pw.out"; then
    echo "[gpu$gpu] $name: pw.in already done, skip"
  else
    echo "[gpu$gpu] $name: pw.x start $(date +%H:%M:%S)"
    ( cd "$dir" && CUDA_VISIBLE_DEVICES="$gpu" "$PW_BIN" -nk 1 -input pw.in > pw.out 2>&1 )
    step_done "$dir/pw.out" \
      || { echo "[gpu$gpu] $name: pw.x FAILED (see $dir/pw.out)"; return 1; }
    echo "[gpu$gpu] $name: pw.x JOB DONE $(date +%H:%M:%S)"
  fi
  # step 2: pw.x 'bands' NSCF along the k-path (set I1 only)
  if [ -f "$dir/pw_bands.in" ]; then
    if step_done "$dir/pw_bands.out"; then
      echo "[gpu$gpu] $name: pw_bands.in already done, skip"
    else
      echo "[gpu$gpu] $name: pw.x bands start $(date +%H:%M:%S)"
      ( cd "$dir" && CUDA_VISIBLE_DEVICES="$gpu" "$PW_BIN" -nk 1 -input pw_bands.in > pw_bands.out 2>&1 )
      step_done "$dir/pw_bands.out" \
        || { echo "[gpu$gpu] $name: pw.x bands FAILED (see $dir/pw_bands.out)"; return 1; }
      echo "[gpu$gpu] $name: pw.x bands JOB DONE $(date +%H:%M:%S)"
    fi
    # steps 3-5: CPU post-processing (bands.x up/dw, projwfc.x)
    local in bin out
    for in in bands_up.in bands_dw.in projwfc.in; do
      [ -f "$dir/$in" ] || continue
      out="${in%.in}.out"
      bin="$BANDS_BIN"; [ "$in" = projwfc.in ] && bin="$PROJWFC_BIN"
      if step_done "$dir/$out"; then
        echo "[gpu$gpu] $name: $in already done, skip"; continue
      fi
      command -v "$bin" >/dev/null 2>&1 \
        || { echo "[gpu$gpu] $name: $bin missing, skipping $in"; return 1; }
      ( cd "$dir" && "$bin" < "$in" > "$out" 2>&1 )
      step_done "$dir/$out" \
        && echo "[gpu$gpu] $name: $in done" \
        || { echo "[gpu$gpu] $name: $in FAILED (see $dir/$out)"; return 1; }
    done
  fi
  echo "[gpu$gpu] $name: chain complete $(date +%H:%M:%S)"
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

list_jobs() {  # all set-I job names, manifest order (I1 chains first)
  python3 - <<'EOF'
import json
for j in json.load(open("manifest_bands.json"))["jobs"]:
    print(j["name"])
EOF
}

echo "== set I: fatband chains + relaxed-water + vdW jobs =="
list_jobs > .queue_bands
run_pool .queue_bands
[ -f STOP ] && { echo "stopped by STOP file"; exit 1; }

echo "== all bands-set jobs finished =="
echo "Summary:"
n_ok=0; n_bad=0
for d in jobs_bands/*/; do
  if job_done "$d"; then n_ok=$((n_ok+1)); else n_bad=$((n_bad+1)); echo "  incomplete: $d"; fi
done
echo "  $n_ok done, $n_bad incomplete"
echo "BANDS COMPLETE"
