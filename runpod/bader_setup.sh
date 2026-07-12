#!/usr/bin/env bash
# Bader charge analysis for all finished SCF jobs.
#   1. Download the Henkelman-group bader Linux binary into tools/.
#   2. Per finished SCF job: pp.x writes the PAW all-electron valence density
#      (plot_num=17) and the valence+core density (plot_num=21) as cubes.
#   3. bader partitions on the total density (-ref) and integrates the valence
#      density -> jobs/<name>/ACF.dat  (Na/Li row gives the donated charge).
# Idempotent: jobs that already have ACF.dat are skipped.
set -uo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT"

PP_BIN="${PP_BIN:-pp.x}"
command -v "$PP_BIN" >/dev/null || { echo "ERROR: $PP_BIN not in PATH"; exit 1; }

BADER=tools/bader
if [ ! -x "$BADER" ]; then
  mkdir -p tools
  echo "downloading Henkelman bader binary..."
  URL="https://theory.cm.utexas.edu/henkelman/code/bader/download/bader_lnx_64.tar.gz"
  curl -fSL --retry 3 -o tools/bader.tar.gz "$URL"
  tar -xzf tools/bader.tar.gz -C tools
  chmod +x "$BADER"
fi
"$BADER" --help >/dev/null 2>&1 || true

write_pp() {  # $1 = dir, $2 = plot_num, $3 = cube name, $4 = pp.in name
  cat > "$1/$4" <<EOF
&INPUTPP
  prefix = 'pw'
  outdir = './out'
  plot_num = $2
  filplot = '$3.pp'
/
&PLOT
  iflag = 3
  output_format = 6
  fileout = '$3.cube'
/
EOF
}

n=0
for d in jobs/*/; do
  name=$(basename "$d")
  case "$name" in nscf_*) continue ;; esac                    # SCF jobs only
  [ -f "$d/pw.out" ] && tail -n 20 "$d/pw.out" | grep -q "JOB DONE" || continue
  [ -f "$d/ACF.dat" ] && { echo "$name: ACF.dat exists, skip"; continue; }
  echo "$name: pp.x + bader"
  write_pp "$d" 17 valence pp_val.in
  write_pp "$d" 21 total   pp_all.in
  ( cd "$d" \
    && "$PP_BIN" -input pp_val.in > pp_val.out 2>&1 \
    && "$PP_BIN" -input pp_all.in > pp_all.out 2>&1 \
    && "$ROOT/$BADER" valence.cube -ref total.cube > bader.out 2>&1 \
    && rm -f valence.pp total.pp ) \
    && { echo "$name: ACF.dat written"; n=$((n+1)); } \
    || echo "$name: FAILED (see $d/pp_*.out, $d/bader.out)"
rm -f "$d"/*.cube "$d"/out/*.wfc* 2>/dev/null
done
echo "bader analysis complete for $n new job(s)"
