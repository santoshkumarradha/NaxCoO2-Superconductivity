#!/bin/bash
# Set D: convergence + sensitivity checks. Run AFTER the main batch (uses its inputs as templates).
# Derives variant inputs from existing jobs via sed; runs them sequentially on GPU 0.
set -u
export PATH=/usr/local/qe/bin:/usr/local/openmpi/bin:/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib:/usr/local/cuda/lib64:/usr/local/fftw/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/24.7/compilers/lib:${LD_LIBRARY_PATH:-}
export OMP_NUM_THREADS=14
export PMIX_MCA_gds=hash
cd "$(dirname "$0")"
mkdir -p checks
run() { d=checks/$1; mkdir -p $d; cp jobs/$2/pw.in $d/; shift 2; for s in "$@"; do sed -i "$s" $d/pw.in; done
  rm -rf /tmp/pmix* /tmp/ompi.* 2>/dev/null; (cd $d && pw.x -nk 1 -input pw.in > pw.out 2>&1); grep -q "JOB DONE" $d/pw.out && echo "OK $d" || echo "FAIL $d"; }
# ecut 60->80, dense k, small smearing on a mid-well point
run ecut80    Na_1x1_c6.9_d0.50 "s/ecutwfc.*/ecutwfc = 80/; s/ecutrho.*/ecutrho = 640/"
run kdense    Na_1x1_c6.9_d0.50 "s/8 8 3/12 12 4/"
run smear01   Na_1x1_c6.9_d0.50 "s/degauss.*/degauss = 0.01/"
# spin off vs on at the two anchor spacings, well center and off-center
for j in Na_1x1_c6.9_d0.00 Na_1x1_c6.9_d0.50 Na_1x1_c9.9_d0.00 Na_1x1_c9.9_d0.75; do
  run nspin1_$j $j "s/nspin.*/nspin = 1/; /starting_magnetization/d"
done
# z_O sensitivity: +-0.06 A around 0.96 at c=6.9, delta=0.5 (fractional zo = zO/c)
for z in 0.90 1.02; do
  d=checks/zO${z}_c6.9; mkdir -p $d; cp jobs/Na_1x1_c6.9_d0.50/pw.in $d/
  python3 - "$d/pw.in" "$z" <<'PY'
import sys
p, z = sys.argv[1], float(sys.argv[2]); c = 6.9
out = []
for l in open(p):
    t = l.split()
    if len(t) == 4 and t[0] == "O":
        zc = float(t[3])
        zc = z/c if zc < 0.5 else 1 - z/c
        l = f"  O  {t[1]}  {t[2]}  {zc:.10f}\n"
    out.append(l)
open(p, "w").writelines(out)
PY
  rm -rf /tmp/pmix* /tmp/ompi.* 2>/dev/null
  (cd $d && pw.x -nk 1 -input pw.in > pw.out 2>&1); grep -q "JOB DONE" $d/pw.out && echo "OK $d" || echo "FAIL $d"
done
echo "SET-D COMPLETE"
