# NaxCoO2 interlayer-2DEG / anharmonic-rattler — RunPod QE package

Self-contained Quantum ESPRESSO job set testing whether Na_xCoO2.yH2O
superconductivity arises from an interlayer 2DEG coupled to an anharmonic Na
double-well mode. Water is replaced by its structural effect: the CoO2
plane-to-plane spacing `c` in {5.5, 6.9, 8.4, 9.9} A is the control parameter.

## Job sets

| Set | What | Jobs |
|-----|------|------|
| A | 1x1 cell E(delta) scans, Na and Li, 4 c x 6 delta each | 48 SCF |
| B | sqrt3 x sqrt3 supercell (x = 1/3, 10 atoms), Na, delta in {0, 0.5} | 8 SCF |
| C | dense-k NSCF (24x24x8) + dos.x at the E(delta) minimum per c: Na 1x1, Li 1x1, Na x=1/3 | 12 NSCF+DOS |

Set C inputs are generated automatically **after** the SCF phase (the minimum
geometry is not known beforehand); `run_all.sh` handles this.

Settings: PBE, PAW (pslibrary), ecutwfc 60 Ry / ecutrho 480 Ry, Marzari-
Vanderbilt smearing 0.02 Ry, nspin=2 with starting Co moment 0.3,
conv_thr 1e-8, verbosity='high' (Gamma eigenvalues are parsed for the alkali
band deformation potential). **Pseudo note:** `Co.pbe-spn-kjpaw_psl.1.0.0.UPF`
is not hosted on pseudopotentials.quantum-espresso.org (verified 2026-07-11);
the closest available pslibrary PAW variant `Co.pbe-spn-kjpaw_psl.0.3.1.UPF`
(same spn semicore flavour) is used instead.

## RunPod steps

1. **Launch a pod** with the NVIDIA NGC Quantum ESPRESSO image:
   `nvcr.io/hpc/quantum_espresso:qe-7.3.1` <-- image tag being verified separately; typically
   `nvcr.io/hpc/quantum_espresso:<tag>`. Any pod with `pw.x`, `dos.x`, `pp.x`
   (GPU builds) in PATH works. Attach a volume (>= 50 GB) mounted at
   `/workspace`. 1-8 GPUs; the launcher uses all of them.

2. **Sync the package up** (from this directory):

       rsync -av --exclude 'jobs/*/out' ./ root@<POD_IP>:/workspace/naxcoo2/ -e "ssh -p <PORT>"

3. **Run everything** on the pod:

       cd /workspace/naxcoo2
       python3 generate_inputs.py        # writes 56 SCF inputs + manifest.json
       bash run_all.sh                   # phases: SCF -> gen NSCF -> NSCF+DOS
       bash bader_setup.sh               # pp.x cubes + Bader ACF.dat per SCF

   `run_all.sh` detects GPUs via `nvidia-smi`, runs one `pw.x -nk 1` per GPU
   (pinned with `CUDA_VISIBLE_DEVICES`) from a flock'd bash job pool, and is
   idempotent — jobs whose `pw.out` ends in `JOB DONE` are skipped, so just
   re-run it after any interruption.

   **Kill switch:** `touch /workspace/naxcoo2/STOP` — each worker exits after
   its current job. Remove STOP and re-run to resume.

4. **Sync results back** (cubes and wavefunctions excluded — only text output
   is needed locally):

       rsync -av --exclude 'out' --exclude '*.cube' --exclude '*.pp' \
             root@<POD_IP>:/workspace/naxcoo2/ ./ -e "ssh -p <PORT>"

5. **Post-process locally** (needs numpy, scipy, matplotlib):

       python3 postprocess.py --selftest   # verify fit + Schroedinger solver
       python3 postprocess.py              # -> results.json + triple_plot.pdf

   `triple_plot.pdf`: alpha(c) (double-well onset where alpha < 0),
   n_2DEG(c) from the Bader alkali charge, and Allen-Dynes Tc(c) with the
   experimental anchors marked (6.9 A non-SC, 9.9 A SC Tc = 4.5 K).

## Rough timing (one A100/H100 per job)

- Set A SCF (4 atoms, 8x8x3, spin): ~3-10 min/job -> 48 jobs / N_GPU.
- Set B SCF (10 atoms, 6x6x3): ~15-40 min/job.
- Set C NSCF (24x24x8): ~10-30 min for 1x1, up to ~1-2 h for the supercell
  (densest-k jobs; the single most expensive part). dos.x and pp.x/bader are
  seconds to minutes.
- Total: a few GPU-hours on 1 GPU, well under an hour of wall time on 8.

Budget guard: everything restarts cleanly, so use spot/community pods freely.
