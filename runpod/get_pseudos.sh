#!/usr/bin/env bash
# Download PAW pseudopotentials (pslibrary via pseudopotentials.quantum-espresso.org).
# NOTE: Co.pbe-spn-kjpaw_psl.1.0.0.UPF is not hosted upstream (404); the
# closest available pslibrary PAW variant, Co.pbe-spn-kjpaw_psl.0.3.1.UPF
# (same spn semicore flavour), is used instead. Verified 2026-07-11.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p pseudo
BASE="https://pseudopotentials.quantum-espresso.org/upf_files"
FILES=(
  Co.pbe-spn-kjpaw_psl.0.3.1.UPF
  O.pbe-n-kjpaw_psl.1.0.0.UPF
  Na.pbe-spn-kjpaw_psl.1.0.0.UPF
  Li.pbe-s-kjpaw_psl.1.0.0.UPF
  # --- extra stage (sets E-H); filenames verified against upstream 2026-07-12,
  #     all return HTTP 200. K/Se use the semicore (spn/dn) PAW flavours to
  #     match the Co.spn / Na.spn convention above (K z_val=9, Se z_val=16,
  #     H z_val=1).
  K.pbe-spn-kjpaw_psl.1.0.0.UPF
  Se.pbe-dn-kjpaw_psl.1.0.0.UPF
  H.pbe-kjpaw_psl.1.0.0.UPF
)
for f in "${FILES[@]}"; do
  if [ -s "pseudo/$f" ]; then
    echo "pseudo/$f already present"
    continue
  fi
  echo "downloading $f"
  curl -fSL --retry 3 -o "pseudo/$f" "$BASE/$f"
  grep -q "PP_HEADER" "pseudo/$f" || { echo "ERROR: $f looks corrupt"; exit 1; }
done
echo "pseudopotentials ready in ./pseudo"
