#!/bin/bash
# Row 1's second half: hash the artifact the driver built, then package it EXACTLY as
# runpreps229.sh packaged the shipped `h1union` (no --contexts on an A/B/C bundle), strictly
# after the build and never alongside one.
set -uo pipefail
cd "$(dirname "$0")/../../.."
H5=experiments/240-driver-fix/work/ka1/240-ka-h1union.h5ad
VCC=experiments/240-driver-fix/work/ka1/240-ka-h1union.vcc
while pgrep -f "make_submission_pm_cd4.py --bundle" > /dev/null; do sleep 30; done
echo "=== build finished $(date -u +%FT%TZ) ==="
[ -f "$H5" ] || { echo "NO ARTIFACT"; exit 1; }
echo -n "h5ad bytes  "; stat -f %z "$H5"
echo -n "h5ad sha256 "; shasum -a 256 "$H5" | awk '{print $1}'
echo "required    8df4dadd5a63b76c7ca3e3d57a511e7e52000b79d016d140897c55dd24c3ecf6 / 12648807800"
echo "=== prep $(date -u +%FT%TZ) ==="
/usr/bin/time -l vcc prep "$H5" -g data/challenge/gene_names.csv \
  --perts data/challenge/pert_counts.csv -o "$VCC" > experiments/240-driver-fix/work/ka1-prep.log 2>&1
echo "prep exit=$?"
if [ -f "$VCC" ]; then
  echo -n "vcc bytes   "; stat -f %z "$VCC"
  echo -n "vcc sha256  "; shasum -a 256 "$VCC" | awk '{print $1}'
  echo "required    f13ce40bca908af7cb85b8e84887af473aae4588410c238529139b7d1f245b1e / 2824785920"
fi
grep -E "maximum resident set size" experiments/240-driver-fix/work/ka1-prep.log | tail -1
echo "=== DONE $(date -u +%FT%TZ) ==="
