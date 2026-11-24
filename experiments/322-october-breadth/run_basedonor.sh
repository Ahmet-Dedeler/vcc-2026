#!/bin/sh
# One (panel, fold) per process: `prepare_main_source` enforces a PEAK-RSS cap and the
# operator's eigenvector matrix is 2.6 GiB resident, so a long-lived process trips it on
# the second fold.  264 used the same rule ("one PLOLO fold per process").
set -e
for p in abc real pess; do
  for f in k562_gwps orion_hct116 orion_hek293t; do
    echo "=== $p / $f"
    uv run python basedonor322.py --panels "$p" --folds "$f"
  done
done
