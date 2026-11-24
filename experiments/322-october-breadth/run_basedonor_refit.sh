#!/bin/sh
set -e
for p in abc real pess; do
  for f in k562_gwps orion_hct116 orion_hek293t; do
    uv run python basedonor322.py --panels "$p" --folds "$f" \
      --kolf2-root 'work/kolf2-refit/{panel}' --arm-suffix '_refit' \
      --out work/basedonor_refit322.json
  done
done
