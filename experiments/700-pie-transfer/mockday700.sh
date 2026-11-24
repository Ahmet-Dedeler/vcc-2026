#!/usr/bin/env bash
# 700: rehearse dumphook700.sh on the A/B/C bundle as a mock D/E/F day (identity from 041's known match), timed, through
# gate602 so it never runs beside a production build. Compares the hook's output dump with 702's offline-built dump
# (same w, seen rule and construction; contexts named by accession here, so A's text is Jurkat CVCL_0065, not E6.1).
#   nohup bash experiments/700-pie-transfer/mockday700.sh > experiments/700-pie-transfer/work/mockday/run.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/700-pie-transfer; W=$X/work/mockday; mkdir -p "$W/out"
DUMP=experiments/495-lodo-cap6/work/dump_h1x2_cd4x15_lodo05hi6.npz
echo "$(date -u +%FT%TZ) queued"
GATE602_QUIET=1 /usr/bin/time -l bash experiments/602-council-gate/gate602.sh 700-mockday 20 -- \
  bash $X/dumphook700.sh "$REPO/data/challenge" "$REPO/$W" "$REPO/$DUMP" "$REPO/$W/identity.json" "$REPO/$W/out" > "$W/dumphook.log" 2>&1
echo "$(date -u +%FT%TZ) hook rc=$?"; tail -5 "$W/dumphook.log"
"$X/.venv/bin/python" - "$W/out/$(basename $DUMP)" "$X/../702-pie-bulk-arm/work/dump_pie702_w025_seen.npz" <<'PY'
import sys, numpy as np
a, b = np.load(sys.argv[1]), np.load(sys.argv[2])
for k in ("perts", "genes", "env", "carriers"):
    print(k, "identical" if np.array_equal(a[k], b[k]) else "DIFFERENT")
x, y = a["bulk"], b["bulk"]
cos = (x * y).sum(1) / np.maximum(np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1), 1e-300)
print("bulk rows cos(hook, 702 offline): median %.6f min %.6f" % (np.median(cos), cos.min()))
PY
echo "$(date -u +%FT%TZ) END"
