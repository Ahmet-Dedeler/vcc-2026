#!/bin/bash
# usage: step650.sh <label> <command string>   (runs in the 650 scratch worktree, never the live tree)
SP=/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/d7cf1fd9-dd97-40e0-ade7-66cf4f59473c/scratchpad
WT=$SP/wt650
LABEL=$1; shift
CMD="$*"
cd "$WT" || exit 99
export PYTHONPATH="$WT/src"
export UV_PROJECT_ENVIRONMENT=~/Code/virtual-cell-challenge-2026/.venv
export UV_NO_SYNC=1
OUT=$SP/650out/$(echo "$LABEL" | tr ' /' '__').log
mkdir -p $SP/650out
T0=$(date -u +%FT%TZ)
/usr/bin/time -l bash -c "$CMD" > "$OUT" 2>&1
RC=$?
T1=$(date -u +%FT%TZ)
python3 - "$LABEL" "$CMD" "$RC" "$T0" "$T1" "$OUT" >> $SP/650log.jsonl <<'PY'
import json, sys
label, cmd, rc, t0, t1, out = sys.argv[1:7]
txt = open(out, errors="replace").read()
peak = [l for l in txt.splitlines() if "maximum resident set size" in l]
print(json.dumps({"label": label, "cmd": cmd, "rc": int(rc), "t0": t0, "t1": t1, "log": out,
                  "peak_rss_bytes": int(peak[-1].split()[0]) if peak else None}))
PY
echo "[$LABEL] rc=$RC ($T0 -> $T1)"
grep -v "^ \+[0-9]\+ \+[a-z]" "$OUT" | grep -v "^$" | tail -${TAILN:-15}
