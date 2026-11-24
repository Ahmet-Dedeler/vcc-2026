#!/usr/bin/env bash
# 2026-09-24 SECOND slot (s2d, 09:45Z 09-23: 433a = 417e + COARSE MDE replaces 426a per experiments/433-coarse-mde-arm/PREREGISTRATION.md; s2c: loop cap by timestamp, 06:15Z 09-23), REVISED 04:30Z 09-23 before any read (docs/NEXT-SLOT-DECISION.md 09-24 revision; replaces
# 417-bulk-sharpen/arm0924_s2.sh, whose waiter pid 27781 was killed while it was still waiting for 417e's capture):
#   417e DE member outside the band (vs 417a) -> no automatic slot 2 (the DE channel must be pinned for the read)
#   otherwise -> 433a (417e + VCC_DE_MDE 0.1,0.1,0.3,coarse), one variable vs 417e whatever 417e's overall
#   (the re-level arms 417h/417g move to 09-25)
# Build-gate failures: no submission (PI's call).
# work/NOGO-slot2 stops it at any point before the submit.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/417-bulk-sharpen/work"; LOG="$W/submit_slot2_0924.log"; W6="$REPO/experiments/426-mde-magnitude/work"
PKG_M="submissions/433a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload-mde010coarse.vcc"; NAME_M="the centrifuge wrote a memoir"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
to410c(){ log "NO SLOT 2: $*"; exit 1; }
log "WAITING for 417e's capture (pid $$)"
n=0; until grep -q "CAPTURED own entry" "$W/submit417e.log" 2>/dev/null; do
  grep -q "NOT SUBMITTING" "$W/submit417e.log" 2>/dev/null && { log "slot 1 did not submit -- slot 2 is the PI's call; exiting"; exit 0; }
  [ -f "$W/NOGO-slot2-0924" ] && { log "NOGO-slot2-0924 present"; exit 1; }
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609241200 ] && { log "417e never captured by 12:00Z 09-24"; exit 1; }; sleep 30; done
ID="$(grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' "$W/submit417e.log" | tail -1 | awk '{print $4}')"
RULE="$(uv run python - "$ID" <<'PY'
import json, sys
def find(x):
    if isinstance(x, dict):
        if "score_avg" in x: return x
        for v in x.values():
            y = find(v)
            if y: return y
a = find(json.load(open(f"data/board-entries/{sys.argv[1]}.json")))
b = find(json.load(open("data/board-entries/T88GsCjp53Bmk8xQcJMS.json")))
keys = {"nmae": "de_wilcoxon_lfc_nmae", "reach": "de_wilcoxon_direction_reach_raw", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
        "jac": "de_wilcoxon_sig_jaccard", "pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm"}
d = {k: a[v] - b[v] for k, v in keys.items()}
dov = a["score_avg"] - b["score_avg"]
leak = abs(d["nmae"]) > 0.010 or abs(d["reach"]) > 0.012 or abs(d["fid"]) > 0.003 or abs(d["jac"]) > 0.001
choice = "leak" if leak else "433a"
print(choice, json.dumps({k: round(v, 6) for k, v in d.items()}), f"d_overall {dov:+.7f} overall {a['score_avg']:.7f}")
PY
)"
log "rule: $RULE"
case "$RULE" in
  433a*) ARM=433a; RLB="$REPO/experiments/428-zero-gene-audit/work/build_433a/run433a.log"; PKG="$PKG_M"; NAME="$NAME_M";;
  *) to410c "$RULE";;
esac
log "slot 2 -> $ARM ($PKG as '$NAME'); waiting for its build gates"
n=0; until grep -q " DONE$" "$RLB" 2>/dev/null; do
  grep -qE "stopping|NO ARTIFACT|GATE never opened" "$RLB" 2>/dev/null && to410c "$ARM build stopped: $(grep -E 'stopping|NO ARTIFACT|GATE never' "$RLB" | tail -1)"
  n=$((n+1)); [ $n -ge 720 ] && to410c "$ARM build not done after 6 h"; sleep 30; done
grep -q "^\[$ARM\] PASS" "$RLB" || to410c "$ARM config diff did not PASS"
grep -q "^\[$ARM\] CLIFF OK" "$RLB" || to410c "$ARM cliff check not OK: $(grep "\[$ARM\] CLIFF" "$RLB" | tail -1)"
grep -q "novelty rc=0" "$RLB" || to410c "$ARM novelty not rc=0"
[ -f "$W/NOGO-slot2-0924" ] && { log "NOGO-slot2-0924 present"; exit 1; }
SHA="$(grep 'package sha256:' "$RLB" | tail -1 | awk '{print $NF}')"; GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ -n "$SHA" ] && [ "$GOT" = "$SHA" ] || to410c "$ARM sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || to410c "$ARM novelty check refused"
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID="$(printf '%s\n' "$OUT" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID" ] && ID="$(printf '%s\n' "$OUT" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID" ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "$ARM slot 2 on 09-24 (count-balanced MDE DE magnitude on 417e): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "$ARM slot 2 on 09-24" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
