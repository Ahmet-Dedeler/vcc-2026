#!/usr/bin/env bash
# 2026-09-23 SECOND slot, rule AMENDED 21:00Z before any read (417 PREREGISTRATION.md addendum):
#   417a DE member outside the noise band (|dnmae|>0.010, |dreach|>0.012, |dfid|>0.003, |djac|>0.001 raw vs 410b) -> 410c
#   417a overall - 410b >= +0.006 -> 417b (q 1.5; needs its build gates, else 410c)
#   otherwise -> 417c (carrier gamma 0 on 417a; needs its build gates, else 410c)
# If slot 1 went to the 410c fallback (417a failed a gate), this exits and leaves slot 2 to the PI.
# work/NOGO-slot2 stops it at any point before the submit.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/417-bulk-sharpen/work"; LOG="$W/submit_slot2_0923.log"
PKG_B="submissions/417b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150.vcc"; NAME_B="my thesis advisor is a cat"
PKG_C="submissions/417c-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-carrierg0.vcc"; NAME_C="the mitochondria unionized"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
to410c(){ log "SLOT 2 -> 410c: $*"; exec bash "$REPO/experiments/410-pertscale-revisit/arm410c.sh"; }
log "WAITING for 417a's capture (pid $$)"
n=0; until grep -q "CAPTURED own entry" "$W/submit417a.log" 2>/dev/null; do
  grep -q "FALLBACK to 410c" "$W/submit417a.log" 2>/dev/null && { log "slot 1 fell back to 410c -- slot 2 is the PI's call; exiting"; exit 0; }
  [ -f "$W/NOGO-slot2" ] && { log "NOGO-slot2 present"; exit 1; }
  n=$((n+1)); [ $n -ge 2880 ] && { log "417a never captured in 24 h"; exit 1; }; sleep 30; done
ID="$(grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' "$W/submit417a.log" | tail -1 | awk '{print $4}')"
RULE="$(uv run python - "$ID" <<'PY'
import json, sys
def find(x):
    if isinstance(x, dict):
        if "score_avg" in x: return x
        for v in x.values():
            y = find(v)
            if y: return y
a = find(json.load(open(f"data/board-entries/{sys.argv[1]}.json")))
b = find(json.load(open("data/board-entries/4WgVi6FIE4iVw9PmRUP3.json")))
keys = {"nmae": "de_wilcoxon_lfc_nmae", "reach": "de_wilcoxon_direction_reach_raw", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
        "jac": "de_wilcoxon_sig_jaccard", "pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm"}
d = {k: a[v] - b[v] for k, v in keys.items()}
dov = a["score_avg"] - b["score_avg"]
leak = abs(d["nmae"]) > 0.010 or abs(d["reach"]) > 0.012 or abs(d["fid"]) > 0.003 or abs(d["jac"]) > 0.001
choice = "410c-leak" if leak else ("417b" if dov >= 0.006 else "417c")
print(choice, json.dumps({k: round(v, 6) for k, v in d.items()}), f"d_overall {dov:+.7f} overall {a['score_avg']:.7f}")
PY
)"
log "rule: $RULE"
case "$RULE" in
  417b*) ARM=417b; RLB="$W/build_b/run417b.log"; PKG="$PKG_B"; NAME="$NAME_B";;
  417c*) ARM=417c; RLB="$W/build_c/run417c.log"; PKG="$PKG_C"; NAME="$NAME_C";;
  *) to410c "$RULE";;
esac
log "slot 2 -> $ARM ($PKG as '$NAME'); waiting for its build gates"
n=0; until grep -q " DONE$" "$RLB" 2>/dev/null; do
  grep -qE "stopping|NO ARTIFACT|GATE never opened" "$RLB" 2>/dev/null && to410c "$ARM build stopped: $(grep -E 'stopping|NO ARTIFACT|GATE never' "$RLB" | tail -1)"
  n=$((n+1)); [ $n -ge 720 ] && to410c "$ARM build not done after 6 h"; sleep 30; done
grep -q "^\[$ARM\] PASS" "$RLB" || to410c "$ARM config diff did not PASS"
grep -q "^\[$ARM\] CLIFF OK" "$RLB" || to410c "$ARM cliff check not OK: $(grep "\[$ARM\] CLIFF" "$RLB" | tail -1)"
grep -q "novelty rc=0" "$RLB" || to410c "$ARM novelty not rc=0"
[ -f "$W/NOGO-slot2" ] && { log "NOGO-slot2 present"; exit 1; }
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
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "$ARM slot 2 on 09-23 ($([ $ARM = 417b ] && echo 'bulk sharpening q1.5 on 410b' || echo 'carrier gamma 0 on 417a')): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "$ARM slot 2 on 09-23" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
