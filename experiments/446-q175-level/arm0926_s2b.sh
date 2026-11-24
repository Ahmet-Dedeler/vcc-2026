#!/usr/bin/env bash
# 2026-09-26 SECOND slot (RE-ARMED 2026-09-25 03:55Z, before any 09-26 entry): 446bja = 446b with context A's block prior
# from Nadig Jurkat (clean4q80ja), read vs 446b whatever 446b scores (450/PREREGISTRATION-446bja.md). The LODO-key rider
# moves to 09-27. Fallback if 446bja's build is missing or fails a gate: the original branch-A rule (446bj2 if 446b - 433k
# >= +0.001 with pds >= -0.0005 and DE pinned, else 433kj). Timestamp-capped (12:00Z 09-26).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/446-q175-level/work"; LOG="$W/submit0926s2.log"; S1="$W/submit446b.log"
REF=T76JWD3tOajk5IYvzwUW   # 433k
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
stop(){ log "NO SLOT 2: $*"; exit 1; }
log "RE-ARMED for 446bja (pid $$); WAITING for 446b's capture"
until grep -q "CAPTURED own entry" "$S1" 2>/dev/null; do
  grep -q "NOT SUBMITTING" "$S1" 2>/dev/null && { log "slot 1 did not submit -- slot 2 is the PI's call; exiting"; exit 0; }
  [ -f "$W/NOGO-0926-s2" ] && stop "NOGO-0926-s2 present"
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609261200 ] && stop "446b never captured by 12:00Z 09-26"; sleep 30; done
ID="$(grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' "$S1" | tail -1 | awk '{print $4}')"
JRL="$REPO/experiments/450-same-line-programme/work/build_446bja/run446bja.log"
JPKG="submissions/446bja-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.vcc"
if grep -q "^\[446bja\] PASS" "$JRL" 2>/dev/null && grep -q "^\[446bja\] CLIFF OK" "$JRL" && grep -q "novelty rc=0" "$JRL" && [ -s "$JPKG" ]; then
  RULE="446bja (preregistered, vs 446b $ID whatever it scores)"
else
  log "446bja not available or not gated -- falling back to the branch-A key-rider rule"
  RULE="$(uv run python - "$ID" "$REF" <<'PY'
import json, sys
def find(x):
    if isinstance(x, dict):
        if "score_avg" in x: return x
        for v in x.values():
            y = find(v)
            if y: return y
a = find(json.load(open(f"data/board-entries/{sys.argv[1]}.json"))); b = find(json.load(open(f"data/board-entries/{sys.argv[2]}.json")))
keys = {"nmae": "de_wilcoxon_lfc_nmae", "reach": "de_wilcoxon_direction_reach_raw", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
        "jac": "de_wilcoxon_sig_jaccard", "pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm"}
d = {k: a[v] - b[v] for k, v in keys.items()}
dov = a["score_avg"] - b["score_avg"]
pinned = abs(d["nmae"]) <= 0.001 and abs(d["fid"]) <= 0.001 and abs(d["jac"]) <= 0.001 and abs(d["reach"]) <= 0.01
pick = "446bj2" if (dov >= 0.001 and d["pds"] >= -0.0005 and pinned) else "433kj"
print(pick, json.dumps({k: round(v, 6) for k, v in d.items()}), f"d_overall {dov:+.7f} pinned {pinned}")
PY
)"
fi
log "rule: $RULE"
case "$RULE" in
  446bja*) TAG=446bja; RL="$JRL"; PKG="$JPKG"; NAME="the vortex mixer wrote a sonnet"
           NOTE="446bja slot 2 on 09-26 (446b + clean4q80ja: A's block prior from Nadig Jurkat)";;
  446bj2*) TAG=446bj2; RL="$REPO/experiments/428-zero-gene-audit/work/build_446bj2/run446bj2.log"; NAME="the petri dish hired an agent"
           PKG="submissions/446bj2-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5.vcc"
           NOTE="446bj2 slot 2 on 09-26 (446b + LODO key at matched rho)";;
  433kj*)  TAG=433kj; RL="$REPO/experiments/428-zero-gene-audit/work/build_433kj/run433kj.log"; NAME="the petri dish hired an agent"
           PKG="submissions/433kj-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stacklodo05load-ctxlevel-mde010coarse-zr5.vcc"
           NOTE="433kj slot 2 on 09-26 (433k + LODO key)";;
  *) stop "rule did not resolve: $RULE";;
esac
grep -q "^\[$TAG\] PASS" "$RL" || stop "$TAG config diff did not PASS"
grep -q "^\[$TAG\] CLIFF OK" "$RL" || stop "$TAG cliff not OK"
grep -q "novelty rc=0" "$RL" || stop "$TAG novelty not rc=0"
[ -f "$W/NOGO-0926-s2" ] && stop "NOGO-0926-s2 present"
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ -n "$SHA" ] && [ "$GOT" = "$SHA" ] || stop "$TAG sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || stop "$TAG novelty check refused"
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID2="$(printf '%s\n' "$OUT" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID2" ] && ID2="$(printf '%s\n' "$OUT" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID2:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID2" ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID2" --note "$NOTE: $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID2" --package "$PKG" --sha256 "$SHA" --note "$NOTE" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID2 -> data/board-entries/$ID2.json" || log "WARN capture of $ID2 failed; re-run scripts/capture_board_entry.py $ID2"
fi
log "END"
