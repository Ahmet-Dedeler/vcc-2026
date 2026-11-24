#!/usr/bin/env bash
# 2026-09-25 SECOND slot (rule fixed 01:30Z 09-24): after 433h's capture,
#   433h reads inside its stop bins vs 433a (nmae <= +0.001 raw; |d pds| <= 0.003 raw) -> 433k (433h with q 1.5 + 417k's re-level);
#   otherwise no automatic slot 2 (PI's call).
# Bins: PREREGISTRATION-0925.md. Waits (timestamp-capped, 12:00Z 09-25) for "CAPTURED own entry" in submit433h.log.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/433-coarse-mde-arm/work"; LOG="$W/submit433k.log"; S1="$W/submit433h.log"
RL="$REPO/experiments/428-zero-gene-audit/work/build_433k/run433k.log"
PKG="submissions/433k-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload-ctxlevel-mde010coarse-zr5.vcc"
NAME="the incubator adopted a goldfish"; REF=FCbBcvGxF1HOFPxPISUC   # 433a
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
stop(){ log "NO SLOT 2: $*"; exit 1; }
log "WAITING for 433h's capture (pid $$)"
until grep -q "CAPTURED own entry" "$S1" 2>/dev/null; do
  grep -q "NOT SUBMITTING" "$S1" 2>/dev/null && { log "slot 1 did not submit -- slot 2 is the PI's call; exiting"; exit 0; }
  [ -f "$W/NOGO-0925-s2" ] && stop "NOGO-0925-s2 present"
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609251200 ] && stop "433h never captured by 12:00Z 09-25"; sleep 30; done
ID="$(grep -o 'CAPTURED own entry [A-Za-z0-9]\{20\}' "$S1" | tail -1 | awk '{print $4}')"
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
ok = d["nmae"] <= 0.001 and abs(d["pds"]) <= 0.003
print(("433k" if ok else "stop"), json.dumps({k: round(v, 6) for k, v in d.items()}), f"d_overall {a['score_avg'] - b['score_avg']:+.7f} overall {a['score_avg']:.7f}")
PY
)"
log "rule: $RULE"
case "$RULE" in 433k*) ;; *) stop "$RULE";; esac
grep -q "^\[433k\] PASS" "$RL" || stop "433k config diff did not PASS"
grep -q "^\[433k\] CLIFF OK" "$RL" || stop "433k cliff not OK"
grep -q "novelty rc=0" "$RL" || stop "433k novelty not rc=0"
[ -f "$W/NOGO-0925-s2" ] && stop "NOGO-0925-s2 present"
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ -n "$SHA" ] && [ "$GOT" = "$SHA" ] || stop "433k sha mismatch got=$GOT want=$SHA"
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || stop "433k novelty check refused"
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID2="$(printf '%s\n' "$OUT" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID2" ] && ID2="$(printf '%s\n' "$OUT" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID2:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID2" ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID2" --note "433k slot 2 on 09-25 (433h + q 1.5 + 417k re-level): $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID2" --package "$PKG" --sha256 "$SHA" --note "433k slot 2 on 09-25" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID2 -> data/board-entries/$ID2.json" || log "WARN capture of $ID2 failed; re-run scripts/capture_board_entry.py $ID2"
fi
log "END"
