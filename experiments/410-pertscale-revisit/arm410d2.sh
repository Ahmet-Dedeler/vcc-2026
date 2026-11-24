#!/usr/bin/env bash
# 410d2 -- wait for run410d2_build.sh to pass its gates (config diff PASS, novelty rc=0, package sha), then submit (2026-09-22 has
# both slots unspent and nothing in flight), record the ledger with the "entry <id>"-anchored parser, capture the record.
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
W="$REPO/experiments/410-pertscale-revisit/work"; RL="$W/build_d2/run410d2.log"; LOG="$W/submit410d2.log"
PKG="submissions/410d2-cap15-partial50-clean4q80-pertscaleinf-ballast-genegainrg.vcc"; NAME="a sourdough starter achieved sentience"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
log "WAITING for 410 build DONE (pid $$)"
n=0; until grep -q " DONE$" "$RL" 2>/dev/null; do
  if grep -qE "stopping|NO ARTIFACT|GATE never opened" "$RL" 2>/dev/null; then log "BUILD STOPPED: $(grep -E 'stopping|NO ARTIFACT|GATE never' "$RL" | tail -1)"; exit 1; fi
  n=$((n+1)); [ $n -ge 720 ] && { log "build not done after 6 h"; exit 1; }; sleep 30; done
grep -q "^\[410d2\] PASS" "$RL" || { log "config diff did not PASS -- not submitting"; exit 1; }
# gates added for 410b: 410 must be terminal (captured) and the PI must have written work/GO-410b after reading H1's ballast rows
n=0; until grep -q "CAPTURED own entry" "$W/submit410c.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 1440 ] && { log "410c never captured -- not submitting"; exit 1; }; sleep 30; done
log "410c is captured; applying the preregistered GO-410d rule (docs/NEXT-SLOT-DECISION.md, 2026-09-23 plan)"
[ -f "$W/NOGO-410d2" ] && { log "NOGO-410d2 present -- not submitting"; exit 1; }
if [ ! -f "$W/GO-410d2" ]; then
  C_ID="$(grep -o 'entry id: [A-Za-z0-9]\{20\}' "$W/submit410c.log" | tail -1 | awk '{print $3}')"
  RULE="$(uv run python - "$C_ID" <<'PY'
import json, sys
def find(x):
    if isinstance(x, dict):
        if "score_avg" in x: return x
        for v in x.values():
            y = find(v)
            if y: return y
c = find(json.load(open(f"data/board-entries/{sys.argv[1]}.json")))
b = find(json.load(open("data/board-entries/4WgVi6FIE4iVw9PmRUP3.json")))
d = {k: c[k] - b[k] for k in ("de_wilcoxon_lfc_nmae", "de_wilcoxon_direction_reach_raw", "pds_cosine")}
ok = abs(d["de_wilcoxon_lfc_nmae"]) <= 0.010 and abs(d["de_wilcoxon_direction_reach_raw"]) <= 0.012 and d["pds_cosine"] >= -0.010
print(("GO" if ok else "NOGO") + " " + json.dumps({k: round(v, 6) for k, v in d.items()}) + f" overall {c['score_avg']:.7f}")
PY
)"
  log "rule: $RULE"
  case "$RULE" in GO*) echo "$(date -u +%FT%TZ) auto $RULE" > "$W/GO-410d2";; *) log "rule says NOGO -- not submitting"; exit 1;; esac
fi
n=0; until [ "$(date -u +%Y%m%d%H%M)" -ge 202609230001 ]; do n=$((n+1)); [ $n -ge 1440 ] && exit 1; sleep 20; done
log "GO-410d2: $(cat "$W/GO-410d2")"
grep -q "novelty rc=0" "$RL" || { log "novelty not rc=0 -- not submitting"; exit 1; }
SHA="$(grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; [ -n "$SHA" ] || { log "no package sha"; exit 1; }
GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ "$GOT" = "$SHA" ] || { log "ERROR sha mismatch got=$GOT want=$SHA -- NOT submitting"; exit 1; }
WHO="$(vcc --json whoami 2>&1)"; log "whoami: $(echo "$WHO" | head -c 300)"
echo "$WHO" | grep -q '"can_submit": true' || { log "ERROR account cannot submit -- NOT submitting"; exit 1; }
uv run python scripts/check_submission_novelty.py check "$PKG" >> "$LOG" 2>&1 || { log "ERROR novelty check refused -- NOT submitting"; exit 1; }
log "SUBMITTING $PKG as '$NAME' (description empty) sha $SHA"
OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?
if [ $RC -ne 0 ] && ! echo "$OUT" | grep -qi "entry\|scor"; then log "submit attempt 1 failed rc=$RC (retrying once in 90 s): $(echo "$OUT" | tail -c 200)"; sleep 90; OUT="$(vcc submit "$PKG" -m "$NAME" --wait 2>&1)"; RC=$?; fi
log "SUBMIT rc=$RC"; printf '%s\n' "$OUT" >> "$LOG"
[ $RC -eq 0 ] && log "SUBMITTED ok" || log "ERROR submit failed rc=$RC"
ID="$(printf '%s\n' "$OUT" | grep -o 'entry [A-Za-z0-9]\{20\}' | head -1 | awk '{print $2}')"
[ -z "$ID" ] && ID="$(printf '%s\n' "$OUT" | grep -o '\b[A-Za-z0-9]\{20\}\b' | grep -v -E '^[0-9]+$' | grep -v -E '^[a-z0-9]+$' | head -1)"
log "entry id: ${ID:-<none>}"
if [ $RC -eq 0 ] && [ -n "$ID" ]; then
  uv run python scripts/check_submission_novelty.py record "$PKG" --entry-id "$ID" --note "410d2 228 per-gene gain on the bulk, render-gated energy match, on 410c: $NAME" >> "$LOG" 2>&1
  uv run python scripts/capture_board_entry.py "$ID" --package "$PKG" --sha256 "$SHA" --note "410d2 render-gated per-gene bulk gain on 410c" --wait --attempts 40 --interval 30 >> "$LOG" 2>&1 \
    && log "CAPTURED own entry $ID -> data/board-entries/$ID.json" || log "WARN capture of $ID failed; re-run scripts/capture_board_entry.py $ID"
fi
log "END"
