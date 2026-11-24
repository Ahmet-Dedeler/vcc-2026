#!/usr/bin/env bash
# After the overnight chain (KA-LOAD PASS + dumps), build the stack LOAD arms for 09-24 (420's ladder, step 1):
#   417e = 417a (q 1.25) + VCC_BULK_LOAD(stack: H1 x2 + CD4 x1.5)   one variable vs 417a
#   417f = 417b (q 1.5)  + VCC_BULK_LOAD(stack)                      one variable vs 417b
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/chain_loadarms417.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
log "WAITING for the overnight chain (pid $$)"
n=0; until grep -qE "END$|STOP|never" "$W/chain_overnight417.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 1440 ] && { log "overnight chain never ended"; exit 1; }; sleep 30; done
grep -q "KA-LOAD PASS" "$W/chain_overnight417.log" || { log "no KA-LOAD PASS -- LOAD arms NOT built"; exit 1; }
D="$W/dump_h1x2_cd4x15.npz"; [ -s "$D" ] || { log "no stack dump -- LOAD arms NOT built"; exit 1; }
B="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
log "stack dump $(shasum -a 256 "$D" | cut -d' ' -f1); builder $B"
bash "$E/run_load_arm417.sh" 417e 1.25 "$D" "$REPO/submissions/417a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125.h5ad.meta.json" \
  the-cell-line-got-a-lawyer 417e-cap15-partial50-clean4q80-pertscale3-ballast-sharpen125-stackload "$B" >> "$LOG" 2>&1
log "417e: $(grep -E '^\[417e\] (PASS|FAIL|CLIFF)' "$W/build_417e/run417e.log" | tr '\n' ' ')"
bash "$E/run_load_arm417.sh" 417f 1.5 "$D" "$REPO/submissions/417b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150.h5ad.meta.json" \
  reviewer-two-was-a-cat-all-along 417f-cap15-partial50-clean4q80-pertscale3-ballast-sharpen150-stackload "$B" >> "$LOG" 2>&1
log "417f: $(grep -E '^\[417f\] (PASS|FAIL|CLIFF)' "$W/build_417f/run417f.log" | tr '\n' ' ')"
log "END"
