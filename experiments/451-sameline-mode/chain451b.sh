#!/usr/bin/env bash
# 451b: apply the per-context-w patch on top of builder 1dcb29a3, run the tests, then the same two known-answer builds as
# 451: KA-1 unset == 446b's output sha; KA-2 VCC_SAMELINE (A = Nadig Jurkat programme, numeric w 0) == 446bja's. Any
# failure restores the byte copy of 1dcb29a3 and stops. The coordinator commits after reading the results.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/451-sameline-mode; LOG=$X/work/chain451b.log; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
OLD=1dcb29a382eebfa3
[ "$(shasum -a 256 $BU | cut -c1-16)" = "$OLD" ] || { log "builder is not 1dcb29a3 before the patch -- stopping"; exit 1; }
cp $BU $X/work/builder-1dcb29a3.py
back(){ cp $X/work/builder-1dcb29a3.py $BU; log "RESTORED 1dcb29a3: $(shasum -a 256 $BU | cut -c1-16) ($1)"; exit 1; }
uv run python $X/patch_builder451b.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: $OLD -> ${NEW:0:16}"
uv run python $X/test_sameline451.py > $X/work/test_sameline451b.log 2>&1 || back "451 unit tests failed"
log "451/451b unit tests: $(tail -1 $X/work/test_sameline451b.log)"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-451b.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-451b.log)"
grep -q "56 passed, 0 failed" $X/work/builder-tests-451b.log || back "builder suite not 56/0"
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t-451b.log 2>&1 || back "$t failed"; done
log "417 unit tests pass"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
R446B="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
R446BJA="$REPO/submissions/446bja-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
( unset VCC_SAMELINE VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $X/run_ka_451.sh ka1b "$R446B" "$NEW" ) >> $X/work/ka1b.nohup 2>&1
grep -q "KA PASS" $X/work/build_ka1b/runka1b.log || back "KA-1b (unset == 446b) failed"
log "KA-1b PASS: unset rebuilds 446b byte for byte"
KA=$X/work/sameline_ka_A.npz; KS=$(shasum -a 256 $KA | cut -d' ' -f1)
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" \
         VCC_SAMELINE="$KA:$KS:0"
  bash $X/run_ka_451.sh ka2b "$R446BJA" "$NEW" ) >> $X/work/ka2b.nohup 2>&1
grep -q "KA PASS" $X/work/build_ka2b/runka2b.log || back "KA-2b (numeric w 0 == 446bja) failed"
log "KA-2b PASS: VCC_SAMELINE numeric w 0 rebuilds 446bja byte for byte"
log "END builder ${NEW:0:16} (uncommitted; coordinator commits)"
