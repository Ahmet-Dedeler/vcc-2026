#!/usr/bin/env bash
# 451: after 446bja is built (chain450ja END), patch the builder with VCC_SAMELINE, run the tests, then two known-answer
# builds: KA-1 unset == 446b's output sha; KA-2 VCC_SAMELINE (A = Nadig Jurkat programme, w 0) == 446bja's output sha.
# Any failure restores the byte copy of b6c3f131 and stops. The coordinator commits after reading the results.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/451-sameline-mode; LOG=$X/work/chain451.log; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q " END$" experiments/450-same-line-programme/work/chain450ja.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain450ja\.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609252000 ] && { log "chain450ja never ended by 20:00Z"; exit 1; }; sleep 30; done
OLD=b6c3f13141e707012588cc354f62758889db9c38c3757f9765ebea37a8a733f8
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$OLD" ] || { log "builder is not b6c3f131 before the patch -- stopping"; exit 1; }
cp $BU $X/work/builder-b6c3f131.py
back(){ cp $X/work/builder-b6c3f131.py $BU; log "RESTORED b6c3f131: $(shasum -a 256 $BU | cut -c1-16) ($1)"; exit 1; }
uv run python $X/patch_builder451.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${OLD:0:16} -> ${NEW:0:16}"
uv run python $X/test_sameline451.py > $X/work/test_sameline451.log 2>&1 || back "451 unit tests failed"
log "451 unit tests: $(tail -1 $X/work/test_sameline451.log)"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1   # the suite's runner (skips the one test stale since 481cd092)
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t.log 2>&1 || back "$t failed"; done
log "417 unit tests pass"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
R446B="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
R446BJA="$REPO/submissions/446bja-cap15-partial50-clean4q80ja-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
[ -s "$R446BJA" ] || back "446bja's meta is missing"
( unset VCC_SAMELINE VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $X/run_ka_451.sh ka1 "$R446B" "$NEW" ) >> $X/work/ka1.nohup 2>&1
grep -q "KA PASS" $X/work/build_ka1/runka1.log || back "KA-1 (unset == 446b) failed"
log "KA-1 PASS: unset rebuilds 446b byte for byte"
KA=$X/work/sameline_ka_A.npz; KS=$(shasum -a 256 $KA | cut -d' ' -f1)
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" \
         VCC_SAMELINE="$KA:$KS:0"
  bash $X/run_ka_451.sh ka2 "$R446BJA" "$NEW" ) >> $X/work/ka2.nohup 2>&1
grep -q "KA PASS" $X/work/build_ka2/runka2.log || back "KA-2 (sameline programme == 446bja) failed"
log "KA-2 PASS: VCC_SAMELINE (A = Nadig Jurkat programme, w 0) rebuilds 446bja byte for byte"
log "END builder ${NEW:0:16} (uncommitted; coordinator commits)"
