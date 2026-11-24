#!/usr/bin/env bash
# 472: the absolute prior gate `VCC_DE_PRIOR_GATE=le<N>` (patch_builder472.py). Serial after chain470b (itself after 471).
#  1. the builder must equal HEAD (471 committed or restored; not pinned); byte copy kept;
#  2. patch; tests: builder suite 56/0, 468, 472, 451, 417;
#  3. KA-le1: 458g's env with VCC_DE_PRIOR_GATE=le1 rebuilds 458g byte for byte (on A/B/C le1 selects q33's weak sets);
#  4. KA-1: every knob unset rebuilds 446b byte for byte.
# Any failure restores the pre-472 builder and stops. Never submits. The coordinator commits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/472-gate-absolute; LOG=$X/work/chain472.log; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain472.pid
log "armed (pid $$): waiting for chain470b to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|.*stopping.*)$" experiments/470-weak-target-prior/work/chain470b.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain470b\.sh|chain471\.sh|run_arm_generic\.sh|run_ka_471\.sh|run_dump471\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271200 ] && { log "470b never ended by 12:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
until [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271200 ] && { log "builder never equal to HEAD (471 uncommitted?) by 12:00Z 09-27 -- stopping"; exit 1; }
  log "builder != HEAD (471's patch not committed yet?) -- waiting"; sleep 300; done
OLD="$(shasum -a 256 $BU | cut -d' ' -f1)"; cp $BU $X/work/builder-pre472.py
back(){ cp $X/work/builder-pre472.py $BU; log "RESTORED ${OLD:0:16}: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
uv run python $X/patch_builder472.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${OLD:0:16} -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
for t in experiments/468-prior-gate/test_gate468.py $X/test_gate472.py experiments/451-sameline-mode/test_sameline451.py; do
  uv run python $t > $X/work/$(basename $t .py).log 2>&1 || back "$t failed"; done
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t.log 2>&1 || back "$t failed"; done
log "unit tests pass (468, 472, 451, 417)"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
R446B="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
R458G="$REPO/submissions/458g-cap15-partial50-arch1dprog-gate33-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
A=experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz; AS=$(shasum -a 256 $A | cut -d' ' -f1)
( unset VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_UNION_SOURCES
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" \
         VCC_SAMELINE="$REPO/$A:$AS:0" VCC_DE_PRIOR_GATE=le1
  bash $X/run_ka_472.sh kale1 "$R458G" "$NEW" ) >> $X/work/kale1.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_kale1/runkale1.log || back "KA-le1 (le1 == 458g) failed"
log "KA-le1 PASS: le1 rebuilds 458g byte for byte ($(/usr/bin/grep -o 'output sha256 [0-9a-f]* vs reference [0-9a-f]*' $X/work/build_kale1/runkale1.log))"
log "le1's gate lines: $(/usr/bin/grep -h '^\[prior_gate\]' $X/work/build_kale1/logs/build.log | tr '\n' ' ' | cut -c1-400)"
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_UNION_SOURCES
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $X/run_ka_472.sh ka1 "$R446B" "$NEW" ) >> $X/work/ka1.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_ka1/runka1.log || back "KA-1 (unset == 446b) failed"
log "KA-1 PASS: unset rebuilds 446b byte for byte"
log "END builder ${NEW:0:16} (uncommitted; the coordinator commits with KA-le1 and KA-1)"
