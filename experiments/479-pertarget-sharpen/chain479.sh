#!/usr/bin/env bash
# 479: per-target bulk sharpening (VCC_BULK_SHARPEN_Q=tercile:/stercile:, patch_builder479.py). Waits for chain471c's
# PROCESS to exit (AGENTS §6.2) and the builder to equal HEAD. Then:
#   tests: builder suite 56/0, 468, 472, 451, 417, 479 (parse + 446bj2's dump terciles);
#   KA-S (scalar q 1.75, 446bj2's env) == 446bj2's output 711ee43b; KA-T (tercile:1.75,1.75,1.75) == 711ee43b;
#   KA-ST (stercile:1.75,1.75,1.75) == 711ee43b.
# Any failure restores the pre-479 builder. Never submits. The arms (P2g, S_P2) are built by the coordinator afterwards.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/479-pertarget-sharpen; LOG=$X/work/chain479.log; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain479.pid
log "armed (pid $$): waiting for chain471c's process to exit"
until ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain471c\.sh|chain475b\.sh|run_arm_generic\.sh|def_day\.sh|run_ka_4/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271600 ] && { log "471c still running at 16:00Z -- stopping"; log "END"; exit 1; }; sleep 60; done
until [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271600 ] && { log "builder never equal to HEAD by 16:00Z -- stopping"; log "END"; exit 1; }
  log "builder != HEAD -- waiting"; sleep 300; done
OLD="$(shasum -a 256 $BU | cut -d' ' -f1)"; cp $BU $X/work/builder-pre479.py
back(){ cp $X/work/builder-pre479.py $BU; log "RESTORED ${OLD:0:16}: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
uv run python $X/patch_builder479.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${OLD:0:16} -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
for t in $X/test_pt479.py experiments/468-prior-gate/test_gate468.py experiments/472-gate-absolute/test_gate472.py experiments/451-sameline-mode/test_sameline451.py; do
  uv run python $t > $X/work/$(basename $t .py).log 2>&1 || back "$t failed"; done
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t.log 2>&1 || back "$t failed"; done
log "unit tests pass (479, 468, 472, 451, 417)"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
[ "${DS:0:16}" = 00ccfcba1b67fa24 ] || back "LODO dump sha ${DS:0:16} unexpected"
R="$REPO/$(ls submissions/446bj2-*.h5ad.meta.json | head -1)"
for tag in kas:1.75 kat:tercile:1.75,1.75,1.75 kast:stercile:1.75,1.75,1.75; do
  T="${tag%%:*}"; QQ="${tag#*:}"
  ( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP
    export KA_Q="$QQ" VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.0962,B=1.1461,C=1.1681"
    bash $X/run_ka_479.sh "$T" "$R" "$NEW" ) >> $X/work/$T.nohup 2>&1
  /usr/bin/grep -q "KA PASS" $X/work/build_$T/run$T.log || back "KA-$T ($QQ == 446bj2) failed"
  log "KA-$T PASS: VCC_BULK_SHARPEN_Q=$QQ rebuilds 446bj2 byte for byte ($(/usr/bin/grep -o 'output sha256 [0-9a-f]* vs reference [0-9a-f]*' $X/work/build_$T/run$T.log))"
done
log "END builder ${NEW:0:16} (uncommitted; the coordinator commits with KA-S/KA-T/KA-ST)"
