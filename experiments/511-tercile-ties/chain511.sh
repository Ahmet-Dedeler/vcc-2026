#!/usr/bin/env bash
# 511: the tie-rule fix in sharpen_q_by_tercile (501 section 6), after the whole 09-27 queue (serial after chain510).
#  1. patch_builder511.py on the builder (expected 24b95dec -> 510f96bd); builder suite 56/0; 479/468/472/451/417 tests;
#     test_ties511.py (old vs new function: q vectors EQUAL on 446bj2's LODO dump and 495's cap-6 dump in both modes; the
#     synthetic tie cases fixed);
#  2. KA-SP2: 479sp2's own configuration (stercile:2.5,2.0,1.75 at its round-2 factors) rebuilds 479sp2's output byte for byte.
# Any failure restores the pre-511 builder. commit_after_chain.sh commits it (spec commit511.spec). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/511-tercile-ties; LOG=$X/work/chain511.log; BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain511.pid
log "armed (pid $$): waiting for chain510's process to exit (no other build), builder == HEAD"
P510="$(cat experiments/510-october-recipe-rehearsal/work/chain510.pid 2>/dev/null)"
until ! ps -p "${P510:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain510 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-8])\.sh|chain496p\.sh|chain497e\.sh|chain510\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610030000 ] && { log "prerequisites not met by 00:00Z 10-03 -- stopping"; log "END"; exit 1; }; sleep 60; done
OLD="$(shasum -a 256 $BU | cut -d' ' -f1)"; cp $BU $X/work/builder-pre511.py
back(){ cp $X/work/builder-pre511.py $BU; log "RESTORED ${OLD:0:16}: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
log "released: builder ${OLD:0:16}"
uv run python $X/patch_builder511.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${OLD:0:16} -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
for t in experiments/479-pertarget-sharpen/test_pt479.py experiments/468-prior-gate/test_gate468.py experiments/472-gate-absolute/test_gate472.py experiments/451-sameline-mode/test_sameline451.py; do
  uv run python $t > $X/work/$(basename $t .py).log 2>&1 || back "$t failed"; done
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t.log 2>&1 || back "$t failed"; done
uv run python $X/test_ties511.py $X/work/builder-pre511.py $BU > $X/work/test_ties511.log 2>&1 || back "test_ties511 failed"
log "unit tests pass (479, 468, 472, 451, 417, ties511: $(tail -1 $X/work/test_ties511.log))"
# 2. KA-SP2: 479sp2's own configuration (the package chain479b built; bulk untouched by the fix on A/B/C)
R="$REPO/$(ls submissions/479sp2-*.h5ad.meta.json 2>/dev/null | head -1)"
[ -s "$R" ] || back "479sp2's meta is missing (chain479b did not build it)"
F2="$(uv run python -c 'import json,sys; s=json.load(open(sys.argv[1]))["config"]["bulk_ctx_scale"]; print(",".join(f"{k}={s[k]:.4f}" for k in "ABC"))' "$R")"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN
  export KA_Q=stercile:2.5,2.0,1.75 VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F2"
  bash $X/run_ka_511.sh kasp2 "$R" "$NEW" ) >> $X/work/kasp2.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_kasp2/runkasp2.log || back "KA-SP2 (479sp2's configuration) failed"
log "KA-SP2 PASS: 479sp2's configuration rebuilds 479sp2 byte for byte ($(/usr/bin/grep -o 'output sha256 [0-9a-f]* vs reference [0-9a-f]*' $X/work/build_kasp2/runkasp2.log))"
log "END builder ${NEW:0:16} (uncommitted; commit_after_chain.sh commits with the KA line)"
