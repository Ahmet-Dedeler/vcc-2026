#!/usr/bin/env bash
# 496p: the uncalled-gene sign hook into production, then its DE-only arm on 446bj2. Serial after chain498 (process exit).
#  1. patch_emitter496.py on src/vcc2026/rc_emit_k.py (must give 1658ae2c400f75ee); patch_builder496.py (emitter pin +
#     VCC_DE_MDE_SIGN); builder suite 56/0 and the 479/468/472/451/417 tests;
#  2. KA-H1: the patched PRODUCTION emitter without a sign source on 435's H1 harness (seed 0) == 435's production run,
#     stats bit for bit (run496.py with RUN426_EMITTER=prod);
#  3. KA-1: 446bj2's env rebuilds 446bj2's output 711ee43b byte for byte;
#  4. the arm 496a = 446bj2 + VCC_DE_MDE_SIGN=clean4:0.05 (bulk untouched: no re-level). Config diff vs 446bj2 exactly
#     de_mde_sign + de_sign_emitter_sha256.
# Any failure in 1-3 restores both pre-496 files. The files stay patched after success; commit_after_chain.sh commits them.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/496-uncalled-sign-hook; LOG=$X/work/chain496p.log; BU=scripts/make_submission_pm_cd4.py; EMP=src/vcc2026/rc_emit_k.py
E=experiments/428-zero-gene-audit
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain496p.pid
EMNEW=1658ae2c400f75ee1659a8ffcaebca569bf3c8afd55ffd1c3006414344ef2791
log "armed (pid $$): waiting for chain498's process to exit (and no other build), builder and emitter == HEAD"
P498="$(cat experiments/498-combined-bulk/work/chain498.pid 2>/dev/null)"
until ! ps -p "${P498:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain498 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain4(7[69]b?|9[5-8])\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] \
      && [ "$(shasum -a 256 $EMP | cut -d' ' -f1)" = "$(git show HEAD:$EMP | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610020200 ] && { log "498 never ended by 02:00Z 10-02 -- stopping"; log "END"; exit 1; }; sleep 60; done
log "chain498 exited; its last lines: $(tail -2 experiments/498-combined-bulk/work/chain498.log | tr '\n' ' ' | cut -c1-300)"
cp $BU $X/work/pre496.builder.py; cp $EMP $X/work/pre496.rc_emit_k.py
back(){ cp $X/work/pre496.builder.py $BU; cp $X/work/pre496.rc_emit_k.py $EMP; log "RESTORED builder $(shasum -a 256 $BU | cut -c1-16) and emitter $(shasum -a 256 $EMP | cut -c1-16) ($1)"; log "END"; exit 1; }
[ "$(shasum -a 256 $EMP | cut -c1-16)" = edb5b8dcbee51968 ] || { log "emitter is not edb5b8dc -- stopping"; log "END"; exit 1; }
# 1
uv run python $X/patch_emitter496.py $EMP >> $LOG 2>&1 || back "emitter patch failed"
[ "$(shasum -a 256 $EMP | cut -d' ' -f1)" = "$EMNEW" ] || back "patched emitter is not ${EMNEW:0:16}"
uv run python $X/patch_builder496.py $BU $EMP >> $LOG 2>&1 || back "builder patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: emitter edb5b8dc -> ${EMNEW:0:16}; builder -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
for t in experiments/479-pertarget-sharpen/test_pt479.py experiments/468-prior-gate/test_gate468.py experiments/472-gate-absolute/test_gate472.py experiments/451-sameline-mode/test_sameline451.py; do
  uv run python $t > $X/work/$(basename $t .py).log 2>&1 || back "$t failed"; done
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t.log 2>&1 || back "$t failed"; done
log "unit tests pass (479, 468, 472, 451, 417)"
# 2
rm -f $X/work/h1/metrics-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.json $X/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz
RUN426_EMITTER=prod /usr/bin/time -l uv run python $X/run496.py base 0 --bulk ps --mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse > $X/work/kah1.log 2>&1 || back "KA-H1 run failed"
uv run python - <<'PY' >> $LOG 2>&1 || back "KA-H1 failed"
import numpy as np, sys
a = np.load("experiments/496-uncalled-sign-hook/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
b = np.load("experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz")
ok = set(a.files) == set(b.files) and all(np.array_equal(a[k], b[k]) for k in b.files)
print("KA-H1 PASS: the patched production emitter without a sign source == 435's production run (stats bit-equal)" if ok else "KA-H1 FAIL")
sys.exit(0 if ok else 1)
PY
# 3
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
R="$REPO/$(ls submissions/446bj2-*.h5ad.meta.json | head -1)"; F0="A=1.0962,B=1.1461,C=1.1681"
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP VCC_DE_MDE_SIGN
  export KA_Q=1.75 VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0"
  bash $X/run_ka_496.sh ka1 "$R" "$NEW" ) >> $X/work/ka1.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_ka1/runka1.log || back "KA-1 (446bj2's env == 446bj2) failed"
log "KA-1 PASS: 446bj2's env rebuilds 446bj2 byte for byte ($(/usr/bin/grep -o 'output sha256 [0-9a-f]* vs reference [0-9a-f]*' $X/work/build_ka1/runka1.log))"
# 4
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0" VCC_DE_MDE_SIGN=clean4:0.05
  bash $E/run_arm_generic.sh 496a 1.75 "$R" hook496a \
    496a-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-as05 "$NEW" "$EMNEW" de_mde_sign,de_sign_emitter_sha256 ) >> $X/work/build496a.nohup 2>&1
log "496a: $(/usr/bin/grep -E '^\[496a\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping' $E/work/build_496a/run496a.log | tr '\n' ' ' | cut -c1-500)"
log "496a de_mde_sign diagnostics: $(/usr/bin/grep -h '\[de_mde_sign\]' $E/work/build_496a/logs/build.log 2>/dev/null | tr '\n' ' ' | cut -c1-300)"
rm -f $E/work/build_496a/hook496a.h5ad
log "END builder ${NEW:0:16} emitter ${EMNEW:0:16} (uncommitted; commit_after_chain.sh commits with KA-H1 and KA-1)"
