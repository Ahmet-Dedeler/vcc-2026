#!/usr/bin/env bash
# 417 window-2 chain (PI, 2026-09-22): after run417b_build.sh restores the builder, apply the three window-2 patches
# (415 LODO + 410 floor-0.5 keys, VCC_BULK_CTX_SCALE, VCC_BULK_DUMP/LOAD, VCC_CARRIER_GAMMA), run every unit test and the
# builder suite (must be 56 passed / 0 failed / 1 skipped), commit the builder by explicit path, pin 417c's builder hash,
# run the carrier known answer (a dump-exit run at VCC_CARRIER_GAMMA=0.65 must print max|carrier-bulk|=0.000e+00), then
# launch 417c's build. Any failure stops the chain and leaves the builder as it was (restored from the committed copy).
set -uo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"; cd "$REPO"
E="$REPO/experiments/417-bulk-sharpen"; W="$E/work"; LOG="$W/window2_chain.log"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
fail(){ log "STOP: $*"; git checkout -- scripts/make_submission_pm_cd4.py 2>/dev/null; log "builder now $(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; exit 1; }
log "WAITING for 417b's builder restore (pid $$)"
n=0; until grep -qE "builder after restore|stopping|NO ARTIFACT|GATE never" "$W/build_b/run417b.log" 2>/dev/null; do n=$((n+1)); [ $n -ge 720 ] && { log "417b never restored in 6 h"; exit 1; }; sleep 20; done
grep -q "builder after restore" "$W/build_b/run417b.log" || log "417b build stopped before restore: $(grep -E 'stopping|NO ARTIFACT|GATE never' "$W/build_b/run417b.log" | tail -1)"
B0="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
[ "$B0" = "400ef6e307896deb8c0adf694547c756000651f77ffb923ab450ef2f01389e31" ] || { log "builder $B0 is not the committed 400ef6e3 -- not patching"; exit 1; }
[ -z "$(git status --short scripts/make_submission_pm_cd4.py)" ] || { log "builder has uncommitted changes -- not patching"; exit 1; }
for p in patch_next417 patch_bulkio417 patch_carrier417; do uv run python "$E/$p.py" --write > "$W/$p.apply.log" 2>&1 || fail "$p did not apply (see $W/$p.apply.log)"; done
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; log "patched builder $B"
for t in test_sharpen417 test_next417 test_bulkio417 test_carrier417; do
  uv run python "$E/$t.py" > "$W/$t.log" 2>&1 || fail "$t FAILED (see $W/$t.log)"; log "$(tail -1 "$W/$t.log")"; done
for t in test_gain410d2 test_gain410d test_ballast410 test_pert_scale; do
  uv run python "experiments/410-pertscale-revisit/$t.py" > "$W/$t.w2.log" 2>&1 || fail "$t FAILED"; log "$(tail -1 "$W/$t.w2.log" | cut -c1-120)"; done
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > "$W/builder-tests-w2.log" 2>&1
grep -q "56 passed, 0 failed, 1 skipped" "$W/builder-tests-w2.log" || fail "builder suite: $(grep -E 'passed.*failed' "$W/builder-tests-w2.log" | tail -1)"
log "builder suite: $(grep -E 'passed.*failed' "$W/builder-tests-w2.log" | tail -1)"
cp "$W/builder-tests-w2.log" "$E/work/builder-tests-417w2.log"
git add scripts/make_submission_pm_cd4.py "$E/work/builder-tests-417w2.log" && git commit -q -m "Builder window 2 (417): 415 LODO + 410 floor-0.5 key artifacts, VCC_BULK_CTX_SCALE, VCC_BULK_DUMP/LOAD, VCC_CARRIER_GAMMA

Builder 400ef6e3 -> ${B:0:8}. All knobs unset is bit-identical (56 passed / 0 failed / 1 skipped). Unit tests PASS:
test_sharpen417, test_next417, test_bulkio417, test_carrier417, test_gain410d2, test_gain410d, test_ballast410,
test_pert_scale. Applied by experiments/417-bulk-sharpen/window2_chain.sh after run417b_build.sh restored the builder.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" || fail "commit failed"
log "committed builder $(git log -1 --format=%h)"
sed -i '' "s/__BUILDER_SHA__/$B/g" "$E/run417c_build.sh"
git add "$E/run417c_build.sh" && git commit -q -m "417c build pinned to builder ${B:0:8}

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
# the carrier known answer: 417a's env + VCC_CARRIER_GAMMA=0.65, dump-exit
export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR=clean4q80 VCC_PERT_SCALE=wshare_lo0.3_hi3 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25 VCC_CARRIER_GAMMA=0.65
unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_BULK_CTX_SCALE VCC_BULK_LOAD
bash "$E/run_dump417.sh" ka417a_g065 "$W/dump_ka417a_g065.npz" "$B" >> "$LOG" 2>&1
KA="$(grep -h '\[carrier-gamma\]' "$W/dump_ka417a_g065/logs/build.log" 2>/dev/null | tail -1)"
log "KA line: $KA"
echo "$KA" | grep -q "max|carrier-bulk|=0.000e+00" || { log "STOP: carrier known answer did not print max|carrier-bulk|=0 -- 417c NOT launched"; exit 1; }
[ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$B" ] || { log "STOP: builder not restored to $B after the dump"; exit 1; }
unset VCC_CARRIER_GAMMA VCC_BULK_DUMP
log "carrier known answer PASS; launching 417c"
nohup bash "$E/run417c_build.sh" > "$W/run417c.nohup" 2>&1 &
log "417c build pid $!"
log "END"
