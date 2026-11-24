#!/usr/bin/env bash
# 468: the builder already carries the VCC_DE_PRIOR_GATE patch (4feae27d -> 693fde15, applied 2026-09-26 16:50Z by
# patch_builder468.py; byte copy of 4feae27d in work/builder-4feae27d.py). This chain runs the unit tests now, waits for
# 467b to END (heavy jobs strictly serial), then two known-answer builds:
#   KA-1: every new knob unset                          == 446b's output sha (the production default is untouched)
#   KA-2: VCC_SAMELINE (Arc's H1 programme, w 0), no gate == 458d's output sha (the programme path is untouched)
# and finally the arm 458g = 458d + VCC_DE_PRIOR_GATE=q33 through run_arm_generic (config diff vs 458d exactly
# de_prior_gate). A test or KA failure restores 4feae27d and stops. Never submits. The coordinator commits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/468-prior-gate; LOG=$X/work/chain468.log; BU=scripts/make_submission_pm_cd4.py; E=experiments/428-zero-gene-audit
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain468.pid
OLD=4feae27d37ffa43687a2eaa96d4ca835d85900fbaac63504dc1f5987fe6cd2e9
NEW=693fde151be7c9994b3aa643480d8dbf759242b2f5a24ed044e688e7e8f0a405
EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
log "armed (pid $$)"
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$NEW" ] || { log "builder is not the 468 patch ${NEW:0:16} -- stopping"; exit 1; }
[ "$(shasum -a 256 $X/work/builder-4feae27d.py | cut -d' ' -f1)" = "$OLD" ] || { log "the byte copy of 4feae27d is missing -- stopping"; exit 1; }
back(){ cp $X/work/builder-4feae27d.py $BU; log "RESTORED 4feae27d: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
uv run python $X/test_gate468.py > $X/work/test_gate468.log 2>&1 || back "468 unit tests failed"
log "468 unit tests: $(tail -1 $X/work/test_gate468.log)"
uv run python experiments/451-sameline-mode/test_sameline451.py > $X/work/test_sameline451.log 2>&1 || back "451 unit tests failed"
log "451 unit tests: $(tail -1 $X/work/test_sameline451.log)"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t.log 2>&1 || back "$t failed"; done
log "417 unit tests pass; waiting for 467b to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ END$" experiments/467-strength-gated-programme/work/chain467b.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain467b\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609262100 ] && { log "467b never ended by 21:00Z -- stopping (builder stays patched; tests passed)"; exit 1; }
  sleep 30; done
log "467b ended; KA-1 starts"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
R446B="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
R458D="$REPO/submissions/458d-cap15-partial50-arch1dprog-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
A=experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz; AS=$(shasum -a 256 $A | cut -d' ' -f1)
[ "$AS" = "6951dbc7498e1f2aeabf405d5cbe94e7b695a7d3c0848a3a9460daa52d87a5ba" ] || back "programme artifact sha $AS unexpected"
[ -s "$R446B" ] && [ -s "$R458D" ] || back "a reference meta is missing"
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $X/run_ka_468.sh ka1 "$R446B" "$NEW" ) >> $X/work/ka1.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_ka1/runka1.log || back "KA-1 (unset == 446b) failed"
log "KA-1 PASS: unset rebuilds 446b byte for byte ($(/usr/bin/grep -o 'output sha256 [0-9a-f]* vs reference [0-9a-f]*' $X/work/build_ka1/runka1.log))"
( unset VCC_DE_PRIOR_GATE VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" \
         VCC_SAMELINE="$REPO/$A:$AS:0"
  bash $X/run_ka_468.sh ka2 "$R458D" "$NEW" ) >> $X/work/ka2.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_ka2/runka2.log || back "KA-2 (programme, no gate == 458d) failed"
log "KA-2 PASS: VCC_SAMELINE without the gate rebuilds 458d byte for byte ($(/usr/bin/grep -o 'output sha256 [0-9a-f]* vs reference [0-9a-f]*' $X/work/build_ka2/runka2.log))"
log "458g build starts"
( export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720" \
         VCC_SAMELINE="$REPO/$A:$AS:0" VCC_DE_PRIOR_GATE=q33
  bash $E/run_arm_generic.sh 458g 1.75 "$R458D" the-petri-dish-joined-a-book-club \
    458g-cap15-partial50-arch1dprog-gate33-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5 "$NEW" "$EM" de_prior_gate ) >> $X/work/build458g.nohup 2>&1
log "458g: $(/usr/bin/grep -E '^\[458g\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping' $E/work/build_458g/run458g.log | tr '\n' ' ' | cut -c1-900)"
log "458g prior gate: $(/usr/bin/grep -hE '^\[prior_gate\]' $E/work/build_458g/driver.log $E/work/build_458g/logs/*.log 2>/dev/null | sort -u | tr '\n' ' ' | cut -c1-600)"
rm -f $E/work/build_458g/the-petri-dish-joined-a-book-club.h5ad
log "END builder ${NEW:0:16} (uncommitted until the coordinator reads KA-1/KA-2)"
