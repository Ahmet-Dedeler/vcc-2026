#!/usr/bin/env bash
# 495: the LODO key at cap 6 on 446bj2 (PREREGISTRATION.md). Serial after chain479b (process exit + its END line).
#  1. register lodo415_lo0.5_hi6 (patch_builder495.py, a pure addition); unit tests (suite 56/0, 479, 468, 472, 451, 417);
#  2. KA-1: 446bj2's env (LODO dump, q 1.75, its factors) rebuilds 446bj2's output 711ee43b byte for byte;
#  3. dump KA: 446bj2's dump env with lodo415_lo0.5_hi3 re-dumps dump_h1x2_cd4x15_lodo05.npz array for array;
#  4. the cap-6 dump (lodo415_lo0.5_hi6);
#  5. arm r1 = 446bj2 + the cap-6 dump at 446bj2's factors (diff exactly bulk_load, bulk_load_sha256) -> ntclevel ->
#     r2 at factors x EV-scale (diff + bulk_ctx_scale) -> ntclevel (rho within 0.5 % of 0.9793).
# Any failure in 1-4 restores the pre-495 builder. The builder stays patched (uncommitted) after success: the coordinator
# commits it with KA-1 and the dump KA. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/495-lodo-cap6; LOG=$X/work/chain495.log; BU=scripts/make_submission_pm_cd4.py; E=experiments/428-zero-gene-audit
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain495.pid
EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
log "armed (pid $$): waiting for chain479b's process to exit (and no other build), then builder == HEAD"
P479B="$(cat experiments/479-pertarget-sharpen/work/chain479b.pid 2>/dev/null)"
until ! ps -p "${P479B:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain479b \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain479b?\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609280600 ] && { log "479b never ended by 06:00Z 09-28 -- stopping"; log "END"; exit 1; }; sleep 60; done
log "chain479b exited; its last lines: $(tail -3 experiments/479-pertarget-sharpen/work/chain479b.log | tr '\n' ' ' | cut -c1-400)"
until [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609280600 ] && { log "builder never equal to HEAD by 06:00Z 09-28 -- stopping"; log "END"; exit 1; }
  log "builder != HEAD (479 uncommitted?) -- waiting"; sleep 300; done
OLD="$(shasum -a 256 $BU | cut -d' ' -f1)"; cp $BU $X/work/builder-pre495.py
back(){ cp $X/work/builder-pre495.py $BU; log "RESTORED ${OLD:0:16}: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
# 1
uv run python $X/patch_builder495.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${OLD:0:16} -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
for t in experiments/479-pertarget-sharpen/test_pt479.py experiments/468-prior-gate/test_gate468.py experiments/472-gate-absolute/test_gate472.py experiments/451-sameline-mode/test_sameline451.py; do
  uv run python $t > $X/work/$(basename $t .py).log 2>&1 || back "$t failed"; done
for t in test_next417 test_sharpen417 test_bulkio417 test_carrier417; do
  uv run python experiments/417-bulk-sharpen/$t.py > $X/work/$t.log 2>&1 || back "$t failed"; done
log "unit tests pass (479, 468, 472, 451, 417)"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15_lodo05.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
[ "${DS:0:16}" = 00ccfcba1b67fa24 ] || back "LODO dump sha ${DS:0:16} unexpected"
R="$REPO/$(ls submissions/446bj2-*.h5ad.meta.json | head -1)"; F0="A=1.0962,B=1.1461,C=1.1681"
# 2
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export KA_Q=1.75 VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$F0"
  bash $X/run_ka_495.sh ka1 "$R" "$NEW" ) >> $X/work/ka1.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_ka1/runka1.log || back "KA-1 (446bj2's env == 446bj2) failed"
log "KA-1 PASS: 446bj2's env rebuilds 446bj2 byte for byte ($(/usr/bin/grep -o 'output sha256 [0-9a-f]* vs reference [0-9a-f]*' $X/work/build_ka1/runka1.log))"
# 3, 4: the dumps (window3's env for dump_h1x2_cd4x15_lodo05; q is a render knob and does not enter the dump)
dump(){ ( unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_BULK_CTX_SCALE \
                VCC_BULK_LOAD VCC_CARRIER_GAMMA VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_DE_MDE VCC_DE_ZERO_REPAIR VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP
          export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR=clean4q80 VCC_PERT_SCALE="$3" \
                 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25 VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51
          bash $X/run_dump495.sh "$1" "$2" "$NEW" ) >> $X/work/dumps.nohup 2>&1; }
dump ka_lodo05 "$REPO/$X/work/dump_ka_lodo05.npz" lodo415_lo0.5_hi3
uv run python - "$REPO/$X/work/dump_ka_lodo05.npz" "$D" <<'PY' >> $LOG 2>&1 || back "dump KA failed"
import sys, numpy as np
a = np.load(sys.argv[1], allow_pickle=True); b = np.load(sys.argv[2], allow_pickle=True)
arrs = [k for k in b.files if k not in ("env", "meta", "provenance")]
same = {k: bool(np.array_equal(a[k], b[k])) for k in arrs}
print("dump KA (lodo415_lo0.5_hi3 == dump_h1x2_cd4x15_lodo05):", same, "PASS" if all(same.values()) else "FAIL")
sys.exit(0 if all(same.values()) else 1)
PY
log "dump KA PASS: lodo415_lo0.5_hi3 re-dumps 446bj2's stack array for array"
rm -f "$REPO/$X/work/dump_ka_lodo05.npz"
N2="$REPO/$X/work/dump_h1x2_cd4x15_lodo05hi6.npz"
dump lodo05hi6 "$N2" lodo415_lo0.5_hi6
[ -s "$N2" ] || back "the cap-6 dump failed"
uv run python - "$N2" "$D" <<'PY' >> $LOG 2>&1
import sys, numpy as np
a = np.load(sys.argv[1], allow_pickle=True); b = np.load(sys.argv[2], allow_pickle=True)
print("cap-6 dump vs 446bj2's:", {k: (a[k].shape, bool(np.array_equal(a[k], b[k]))) for k in b.files if k not in ("env", "meta", "provenance")})
for k in ("bulk", "effects"):
    if k in b.files:
        x, y = np.asarray(a[k], float), np.asarray(b[k], float)
        na, nb = np.linalg.norm(x, axis=1), np.linalg.norm(y, axis=1); r = na / np.where(nb > 0, nb, 1)
        print(f"{k}: rows changed {int((np.abs(x - y).max(1) > 0).sum())}; row-norm ratio min {r.min():.4f} median {np.median(r):.4f} max {r.max():.4f}; total energy ratio {float((na**2).sum() / (nb**2).sum()):.4f}")
PY
DS2=$(shasum -a 256 "$N2" | cut -d' ' -f1); log "cap-6 dump ${DS2:0:16}: $(tail -2 $LOG | tr '\n' ' ' | cut -c1-400)"
# 5
build(){ local tag=$1 fact=$2 name=$3 stem=$4 expect=$5
  ( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP
    export VCC_BULK_LOAD="$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="$fact"
    bash $E/run_arm_generic.sh "$tag" 1.75 "$R" "$name" "$stem" "$NEW" "$EM" "$expect" ) >> $X/work/build_$tag.nohup 2>&1
  log "$tag: $(/usr/bin/grep -E "^\[$tag\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping" $E/work/build_$tag/run$tag.log | tr '\n' ' ' | cut -c1-500)"
  /usr/bin/grep -q "^\[$tag\] PASS" $E/work/build_$tag/run$tag.log; }
level(){ uv run python experiments/417-bulk-sharpen/ntclevel417.py "$1" --bundle data/challenge --contexts A,B,C --sd 0.01 >> $LOG 2>&1; }
NM=the-water-bath-filed-for-divorce; S1=495r1-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevel-mde010coarse-zr5
build 495r1 "$F0" "$NM-r1" "$S1" bulk_load,bulk_load_sha256 || { log "495r1 build failed"; log "END builder ${NEW:0:16} (uncommitted; KA-1 + dump KA passed)"; exit 1; }
H1="$E/work/build_495r1/$NM-r1.h5ad"; level "$H1" || { log "495r1 ntclevel failed"; log "END"; exit 1; }
EV="$(/usr/bin/grep -o 'VCC_BULK_CTX_SCALE (relative to THIS artifact) = [A-Z0-9=.,]*' $LOG | tail -1 | awk -F'= ' '{print $2}')"
F1="$(uv run python -c "import sys; f0={k:float(v) for k,v in (p.split('=') for p in sys.argv[1].split(','))}; ev={k:float(v) for k,v in (p.split('=') for p in sys.argv[2].split(','))}; print(','.join(f'{k}={f0[k]*ev[k]:.4f}' for k in 'ABC'))" "$F0" "$EV")"
log "495r1 rho: $(uv run python scripts/def_day_helpers.py rho "${H1%.h5ad}.ntclevel.json" 2>&1 | tail -1); EV $EV -> round-2 factors $F1"
rm -f "$H1" "$REPO/submissions/$S1.vcc"
S2=495-cap15-partial50-clean4q80-lodo6-ballast-sharpen175-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5
build 495 "$F1" "$NM" "$S2" bulk_ctx_scale,bulk_load,bulk_load_sha256 || { log "495 build failed"; log "END"; exit 1; }
H2="$E/work/build_495/$NM.h5ad"; level "$H2"
log "495 rho: $(uv run python scripts/def_day_helpers.py rho "${H2%.h5ad}.ntclevel.json" 2>&1 | tail -1) (must be within 0.5 % of 0.9793)"
rm -f "$H2"
log "END builder ${NEW:0:16} (uncommitted; the coordinator commits with KA-1 and the dump KA)"
