#!/usr/bin/env bash
# 471: the denoised H1 bank (h1_2025_full) in place of h1_2025 in the stack's H1 union. Serial after chain470.
#  1. bctrl471.py --write: data/donors/h1_2025_full/b_ctrl.npy (KA: the recomputed delta equals 458d's bit for bit);
#  2. the builder registers h1_2025_full (patch_builder471.py, a pure addition); unit tests (suite 56/0, 468, 451);
#  3. KA-1: every knob unset rebuilds 446b byte for byte;
#  4. dump KA: the default union re-dumps 446b's stack (dump_h1x2_cd4x15.npz) array for array;
#  5. the new dump with VCC_UNION_SOURCES=gse295214_h1,h1_2025_full;
#  6. arm 471 = 446b + VCC_BULK_LOAD=<new dump> via run_arm_generic (config diff vs 446b exactly bulk_load,
#     bulk_load_sha256; 446b's re-level factors kept; the cliff check must pass).
# Needs HEAD's builder = 693fde15 (468 committed) at the start. Any failure restores that builder and stops. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/471-h1full-union; LOG=$X/work/chain471.log; BU=scripts/make_submission_pm_cd4.py; E=experiments/428-zero-gene-audit
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain471.pid
B0=693fde151be7c9994b3aa643480d8dbf759242b2f5a24ed044e688e7e8f0a405
EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
log "armed (pid $$): waiting for chain470 to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|.*stopping)$" experiments/470-weak-target-prior/work/chain470.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain470\.sh|chain469\.sh|chain468\.sh|def_day\.sh|run_arm_generic\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271200 ] && { log "470 never ended by 12:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
BF="$(shasum -a 256 $BU | cut -d' ' -f1)"; BH="$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)"
[ "$BF" = "$B0" ] && [ "$BH" = "$B0" ] || { log "builder ${BF:0:16} / HEAD ${BH:0:16} is not the committed 693fde15 -- stopping"; exit 1; }
cp $BU $X/work/builder-693fde15.py
back(){ cp $X/work/builder-693fde15.py $BU; log "RESTORED 693fde15: $(shasum -a 256 $BU | cut -c1-16) ($1)"; log "END"; exit 1; }
# 1
if [ ! -s data/donors/h1_2025_full/b_ctrl.npy ]; then
  /usr/bin/time -l uv run python $X/bctrl471.py --write > $X/work/bctrl471.log 2>&1 || { log "bctrl471 failed: $(tail -2 $X/work/bctrl471.log | tr '\n' ' ')"; exit 1; }
fi
log "b_ctrl: $(/usr/bin/grep -E 'KA-1|KA-2|wrote' $X/work/bctrl471.log | tr '\n' ' ' | cut -c1-400)"
# 2
uv run python $X/patch_builder471.py $BU >> $LOG 2>&1 || back "patch failed"
NEW="$(shasum -a 256 $BU | cut -d' ' -f1)"; log "patched: ${B0:0:16} -> ${NEW:0:16}"
uv run python experiments/410-pertscale-revisit/run_builder_tests.py > $X/work/builder-tests-after.log 2>&1
log "builder suite: $(tail -1 $X/work/builder-tests-after.log)"
/usr/bin/grep -q "56 passed, 0 failed" $X/work/builder-tests-after.log || back "builder suite not 56/0"
uv run python experiments/468-prior-gate/test_gate468.py > $X/work/test_gate468.log 2>&1 || back "468 tests failed"
uv run python experiments/451-sameline-mode/test_sameline451.py > $X/work/test_sameline451.log 2>&1 || back "451 tests failed"
log "468 and 451 unit tests pass"
D="$REPO/experiments/417-bulk-sharpen/work/dump_h1x2_cd4x15.npz"; DS=$(shasum -a 256 "$D" | cut -d' ' -f1)
R446B="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
# 3
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES VCC_DE_MAG VCC_RC_AUDIT_DUMP
  export VCC_BULK_LOAD="$D:$DS" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $X/run_ka_471.sh ka1 "$R446B" "$NEW" ) >> $X/work/ka1.nohup 2>&1
/usr/bin/grep -q "KA PASS" $X/work/build_ka1/runka1.log || back "KA-1 (unset == 446b) failed"
log "KA-1 PASS: unset rebuilds 446b byte for byte"
# 4, 5: the dumps (417's dump env for dump_h1x2_cd4x15; q is a render knob and does not enter the dump)
dump(){ ( unset VCC_NORMREL_PM VCC_RC_ORDER VCC_CTX_A_PRIOR VCC_CTX_A_PRIOR_BETA VCC_SIGN_REPAIR VCC_BULK_GENE_GAIN VCC_BULK_CTX_SCALE \
                VCC_BULK_LOAD VCC_CARRIER_GAMMA VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_DE_MDE VCC_DE_ZERO_REPAIR VCC_UNION_SOURCES
          export VCC_NORMREL_EXP=-1.5 VCC_CD4_PARTIAL_G=0.5 VCC_NORMREL_MODE=cap VCC_DE_SIGN_PRIOR=clean4q80 VCC_PERT_SCALE=wshare_lo0.3_hi3 \
                 VCC_PERT_SCALE_BALLAST=1 VCC_BULK_SHARPEN_Q=1.25 VCC_H1_UNION_G=0.56 VCC_CD4_TARGET_F=0.51
          [ -n "$3" ] && export VCC_UNION_SOURCES="$3"
          bash $X/run_dump471.sh "$1" "$2" "$NEW" ) >> $X/work/dumps.nohup 2>&1; }
dump ka_h1x2_cd4x15 "$REPO/$X/work/dump_ka_h1x2_cd4x15.npz" ""
uv run python - "$REPO/$X/work/dump_ka_h1x2_cd4x15.npz" "$D" <<'PY' >> $LOG 2>&1 || back "dump KA failed"
import sys, numpy as np
a = np.load(sys.argv[1], allow_pickle=True); b = np.load(sys.argv[2], allow_pickle=True)
arrs = [k for k in b.files if k not in ("env", "meta", "provenance")]
same = {k: bool(np.array_equal(a[k], b[k])) for k in arrs}
print("dump KA (default union == dump_h1x2_cd4x15):", same, "PASS" if all(same.values()) else "FAIL")
sys.exit(0 if all(same.values()) else 1)
PY
log "dump KA PASS: the default union re-dumps 446b's stack array for array"
N2="$REPO/$X/work/dump_h1fx2_cd4x15.npz"
dump h1fx2_cd4x15 "$N2" "gse295214_h1,h1_2025_full"
[ -s "$N2" ] || back "the h1_2025_full dump failed"
uv run python - "$N2" "$D" <<'PY' >> $LOG 2>&1
import sys, numpy as np
a = np.load(sys.argv[1], allow_pickle=True); b = np.load(sys.argv[2], allow_pickle=True)
print("new dump vs 446b's:", {k: (a[k].shape, bool(np.array_equal(a[k], b[k]))) for k in b.files if k not in ("env", "meta", "provenance")})
if "bulk" in b.files:
    d = np.abs(np.asarray(a["bulk"], float) - np.asarray(b["bulk"], float)); rows = (d.max(1) > 0) if d.ndim == 2 else None
    print("bulk rows changed:", None if rows is None else int(rows.sum()), "max |d|", float(d.max()))
PY
DS2=$(shasum -a 256 "$N2" | cut -d' ' -f1); log "new dump ${DS2:0:16}: $(tail -2 $LOG | tr '\n' ' ' | cut -c1-300)"
# 6
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES
  export VCC_BULK_LOAD="$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1004,B=1.1517,C=1.1720"
  bash $E/run_arm_generic.sh 471 1.75 "$R446B" the-hemocytometer-lost-count \
    471-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackh1fload-ctxlevel-mde010coarse-zr5 "$NEW" "$EM" bulk_load,bulk_load_sha256 ) >> $X/work/build471.nohup 2>&1
log "471: $(/usr/bin/grep -E '^\[471\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping' $E/work/build_471/run471.log | tr '\n' ' ' | cut -c1-900)"
rm -f $E/work/build_471/the-hemocytometer-lost-count.h5ad
log "END builder ${NEW:0:16} (uncommitted; the coordinator commits with KA-1 and the dump KA)"
