#!/usr/bin/env bash
# 471b: 471's dump (h1_2025_full in the H1 union) with 446b's re-level factors divided by the dump's global scale 0.988146
# (energy-preserving key; see PREREGISTRATION addendum). Serial after chain472. Then ntclevel417 on the artifact (rho must
# match 446b's level within 0.5 %). Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/471-h1full-union; LOG=$X/work/chain471b.log; E=experiments/428-zero-gene-audit
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain471b.pid
log "armed (pid $$): waiting for chain472 to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ END( builder .*)?$" experiments/472-gate-absolute/work/chain472.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain472\.sh|chain470b\.sh|chain471\.sh|run_arm_generic\.sh|run_ka_47[12]\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271400 ] && { log "472 never ended by 14:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
until [ "$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)" = "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271400 ] && { log "builder never equal to HEAD by 14:00Z 09-27 -- stopping"; exit 1; }
  log "builder != HEAD (472 uncommitted?) -- waiting"; sleep 300; done
B="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
EM=edb5b8dcbee519683b9d50c3556d4470ca6290594d4d1611be3d55c4fb7d41c7
N2="$REPO/$X/work/dump_h1fx2_cd4x15.npz"; DS2=$(shasum -a 256 "$N2" | cut -d' ' -f1)
[ "${DS2:0:16}" = 30086b5bfc6a5d07 ] || { log "dump sha ${DS2:0:16} is not chain471's 30086b5bfc6a5d07 -- stopping"; exit 1; }
log "dump ${DS2:0:16} builder ${B:0:16}"
R446B="$REPO/submissions/446b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackload-ctxlevel-mde010coarse-zr5.h5ad.meta.json"
( unset VCC_SAMELINE VCC_DE_PRIOR_GATE VCC_UNION_SOURCES
  export VCC_BULK_LOAD="$N2:$DS2" VCC_DE_MDE=0.1,0.1,0.3,coarse VCC_DE_ZERO_REPAIR=-5 VCC_BULK_CTX_SCALE="A=1.1136,B=1.165516,C=1.186059"
  bash $E/run_arm_generic.sh 471b 1.75 "$R446B" the-ph-meter-joined-a-union \
    471b-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stackh1fload-ctxlevelr-mde010coarse-zr5 "$B" "$EM" bulk_load,bulk_load_sha256,bulk_ctx_scale ) >> $X/work/build471b.nohup 2>&1
log "471b: $(/usr/bin/grep -E '^\[471b\] (PASS|FAIL|CLIFF)|audit|package sha|novelty|DONE$|stopping' $E/work/build_471b/run471b.log | tr '\n' ' ' | cut -c1-900)"
H=$E/work/build_471b/the-ph-meter-joined-a-union.h5ad
if [ -s "$H" ]; then
  uv run python experiments/417-bulk-sharpen/ntclevel417.py "$H" --bundle data/challenge --contexts A,B,C --sd 0.01 >> $LOG 2>&1
  log "471b rho: $(uv run python scripts/def_day_helpers.py rho "${H%.h5ad}.ntclevel.json" 2>&1 | tail -1)"
fi
rm -f "$H"
log "END"
