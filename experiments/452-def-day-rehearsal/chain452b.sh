#!/usr/bin/env bash
# 452b: def_day.sh on ef-bundle-r with the same-line call FORCED (D=K562 300/300, E=HCT116 300/300, F=A549 6/300 rows only):
# heavy blending at w 2, :pt:noprog. Checks the loop converges when the blend moves every target of two contexts. After 452.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/452-def-day-rehearsal; LOG=$X/work/chain452b.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -qE " END$|stopping|never ended" $X/work/chain452.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain452\.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609260000 ] && { log "452 never ended by 00:00Z"; exit 1; }; sleep 30; done
grep -q "VCC_SAMELINE" scripts/make_submission_pm_cd4.py || { log "builder has no VCC_SAMELINE -- stopping"; exit 1; }
log "def_day forced rehearsal start"
DEF_DAY_FORCE="D=K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE,E=HCT116_LARGE_INTESTINE,F=A549_LUNG" \
  bash scripts/def_day.sh experiments/238-essential-free-panel/work/ef-bundle-r $X/work/efr_forced def452b 1.75 2 pt:noprog 4 >> $X/work/def_day_b.nohup 2>&1
log "def_day rc=$?: $(tail -3 $X/work/efr_forced/def_day.log | tr '\n' ' ' | cut -c1-400)"
rm -f $X/work/efr_forced/def452b_final.h5ad $X/work/efr_forced/def452b-lvl*.h5ad
log "END"
