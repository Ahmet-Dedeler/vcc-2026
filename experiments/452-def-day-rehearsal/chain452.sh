#!/usr/bin/env bash
# 452: the first end-to-end rehearsal of scripts/def_day.sh on ef-bundle-r with its REAL identity file path (fresh work dir:
# the whitener is re-fitted). Expected: D = Jurkat is matched (programme off by the default flags pt:noprog; 9 covered
# targets blended at w 2), E/F unchanged; the secant loop converges in <= 3 levelled builds. Serial, after 450d.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/452-def-day-rehearsal; LOG=$X/work/chain452.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
until grep -q "QUEUE DONE\|failed\|never ended" experiments/450-same-line-programme/work/drive450d.log 2>/dev/null && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /drive450d\.sh/' | grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609252300 ] && { log "450d never ended by 23:00Z"; exit 1; }; sleep 30; done
grep -q "VCC_SAMELINE" scripts/make_submission_pm_cd4.py || { log "builder has no VCC_SAMELINE (451 patch not applied) -- stopping"; exit 1; }
log "def_day rehearsal start"
bash scripts/def_day.sh experiments/238-essential-free-panel/work/ef-bundle-r $X/work/efr def452 1.75 2 pt:noprog 4 >> $X/work/def_day.nohup 2>&1
log "def_day rc=$?: $(tail -3 $X/work/efr/def_day.log | tr '\n' ' ' | cut -c1-400)"
rm -f $X/work/efr/def452_final.h5ad $X/work/efr/def452-lvl*.h5ad
log "END"
