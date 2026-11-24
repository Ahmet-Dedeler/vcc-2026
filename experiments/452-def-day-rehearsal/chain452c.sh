#!/usr/bin/env bash
# 452c: 452b with the same-line weight on `auto` (builder >= 4feae27d, 451b): the forced calls D=K562 (k562_gwps, strong ->
# w 2), E=HCT116 (orion_hct116, strong -> w 2), F=A549 (jiang_A549 + myllia_a549, weak -> w 0.5). The October default path
# (`auto`) has never run end to end through def_day.sh. Expected: D and E exactly as 452b (same w), F at w 0.5 with its
# 6 rows. Waits for 456's chain to END (builder serialisation).
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/452-def-day-rehearsal; LOG=$X/work/chain452c.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $X/work/chain452c.pid
log "armed (pid $$): waiting for 456's chain to END"
until /usr/bin/grep -qE " END$|stopping" experiments/456-zs-riders/work/chain456.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain456\.sh|run_arm_generic\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609262000 ] && { log "456 never ended by 20:00Z -- stopping"; exit 1; }; sleep 60; done
/usr/bin/grep -q "VCC_SAMELINE" scripts/make_submission_pm_cd4.py || { log "builder has no VCC_SAMELINE -- stopping"; exit 1; }
log "def_day forced rehearsal with auto w start"
DEF_DAY_FORCE="D=K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE,E=HCT116_LARGE_INTESTINE,F=A549_LUNG" \
  bash scripts/def_day.sh experiments/238-essential-free-panel/work/ef-bundle-r $X/work/efr_forced_auto def452c 1.75 auto pt:noprog 4 >> $X/work/def_day_c.nohup 2>&1
log "def_day rc=$?: $(tail -3 $X/work/efr_forced_auto/def_day.log | tr '\n' ' ' | cut -c1-400)"
uv run python - >> "$LOG" 2>&1 <<'PY'
import json
m = json.load(open("experiments/452-def-day-rehearsal/work/efr_forced_auto/def452c_final.h5ad.meta.json"))
sl = m["config"].get("sameline", {})
print("sameline w:", {c: v.get("w") for c, v in sl.get("contexts", {}).items()}, "| flag w:", sl.get("w"))
print("bulk_ctx_scale:", m["config"].get("bulk_ctx_scale"))
PY
rm -f $X/work/efr_forced_auto/def452c_final.h5ad $X/work/efr_forced_auto/def452c-lvl*.h5ad
log "END"
