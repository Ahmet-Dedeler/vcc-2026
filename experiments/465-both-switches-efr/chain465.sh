#!/usr/bin/env bash
# 465: the October day with BOTH new switches on (DEF_DAY_LODO=1, DEF_DAY_PROGRAMME_SOURCE=h1_2025_full, flags pt) on a
# panel that is NOT A/B/C's: 238's essential-free realistic draw `ef-bundle-r` (efr; A/B/C controls relabelled D/E/F,
# real identity). 460 and 461 rehearsed each switch on the A/B/C panel; this exercises their panel-dependent code (the
# LODO fit on a new panel's donors/CD4 banks/operator; the programme source with a new panel excluded) end to end.
# A rehearsal only: never prepped for real, never submitted. Serial after 462l (drive462). 452's efr operator
# (75773d59) is cloned in so the whitener stage verifies and skips the fit.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/465-both-switches-efr; W=$X/work; LOG=$W/chain465.log; R452=experiments/452-def-day-rehearsal/work/efr
BUNDLE=experiments/238-essential-free-panel/work/ef-bundle-r; N=465efr
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/chain465.pid"
log "armed (pid $$): waiting for 462l's runs (drive462) to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|KA run failed|KA FAIL -- stopping|.*stopping)$" experiments/462-broad16-lab-signature/work/drive462.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /drive462\.sh|chain458e\.sh|chain461\.sh|run_arm_generic\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609272000 ] && { log "462l never ended by 20:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
[ "$BF" = "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)" ] || { log "builder != HEAD -- stopping"; exit 1; }
[ "$(shasum -a 256 scripts/def_day.sh | cut -d' ' -f1)" = "$(git show HEAD:scripts/def_day.sh | shasum -a 256 | cut -d' ' -f1)" ] || { log "def_day.sh != HEAD -- stopping"; exit 1; }
[ "$(shasum -a 256 "$R452/whitener_union_pruned_def.npz" | cut -c1-16)" = 75773d59131e0da1 ] || { log "452's efr operator is not 75773d59 -- stopping"; exit 1; }
mkdir -p "$W/run"; cp -c "$R452/whitener_union_pruned_def.npz" "$R452/whitener_union_pruned_def.json" "$W/run/"
log "starting DEF_DAY_LODO=1 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full def_day.sh on $BUNDLE (flags pt)"
DEF_DAY_LODO=1 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full bash scripts/def_day.sh "$BUNDLE" "$W/run" "$N" 1.75 auto pt 4 > "$W/def_day.nohup" 2>&1
log "def_day: $(/usr/bin/grep -E "LODO key ON|programme-source|same-line mode|CONVERGED at|NOT converged|final:|DONE|STOP" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-900)"
log "END"
