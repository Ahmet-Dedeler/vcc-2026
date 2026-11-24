#!/usr/bin/env bash
# R7a access, v3 (R7 owner [d80a10], 10-16 ~10:35Z; coordinator [2d6f46] OK 10:33Z): resume after ADMISSION-ONLY DEVIATION 2
# (T4:comp551 through 552's rss_lift552.py; ledger R7 section, driver pin 423b8bfe). Starts r7a551.py access only when
#   (1) no bash process runs experiments/707-strong-tier-direction/build707.sh (10-17 s2's package; the coordinator's order:
#       T4 after build707 exits), matched as interpreter + path (AGENTS 6.2), and
#   (2) def_pipeline.other_builders() is empty.
# Then 551's own gate admits the step (22 GiB free+inactive, pressure 55 %, swap 400 MiB, the window). Resumable: only
# T4:comp551, T4's operator deletion and the aggregate remain. Its name matches gate602's r7a_access_* yield pattern.
# Start cap: the ledger's window extension minus 30 min. Terminal lines, always: "<ts> END rc=<rc>" appended to BOTH
# experiments/551-r6-c10b/work/r7a/access.log (the coordinator's GO and the lab's HOLD key on it) and access_wrapper_1015.log.
#   nohup bash scripts/detach.sh bash experiments/552-r7-reserve/work/r7a_access_1016c.sh > experiments/552-r7-reserve/work/r7a_access_1016c.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
WL=experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log; AL=experiments/551-r6-c10b/work/r7a/access.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$WL"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
r7cap(){ local w c
  w="$(awk '/^## Reserve R7/{s=1;next} s&&/^## /{s=0} s' docs/CONFIRMATION-LEDGER.md \
       | /usr/bin/grep -oE 'R7a WINDOW END EXTENDED TO [0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}Z' | tail -1 | awk '{print $NF}')"
  c="$(date -j -u -v-30M -f '%Y-%m-%dT%H:%MZ' "$w" +%Y%m%d%H%M 2>/dev/null)"; echo "${c:-202610161730}"; }
RC=9
fin(){ local ts; ts="$(date -u +%FT%TZ)"; printf '%s END rc=%s\n' "$ts" "$RC" >> "$AL"; printf '%s END rc=%s\n' "$ts" "$RC" >> "$WL"; }
trap fin EXIT
log "v3 armed (pid $$): waiting for build707.sh to exit and no def_pipeline build (start cap $(r7cap))"
until ! alive experiments/707-strong-tier-direction/build707.sh \
      && ! uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' 2>/dev/null; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$(r7cap)" ] && { log "start cap passed before the conditions held"; RC=9; exit 9; }
  sleep 60
done
[ -e experiments/691-release-watcher/work/RELEASE-DETECTED ] && { log "skipped: release detected"; RC=8; exit 8; }
log "conditions hold: resuming r7a551.py access (T4:comp551 through rss_lift552.py)"
printf '\n===== %s resume v3 (DEVIATION 2) =====\n' "$(date -u +%FT%TZ)" >> "$AL"
uv run python experiments/551-r6-c10b/r7a551.py access >> "$AL" 2>&1
RC=$?
log "access rc=$RC: $(tail -3 "$AL" | tr '\n' ' ' | cut -c1-400)"
