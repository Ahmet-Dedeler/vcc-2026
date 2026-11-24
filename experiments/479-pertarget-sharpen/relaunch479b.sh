#!/usr/bin/env bash
# Starts chain479b only once chain495 is past its own start gate. chain495 waits for "no bash process named chain479b.sh"
# (and nothing else of the 479 family), so a chain479b waiting beside it would deadlock both. This file's name is chosen
# so it matches none of the chains' gate regexes. Exits after launching.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
L=experiments/479-pertarget-sharpen/work/chain479b.log
P495="$(cat experiments/495-lodo-cap6/work/chain495.pid 2>/dev/null)"
until /usr/bin/grep -q "patched:" experiments/495-lodo-cap6/work/chain495.log 2>/dev/null \
      || ! ps -p "${P495:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain495; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609281800 ] && { printf '%s %s\n' "$(date -u +%FT%TZ)" "relauncher: chain495 never started by 18:00Z 09-28" >> "$L"; exit 1; }
  sleep 30; done
printf '%s %s\n' "$(date -u +%FT%TZ)" "relauncher: chain495 is past its gate; starting chain479b" >> "$L"
nohup bash experiments/479-pertarget-sharpen/chain479b.sh > experiments/479-pertarget-sharpen/work/chain479b.nohup 2>&1 &
