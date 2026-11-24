#!/bin/bash
# wait for gates633j (its gate602 ticket exits), then run dose633j through gate602. Cap: 6 h.
set -u
R=~/Code/virtual-cell-challenge-2026
X=$R/experiments/633-updose-folds/work
cd $R
END_TS=$(( $(date +%s) + 21600 ))
until /usr/bin/grep -Eq "wrote work/gate-nadig_jurkat|KA FAIL|Traceback" $X/gates633j.log 2>/dev/null; do
  [ "$(date +%s)" -gt "$END_TS" ] && { echo "$(date -u +%FT%TZ) cap"; echo "$(date -u +%FT%TZ) END"; exit 1; }
  sleep 30
done
if ! /usr/bin/grep -q "wrote work/gate-nadig_jurkat" $X/gates633j.log; then echo "$(date -u +%FT%TZ) gate build failed"; echo "$(date -u +%FT%TZ) END"; exit 1; fi
bash experiments/602-council-gate/gate602.sh 633j-dose 6 -- /usr/bin/time -l uv run python experiments/633-updose-folds/dose633j.py > $X/dose633j.log 2>&1
echo "$(date -u +%FT%TZ) dose633j rc $?"
echo "$(date -u +%FT%TZ) END"
