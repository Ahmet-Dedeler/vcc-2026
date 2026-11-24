#!/usr/bin/env bash
# kill the evidence S3 sync (pid 86917) at 23:56Z so the DN upload (00:01Z) has the line; re-run the same sync after "DN upload done"
while ps -p 86917 > /dev/null 2>&1; do
  if [ "$(date -u +%Y%m%d%H%M)" -ge 202610142356 ]; then pkill -P 86917 2>/dev/null; kill 86917 2>/dev/null; echo "$(date -u +%FT%TZ) killed sync 86917 (upload window)" >> experiments/700-pie-transfer/work/ev_guard700.log; exit 0; fi
  sleep 15
done
echo "$(date -u +%FT%TZ) sync 86917 finished" >> experiments/700-pie-transfer/work/ev_guard700.log
