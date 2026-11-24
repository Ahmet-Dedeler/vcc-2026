#!/usr/bin/env bash
# 700: pause the PIE downloader (SIGSTOP) at 10-15 00:25Z if it is still running, so the coordinator's DN upload
# window (~00:30Z-01:15Z) has the line. Exits when the downloader exits or after acting.
P=21868; LOG=~/Code/virtual-cell-challenge-2026/experiments/700-pie-transfer/work/dl_guard700.log
while ps -p $P > /dev/null 2>&1; do
  if [ "$(date -u +%Y%m%d%H%M)" -ge 202610150025 ]; then kill -STOP $P && echo "$(date -u +%FT%TZ) SIGSTOP $P (upload window)" >> $LOG; exit 0; fi
  sleep 30
done
echo "$(date -u +%FT%TZ) downloader exited; guard done" >> $LOG
