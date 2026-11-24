#!/bin/bash
# 600 [4d0207]: tier 2b = experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad (16.8 GB; AGENTS §10 says keep it;
# public but slow to re-download). One `aws s3 cp` (multipart, 16 MB parts) at 4 concurrent requests and a 3 MB/s cap
# (the first harness pass showed that 20 starved connections time out). Hard stop 22:30Z, before the 00:01Z uploads.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
W=experiments/600-breakthrough-council/work; LOG=$W/backup-s3-1002-tier2.log
BUCKET=<your-s3-bucket>; PREFIX=backup-20261002; STOP_AT=${STOP_AT:-202610022230}
ts(){ date -u +%FT%TZ; }
echo "$(ts) START tier2b experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad (pid $$, 4 concurrent, cap 3 MB/s, stop $STOP_AT)" >> $LOG
CFG=$(mktemp -t awscfg600r)
sed -e 's/^    max_concurrent_requests = 20$/    max_concurrent_requests = 4/' -e 's/^    multipart_chunksize = 16MB$/    multipart_chunksize = 16MB\
    max_bandwidth = 3MB\/s/' ~/.aws/config > "$CFG"
grep -q "max_bandwidth = 3MB/s" "$CFG" && grep -q "max_concurrent_requests = 4" "$CFG" || { echo "$(ts) FAIL config; END tier2b rc=9" >> $LOG; exit 1; }
export AWS_CONFIG_FILE="$CFG"
left=$(( $(date -j -u -f %Y%m%d%H%M "$STOP_AT" +%s) - $(date -u +%s) ))
timeout "$left" aws s3 cp experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad "s3://$BUCKET/$PREFIX/experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad" --storage-class STANDARD_IA --only-show-errors >> $LOG 2>&1
rc=$?
rm -f "$CFG"
echo "$(ts) END tier2b rc=$rc" >> $LOG
