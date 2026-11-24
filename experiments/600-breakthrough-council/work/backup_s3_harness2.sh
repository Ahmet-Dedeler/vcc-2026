#!/bin/bash
# 600 [4d0207]: retry of tier 2a (data/harness). The first pass (1 MB/s cap over 20 concurrent part uploads) starved
# connections: S3 RequestTimeout on UploadPart, and 8 of 18 files never landed. This pass uses 4 concurrent requests at
# 3 MB/s (tier 1 is done, so the uplink is free). `aws s3 sync` skips files already there. Hard stop 22:30Z.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
W=experiments/600-breakthrough-council/work; LOG=$W/backup-s3-1002-tier2.log
BUCKET=<your-s3-bucket>; PREFIX=backup-20261002; STOP_AT=${STOP_AT:-202610022230}
ts(){ date -u +%FT%TZ; }
echo "$(ts) START tier2a-retry data/harness (pid $$, 4 concurrent, cap 3 MB/s, stop $STOP_AT)" >> $LOG
CFG=$(mktemp -t awscfg600r)
sed -e 's/^    max_concurrent_requests = 20$/    max_concurrent_requests = 4/' -e 's/^    multipart_chunksize = 16MB$/    multipart_chunksize = 16MB\
    max_bandwidth = 3MB\/s/' ~/.aws/config > "$CFG"
grep -q "max_bandwidth = 3MB/s" "$CFG" && grep -q "max_concurrent_requests = 4" "$CFG" || { echo "$(ts) FAIL config; END tier2a-retry rc=9" >> $LOG; exit 1; }
export AWS_CONFIG_FILE="$CFG"
left=$(( $(date -j -u -f %Y%m%d%H%M "$STOP_AT" +%s) - $(date -u +%s) ))
timeout "$left" aws s3 sync data/harness "s3://$BUCKET/$PREFIX/data/harness" --storage-class STANDARD_IA --exclude '.DS_Store' --only-show-errors >> $LOG 2>&1
rc=$?
rm -f "$CFG"
echo "$(ts) END tier2a-retry rc=$rc" >> $LOG
