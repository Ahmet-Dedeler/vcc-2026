#!/bin/bash
# 600 [4d0207]: tier-2a of the S3 backup = data/harness (12 GB of local scoring surfaces: h1-test, h1-full, the Orion
# harnesses and their baseline predictions; rebuildable, but they took real effort). Runs BESIDE tier 1 (backup_s3.sh is
# live and is not edited), capped at 1 MB/s through its own private AWS config. Hard stop at 22:30Z; `aws s3 sync` resumes.
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
W=experiments/600-breakthrough-council/work; LOG=$W/backup-s3-1002-tier2.log
BUCKET=<your-s3-bucket>; PREFIX=backup-20261002; STOP_AT=${STOP_AT:-202610022230}
ts(){ date -u +%FT%TZ; }
echo "$(ts) START tier2a data/harness (pid $$, cap 1 MB/s, stop $STOP_AT)" >> $LOG
CFG=$(mktemp -t awscfg600h)
sed 's/^    multipart_chunksize = 16MB$/    multipart_chunksize = 16MB\
    max_bandwidth = 1MB\/s/' ~/.aws/config > "$CFG"
grep -q "max_bandwidth = 1MB/s" "$CFG" || { echo "$(ts) FAIL bandwidth cap; END" >> $LOG; exit 1; }
export AWS_CONFIG_FILE="$CFG"
MAN=$W/backup-s3-manifest-tier2a.txt
[ -s "$MAN" ] || { find data/harness -type f ! -name .DS_Store -print0 | xargs -0 shasum -a 256 > "$MAN.tmp" && mv "$MAN.tmp" "$MAN"; echo "$(ts) manifest $(wc -l < $MAN | tr -d ' ') files" >> $LOG; }
left=$(( $(date -j -u -f %Y%m%d%H%M "$STOP_AT" +%s) - $(date -u +%s) ))
timeout "$left" aws s3 sync data/harness "s3://$BUCKET/$PREFIX/data/harness" --storage-class STANDARD_IA --exclude '.DS_Store' --only-show-errors >> $LOG 2>&1
rc=$?
rm -f "$CFG"
echo "$(ts) done data/harness rc=$rc" >> $LOG
echo "$(ts) END tier2a rc=$rc" >> $LOG
