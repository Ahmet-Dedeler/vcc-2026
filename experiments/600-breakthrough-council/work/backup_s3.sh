#!/bin/bash
# 600 [4d0207]: off-machine backup of the data the October recipe cannot rebuild quickly (SYNTHESIS §2 item 6).
# No backup existed (Time Machine has no destination, GitHub holds only code and small files). This copies the tier-1 set
# to the project's private S3 bucket (SSE AES256, public access blocked, no lifecycle rule on this prefix) as STANDARD_IA.
# - Bandwidth is capped at 2 MB/s through a PRIVATE AWS config file (the global ~/.aws/config is not touched).
# - Hard stop at STOP_AT (default 22:30Z), well before the 00:01Z slot uploads. `aws s3 sync` is idempotent, so a re-run
#   resumes where it stopped.
# - Writes a sha256 manifest first (git-tracked), so a restore can be verified file by file.
# usage: nohup bash experiments/600-breakthrough-council/work/backup_s3.sh [tier1|tier2] >/dev/null 2>&1 &
set -u
cd ~/Code/virtual-cell-challenge-2026 || exit 1
W=experiments/600-breakthrough-council/work
LOG=$W/backup-s3-1002.log
TIER=${1:-tier1}
BUCKET=<your-s3-bucket>
PREFIX=backup-20261002
STOP_AT=${STOP_AT:-202610022230}
ts() { date -u +%FT%TZ; }
echo "$$" > $W/backup_s3.pid
echo "$(ts) START $TIER (pid $$, stop at $STOP_AT UTC, cap 2 MB/s) -> s3://$BUCKET/$PREFIX/" >> $LOG

if [ "$TIER" = tier1 ]; then
  SRCS="data/challenge data/external/perturbmatch data/external/cd4t/cloud-compact/cd4_compact data/external/ccle_broad_2025
experiments/223-silent-zero-audit/work/tables/dld1_lowmoi experiments/166-donor-acquisition/donors/nadig_hepg2
experiments/225-source-sweep/work/tables experiments/120-myllia-a549/work/compact-donor
experiments/402-ctxA-jurkat-knn/work/knn25o400_A_jurkat50_prior_challenge.npz
experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.npz data/identity_refs data/donors"
else
  SRCS="data/harness experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
fi

CFG=$(mktemp -t awscfg600)
sed 's/^    multipart_chunksize = 16MB$/    multipart_chunksize = 16MB\
    max_bandwidth = 2MB\/s/' ~/.aws/config > "$CFG"
grep -q "max_bandwidth = 2MB/s" "$CFG" || { echo "$(ts) FAIL could not write the bandwidth cap; END" >> $LOG; exit 1; }
export AWS_CONFIG_FILE="$CFG"

MAN=$W/backup-s3-manifest-$TIER.txt
if [ ! -s "$MAN" ]; then
  echo "$(ts) manifest: hashing" >> $LOG
  for s in $SRCS; do
    if [ -d "$s" ]; then find "$s" -type f ! -name '.DS_Store' -print0 | xargs -0 shasum -a 256
    elif [ -f "$s" ]; then shasum -a 256 "$s"; fi
  done > "$MAN.tmp" && mv "$MAN.tmp" "$MAN"
  echo "$(ts) manifest: $(wc -l < "$MAN" | tr -d ' ') files -> $MAN" >> $LOG
fi

rc_all=0
for s in $SRCS; do
  now=$(date -u +%Y%m%d%H%M)
  if [ "$now" -ge "$STOP_AT" ]; then echo "$(ts) STOP_AT reached before $s; re-run later to resume" >> $LOG; rc_all=3; break; fi
  left=$(( ( $(date -j -u -f %Y%m%d%H%M "$STOP_AT" +%s) - $(date -u +%s) ) ))
  echo "$(ts) sync $s (time left ${left}s)" >> $LOG
  if [ -d "$s" ]; then
    timeout "$left" aws s3 sync "$s" "s3://$BUCKET/$PREFIX/$s" --storage-class STANDARD_IA --exclude '.DS_Store' --only-show-errors >> $LOG 2>&1
  else
    timeout "$left" aws s3 cp "$s" "s3://$BUCKET/$PREFIX/$s" --storage-class STANDARD_IA --only-show-errors >> $LOG 2>&1
  fi
  rc=$?
  echo "$(ts) done $s rc=$rc" >> $LOG
  [ $rc -ne 0 ] && rc_all=$rc
done
rm -f "$CFG"
n=$(aws s3 ls "s3://$BUCKET/$PREFIX/" --recursive --summarize 2>/dev/null | awk '/Total Objects/ {o=$3} /Total Size/ {s=$3} END {printf "%s objects, %.2f GB", o, s/1e9}')
echo "$(ts) remote now: $n" >> $LOG
echo "$(ts) END $TIER rc=$rc_all" >> $LOG
