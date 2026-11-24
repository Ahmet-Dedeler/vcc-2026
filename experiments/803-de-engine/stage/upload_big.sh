#!/bin/bash
# 803 staging of the two big inputs (stops itself before the 23:30Z network window)
cd ~/Code/virtual-cell-challenge-2026
log(){ echo "$(date -u +%FT%TZ) $*"; }
log start
( sleep $(( $(date -u -j -f %Y-%m-%dT%H:%M:%SZ 2026-10-15T23:29:00Z +%s) - $(date -u +%s) )); log "deadline: killing uploads"; pkill -f "aws s3 cp - s3://<your-s3-bucket>/input/lab/repo/experiments/803-de-engine" ; pkill -f "aws s3 cp experiments/323-current-family-validation/work/whitener" ) &
zstd -T6 -6 -c experiments/323-current-family-validation/work/reference.npz | aws s3 cp - s3://<your-s3-bucket>/input/lab/repo/experiments/803-de-engine/stage/reference323.npz.zst --only-show-errors && log "reference323 zst uploaded" &
aws s3 cp experiments/323-current-family-validation/work/whitener_h1_pruned.npz s3://<your-s3-bucket>/input/lab/repo/experiments/323-current-family-validation/work/whitener_h1_pruned.npz --only-show-errors && log "whitener uploaded" &
wait %2 %3
log END
