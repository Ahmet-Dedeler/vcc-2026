#!/bin/bash
cd ~/Code/virtual-cell-challenge-2026
echo "$(date -u +%FT%TZ) start"
( sleep $(( $(date -u -j -f %Y-%m-%dT%H:%M:%SZ 2026-10-15T23:29:00Z +%s) - $(date -u +%s) )); pkill -f "aws s3 cp /private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/dfcab675-ce16-4481-bec4-bacb2383dec4/scratchpad/small803" ) &
aws s3 cp /private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/dfcab675-ce16-4481-bec4-bacb2383dec4/scratchpad/small803.tar.zst s3://<your-s3-bucket>/input/lab/repo/experiments/803-de-engine/stage/small803.tar.zst --only-show-errors && echo "$(date -u +%FT%TZ) small uploaded"
echo "$(date -u +%FT%TZ) END"
