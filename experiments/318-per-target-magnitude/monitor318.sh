#!/usr/bin/env bash
# 318 — poll S3 for each run's agg_results.csv, pull it down, then TERMINATE that worker.
# Same shape as 317's monitor: the worker syncs `out/<name>/` to S3 immediately after its
# `cell-eval2 run`, before the ~4 h `prep-real-bundle` step this experiment does not need.
# Terminating on arrival is what keeps the sweep inside its $12 cap. (bash 3.2: no assoc arrays)
set -uo pipefail
cd ~/Code/virtual-cell-challenge-2026
B=<your-s3-bucket>
OUT=experiments/318-per-target-magnitude/cloud-output
PAIRS="s318-ptrule:ptrule s318-ptshuf1:ptshuf1 s318-ptshuf2:ptshuf2 s318-ptshuf3:ptshuf3 \
s318-ptamp6:ptamp6 s318-ptunif:ptunif"
N=$(echo $PAIRS | wc -w | tr -d ' ')
for i in $(seq 1 260); do
  DONE=0; REMAIN=""
  for P in $PAIRS; do
    R="${P%%:*}"; A="${P##*:}"
    if [ -f "$OUT/$R/$A/run/agg_results.csv" ]; then DONE=$((DONE+1)); continue; fi
    if aws s3 ls "s3://$B/output/cloud-score/$R/$A/run/agg_results.csv" >/dev/null 2>&1; then
      mkdir -p "$OUT/$R"
      aws s3 sync "s3://$B/output/cloud-score/$R/" "$OUT/$R/" --exclude '*.h5ad' --only-show-errors
      ID=$(aws ec2 describe-instances --region us-east-1 \
        --filters "Name=tag:Name,Values=vcc-cloud-score-$R" \
                  Name=instance-state-name,Values=running,pending \
        --query 'Reservations[].Instances[].InstanceId' --output text 2>/dev/null)
      if [ -n "$ID" ] && [ "$ID" != "None" ]; then
        aws ec2 terminate-instances --region us-east-1 --instance-ids $ID >/dev/null 2>&1 \
          && echo "$(date -u +%FT%TZ) $R: collected, terminated $ID"
      else
        echo "$(date -u +%FT%TZ) $R: collected (instance already gone)"
      fi
      DONE=$((DONE+1))
    else
      REMAIN="$REMAIN $R"
    fi
  done
  if [ "$DONE" -ge "$N" ]; then echo "$(date -u +%FT%TZ) ALL $N COLLECTED"; break; fi
  [ $((i % 10)) -eq 1 ] && echo "$(date -u +%FT%TZ) $DONE/$N done; waiting on:$REMAIN"
  sleep 60
done
echo "$(date -u +%FT%TZ) monitor exit"
