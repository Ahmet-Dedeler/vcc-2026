#!/usr/bin/env bash
# 317 — poll S3 for each run's agg_results.csv, pull it down, then TERMINATE that worker.
# The worker syncs `out/<name>/` to S3 immediately after its `cell-eval2 run`, before the
# ~4 h `prep-real-bundle` step this experiment does not need. Terminating on arrival is what
# keeps the sweep inside its $25 cap. (No associative arrays: macOS ships bash 3.2.)
set -uo pipefail
cd ~/Code/virtual-cell-challenge-2026
B=<your-s3-bucket>
OUT=experiments/317-operating-point/cloud-output
PAIRS="ka1-317:base s317-od1:od1 s317-od8:od8 s317-ks010:ks010 s317-ts070:ts070 \
s317-kt16000:kt16000 s317-amp070:amp070 s317-amp085:amp085 s317-amp120:amp120 \
s317-amp145:amp145 s317-amp200:amp200"
N=$(echo $PAIRS | wc -w | tr -d ' ')
for i in $(seq 1 200); do
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
