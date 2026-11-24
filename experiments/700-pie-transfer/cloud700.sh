#!/usr/bin/env bash
# 700: launch ONE disposable EC2 worker that builds PIE's evidence cache (too big for this Mac: ~55 GB of accumulators) and
# runs PIE on the k562 / h1test / vcc26abc queries. Modelled on scripts/cloud_score.sh: same bucket, IAM profile, AL2023
# AMI, terminate-on-shutdown, a user-data dead-man `shutdown -h +MAX_MINUTES`, worker sha-checked before it runs.
# usage: bash experiments/700-pie-transfer/cloud700.sh            (prints the run id, instance id and output prefix)
#        bash experiments/700-pie-transfer/cloud700.sh --collect <run-id>   (pulls outputs into work/cloud/<run-id>/)
set -Eeuo pipefail
REPO=~/Code/virtual-cell-challenge-2026; X=$REPO/experiments/700-pie-transfer
REGION=us-east-1; BUCKET=<your-s3-bucket>; IAM_PROFILE=vcc-scoring-20260822
INSTANCE_TYPE="${INSTANCE_TYPE:-r7i.8xlarge}"; ROOT_GB=300; MAX_MINUTES="${MAX_MINUTES:-300}"
if [ "${1:-}" = "--collect" ]; then
  RID="$2"; D="$X/work/cloud/$RID"; mkdir -p "$D"
  aws s3 cp "s3://$BUCKET/output/700-pie/$RID/STATUS" "$D/STATUS" --region $REGION --only-show-errors || true
  aws s3 cp "s3://$BUCKET/output/700-pie/$RID/worker.log" "$D/worker.log" --region $REGION --only-show-errors || true
  aws s3 sync "s3://$BUCKET/output/700-pie/$RID/out/" "$D/out/" --region $REGION --only-show-errors || true
  cat "$D/STATUS" 2>/dev/null; tail -5 "$D/worker.log" 2>/dev/null; exit 0
fi
RID="$(date -u +%Y%m%dT%H%M%SZ)"; PREFIX="input/700-pie/$RID"; OUT="output/700-pie/$RID"
ST="$X/work/cloud/$RID/stage/inputs"; mkdir -p "$ST/prep"
cp "$X/infer700.py" "$ST/"
cp -R "$REPO/data/external/pie/prep/vcc26abc" "$REPO/data/external/pie/prep/h1test" "$ST/prep/"
cp -R "$REPO/data/external/pie/local_sources/context_text" "$ST/context_text"
cp "$X/work/query_k562.json" "$X/work/query_h1test.json" "$X/work/query_abc.json" "$ST/"
( cd "$X/work/cloud/$RID/stage" && tar czf ../inputs.tar.gz inputs )
WSHA="$(shasum -a 256 "$X/cloud700_worker.sh" | cut -d' ' -f1)"
aws s3 cp "$X/cloud700_worker.sh" "s3://$BUCKET/$PREFIX/worker.sh" --region $REGION --only-show-errors
aws s3 cp "$X/work/cloud/$RID/inputs.tar.gz" "s3://$BUCKET/$PREFIX/inputs.tar.gz" --region $REGION --only-show-errors
AMI="$(aws ssm get-parameter --region $REGION --name /aws/service/ami-amazon-linux-latest/al2023-ami-kernel-default-x86_64 \
  --query 'Parameter.Value' --output text)"
USER_DATA="$(cat <<EOF
#!/bin/bash
shutdown -h +${MAX_MINUTES}
export AWS_REGION=${REGION}
aws s3 cp s3://${BUCKET}/${PREFIX}/worker.sh /root/worker.sh --only-show-errors
printf '${WSHA}  /root/worker.sh\n' | sha256sum --check - || shutdown -h now
systemd-run --unit=vcc-pie700 --collect --setenv=AWS_REGION=${REGION} --setenv=BUCKET=${BUCKET} \
  --setenv=PREFIX=${PREFIX} --setenv=OUT=${OUT} /bin/bash /root/worker.sh
EOF
)"
IID="$(aws ec2 run-instances --region $REGION --image-id "$AMI" --instance-type "$INSTANCE_TYPE" --count 1 \
  --iam-instance-profile "Name=${IAM_PROFILE}" --instance-initiated-shutdown-behavior terminate \
  --block-device-mappings "[{\"DeviceName\":\"/dev/xvda\",\"Ebs\":{\"VolumeSize\":${ROOT_GB},\"VolumeType\":\"gp3\",\"Throughput\":500,\"Iops\":6000,\"DeleteOnTermination\":true}}]" \
  --metadata-options "HttpTokens=required,HttpEndpoint=enabled" \
  --tag-specifications "ResourceType=instance,Tags=[{Key=Name,Value=vcc-pie700-$RID},{Key=project,Value=vcc2026},{Key=purpose,Value=pie-evidence-infer},{Key=reap-after-minutes,Value=${MAX_MINUTES}}]" \
  --user-data "$USER_DATA" --query 'Instances[0].InstanceId' --output text)"
printf '{"run_id":"%s","instance_id":"%s","instance_type":"%s","ami":"%s","max_minutes":%s,"worker_sha256":"%s","output":"s3://%s/%s/"}\n' \
  "$RID" "$IID" "$INSTANCE_TYPE" "$AMI" "$MAX_MINUTES" "$WSHA" "$BUCKET" "$OUT" | tee "$X/work/cloud/$RID/launch.json"
