#!/usr/bin/env bash
# Arm a4 — the ordering channel with n_pred held fixed.
#
# a3 is `--tail-spread 0.25 --keep-top 2000 --keep-spread 0.0`. It does two things at once:
# it promotes the top-2,000 up the confidence order (worth +0.0058 on the repaired harness)
# AND, because only 37.9% of them were already significant, it can grow the call set by up
# to 31%. Inverting the board's three fid points shows coverage is already saturated at v11,
# so that growth costs pure precision -- about -0.005 overall. The two roughly cancel.
#
# a4 raises the tail just enough to shed the calls the promotion adds, so `n_pred` lands back
# at v11's and the ordering effect is measured on its own.
#
#   TAIL=0.38 bash launch_a4.sh
set -Eeuo pipefail
TAIL="${TAIL:?set TAIL, computed from the probe's measured n_pred gap}"
ARM="${ARM:-a4-order-neutral}"
B="${BUCKET:-<your-s3-bucket>}"
SRC="${SRC_ID:-20260827T163833Z}"
I="${INSTANCE_ID:-i-08300b395f88f6f55}"

cid=$(aws ssm send-command --instance-ids "$I" --document-name AWS-RunShellScript \
  --parameters "commands=[\"set -Eeuo pipefail\",\"aws s3 cp s3://${B}/input/exp094/${SRC}/fastbuild.sh /root/fast094-${ARM}.sh --only-show-errors\",\"chmod +x /root/fast094-${ARM}.sh\",\"systemd-run --unit=vcc094-${ARM} --setenv=ARM=${ARM} --setenv=BUCKET=${B} --setenv=SRC_ID=${SRC} --setenv=WHITENER=shipped --setenv=EXTRA_ARGS='--tail-spread ${TAIL} --keep-top 2000 --keep-spread 0.0' /bin/bash /root/fast094-${ARM}.sh\",\"sleep 3; systemctl is-active vcc094-${ARM} || true\"]" \
  --comment "exp094 a4 n_pred-neutral ordering arm" --query 'Command.CommandId' --output text)
sleep 14
aws ssm get-command-invocation --command-id "$cid" --instance-id "$I" \
  --query '[Status,StandardOutputContent]' --output text
echo "arm=${ARM} tail=${TAIL}"
