#!/usr/bin/env bash
# Arm a7 — over-dispersed emission, the one lever experiment 095 found that moves the order.
#
#   SPREAD=0.25 PHI=2.0 bash launch_nb.sh
#
# `spread` is the fraction of each gene's panel total that is scattered at all; `phi` is the
# gamma-Poisson dispersion of that scatter (Var = mu + phi*mu^2). Both are free on `pds`,
# `mse` and `nmae` -- the exact-total repair holds each gene's panel total whatever the
# scatter, board-verified in v14 where `pds` moved by exactly 0.000000.
#
# The operating point must be chosen so `n_pred` stays near production's ~5,961 on the
# challenge. Over-dispersion at full spread collapses it to ~1,600, and the challenge's
# coverage term -- saturated at v11 -- charges about -0.027 for that, more than reach gains.
set -Eeuo pipefail
PHI="${PHI:?set PHI}"
SPREAD="${SPREAD:-0.25}"
ARM="${ARM:-a7-nb-phi${PHI}}"
B="${BUCKET:-<your-s3-bucket>}"
SRC="${SRC_ID:-20260827T183341Z}"
I="${INSTANCE_ID:-i-08300b395f88f6f55}"
WH="${WHITENER:-pruned}"

cid=$(aws ssm send-command --instance-ids "$I" --document-name AWS-RunShellScript \
  --parameters "commands=[\"set -Eeuo pipefail\",\"aws s3 cp s3://${B}/input/exp094/${SRC}/fastbuild.sh /root/fast094-${ARM}.sh --only-show-errors\",\"chmod +x /root/fast094-${ARM}.sh\",\"systemd-run --unit=vcc094-${ARM} --setenv=ARM=${ARM} --setenv=BUCKET=${B} --setenv=SRC_ID=${SRC} --setenv=WHITENER=${WH} --setenv=EXTRA_ARGS='--tail-spread ${SPREAD} --overdispersion ${PHI}' /bin/bash /root/fast094-${ARM}.sh\",\"sleep 3; systemctl is-active vcc094-${ARM} || true\"]" \
  --comment "exp094 over-dispersed emission arm" --query 'Command.CommandId' --output text)
sleep 14
aws ssm get-command-invocation --command-id "$cid" --instance-id "$I" \
  --query '[Status,StandardOutputContent]' --output text
echo "arm=${ARM} spread=${SPREAD} phi=${PHI} whitener=${WH}"
