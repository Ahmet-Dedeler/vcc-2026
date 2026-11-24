#!/usr/bin/env bash
# Generic exp094-tree launcher that also threads VCC_ALPHA, which launch_nb.sh does not.
#
# VCC_ALPHA (default 1.0796) multiplies the effect in LOG space before whitening. Experiment
# 096 measured that our emitted cells are 0.67x as dispersed as the reference's, so we forfeit
# most of the jackknife correction the expression-error member allows; the metric's optimum is
# therefore the rho cliff at s* = sqrt(claim/B), NOT the least-squares X/A.
#
#   ARM=a11-alpha0838 SPREAD=0.15 PHI=16.0 ALPHA=0.838 bash launch_alpha.sh
#
# TARGET_F (optional) sets VCC_CD4_TARGET_F, which makes the builder solve CD4_FRACTION so
# that the post-whitening residual lands at experiment 099's `f` -- the added residual's
# median row norm over the whitened base's. The two knobs are NOT the same scale (099
# §5.4), so hand-picking a fraction is the error that section warns about. Unset threads
# an empty string, which the builder reads as "not requested" and ignores.
#
#   ARM=a15-cd4f010 SPREAD=1.0 PHI=4.0 ALPHA=0.982 TARGET_F=0.10 \
#     EXTRA_ARGS="--keep-top 8000 --keep-spread 0.25 --tail-spread 1.0 --overdispersion 4.0" \
#     bash launch_alpha.sh
set -Eeuo pipefail
PHI="${PHI:?set PHI}"; SPREAD="${SPREAD:?set SPREAD}"; ALPHA="${ALPHA:?set ALPHA}"
ARM="${ARM:?set ARM}"
B="${BUCKET:-<your-s3-bucket>}"
SRC="${SRC_ID:-20260827T183341Z}"
I="${INSTANCE_ID:-i-08300b395f88f6f55}"
WH="${WHITENER:-pruned}"
EXTRA="${EXTRA_ARGS:---tail-spread ${SPREAD} --overdispersion ${PHI}}"
TF="${TARGET_F:-}"
cid=$(aws ssm send-command --instance-ids "$I" --document-name AWS-RunShellScript \
  --parameters "commands=[\"set -Eeuo pipefail\",\"systemctl reset-failed vcc096-${ARM} 2>/dev/null || true\",\"aws s3 cp s3://${B}/input/exp094/${SRC}/fastbuild.sh /root/fast096-${ARM}.sh --only-show-errors\",\"chmod +x /root/fast096-${ARM}.sh\",\"systemd-run --unit=vcc096-${ARM} --setenv=ARM=${ARM} --setenv=BUCKET=${B} --setenv=SRC_ID=${SRC} --setenv=WHITENER=${WH} --setenv=VCC_ALPHA=${ALPHA} --setenv=VCC_CD4_TARGET_F=${TF} --setenv=EXTRA_ARGS='${EXTRA}' /bin/bash /root/fast096-${ARM}.sh\",\"sleep 3; systemctl is-active vcc096-${ARM} || true\"]" \
  --comment "exp096 amplitude-to-cliff arm" --query 'Command.CommandId' --output text)
sleep 14
aws ssm get-command-invocation --command-id "$cid" --instance-id "$I" \
  --query '[Status,StandardOutputContent]' --output text
echo "arm=${ARM} spread=${SPREAD} phi=${PHI} VCC_ALPHA=${ALPHA} VCC_CD4_TARGET_F=${TF:-<unset>} whitener=${WH}"
