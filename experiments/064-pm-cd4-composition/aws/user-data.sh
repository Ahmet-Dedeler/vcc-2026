#!/bin/bash
set -euxo pipefail

aws s3 cp \
  s3://<your-s3-bucket>/input/exp064-pm-cd4-composition/run.sh \
  /root/exp064-run.sh
chmod +x /root/exp064-run.sh
systemd-run \
  --unit=vcc-exp064 \
  --description="VCC experiment 064 frozen composition" \
  --collect \
  /root/exp064-run.sh

