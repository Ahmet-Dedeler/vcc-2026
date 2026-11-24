#!/bin/bash
# 803 instance setup (run from the repo root on EC2, after labworker's sync): unpack the staged inputs, re-splice ref_tv,
# and check every input against the Mac's sha256 lists. Any mismatch aborts the batch.
set -euo pipefail
X=experiments/803-de-engine
command -v zstd >/dev/null || dnf install -y -q zstd >/dev/null
zstd -d -q -f $X/stage/reference323.npz.zst -o experiments/323-current-family-validation/work/reference.npz
tar --zstd -xf $X/stage/small803.tar.zst
uv run python $X/stage/splice803.py splice experiments/323-current-family-validation/work/reference.npz \
  $X/stage/ref_tv_rest.bin $X/stage/ref_tv_rest.json experiments/618-h1-trainval-surface/work/ref_tv/reference.npz
sha256sum -c $X/stage/SHA256SUMS.inputs
sha256sum -c $X/stage/SHA256SUMS.small
rm -f $X/stage/reference323.npz.zst $X/stage/small803.tar.zst $X/stage/ref_tv_rest.bin
echo "setup803 OK"
