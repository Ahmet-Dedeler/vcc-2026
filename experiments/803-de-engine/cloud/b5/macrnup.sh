#!/bin/bash
# 803 diagnosis: 630's rnup (UPRN, + DN) run on Linux over the MAC's own dumped h1-test s0 blocks (630/work/blocks,
# 435's construction). If the stats equal 630's Mac stats bit for bit, the UP pass is platform-deterministic and the
# Linux-vs-Mac UPRN gap comes from its input blocks (the ±1-count emission differences); otherwise the pass itself differs.
set -euo pipefail
X=experiments/803-de-engine
S=$X/work/scratch/macsrc
tar -xf $X/stage2/macblocks803.tar
sha256sum -c $X/stage/SHA256SUMS.macblocks > /dev/null && echo "mac blocks sha OK"
mkdir -p $S/blocks
mv experiments/630-rank-neutral-up/work/blocks/base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_*.npz $S/blocks/
cp experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz $S/
uv run python $X/emit803.py rnup h1test $X/work/scratch/macrnup 0 b_k100_q20r3p10 b_k100_q20r3p10_D200 \
  --src $S --name mac --tag 'base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse'
