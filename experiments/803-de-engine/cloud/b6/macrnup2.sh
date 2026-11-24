#!/bin/bash
# 803 diagnosis (b6): Linux stats + Linux rnup630 (UPRN, + DN) on the Mac's own h1-test s0 blocks.
set -euo pipefail
X=experiments/803-de-engine
S=$X/work/scratch/macsrc
tar -xf $X/stage2/macblocks803.tar 2>/dev/null
sha256sum -c $X/stage/SHA256SUMS.macblocks > /dev/null && echo "mac blocks sha OK"
mkdir -p $S/blocks
mv experiments/630-rank-neutral-up/work/blocks/base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_*.npz $S/blocks/
uv run python $X/diag803.py $S experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse_dnk250zs.npz $X/work/out/diag/zs_stats_linux_vs_mac.json
uv run python $X/emit803.py rnup h1test $X/work/scratch/macrnup 0 b_k100_q20r3p10 b_k100_q20r3p10_D200 \
  --src $S --name mac --tag 'base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse'
W=experiments/630-rank-neutral-up/work/h1; M=$X/work/scratch/macrnup
uv run python $X/de803.py ka --name mac_blocks_rnup_upr --new $M/stats-rnup_mac_s0_b_k100_q20r3p10.npz --ref $W/stats-rnup_s0_b_k100_q20r3p10.npz || true
uv run python $X/de803.py ka --name mac_blocks_rnup_updn --new $M/stats-rnup_mac_s0_b_k100_q20r3p10_D200.npz --ref $W/stats-rnup_s0_b_k100_q20r3p10_D200.npz || true
cp $M/rnup630_mac_s0_*.json $X/work/out/diag/ ; cp $M/manifest-rnup.json $X/work/out/diag/
