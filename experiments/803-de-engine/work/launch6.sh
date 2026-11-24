#!/bin/bash
cd ~/Code/virtual-cell-challenge-2026
L=experiments/800-checker/cloud/labcloud.sh
bash $L launch --name de803-b6 --jobs experiments/803-de-engine/cloud/b6/jobs.txt --type r7i.2xlarge --minutes 120 --parallel 1 --sync experiments/323-current-family-validation/work --sync experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz --sync experiments/528-callprop-block/work/keys528.npz --sync data/donors/k562_gwps --sync data/donors/orion_hct116 --sync data/donors/orion_hek293t --sync experiments/803-de-engine/stage --sync experiments/803-de-engine/stage2 --collect experiments/803-de-engine/work/out --code-extra experiments/803-de-engine/arms --code-extra experiments/803-de-engine/stage/SHA256SUMS.inputs --code-extra experiments/803-de-engine/stage/SHA256SUMS.small --code-extra experiments/803-de-engine/stage/SHA256SUMS.macblocks --code-extra experiments/803-de-engine/cloud/b6
R=$(grep -o '[0-9T]*Z-de803-b6' experiments/800-checker/cloud/runs.jsonl | tail -1)
bash $L wait $R > experiments/803-de-engine/work/wait-b6.log 2>&1
