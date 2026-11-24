#!/bin/bash
cd ~/Code/virtual-cell-challenge-2026
L=experiments/800-checker/cloud/labcloud.sh
bash $L launch --name de803-b3 --jobs experiments/803-de-engine/cloud/b3/jobs.txt --type r7i.4xlarge --minutes 150 --parallel 1 --sync data/donors/k562_gwps --sync data/donors/orion_hct116 --sync data/donors/orion_hek293t --sync experiments/323-current-family-validation/work --sync experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz --sync experiments/528-callprop-block/work/keys528.npz --sync experiments/803-de-engine/stage --collect experiments/803-de-engine/work/out --code-extra experiments/803-de-engine/arms --code-extra experiments/803-de-engine/stage/SHA256SUMS.inputs --code-extra experiments/803-de-engine/stage/SHA256SUMS.small --code-extra experiments/803-de-engine/cloud/b3
R3=$(grep -o '[0-9T]*Z-de803-b3' experiments/800-checker/cloud/runs.jsonl | tail -1)
bash $L launch --name de803-b4 --jobs experiments/803-de-engine/cloud/b4/jobs.txt --type r7i.8xlarge --minutes 150 --parallel 1 --sync data/donors/k562_gwps --sync data/donors/orion_hct116 --sync data/donors/orion_hek293t --sync experiments/323-current-family-validation/work --sync experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz --sync experiments/528-callprop-block/work/keys528.npz --sync experiments/803-de-engine/stage --collect experiments/803-de-engine/work/out --code-extra experiments/803-de-engine/arms --code-extra experiments/803-de-engine/stage/SHA256SUMS.inputs --code-extra experiments/803-de-engine/stage/SHA256SUMS.small --code-extra experiments/803-de-engine/cloud/b4
R4=$(grep -o '[0-9T]*Z-de803-b4' experiments/800-checker/cloud/runs.jsonl | tail -1)
bash $L wait $R3 > experiments/803-de-engine/work/wait-b3.log 2>&1 &
bash $L wait $R4 > experiments/803-de-engine/work/wait-b4.log 2>&1 &
wait
