#!/bin/bash
cd ~/Code/virtual-cell-challenge-2026
L=experiments/800-checker/cloud/labcloud.sh
bash $L launch --name de803-b7 --jobs experiments/803-de-engine/cloud/b7/jobs.txt --type r7i.large --minutes 40 --parallel 1 --sync experiments/803-de-engine/outmerged --collect experiments/803-de-engine/work/out --code-extra experiments/803-de-engine/arms --code-extra experiments/803-de-engine/cloud/b7 --code-extra experiments/618-h1-trainval-surface/work/gate618.json --code-extra experiments/619-board-multipliers-618/work/score619.json
R=$(grep -o '[0-9T]*Z-de803-b7' experiments/800-checker/cloud/runs.jsonl | tail -1)
bash $L wait $R > experiments/803-de-engine/work/wait-b7.log 2>&1
