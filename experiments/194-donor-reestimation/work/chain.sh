set -e
cd ~/Code/virtual-cell-challenge-2026/experiments/194-donor-reestimation
for L in orion_hct116 orion_hek293t; do
  uv run python accumulate.py --donor $L --stage obs   >> work/chain.log 2>&1
  uv run python accumulate.py --donor $L --stage sums --targets restricted >> work/chain.log 2>&1
done
echo CHAIN-DONE >> work/chain.log
