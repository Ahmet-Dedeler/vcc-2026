set -e
cd ~/Code/virtual-cell-challenge-2026/experiments/194-donor-reestimation
uv run python accumulate.py --donor k562_essential --stage obs >> work/chain.log 2>&1
uv run python accumulate.py --donor k562_essential --stage sums --targets restricted >> work/chain.log 2>&1
echo ESSENTIAL-DONE >> work/chain.log
