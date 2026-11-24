set -e
cd ~/Code/virtual-cell-challenge-2026/experiments/194-donor-reestimation
uv run python accumulate.py --donor orion_hct116 --stage sums --targets restricted >> work/chain.log 2>&1
uv run python accumulate.py --donor orion_hek293t --stage obs >> work/chain.log 2>&1
uv run python accumulate.py --donor orion_hek293t --stage sums --targets restricted >> work/chain.log 2>&1
echo CHAIN-DONE >> work/chain.log
