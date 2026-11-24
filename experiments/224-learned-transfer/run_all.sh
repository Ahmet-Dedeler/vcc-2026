#!/bin/zsh
cd ~/Code/virtual-cell-challenge-2026/experiments/224-learned-transfer
for s in S5 S7 S8 P4 P5 P6; do
  echo "=== $s ==="
  uv run python ladder224.py --surface $s --tag $s > work/ladder-$s.log 2>&1 || echo "FAILED $s"
  tail -1 work/ladder-$s.log
done
for s in S7 P5; do
  echo "=== $s-perm ==="
  uv run python ladder224.py --surface $s --tag $s-perm --permute-labels > work/ladder-$s-perm.log 2>&1 || echo "FAILED $s-perm"
  tail -1 work/ladder-$s-perm.log
done
echo ALLDONE
