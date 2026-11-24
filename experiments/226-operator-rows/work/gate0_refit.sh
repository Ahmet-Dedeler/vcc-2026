set -x
cd ~/Code/virtual-cell-challenge-2026
D=experiments/226-operator-rows/work
uv run python scripts/build_union_whitener.py \
  --basis-prune essential --max-rss-gib 11 \
  --scratch $D/scratch --out $D/refit_pruned.npz \
  --expect-sha256 e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123
echo "EXIT_PRUNED=$?"
uv run python scripts/build_union_whitener.py \
  --basis-prune none --max-rss-gib 11 \
  --scratch $D/scratch --out $D/refit_union.npz \
  --expect-sha256 5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3
echo "EXIT_UNION=$?"
