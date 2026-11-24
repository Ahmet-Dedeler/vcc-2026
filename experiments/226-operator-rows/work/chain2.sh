cd ~/Code/virtual-cell-challenge-2026
while pgrep -f "python.*rows226.py" > /dev/null; do sleep 20; done
uv run python experiments/226-operator-rows/rows226.py \
  --cap 18000 --fracs 1.0,0.75,0.5 --seeds 1 \
  --out work/rows226_caplift.json > experiments/226-operator-rows/work/rows226_caplift.log 2>&1
echo "CAPLIFT_EXIT=$?" >> experiments/226-operator-rows/work/rows226_caplift.log
