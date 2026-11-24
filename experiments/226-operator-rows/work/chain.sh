cd ~/Code/virtual-cell-challenge-2026
while pgrep -f "python.*rows226.py" > /dev/null; do sleep 20; done
uv run python experiments/226-operator-rows/gamma226.py > experiments/226-operator-rows/work/gamma226.log 2>&1
echo "GAMMA_EXIT=$?" >> experiments/226-operator-rows/work/gamma226.log
