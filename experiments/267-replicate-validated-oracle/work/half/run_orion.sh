#!/bin/zsh
cd ~/Code/virtual-cell-challenge-2026/experiments/267-replicate-validated-oracle
export PYTHONDONTWRITEBYTECODE=1
uv run python half267.py orion_hct116 > work/half/half267-hct116.log 2>&1
( uv run python truth267.py orion_hct116 > work/half/truth267-hct116.nohup 2>&1 && uv run python sweep267.py orion_hct116 > work/sweep267-orion_hct116.nohup 2>&1 ) &
uv run python half267.py orion_hek293t > work/half/half267-hek293t.log 2>&1
wait
uv run python truth267.py orion_hek293t > work/half/truth267-hek293t.nohup 2>&1
uv run python sweep267.py orion_hek293t > work/sweep267-orion_hek293t.nohup 2>&1
