set -e
cd ~/Code/virtual-cell-challenge-2026/experiments/194-donor-reestimation
uv run python score194.py --variants v0 "kd_post@0.1" "kd_post@0.25" "kd_post@0.4" "kd_perm@0.25" "kd_perm_depth@0.25" "kd_gate80@0.25" "kdrescale@0.25" --folds k562_gwps orion_hct116 orion_hek293t kolf2 --swap k562_gwps orion_hct116 orion_hek293t --tag dose --out work/scores-dose.json >> work/score-dose.log 2>&1
uv run python score194.py --variants v0 kd_post kdrescale debias cellmean --folds k562_gwps orion_hct116 orion_hek293t --swap k562_gwps orion_hct116 orion_hek293t --truth variant --tag all3 --out work/scores-truthvar.json >> work/score-truthvar.log 2>&1
echo FINAL-DONE >> work/score-dose.log
