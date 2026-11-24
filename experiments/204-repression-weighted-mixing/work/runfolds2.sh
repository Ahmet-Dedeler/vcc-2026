#!/bin/zsh
set -u
cd ~/Code/virtual-cell-challenge-2026/experiments/204-repression-weighted-mixing
ARMS="equal a_self deatten a_slope minvar deatten_slope a_inst minvar_self a_self_pt deatten_pt perm_a_1 perm_a_2 perm_a_3 perm_a_4 perm_a_5 share_4_1_1 share_2_1_1 share_1_2_2 k05 k07 k085 k115 k13 k15 h05 h07 h085 h115 h13 h15"
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do
  echo "===== fold $f ====="
  uv run python score204.py --arms ${=ARMS} --folds $f || echo "FOLD $f FAILED"
done
echo ALLDONE
