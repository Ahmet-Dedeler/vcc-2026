#!/bin/zsh
cd ~/Code/virtual-cell-challenge-2026/experiments/268-broad16-context-residuals/work
for L in T3M5 C4I SKGII KMRC20 ONE58 SNU8 HKGZCC U343 RPMI7951 HCT15 KYSE140 KYSE450; do
  case $L in
    KYSE450) GE=KYSE450_all_GE_filtered_feature_bc_matrix;;
    HCT15|HKGZCC|KYSE140|RPMI7951|U343) GE=${L}_combined_GE_filtered_feature_bc_matrix;;
    *) GE=${L}_GE_filtered_feature_bc_matrix;;
  esac
  EXCLUDE_ARM=1 /usr/bin/time -l uv run python build_halves2.py $L sc/$GE 2>&1 | grep -E "^\{|maximum resident|Error|Traceback|assert"
done
echo ALLDONE
