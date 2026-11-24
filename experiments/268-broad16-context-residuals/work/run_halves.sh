#!/bin/zsh
cd ~/Code/virtual-cell-challenge-2026/experiments/268-broad16-context-residuals/work
for L in "$@"; do
  case $L in
    KYSE450) GE=KYSE450_all_GE_filtered_feature_bc_matrix; GU=KYSE450_all_CRISPR.filt;;
    HCT15|HKGZCC|KYSE140|RPMI7951|U343) GE=${L}_combined_GE_filtered_feature_bc_matrix; GU=${L}_CRISPR.filt;;
    *) GE=${L}_GE_filtered_feature_bc_matrix; GU=${L}_CRISPR.filt;;
  esac
  /usr/bin/time -l uv run python build_halves.py $L sc/$GE sc/$GU 3 0.1 2>&1 | grep -E "^\{|h5|maximum resident|Error|Traceback|assert" 
done
