#!/bin/sh
# KA-5: 233 §5's own Phase-3 rows, re-run through 322's copy of `oct233.run()`.
# 233 published `rpe1` on the two disjoint draws (real +10.49 % rel cos / +0.00271 pds;
# pess +17.66 % / +0.00245).  Reproducing them checks the ARM path, not only BASE.
set -e
uv run python breadth322.py --panels real pess --sources rpe1 --out work/ka5_322.json
