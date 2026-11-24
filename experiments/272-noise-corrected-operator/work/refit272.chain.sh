#!/bin/bash
# construction B, strictly serial, ONE fit per process (E062's guard reads the process's peak RSS).
cd ~/Code/virtual-cell-challenge-2026/experiments/272-noise-corrected-operator
run() { PYTHONDONTWRITEBYTECODE=1 /usr/bin/time -l uv run python refit272.py "$@" > "work/refit272-$1-$2${3:+-$3}.nohup" 2>&1; }
run fit-full orion_hct116
run fit-half orion_hct116 A
run fit-half orion_hct116 B
run measure orion_hct116
touch work/refit272-orion_hct116.done
run fit-full k562_gwps
run fit-half k562_gwps A
run fit-half k562_gwps B
run measure k562_gwps
touch work/refit272-k562_gwps.done
run cleanup orion_hct116
run cleanup k562_gwps
touch work/refit272.chain.done
