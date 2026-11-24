#!/bin/bash
# 630 on the hooks that ship: wait for 21dd4c's chain612e dumps (c4kgd, idkg3d; seed 0) to log KA PASS and move their
# blocks, then run the v3 UP pass on each (through gate602) and price it against its own zs base.
# Wait cap: 10 h from start (never an iteration count).
set -u
R=~/Code/virtual-cell-challenge-2026
W=$R/experiments/630-rank-neutral-up/work
D=$R/experiments/612-key-gated-hook/work
LOG=$D/logs/chain612e.log
cd $R
END_TS=$(( $(date +%s) + 36000 ))
ready() {
  [ -f "$LOG" ] || return 1
  /usr/bin/grep -Eqi "fail" "$LOG" && { echo "$(date -u +%FT%TZ) chain612e logged a FAIL; not stacking"; echo "$(date -u +%FT%TZ) END"; exit 1; }
  n=$(ls $D/h1_c4kgd/blocks 2>/dev/null | wc -l | tr -d ' '); m=$(ls $D/h1_idkg3d/blocks 2>/dev/null | wc -l | tr -d ' ')
  [ "$n" = "85" ] && [ "$m" = "85" ] && [ "$(/usr/bin/grep -Eci 'pass' "$LOG")" -ge 2 ]
}
until ready; do
  [ "$(date +%s)" -gt "$END_TS" ] && { echo "$(date -u +%FT%TZ) cap reached, dumps not ready"; echo "$(date -u +%FT%TZ) END"; exit 1; }
  sleep 60
done
echo "$(date -u +%FT%TZ) dumps ready; KA lines: $(/usr/bin/grep -Ei 'pass' "$LOG" | tr '\n' ' ' | cut -c1-300)"
G=experiments/602-council-gate/gate602.sh     # one gate602 ticket PER RUN (v3 rule)
bash $G 630-stack-c4kg 6 -- uv run python experiments/630-rank-neutral-up/rnup630.py 0 b_k100_q20r3p10 --src $D/h1_c4kgd --name c4kg > $W/rnup-c4kg-s0.log 2>&1
echo "$(date -u +%FT%TZ) c4kg rc $?"
bash $G 630-stack-idkg3 6 -- uv run python experiments/630-rank-neutral-up/rnup630.py 0 b_k100_q20r3p10 --src $D/h1_idkg3d --name idkg3 > $W/rnup-idkg3-s0.log 2>&1
echo "$(date -u +%FT%TZ) idkg3 rc $?"
bash $G 630-stack-price 4 -- bash -c "uv run python experiments/630-rank-neutral-up/price630.py --base $D/h1_c4kgd --name c4kg b_k100_q20r3p10 > $W/price630-c4kg.log 2>&1; uv run python experiments/630-rank-neutral-up/price630.py --base $D/h1_idkg3d --name idkg3 b_k100_q20r3p10 > $W/price630-idkg3.log 2>&1"
echo "$(date -u +%FT%TZ) stack rc $?"
echo "$(date -u +%FT%TZ) END"
