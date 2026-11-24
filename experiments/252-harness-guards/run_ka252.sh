#!/bin/zsh
# 252 -- the known-answer suite, strictly serial, every step behind the memory gate and under
# /usr/bin/time -l.  Logs in work/<step>.log, the ledger in work/ka-status.txt.
#
#     nohup zsh run_ka252.sh > work/run_ka252.nohup 2>&1 &
set -u
cd ~/Code/virtual-cell-challenge-2026/experiments/252-harness-guards
STATUS=work/ka-status.txt
: > "$STATUS"

run() {
  local name=$1; shift
  uv run python gate252.py --wait >> work/gate.log 2>&1
  echo "[$(date -u +%H:%M:%SZ)] START $name  ($*)" >> "$STATUS"
  /usr/bin/time -l uv run python "$@" > "work/$name.log" 2>&1
  local rc=$?
  local rss
  rss=$(grep "maximum resident" "work/$name.log" | awk '{print $1}')
  local secs
  secs=$(grep " real " "work/$name.log" | awk '{print $1}')
  echo "[$(date -u +%H:%M:%SZ)] END   $name  exit=$rc  real=${secs}s  maxrss_gib=$(python3 -c "print(round(${rss:-0}/2**30,2))")" >> "$STATUS"
}

run ka216 ka216_252.py
run ka218 ka218_252.py
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do run "ka246-$f" ka246_252.py --fold "$f"; done
run ka_cache ka_cache252.py
run repro222 repro222_copy.py
for f in k562_gwps orion_hct116 orion_hek293t kolf2; do run "repro249-$f" repro249_copy.py --fold "$f"; done
run repro249-summary repro249_copy.py --summary
run ka233 ka233_252.py
echo "[$(date -u +%H:%M:%SZ)] ALL DONE" >> "$STATUS"
