#!/bin/bash
# 618 stages A + B, ONE serial chain through gate602 v6 (5197f00b) with quiet admission (GATE602_QUIET=1) for the three big
# steps, so no two quiet jobs ever overlap (4d0207's rule). Replaces chain618a2 / chain618b2 (their polite waits could not
# be met on this machine: the compressor holds ~16 GiB). Order: ka operator (quiet 14, measured 13.7 GB) -> ka compose ->
# ka618 (KA must PASS) -> tv cd4 / freeze / commit / stage -> tv operator (quiet 14) -> tv compose -> ref_tv (quiet 17:
# extrapolated, h1-test's measured 11.87 GB x 1.31 cells / nnz = 15.6, plus margin) -> diag618.
# Writes the lines chain618c / chain614f wait for into chain618a.log ("ka618 rc=", "stage A done") and chain618b.log
# ("ref_tv rc=0", "diag618"). Every exit path writes END to both logs.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; LA=$X/work/logs/chain618a.log; LB=$X/work/logs/chain618b.log
loga(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LA"; }
logb(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LB"; }
trap 'loga "END"; logb "END"' EXIT
for p in "$@"; do while ps -p "$p" >/dev/null 2>&1; do sleep 2; done; done   # the old chains' pids, killed by the caller
loga "chain618ab4 start (gate602 v6 quiet admission; resume at ka operator)"
logb "chain618ab4 start (ref_tv after stage A)"
step(){ local p=$1 s=$2 gb=$3 how=$4
  case $how in
    none) /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1 ;;
    gate) bash experiments/602-council-gate/gate602.sh "618-$p-$s" $gb -- /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1 ;;
    quiet) GATE602_QUIET=1 bash experiments/602-council-gate/gate602.sh "618-$p-$s" $gb -- /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1 ;;
  esac
  local rc=$?; loga "$p $s rc=$rc"; [ $rc -eq 0 ] || exit 1; }
commit_inputs(){ local p=$1
  git add -f $X/work/$p/inputs.json $X/work/$p/bundle/gene_names.csv $X/work/$p/bundle/pert_counts.csv && \
  git commit -q -m "618 $p: frozen inputs (run618 freeze) committed before any measurement" -- $X/work/$p/inputs.json $X/work/$p/bundle/gene_names.csv $X/work/$p/bundle/pert_counts.csv
  loga "$p inputs committed rc=$? $(git log -1 --format=%h -- $X/work/$p/inputs.json)"; }
step ka operator 14 quiet
step ka compose 6 gate
uv run python $X/ka618.py >> $LA 2>&1; rc=$?; loga "ka618 rc=$rc"; [ $rc -eq 0 ] || exit 1
step tv cd4 0 none
step tv freeze 0 none
commit_inputs tv
step tv stage 8 gate
step tv operator 14 quiet
step tv compose 6 gate
loga "stage A done"
GATE602_QUIET=1 bash experiments/602-council-gate/gate602.sh 618-ref 17 -- /usr/bin/time -l uv run python experiments/614-k562-real-emitter/ref614.py $X/work/tv-truth.h5 $X/work/targets_tv.txt $X/work/ref_tv > $X/work/logs/ref_tv.log 2>&1
rc=$?; logb "ref_tv rc=$rc"; [ $rc -eq 0 ] || exit 1
uv run python $X/diag618.py >> $LB 2>&1; logb "diag618 rc=$?"
