#!/bin/bash
# 618 stage A, RESUMED at the KA panel's operator (19:52Z): the 14 GB operator ticket headed gate602's FIFO queue for 20 min
# and blocked the council's smaller jobs, so big steps now enter through polite618.sh. Same steps and order as chain618a
# (ka panel first; ka618 must PASS before any tv step); logs to the same chain618a.log. Every exit path writes END.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; L=$X/work/logs/chain618a.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
log "chain618a2 start (resume at ka operator)"
step(){ local p=$1 s=$2 gb=$3 how=$4
  case $how in
    none) /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1 ;;
    gate) bash experiments/602-council-gate/gate602.sh "618-$p-$s" $gb -- /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1 ;;
    polite) bash $X/polite618.sh "618-$p-$s" $gb -- /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1 ;;
  esac
  local rc=$?; log "$p $s rc=$rc"; [ $rc -eq 0 ] || exit 1; }
commit_inputs(){ local p=$1
  git add -f $X/work/$p/inputs.json $X/work/$p/bundle/gene_names.csv $X/work/$p/bundle/pert_counts.csv && \
  git commit -q -m "618 $p: frozen inputs (run618 freeze) committed before any measurement" -- $X/work/$p/inputs.json $X/work/$p/bundle/gene_names.csv $X/work/$p/bundle/pert_counts.csv
  log "$p inputs committed rc=$? $(git log -1 --format=%h -- $X/work/$p/inputs.json)"; }
step ka operator 14 polite
step ka compose 6 gate
uv run python $X/ka618.py >> $L 2>&1; rc=$?; log "ka618 rc=$rc"; [ $rc -eq 0 ] || exit 1
step tv cd4 0 none
step tv freeze 0 none
commit_inputs tv
step tv stage 8 gate
step tv operator 14 polite
step tv compose 6 gate
log "stage A done"
