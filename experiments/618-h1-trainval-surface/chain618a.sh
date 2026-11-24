#!/bin/bash
# 618 stage A (PLAN.md): the faithful production-key pipeline. KA panel first (h1-test's 85 targets must reproduce 323's
# cap15 bit for bit, ka618.py), then the trainval panel. One gate602 ticket per heavy step. Every exit path writes END.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/618-h1-trainval-surface; L=$X/work/logs/chain618a.log
mkdir -p $X/work/logs
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
log "chain618a start"
step(){ local p=$1 s=$2 gb=$3
  if [ "$gb" = "0" ]; then
    /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1
  else
    bash experiments/602-council-gate/gate602.sh "618-$p-$s" $gb -- /usr/bin/time -l uv run python $X/run618.py $p $s > $X/work/logs/$p-$s.log 2>&1
  fi
  local rc=$?; log "$p $s rc=$rc"; [ $rc -eq 0 ] || exit 1; }
commit_inputs(){ local p=$1
  git add -f $X/work/$p/inputs.json $X/work/$p/bundle/gene_names.csv $X/work/$p/bundle/pert_counts.csv && \
  git commit -q -m "618 $p: frozen inputs (run618 freeze) committed before any measurement" -- $X/work/$p/inputs.json $X/work/$p/bundle/gene_names.csv $X/work/$p/bundle/pert_counts.csv
  log "$p inputs committed rc=$? $(git log -1 --format=%h -- $X/work/$p/inputs.json)"; }
for p in ka tv; do
  step $p cd4 0
  step $p freeze 0
  commit_inputs $p
  step $p stage 8
  step $p operator 14
  step $p compose 6
  if [ "$p" = "ka" ]; then
    uv run python $X/ka618.py >> $L 2>&1; rc=$?; log "ka618 rc=$rc"; [ $rc -eq 0 ] || exit 1
  fi
done
log "stage A done"
