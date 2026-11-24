#!/bin/bash
# 614 faithful-key re-run (TRUSTGATE.md 3d851ae7; coordinator 20:25Z). Queued behind 618's KA / ref_tv, any live
# 655 / 656 gate jobs, and 618's magnitude stacks (the STACKS-DONE sentinel): fkey614.py (bundle, cd4, freeze, commit, stage, operator, compose), then the same 12 gate runs with
# only the key changed (run614 k562f / k562fcp), then gate614f.py. One gate602 ticket per heavy run; END on every exit.
set -u
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO"
X=experiments/614-k562-real-emitter; L=$X/work/logs/chain614f.log
Y=experiments/618-h1-trainval-surface
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$L"; }
trap 'log "END"' EXIT
MDE="--mde-c 0.1 --mde-cap 0.1 --mde-min-count 0.3 --mde-coarse"
log "chain614f start; waiting for 618's KA and ref_tv, then for live 655/656 gate jobs"
until /usr/bin/grep -q "ka618 rc=" $Y/work/logs/chain618a.log || ! ps -p $(cat $Y/work/chain618a.pid) >/dev/null 2>&1; do sleep 60; done
/usr/bin/grep -q "ka618 rc=0" $Y/work/logs/chain618a.log || { log "ka618 did not pass: the faithful-key code path is unvalidated; stopping"; exit 1; }
P=$(cat $Y/work/chain618b.pid); while ps -p "$P" >/dev/null 2>&1; do sleep 60; done
busy(){ for f in experiments/602-council-gate/work/queue/* experiments/602-council-gate/work/slot*.lock/owner; do
          [ -f "$f" ] && head -1 "$f" | /usr/bin/grep -Eq "(^| )65[56]" && return 0; done; return 1; }
while busy; do sleep 120; done
# coordinator 10-02 ~20:30Z: 618's two magnitude stacks (after its gate and the 613d re-check) go BEFORE this re-run.
# 618's stacks chain writes work/STACKS-DONE when it ends (pass, fail or skip), so this waits for that sentinel.
log "waiting for 618's STACKS-DONE sentinel"
until [ -f $Y/work/STACKS-DONE ]; do sleep 120; done
while busy; do sleep 120; done
log "queue clear of 618 KA / ref_tv and 655/656; building the faithful key"
step(){ local s=$1 gb=$2
  if [ "$gb" = "0" ]; then /usr/bin/time -l uv run python $X/fkey614.py $s > $X/work/logs/fkey-$s.log 2>&1
  else bash experiments/602-council-gate/gate602.sh "614f-$s" $gb -- /usr/bin/time -l uv run python $X/fkey614.py $s > $X/work/logs/fkey-$s.log 2>&1; fi
  local rc=$?; log "fkey $s rc=$rc"; [ $rc -eq 0 ] || exit 1; }
step bundle 0
step cd4 0
step freeze 0
git add -f $X/work/fkey/inputs.json $X/work/fkey/bundle/gene_names.csv $X/work/fkey/bundle/pert_counts.csv && \
  git commit -q -m "614 fkey: frozen inputs (fkey614 freeze) committed before any measurement" -- $X/work/fkey/inputs.json $X/work/fkey/bundle/gene_names.csv $X/work/fkey/bundle/pert_counts.csv
log "fkey inputs committed rc=$?"
step stage 8
step operator 14
step compose 6
run(){ local mode=$1 nm=$2 tag=$3; shift 3
  log "run $tag"
  bash experiments/602-council-gate/gate602.sh "614-$tag" 8 -- uv run python $X/run614.py $mode $nm base "$@" > $X/work/logs/$tag.log 2>&1
  log "done $tag rc=$?"; }
for s in 0 1 2; do
  run k562f g1f g1f-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
  run k562f g2nomdef g2nomdef-s$s $s --bulk base
  run k562f g2mdef g2mdef-s$s $s --bulk base $MDE
  run k562fcp g3cpf g3cpf-s$s $s --bulk base $MDE --zero-repair -5 --zs dn_k250_zs
done
uv run python $X/gate614f.py >> $L 2>&1; log "gate614f rc=$?"
