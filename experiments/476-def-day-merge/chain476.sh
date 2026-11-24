#!/usr/bin/env bash
# 476: 485's October day-of fixes merged onto 474 (work/def_day.resolved.sh; 485's three other patches apply clean), then
# four rehearsals with known answers. Serial after chain495 (process exit, no other build, builder == HEAD).
#   apply: scripts/def_day.sh <- work/def_day.resolved.sh; patch -p1 def_day_helpers.py, def_pipeline.py (dump-assert
#          skip), relevel_iter.py; each applied file must equal repro/*.new, and repro/repro_fixes.py must pass;
#          run_arm_generic.sh's driver pin list gains the new def_pipeline sha (atomic rename).
#   (a) efp strict: stops at setup with the CD4-gate message (474's (a));
#   (b) 465 retraced on efr (DEF_DAY_LODO=1, programme h1_2025_full, flags pt): final output == 465's 6113937f, postbuild
#       PASS, and the dump stage now exits 0;
#   (c) efp relaxed: final output == 474's (c) ac49d319 (D 0.9831 E 1.0324 F 1.0663), postbuild PASS;
#   (d) in-place rerun of (c) in the same work dir with q 1.5 and max_levels 1, (c)'s h5ads still on disk: must either
#       STOP loudly or finish DONE with a fresh build (meta q 1.5, output != ac49d319, every rho <= 1). The pre-485 script
#       would have shipped (c)'s h5ad under q 1.5's factor string.
# Any failure restores the five pre-476 files and stops. Never submits. The coordinator commits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/476-def-day-merge; W=$X/work; LOG=$W/chain476.log; P=experiments/485-cloud-october-review/patches; T=$X/repro
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; DP=scripts/def_pipeline.py; RI=experiments/448-def-rehearsal-q175/relevel_iter.py
RAG=experiments/428-zero-gene-audit/run_arm_generic.sh; BU=scripts/make_submission_pm_cd4.py
EFP=experiments/238-essential-free-panel/work/ef-bundle; EFR=experiments/238-essential-free-panel/work/ef-bundle-r
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain476.pid
log "armed (pid $$): waiting for chain495's process to exit (and no other build), then builder == HEAD"
P495="$(cat experiments/495-lodo-cap6/work/chain495.pid 2>/dev/null)"
until ! ps -p "${P495:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain495 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain479b?\.sh|chain495\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609300200 ] && { log "495 never ended by 02:00Z 09-30 -- stopping"; log "END"; exit 1; }; sleep 60; done
log "chain495 exited; its last lines: $(tail -3 experiments/495-lodo-cap6/work/chain495.log | tr '\n' ' ' | cut -c1-400)"
until [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609300200 ] && { log "builder never equal to HEAD by 02:00Z 09-30 -- stopping"; log "END"; exit 1; }
  log "builder != HEAD (495 uncommitted?) -- waiting"; sleep 300; done
for f in $DD $HP $DP $RI $RAG; do [ "$(shasum -a 256 $f | cut -d' ' -f1)" = "$(git show HEAD:$f | shasum -a 256 | cut -d' ' -f1)" ] || { log "$f differs from HEAD -- stopping"; log "END"; exit 1; }; done
for f in $DD $HP $DP $RI $RAG; do cp $f $W/pre476.$(basename $f); done
back(){ cp $W/pre476.def_day.sh $DD; cp $W/pre476.def_day_helpers.py $HP; cp $W/pre476.def_pipeline.py $DP; cp $W/pre476.relevel_iter.py $RI
        cp $W/pre476.run_arm_generic.sh $RAG.tmp476 && mv $RAG.tmp476 $RAG
        log "RESTORED the five pre-476 files ($1)"; log "END"; exit 1; }
memwait(){ local n=0 fi; while :; do
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && return 0
  n=$((n+1)); [ $n -ge 120 ] && return 1; log "waiting for memory: free+inactive $fi GiB < 30"; sleep 60; done; }
# --- apply
cp $W/def_day.resolved.sh $DD
for f in def_day_helpers.py def_pipeline.py relevel_iter.py; do patch -p1 -s < $P/$f.patch >> $LOG 2>&1 || back "$f patch failed"; done
cmp -s $DD $T/def_day.sh.new && cmp -s $HP $T/def_day_helpers.py.new && cmp -s $DP $T/def_pipeline.py.new && cmp -s $RI $T/relevel_iter.py.new \
  || back "an applied file differs from its repro/*.new copy"
bash -n $DD || back "def_day.sh syntax"
(cd $T && uv run python repro_fixes.py > last-fixes-chain.log 2>&1) || back "repro_fixes.py failed: $(tail -3 $T/last-fixes-chain.log | tr '\n' ' ')"
ND="$(shasum -a 256 $DP | cut -d' ' -f1)"
uv run python - "$RAG" "$ND" <<'PY' >> $LOG 2>&1 || back "run_arm_generic pin-list update failed"
import sys
from pathlib import Path
p, nd = Path(sys.argv[1]), sys.argv[2]; s = p.read_text()
old = "|62b5110dbec32a41a9357ec6aa68329cd3888599c81a974c15fe9323222ccb4b) ;;"
assert s.count(old) == 1, s.count(old)
s = s.replace(old, f"|62b5110dbec32a41a9357ec6aa68329cd3888599c81a974c15fe9323222ccb4b|{nd}) ;;  # 476: + 485's dump-assert skip")
tmp = p.with_name(p.name + ".tmp476"); tmp.write_text(s); tmp.rename(p)
print("run_arm_generic.sh pin list +", nd[:16])
PY
log "applied: def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16), driver ${ND:0:16}, relevel $(shasum -a 256 $RI | cut -c1-16); repro_fixes PASS; run_arm_generic pins + ${ND:0:16}"
# (a) efp strict
rm -rf $W/efp_strict
bash $DD $EFP $W/efp_strict 476efps 1.75 auto pt:noprog 4 > $W/efp_strict.nohup 2>&1; RC=$?
if [ $RC -ne 0 ] && /usr/bin/grep -q "STOP: CD4 extraction gate fired" $W/efp_strict/def_day.log; then
  log "(a) PASS: efp strict stops at setup: $(/usr/bin/grep 'STOP:' $W/efp_strict/def_day.log | head -1 | cut -c1-200)"
else back "(a) efp strict did not stop with the CD4-gate message (rc $RC): $(tail -2 $W/efp_strict/def_day.log | tr '\n' ' ')"; fi
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || back "builder not restored after (a)"
# (b) 465 retraced on efr
memwait || back "memory never reached 30 GiB before (b)"
DEF_DAY_LODO=1 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full bash $DD $EFR $W/efr_retrace 465efr 1.75 auto pt 4 > $W/efr_retrace.nohup 2>&1
/usr/bin/grep -q "DONE. Next" $W/efr_retrace/def_day.log || back "(b) retrace did not finish: $(tail -2 $W/efr_retrace/def_day.log | tr '\n' ' ')"
/usr/bin/grep -q "postbuild PASS" $W/efr_retrace/def_day.log || back "(b) postbuild did not PASS"
GOT="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $W/efr_retrace/465efr_final.h5ad.meta.json)"
[ "$GOT" = 6113937f1acafb34c8cc323c7d0a35d4fabd7fc505c7e52404745e26263a949b ] || back "(b) final output ${GOT:0:16} != 465's 6113937f"
/usr/bin/grep -q "skip" $W/efr_retrace/driver-dump.log && log "(b) driver-dump: $(/usr/bin/grep -i 'skip' $W/efr_retrace/driver-dump.log | head -1 | cut -c1-200)"
log "(b) PASS: 465 retraced byte for byte (${GOT:0:16}); $(/usr/bin/grep -E 'CONVERGED at|NOT converged|postbuild PASS' $W/efr_retrace/def_day.log | tr '\n' ' ' | cut -c1-300)"
rm -f $W/efr_retrace/*.h5ad
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || back "builder not restored after (b)"
# (c) efp relaxed
memwait || back "memory never reached 30 GiB before (c)"
DEF_DAY_CD4_THRESHOLDS=relaxed bash $DD $EFP $W/efp_relaxed 474efp 1.75 auto pt:noprog 4 > $W/efp_relaxed.nohup 2>&1
/usr/bin/grep -q "DONE. Next" $W/efp_relaxed/def_day.log || back "(c) efp relaxed did not finish: $(tail -2 $W/efp_relaxed/def_day.log | tr '\n' ' ')"
/usr/bin/grep -q "postbuild PASS" $W/efp_relaxed/def_day.log || back "(c) postbuild did not PASS"
GOTC="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $W/efp_relaxed/474efp_final.h5ad.meta.json)"
[ "$GOTC" = ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3 ] || back "(c) final output ${GOTC:0:16} != 474 (c)'s ac49d319"
log "(c) PASS: 474's efp relaxed retraced byte for byte (${GOTC:0:16}); $(/usr/bin/grep -E 'CONVERGED at|NOT converged|postbuild PASS' $W/efp_relaxed/def_day.log | tr '\n' ' ' | cut -c1-300)"
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || back "builder not restored after (c)"
# (d) in-place rerun of (c) at q 1.5, max_levels 1, with (c)'s h5ads still on disk
memwait || back "memory never reached 30 GiB before (d)"
log "(d) in-place rerun: $(ls $W/efp_relaxed/*.h5ad 2>/dev/null | wc -l | tr -d ' ') h5ad(s) from (c) on disk"
mv $W/efp_relaxed/def_day.log $W/efp_relaxed.def_day.c.log   # (d)'s log must start empty: (c)'s DONE line would satisfy the checks
DEF_DAY_CD4_THRESHOLDS=relaxed bash $DD $EFP $W/efp_relaxed 474efp 1.5 auto pt:noprog 1 > $W/efp_rerun.nohup 2>&1; RC=$?
if [ $RC -ne 0 ]; then
  /usr/bin/grep -q "STOP:" $W/efp_relaxed/def_day.log || back "(d) rerun failed rc $RC without a STOP line"
  log "(d) PASS (fail-closed): the in-place rerun stopped loudly rc $RC: $(/usr/bin/grep 'STOP:' $W/efp_relaxed/def_day.log | tail -1 | cut -c1-240)"
else
  /usr/bin/grep -q "DONE. Next" $W/efp_relaxed/def_day.log || back "(d) rc 0 without DONE"
  uv run python - $W/efp_relaxed/474efp_final.h5ad.meta.json <<'PY' >> $LOG 2>&1 || back "(d) the rerun's final is not a fresh q-1.5 build"
import json, sys
m = json.load(open(sys.argv[1])); c = m["config"]; sha = m["output"]["sha256"]
ok = abs(float(c["bulk_sharpen_q"]) - 1.5) < 1e-9 and sha != "ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3"
print(f"(d) final meta q {c['bulk_sharpen_q']} scale {c.get('bulk_ctx_scale')} output {sha[:16]}:", "fresh" if ok else "STALE")
sys.exit(0 if ok else 1)
PY
  FIN="$(/usr/bin/grep -o 'final: [^ ]* factors [^ ]*' $W/efp_relaxed/def_day.log | tail -1)"
  RHO="$(uv run python scripts/def_day_helpers.py rho "$W/efp_relaxed/$(echo "$FIN" | awk '{print $2}').ntclevel.json" 2>&1 | tail -1)"
  log "(d) PASS (fresh build): $FIN; rho $RHO; $(/usr/bin/grep -E 'CONVERGED at|NOT converged|shipping the nearest' $W/efp_relaxed/def_day.log | tail -2 | tr '\n' ' ' | cut -c1-300)"
fi
rm -f $W/efp_relaxed/*.h5ad
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || back "builder not restored after (d)"
log "END all PASS: def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16), driver ${ND:0:16}, relevel $(shasum -a 256 $RI | cut -c1-16), run_arm_generic $(shasum -a 256 $RAG | cut -c1-16) (uncommitted; the coordinator commits)"
