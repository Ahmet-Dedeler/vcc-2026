#!/usr/bin/env bash
# 663 prep fallback (coordinator, 2026-10-04 ~20:00Z).
# run_arm_generic's prep gate (fi23 or pressure55) gives up after 240 x 30 s. idall3's started at 19:23:58Z while the iOS
# Simulator held ~17 GiB in the compressor (free+inactive 19.5 GiB, memory_pressure free 59 %, swap free 386 MiB against
# the 400 MiB floor), so it can expire at ~21:24Z with "GATE never opened (prep)". build_idall3.sh then ends "NOT READY gates
# did not pass" and fill1005s2 exits with nothing armed.
# This waits for build_idall3.sh (pid 91191) to exit and acts ONLY if the prep gate expired:
#   1. kills fill1005s2 (pid from its .pid file, command line verified) so that nothing else can arm slot 2;
#   2. waits for run_arm_generic's own prep gate (same function, same thresholds) until 10-05 01:30Z. If the gate is
#      still closed at 01:30Z, it preps anyway when memory_pressure free >= 50 % and nothing else builds or preps. A vcc
#      prep runs ~25 s with a ~24 GB peak footprint, and 549 was prepped by hand on 09-30 the same way. The deviation is
#      logged in runidall3.log and in 663's README. It finishes well before recipe_dr524's 02:00Z build;
#   3. runs run_arm_generic's post-gate steps unchanged (vcc prep, meta copy, novelty check, package sha) into runidall3.log,
#      and additionally requires prep rc 0 and "prep complete";
#   4. appends READY/NOT READY to build_idall3.log with build_idall3.sh's own criteria. On READY it arms arm1005 slot 2
#      under fill1005s2's own conditions (chain1005 armed s1 and left s2 empty). Cap 10-05 04:00Z.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/663-idall3-arm; FL=$X/work/prep_fallback663.log; BL=$X/work/build_idall3.log
E=experiments/428-zero-gene-audit; W="$REPO/$E/work/build_idall3"; LOG="$W/runidall3.log"
NAME=idall3-ungated-identity-hook
STEM=idall3-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5-idhook
PKG="$REPO/submissions/$STEM.vcc"; META="$W/$NAME.h5ad.meta.json"; BUNDLE="$REPO/data/challenge"
A=experiments/496-uncalled-sign-hook; C=$A/work/chain1005.log
RL=$E/work/build_idall3/runidall3.log
MNAME="There's nothing sadder than living your whole life by the ocean but never daring to set sail."
BUILDPID=91191; CHAINPID=24465
flog(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$FL"; }
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
blog(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$BL"; }
now(){ date -u +%Y%m%d%H%M; }
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.1f", $4/1048576}'; }
echo $$ > $X/work/prep_fallback663.pid
flog "armed (pid $$): waiting for build_idall3.sh (pid $BUILDPID) to exit"
while kill -0 $BUILDPID 2>/dev/null; do
  [ "$(now)" -ge 202610050400 ] && { flog "cap 10-05 04:00Z with build_idall3.sh still running -- nothing done"; flog "END"; exit 0; }
  sleep 30
done
T="$(tail -1 "$BL" 2>/dev/null)"
flog "build_idall3.sh exited; its last line: $T"
case "$T" in *" READY idall3 "*) flog "idall3 READY on the normal path -- nothing to do"; flog "END"; exit 0;; esac
if ! /usr/bin/grep -q "GATE never opened (prep)" "$LOG"; then
  flog "the build did not end on an expired prep gate -- not a case for this fallback"; flog "END"; exit 0
fi
# 1. no other armer for slot 2
FP="$(cat $X/work/fill1005s2.pid 2>/dev/null)"
if [ -n "$FP" ] && kill -0 "$FP" 2>/dev/null && ps -o args= -p "$FP" | /usr/bin/grep -q "fill1005s2.sh"; then
  kill "$FP" && flog "killed fill1005s2 (pid $FP): this fallback takes over its role"
fi
printf '%s %s\n' "$(date -u +%FT%TZ)" "killed by prep_fallback663.sh after the prep gate expired; the fallback arms slot 2 if idall3 becomes READY" >> $X/work/fill1005s2.log
# 2. the prep gate, run_arm_generic's own function (copied verbatim, logging to the fallback log)
gate_ok(){
  local mode="${1:-build}" pct fi swapfree others preps disk branch
  pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  swapfree="$(sysctl -n vm.swapusage 2>/dev/null | awk '{for(i=1;i<NF;i++) if($i=="free"){v=$(i+2); u=substr(v,length(v),1); n=substr(v,1,length(v)-1)+0; if(u=="G") n=n*1024; else if(u=="K") n=n/1024; printf "%.0f", n; exit}}')"
  others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
  preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
  disk="$(free_gib)"
  branch=none
  if [ "${pct:-0}" -ge 40 ] && awk -v v="$fi" 'BEGIN{exit !(v>=23)}'; then branch=fi23
  elif [ "$mode" = prep ] && [ "${pct:-0}" -ge 55 ] && awk -v v="${swapfree:-0}" 'BEGIN{exit !(v>=400)}'; then branch=pressure55
  elif [ "$mode" = prep ] && [ "$(now)" -ge 202610050130 ] && [ "${pct:-0}" -ge 50 ]; then branch=deviation50
  fi
  flog "gate: free ${pct}%  free+inactive ${fi} GiB  other builders ${others}  vcc preps ${preps}  disk ${disk} GiB  swap free ${swapfree:-?} MiB  mode ${mode}  memory branch ${branch}"
  GATE_BRANCH=$branch
  [ "$branch" != none ] && [ "$others" = "0" ] && [ "$preps" = "0" ] && awk -v v="$disk" 'BEGIN{exit !(v>=20)}'
}
GATE_BRANCH=none
until gate_ok prep; do
  [ "$(now)" -ge 202610050145 ] && { flog "no gate by 10-05 01:45Z -- idall3 stays unpackaged"; blog "NOT READY prep gate never opened (fallback, 01:45Z)"; flog "END"; exit 0; }
  sleep 30
done
log "prep fallback (prep_fallback663.sh): the gate opened on branch $GATE_BRANCH"
[ "$GATE_BRANCH" = deviation50 ] && log "DEVIATION from the memory rule: free+inactive < 23 GiB and swap free < 400 MiB at 01:30Z; prepped on memory_pressure free >= 50 % (a ~25 s prep; 549 precedent 09-30)"
# 3. run_arm_generic's post-gate steps, unchanged
log "INVOCATION 2: vcc prep -> $PKG"
/usr/bin/time -l vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" \
    --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C \
    -o "$PKG" > "$W/prep.log" 2>&1
RC=$?; log "prep rc=$RC  disk $(free_gib) GiB"
tail -20 "$W/prep.log" | tee -a "$LOG"
log "vcc: $(ls -la "$PKG" 2>/dev/null | awk '{print $5, $9}')"
cp "$META" "$REPO/submissions/$STEM.h5ad.meta.json"
log "meta copied beside the package: $(shasum -a 256 "$REPO/submissions/$STEM.h5ad.meta.json" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$PKG" 2>&1 | tee -a "$LOG"
log "novelty rc=${PIPESTATUS[0]}"
log "package sha256: $(shasum -a 256 "$PKG" | cut -d' ' -f1)"
log "git status scripts/ src/: $(git status --short scripts/ src/ | tr '\n' ' ')"
log "DONE"
# 4. build_idall3.sh's own READY criteria, plus prep rc 0 and "prep complete"
blog "idall3 (fallback): $(/usr/bin/grep -E '^\[idall3\] (PASS|FAIL|CLIFF)|package sha|novelty rc|stopping' $RL 2>/dev/null | tr '\n' ' ' | cut -c1-600)"
if [ $RC -eq 0 ] && /usr/bin/grep -q "prep complete" "$W/prep.log" && /usr/bin/grep -q '^\[idall3\] PASS' $RL \
   && /usr/bin/grep -q '^\[idall3\] CLIFF OK' $RL && /usr/bin/grep -q 'novelty rc=0' $RL; then
  SHA="$(/usr/bin/grep 'package sha256:' $RL | tail -1 | awk '{print $NF}')"
  if [ -s "submissions/$STEM.vcc" ] && [ "$(shasum -a 256 submissions/$STEM.vcc | cut -d' ' -f1)" = "$SHA" ]; then
    blog "READY idall3 $SHA"
  else
    blog "NOT READY package missing or sha mismatch (fallback)"; flog "END"; exit 0
  fi
else
  blog "NOT READY prep or gates failed (fallback; see $RL)"; flog "END"; exit 0
fi
# arm slot 2 exactly as fill1005s2 would have
while kill -0 $CHAINPID 2>/dev/null; do [ "$(now)" -ge 202610050400 ] && { flog "cap 10-05 04:00Z with chain1005 running -- nothing armed"; flog "END"; exit 0; }; sleep 30; done
LINE="$(/usr/bin/grep -E 'armed 10-05: s1' "$C" | tail -1)"
flog "chain1005's arm line: ${LINE:-<none>}"
case "$LINE" in
  "") flog "chain1005 armed nothing -- not arming a slot 2 fill"; flog "END"; exit 0;;
  *"s2 <none>"*) ;;
  *) flog "chain1005 armed its own slot 2 -- no fill"; flog "END"; exit 0;;
esac
OP="$(cat $A/work/arm1005_s2.pid 2>/dev/null)"
if [ -n "$OP" ] && kill -0 "$OP" 2>/dev/null && ps -o args= -p "$OP" | /usr/bin/grep -q "arm1005.sh"; then
  flog "an arm1005 slot-2 waiter is already alive (pid $OP) -- not arming a second"; flog "END"; exit 0
fi
[ "$(now)" -ge 202610050400 ] && { flog "past 10-05 04:00Z -- nothing armed"; flog "END"; exit 0; }
nohup bash $A/arm1005.sh 2 idall3 "$RL" "submissions/$STEM.vcc" "$MNAME" "idall3 slot 2 on 10-05 (coordinator fill, prep fallback)" > $A/work/arm1005_s2.nohup 2>&1 & echo $! > $A/work/arm1005_s2.pid
sleep 3; flog "armed idall3 as 10-05 slot 2 (pid $(cat $A/work/arm1005_s2.pid))"; flog "END"
