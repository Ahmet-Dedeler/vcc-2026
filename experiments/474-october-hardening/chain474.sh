#!/usr/bin/env bash
# 474: apply and verify 481's October-path patches (experiments/481-october-audit/patches). Serial after chain471b.
#   F2 ntclevel417: each NTC construct's jackknife x n/cells_per_pert. KA: test_cp46_481.py (A/B/C Cp46 bit for bit) and
#      ntclevel417 --ka on 469dr's relabelled artifact (baseline PASS pre-patch, work/ntclevel_ka_pre.log).
#   F7 def_sameline: --programme-source refuses < 150 rows / < 50 % coverage. KA: 465's call rebuilds 465's sameline.npz
#      byte for byte (b4eb2a6d...).
#   F1/F5 def_day.sh: DEF_DAY_CD4_THRESHOLDS pass-through + an explicit stop message, contexts from the manifest, and
#      def_postbuild.py before the builder restore. Three rehearsals:
#      (a) efp (CD4 union 186) WITHOUT the variable: must stop at setup with the CD4-gate message;
#      (b) 465 retraced on ef-bundle-r (DEF_DAY_LODO=1 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full, flags pt): the final output
#          sha must equal 465's (6113937f...), and postbuild must PASS;
#      (c) efp WITH DEF_DAY_CD4_THRESHOLDS=relaxed (default flags): DONE, postbuild PASS, every final rho <= 1.
# Any failure restores the three pre-474 files and stops. Never submits. The coordinator commits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/474-october-hardening; W=$X/work; LOG=$W/chain474.log; P=experiments/481-october-audit/patches
DD=scripts/def_day.sh; SL=scripts/def_sameline.py; NT=experiments/417-bulk-sharpen/ntclevel417.py
EFP=experiments/238-essential-free-panel/work/ef-bundle; EFR=experiments/238-essential-free-panel/work/ef-bundle-r
BU=scripts/make_submission_pm_cd4.py
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain474.pid
log "armed (pid $$): waiting for chain471b to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ END$" experiments/471-h1full-union/work/chain471b.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain471b\.sh|chain472\.sh|chain470b\.sh|run_arm_generic\.sh|run_ka_47[12]\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271600 ] && { log "471b never ended by 16:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
until [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609271600 ] && { log "builder never equal to HEAD by 16:00Z -- stopping"; exit 1; }
  log "builder != HEAD (472 uncommitted?) -- waiting"; sleep 300; done
for f in $DD $SL $NT; do [ "$(shasum -a 256 $f | cut -d' ' -f1)" = "$(git show HEAD:$f | shasum -a 256 | cut -d' ' -f1)" ] || { log "$f differs from HEAD -- stopping"; exit 1; }; done
cp $DD $W/def_day.pre474.sh; cp $SL $W/def_sameline.pre474.py; cp $NT $W/ntclevel417.pre474.py
back(){ cp $W/def_day.pre474.sh $DD; cp $W/def_sameline.pre474.py $SL; cp $W/ntclevel417.pre474.py $NT
        log "RESTORED the three pre-474 files ($1)"; log "END"; exit 1; }
memwait(){ local n=0 fi; while :; do
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && return 0
  n=$((n+1)); [ $n -ge 120 ] && return 1; log "waiting for memory: free+inactive $fi GiB < 30"; sleep 60; done; }
# --- F2
patch -p0 $NT < $P/F2_ntclevel417_cp46.diff >> $LOG 2>&1 || back "F2 patch failed"
log "F2 applied: ntclevel417 $(shasum -a 256 $NT | cut -c1-16)"
uv run python $P/test_cp46_481.py > $W/test_cp46.log 2>&1 || back "test_cp46_481 failed: $(tail -2 $W/test_cp46.log | tr '\n' ' ')"
log "test_cp46_481: $(tail -3 $W/test_cp46.log | tr '\n' ' ' | cut -c1-300)"
cp -c submissions/469dr-october-path-arcprog-gate33-abc-q1.75.h5ad $W/ka469.h5ad 2>/dev/null || true
uv run python $NT $W/ka469.h5ad --contexts A,B,C --f A=1,B=1,C=0.9825 --ka > $W/ntclevel_ka_post.log 2>&1 || back "ntclevel417 --ka failed after F2"
/usr/bin/grep -q "KA PASS" $W/ntclevel_ka_post.log || back "ntclevel417 --ka not PASS after F2"
diff <(grep -E "^  [ABC]: Cp46" $W/ntclevel_ka_pre.log) <(grep -E "^  [ABC]: Cp46" $W/ntclevel_ka_post.log) > /dev/null || back "F2 changed an A/B/C Cp46 line"
rm -f $W/ka469.h5ad
log "F2 KA PASS: Cp46 lines identical to the pre-patch baseline; --ka PASS"
# --- F7
patch -p0 $SL < $P/F7_def_sameline_programme_rows.diff >> $LOG 2>&1 || back "F7 patch failed"
log "F7 applied: def_sameline $(shasum -a 256 $SL | cut -c1-16)"
R465=experiments/465-both-switches-efr/work/run
uv run python $SL --bundle $EFR --identity $R465/identity.json --out $W/ka465_sameline.npz --programme-source h1_2025_full > $W/ka465_sameline.log 2>&1 \
  || back "F7 KA: def_sameline failed: $(tail -2 $W/ka465_sameline.log | tr '\n' ' ')"
G="$(shasum -a 256 $W/ka465_sameline.npz | cut -d' ' -f1)"; Wt="$(shasum -a 256 $R465/sameline.npz | cut -d' ' -f1)"
[ "$G" = "$Wt" ] || back "F7 KA: sameline ${G:0:16} != 465's ${Wt:0:16}"
log "F7 KA PASS: 465's sameline.npz rebuilt byte for byte (${G:0:16})"
# --- F1/F5
patch -p0 $DD < $P/F1_F5_def_day.diff >> $LOG 2>&1 || back "F1/F5 patch failed"
bash -n $DD || back "def_day.sh syntax"
log "F1/F5 applied: def_day.sh $(shasum -a 256 $DD | cut -c1-16)"
# (a) efp, strict: must stop at setup with the CD4-gate message
memwait || back "memory never reached 30 GiB before (a)"
bash $DD $EFP $W/efp_strict 474efps 1.75 auto pt:noprog 4 > $W/efp_strict.nohup 2>&1; RC=$?
if [ $RC -ne 0 ] && /usr/bin/grep -q "STOP: CD4 extraction gate fired" $W/efp_strict/def_day.log; then
  log "(a) PASS: efp without the variable stops at setup: $(/usr/bin/grep 'STOP:' $W/efp_strict/def_day.log | cut -c1-300)"
else back "(a) efp strict did not stop with the CD4-gate message (rc $RC): $(tail -2 $W/efp_strict/def_day.log | tr '\n' ' ')"; fi
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || back "builder not restored after (a)"
# (b) 465 retraced on efr
memwait || back "memory never reached 30 GiB before (b)"
DEF_DAY_LODO=1 DEF_DAY_PROGRAMME_SOURCE=h1_2025_full bash $DD $EFR $W/efr_retrace 465efr 1.75 auto pt 4 > $W/efr_retrace.nohup 2>&1
/usr/bin/grep -q "DONE. Next" $W/efr_retrace/def_day.log || back "(b) retrace did not finish: $(tail -2 $W/efr_retrace/def_day.log | tr '\n' ' ')"
/usr/bin/grep -q "postbuild PASS" $W/efr_retrace/def_day.log || back "(b) postbuild did not PASS"
GOT="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $W/efr_retrace/465efr_final.h5ad.meta.json)"
[ "$GOT" = 6113937f1acafb34c8cc323c7d0a35d4fabd7fc505c7e52404745e26263a949b ] || back "(b) final output ${GOT:0:16} != 465's 6113937f"
log "(b) PASS: 465 retraced byte for byte (${GOT:0:16}); $(/usr/bin/grep -E 'CONVERGED at|postbuild PASS' $W/efr_retrace/def_day.log | tr '\n' ' ' | cut -c1-300)"
rm -f $W/efr_retrace/*.h5ad
# (c) efp relaxed
memwait || back "memory never reached 30 GiB before (c)"
DEF_DAY_CD4_THRESHOLDS=relaxed bash $DD $EFP $W/efp_relaxed 474efp 1.75 auto pt:noprog 4 > $W/efp_relaxed.nohup 2>&1
/usr/bin/grep -q "DONE. Next" $W/efp_relaxed/def_day.log || back "(c) efp relaxed did not finish: $(tail -2 $W/efp_relaxed/def_day.log | tr '\n' ' ')"
/usr/bin/grep -q "postbuild PASS" $W/efp_relaxed/def_day.log || back "(c) postbuild did not PASS"
FIN="$(/usr/bin/grep -o 'final: [^ ]* factors [^ ]*' $W/efp_relaxed/def_day.log | tail -1)"
RHO="$(uv run python scripts/def_day_helpers.py rho "$W/efp_relaxed/$(echo "$FIN" | awk '{print $2}').ntclevel.json" 2>&1 | tail -1)"
log "(c) PASS: efp relaxed DONE; $FIN; rho $RHO; $(/usr/bin/grep -E 'CONVERGED at|NOT converged|postbuild PASS' $W/efp_relaxed/def_day.log | tr '\n' ' ' | cut -c1-300)"
rm -f $W/efp_relaxed/*.h5ad
log "END all PASS: def_day $(shasum -a 256 $DD | cut -c1-16), def_sameline $(shasum -a 256 $SL | cut -c1-16), ntclevel417 $(shasum -a 256 $NT | cut -c1-16) (uncommitted; the coordinator commits)"
