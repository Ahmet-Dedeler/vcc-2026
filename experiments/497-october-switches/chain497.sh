#!/usr/bin/env bash
# 497: the October pass-throughs DEF_DAY_LODO_HI_PCT and DEF_DAY_SHARPEN_PT (patch_defday497.py) on top of 476, then:
#   tests: bash -n; test_helpers497.py (9 cases); 485's repro_fixes against the 497 copies (repro/, rc 0 at prep time);
#   (a) inert KA: efp relaxed with the 497 script and neither knob == 474/476 (c)'s ac49d319 byte for byte;
#   (b) 497dr: 460dr's recipe (180's def-bundle, 453's operator cloned, DEF_DAY_LODO=1, flags pt:noprog) + the percentile
#       cap (DEF_DAY_LODO_HI_PCT=96.64697522769474 == cap 6 on A/B/C) + S_P2 (DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75);
#       relabelled to A/B/C, prepped, novelty-checked; gate497.py vs 460dr.
# Serial after chain476 (process exit) and after the coordinator commits 476 (def_day.sh and the helpers equal 476's
# versions at HEAD). Any failure in apply/tests/(a) restores the two pre-497 files. Never submits; the coordinator commits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/497-october-switches; W=$X/work; LOG=$W/chain497.log; R453=experiments/453-abc-dress-rehearsal/work/run
DD=scripts/def_day.sh; HP=scripts/def_day_helpers.py; BU=scripts/make_submission_pm_cd4.py
EFP=experiments/238-essential-free-panel/work/ef-bundle; BUNDLE=experiments/180-def-readiness/work/def-bundle
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain497.pid
log "armed (pid $$): waiting for chain476's process to exit and 476 committed"
P476="$(cat experiments/476-def-day-merge/work/chain476.pid 2>/dev/null)"
until ! ps -p "${P476:-1}" -o args= 2>/dev/null | /usr/bin/grep -q chain476 \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain47[69]b?\.sh|chain495\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4/' | /usr/bin/grep -q . \
      && [ "$(shasum -a 256 $DD | cut -d' ' -f1)" = "$(shasum -a 256 experiments/476-def-day-merge/work/def_day.resolved.sh | cut -d' ' -f1)" ] \
      && [ "$(shasum -a 256 $DD | cut -d' ' -f1)" = "$(git show HEAD:$DD | shasum -a 256 | cut -d' ' -f1)" ] \
      && [ "$(shasum -a 256 $HP | cut -d' ' -f1)" = "$(git show HEAD:$HP | shasum -a 256 | cut -d' ' -f1)" ] \
      && [ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609301800 ] && { log "476 not landed by 18:00Z 09-30 -- stopping"; log "END"; exit 1; }; sleep 60; done
log "chain476 exited and 476 is at HEAD; its last lines: $(tail -2 experiments/476-def-day-merge/work/chain476.log | tr '\n' ' ' | cut -c1-300)"
cp $DD $W/pre497.def_day.sh; cp $HP $W/pre497.def_day_helpers.py
back(){ cp $W/pre497.def_day.sh $DD; cp $W/pre497.def_day_helpers.py $HP; log "RESTORED the two pre-497 files ($1)"; log "END"; exit 1; }
# --- apply + tests
uv run python $X/patch_defday497.py $DD $HP >> $LOG 2>&1 || back "patch failed"
cmp -s $DD $X/repro/def_day.sh.new && cmp -s $HP $X/repro/def_day_helpers.py.new || back "applied files differ from repro/*.new"
bash -n $DD || back "def_day.sh syntax"
uv run python $X/test_helpers497.py $HP > $W/test_helpers497.log 2>&1 || back "test_helpers497 failed"
(cd $X/repro && uv run python repro_fixes.py > last-fixes-chain.log 2>&1) || back "repro_fixes on the 497 files failed"
log "applied: def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16); test_helpers497 9/9; repro_fixes PASS"
memwait(){ local n=0 fi; while :; do
  fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
  awk -v v="$fi" 'BEGIN{exit !(v>=30)}' && return 0
  n=$((n+1)); [ $n -ge 120 ] && return 1; log "waiting for memory: free+inactive $fi GiB < 30"; sleep 60; done; }
# --- (a) inert KA
memwait || back "memory never reached 30 GiB before (a)"
rm -rf $W/efp_inert
DEF_DAY_CD4_THRESHOLDS=relaxed bash $DD $EFP $W/efp_inert 474efp 1.75 auto pt:noprog 4 > $W/efp_inert.nohup 2>&1
/usr/bin/grep -q "DONE. Next" $W/efp_inert/def_day.log || back "(a) did not finish: $(tail -2 $W/efp_inert/def_day.log | tr '\n' ' ')"
GOTA="$(uv run python -c 'import json,sys; print(json.load(open(sys.argv[1]))["output"]["sha256"])' $W/efp_inert/474efp_final.h5ad.meta.json)"
[ "$GOTA" = ac49d319f5982507003b1f2070660d2d03ee351fd13c9fdd30df8f5fb08ec0b3 ] || back "(a) final ${GOTA:0:16} != ac49d319"
rm -f $W/efp_inert/*.h5ad
[ "$(shasum -a 256 $BU | cut -d' ' -f1)" = "$(git show HEAD:$BU | shasum -a 256 | cut -d' ' -f1)" ] || back "builder not restored after (a)"
log "(a) PASS: the 497 script with neither knob retraces efp relaxed byte for byte (${GOTA:0:16})"
# --- (b) 497dr
[ "$(shasum -a 256 "$R453/whitener_union_pruned_def.npz" | cut -d' ' -f1)" = 66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe ] \
  || { log "453's day operator is not 66610be7 -- stopping (497 files stay applied: (a) passed)"; log "END"; exit 1; }
memwait || { log "memory never reached 30 GiB before (b)"; log "END"; exit 1; }
mkdir -p "$W/run"; cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "(b) 497dr: DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474 DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75 on $BUNDLE"
DEF_DAY_LODO=1 DEF_DAY_LODO_HI_PCT=96.64697522769474 DEF_DAY_SHARPEN_PT=stercile:2.5,2.0,1.75 \
  bash $DD "$BUNDLE" "$W/run" 497dr 1.75 auto pt:noprog 4 > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || { log "(b) def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ')"; log "END"; exit 1; }
log "(b) def_day: $(/usr/bin/grep -E "LODO key ON|per-target sharpening|CONVERGED at|NOT converged|final:|postbuild PASS" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-600)"
STEM="submissions/497dr-october-path-lodopct-spt-abc"
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/497dr_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || { log "relabel failed"; log "END"; exit 1; }
FREE="$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$NF); s+=$NF} END {printf "%d", s*16384/2^30}')"
[ "$FREE" -ge 23 ] || { log "free+inactive ${FREE} GiB < 23 -- not prepping; rerun prep by hand"; log "END"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || { log "prep dry run failed"; log "END"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || { log "prep failed"; log "END"; exit 1; }
log "package $STEM.vcc sha256 $(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
uv run python $X/gate497.py >> "$LOG" 2>&1; log "gate497 rc=$?"
rm -f "$W/run/"*.h5ad
log "END def_day $(shasum -a 256 $DD | cut -c1-16), helpers $(shasum -a 256 $HP | cut -c1-16) (uncommitted; the coordinator commits with (a))"
