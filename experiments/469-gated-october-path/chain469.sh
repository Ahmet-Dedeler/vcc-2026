#!/usr/bin/env bash
# 469: the October path with Arc's denoised H1 programme AND the strength gate, end to end on A/B/C:
# `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=q33 def_day.sh` with sameline_flags `pt` on 180's relabelled
# bundle (q 1.75, auto). It is 461dr with the gate on: the October-path twin of 458g, and the end-to-end rehearsal of the
# DEF_DAY_PRIOR_GATE switch. 453's re-saved day operator (66610be7) is cloned in, as for 460/461. Waits for chain468 to
# END (strictly serial) and for the 468 builder (693fde15) to be committed; def_day.sh must be the gate-switch bytes
# (5f9bfa26). Relabel -> prep -> novelty. Never submits.
#   bash experiments/469-gated-october-path/chain469.sh
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/469-gated-october-path; W=$X/work; LOG=$W/chain469.log; R453=experiments/453-abc-dress-rehearsal/work/run
BUNDLE=experiments/180-def-readiness/work/def-bundle; N=469dr; Q=1.75; FLAGS=pt; GATE=q33
C468=experiments/468-prior-gate/work/chain468.log
B468=693fde151be7c9994b3aa643480d8dbf759242b2f5a24ed044e688e7e8f0a405
DD468=5f9bfa26823f355856a720756864b0a3c6e50d273102c3b6db91f08029131491
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/chain469.pid"
log "armed (pid $$): waiting for chain468 to END and the 468 builder to be committed"
until /usr/bin/grep -qE "^[0-9TZ:-]+ END( builder .*)?$" "$C468" 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain468\.sh|run_arm_generic\.sh|run_ka_468\.sh|def_day\.sh/' | /usr/bin/grep -q . \
      && [ "$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)" = "$B468" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609270030 ] && { log "chain468 END + committed builder not seen by 00:30Z 09-27 -- stopping"; exit 1; }
  sleep 60; done
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
[ "$BF" = "$B468" ] || { log "builder ${BF:0:16} is not 693fde15 (pinned, or restored after a failed KA) -- stopping"; exit 1; }
DD="$(shasum -a 256 scripts/def_day.sh | cut -d' ' -f1)"
[ "$DD" = "$DD468" ] || { log "def_day.sh ${DD:0:16} is not the gate-switch bytes 5f9bfa26 -- stopping"; exit 1; }
/usr/bin/grep -q "DEF_DAY_PRIOR_GATE" scripts/def_day.sh || { log "def_day.sh has no DEF_DAY_PRIOR_GATE -- stopping"; exit 1; }
[ "$(shasum -a 256 data/donors/h1_2025_full/delta.npy | cut -d' ' -f1)" = 26662c527164ddd434ae19a13076e915d1f82de9683c75e27e0af1f5b28591c0 ] \
  || { log "h1_2025_full delta sha unexpected -- stopping"; exit 1; }
[ "$(shasum -a 256 "$R453/whitener_union_pruned_def.npz" | cut -d' ' -f1)" = 66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe ] \
  || { log "453's day operator is not 66610be7 -- stopping"; exit 1; }
cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "builder ${BF:0:16} = HEAD; def_day.sh ${DD:0:16}; starting DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=$GATE def_day.sh (flags $FLAGS)"
DEF_DAY_PROGRAMME_SOURCE=h1_2025_full DEF_DAY_PRIOR_GATE=$GATE bash scripts/def_day.sh "$BUNDLE" "$W/run" "$N" "$Q" auto "$FLAGS" 4 > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || { log "def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ')"; exit 1; }
log "def_day: $(/usr/bin/grep -E "programme-source|prior gate|CONVERGED at|NOT converged|final:" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-500)"
log "prior gate in the builds: $(/usr/bin/grep -h '^\[prior_gate\]' "$W"/run/logs/build*.log 2>/dev/null | sort -u | tr '\n' ' ' | cut -c1-500)"
STEM="submissions/469dr-october-path-arcprog-gate33-abc-q${Q}"
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/${N}_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || { log "relabel failed"; exit 1; }
FREE="$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$NF); s+=$NF} END {printf "%d", s*16384/2^30}')"
[ "$FREE" -ge 23 ] || { log "free+inactive ${FREE} GiB < 23 -- not prepping; rerun prep by hand"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || { log "prep dry run failed"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || { log "prep failed"; exit 1; }
log "package $STEM.vcc sha256 $(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
log "END"
