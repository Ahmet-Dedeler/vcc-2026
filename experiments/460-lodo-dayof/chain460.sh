#!/usr/bin/env bash
# 460: the October path WITH the LODO stack key (DEF_DAY_LODO=1) end to end on A/B/C, i.e. 453 + the LODO key. It is the
# October-path twin of 446bj2, board-scorable against 453. `scripts/def_day.sh` (HEAD) runs on 180's `def-bundle` (the
# A/B/C controls and panel, relabelled A->D, B->E, C->F) with 453's re-saved day operator cloned in (66610be7; a fresh
# refit on this panel equals A/B/C's operator and trips the driver's panel-changed guard, 453 attempt 1). Then
# relabel453.py maps the final artifact back to A/B/C, then vcc prep (dry run, real) and the novelty check. Never submits.
#   bash experiments/460-lodo-dayof/chain460.sh
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/460-lodo-dayof; W=$X/work; LOG=$W/chain460.log; R453=experiments/453-abc-dress-rehearsal/work/run
BUNDLE=experiments/180-def-readiness/work/def-bundle; N=460dr; Q=1.75; FLAGS=pt:noprog
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/chain460.pid"
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
BH="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
[ "$BF" = "$BH" ] || { log "builder ${BF:0:16} differs from HEAD's ${BH:0:16} -- stopping"; exit 1; }
DD="$(shasum -a 256 scripts/def_day.sh | cut -d' ' -f1)"; DH="$(git show HEAD:scripts/def_day.sh | shasum -a 256 | cut -d' ' -f1)"
[ "$DD" = "$DH" ] || { log "def_day.sh differs from HEAD -- stopping"; exit 1; }
[ "$(shasum -a 256 "$R453/whitener_union_pruned_def.npz" | cut -d' ' -f1)" = 66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe ] \
  || { log "453's day operator is not 66610be7 -- stopping"; exit 1; }
mkdir -p "$W/run"
cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "builder ${BF:0:16} = HEAD; def_day.sh ${DD:0:16} = HEAD; operator 66610be7 cloned; starting DEF_DAY_LODO=1 def_day.sh on $BUNDLE"
DEF_DAY_LODO=1 bash scripts/def_day.sh "$BUNDLE" "$W/run" "$N" "$Q" auto "$FLAGS" 4 > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || { log "def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ')"; exit 1; }
log "def_day: $(/usr/bin/grep -E "LODO key ON|CONVERGED at|NOT converged|final:" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-400)"
STEM="submissions/460dr-october-path-lodo-abc-q${Q}"
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
