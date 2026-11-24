#!/usr/bin/env bash
# 457: after 452c ends, build the CD4 x1.5 twin of 453 through def_day_cd4x15.sh on 180's relabelled bundle, relabel it
# to A/B/C, prep it, check novelty. The work dir is pre-seeded with 453's re-saved operator (66610be7, arrays = e2f65ecf) so
# the driver skips the refit and its "panel did not change" guard, exactly as 453's attempt 3 did. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/457-cd4-half-def; W=$X/work; LOG=$W/chain457.log; BUNDLE=experiments/180-def-readiness/work/def-bundle; N=457dr
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/chain457.pid
log "armed (pid $$): waiting for 452c's chain to END"
until /usr/bin/grep -qE " END$|stopping" experiments/452-def-day-rehearsal/work/chain452c.log 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain452c\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609262200 ] && { log "452c never ended by 22:00Z -- stopping"; exit 1; }; sleep 60; done
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"; BH="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
[ "$BF" = "$BH" ] || { log "builder ${BF:0:16} != HEAD ${BH:0:16} -- stopping"; exit 1; }
mkdir -p $W/run; S3=experiments/453-abc-dress-rehearsal/work/run
cp -c $S3/whitener_union_pruned_def.npz $W/run/ && cp $S3/whitener_union_pruned_def.json $W/run/ || { log "seed copy failed"; exit 1; }
log "seeded operator $(shasum -a 256 $W/run/whitener_union_pruned_def.npz | cut -c1-16); starting def_day_cd4x15.sh"
bash $X/def_day_cd4x15.sh "$BUNDLE" "$W/run" "$N" 1.75 auto pt:noprog 4 > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || { log "def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ')"; exit 1; }
log "dump CD4 check: $(/usr/bin/grep -o '\[cd4\] sizing=[a-z]* g=[0-9.]*' $W/run/driver-dump.log | head -1)"
log "def_day: $(/usr/bin/grep -E "CONVERGED at|NOT converged|final:" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-300)"
STEM="submissions/457dr-october-path-abc-cd4x15-q1.75"
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/${N}_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || { log "relabel failed"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 || { log "prep dry run failed"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 || { log "prep failed"; exit 1; }
log "package $STEM.vcc sha256 $(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
rm -f "$W/run/${N}_final.h5ad" "$W/run/${N}-lvl"*.h5ad
log "END"
