#!/usr/bin/env bash
# 453: the October path end to end on A/B/C, board-scorable. `scripts/def_day.sh` (unmodified, its default D/E/F path)
# on 180's `def-bundle` (= the A/B/C controls and panel, relabelled A->D, B->E, C->F), then relabel453.py maps the final
# artifact back to A/B/C, then vcc prep (dry run, real) and the novelty check. It never submits.
# Waits for chain0926 to finish (it builds with the production builder) and for work/args.env (Q, FLAGS), which the
# coordinator writes after the 09-26 reads.   bash experiments/453-abc-dress-rehearsal/chain453.sh
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/453-abc-dress-rehearsal; W=$X/work; LOG=$W/chain453.log; C926=experiments/446-q175-level/work/chain0926.log
BUNDLE=experiments/180-def-readiness/work/def-bundle; N=453dr
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/chain453.pid"
log "armed (pid $$): waiting for chain0926 to end and for $W/args.env"
until /usr/bin/grep -qE " END$|stopping" "$C926" 2>/dev/null && [ -s "$W/args.env" ]; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609262000 ] && { log "inputs not ready by 20:00Z 09-26 -- stopping"; exit 1; }
  sleep 60
done
# shellcheck disable=SC1090
source "$W/args.env"
: "${Q:?args.env must set Q}" "${FLAGS:?args.env must set FLAGS}"
log "inputs: chain0926 $(tail -1 "$C926" | cut -c1-80); Q $Q FLAGS $FLAGS"
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
BH="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
[ "$BF" = "$BH" ] || { log "builder ${BF:0:16} differs from HEAD's ${BH:0:16} -- stopping"; exit 1; }
DD="$(shasum -a 256 scripts/def_day.sh | cut -d' ' -f1)"; DH="$(git show HEAD:scripts/def_day.sh | shasum -a 256 | cut -d' ' -f1)"
[ "$DD" = "$DH" ] || { log "def_day.sh differs from HEAD -- stopping"; exit 1; }
log "builder ${BF:0:16} = HEAD; def_day.sh ${DD:0:16} = HEAD; starting def_day.sh on $BUNDLE"
bash scripts/def_day.sh "$BUNDLE" "$W/run" "$N" "$Q" auto "$FLAGS" 4 > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || { log "def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ')"; exit 1; }
log "def_day: $(/usr/bin/grep -E "CONVERGED at|NOT converged|final:" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-300)"
STEM="submissions/453dr-october-path-abc-q${Q}"
uv run python "$X/relabel453.py" "$W/run/${N}_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || { log "relabel failed"; exit 1; }
FREE="$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$NF); s+=$NF} END {printf "%d", s*16384/2^30}')"
[ "$FREE" -ge 23 ] || { log "free+inactive ${FREE} GiB < 23 -- not prepping; rerun prep by hand"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || { log "prep dry run failed"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || { log "prep failed"; exit 1; }
log "package $STEM.vcc sha256 $(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
log "END"
