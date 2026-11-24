#!/usr/bin/env bash
# 2026-09-28: finish chain497's (b) by hand after def_day's postbuild was fixed and re-run (see run/def_day.log). Same steps
# and log lines as chain497.sh lines 60-75: relabel, vcc prep (dry run, then real), novelty, gate497. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/497-october-switches; W=$X/work; LOG=$W/chain497.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
STEM="submissions/497dr-october-path-lodopct-spt-abc"
uv run python experiments/453-abc-dress-rehearsal/relabel453.py "$W/run/497dr_final.h5ad" "$STEM.h5ad" >> "$LOG" 2>&1 || { log "relabel failed"; log "END"; exit 1; }
until [ "$(vm_stat | awk '/Pages free|Pages inactive/ {gsub("\\.","",$NF); s+=$NF} END {printf "%d", s*16384/2^30}')" -ge 23 ]; do sleep 60; done
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" --dry-run >> "$LOG" 2>&1 \
  || { log "prep dry run failed"; log "END"; exit 1; }
vcc prep "$STEM.h5ad" -g data/challenge/gene_names.csv --perts data/challenge/pert_counts.csv -o "$STEM.vcc" >> "$LOG" 2>&1 \
  || { log "prep failed"; log "END"; exit 1; }
log "package $STEM.vcc sha256 $(shasum -a 256 "$STEM.vcc" | cut -d' ' -f1)"
uv run python scripts/check_submission_novelty.py check "$STEM.vcc" >> "$LOG" 2>&1; log "novelty rc=$?"
uv run python $X/gate497.py --write >> "$LOG" 2>&1; log "gate497 rc=$?"
log "END finish497dr (manual continuation; run/*.h5ad kept until the gate log exists)"
