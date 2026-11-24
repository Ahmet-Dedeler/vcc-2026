#!/usr/bin/env bash
# 461: the October path with Arc's own (denoised) H1 programme as every context's block prior, end to end on A/B/C:
# `DEF_DAY_PROGRAMME_SOURCE=h1_2025_full def_day.sh` with sameline_flags `pt` on 180's relabelled bundle (q 1.75, auto),
# i.e. 453dr with the programme switched on. The October-path twin of 458d (vs 446b). 453's re-saved day operator
# (66610be7) is cloned in, as for 460. Waits for 458d's build (strictly serial) and fails closed unless def_day.sh at HEAD
# carries the DEF_DAY_PROGRAMME_SOURCE pass-through. Relabel -> prep -> novelty. Never submits.
#   bash experiments/461-arcprog-october-path/chain461.sh
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/461-arcprog-october-path; W=$X/work; LOG=$W/chain461.log; R453=experiments/453-abc-dress-rehearsal/work/run
BUNDLE=experiments/180-def-readiness/work/def-bundle; N=461dr; Q=1.75; FLAGS=pt
C458D=experiments/458-arc-lab-programme/work/chain458d.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > "$W/chain461.pid"
log "armed (pid $$): waiting for 458d's chain to END"
until /usr/bin/grep -qE "^[0-9TZ:-]+ (END|458d never.*|builder .* stopping|artifact sha .* stopping)$" "$C458D" 2>/dev/null \
      && ! ps -axo comm,args | awk '$1 ~ /^bash$/ && /chain458d\.sh|run_arm_generic\.sh|chain460\.sh|def_day\.sh/' | /usr/bin/grep -q .; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202609270000 ] && { log "458d never ended by 00:00Z 09-27 -- stopping"; exit 1; }; sleep 60; done
BF="$(shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1)"
BH="$(git show HEAD:scripts/make_submission_pm_cd4.py | shasum -a 256 | cut -d' ' -f1)"
[ "$BF" = "$BH" ] || { log "builder ${BF:0:16} differs from HEAD's ${BH:0:16} -- stopping"; exit 1; }
DD="$(shasum -a 256 scripts/def_day.sh | cut -d' ' -f1)"; DH="$(git show HEAD:scripts/def_day.sh | shasum -a 256 | cut -d' ' -f1)"
[ "$DD" = "$DH" ] || { log "def_day.sh differs from HEAD -- stopping"; exit 1; }
/usr/bin/grep -q "DEF_DAY_PROGRAMME_SOURCE" scripts/def_day.sh || { log "def_day.sh has no DEF_DAY_PROGRAMME_SOURCE pass-through -- stopping"; exit 1; }
[ -s data/donors/h1_2025_full/delta.npy ] || { log "no h1_2025_full bank -- stopping"; exit 1; }
[ "$(shasum -a 256 data/donors/h1_2025_full/delta.npy | cut -d' ' -f1)" = 26662c527164ddd434ae19a13076e915d1f82de9683c75e27e0af1f5b28591c0 ] \
  || { log "h1_2025_full delta sha unexpected -- stopping"; exit 1; }
[ "$(shasum -a 256 "$R453/whitener_union_pruned_def.npz" | cut -d' ' -f1)" = 66610be766ff44b266e076dbf5434833b79730568d030e4e4a22f83c1ffaaebe ] \
  || { log "453's day operator is not 66610be7 -- stopping"; exit 1; }
mkdir -p "$W/run"
cp -c "$R453/whitener_union_pruned_def.npz" "$R453/whitener_union_pruned_def.json" "$W/run/"
log "builder ${BF:0:16} = HEAD; def_day.sh ${DD:0:16} = HEAD (pass-through present); starting DEF_DAY_PROGRAMME_SOURCE=h1_2025_full def_day.sh (flags $FLAGS)"
DEF_DAY_PROGRAMME_SOURCE=h1_2025_full bash scripts/def_day.sh "$BUNDLE" "$W/run" "$N" "$Q" auto "$FLAGS" 4 > "$W/def_day.nohup" 2>&1
/usr/bin/grep -q "DONE. Next" "$W/run/def_day.log" || { log "def_day did not finish: $(tail -2 "$W/run/def_day.log" | tr '\n' ' ')"; exit 1; }
log "def_day: $(/usr/bin/grep -E "programme-source|CONVERGED at|NOT converged|final:" "$W/run/def_day.log" | tr '\n' ' ' | cut -c1-400)"
STEM="submissions/461dr-october-path-arcprog-abc-q${Q}"
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
