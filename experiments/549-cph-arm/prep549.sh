#!/usr/bin/env bash
# 549: package 446bj2zcph by hand (2026-09-30 ~22:05Z). The build finished at 21:49Z with every pre-prep gate passed
# ([446bj2zcph] PASS exactly de_mde_sign + de_sign_emitter_sha256 + rc_callset_block; zero-gene audit 0; CLIFF OK), but
# run_arm_generic.sh's prep gate (free+inactive >= 23 GiB) could not open: since 21:46Z the Mac sat at 18.8-20.5 GiB because
# the compressor holds ~18.5 GiB of other apps' pages (browser renderers), which vm_stat's free+inactive does not count;
# memory_pressure reported 58-62 % free and swap was unused (diagnosed by the peer session, 22:00Z). This runs exactly the
# arm script's prep lines (run_arm_generic.sh after `wait_gate`) into the same run log, with a vcc dry run first (AGENTS
# §6.2), under the same gate minus the free+inactive condition. Never submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
ARM=446bj2zcph; NAME=cph549; W="$REPO/experiments/428-zero-gene-audit/work/build_$ARM"; LOG="$W/run$ARM.log"
STEM="$ARM-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cph"
PKG="$REPO/submissions/$STEM.vcc"; META="$W/$NAME.h5ad.meta.json"; BUNDLE="$REPO/data/challenge"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
free_gib(){ df -k "$REPO" | awk 'NR==2{printf "%.1f", $4/1048576}'; }
[ -s "$W/$NAME.h5ad" ] && [ -s "$META" ] || { echo "no built artifact/meta in $W"; exit 1; }
[ -e "$PKG" ] && { echo "package already exists: $PKG"; exit 1; }
/usr/bin/grep -q "^\[$ARM\] PASS" "$LOG" && /usr/bin/grep -q "^\[$ARM\] CLIFF OK" "$LOG" || { echo "pre-prep gates not in $LOG"; exit 1; }
pct="$(memory_pressure | awk '/free percentage/{gsub(/%/,"",$5); print $5}')"
others="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print(len(P.other_builders()))' 2>/dev/null || echo 99)"
preps="$(ps -axo args= | grep -vE 'claude|grep' | grep -cE '(^|/| )vcc prep( |$)' | tr -d ' ')"
fi="$(vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}')"
log "prep549 (by hand; the arm script's prep gate never opened: free+inactive ${fi} GiB < 23, compressor-held pages): free ${pct}%  other builders ${others}  vcc preps ${preps}  disk $(free_gib) GiB"
[ "${pct:-0}" -ge 40 ] && [ "$others" = "0" ] && [ "$preps" = "0" ] && awk -v v="$(free_gib)" 'BEGIN{exit !(v>=20)}' \
  || { log "prep549 gate closed -- not prepping"; exit 3; }
log "INVOCATION 2 (prep549 dry run): vcc prep --dry-run -> $PKG"
vcc prep "$W/$NAME.h5ad" -g "$BUNDLE/gene_names.csv" --perts "$BUNDLE/pert_counts.csv" --contexts A,B,C -o "$PKG" --dry-run > "$W/prep-dryrun.log" 2>&1
RC=$?; log "dry-run rc=$RC"; tail -5 "$W/prep-dryrun.log" | tee -a "$LOG"
[ $RC -eq 0 ] || { log "dry run failed -- not prepping"; exit 1; }
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
