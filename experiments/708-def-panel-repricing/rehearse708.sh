#!/usr/bin/env bash
# 708 rehearsal (PREREGISTRATION "Rehearsal before 10-21"; the coordinator's condition: it must reproduce 438's abcdef
# outputs exactly, or 708 is OFF). Runs once the science chain (sci1016.sh) has written a NEW "END" and its process is
# gone (step 7 of that chain, as a separate waiter so the live chain is never edited), then strictly serially, each step
# behind gate602 (FIFO; it refuses beside a builder):
#   A. dayof438.py on 180's relabelled A/B/C bundle (the D/E/F code path: contexts D/E/F, --synthetic) with 438's abcdef
#      inputs (A/B/C production operator, A/B/C CD4 banks, 062's held fold operators, board q 1.25), graded against
#      work/expected_abcdef.txt by kaA708.py (AMENDMENT 5: like-labelled; claim-independent cells exact, rho by the claim
#      ratio, RECOMMENDATION identical).
#   B. dump708.py on the same inputs, one process per fold, with --ka438 = 438's stored abcdef comps (|d| < 1e-12 or it
#      aborts) and --designs h2.5c1.5 (abcdef's Part A winner, for AMENDMENT 1's confirmation).
#   C. score708.py on the rehearsal tag (informational on A/B/C's panel; it proves the day's scoring path end to end).
# Terminal line, always: "<ts> REHEARSAL PASS" | "<ts> REHEARSAL FAIL <why>", then "<ts> END".
#   nohup bash scripts/detach.sh bash experiments/708-def-panel-repricing/rehearse708.sh > experiments/708-def-panel-repricing/work/rehearse708.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/708-def-panel-repricing; W=$X/work; mkdir -p "$W"; LOG=$W/rehearse708.log
T0="$(date -u +%FT%TZ)"; CAP="${R708_CAP:-202610201200}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
newline(){ awk -v t0="$1" -v re="$3" 'substr($0, 1, 20) > t0 && $0 ~ re {f = 1} END {exit f ? 0 : 1}' "$2" 2>/dev/null; }
TERM_=0; trap 'rc=$?; [ "$TERM_" = 1 ] || log "REHEARSAL FAIL unexpected exit rc $rc"; log "END"' EXIT
fail(){ TERM_=1; log "REHEARSAL FAIL $1"; exit 1; }
released(){ [ -e experiments/691-release-watcher/work/RELEASE-DETECTED ]; }   # coordinator 10-16: after the D/E/F release no new A/B/C build may START (it would delay R0)
relguard(){ released && { TERM_=1; log "REHEARSAL $1 skipped: release detected (708 runs on the real bundle after R0 instead)"; exit 0; }; return 0; }
echo $$ > "$W/rehearse708.pid"
# 10-16 re-arm (AMENDMENT 5): after the production queue (build711, R-NEW via chain690j/run690, build710) is clear, over two
# checks 60 s apart; R-ABC2's chain (rabc2_695) is not waited for here (rehearse708b yields to it)
quiet(){ ! alive experiments/710-h1-dose-board/build710.sh && ! alive experiments/711-kup-dose/build711.sh && ! alive chain690j \
         && ! alive run690.sh && ! alive after690.sh && ! alive scripts/def_day.sh && ! alive run_arm_generic.sh; }
log "armed (pid $$, T0 $T0, cap $CAP): waiting for the production queue (build711, R-NEW, build710) to clear"
until quiet && sleep 60 && quiet; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && fail "the production queue never cleared by $CAP"
  sleep 60
done
log "queue clear"
B=experiments/180-def-readiness/work/def-bundle
OP=data/donors/whitener_union_pruned.npz
CD4=data/external/cd4t/cloud-compact/cd4_compact
HO=experiments/062-cd4-post-whiten/cloud-output/whiteners
TAG=r708abcdef
G(){ bash experiments/602-council-gate/gate602.sh "$1" "$2" -- "${@:3}"; }
# ---- A (AMENDMENT 5: graded by kaA708.py, like-labelled; the 12:24Z output is reused if its bytes are the graded ones)
relguard A
GRADED=7f054544802225f8306762ddf581303f39385d2327fd0f293f65fa9f74a90689
if [ -s "$W/dayof438-$TAG.txt" ] && [ "$(shasum -a 256 "$W/dayof438-$TAG.txt" | cut -d' ' -f1)" = "$GRADED" ]; then
  log "A: reusing the graded Part A output (sha256 ${GRADED:0:16}, AMENDMENT 5, commit 20cb2bc8)"
else
  log "A: dayof438 on $B (tag $TAG, held fold ops, q 1.25)"
  GATE602_QUIET=1 G 708-dayof438 9 uv run python experiments/438-october-stack-repricing/dayof438.py --bundle "$B" --op "$OP" \
    --cd4-root "$CD4" --tag "$TAG" --fold-ops held --held-ops "$HO" --synthetic --board-q 1.25 > "$W/dayof438-$TAG.txt" 2>&1 \
    || fail "dayof438 rc $? (see $W/dayof438-$TAG.txt)"
fi
uv run python $X/kaA708.py "$W/dayof438-$TAG.txt" --out="$W/kaA708-$TAG.json" >> "$LOG" 2>&1 \
  || fail "A: PART A KA FAIL (Amendment 5; see $W/kaA708-$TAG.json)"
log "A PASS: PART A KA PASS (Amendment 5): unlevelled cells identical, rho = ref x claim ratio, RECOMMENDATION identical"
# ---- B
for f in k562_gwps orion_hct116 orion_hek293t; do
relguard B
  log "B: dump708 fold $f"
  GATE602_QUIET=1 G "708-dump-$f" 9 uv run python $X/dump708.py --tag "$TAG" --bundle "$B" --op "$OP" --cd4-root "$CD4" \
    --fold-ops held --held-ops "$HO" --synthetic --designs h2.5c1.5 --ka438 experiments/438-october-stack-repricing/work/abcdef "$f" \
    > "$W/dump708-$TAG-$f.log" 2>&1 || fail "dump708 $f rc $? (see $W/dump708-$TAG-$f.log)"
  log "B $f: $(/usr/bin/grep -E 'KA ' "$W/dump708-$TAG-$f.log" | tr '\n' ' ' | cut -c1-300)"
done
log "B PASS: the 438 KA reproduced on 3/3 folds"
# ---- C
relguard C
GATE602_QUIET=1 G 708-score 6 uv run python $X/score708.py "$TAG" > "$W/score708-$TAG.txt" 2>&1 || fail "score708 rc $?"
log "C: $(/usr/bin/grep -E '^(r |s |rs |A:)' "$W/score708-$TAG.txt" | tr '\n' ' ' | cut -c1-600)"
TERM_=1; log "REHEARSAL PASS"
