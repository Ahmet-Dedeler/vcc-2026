#!/usr/bin/env bash
# 708 DAY-OF (10-22, AFTER R0 is submitted): re-price the panel-dependent switches on D/E/F's real panel.
#   nohup bash scripts/detach.sh bash experiments/708-def-panel-repricing/day708.sh <bundle> <def_day work dir> <tag> \
#       [--held-ops <dir with held_<line>.npz: skip the fits>] > experiments/708-def-panel-repricing/work/day708-<tag>.nohup 2>&1 &
# Env (the A/B/C board reads, from BOARD-LEDGER): R705=CLEARED|HARM|NONE, R707=TRANSFERS|WEAK|NULL|HARM|NONE.
# Steps, strictly serial, each behind gate602:
#   1. three held-out fold operators for THIS panel (build_union_whitener.py --donors <the other two production donors>
#      --basis-prune essential, ~5 min / 13-14 GB each), unless --held-ops is given;
#   2. Part A: dayof438.py (--fold-ops held, board q 1.25), parsed for its RECOMMENDATION;
#   3. Part B: dump708.py on the 3 folds (--designs = Part A's winner if it moved), floor708.py's floor-kept targets
#      (AMENDMENT 4), then score708.py on all and on floor-kept targets;
#   4. the decisions (PREREGISTRATION + AMENDMENTS 1-4; an action needs both versions), one line per switch, into
#      work/decision708-<tag>.txt.
# Inputs from the day's def_day run: <work>/whitener_union_pruned_def.npz and <work>/cd4_compact.
# Terminal line, always: "<ts> DECISIONS WRITTEN" | "<ts> FAIL <why>", then "<ts> END". It never builds or submits.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
B="${1:?bundle}"; DW="${2:?def_day work dir}"; TAG="${3:?tag}"; shift 3
HO=""; [ "${1:-}" = "--held-ops" ] && HO="${2:?held ops dir}"
X=experiments/708-def-panel-repricing; W=$X/work; mkdir -p "$W"; LOG=$W/day708-$TAG.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
TERM_=0; trap 'rc=$?; [ "$TERM_" = 1 ] || log "FAIL unexpected exit rc $rc"; log "END"' EXIT
fail(){ TERM_=1; log "FAIL $1"; exit 1; }
G(){ bash experiments/602-council-gate/gate602.sh "$1" "$2" -- "${@:3}"; }
OP="$DW/whitener_union_pruned_def.npz"; CD4="$DW/cd4_compact"
[ -s "$OP" ] || fail "no day operator at $OP"; [ -d "$CD4" ] || fail "no day CD4 banks at $CD4"
SYN=""; python3 -c 'import json,sys; m=json.load(open(sys.argv[1])); sys.exit(0 if m.get("partition","").startswith("synthetic") else 1)' "$B/manifest.json" && SYN="--synthetic"
log "armed (pid $$): bundle $B ${SYN:-(real)}; day op $(shasum -a 256 "$OP" | cut -c1-16); R705=${R705:-NONE} R707=${R707:-NONE}"
# ---- 1. held-out fold operators
if [ -z "$HO" ]; then
  HO="$W/heldops-$TAG"; mkdir -p "$HO"
  for h in k562_gwps orion_hct116 orion_hek293t; do
    case $h in k562_gwps) D=orion_hct116,orion_hek293t;; orion_hct116) D=k562_gwps,orion_hek293t;; orion_hek293t) D=k562_gwps,orion_hct116;; esac
    [ -s "$HO/held_$h.npz" ] && { log "1 held_$h exists: $(shasum -a 256 "$HO/held_$h.npz" | cut -c1-16)"; continue; }
    log "1 fitting held_$h (donors $D)"
    GATE602_QUIET=1 G "708-fit-$h" 14 uv run python scripts/build_union_whitener.py --bundle "$B" --donors "$D" \
      --basis-prune essential --out "$HO/held_$h.npz" --scratch "$HO/.scratch_$h" --max-rss-gib 14 > "$HO/fit_$h.log" 2>&1 \
      || fail "fit held_$h rc $? (see $HO/fit_$h.log)"
    log "1 held_$h: $(shasum -a 256 "$HO/held_$h.npz" | cut -c1-16)"
  done
fi
# ---- 2. Part A
log "2 Part A: dayof438 (tag $TAG)"
GATE602_QUIET=1 G 708-dayof438 9 uv run python experiments/438-october-stack-repricing/dayof438.py --bundle "$B" --op "$OP" \
  --cd4-root "$CD4" --tag "$TAG" --fold-ops held --held-ops "$HO" $SYN --board-q 1.25 > "$W/dayof438-$TAG.txt" 2>&1 \
  || fail "dayof438 rc $? (see $W/dayof438-$TAG.txt)"
REC="$(/usr/bin/grep -E '^RECOMMENDATION' "$W/dayof438-$TAG.txt" | tail -1)"
[ -n "$REC" ] || fail "dayof438 printed no RECOMMENDATION"
HX="$(echo "$REC" | sed -nE 's/.*H1 x([0-9.]+), CD4 x([0-9.]+).*/\1/p')"; CX="$(echo "$REC" | sed -nE 's/.*H1 x([0-9.]+), CD4 x([0-9.]+).*/\2/p')"
[ -n "$HX" ] && [ -n "$CX" ] || fail "cannot parse: $REC"
DES=""; [ "$HX" = 2 ] && [ "$CX" = 1 ] || DES="h${HX}c${CX}"
log "2 Part A: $REC -> design ${DES:-incumbent (h2c1)}"
# ---- 3. Part B
for f in k562_gwps orion_hct116 orion_hek293t; do
  GATE602_QUIET=1 G "708-dump-$f" 9 uv run python $X/dump708.py --tag "$TAG" --bundle "$B" --op "$OP" --cd4-root "$CD4" \
    --fold-ops held --held-ops "$HO" $SYN ${DES:+--designs "$DES"} "$f" > "$W/dump708-$TAG-$f.log" 2>&1 \
    || fail "dump708 $f rc $? (see $W/dump708-$TAG-$f.log)"
  [ -s "$W/dump708-$TAG-$f.json" ] && [ -s "$W/dump708-$TAG-$f.npz" ] || fail "dump708 $f wrote no npz/json (10-16 json fix: json708.default)"
  log "3 dump708 $f: $(/usr/bin/grep -E 'KA per-donor' "$W/dump708-$TAG-$f.log" | tail -1 | cut -c1-160)"
done
# 708 AMENDMENT 4: Arc's basal floor from the bundle's OWN controls (801's rule), then each fold's floor-kept targets; on 180's
# relabelled A/B/C bundle (the rehearsal) the floor must reproduce 801's A/B/C threshold
FKA=""; [ "$B" = "experiments/180-def-readiness/work/def-bundle" ] && FKA="--expect-thr 13.412203349132618"
GATE602_QUIET=1 G 708-floor 2 uv run python $X/floor708.py --bundle "$B" --tag "$TAG" $FKA > "$W/floor708-$TAG.log" 2>&1 \
  || fail "floor708 rc $? (see $W/floor708-$TAG.log)"
[ -s "$W/floormask708-$TAG.json" ] || fail "floor708 wrote no floormask708-$TAG.json"
log "3 floor: $(/usr/bin/grep -E '^(floor F|KA:|k562|orion)' "$W/floor708-$TAG.log" | tr '\n' ' ' | cut -c1-300)"
GATE602_QUIET=1 G 708-score 6 uv run python $X/score708.py "$TAG" --floor-json "$W/floormask708-$TAG.json" > "$W/score708-$TAG.txt" 2>&1 \
  || fail "score708 rc $?"
log "3 score708: $(/usr/bin/grep -E '^(r |s |rs |A:)' "$W/score708-$TAG.txt" | tr '\n' ' ' | cut -c1-600)"
log "3 score708 floor-kept: $(/usr/bin/grep -E '^floor (r |s |rs |A:|s\|r)' "$W/score708-$TAG.txt" | tr '\n' ' ' | cut -c1-600)"
log "3 $(/usr/bin/grep -E '^fold cap hi' "$W/score708-$TAG.txt" | cut -c1-200)"
# ---- 4. decisions
# AMENDMENT 5: on a real bundle decide708 requires the rehearsal on file; a synthetic (rehearsal) bundle is exempt
R708_REHEARSAL="$([ -n "$SYN" ] && echo 1 || echo 0)" \
uv run python $X/decide708.py "$TAG" "$W/dayof438-$TAG.txt" "${DES:-}" "${R705:-NONE}" "${R707:-NONE}" > "$W/decision708-$TAG.txt" 2>&1 \
  || fail "decide708 rc $?"
log "4 decisions: $(tr '\n' ' ' < "$W/decision708-$TAG.txt" | cut -c1-700)"
TERM_=1; log "DECISIONS WRITTEN ($W/decision708-$TAG.txt)"
