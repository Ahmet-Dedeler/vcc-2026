#!/usr/bin/env bash
# 690j (10-16, readiness [4592ee]): prove def_lodo_key's RSS-guard lift (scripts/def_lodo_key.py, DEF_DAY_LODO_MAX_RSS_GIB,
# default 16) on its own known answers, then rerun R-NEW, which died at 07:27:51Z on the old 8.0 GiB guard ("cd4_rest
# centering: peak RSS 8.530 GiB reached 8.000 GiB"). The coordinator orders the builder queue, so this starts only after:
#   - $W/GO-690j exists (the coordinator touches it when the slot is free), AND
#   - no run690 / after690 / R7a wrapper / def_day / run_arm / run_ka chain alive under bash, no builder or operator fit,
#     free+inactive >= 20 GiB AND memory_pressure free >= 40 % (AGENTS §9); cap CAP690J (default 10-19 12:00Z)
# Then:
#   1. def_lodo_key --ka (A/B/C, production operator e2f65ecf, rules abc, vs 420's 8c82bef3) and --ka-pct (vs 495's
#      dc374281), each on a FRESH lib438 tag so they recompute through the lifted guard; both must print PASS (rc 0)
#   2. lib438's cache of the failed R-NEW (work/lodo_new_def) is moved aside, then this process becomes
#      `ONLY690=rnew run690.sh` (R-NEW alone; R-KA already proved this tree). run690.log takes over.
# Log: work/chain690j.log; verdict lines "690j KA PASS|FAIL ..."; terminal "<ts> END" or "<ts> exec run690 ONLY690=rnew".
# usage: nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/chain690j.sh >> experiments/690-def-readiness-1015/work/chain690j.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; W=$E/work
if [ -z "${CHAIN690J_FROZEN:-}" ]; then
  FZ="$W/.frozen-$$-chain690j.sh"; cp "$REPO/$E/chain690j.sh" "$FZ" || exit 1
  CHAIN690J_FROZEN=1 exec bash "$FZ" "$@"
fi
LOG=$W/chain690j.log; CAP="${CAP690J:-202610191200}"; KO=$W/ka690j; TS="$(date -u +%Y%m%dT%H%M%S)"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'rm -f "$W/.frozen-$$-chain690j.sh"; log "END"' EXIT
alive(){ ps -axo comm=,args= | awk -v p="$1" '$1 ~ /bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
fi_gib(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(int(d.free_gib()))" 2>/dev/null || echo 0; }
n_builders(){ uv run python -c "import sys; sys.path.insert(0, 'scripts'); import def_pipeline as d; print(len(d.other_builders()))" 2>/dev/null || echo 99; }
gate(){
  [ -e "$W/GO-690j" ] || { echo "no GO-690j (the coordinator's)"; return 1; }
  # r7a_after690.sh: science's R7a waiter (alive until R7a's wrapper ends; it waits only for run690's END, never for this
  # chain). Without it, the 60-70 s between run690's END and the wrapper's launch could let R-NEW's rerun jump R7a.
  local c; for c in run690.sh after690.sh r7a_after690.sh r7a_access_1015b.sh scripts/def_day.sh run_arm_generic.sh run_ka_527.sh; do
    alive "$c" && { echo "$c alive"; return 1; }; done
  [ "$(n_builders)" = 0 ] || { echo "a builder / operator fit is running"; return 1; }
  local f mp; f="$(fi_gib)"; mp="$(memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}')"
  [ "$f" -ge 20 ] || { echo "free+inactive $f GiB < 20"; return 1; }
  [ "${mp:-0}" -ge 40 ] || { echo "memory_pressure free ${mp:-?} % < 40"; return 1; }
  echo "GO, nothing alive, free+inactive $f GiB, pressure free $mp %"
}
echo $$ > "$W/chain690j.pid"
log "chain690j armed (pid $$): KA of def_lodo_key's RSS-guard lift, then ONLY690=rnew; cap $CAP"
last=""
until r="$(gate)"; do
  [ "$r" != "$last" ] && log "waiting: $r"; last="$r"
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "690j STOP: not started by cap $CAP ($r)"; exit 3; }
  sleep 60
done
log "gate PASS: $r"
mkdir -p "$KO"
OP=data/donors/whitener_union_pruned.npz; CD4=data/external/cd4t/cloud-compact/cd4_compact
/usr/bin/time -l uv run python scripts/def_lodo_key.py --bundle data/challenge --op "$OP" --cd4-root "$CD4" --rules abc --ka \
  --tag "ka690j_${TS}" --out "$KO/ka.npz" > "$KO/ka.log" 2>&1; R1=$?
/usr/bin/time -l uv run python scripts/def_lodo_key.py --bundle data/challenge --op "$OP" --cd4-root "$CD4" --rules abc --ka-pct \
  --hi-pct 96.64697522769474 --tag "ka690jp_${TS}" --out "$KO/ka_pct.npz" > "$KO/ka_pct.log" 2>&1; R2=$?
peak(){ awk '/maximum resident set size/{printf "%.2f GiB", $1/2^30}' "$1"; }
L1="$(/usr/bin/grep -E '^KA (PASS|FAIL)' "$KO/ka.log" | tail -1 | cut -c1-200)"; L2="$(/usr/bin/grep -E '^KA-pct (PASS|FAIL)' "$KO/ka_pct.log" | tail -1 | cut -c1-200)"
G1="$(/usr/bin/grep -E 'RSS guard lifted' "$KO/ka.log" | tail -1)"
log "KA rc $R1 peak $(peak "$KO/ka.log"): ${L1:-no KA line}"
log "KA-pct rc $R2 peak $(peak "$KO/ka_pct.log"): ${L2:-no KA-pct line}"
log "guard: ${G1:-no lift line}"
if [ "$R1" = 0 ] && [ "$R2" = 0 ] && [ "${L1#KA PASS}" != "$L1" ] && [ "${L2#KA-pct PASS}" != "$L2" ] && [ -n "$G1" ]; then
  log "690j KA PASS: def_lodo_key with the lifted guard reproduces 8c82bef3 and dc374281"
else
  log "690j KA FAIL: R-NEW not rerun (see $KO/ka.log, $KO/ka_pct.log)"; exit 1
fi
C=experiments/438-october-stack-repricing/work/lodo_new_def
[ -e "$C" ] && { mv "$C" "$C.failed-0727Z-$TS"; log "moved the failed R-NEW's lib438 cache aside: $C"; }
trap - EXIT; rm -f "$W/.frozen-$$-chain690j.sh"
log "exec run690 ONLY690=rnew (CAP690 202610210000); its log: $W/run690.log"
ONLY690=rnew CAP690=202610210000 exec bash "$REPO/$E/run690.sh"
