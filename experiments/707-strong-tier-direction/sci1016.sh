#!/usr/bin/env bash
# The science lane's builder chain after the 10-15 restart, in the coordinator [2d6f46]'s order (10-15 ~21:30Z):
#   1. wait: [17adab]'s rabc694 (R-ABC, 10-16 s1) has written a NEW terminal line (" READY rabc " / " NOT READY ") after this
#      chain started, and its process is gone
#   2. build701.sh            K-UP ku30 (10-16 s2)
#   3. r7a_access_1015b.sh    R7a, the sealed reserve's single access (resumes; its own gates hold it off any builder)
#   4. build705.sh            lstr705, the LODO strength rider safety read (10-17 s1)
#   5. wait: readiness's run690 (R-NEW + R-KA; launched by after690.sh once the coordinator's 690 HOLD is gone, which stays
#      until build705 exits) has written a NEW "END" after step 4 finished, and neither after690.sh nor run690.sh is alive
#   6. build707.sh            sw707, the strong-tier direction read (10-17 s2 / 10-18 s1)
# A step's NOT READY does not stop the chain: each later step is an independent read with its own gates and terminal lines.
# Launch (detached, so an app quit cannot kill it):
#   nohup bash scripts/detach.sh bash experiments/707-strong-tier-direction/sci1016.sh > experiments/707-strong-tier-direction/work/sci1016.nohup 2>&1 &
# Release guard (coordinator 10-16): before starting each step, RELEASE-DETECTED ends the chain ("skipped: release detected").
# Terminal line, always: "<ts> END". Process matching uses bash + script path (never pgrep -f; AGENTS 6.2).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/707-strong-tier-direction/work; mkdir -p "$W"; LOG=$W/sci1016.log
T0="${SCI_T0:-$(date -u +%FT%TZ)}"   # a re-arm passes the original T0 so an upstream line written meanwhile still counts
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk -v p="$1" '$1 ~ /^bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
# a line written after time $1 (ISO, compared as text) that matches the ERE $3, in file $2
newline(){ awk -v t0="$1" -v re="$3" 'substr($0, 1, 20) > t0 && $0 ~ re {f = 1} END {exit f ? 0 : 1}' "$2" 2>/dev/null; }
capped(){ [ "$(date -u +%Y%m%d%H%M)" -ge "$1" ]; }
released(){ [ -e experiments/691-release-watcher/work/RELEASE-DETECTED ]; }   # coordinator 10-16: after the D/E/F release no new A/B/C build may START (it would delay R0)
relguard(){ released && { log "$1 skipped: release detected"; exit 0; }; return 0; }
echo $$ > "$W/sci1016.pid"
trap 'log "END"' EXIT
log "armed (pid $$, T0 $T0)"

# ---- 1. R-ABC (rabc694) must finish first
RL=experiments/694-rabc/work/rabc694.log
until newline "$T0" "$RL" " (READY rabc|NOT READY) " && ! alive rabc694.sh; do
  capped 202610161000 && { log "rabc694 wrote no new terminal line by 10-16 10:00Z: continuing (build701's own gates still hold off any builder)"; break; }
  sleep 60
done
log "1 rabc694 done: $(tail -1 "$RL" | cut -c1-220)"

# ---- 2. K-UP
relguard "2 build701"
log "2 START build701"
bash experiments/701-kup-callset/build701.sh > experiments/701-kup-callset/work/build701.nohup 2>&1
log "2 build701 exited rc=$?: $(tail -2 experiments/701-kup-callset/work/build701.log | tr '\n' ' ' | cut -c1-300)"

# ---- 3. R7a (single access; resumable)
relguard "3 R7a"
log "3 START R7a wrapper"
bash experiments/552-r7-reserve/work/r7a_access_1015b.sh > experiments/552-r7-reserve/work/r7a_access_1015b.nohup 2>&1
log "3 R7a wrapper exited rc=$?: $(tail -2 experiments/551-r6-c10b/work/r7a/access_wrapper_1015.log | tr '\n' ' ' | cut -c1-300)"

# ---- 4. lstr705
relguard "4 build705"
log "4 START build705"
bash experiments/705-lodo-strength-safety/build705.sh > experiments/705-lodo-strength-safety/work/build705.nohup 2>&1
log "4 build705 exited rc=$?: $(tail -2 experiments/705-lodo-strength-safety/work/build_lstr705.log | tr '\n' ' ' | cut -c1-300)"
T4="$(date -u +%FT%TZ)"

# ---- 5. [17adab]'s run690 (R-NEW + R-KA) runs next; wait for its NEW terminal line
PL=experiments/690-def-readiness-1015/work/run690.log
until newline "$T4" "$PL" " END$" && ! alive after690.sh && ! alive run690.sh; do
  capped 202610171600 && { log "run690 wrote no END after build705 by 10-17 16:00Z: build707 NOT started (coordinator decides)"; exit 0; }
  sleep 60
done
log "5 run690 done: $(tail -1 "$PL" | cut -c1-200)"

# ---- 6. sw707
relguard "6 build707"
log "6 START build707"
bash experiments/707-strong-tier-direction/build707.sh > experiments/707-strong-tier-direction/work/build707.nohup 2>&1
log "6 build707 exited rc=$?: $(tail -2 experiments/707-strong-tier-direction/work/build_sw707.log | tr '\n' ' ' | cut -c1-300)"
