#!/usr/bin/env bash
# One-shot (coordinator, 2026-10-04): once chain1005 (pid 24465) exits with 10-05 slot 2 empty (rule1005's s2 list is
# exhausted: 498p2g and 460dr are scored), and idall3's build log ends READY, arm 496/arm1005.sh slot 2 with idall3.
# arm1005 slot 2 waits for slot 1 (496a) to be captured. Nothing is armed if idall3 is not READY by 10-05 04:00Z.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/663-idall3-arm; W=$X/work; LOG=$W/fill1005s2.log; C=experiments/496-uncalled-sign-hook/work/chain1005.log
A=experiments/496-uncalled-sign-hook; BL=$W/build_idall3.log
STEM=idall3-cap15-partial50-clean4q80-lodo6-ballast-sharpenp2g-stacklodo05hi6load-ctxlevelr-mde010coarse-zr5-idhook
RL=experiments/428-zero-gene-audit/work/build_idall3/runidall3.log
NAME="There's nothing sadder than living your whole life by the ocean but never daring to set sail."
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/fill1005s2.pid
log "armed (pid $$): waiting for chain1005 (pid 24465) to exit, and for idall3 READY"
while kill -0 24465 2>/dev/null; do [ "$(date -u +%Y%m%d%H%M)" -ge 202610050400 ] && { log "cap 10-05 04:00Z with chain1005 running -- nothing armed"; log "END"; exit 0; }; sleep 30; done
LINE="$(/usr/bin/grep -E 'armed 10-05: s1' "$C" | tail -1)"
log "chain1005 exited; its arm line: ${LINE:-<none>}"
case "$LINE" in
  "") log "chain1005 armed nothing -- not arming a slot 2 fill"; log "END"; exit 0;;
  *"s2 <none>"*) ;;
  *) log "chain1005 armed its own slot 2 -- no fill"; log "END"; exit 0;;
esac
until tail -1 "$BL" 2>/dev/null | /usr/bin/grep -qE " (READY|NOT READY) "; do
  [ "$(date -u +%Y%m%d%H%M)" -ge 202610050400 ] && { log "idall3 not READY by 10-05 04:00Z -- nothing armed"; log "END"; exit 0; }
  sleep 60; done
T="$(tail -1 "$BL")"
case "$T" in *" READY idall3 "*) ;; *) log "idall3 build ended: $T -- nothing armed"; log "END"; exit 0;; esac
nohup bash $A/arm1005.sh 2 idall3 "$RL" "submissions/$STEM.vcc" "$NAME" "idall3 slot 2 on 10-05 (coordinator fill)" > $A/work/arm1005_s2.nohup 2>&1 & echo $! > $A/work/arm1005_s2.pid
sleep 3; log "armed idall3 as 10-05 slot 2 (pid $(cat $A/work/arm1005_s2.pid))"; log "END"
