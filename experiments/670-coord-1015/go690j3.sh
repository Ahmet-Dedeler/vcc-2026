#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-16): go690j v3 = v2 (launched 10:31:33Z) + science's run701c/run701d must have exited first, so the
# October K rule's condition (i) (reach701.json) and the K 5,000 price exist before R-NEW's 4.5 h run. Releases [4592ee]'s chain690j only
# after build707 has EXITED, R7a's T4 resume has ended (R7a's gate needs "builders 0, chains 0", so an R-NEW started at build707's exit would
# have pushed it past its 18:00Z window), and run701c/run701d have exited (cap 20:00Z).
# Order: build707 → R7a T4 resume → run701c/d (via gate602) → GO-690j → R-NEW → 708 → 710 → 711.
# R7a has ended when its access.log gains an "END rc=" line, or when work/R7A-RELEASE exists (science or the coordinator writes it if the
# resume is abandoned), or at 10-16 18:15Z (the window has closed by then).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
W=experiments/670-coord-1015/work; LOG=$W/go690j.log; GO=experiments/690-def-readiness-1015/work/GO-690j
ACC=experiments/551-r6-c10b/work/r7a/access.log
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /experiments\/707-strong-tier-direction\/build707\.sh/' | /usr/bin/grep -q .; }
cap(){ [ "$(date -u +%Y%m%d%H%M)" -ge 202610181200 ]; }
r7a_done(){ /usr/bin/grep -q "END rc=" "$ACC" 2>/dev/null || [ -f "$W/R7A-RELEASE" ] || [ "$(date -u +%Y%m%d%H%M)" -ge 202610161815 ]; }
echo $$ > $W/go690j3.pid
log "v3 armed (pid $$): GO-690j once build707.sh has exited, R7a's T4 resume has ended, and run701c/run701d have exited"
r701(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /experiments\/701-kup-callset\/run701[cd]\.sh/' | /usr/bin/grep -q .; }
cap701(){ [ "$(date -u +%Y%m%d%H%M)" -ge 202610162000 ]; }
while alive; do cap && { log "cap while build707 alive"; log "END"; exit 1; }; sleep 60; done
log "build707 not running"
until r7a_done; do cap && { log "cap while R7a open"; log "END"; exit 1; }; sleep 60; done
log "R7a ended: $(/usr/bin/grep "END rc=" "$ACC" 2>/dev/null | tail -1 | cut -c1-120)$([ -f "$W/R7A-RELEASE" ] && echo ' R7A-RELEASE present')"
while r701; do cap701 && { log "run701c/d still alive at 20:00Z: releasing R-NEW anyway"; break; }; sleep 60; done
log "run701c/d not running: $(tail -1 experiments/701-kup-callset/work/run701d.log 2>/dev/null | cut -c1-100)"
touch "$GO" && log "GO-690j written"
log "END"
