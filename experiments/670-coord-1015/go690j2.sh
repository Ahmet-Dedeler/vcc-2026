#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-16 10:40Z): go690j v2. Releases [4592ee]'s chain690j (lifted-guard proof + the R-NEW rerun) only
# after build707 has EXITED and R7a's T4 resume has ended. R7a's own gate needs "builders 0, chains 0", so a 4.5 h R-NEW started at
# build707's exit would have pushed R7a's resume past its 18:00Z window. Order: build707 → R7a T4 resume → GO-690j → R-NEW → 708 → 710 → 711.
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
echo $$ > $W/go690j2.pid
log "v2 armed (pid $$): GO-690j once build707.sh has exited AND R7a's T4 resume has ended"
while alive; do cap && { log "cap while build707 alive"; log "END"; exit 1; }; sleep 60; done
log "build707 not running"
until r7a_done; do cap && { log "cap while R7a open"; log "END"; exit 1; }; sleep 60; done
log "R7a ended: $(/usr/bin/grep "END rc=" "$ACC" 2>/dev/null | tail -1 | cut -c1-120)$([ -f "$W/R7A-RELEASE" ] && echo ' R7A-RELEASE present')"
touch "$GO" && log "GO-690j written"
log "END"
