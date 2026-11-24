#!/usr/bin/env bash
# 670 (coordinator [2d6f46], 2026-10-16): release [4592ee]'s chain690j (proof of the lifted LODO RSS guard + the R-NEW rerun) only after
# build707 has been seen alive and has then EXITED, so the 4.5 h rerun cannot jump R7a (window ends 18:00Z) or split build707's levels.
# Order: R-KA → R7a → build707 → GO-690j → (R-NEW rerun; rehearse708 / build710 / build711 follow on their own gates). Cap 10-18 12:00Z.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
LOG=experiments/670-coord-1015/work/go690j.log; GO=experiments/690-def-readiness-1015/work/GO-690j
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
alive(){ ps -axo comm,args | awk '$1 ~ /^bash$/ && /experiments\/707-strong-tier-direction\/build707\.sh/' | /usr/bin/grep -q .; }
cap(){ [ "$(date -u +%Y%m%d%H%M)" -ge 202610181200 ]; }
log "armed (pid $$): GO-690j once build707.sh has run and exited"
until alive; do cap && { log "cap before build707 started"; log "END"; exit 1; }; sleep 60; done
log "build707 alive"
while alive; do cap && { log "cap while build707 alive"; log "END"; exit 1; }; sleep 60; done
touch "$GO" && log "build707 exited: GO-690j written"
log "END"
