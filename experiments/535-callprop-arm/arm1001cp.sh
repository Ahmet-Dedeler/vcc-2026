#!/usr/bin/env bash
# 10-01 slot 2 for 446bj2zcp (rule amendment 2026-09-30 06:38Z, before any 10-01 entry exists): after the 09-30 HARM read the
# 10-01 rule (464/rule1001b.py) armed only s1 = 460dr, so slot 2 was empty. This waits for build535's PROCESS to exit, requires
# its arm's gates ([446bj2zcp] PASS, CLIFF OK, novelty rc=0, a package sha equal to the file's), then arms 464's arm1001.sh as
# slot 2 (it waits for slot 1's capture on 10-01 and has 536's retry/resume loop). If build535 fails or is not done by
# 10-01 00:00Z, nothing is armed. Never builds, never submits by itself. Named arm1001cp.sh (no chain5NN trap).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/535-callprop-arm; W=$X/work; LOG=$W/arm1001cp.log; E=experiments/428-zero-gene-audit/work
TAG=446bj2zcp; RL=$E/build_$TAG/run$TAG.log
PKG=submissions/$TAG-cap15-partial50-clean4q80-pertscale3-ballast-sharpen175-stacklodo05load-ctxlevelm-mde010coarse-zr5-zs250-cp.vcc
NAME="I see problems and I solve them."
# 2026-09-30 10:55Z: optional args name the build to wait for (pid file, log); default build535. Re-armed on rebuild535
# after build535's arm never got its package (prep gate closed on disk, then the prediction was removed).
PIDF="${1:-$W/build535.pid}"; BLOG="${2:-$W/build535.log}"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
echo $$ > $W/arm1001cp.pid
log "armed (pid $$): waiting for the process in $PIDF to exit (cap 10-01 00:00Z)"
alive(){ local p; p="$(cat "$PIDF" 2>/dev/null)"; [ -n "$p" ] && ps -p "$p" -o args= 2>/dev/null | /usr/bin/grep -q build535; }
until ! alive; do [ "$(date -u +%Y%m%d%H%M)" -ge 202610010000 ] && { log "$PIDF's build still running at 10-01 00:00Z -- no slot 2"; log "END"; exit 0; }; sleep 60; done
log "$(basename "$PIDF" .pid) exited; its last lines: $(tail -2 "$BLOG" | tr '\n' ' ' | cut -c1-300)"
[ -s "$PKG" ] && /usr/bin/grep -q "^\[$TAG\] PASS" "$RL" 2>/dev/null && /usr/bin/grep -q "^\[$TAG\] CLIFF OK" "$RL" && /usr/bin/grep -q "novelty rc=0" "$RL" \
  || { log "446bj2zcp not ready (package/gates/novelty) -- no slot 2"; log "END"; exit 0; }
WANT="$(/usr/bin/grep 'package sha256:' "$RL" | tail -1 | awk '{print $NF}')"; GOT="$(shasum -a 256 "$PKG" | cut -d' ' -f1)"
[ -n "$WANT" ] && [ "$WANT" = "$GOT" ] || { log "package sha ${GOT:0:16} != run log ${WANT:0:16} -- no slot 2"; log "END"; exit 0; }
uv run python scripts/check_submission_novelty.py check "$PKG" > /dev/null 2>&1 || { log "novelty re-check refused -- no slot 2"; log "END"; exit 0; }
nohup bash experiments/464-slot-1001/arm1001.sh 2 "$TAG" "$RL" "$PKG" "$NAME" "$TAG slot 2 on 10-01 (535; 09-30 06:38Z amendment)" > $W/arm1001cp_s2.nohup 2>&1 & echo $! > $W/arm1001cp_s2.pid
sleep 5
log "armed 10-01 slot 2: $TAG (pid $(cat $W/arm1001cp_s2.pid))"
log "END"
