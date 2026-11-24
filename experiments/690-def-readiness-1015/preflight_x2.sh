#!/usr/bin/env bash
# 690 preflight_x2 (readiness [401d9b], 2026-10-22; the coordinator's GO-X2 condition 3, SIXA 344c4c31): once the lab commits
# 811's banks extended to D/E/F's panel and re-pins BANKS_MANIFEST_SHA in hook811.py, run X2's hook (dumphook811c0.sh) on r0's
# stack dump (ada1baa2) with the D/E/F bundle and report whether it now runs. Read-only for every tracked file.
#   1. waits until hook811.py's pinned BANKS_MANIFEST_SHA differs from the 10-22 pin 37e34335…, that hook811.py, the banks manifest
#      and dumphook811c0.sh are committed and clean, and the manifest on disk has the pinned sha (the lab's commit has landed);
#   2. runs the hook into a fresh dir experiments/690-def-readiness-1015/work/def-day/hook811c0_pre-<ts>/ and logs rc, the output
#      sha, the hook's own info line, and the diff of the lab's commit touching 811 (paths only: nothing but banks + the pin may move);
#   3. "PREFLIGHT PASS: …" (rc 0, output written) or "PREFLIGHT FAIL: …".
# Cap PRE_CAP (default 202610240000). Terminal line, always: "<ts> END".
#   nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/preflight_x2.sh > experiments/690-def-readiness-1015/work/preflight_x2.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; X=experiments/811-lolo-model
LOG="${PRE_LOG:-$E/work/preflight_x2.log}"; CAP="${PRE_CAP:-202610240000}"
OLD=37e343351000b9b65ceb4adc6a18004563219908cf03b805ed8b27215e57d83b; R0W=experiments/692-def-r0/work/r0
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'log "END"' EXIT
pin(){ sed -nE 's/^BANKS_MANIFEST_SHA = "([0-9a-f]{64})".*/\1/p' "$X/hook811.py"; }
log "preflight_x2 armed (pid $$): waits for a re-pinned hook811.py (pin != ${OLD:0:16}), committed + clean; cap $CAP"
last=""
while :; do
  P="$(pin)"
  if [ -n "$P" ] && [ "$P" != "$OLD" ] && [ -z "$(git status --porcelain -- "$X/hook811.py" "$X/dumphook811c0.sh" "$X/work/banks/manifest.json")" ] \
     && git ls-files --error-unmatch "$X/work/banks/manifest.json" >/dev/null 2>&1 \
     && [ "$(shasum -a 256 "$X/work/banks/manifest.json" | cut -d' ' -f1)" = "$P" ]; then break; fi
  [ "$last" = w ] || log "waiting: pin ${P:0:16} (old ${OLD:0:16}), hook / manifest committed + clean + matching"; last=w
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "PREFLIGHT FAIL: no re-pinned hook by the cap"; exit 1; }
  sleep 120
done
C="$(git log -1 --format='%h %s' -- "$X/hook811.py" | cut -c1-200)"
log "re-pinned: BANKS_MANIFEST_SHA ${P:0:16} (commit $C)"
log "811 paths changed since the 10-22 dossier ad212d84: $(git diff --name-only ad212d84 HEAD -- "$X" | grep -v '^experiments/811-lolo-model/work/run-' | tr '\n' ' ' | cut -c1-600)"
REF="$E/work/def-day/hook811c0_pre-$(date -u +%H%M%S)"; mkdir -p "$REF"
/usr/bin/time -l bash "$X/dumphook811c0.sh" data/challenge-def "$R0W" "$R0W/dump_h1x2_cd4x1_s0.25.npz" "$R0W/identity.json" "$REF" > "$REF/hook.log" 2>&1; rc=$?
OUT="$REF/dump_h1x2_cd4x1_s0.25.npz"
log "hook rc $rc; $(/usr/bin/grep -E 'wrote|Error|refuse' "$REF/hook.log" | tail -2 | tr '\n' ' ' | cut -c1-300); peak $(/usr/bin/grep 'maximum resident' "$REF/hook.log" | tr -s ' ')"
if [ "$rc" = 0 ] && [ -s "$OUT" ]; then
  log "PREFLIGHT PASS: dumphook811c0 ($(shasum -a 256 "$X/dumphook811c0.sh" | cut -c1-8), hook811 $(shasum -a 256 "$X/hook811.py" | cut -c1-8), banks ${P:0:8}) ran on r0's stack dump ada1baa2 for D/E/F: $OUT ($(shasum -a 256 "$OUT" | cut -c1-16))"
else
  log "PREFLIGHT FAIL: rc $rc (see $REF/hook.log)"; exit 1
fi
