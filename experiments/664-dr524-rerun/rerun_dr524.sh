#!/usr/bin/env bash
# 664 (coordinator, 2026-10-05 ~02:10Z): re-run 524's October-recipe dress rehearsal after a memory refusal.
# recipe_dr524.sh released at 02:00:41Z (its gate: free+inactive >= 25 GiB). 31 s later, def_pipeline's own pre-build check read
# free+inactive(+speculative+purgeable) 22.5 GiB < 23 and refused ("setup did not complete"); the builder was restored
# (4daae237 == HEAD). Two iOS Simulators (GlowUp) keep free+inactive at 19-26 GiB, and setup's identity/key stages (2.7 / 4.1 GiB
# peaks) eat the 2 GiB margin the 25 GiB gate leaves. A full def_day holds ~20 GiB for ~24 min (497 efr: build 20.2, prep 22.2).
# This wrapper re-runs recipe_dr524.sh UNCHANGED (same script, same defaults; it moves the failed run dir aside itself and
# re-resolves the recipe, with mde_sign still PENDING -> OFF as 658 assumes; the coordinator decides mde_sign only after this
# rehearsal exits). It releases only when free+inactive >= 27 GiB on two readings 60 s apart, with no other build or prep
# running. It retries ONLY on that same memory refusal, after a 5-minute pause, and stops on any other failure. Cap: 10-05
# 19:00Z (a ~40 min build then still ends before chain1006 decides at 20:00Z; chain1006 waits for the pid file's process
# anyway, until 23:30Z).
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=experiments/664-dr524-rerun; mkdir -p $X/work; LOG=$X/work/rerun_dr524.log
D=experiments/524-october-recipe-dr; DL=$D/work/recipe_dr524.log; CAP=202610051900
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%.1f",(f+i)*16384/2^30}'; }
busy(){ uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; raise SystemExit(0 if P.other_builders() else 1)' \
        || ps -axo args= | grep -vE 'claude|grep' | grep -qE '(^|/| )vcc prep( |$)'; }
echo $$ > $X/work/rerun_dr524.pid
log "armed (pid $$): release at free+inactive >= 27 GiB twice 60 s apart, no build/prep running; retry only on def_pipeline's memory refusal; cap $CAP"
N=0
while [ "$(date -u +%Y%m%d%H%M)" -lt $CAP ]; do
  a="$(fi_gib)"
  if awk -v v="$a" 'BEGIN{exit !(v>=27)}' && ! busy; then
    sleep 60; b="$(fi_gib)"
    if awk -v v="$b" 'BEGIN{exit !(v>=27)}' && ! busy; then
      N=$((N + 1)); log "attempt $N: free+inactive $a then $b GiB -> bash $D/recipe_dr524.sh"
      bash $D/recipe_dr524.sh > $X/work/attempt$N.nohup 2>&1; RC=$?
      T="$(tail -3 $DL | tr '\n' ' ' | cut -c1-400)"
      log "attempt $N: rc=$RC; recipe_dr524.log tail: $T"
      if [ $RC -eq 0 ]; then log "rehearsal finished rc 0"; log "END"; exit 0; fi
      if /usr/bin/grep -q "refusing" $D/work/run/driver-setup.log 2>/dev/null && /usr/bin/grep -q "setup did not complete" $D/work/run/def_day.log 2>/dev/null; then
        log "attempt $N: def_pipeline's memory refusal again -- retrying after 5 min"; sleep 300; continue
      fi
      log "attempt $N failed for another reason -- stopping for diagnosis"; log "END"; exit 1
    fi
  fi
  sleep 60
done
log "cap $CAP reached without a finished rehearsal ($N attempts)"; log "END"; exit 1
