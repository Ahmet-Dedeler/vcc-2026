#!/usr/bin/env bash
# 432 chain (RESUME 2, 2026-09-23 ~11:35Z): the builder is free (PI, 11:28Z); no PI signal is waited for any more. Runs
# build432.sh's phases strictly in order. The main path (a-d) runs first; the diagnostic dumps and the negative
# controls last. Critical phases stop the chain on failure; diagnostic phases log and continue.
# RESUMABLE: a phase whose "[<phase>] phase done" line is already in work/state432.log is skipped (a phase with no done
# line is re-run from scratch). Do not edit while running.
# No `pgrep -f` anywhere (AGENTS 6.2 trap: it matches agent prompts); build432.sh's gate uses def_pipeline.other_builders(),
# which requires a launcher executable (python/uv/time) with the builder's file name in its arguments.
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
X=$REPO/experiments/432-def-rehearsal-new-knobs; W=$X/work; LOG=$W/chain432.log; STATE=$W/state432.log
HEADB=2c90d9e1cc2c19fd427da35ffe6a6d1306d3097382ff7a343233423cb710573f
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" | tee -a "$LOG"; }
bsha(){ shasum -a 256 scripts/make_submission_pm_cd4.py | cut -d' ' -f1; }
log "START (RESUME 2; pid $$)"
# a reboot mid-build leaves one of OUR pins on disk: restore it from its own record, never from git
if [ "$(bsha)" != "$HEADB" ]; then
  cur="$(bsha)"; rec=""
  for r in "$W"/*/pin-artifacts.json; do
    [ -f "$r" ] && grep -q "\"builder_sha256_after\": \"$cur\"" "$r" && rec="$r"
  done
  if [ -n "$rec" ]; then log "builder $cur is OUR pin ($rec) -- restoring"
    uv run python scripts/def_pin_panel_artifacts.py --restore --record "$rec" 2>&1 | tee -a "$LOG"
  fi
  [ "$(bsha)" = "$HEADB" ] || { log "builder $(bsha) != HEAD $HEADB and not our pin -- STOP"; exit 1; }
fi
log "builder at START: $(bsha)"
done_already(){ grep -qE " \[$1\] phase done$" "$STATE" 2>/dev/null; }
run(){ local crit="$1"; shift
       if done_already "$1"; then log "phase $* already done (state432.log) -- skipped"; return 0; fi
       log "phase $* START"; bash "$X/build432.sh" "$@"; local rc=$?; log "phase $* rc=$rc"
       if [ $rc -ne 0 ] && [ "$crit" = 1 ]; then log "critical phase $* failed -- STOP"; exit 1; fi; }
# (a) base path: operator re-fit, 410 key for the panel, 433h minus LOAD/CTX
run 1 whitener
run 1 key
run 1 b0
run 0 level_b0
run 0 audit_b0
run 0 post_b0
run 0 prep_b0
# (b) the stack dump on D/E/F, LOADed; (c) its re-level string; (d) + zs
run 1 dump_stack
run 1 b1
run 1 level_b1
run 0 audit_b1
run 0 post_b1
run 1 b2
run 0 level_b2
run 0 audit_b2
run 0 post_b2
run 0 prep_b2
# diagnostics: 417's literal CD4 x1.5 recipe vs H1 x2 alone (the carrying-sizing trap), and the three refusals
run 0 dump_lit
run 0 dump_h1
run 0 neg_abcdump
run 0 neg_abckey
run 0 neg_abcctx
log "builder at END: $(bsha)"
log "END"
