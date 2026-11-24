#!/bin/bash
# gate602.sh -- shared memory gate for the 600-639 council's heavy jobs (> ~4 GB), so four sessions never stack
# 7-25 GB jobs on one 64 GB Mac (a job that runs out of memory here dies without a traceback; AGENTS 9).
#
# usage:  bash experiments/602-council-gate/gate602.sh <label> <max_gb> -- <command ...>
#   <label>   short id for the log, e.g. 611-h1full-cpm-s0
#   <max_gb>  the job's expected peak RSS in GB (integer); a slot opens only if free+inactive >= max(20, max_gb + 12) GiB
#             (the coordinator's peer rule: <= 8 GB jobs at >= 20 GiB)
#
# Rules (from AGENTS 9, the parallel-session protocol and the coordinator's 10-02 windows):
#   - FIFO (v2, 2026-10-02 ~12:10Z, after [21dd4c] found newcomers jumping the queue): every job takes a ticket at WAIT
#     (an atomic counter under a mkdir lock) and may start only when no LIVE job holds a smaller ticket in work/queue/;
#   - at most SLOTS602 (default 2) gated jobs at once (mkdir locks, atomic; stale locks of dead pids are cleared);
#   - v3 (2026-10-02 ~14:20Z, [21dd4c]'s throughput report): up to SLOTS602_MAX (default 3) when memory sits idle. A slot
#     above SLOTS602 opens only if memory_pressure free >= 60 % (on top of every check below). It goes in ticket order among
#     v3 waiters (their queue file says "v3"), so pre-v3 waiters, which can only see slots 1..2, keep their priority there;
#   - CHAINS take one ticket per RUN (call gate602.sh around each step, as 612's chain612b/c/d do), never one ticket around
#     a whole chain script: one ticket holds its slot for the chain's whole duration and starves the FIFO;
#   - start only if free+inactive >= max(20, max_gb + 12) GiB AND memory_pressure free >= 55 % AND swap free >= 400 MiB
#     AND no production builder is running (def_pipeline.other_builders(), as recipe_dr524 checks);
#   - never inside the coordinator's reserved windows: 10-04 05:00Z .. 10-05 18:00Z, 10-06 06:00Z .. 18:00Z,
#     10-07 06:00Z .. 12:00Z (jobs already running are not killed; only new starts wait);
#   - v4 (2026-10-02 ~18:25Z, [c941b4]): DISK FLOOR. Start only if /System/Volumes/Data has >= the floor free: 70 GiB until
#     10-06 18:00Z (recipe_dr524 builds at 10-05 02:00Z and needs >= 35 GiB at its prep gate plus ~15 GB for its h5ad and
#     package), 60 GiB after that (AGENTS 10's final-phase floor). DISK602_FLOOR_GIB overrides both. A gated job writes work
#     files, so the floor is checked at every start; jobs already running are not touched;
#   - v5 (2026-10-02 ~19:30Z, [21dd4c]'s head-of-line report): BACKFILL. A waiter may start ahead of an earlier ticket only
#     when that ticket cannot fit memory now (its need > free+inactive) and this one can, and only while the skipped
#     ticket has waited < BACKFILL602_MAX_S (default 1800 s) since it joined: after that, strict FIFO again, so a big job
#     is never starved. Tickets record their need on a second line ("need=<GiB>"); a ticket without one (a pre-v5 waiter)
#     is never skipped. START lines name the tickets a job backfilled past;
#   - v6 (2026-10-02 ~21:55Z, [21dd4c]'s quiet-machine report): QUIET ADMISSION, opt-in per call with GATE602_QUIET=1. The
#     compressor can hold ~16 GiB of other apps' pages, so free+inactive sits near 20 GiB even on an idle machine and a
#     14 GB job (need 26) never starts. An opted-in job may start at free+inactive >= max(20, max_gb + 6) when
#     memory_pressure free >= 60 %, at most one gated job holds a slot, and no pre-v6 waiter is queued. Every start now goes
#     through one admission lock (work/admit.lock). A quiet start holds ALL other starts for QUIET602_HOLD_S (default 300 s)
#     while it ramps up (work/quiet_hold), so free+inactive catches up before anyone else starts. Opt in only with a
#     MEASURED peak (max_gb from a previous run's /usr/bin/time -l);
#   - logs WAIT / START / END with rc and peak RSS (/usr/bin/time -l) to experiments/602-council-gate/work/gate602.log.
# The script name deliberately matches none of the build gates' patterns (chain4/5.., run_arm_generic, def_day, run_ka_).
set -u
REPO=~/Code/virtual-cell-challenge-2026
W=${GATE602_W:-$REPO/experiments/602-council-gate/work}; mkdir -p "$W/queue"; LOG=$W/gate602.log   # GATE602_W: tests only
SLOTS=${SLOTS602:-2}; SLOTS_MAX=${SLOTS602_MAX:-3}; [ "$SLOTS_MAX" -lt "$SLOTS" ] && SLOTS_MAX=$SLOTS
[ $# -ge 4 ] && [ "$3" = "--" ] || { echo "usage: gate602.sh <label> <max_gb> -- <command ...>" >&2; exit 2; }
LABEL=$1; MAXGB=$2; shift 3
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
fi_gib(){ vm_stat | awk '/Pages free/{f=$3}/Pages inactive/{i=$3}END{gsub(/\./,"",f);gsub(/\./,"",i);printf "%d",(f+i)*16384/2^30}'; }
pfree(){ memory_pressure 2>/dev/null | awk -F': ' '/free percentage/{gsub(/%/,"",$2); print int($2)}'; }
# 670 (10-16): no swap file at all (total = 0.00M, e.g. right after a reboot) means nothing is swapped, not "swap exhausted"; pass it.
swapfree_mb(){ sysctl -n vm.swapusage | awk '{for(i=1;i<=NF;i++){ if($i=="total") t=$(i+2); if($i=="free") {v=$(i+2); u=substr(v,length(v)); n=substr(v,1,length(v)-1); if(u=="G") n*=1024}} } END{ if(t=="0.00M") print 999999; else print int(n) }'; }
disk_gib(){ df -k /System/Volumes/Data | awk 'NR==2{printf "%d", $4/1048576}'; }
disk_floor(){ if [ -n "${DISK602_FLOOR_GIB:-}" ]; then echo "$DISK602_FLOOR_GIB"
  elif [ "$(date -u +%Y%m%d%H%M)" -lt 202610061800 ]; then echo 70; else echo 60; fi; }
builder_running(){ (cd "$REPO" && uv run python -c 'import sys; sys.path.insert(0, "scripts"); import def_pipeline; raise SystemExit(0 if def_pipeline.other_builders() else 1)') 2>/dev/null && return 0
  # 670 (10-16): also yield while a PRODUCTION chain is alive, waiting or running. R-NEW (run690) sat at its own 26 GiB gate while 8 GB
  # jobs kept starting between builds. Deliberately narrow: watch691 / build7NN waiters live for days and must not block this gate.
  # 670 (10-16 09:20Z): also yields to a sealed reserve's access wrapper (R7a: r7a_access_1015b.sh), whose own gate is not FIFO with ours.
  ps -axo comm,args | awk '$1 ~ /^bash$/ && /(def_day|run_arm_generic|run_ka_527|run690|rabc694|r7a_access_[0-9a-z]+)\.sh/' | /usr/bin/grep -q .; }
blackout(){ local t; t=$(date -u +%Y%m%d%H%M)
  { [ "$t" -ge 202610040500 ] && [ "$t" -lt 202610051800 ]; } || { [ "$t" -ge 202610060600 ] && [ "$t" -lt 202610061800 ]; } \
  || { [ "$t" -ge 202610070600 ] && [ "$t" -lt 202610071200 ]; }; }
# --- FIFO ticket: atomic increment of work/ticket602 under a mkdir lock
until mkdir "$W/ticket.lock" 2>/dev/null; do sleep 0.2; done
T=$(( $(cat "$W/ticket602" 2>/dev/null || echo 0) + 1 )); echo "$T" > "$W/ticket602"; rmdir "$W/ticket.lock"
NEED=$((MAXGB + 12)); [ "$NEED" -lt 20 ] && NEED=20
NEEDQ=$((MAXGB + 6)); [ "$NEEDQ" -lt 20 ] && NEEDQ=20; QOK=${GATE602_QUIET:-0}      # v6 quiet admission (opt-in)
QF=$(printf '%s/queue/%08d.%s' "$W" "$T" "$$")
{ printf '%s v3\nneed=%s\n' "$LABEL" "$NEED"; [ "$QOK" = 1 ] && printf 'needq=%s\n' "$NEEDQ"; printf 'gate=v6\n'; } > "$QF"
skippable(){ local n q a   # v5 backfill: an earlier ticket that cannot fit memory now, while this one can (FI_NOW >= MYNEED)
  n=$(sed -n 's/^need=//p' "$1" 2>/dev/null); [ -n "$n" ] || return 1                 # pre-v5 ticket: never skipped
  q=$(sed -n 's/^needq=//p' "$1" 2>/dev/null); [ -n "$q" ] && [ "$QCOND" = 1 ] && n=$q  # v6: an opted-in ticket, machine quiet
  a=$(( $(date +%s) - $(stat -f %m "$1" 2>/dev/null || date +%s) ))
  [ "$a" -lt "${BACKFILL602_MAX_S:-1800}" ] && [ "$n" -gt "$FI_NOW" ] && [ "$MYNEED" -le "$FI_NOW" ]; }
first_in_line(){ local f t p; SKIPPED=""
  for f in "$W"/queue/*; do [ -e "$f" ] || continue; t=$(basename "$f"); p=${t#*.}; t=${t%%.*}
    if ! kill -0 "$p" 2>/dev/null; then rm -f "$f"; continue; fi        # dead waiter: drop its ticket
    if [ "$((10#$t))" -lt "$T" ]; then skippable "$f" && { SKIPPED="$SKIPPED $((10#$t))"; continue; }; return 1; fi; done; return 0; }
first_v3(){ local f t p; SKIPPED=""   # first among the waiters that can use a slot above SLOTS (v3 queue files)
  for f in "$W"/queue/*; do [ -e "$f" ] || continue; t=$(basename "$f"); p=${t#*.}; t=${t%%.*}
    if ! kill -0 "$p" 2>/dev/null; then rm -f "$f"; continue; fi
    if [ "$((10#$t))" -lt "$T" ] && grep -q " v3$" "$f" 2>/dev/null; then
      skippable "$f" && { SKIPPED="$SKIPPED $((10#$t))"; continue; }; return 1; fi; done; return 0; }
try_slots(){ local lo=$1 hi=$2 o; i=$lo; while [ $i -le "$hi" ]; do
    if mkdir "$W/slot$i.lock" 2>/dev/null; then SLOT=$i; echo "$$ $LABEL" > "$W/slot$i.lock/owner"; return 0; fi
    o=$(cut -d' ' -f1 "$W/slot$i.lock/owner" 2>/dev/null); if [ -n "$o" ] && ! kill -0 "$o" 2>/dev/null; then rm -rf "$W/slot$i.lock"; log "cleared stale slot$i (pid $o)"; continue; fi
    i=$((i + 1)); done; return 1; }
held_slots(){ local d o n=0; for d in "$W"/slot*.lock; do [ -d "$d" ] || continue
    o=$(cut -d' ' -f1 "$d/owner" 2>/dev/null); [ -n "$o" ] && kill -0 "$o" 2>/dev/null && n=$((n + 1)); done; echo "$n"; }
hold_active(){ local e; e=$(cat "$W/quiet_hold" 2>/dev/null); [ -n "$e" ] && [ "$(date +%s)" -lt "$e" ]; }
pre_v6_queued(){ local f p; for f in "$W"/queue/*; do [ -e "$f" ] || continue; p=$(basename "$f"); p=${p#*.}
    kill -0 "$p" 2>/dev/null && ! grep -q '^gate=v6$' "$f" 2>/dev/null && return 0; done; return 1; }
admit(){   # v6: every start goes through one admission lock, so a quiet start never races another start
  local o r=1; until mkdir "$W/admit.lock" 2>/dev/null; do
    o=$(cat "$W/admit.lock/owner" 2>/dev/null)
    if [ -n "$o" ] && ! kill -0 "$o" 2>/dev/null; then rm -rf "$W/admit.lock"
    elif [ -z "$o" ] && [ $(( $(date +%s) - $(stat -f %m "$W/admit.lock" 2>/dev/null || date +%s) )) -gt 10 ]; then rm -rf "$W/admit.lock"; fi
    sleep 0.2; done
  echo $$ > "$W/admit.lock/owner"
  if ! hold_active; then case "$1" in
    quiet)  if [ "$(held_slots)" -le 1 ] && ! pre_v6_queued && try_slots 1 "$SLOTS_MAX"; then
              echo $(( $(date +%s) + ${QUIET602_HOLD_S:-300} )) > "$W/quiet_hold"; QSTART=1; r=0; fi;;
    normal) if try_slots 1 "$SLOTS"; then r=0
            elif [ "$PF" -ge 60 ] && [ "$SLOTS_MAX" -gt "$SLOTS" ] && try_slots $((SLOTS + 1)) "$SLOTS_MAX"; then r=0; fi;;
    extra)  try_slots $((SLOTS + 1)) "$SLOTS_MAX" && r=0;;
  esac; fi
  rm -rf "$W/admit.lock"; return $r; }
SLOT=""; WAITED=0; SKIPPED=""; FI_NOW=0; PF=0; QCOND=0; MYNEED=$NEED; QSTART=0
trap 'rm -f "$QF"; [ "$(cat "$W/admit.lock/owner" 2>/dev/null)" = "$$" ] && rm -rf "$W/admit.lock"' EXIT
log "WAIT $LABEL ticket $T (max ${MAXGB} GB, needs free+inactive >= ${NEED} GiB$([ "$QOK" = 1 ] && echo ", quiet opt-in at ${NEEDQ}"), disk >= $(disk_floor) GiB) pid $$"
while :; do
  FI_NOW=$(fi_gib); PF=$(pfree); QCOND=0; MYNEED=$NEED
  [ "$PF" -ge 60 ] && [ "$(held_slots)" -le 1 ] && QCOND=1
  [ "$QOK" = 1 ] && [ "$QCOND" = 1 ] && MYNEED=$NEEDQ
  if ! blackout && ! hold_active && [ "$(disk_gib)" -ge "$(disk_floor)" ] && [ "$FI_NOW" -ge "$MYNEED" ] && [ "$(swapfree_mb)" -ge 400 ]; then
    if first_in_line && [ "$PF" -ge 55 ] && ! builder_running; then
      if [ "$FI_NOW" -lt "$NEED" ]; then admit quiet; else admit normal; fi
    elif [ "$PF" -ge 60 ] && [ "$SLOTS_MAX" -gt "$SLOTS" ] && [ "$FI_NOW" -ge "$NEED" ] && first_v3 && ! builder_running; then
      admit extra   # an idle-memory slot, in ticket order among v3 waiters
    fi
    [ -n "$SLOT" ] && break
  fi
  [ $((WAITED % 600)) -eq 0 ] && [ $WAITED -gt 0 ] && log "still waiting $LABEL (ticket $T): first_in_line $(first_in_line && echo yes || echo no), fi $(fi_gib) GiB, pressure $(pfree) %, swapfree $(swapfree_mb) MB, disk $(disk_gib)/$(disk_floor) GiB, quiet opt-in/cond $QOK/$QCOND, hold $(hold_active && echo on || echo off), blackout $(blackout && echo yes || echo no)"
  sleep 10; WAITED=$((WAITED + 10))
done
rm -f "$QF"
trap 'rm -rf "$W/slot$SLOT.lock"' EXIT
QNOTE=""; [ "$QSTART" = 1 ] && QNOTE=" QUIET start (fi $FI_NOW < need $NEED, needq $NEEDQ; other starts held until $(date -u -r "$(cat "$W/quiet_hold")" +%TZ))"
log "START $LABEL ticket $T slot $SLOT (fi $(fi_gib) GiB, pressure $(pfree) %, disk $(disk_gib) GiB)${SKIPPED:+ backfilled past ticket(s)$SKIPPED}$QNOTE: $*"
TF=$W/time.$$.txt
/usr/bin/time -l "$@" 2> >(tee "$TF" >&2); RC=$?
# the job may print its own /usr/bin/time block too (a nested time): keep the largest value, the gate's outer one
PEAK=$(awk '/maximum resident set size/{if ($1 > v) v = $1} END{if (v != "") printf "%.2f", v/1e9}' "$TF" 2>/dev/null); rm -f "$TF"
log "END $LABEL rc=$RC peakRSS=${PEAK:-?} GB"
exit $RC
