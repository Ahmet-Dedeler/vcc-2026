#!/bin/bash
# test_gate602.sh -- behaviour tests for gate602.sh (v4 disk floor, v5 backfill) in an isolated GATE602_W, with PATH shims
# for vm_stat / memory_pressure / sysctl / df / uv and the minute clock, so memory, pressure, swap, disk, builders and
# blackouts are deterministic. Nothing here touches the live queue in experiments/602-council-gate/work.
#   usage: bash experiments/602-council-gate/test_gate602.sh [path/to/gate602.sh]     (~5 min; exit 0 = all PASS)
#   T1 backfill: a small job starts ahead of an earlier big ticket that cannot fit; the big one starts once memory frees
#   T2 starvation bound: once the big ticket has waited >= BACKFILL602_MAX_S, nobody skips it (strict FIFO)
#   T3 a pre-v5 ticket (no need= line) is never skipped; a dead waiter's ticket is dropped
#   T4 FIFO among tickets that fit: with the only slot held, C (earlier) starts before D, and D does not "backfill"
#   T5 the v4 disk floor: no start below the floor, start once disk frees
#   T6-T11 v6 quiet admission: opted-in start on a quiet machine + hold; no opt-in, no start; the hold blocks every other
#          start until it expires; none with 2 slots held, at pressure < 60 %, or while a pre-v6 waiter is queued
set -u
G=${1:-~/Code/virtual-cell-challenge-2026/experiments/602-council-gate/gate602.sh}
TD=$(mktemp -d -t g602test); SHIM=$TD/shim; mkdir -p "$SHIM/bin"
PIDS=""; PASS=0; FAIL=0
cleanup(){ for p in $PIDS; do kill "$p" 2>/dev/null; done; pkill -f "$TD/" 2>/dev/null; rm -rf "$TD"; }
trap cleanup EXIT
ok(){ if [ "$1" = 0 ]; then PASS=$((PASS + 1)); echo "PASS $2"; else FAIL=$((FAIL + 1)); echo "FAIL $2"; fi; }
cat > "$SHIM/bin/vm_stat" <<'S'
#!/bin/bash
p=$(( $(cat "$SHIM/fi_gib") * 65536 ))
printf 'Mach Virtual Memory Statistics: (page size of 16384 bytes)\nPages free:     %s.\nPages inactive:     0.\n' "$p"
S
cat > "$SHIM/bin/memory_pressure" <<'S'
#!/bin/bash
echo "System-wide memory free percentage: $(cat "$SHIM/pf")%"
S
cat > "$SHIM/bin/sysctl" <<'S'
#!/bin/bash
[ "$*" = "-n vm.swapusage" ] && { echo "total = 2048.00M  used = 0.00M  free = 2048.00M  (encrypted)"; exit 0; }
exec /usr/sbin/sysctl "$@"
S
cat > "$SHIM/bin/df" <<'S'
#!/bin/bash
echo "Filesystem 1024-blocks Used Available Capacity iused ifree %iused Mounted on"
echo "/dev/fake 4000000000 1 $(( $(cat "$SHIM/disk_gib") * 1048576 )) 1% 1 1 1% /System/Volumes/Data"
S
printf '#!/bin/bash\nexit 1\n' > "$SHIM/bin/uv"                       # builder_running(): no production builder
cat > "$SHIM/bin/date" <<'S'
#!/bin/bash
for a in "$@"; do [ "$a" = "+%Y%m%d%H%M" ] && { echo 202610021930; exit 0; }; done   # outside every blackout
exec /bin/date "$@"
S
chmod +x "$SHIM"/bin/*
export SHIM PATH="$SHIM/bin:$PATH"
set_env(){ echo "$1" > "$SHIM/fi_gib"; echo "$2" > "$SHIM/pf"; echo "$3" > "$SHIM/disk_gib"; }
gate(){ local w=$1 label=$2 gb=$3 dur=$4; shift 4
  env "$@" GATE602_W="$w" bash "$G" "$label" "$gb" -- sleep "$dur" > /dev/null 2>&1 & PIDS="$PIDS $!"; }
until_log(){ local w=$1 pat=$2 cap=$3 t=0; while [ $t -lt "$cap" ]; do /usr/bin/grep -qE "$pat" "$w/gate602.log" 2>/dev/null && return 0; sleep 1; t=$((t + 1)); done; return 1; }
line_no(){ /usr/bin/grep -nE "$2" "$1/gate602.log" | head -1 | cut -d: -f1; }

# T1 backfill
W=$TD/w1; mkdir -p "$W"; set_env 22 70 500
gate "$W" big1 14 3; sleep 2; gate "$W" small1 2 2
until_log "$W" "START small1 .*backfilled past ticket\(s\) 1:" 30; ok $? "T1 small job backfills past the big ticket (fi 22 < need 26)"
/usr/bin/grep -q "START big1" "$W/gate602.log"; [ $? -ne 0 ]; ok $? "T1 the big ticket has not started while it cannot fit"
set_env 30 70 500; until_log "$W" "START big1" 30; ok $? "T1 the big ticket starts once fi reaches its need"
# T2 starvation bound
W=$TD/w2; mkdir -p "$W"; set_env 22 70 500
gate "$W" big2 14 3; sleep 4; gate "$W" small2 2 2 BACKFILL602_MAX_S=2
until_log "$W" "START small2" 20; [ $? -ne 0 ]; ok $? "T2 no backfill past a ticket older than BACKFILL602_MAX_S"
set_env 30 70 500; until_log "$W" "START small2" 40
b=$(line_no "$W" "START big2"); s=$(line_no "$W" "START small2"); [ -n "$b" ] && [ -n "$s" ] && [ "$b" -lt "$s" ]
ok $? "T2 strict FIFO: the big ticket starts first once it fits, then the small one"
# T3 pre-v5 ticket never skipped; a dead waiter is dropped
W=$TD/w3; mkdir -p "$W/queue"; set_env 22 70 500
sleep 300 & OLD=$!; PIDS="$PIDS $OLD"; echo 1 > "$W/ticket602"; echo "old v3" > "$W/queue/00000001.$OLD"
gate "$W" small3 2 1
until_log "$W" "START small3" 20; [ $? -ne 0 ]; ok $? "T3 a pre-v5 ticket (no need= line) is never skipped"
kill "$OLD" 2>/dev/null; until_log "$W" "START small3" 30; ok $? "T3 the dead waiter's ticket is dropped and the job starts"
/usr/bin/grep -q "START small3 .*backfilled" "$W/gate602.log"; [ $? -ne 0 ]; ok $? "T3 that start is not logged as a backfill"
# T4 FIFO among tickets that fit (one slot, held)
W=$TD/w4; mkdir -p "$W/slot1.lock"; set_env 30 70 500
sleep 300 & HOLD=$!; PIDS="$PIDS $HOLD"; echo "$HOLD holder" > "$W/slot1.lock/owner"
gate "$W" c4 2 2 SLOTS602=1 SLOTS602_MAX=1; sleep 2; gate "$W" d4 2 1 SLOTS602=1 SLOTS602_MAX=1
sleep 15; /usr/bin/grep -qE "START (c4|d4)" "$W/gate602.log"; [ $? -ne 0 ]; ok $? "T4 nothing starts while the only slot is held"
kill "$HOLD" 2>/dev/null; until_log "$W" "START d4" 60
c=$(line_no "$W" "START c4"); d=$(line_no "$W" "START d4"); [ -n "$c" ] && [ -n "$d" ] && [ "$c" -lt "$d" ]; ok $? "T4 C (earlier ticket) starts before D"
/usr/bin/grep -q "START d4 .*backfilled" "$W/gate602.log"; [ $? -ne 0 ]; ok $? "T4 D never backfills past a ticket that fits"
# T5 v4 disk floor
W=$TD/w5; mkdir -p "$W"; set_env 30 70 10
gate "$W" disk5 2 1
until_log "$W" "START disk5" 15; [ $? -ne 0 ]; ok $? "T5 no start below the disk floor (10 < 70 GiB)"
set_env 30 70 500; until_log "$W" "START disk5" 30; ok $? "T5 starts once disk is above the floor"
# v6 quiet admission (opt-in GATE602_QUIET=1): needq = max(20, max_gb + 6) when pressure >= 60 % and <= 1 slot held
hold_slot(){ local w=$1 i=$2; mkdir -p "$w/slot$i.lock"; sleep 300 & PIDS="$PIDS $!"; echo "$! holder" > "$w/slot$i.lock/owner"; }
# T6 quiet admission
W=$TD/w6; mkdir -p "$W"; set_env 21 65 500; hold_slot "$W" 1
gate "$W" big6 14 20 GATE602_QUIET=1
until_log "$W" "START big6 .*QUIET start \(fi 21 < need 26, needq 20" 25; ok $? "T6 an opted-in 14 GB job starts on a quiet machine (fi 21, 1 slot held, pressure 65 %)"
e=$(cat "$W/quiet_hold" 2>/dev/null); [ -n "$e" ] && [ "$e" -gt "$(/bin/date +%s)" ]; ok $? "T6 the quiet start set a hold"
# T7 opt-in required
W=$TD/w7; mkdir -p "$W"; set_env 21 65 500; hold_slot "$W" 1
gate "$W" big7 14 5
until_log "$W" "START big7" 15; [ $? -ne 0 ]; ok $? "T7 without GATE602_QUIET the same job does not start"
# T8 the hold blocks every other start until it expires
W=$TD/w8; mkdir -p "$W"; set_env 21 65 500; hold_slot "$W" 1
gate "$W" big8 14 60 GATE602_QUIET=1 QUIET602_HOLD_S=20
until_log "$W" "START big8 .*QUIET" 25; ok $? "T8 the quiet job starts"
gate "$W" small8 2 1
until_log "$W" "START small8" 10; [ $? -ne 0 ]; ok $? "T8 a small job does not start during the hold"
until_log "$W" "START small8" 40; ok $? "T8 the small job starts after the hold expires"
h=$(cat "$W/quiet_hold"); st=$(/usr/bin/grep "START small8" "$W/gate602.log" | head -1 | cut -c1-19)
ss=$(/bin/date -j -u -f "%Y-%m-%dT%H:%M:%S" "$st" +%s 2>/dev/null); [ -n "$ss" ] && [ "$ss" -ge "$h" ]; ok $? "T8 ... and not before the hold's expiry"
# T9 busy machine: two slots held
W=$TD/w9; mkdir -p "$W"; set_env 21 65 500; hold_slot "$W" 1; hold_slot "$W" 2
gate "$W" big9 14 5 GATE602_QUIET=1
until_log "$W" "START big9" 15; [ $? -ne 0 ]; ok $? "T9 no quiet start while 2 slots are held"
# T10 pressure below 60 %
W=$TD/w10; mkdir -p "$W"; set_env 21 58 500; hold_slot "$W" 1
gate "$W" big10 14 5 GATE602_QUIET=1
until_log "$W" "START big10" 15; [ $? -ne 0 ]; ok $? "T10 no quiet start at pressure 58 %"
# T11 a pre-v6 waiter queued (it cannot see the hold): no quiet start until it is gone
W=$TD/w11; mkdir -p "$W/queue"; set_env 21 65 500; hold_slot "$W" 1
sleep 300 & OLD=$!; PIDS="$PIDS $OLD"; printf 'old v3\nneed=20\n' > "$W/queue/00000099.$OLD"
gate "$W" big11 14 5 GATE602_QUIET=1
until_log "$W" "START big11" 15; [ $? -ne 0 ]; ok $? "T11 no quiet start while a pre-v6 waiter is queued"
kill "$OLD" 2>/dev/null; until_log "$W" "START big11 .*QUIET" 30; ok $? "T11 the quiet start happens once that waiter is gone"
echo; echo "$PASS PASS, $FAIL FAIL"; [ "$FAIL" -eq 0 ]
