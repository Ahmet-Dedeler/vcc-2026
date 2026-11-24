#!/usr/bin/env bash
# The slot-1 wait in arm1002.sh:20-23, arm1003.sh:20-23, arm1004.sh:20-23.
# The noon cap sits inside an until-loop whose exit condition is "clock >= 00:01",
# which is nine hours before noon. The body cannot observe noon, and a late start
# skips the body entirely.
set -u
loop() { # now slot cap  -> prints what the real loop does (no sleep)
  local now="$1" slot="$2" cap="$3" body=0 capped=0
  if [ "$now" -ge "$slot" ]; then
    printf 'now=%s slot=%s cap=%s  until-already-true  cap_checked=no  WOULD_SUBMIT=yes\n' "$now" "$slot" "$cap"
    return
  fi
  # one iteration of the body at this instant
  body=1
  if [ "$now" -ge "$cap" ]; then capped=1; fi
  printf 'now=%s slot=%s cap=%s  body=%s cap_would_stop=%s  (loop continues; cap is %s, exit is %s)\n' \
    "$now" "$slot" "$cap" "$body" "$capped" "$cap" "$slot"
}
echo "=== on-time arming (10-01 20:00 through 10-02 00:01), cap 10-02 12:00 ==="
loop 202610012000 202610020001 202610021200
loop 202610012359 202610020001 202610021200
loop 202610020001 202610020001 202610021200
echo "=== clock jump / late wake: past noon, and the next day ==="
loop 202610021500 202610020001 202610021200
loop 202610030900 202610020001 202610021200
echo "=== same shape for 10-03 and 10-04 ==="
loop 202610031500 202610030001 202610031200
loop 202610051200 202610040001 202610041200
echo "=== slot 2's cap IS reachable: its until waits for a log line, not for a time ==="
# faithful fragment from arm1002.sh:25-29, with a fixture log and a frozen clock
D=$(mktemp -d)
printf '2026-10-02T00:05:00Z ARMED slot 1: 479p2g\n2026-10-02T00:20:00Z SUBMITTED ok\n2026-10-02T00:21:00Z WARN capture failed\n' > "$D/s1.log"
now=202610021500
if ! grep -q "CAPTURED own entry" "$D/s1.log"; then
  if grep -q "NOT SUBMITTING" "$D/s1.log"; then echo "slot2: slot1 refused"; 
  elif [ "$now" -ge 202610021200 ]; then echo "slot2: CAPTURED absent and now=$now >= noon -> STOP (cap fires)"; fi
fi
rm -rf "$D"
