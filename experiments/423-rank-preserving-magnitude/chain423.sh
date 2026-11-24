#!/bin/bash
# wait for a driver log to say "queue done", then run the next queue (serial H1 arms)
cd "$(dirname "$0")" || exit 1
while ! grep -q "queue done\|FAILED" "$1"; do sleep 30; done
grep -q "FAILED" "$1" && { echo "previous queue failed; not chaining"; exit 1; }
exec bash drive423.sh "$2"
