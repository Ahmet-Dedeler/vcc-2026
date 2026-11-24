#!/usr/bin/env bash
# wait for pid $1 to exit (polling kill -0; the wait is capped at $2 hours), then run chain438.sh with the rest.
pid=$1; cap_h=$2; shift 2
deadline=$(( $(date +%s) + cap_h * 3600 ))
while kill -0 "$pid" 2>/dev/null; do [ "$(date +%s)" -ge "$deadline" ] && { echo "$(date -u +%FT%TZ) after438: cap reached waiting for $pid" >> "$(dirname "$0")/work/chain438.log"; exit 3; }; sleep 20; done
exec "$(dirname "$0")/${CHAIN:-chain438.sh}" "$@"
