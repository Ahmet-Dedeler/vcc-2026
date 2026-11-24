#!/usr/bin/env bash
# wait for pid $1 (capped at $2 hours), then exec the command in the remaining args.
pid=$1; cap_h=$2; shift 2
deadline=$(( $(date +%s) + cap_h * 3600 ))
while kill -0 "$pid" 2>/dev/null; do [ "$(date +%s)" -ge "$deadline" ] && exit 3; sleep 20; done
exec "$@"
