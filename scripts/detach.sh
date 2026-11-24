#!/usr/bin/env bash
# Run a command in a NEW SESSION (setsid) so that quitting or restarting the launching app (Claude Desktop), which signals its
# whole process group, cannot kill it. nohup alone only ignores SIGHUP. 2026-10-15: every nohup'd waiter, caffeinate and snapshot
# loop launched from an agent shell at 13:33Z was gone by 21:19Z with no reboot (uptime 8 h).
#   nohup bash scripts/detach.sh <cmd> [args...] > <log> 2>&1 &
exec python3 -c 'import os, sys; os.setsid(); os.execvp(sys.argv[1], sys.argv[1:])' "$@"
