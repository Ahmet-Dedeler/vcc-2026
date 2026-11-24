#!/usr/bin/env bash
# The replacement predicate from patches/f4-relauncher.patch, without ps:
# start only after this relauncher has seen chain495's pid alive, and that pid
# has written past-gate or has since exited.
set -u
decide() { # $1=seen pid (empty if never observed)  $2=past-gate contents  $3=1 if that pid is still alive
  local SEEN="$1" GATE="$2" ALIVE="$3"
  if [ -n "$SEEN" ] && [ "$GATE" = "$SEEN" ]; then echo START; return; fi
  if [ -n "$SEEN" ] && [ "$ALIVE" = 0 ]; then echo START; return; fi
  echo WAIT
}
printf 'not started (no pid, stale gate ignored) : %s\n' "$(decide '' oldpid 0)"
printf 'running, gate not written yet            : %s\n' "$(decide 123 '' 1)"
printf 'running, past-gate is this pid           : %s\n' "$(decide 123 123 1)"
printf 'seen, then exited before patch (back)    : %s\n' "$(decide 123 '' 0)"
printf 'stale gate, a new pid is what we saw     : %s\n' "$(decide 200 111 1)"
printf 'never saw it, gate equals some old pid   : %s\n' "$(decide '' 111 0)"
