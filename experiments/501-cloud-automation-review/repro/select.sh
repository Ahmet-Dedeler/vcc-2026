#!/usr/bin/env bash
# chain1002/1003 pick the first two distinct READY tags. chain1004 skips s1's tag in s2
# and promotes s2 to s1 if s1's list produced nothing. Fixture: a ready-set, no packages.
set -u
pick12() { # ready flags in order, one per candidate name passed after --
  local A1="" A2="" c
  for c in "$@"; do
    # name:0 or name:1
    local n="${c%%:*}" ok="${c#*:}"
    if [ "$ok" = 1 ]; then [ -z "$A1" ] && A1=$n || { [ -z "$A2" ] && A2=$n; }; fi
  done
  echo "s1=${A1:-<none>} s2=${A2:-<none>} equal=$([ "$A1" = "$A2" ] && [ -n "$A1" ] && echo YES || echo no)"
}
echo "chain1002 shape (first two ready):"
pick12 479p2g:1 479sp2:1 495:1 471c:1
pick12 479p2g:0 479sp2:0 495:1 471c:1
pick12 479p2g:0 479sp2:1 495:0 471c:0
pick1004() { # $1 = L1 words, $2 = L2 words, then ready names
  local L1="$1" L2="$2" A1="" A2="" c ready="$3"
  isready() { case " $ready " in *" $1 "*) return 0;; *) return 1;; esac; }
  for c in $L1; do [ "$c" = "-" ] && continue; if isready "$c"; then A1=$c; break; fi; done
  for c in $L2; do [ "$c" = "-" ] && continue; [ "$c" = "$A1" ] && continue; if isready "$c"; then A2=$c; break; fi; done
  [ -z "$A1" ] && [ -n "$A2" ] && { A1=$A2; A2=""; }
  echo "L1=[$L1] L2=[$L2] ready=[$ready] -> s1=${A1:-<none>} s2=${A2:-<none>}"
}
echo "chain1004 shape:"
pick1004 "498p2g 460dr 461dr" "497dr 460dr 461dr" "498p2g 497dr"
pick1004 "498p2g 460dr" "497dr 460dr" "460dr"
pick1004 "498sp2 460dr" "498sp2 460dr" "498sp2 460dr"
pick1004 "460dr 461dr" "460dr 461dr" "461dr"
