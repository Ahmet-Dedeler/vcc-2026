#!/bin/bash
# Coordinator check before apply_october538.py. Match interpreter + an argv word that
# IS one of the nine files or recipe_dr524.sh. A mention buried in a different argument
# does not count. Do not use pgrep -f on a bare name: this script's own argv is
# procs_clear550.sh, so it does not match itself.
set -uo pipefail
hits="$(ps -axo pid=,comm=,args= | awk '
  {
    comm = $2
    if (comm !~ /^(bash|python[0-9.]*|uv)$/) next
    for (i = 3; i <= NF; i++) {
      if ($i ~ /(^|\/)(recipe_dr524\.sh|def_day\.sh|def_day_helpers\.py|relevel_iter\.py|def_sameline\.py|october_recipe\.py|october-decisions\.json|gate524\.py|def_postbuild\.py|def_pipeline\.py)$/) {
        print
        hit = 1
        break
      }
    }
  }
  END { exit (hit ? 0 : 1) }
')"
rc=$?
if [ "$rc" -eq 0 ]; then
  echo "BUSY: an interpreter is still executing recipe_dr524.sh or one of the nine files"
  echo "$hits"
  exit 2
fi
echo "CLEAR: no bash/python/uv argv word is recipe_dr524.sh or one of the nine files"
exit 0
