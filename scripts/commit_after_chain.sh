#!/usr/bin/env bash
# Commit a chain's production edits once the chain has exited with every required success line in its log, so the next
# chain (which waits for "file == HEAD") does not depend on a live coordinator session (sessions die; AGENTS §6.2).
#   bash scripts/commit_after_chain.sh <spec file>
# The spec file keeps every path OFF this process's command line: chains wait for "no bash process whose args contain
# run_arm_generic.sh / def_day.sh", and a committer with those paths in its argv deadlocked them (2026-09-27, AGENTS §6.2).
# Spec lines (KEY=value; REQ and FILE repeat):
#   PIDF=<chain.pid>  CLOG=<chain.log>  DEADLINE=<YYYYmmddHHMM>  MSGF=<message file>
#   EXPECT=<path=sha16 path=sha16 ...>   (optional: each file must have exactly that sha prefix)
#   REQ=<extended regex that must match a line of CLOG>
#   FILE=<path to commit>
# Waits for the chain's PROCESS to exit (fire on exit, then read the log). Refuses if the log has "RESTORED", if any REQ
# is missing, or if an EXPECT sha differs. Commits only the FILE paths (git commit --only), retrying while another process
# holds .git/index.lock. The message file's text gets each matched REQ line appended.
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
SPEC="${1:?spec}"
PIDF=""; CLOG=""; DEADLINE=""; MSGF=""; EXPECT=""; REQ=(); FILES=()
while IFS= read -r line || [ -n "$line" ]; do
  case "$line" in ''|\#*) continue;; esac
  k="${line%%=*}"; v="${line#*=}"
  case "$k" in
    PIDF) PIDF="$v";; CLOG) CLOG="$v";; DEADLINE) DEADLINE="$v";; MSGF) MSGF="$v";; EXPECT) EXPECT="$v";;
    REQ) REQ+=("$v");; FILE) FILES+=("$v");;
    *) echo "unknown spec key $k"; exit 2;;
  esac
done < "$SPEC"
[ -n "$PIDF" ] && [ -n "$CLOG" ] && [ -n "$DEADLINE" ] && [ -n "$MSGF" ] && [ ${#REQ[@]} -gt 0 ] && [ ${#FILES[@]} -gt 0 ] \
  || { echo "incomplete spec $SPEC"; exit 2; }
OUT="${CLOG%.log}.commit.log"
say(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$OUT"; }
NAME="$(basename "${PIDF%.pid}")"
say "armed (pid $$) from spec $SPEC: waiting for $NAME (pid file $PIDF) to exit; files ${FILES[*]}"
# the pid file is re-read every minute, so a chain that is killed and re-armed (new pid) is followed, not skipped
while P="$(cat "$PIDF" 2>/dev/null)"; ps -p "${P:-1}" -o args= 2>/dev/null | /usr/bin/grep -q "$NAME"; do
  [ "$(date -u +%Y%m%d%H%M)" -ge "$DEADLINE" ] && { say "chain still running at the deadline -- not committing"; exit 1; }
  sleep 60; done
sleep 5
/usr/bin/grep -q "RESTORED" "$CLOG" && { say "the chain restored its files -- not committing"; exit 1; }
MSG="$(cat "$MSGF")"
for r in "${REQ[@]}"; do
  L="$(/usr/bin/grep -E "$r" "$CLOG" | tail -1)"
  [ -n "$L" ] || { say "required line missing: $r -- not committing"; exit 1; }
  MSG="$MSG
$L"; done
NDIFF=0
for p in "${FILES[@]}"; do
  [ "$(shasum -a 256 "$p" | cut -d' ' -f1)" != "$(git show HEAD:"$p" 2>/dev/null | shasum -a 256 | cut -d' ' -f1)" ] && NDIFF=$((NDIFF + 1))
done
[ "$NDIFF" -gt 0 ] || { say "every file equals HEAD -- nothing to commit (already committed?)"; exit 0; }
for e in $EXPECT; do
  f="${e%%=*}"; want="${e#*=}"; got="$(shasum -a 256 "$f" | cut -c1-16)"
  [ "$got" = "$want" ] || { say "$f is $got, expected $want -- not committing"; exit 1; }; done
MSG="$MSG

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
for i in $(seq 1 20); do
  [ -e .git/index.lock ] && { sleep 30; continue; }
  if git commit -q --only -m "$MSG" -- "${FILES[@]}" >> "$OUT" 2>&1; then say "COMMITTED $(git log --oneline -1 | cut -c1-120)"; exit 0; fi
  sleep 30; done
say "commit failed after 20 attempts"; exit 1
