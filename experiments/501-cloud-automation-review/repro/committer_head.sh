#!/usr/bin/env bash
# The loop at scripts/commit_after_chain.sh:45-47 exits 0 on the FIRST file that
# equals HEAD, before git commit, even when a later file differs.
set -u
D=$(mktemp -d)
cd "$D" || exit 1
git init -q
git config user.email t@example.com
git config user.name t
echo one > a.txt
echo one > b.txt
git add a.txt b.txt
git commit -q -m init
echo TWO > b.txt   # a.txt still equals HEAD; b.txt does not
FILES=(a.txt b.txt)
action=commit
for p in "${FILES[@]}"; do
  [ "$(shasum -a 256 "$p" | cut -d' ' -f1)" != "$(git show HEAD:"$p" | shasum -a 256 | cut -d' ' -f1)" ] \
    || { echo "EXIT 0 at $p (equals HEAD) — commit not reached"; action=skipped; break; }
done
if [ "$action" = commit ]; then
  git commit -q --only -m would -- "${FILES[@]}" && echo "committed"
fi
echo "HEAD after the loop: $(git log --oneline)"
echo "worktree b.txt=$(cat b.txt)  (uncommitted if the loop exited early)"
echo "status: $(git status --short)"
rm -rf "$D"
