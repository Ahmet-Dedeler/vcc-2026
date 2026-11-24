#!/usr/bin/env bash
# Hop the -m strings from chain1002/1003/1004's name() through the arm's "$NAME"
# and into a fake `vcc submit -m`. Non-interactive bash, same quoting as the chains.
set -u
fake=$(mktemp)
cat > "$fake" <<'EOS'
#!/usr/bin/env bash
NAME="$5"; NOTE="$6"
printf 'NAME=%s\n' "$NAME"
printf 'NOTE=%s\n' "$NOTE"
# what execve would see for: vcc submit pkg -m "$NAME" --wait
printf 'ARGV'
printf ' <%s>' vcc submit pkg.vcc -m "$NAME" --wait
printf '\n'
EOS
chmod +x "$fake"
run() { # $1 chain file, $2 tag
  local src="$1" tag="$2"
  # pull the echo argument for this tag out of name()'s case, then pass it the way the chain does
  local msg
  msg=$(awk -v tag="$tag" '$0 ~ "^  " tag "\\)" { if (match($0, /echo ".*"/)) { s=substr($0, RSTART+6, RLENGTH-7); print s } }' "$src")
  # arm1002.sh: $1 slot, $2 tag, $3 log, $4 pkg, $5 name, $6 note
  bash "$fake" 1 "$tag" run.log pkg.vcc "$msg" "$tag slot" 
  echo "---"
}
echo "=== 479sp2 apostrophe + parentheses ==="
run /workspace/experiments/479-pertarget-sharpen/chain1002.sh 479sp2
echo "=== 495 em dash + !! ==="
run /workspace/experiments/479-pertarget-sharpen/chain1002.sh 495
echo "=== 471c apostrophe + parentheses ==="
run /workspace/experiments/479-pertarget-sharpen/chain1002.sh 471c
echo "=== 498sp2 comma + apostrophe ==="
run /workspace/experiments/498-combined-bulk/chain1004.sh 498sp2
echo "=== history expansion is off in this non-interactive bash (histexpand) ==="
set -o | awk '$1=="histexpand"{print}'
# prove !! is not expanded by running the chain's own echo via bash -c
bash -c 'echo "taking action is the best action—anytime. dont let it stay in your mind, go grind!!"'
rm -f "$fake"
