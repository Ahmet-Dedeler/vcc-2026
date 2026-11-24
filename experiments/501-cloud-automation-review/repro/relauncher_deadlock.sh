#!/usr/bin/env bash
# Two facts, both without the real chains:
#  1. chain495's start gate matches a live chain479b.sh, and chain479b's gate matches a live chain495.sh.
#  2. relaunch479b.sh treats "chain495 is not running" as "past the gate" and launches chain479b.
set -u
echo "=== awk self/cross match on synthetic ps lines (the gates' exact regexes) ==="
line495='bash bash experiments/495-lodo-cap6/chain495.sh'
line479b='bash bash experiments/479-pertarget-sharpen/chain479b.sh'
line479='bash bash experiments/479-pertarget-sharpen/chain479.sh'
linerel='bash bash experiments/479-pertarget-sharpen/relaunch479b.sh'
linecommit='bash bash scripts/commit_after_chain.sh experiments/495-lodo-cap6/commit495.spec'
linearm='bash bash experiments/479-pertarget-sharpen/arm1002.sh 1 479p2g'
gate495='chain479b?\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4'
gate479b='chain479\.sh|chain495\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4'
match() { printf '%s\n' "$1" | awk -v re="$2" '$1 ~ /^bash$/ && $0 ~ re { found=1 } END { exit !found }'; }
for pair in "495-gate|479b-proc" "495-gate|495-proc" "495-gate|relauncher" "495-gate|commit495" \
            "479b-gate|495-proc" "479b-gate|479b-proc" "479b-gate|479-proc" "479b-gate|arm1002"; do
  gate="${pair%%|*}"; proc="${pair#*|}"
  case "$gate" in 495-gate) re="$gate495";; *) re="$gate479b";; esac
  case "$proc" in
    479b-proc) ln="$line479b";; 495-proc) ln="$line495";; 479-proc) ln="$line479";;
    relauncher) ln="$linerel";; commit495) ln="$linecommit";; arm1002) ln="$linearm";; esac
  if match "$ln" "$re"; then bit=MATCH; else bit=no; fi
  printf '  %-10s vs %-12s %s\n' "$gate" "$proc" "$bit"
done

echo "=== relauncher condition: 'not running' is success (relaunch479b.sh:9-10) ==="
# until grep patched || ! ps chain495; do ...; done
# Evaluate the until-predicate for three states.
pred() { # $1=patched_in_log (0/1)  $2=process_running (0/1)
  local patched="$1" running="$2"
  if [ "$patched" = 1 ] || [ "$running" = 0 ]; then echo START_479b; else echo KEEP_WAITING; fi
}
echo "  495 running, no patched yet : $(pred 0 1)   (correct)"
echo "  495 running, patched        : $(pred 1 1)   (correct)"
echo "  495 not started, no patched : $(pred 0 0)   (BUG: launches 479b)"
echo "  495 not started, STALE patched in the appended log : $(pred 1 0)   (BUG: launches 479b)"

echo "=== live circular wait (fake scripts, killed after 3s) ==="
D=$(mktemp -d)
cat > "$D/chain495.sh" <<'EOS'
#!/usr/bin/env bash
echo "495 gate: waiting until no chain479b.sh"
while ps -axo args= | awk '/[c]hain479b\.sh/ {found=1} END{exit !found}'; do sleep 0.2; done
echo "495 RELEASED"
EOS
cat > "$D/chain479b.sh" <<'EOS'
#!/usr/bin/env bash
echo "479b gate: waiting until no chain495.sh"
while ps -axo args= | awk '/[c]hain495\.sh/ {found=1} END{exit !found}'; do sleep 0.2; done
echo "479b RELEASED"
EOS
chmod +x "$D/chain495.sh" "$D/chain479b.sh"
bash "$D/chain495.sh" > "$D/495.out" 2>&1 & p1=$!
bash "$D/chain479b.sh" > "$D/479b.out" 2>&1 & p2=$!
sleep 3
echo "  after 3s, 495 pid $p1 still $(kill -0 $p1 2>/dev/null && echo RUNNING || echo dead)"
echo "  after 3s, 479b pid $p2 still $(kill -0 $p2 2>/dev/null && echo RUNNING || echo dead)"
echo "  495 log: $(cat "$D/495.out")"
echo "  479b log: $(cat "$D/479b.out")"
kill "$p1" "$p2" 2>/dev/null
wait "$p1" "$p2" 2>/dev/null
rm -rf "$D"
