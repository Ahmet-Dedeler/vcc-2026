#!/bin/bash
# 591 — stubbed-reader tests of the patched run_arm_generic.sh gate (never touches the live script).
# Runs under macOS /bin/bash 3.2 with BSD awk/grep. Every memory/process/disk reader the gate calls (memory_pressure,
# vm_stat, sysctl, uv, ps, df) is a stub on PATH driven by STUB_* variables. The gate functions are cut out of the
# ORIGINAL and the PATCHED script text (from `log(){` up to `restore(){`) and evaluated in a subshell per grid point.
#   R1  build mode: patched exit code == original's at every grid point (identical behaviour)
#   R2  prep mode: patched == original OR the pressure branch (pct >= 55 AND swap >= 400 MiB), all else unchanged;
#       in particular identical wherever free+inactive >= 23 and pressure >= 40
#   R3  the gate line keeps the "gate: free " prefix and names the mode and the expected branch
#   R4  call sites: "wait_gate build" precedes invocation 1, "wait_gate prep" precedes vcc prep; bash -n passes
#   R5  live readers (no stubs): the patched gate parses this Mac's real sysctl/vm_stat/memory_pressure output
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"; REPO_ROOT="$(cd "$HERE/../.." && pwd)"
TARGET="${TARGET591:-$REPO_ROOT/experiments/428-zero-gene-audit/run_arm_generic.sh}"
T="$(mktemp -d "$HERE/test591.XXXXXX")"; trap 'rm -rf "$T"' EXIT
cd "$REPO_ROOT" || exit 2
ORIG_COPY="$HERE/original_run_arm_generic_c58b9553.sh"   # byte copy of the pre-591 script (sha c58b95533782c597)
[ "$(shasum -a 256 "$ORIG_COPY" | cut -c1-16)" = c58b95533782c597 ] || { echo "FAIL: the saved original copy changed"; exit 1; }
if grep -q "memory branch" "$TARGET"; then   # after --apply: test the live (patched) file against the saved original
  cp "$TARGET" "$T/patched.sh"; echo "(target already patched: testing the live file against the saved original)"
else
  uv run python -c "import sys; sys.path.insert(0,'experiments/591-prep-gate'); import patch_gate591 as P; sys.stdout.write(P.patched(open('$TARGET').read()))" > "$T/patched.sh" || { echo "FAIL: patch did not apply"; exit 1; }
fi
cp "$ORIG_COPY" "$T/original.sh"
cut_block(){ awk '/^log\(\)\{/{p=1} /^restore\(\)\{/{p=0} p' "$1"; }
cut_block "$T/original.sh" > "$T/orig_block.sh"; cut_block "$T/patched.sh" > "$T/pat_block.sh"
[ -s "$T/orig_block.sh" ] && [ -s "$T/pat_block.sh" ] || { echo "FAIL: could not cut the gate block"; exit 1; }

mkdir -p "$T/stubs"
cat > "$T/stubs/memory_pressure" <<'EOS'
#!/bin/bash
echo "The system has 68719476736 (4194304 pages with a page size of 16384)."
echo "System-wide memory free percentage: ${STUB_PCT}%"
EOS
cat > "$T/stubs/vm_stat" <<'EOS'
#!/bin/bash
inact=$(awk -v g="$STUB_FI" 'BEGIN{printf "%d", g*65536}')
echo "Mach Virtual Memory Statistics: (page size of 16384 bytes)"
echo "Pages free:                                       0."
echo "Pages active:                                1000000."
echo "Pages inactive:                         $inact."
EOS
cat > "$T/stubs/sysctl" <<'EOS'
#!/bin/bash
echo "total = 4096.00M  used = 1000.00M  free = ${STUB_SWAP}  (encrypted)"
EOS
cat > "$T/stubs/uv" <<'EOS'
#!/bin/bash
echo "${STUB_OTHERS:-0}"
EOS
cat > "$T/stubs/ps" <<'EOS'
#!/bin/bash
echo "/bin/bash some_chain.sh"
[ "${STUB_PREP:-0}" = 1 ] && echo "/Users/x/.local/bin/vcc prep pred.h5ad -g gene_names.csv -o p.vcc"
exit 0
EOS
cat > "$T/stubs/df" <<'EOS'
#!/bin/bash
echo "Filesystem 1024-blocks Used Available Capacity iused ifree %iused Mounted on"
awk -v g="$STUB_DISK" 'BEGIN{printf "/dev/disk3s5 1 1 %d 99%% 1 1 1%% /System/Volumes/Data\n", g*1048576}'
EOS
chmod +x "$T/stubs/"*

run_gate(){ # $1 = block file, $2 = mode ("" = no argument, as the original calls it); prints "rc|last log line"
  ( export PATH="$T/stubs:$PATH" REPO="$T" LOG="$T/gate.log"; : > "$LOG"
    eval "$(cat "$1")"
    if [ -n "$2" ]; then gate_ok "$2" >/dev/null 2>&1; else gate_ok >/dev/null 2>&1; fi; rc=$?
    echo "$rc|$(tail -1 "$LOG")" )
}

fails=0; n=0
for pct in 39 40 54 55 81; do
 for fi in 22.9 23.0 30.0; do
  for swap in 0.00M 399.00M 400.00M 1.50G; do
   for base in ok others1 prep1 disk19 disk20; do
    case $base in
      ok) o=0; p=0; d=50;; others1) o=1; p=0; d=50;; prep1) o=0; p=1; d=50;; disk19) o=0; p=0; d=19.9;; disk20) o=0; p=0; d=20;;
    esac
    export STUB_PCT=$pct STUB_FI=$fi STUB_SWAP=$swap STUB_OTHERS=$o STUB_PREP=$p STUB_DISK=$d
    swapmib=$(awk -v s="$swap" 'BEGIN{u=substr(s,length(s),1); v=substr(s,1,length(s)-1)+0; if(u=="G") v=v*1024; print v}')
    base_ok=$(awk -v o=$o -v p=$p -v d=$d 'BEGIN{print (o==0 && p==0 && d>=20)?1:0}')
    fi_br=$(awk -v c=$pct -v f=$fi 'BEGIN{print (c>=40 && f>=23)?1:0}')
    pr_br=$(awk -v c=$pct -v s=$swapmib 'BEGIN{print (c>=55 && s>=400)?1:0}')
    orig="$(run_gate "$T/orig_block.sh" "")"; orc=${orig%%|*}
    for mode in build prep; do
      n=$((n+1))
      out="$(run_gate "$T/pat_block.sh" "$mode")"; rc=${out%%|*}; line=${out#*|}
      if [ "$mode" = build ]; then want_branch=$([ $fi_br = 1 ] && echo fi23 || echo none)
      else want_branch=$([ $fi_br = 1 ] && echo fi23 || { [ $pr_br = 1 ] && echo pressure55 || echo none; }); fi
      want_ok=$([ $base_ok = 1 ] && [ $want_branch != none ] && echo 0 || echo 1)
      bad=""
      [ "$rc" = "$want_ok" ] || bad="$bad rc=$rc want=$want_ok"
      [ "$mode" = build ] && [ "$rc" != "$orc" ] && bad="$bad R1(build!=orig $orc)"
      if [ "$mode" = prep ] && [ "$rc" != "$orc" ] && ! { [ $pr_br = 1 ] && [ $base_ok = 1 ] && [ $fi_br = 0 ]; }; then
        bad="$bad R2(prep!=orig $orc outside the pressure branch)"; fi
      case "$line" in *" gate: free ${pct}%"*"mode ${mode}  memory branch ${want_branch}") ;; *) bad="$bad R3(line: $line)";; esac
      if [ -n "$bad" ]; then fails=$((fails+1)); [ $fails -le 10 ] && echo "FAIL pct=$pct fi=$fi swap=$swap base=$base mode=$mode:$bad"; fi
    done
   done
  done
 done
done
echo "R1-R3: $n patched gate calls checked against the original and the spec, $fails failure(s)"

# R4: call sites and syntax
r4=0
/bin/bash -n "$T/patched.sh" || r4=1
[ "$(grep -E '^wait_gate (build|prep)$' "$T/patched.sh" | tr '\n' ' ')" = "wait_gate build wait_gate prep " ] || { echo "R4 call-site order wrong"; r4=1; }
[ "$(grep -cE '^wait_gate$' "$T/patched.sh")" = 0 ] || { echo "R4 a bare wait_gate call remains"; r4=1; }
grep -n "^wait_gate build$" "$T/patched.sh" | head -1 | cut -d: -f1 > "$T/l1"; grep -n 'DRIVER invocation 1' "$T/patched.sh" | head -1 | cut -d: -f1 > "$T/l2"
grep -n "^wait_gate prep$" "$T/patched.sh" | head -1 | cut -d: -f1 > "$T/l3"; grep -n 'INVOCATION 2: vcc prep' "$T/patched.sh" | head -1 | cut -d: -f1 > "$T/l4"
[ "$(( $(cat "$T/l1") + 1 ))" = "$(cat "$T/l2")" ] && [ "$(( $(cat "$T/l3") + 1 ))" = "$(cat "$T/l4")" ] || r4=1
nd=$(diff "$T/original.sh" "$T/patched.sh" | grep -c '^[<>]')
echo "R4: bash -n + call sites $( [ $r4 = 0 ] && echo PASS || echo FAIL ); diff lines original->patched: $nd"
[ $r4 = 0 ] || fails=$((fails+1))

# R5: live readers (informational + numeric parse)
live="$( ( REPO="$REPO_ROOT" LOG="$T/live.log"; : > "$LOG"; eval "$(cat "$T/pat_block.sh")"; gate_ok prep >/dev/null 2>&1; echo "rc=$?" ) )"
echo "R5 live: $(tail -1 "$T/live.log")  $live"
tail -1 "$T/live.log" | grep -qE "swap free [0-9]+ MiB  mode prep  memory branch (fi23|pressure55|none)" || { echo "R5 FAIL"; fails=$((fails+1)); }

[ $fails = 0 ] && echo "ALL PASS" || { echo "FAILED ($fails)"; exit 1; }
