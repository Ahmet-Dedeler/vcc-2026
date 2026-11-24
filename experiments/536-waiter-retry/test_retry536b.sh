#!/bin/bash
# test_retry536b.sh -- 536b on a minimal waiter built from 536's own block, with a vcc shim and a fake vcc state.json
# (HOME override), VCC_RETRY_BASE_S=1 so the backoff runs in seconds. ~2 min. exit 0 = all PASS.
#   R1 resume path: 3 failures then success; sleeps 1, 2, 4 s; every retry resumes the pending entry
#   R2 cap: always failing, cap ~1-2 min ahead -> stops with "no time left before the cap", never sleeps across it
#   R3 an entry existed and no upload is pending -> no retry (double-submit guard), 1 call
#   R4 a failure before any entry -> one fresh retry, then success
#   R5 the patcher refuses a file a live bash process is executing; R6 a second patch is a no-op
set -u
R=~/Code/virtual-cell-challenge-2026; X=$R/experiments/536-waiter-retry
TD=$(mktemp -d -t r536b); PASS=0; FAIL=0
trap 'rm -rf "$TD"' EXIT
ok(){ if [ "$1" = 0 ]; then PASS=$((PASS + 1)); echo "PASS $2"; else FAIL=$((FAIL + 1)); echo "FAIL $2"; fi; }
mkdir -p "$TD/shim" "$TD/home/.config/vcc"
cat > "$TD/shim/vcc" <<'S'
#!/bin/bash
n=$(( $(cat "$SHIMD/calls" 2>/dev/null || echo 0) + 1 )); echo $n > "$SHIMD/calls"; echo "$(date +%s) $*" >> "$SHIMD/argv"
if [ "$n" -gt "${FAILS:-0}" ]; then echo "submitted, entry ABCDEFGHIJKLMNOPQRST"; exit 0; fi
echo "${FAILTEXT:-Upload failed after 6 retries: [Errno 8] nodename nor servname provided}"; exit 1
S
chmod +x "$TD/shim/vcc"
mk_waiter(){ # $1 = path, $2 = cap (YYYYmmddHHMM): a minimal waiter carrying 536's exact block
  uv run --project "$R" python - "$1" "$2" <<'PY'
import sys
sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/experiments/536-waiter-retry")
import importlib.util
spec = importlib.util.spec_from_file_location("p536", "~/Code/virtual-cell-challenge-2026/experiments/536-waiter-retry/patch_retry536.py")
src = open(spec.origin).read()
block = src.split("BLOCK = r'''", 1)[1].split("'''", 1)[0]
import patch_retry536b as P
pend = block.split('OUT="$(vcc submit', 1)[0]
body = ('#!/bin/bash\nset -uo pipefail\nlog(){ printf \'%s %s\\n\' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }\n'
        + pend + P.OLD.replace("__CAP__", sys.argv[2]) + '\nlog "SUBMIT rc=$RC"\n'
        + '[ "$(date -u +%Y%m%d%H%M)" -ge ' + sys.argv[2] + ' ] && stop "clock is past 12:00Z"\n')
open(sys.argv[1], "w").write(body)
PY
}
run(){ # $1 = case dir; env: FAILS, FAILTEXT, PENDING (0/1)
  local d=$1; mkdir -p "$d"; rm -f "$TD/home/.config/vcc/state.json"
  [ "${PENDING:-0}" = 1 ] && printf '{"profiles": {"default": {"pending_uploads": {"PEND123": {"local_path": "%s"}}}}}' "$TD/pkg.vcc" > "$TD/home/.config/vcc/state.json"
  env SHIMD="$d" FAILS="${FAILS:-0}" FAILTEXT="${FAILTEXT:-Upload failed after 6 retries: [Errno 8] nodename nor servname provided}" \
    HOME="$TD/home" PATH="$TD/shim:$PATH" PKG="$TD/pkg.vcc" NAME=joke LOG="$d/log" VCC_RETRY_BASE_S=1 bash "$d/arm9999.sh"; }
touch "$TD/pkg.vcc"
CAP_FAR=$(date -u -v+6H +%Y%m%d%H%M)
# R1
D=$TD/r1; mkdir -p $D; mk_waiter $D/arm9999.sh $CAP_FAR; uv run --project $R python $X/patch_retry536b.py $D/arm9999.sh > /dev/null
FAILS=3 PENDING=1 run $D
[ "$(cat $D/calls)" = 4 ] && grep -q "SUBMIT rc=0" $D/log; ok $? "R1 3 failures then success: 4 calls, rc 0"
grep -q "in 1 s" $D/log && grep -q "in 2 s" $D/log && grep -q "in 4 s" $D/log; ok $? "R1 backoff 1, 2, 4 s (base 1)"
[ "$(grep -c -- "--resume PEND123" $D/argv)" = 3 ]; ok $? "R1 every retry resumes the pending entry"
# R2
D=$TD/r2; mkdir -p $D; CAP_NEAR=$(date -u -v+2M +%Y%m%d%H%M); mk_waiter $D/arm9999.sh $CAP_NEAR
uv run --project $R python $X/patch_retry536b.py $D/arm9999.sh > /dev/null
T0=$(date +%s); FAILS=100000 PENDING=1 run $D; T1=$(date +%s)
CAPS=$(date -j -u -f %Y%m%d%H%M $CAP_NEAR +%s)
grep -q "no time left before the cap: stopping" $D/log && grep -q "SUBMIT rc=1" $D/log; ok $? "R2 stops with 'no time left before the cap'"
last=$(tail -1 $D/argv | cut -d' ' -f1); [ "$last" -lt "$CAPS" ] && [ "$T1" -lt "$CAPS" ]; ok $? "R2 no call and no sleep reaches the cap ($((T1 - T0)) s)"
# R3
D=$TD/r3; mkdir -p $D; mk_waiter $D/arm9999.sh $CAP_FAR; uv run --project $R python $X/patch_retry536b.py $D/arm9999.sh > /dev/null
FAILS=100000 FAILTEXT="error: entry XYZ is scoring" PENDING=0 run $D
[ "$(cat $D/calls)" = 1 ] && grep -q "not retrying (a fresh submit could double-submit)" $D/log; ok $? "R3 entry without a pending upload: no retry"
# R4
D=$TD/r4; mkdir -p $D; mk_waiter $D/arm9999.sh $CAP_FAR; uv run --project $R python $X/patch_retry536b.py $D/arm9999.sh > /dev/null
FAILS=1 FAILTEXT="connection reset by peer" PENDING=0 run $D
[ "$(cat $D/calls)" = 2 ] && grep -q "fresh retry in 1 s" $D/log && grep -q "SUBMIT rc=0" $D/log; ok $? "R4 failure before any entry: one fresh retry, then rc 0"
# R5 / R6
D=$TD/r5; mkdir -p $D; mk_waiter $D/arm9999.sh $CAP_FAR; printf 'sleep 20\nexit 0\n' | cat - $D/arm9999.sh > $D/live.sh; mv $D/live.sh $D/arm9999.sh
( cd $TD && bash r5/arm9999.sh > /dev/null 2>&1 & echo $! > $D/pid ); sleep 1
uv run --project $R python $X/patch_retry536b.py $D/arm9999.sh 2>&1 | grep -q "live bash process is executing it"; ok $? "R5 refuses a file a live bash process is executing"
kill $(cat $D/pid) 2>/dev/null; sleep 1
uv run --project $R python $X/patch_retry536b.py $TD/r1/arm9999.sh | grep -q "already 536b"; ok $? "R6 a second patch is a no-op"
echo; echo "$PASS PASS, $FAIL FAIL"; [ "$FAIL" -eq 0 ]
