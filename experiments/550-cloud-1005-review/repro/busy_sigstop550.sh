#!/bin/bash
# 550: build549.sh busy() false positive, and SIGSTOP left behind when the pauser is SIGKILL'd.
# Cleans up every process it starts. Does not touch the live tree.
set -uo pipefail
ROOT="$(CDPATH= cd "$(dirname "$0")/../../.." && pwd)"
cd "$ROOT"
fail=0

echo "=== busy() pattern from build549.sh:25 ==="
# The pattern matches any bash whose argv CONTAINS the name, not an argv word that ends in it.
PAT='$1 ~ /^bash$/ && /run_arm_generic\.sh|def_day\.sh|build527\.sh/'
self_line="bash ${ROOT}/experiments/549-cph-arm/build549.sh"
worker_line="bash -c echo run_arm_generic.sh is in the mandate"
real_line="bash ${ROOT}/experiments/428-zero-gene-audit/run_arm_generic.sh TAG"
python_line="python scripts/make_submission_pm_cd4.py --bundle x"
echo "$self_line" | awk "$PAT" | grep -q . && echo "SELF MATCH (bad)" || echo "self: no match"
echo "$worker_line" | awk "$PAT" | grep -q . && echo "worker argv mention: MATCH" || { echo "worker: no match (unexpected)"; fail=1; }
echo "$real_line" | awk "$PAT" | grep -q . && echo "real run_arm: MATCH" || { echo "real: no match"; fail=1; }
echo "$python_line" | awk "$PAT" | grep -q . && echo "python builder matched by the bash awk (bad)" || echo "python builder: awk misses it (other_builders must catch it)"

# live: a bash whose args mention the name, and a python whose args mention the builder
bash -c 'sleep 30; true' run_arm_generic.sh &
BPID=$!
# $0 is run_arm_generic.sh, so the bash argv contains that word. The trailing `true` keeps bash from exec'ing sleep.
sleep 0.3
echo "live bash pid $BPID args: $(ps -p $BPID -o args=)"
ps -axo comm,args | awk '$1 ~ /^bash$/ && /run_arm_generic\.sh|def_day\.sh/' | grep -q . && echo "live busy() would return busy" || { echo "live busy() missed the sleeper"; fail=1; }

uv run python -c 'import time; time.sleep(30)' scripts/make_submission_pm_cd4.py >/tmp/550py.out 2>&1 &
PPID_PY=$!
sleep 0.6
HITS="$(uv run python -c 'import sys; sys.path.insert(0,"scripts"); import def_pipeline as P; print("\n".join(P.other_builders()))')"
echo "other_builders:"
echo "$HITS"
echo "$HITS" | grep -q make_submission_pm_cd4 || { echo "other_builders missed the python sleeper"; fail=1; }
kill "$BPID" "$PPID_PY" 2>/dev/null || true
wait "$BPID" "$PPID_PY" 2>/dev/null || true

echo "=== SIGSTOP: EXIT on SIGTERM resumes; SIGKILL leaves T ==="
WORKDIR="$(mktemp -d)"
cat > "$WORKDIR/child.sh" <<'EOF'
#!/bin/bash
sleep 60
EOF
chmod +x "$WORKDIR/child.sh"
cat > "$WORKDIR/pauser.sh" <<'EOF'
#!/bin/bash
DP="$1"
trap 'kill -CONT "$DP" 2>/dev/null; echo RESUMED >> "$2"' EXIT
kill -STOP "$DP"
echo STOPPED >> "$2"
# no cap, same shape as build549.sh:43
while kill -0 "$DP" 2>/dev/null; do sleep 0.2; done
EOF
chmod +x "$WORKDIR/pauser.sh"
"$WORKDIR/child.sh" &
CPID=$!
sleep 0.2
"$WORKDIR/pauser.sh" "$CPID" "$WORKDIR/log" &
PAUSER=$!
sleep 0.5
STATE="$(ps -p "$CPID" -o stat= | tr -d ' ')"
echo "after STOP, child stat=$STATE"
# SIGTERM runs EXIT
kill -TERM "$PAUSER"
sleep 0.4
STATE2="$(ps -p "$CPID" -o stat= | tr -d ' ')"
echo "after SIGTERM pauser, child stat=$STATE2 log=$(tr '\n' ' ' < "$WORKDIR/log")"
echo "$STATE2" | grep -q T && { echo "SIGTERM left the child stopped"; fail=1; }

# stop again, then SIGKILL the pauser
kill -STOP "$CPID"
: > "$WORKDIR/log"
"$WORKDIR/pauser.sh" "$CPID" "$WORKDIR/log" &
PAUSER=$!
sleep 0.5
kill -KILL "$PAUSER" 2>/dev/null || true
sleep 0.3
STATE3="$(ps -p "$CPID" -o stat= | tr -d ' ')"
echo "after SIGKILL pauser, child stat=$STATE3 (T means the driver is left stopped)"
echo "$STATE3" | grep -q T || { echo "expected T after SIGKILL"; fail=1; }
# cleanup: the review must not leave a stopped process
kill -CONT "$CPID" 2>/dev/null || true
kill "$CPID" 2>/dev/null || true
wait "$CPID" 2>/dev/null || true
rm -rf "$WORKDIR"
if [ "$fail" -ne 0 ]; then
  echo "BUSY/SIGSTOP FAIL"
  exit 1
fi
echo "BUSY/SIGSTOP PASS"
