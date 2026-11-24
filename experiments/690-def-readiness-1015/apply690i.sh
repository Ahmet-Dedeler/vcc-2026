#!/usr/bin/env bash
# 690i applier (10-15 22:10Z, readiness [4592ee]). Applies patch690i (DEF_DAY_H1_STACK + the CD4 grid; inert unset) in the
# window when no def_day run can be hurt: after R-ABC (rabc694) has EXITED and while no def_day.sh / run690 is alive.
# Editing def_day.sh under a live bash run corrupts that run. 708 needs the knob only for a post-R0 final build.
#   1. wait: rabc694 not alive, no def_day.sh / run690.sh under bash, no NOGO file ($W/NOGO-690i), the four target files
#      clean in git (so the pathspec commit carries only this patch); cap APPLY690I_CAP (default 10-16 18:00Z)
#   2. record the live recipe's env (plain `command` and --allow-pending), byte-copy the four files
#   3. patch690i --root <repo>; then on the LIVE tree: bash -n, py_compile, the recipe's env unchanged (both forms),
#      scripts/test_october_recipe.py passes, the helpers' old dump names unchanged
#   4. all pass -> commit the four paths (pathspec); any failure -> restore the byte copies, verify the shas, no commit
# The byte-level known answer comes after: run690's R-KA rebuilds 524dr2's board-scored prediction through the patched
# def_day with every new knob unset. Log: work/apply690i.log; verdict "APPLY690I DONE <commit>" | "APPLY690I FAIL ..."; then END.
# usage: nohup bash scripts/detach.sh bash experiments/690-def-readiness-1015/apply690i.sh >> experiments/690-def-readiness-1015/work/apply690i.nohup 2>&1 &
set -uo pipefail
REPO=~/Code/virtual-cell-challenge-2026; cd "$REPO" || exit 1
E=experiments/690-def-readiness-1015; W=$E/work
if [ -z "${APPLY690I_FROZEN:-}" ]; then
  FZ="$W/.frozen-$$-apply690i.sh"; cp "$REPO/$E/apply690i.sh" "$FZ" || exit 1
  APPLY690I_FROZEN=1 exec bash "$FZ" "$@"
fi
LOG=$W/apply690i.log; CAP="${APPLY690I_CAP:-202610161800}"; BK=$W/690i-backup
FILES="scripts/def_day.sh scripts/def_day_helpers.py scripts/october_recipe.py docs/october-decisions.json"
log(){ printf '%s %s\n' "$(date -u +%FT%TZ)" "$*" >> "$LOG"; }
trap 'rm -f "$W/.frozen-$$-apply690i.sh"; log "END"' EXIT
alive(){ ps -axo comm=,args= | awk -v p="$1" '$1 ~ /bash$/ && index($0, p) {n++} END {exit n ? 0 : 1}'; }
sha(){ shasum -a 256 "$1" | cut -d' ' -f1; }
envs(){ uv run python scripts/october_recipe.py command --bundle B --work W --name N --json "$@" 2>/dev/null \
          | python3 -c 'import json,sys; print(json.dumps(json.load(sys.stdin)["env"], sort_keys=True))' 2>/dev/null || echo "rc!=0"; }
echo $$ > "$W/apply690i.pid"
log "apply690i armed (pid $$): waits for rabc694 to exit, no def_day.sh / run690.sh alive, no $W/NOGO-690i; cap $CAP"
last=""
while :; do
  why=""
  alive rabc694.sh && why="$why rabc694"
  alive scripts/def_day.sh && why="$why def_day"
  alive run690.sh && why="$why run690"
  [ -e "$W/NOGO-690i" ] && why="$why NOGO"
  # shellcheck disable=SC2086
  [ -n "$(git status --porcelain -- $FILES)" ] && why="$why dirty:$(git status --porcelain -- $FILES | awk '{print $2}' | tr '\n' ',')"
  [ -z "$why" ] && break
  [ "$why" != "$last" ] && log "waiting for:$why"; last="$why"
  [ "$(date -u +%Y%m%d%H%M)" -ge "$CAP" ] && { log "APPLY690I FAIL: not applied by cap $CAP (for:$why); apply by hand in a safe window"; exit 3; }
  sleep 60
done
grep -q "690i: the doses come from 438's grid" scripts/def_day.sh && { log "APPLY690I FAIL: already applied"; exit 1; }
PRE1="$(envs)"; PRE2="$(envs --allow-pending)"
log "pre-apply recipe env: $PRE1"
rm -rf "$BK"; mkdir -p "$BK"
for f in $FILES; do cp -p "$f" "$BK/$(echo "$f" | tr '/' '_')"; done
PRESHA="$(for f in $FILES; do sha "$f"; done | tr '\n' ' ')"
restore(){
  for f in $FILES; do cp -p "$BK/$(echo "$f" | tr '/' '_')" "$f"; done
  local now; now="$(for f in $FILES; do sha "$f"; done | tr '\n' ' ')"
  [ "$now" = "$PRESHA" ] && log "restored the four files byte for byte" || log "RESTORE MISMATCH: diagnose by hand ($now vs $PRESHA)"
}
if ! uv run python $E/patch690i.py --root "$REPO" >> "$LOG" 2>&1; then restore; log "APPLY690I FAIL: patch690i refused"; exit 1; fi
bad=""
bash -n scripts/def_day.sh || bad="$bad bash-n"
uv run python -m py_compile scripts/def_day_helpers.py scripts/october_recipe.py || bad="$bad py_compile"
[ "$(envs)" = "$PRE1" ] || bad="$bad env-plain($(envs))"
[ "$(envs --allow-pending)" = "$PRE2" ] || bad="$bad env-allow-pending"
uv run python scripts/test_october_recipe.py >> "$LOG" 2>&1 || bad="$bad test_october_recipe"
uv run python - <<'PY' >> "$LOG" 2>&1 || bad="$bad helpers-names"
import importlib.util
s = importlib.util.spec_from_file_location("h", "scripts/def_day_helpers.py"); h = importlib.util.module_from_spec(s); s.loader.exec_module(h)
assert h._dump_name530("x1", "-") == "dump_h1x2_cd4x1.npz" and h._dump_name530("x1.5", "-") == "dump_h1x2_cd4x15.npz"
assert h._dump_name530("x1.5", "0.25") == "dump_h1x2_cd4x15_s0.25.npz" and h._dump_name530("x1", "-0.25") == "dump_h1x2_cd4x1_sm0.25.npz"
print("helpers: old dump names unchanged")
PY
if [ -n "$bad" ]; then restore; log "APPLY690I FAIL: post-apply checks:$bad (restored; nothing committed)"; exit 1; fi
POST="$(for f in $FILES; do echo "$f $(sha "$f" | cut -c1-16)"; done | tr '\n' ';')"
# shellcheck disable=SC2086
git commit -q -m "690i applied by apply690i.sh: DEF_DAY_H1_STACK (x1|x1.5|x2|x2.5|x3, default x2) + DEF_DAY_CD4_STACK widened to 438's grid (x0.75..x2), dump subshell only, for 708's post-R0 re-pricing. Inert unset: test690i 13/13 (dump-subshell env identical to the unpatched script on all 12 shared combinations), live recipe env unchanged (plain + --allow-pending), test_october_recipe PASS. Byte-level KA to come: run690 R-KA (524dr2's board-scored prediction through the patched def_day). shas: $POST

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>" -- $FILES >> "$LOG" 2>&1 || { restore; log "APPLY690I FAIL: git commit failed (restored)"; exit 1; }
log "APPLY690I DONE $(git log -1 --format=%h): $POST"
