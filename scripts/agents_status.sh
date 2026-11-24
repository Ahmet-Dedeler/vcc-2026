#!/usr/bin/env bash
# One line per subagent that was ever launched in the last two days: alive?, elapsed, model, last words.
cd "$(dirname "${BASH_SOURCE[0]}")/.."
for d in experiments/2[3-9][0-9]-*/; do
  n=$(basename "$d"); [ -f "$d/claude-run.pid" ] || continue
  pid=$(cat "$d/claude-run.pid")
  if kill -0 "$pid" 2>/dev/null; then st="RUNNING $(ps -p "$pid" -o etime= | tr -d ' ')"; else st="done"; fi
  [ -f "$d/RESULTS.md" ] && st="$st +RESULTS"
  m=$(grep -o '"model":"[^"]*"' "$d/claude-run.jsonl" 2>/dev/null | head -1 | cut -d'"' -f4)
  last=$(python3 - "$d/claude-run.jsonl" <<'PY'
import json,sys
t=""
for l in open(sys.argv[1]):
    if l.startswith('{"type":"assistant"'):
        try:
            for c in json.loads(l)["message"]["content"]:
                if c.get("type")=="text" and c["text"].strip(): t=c["text"].strip().replace("\n"," ")
        except Exception: pass
print(t[:110])
PY
)
  printf "%-28s %-22s %-18s %s\n" "$n" "$st" "${m:-?}" "$last"
done
