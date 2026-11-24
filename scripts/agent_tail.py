#!/usr/bin/env python3
"""Follow a headless subagent's transcript as plain text.

    python3 scripts/agent_tail.py 246-direction-hypotheses          # history so far, then follow live
    python3 scripts/agent_tail.py 246-direction-hypotheses --all    # also show tool results (truncated)
    python3 scripts/agent_tail.py 246-direction-hypotheses --no-follow

Prints the agent's own words, every tool call it makes, and its final result. Ctrl-C to stop.
"""
import json, os, sys, time

name = next((a for a in sys.argv[1:] if not a.startswith("--")), None)
if not name:
    print(__doc__); sys.exit(2)
show_all = "--all" in sys.argv
follow = "--no-follow" not in sys.argv
root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
path = os.path.join(root, "experiments", name, "claude-run.jsonl")

def fmt_tool(c):
    n, inp = c.get("name", "?"), c.get("input", {}) or {}
    if n == "Bash":
        return f"🔧 Bash: {inp.get('description') or (inp.get('command') or '')[:140]}"
    if n in ("Read", "Write", "Edit", "MultiEdit", "NotebookEdit"):
        return f"🔧 {n}: {os.path.relpath(inp.get('file_path', '?'), root) if inp.get('file_path') else '?'}"
    if n in ("Grep", "Glob"):
        return f"🔧 {n}: {inp.get('pattern', '')}  {inp.get('path', '')}"
    return f"🔧 {n}: {json.dumps(inp)[:140]}"

def emit(line):
    try:
        d = json.loads(line)
    except Exception:
        return
    t = d.get("type")
    if t == "assistant":
        for c in d.get("message", {}).get("content", []) or []:
            if c.get("type") == "text" and c["text"].strip():
                print("\n🧠 " + c["text"].strip() + "\n")
            elif c.get("type") == "tool_use":
                print(fmt_tool(c))
    elif t == "user" and show_all:
        for c in d.get("message", {}).get("content", []) or []:
            if isinstance(c, dict) and c.get("type") == "tool_result":
                s = c.get("content")
                if isinstance(s, list):
                    s = " ".join(x.get("text", "") for x in s if isinstance(x, dict))
                s = str(s or "").strip().replace("\n", " ⏎ ")
                print("   ↳ " + s[:400])
    elif t == "result":
        print("\n✅ FINISHED (" + str(d.get("subtype")) + "): " + str(d.get("result", ""))[:3000] + "\n")
    sys.stdout.flush()

pos = 0
buf = ""
while True:
    try:
        with open(path, "r") as f:
            f.seek(pos)
            chunk = f.read()
            pos = f.tell()
    except FileNotFoundError:
        print(f"(no transcript yet at {path})"); time.sleep(2); continue
    buf += chunk
    while "\n" in buf:
        line, buf = buf.split("\n", 1)
        emit(line)
    if not follow:
        break
    time.sleep(1)
