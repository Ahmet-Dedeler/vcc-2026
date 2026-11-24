import json, glob, os
files = sorted(glob.glob("~/.claude/projects/-Users-ahmet-Code-virtual-cell-challenge-2026/*.jsonl"))
for f in files:
    if "76a42756" in f: continue
    for ln, line in enumerate(open(f, errors="replace")):
        try: rec = json.loads(line)
        except Exception: continue
        msg = rec.get("message") or {}
        content = msg.get("content")
        if not isinstance(content, list): continue
        for c in content:
            if not isinstance(c, dict): continue
            if c.get("type") != "tool_use": continue
            s = json.dumps(c.get("input", {}))
            if "whitener_union" in s or "Whitener.fit" in s or "whiten-basis-perts" in s or "fit_union" in s:
                print("="*100)
                print(f"FILE {os.path.basename(f)} line {ln} ts={rec.get('timestamp')} tool={c.get('name')}")
                print(s[:4000])
