import json, glob, os
files = sorted(glob.glob("~/.claude/projects/-Users-ahmet-Code-virtual-cell-challenge-2026/*.jsonl"))
hits=[]
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
            if "whitener_union" in s:
                hits.append((rec.get("timestamp"), os.path.basename(f), ln, c.get("name"), len(s), s))
hits.sort(key=lambda x: (x[0] or ""))
for ts,f,ln,name,L,s in hits:
    print(f"{ts}  {f[:8]}  line={ln:6d} tool={name:12s} len={L}")
