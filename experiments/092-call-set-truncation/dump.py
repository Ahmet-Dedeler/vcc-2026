import json, sys
f="~/.claude/projects/-Users-ahmet-Code-virtual-cell-challenge-2026/d207ce74-fe80-4325-90f0-b719244ed1ee.jsonl"
lo,hi=int(sys.argv[1]),int(sys.argv[2])
lines=open(f,errors="replace").readlines()
for ln in range(lo,hi):
    rec=json.loads(lines[ln]); msg=rec.get("message") or {}
    content=msg.get("content")
    if not isinstance(content,list): continue
    for c in content:
        if not isinstance(c,dict): continue
        t=c.get("type")
        if t=="tool_use":
            print(f"\n########## line {ln} {rec.get('timestamp')} TOOL_USE {c.get('name')}")
            inp=c.get("input",{})
            print(inp.get("command") or inp.get("content") or json.dumps(inp)[:3000])
        elif t=="tool_result":
            out=c.get("content")
            if isinstance(out,list): out=" ".join(x.get("text","") for x in out if isinstance(x,dict))
            print(f"\n---------- line {ln} RESULT:\n{str(out)[:2500]}")
        elif t=="text" and c.get("text","").strip():
            print(f"\n>>>>> line {ln} ASSISTANT: {c['text'][:1200]}")
