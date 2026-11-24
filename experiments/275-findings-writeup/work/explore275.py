"""Explore saved JSON files: print key structure to a given depth, or paths matching regexes.
usage: explore275.py keys <file> [depth]
       explore275.py grep <file> <regex> [regex...]
Read-only. Used by 275 to re-verify numbers quoted in docs/FINDINGS-2026-09.md."""
import json, sys, re

def load(p):
    with open(p) as f:
        return json.load(f)

def keys(o, depth, prefix=""):
    if depth < 0:
        return
    if isinstance(o, dict):
        for k, v in o.items():
            desc = type(v).__name__
            if isinstance(v, dict):
                desc += f"[{len(v)}]"
            elif isinstance(v, list):
                desc += f"[{len(v)}]"
            else:
                desc += f"={str(v)[:40]}"
            print(f"{prefix}{k}: {desc}")
            if isinstance(v, (dict, list)):
                keys(v, depth - 1, prefix + "  ")
    elif isinstance(o, list) and o and isinstance(o[0], dict):
        print(f"{prefix}[0]:")
        keys(o[0], depth - 1, prefix + "  ")

def walk(o, path, out):
    if isinstance(o, dict):
        for k, v in o.items():
            walk(v, path + [str(k)], out)
    elif isinstance(o, list):
        if all(not isinstance(x, (dict, list)) for x in o):
            out.append(("/".join(path), o))
        else:
            for i, v in enumerate(o):
                walk(v, path + [f"[{i}]"], out)
    else:
        out.append(("/".join(path), o))

mode, fn = sys.argv[1], sys.argv[2]
o = load(fn)
if mode == "keys":
    keys(o, int(sys.argv[3]) if len(sys.argv) > 3 else 1)
else:
    pats = [re.compile(p) for p in sys.argv[3:]]
    out = []
    walk(o, [], out)
    for p, v in out:
        if all(r.search(p) for r in pats):
            s = str(v)
            print(f"{p}: {s[:160]}")
