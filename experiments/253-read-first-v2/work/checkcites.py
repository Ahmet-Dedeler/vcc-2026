"""Every (KEY:n[-m][,n[-m]]) citation in READ-FIRST.md must resolve to existing lines."""
import re, glob, sys
p = "experiments/253-read-first-v2/READ-FIRST.md"
s = open(p, encoding="utf-8").read()
files = {"A": "AGENTS.md", "BL": "docs/BOARD-LEDGER.md", "RB": "docs/DEF-RUNBOOK.md", "DP": "scripts/def_pipeline.py",
         "MB": "scripts/make_submission_pm_cd4.py", "SL": "submissions/SUBMITTED_LEDGER.json",
         "LB": "experiments/253-read-first-v2/work/ranks.json", "VS": "experiments/253-read-first-v2/work/vcc-status.log",
         "P237": "experiments/237-gate-validation/PREREGISTRATION.md", "AD": "experiments/237-gate-validation/ADDENDUM-243.md",
         "P249": "experiments/249-cd4-deflate/PREREGISTRATION.md"}
def resolve(key):
    if key in files: return files[key]
    m = re.match(r"^(BR|R)(\d+)$", key)
    if m:
        pat = f"experiments/{m.group(2)}-*/" + ("BOARD-RESULT*.md" if m.group(1) == "BR" else "RESULTS.md")
        g = sorted(glob.glob(pat)); return g[0] if g else None
    m = re.match(r"^E(\d+)/(.+)$", key)
    if m:
        g = sorted(glob.glob(f"experiments/{m.group(1)}-*/{m.group(2)}")); return g[0] if g else None
    return None
cites = re.findall(r"\(((?:[A-Z]+\d*(?:/[\w./-]+)?):\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*(?:; ?[A-Z]+\d*(?:/[\w./-]+)?:\d+(?:-\d+)?(?:,\d+(?:-\d+)?)*)*)\)", s)
bad, n = [], 0
for c in cites:
    for part in re.split(r"; ?", c):
        key, rngs = part.split(":", 1)
        f = resolve(key)
        if f is None: bad.append((part, "NO FILE")); continue
        nlines = sum(1 for _ in open(f, encoding="utf-8", errors="replace"))
        for r in rngs.split(","):
            a, b = (r.split("-") + [None])[:2]; a = int(a); b = int(b) if b else a
            n += 1
            if not (1 <= a <= b <= nlines): bad.append((part, f"{f} has {nlines} lines"))
print("citations checked:", n, "bad:", len(bad))
for b in bad: print("  BAD", b)
# bare cites like SL:31 inside text without parens
extra = re.findall(r"\b(?:SL|LB|VS|RB|BL|A|DP|MB|R\d+|BR\d+|P237|P249|AD):\d+", s)
print("total cite tokens (incl. inline):", len(extra))
