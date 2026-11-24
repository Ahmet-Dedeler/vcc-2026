"""804 — the known answers' board side: each pair's six SCALED members from docs/BOARD-LEDGER.md rows (read, never written),
their deltas and the overall delta. Light. Output: work/board_pairs.json.
    uv run python experiments/804-known-answers/board_pairs804.py
"""
import json, re
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"; OUT.mkdir(exist_ok=True)
M = ("pds", "mse", "nmae", "fid", "reach", "jac")
rows = {}
for ln in (ROOT / "docs/BOARD-LEDGER.md").read_text().splitlines():
    c = [x.strip() for x in ln.strip().strip("|").split("|")]
    if len(c) != 10 or not re.fullmatch(r"[A-Za-z0-9]{20}", c[9]):
        continue
    name = re.sub(r"\*", "", c[0]).split(" = ")[0].split(" ")[0]
    try:
        vals = [float("nan") if v.strip("* ") in ("—", "-", "") else float(v.replace("−", "-").replace("*", "")) for v in c[2:9]]
    except ValueError:
        continue
    rows.setdefault(name, dict(entry=c[9], overall=vals[0], **dict(zip(M, vals[1:]))))
# (cand, base, channel/class, what changed)  -- the DESIGN 9 table plus extras
PAIRS = [
    ("495", "446bj2", "bulk-key", "LODO key cap 3 -> cap 6 (== pct on A/B/C)"),
    ("498p2g", "495", "bulk-shape", "tercile sharpening 2.5/2.0/1.75 on cap 6"),
    ("479p2g", "446bj2", "bulk-shape", "tercile sharpening on cap 3"),
    ("479sp2", "446bj2", "bulk-shape", "strength-tercile sharpening on cap 3"),
    ("702", "ALLSW", "bulk-external", "PIE blend w 0.25 on seen targets"),
    ("446bj2z", "446bj2", "de-mag", "zs (dn_k250_zs)"),
    ("496a", "446bj2", "de-sign", "CLEAN4 uncalled-gene hook 0.05"),
    ("DN", "630arm", "de-mag", "DOWN magnitude pass d_m200"),
    ("630arm", "446bj2z", "de-mag", "UPRN rank-neutral UP magnitude"),
    ("idall3", "498p2g", "de-sign", "ungated identity hook (vs no hook; vs plain hook = idall3 - 496a's delta)"),
    ("446bj2zk20", "446bj2z", "de-callset", "K 2,000"),
    ("446bj2zk16", "446bj2z", "de-callset", "K 1,600"),
    ("458d", "446b", "de-prior", "Arc H1 programme as block prior"),
    ("471c", "446b", "bulk-union", "h1_2025_full union swap"),
    ("460dr", "453dr", "bulk-key", "day-of LODO key on the October path"),
    ("ALLSW", "498p2g", "de-mixed", "zs + CLEAN4 hook on cap6+P2g"),
]
out = {}
for c, b, cls, what in PAIRS:
    if c not in rows or b not in rows:
        out[f"{c}-vs-{b}"] = dict(missing=[x for x in (c, b) if x not in rows]); continue
    d = {m: rows[c][m] - rows[b][m] for m in M}
    d = {m: (v if v == v else 0.0) for m, v in d.items()}      # a member the ledger left blank ("—": DE flat) counts 0
    out[f"{c}-vs-{b}"] = dict(cand=c, base=b, cls=cls, what=what, d_scaled=d, d_overall=rows[c]["overall"] - rows[b]["overall"],
                              d_overall_from_members=sum(d.values()) / 6, entries=[rows[c]["entry"], rows[b]["entry"]])
(OUT / "board_pairs.json").write_text(json.dumps(dict(rows=rows, pairs=out), indent=1))
for k, v in out.items():
    if "missing" in v: print(f"{k:28s} MISSING {v['missing']}"); continue
    print(f"{k:28s} {v['cls']:13s} overall {v['d_overall']:+.7f} (members {v['d_overall_from_members']:+.7f}) " +
          " ".join(f"{m} {v['d_scaled'][m]:+.5f}" for m in M))
