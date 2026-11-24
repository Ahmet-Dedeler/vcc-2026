#!/usr/bin/env python3
"""198 — high-precision flags. Only literal constructions that the invalidating
findings condemn, not keyword co-occurrence."""
from __future__ import annotations
import json, re
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
rows = json.loads((OUT / "enumeration.json").read_text())["rows"]

HARD = {
    # literal use of the instrument audit 181 condemned
    "H_n_real_target": re.compile(r"n_real_target"),
    # the refuted 2,670 inference
    "H_2670": re.compile(r"2,670|n_real\s*[≈~=]\s*2670|\b2670\b"),
    # debench Surface used to size a pds arm  (188 Phase 3)
    "H_surface_pds": re.compile(r"Surface\s*\([^)]*\)[^\n]{0,200}pds|\.pds\(\)|pds\(\)\s"),
    # a DE pool built by top-N magnitude  (188 Phase 1)
    "H_topN_magnitude": re.compile(
        r"top[\s\-]?\d{2,5}\s+(?:genes\s+)?by\s+\|?(?:lfc|pred|effect|delta|δ)"
        r"|top[\s\-]?\d{2,5}\s+\|lfc\||magnitude[\s\-]selected|by \|pred_lfc\|", re.I),
    # rule 13's emission toll applied
    "H_toll_0146": re.compile(r"0\.0146"),
    # explicit absolute pds/reach/fid/jac/n_real quoted ON an orion surface line
    "H_orion_absolute": None,   # computed below, line-level
    # h1-test reach on an emission-law change
    "H_emission_reach": None,   # line-level
}

EMISSION_LAW = re.compile(
    r"resample|emit_counts|emission law|emitter|dispersion|overdispersion|"
    r"library size|per-cell depth|NB emis|quasi|spread", re.I)
ORION = re.compile(r"orion|hct116|hek293t|HCT-116|HEK-?293", re.I)
H1T = re.compile(r"h1[\s\-_]?test", re.I)
ABS_MEMBER = re.compile(r"\b(pds|reach|fid|jac|n_real|cosine)\b", re.I)
NUM = re.compile(r"(?<![\w.])\d\.\d{3,}|(?<![\w.])0\.\d{2,}")

out = []
for r in rows:
    p = ROOT / r["path"]
    text = p.read_text(errors="replace")
    lines = text.split("\n")
    hits = {}
    for k, pat in HARD.items():
        if pat is None:
            continue
        m = [(i + 1, lines[i].strip()[:220]) for i in range(len(lines)) if pat.search(lines[i])]
        if m:
            hits[k] = m[:6]

    # line-level: an absolute member number quoted on the same line as an Orion surface
    orion_abs = []
    for i, ln in enumerate(lines):
        if ORION.search(ln) and ABS_MEMBER.search(ln) and NUM.search(ln) and not H1T.search(ln):
            # exclude pure delta language
            if re.search(r"\b(delta|Δ|d ?pds|change|vs\.?|difference|improve)", ln, re.I):
                continue
            orion_abs.append((i + 1, ln.strip()[:220]))
    if orion_abs:
        hits["H_orion_absolute"] = orion_abs[:6]

    # h1-test reach quoted on an emission-law arm
    em_reach = []
    for i, ln in enumerate(lines):
        if re.search(r"\breach\b", ln, re.I) and EMISSION_LAW.search(ln) and NUM.search(ln):
            em_reach.append((i + 1, ln.strip()[:220]))
    if em_reach and ("h1-test" in r["surfaces"] or "reach" in r["endpoints"]):
        hits["H_emission_reach"] = em_reach[:6]

    if hits:
        out.append(dict(exp=r["exp"], title=r["title"], verdict=r["verdict"],
                        primary_endpoint=r["primary_endpoint"],
                        primary_surface=r["primary_surface"],
                        hard=hits, path=r["path"]))

(OUT / "hardflags.json").write_text(json.dumps(out, indent=1))

c = Counter()
byclass = defaultdict(list)
for o in out:
    for k in o["hard"]:
        c[k] += 1
        byclass[k].append(o["exp"])
print(f"{len(out)} / {len(rows)} experiments carry >=1 hard flag\n")
for k, v in c.most_common():
    print(f"{k:22s} {v:3d}  {sorted(byclass[k])}")
    print()
