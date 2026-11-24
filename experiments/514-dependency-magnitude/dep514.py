"""514: does line-specific dependency (DepMap Chronos, 24Q4 subset from 360) predict how much stronger a target's
transcriptional response is in one line than another? Public data only; no challenge truth.

For targets screened in two donor lines X, Y (same-target rows in data/donors/<bank>/delta.npy, natural-log pseudobulk
effects at TS 5e4, the scorer's units), energy E = sum over non-target genes of delta^2 (shared gene axis). We regress
log(E_X / E_Y) on (chronos_X - chronos_Y). A strong negative slope means: more essential in X -> bigger response in X.
Pairs with DepMap models: K562 (ACH-000551), HCT116 (ACH-000971); HEK293T is not a DepMap model.
Also reported: within one line, log E vs chronos (does dependency predict response size at all, on non-essential targets?).
    uv run python experiments/514-dependency-magnitude/dep514.py
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEP = pd.read_csv(ROOT / "experiments/360-depmap-three-lines/work/CRISPRGeneEffect_subset.csv", index_col=0)
DEP.columns = [c.split(" (")[0] for c in DEP.columns]
IDS = {"k562_gwps": "ACH-000551", "orion_hct116": "ACH-000971", "nadig_hepg2": "ACH-000739", "nadig_jurkat": "ACH-000995"}


def bank(name):
    b = ROOT / "data/donors" / name
    perts = [l.strip() for l in open(b / "perts.txt")]
    genes = [l.strip() for l in open(b / "genes.txt")]
    return perts, genes, np.load(b / "delta.npy", mmap_mode="r")


def energies(name, targets=None):
    perts, genes, M = bank(name)
    gi = {g: i for i, g in enumerate(genes)}
    out = {}
    for i, p in enumerate(perts):
        if targets is not None and p not in targets:
            continue
        row = np.asarray(M[i], float)
        if p in gi:
            row = row.copy(); row[gi[p]] = 0.0
        out[p] = float(np.nansum(row ** 2))
    return out, genes


def chron(line, genes):
    row = DEP.loc[IDS[line]] if IDS[line] in DEP.index else None
    return {g: float(row[g]) for g in genes if row is not None and g in row.index and np.isfinite(row[g])}


print("DepMap rows present:", {k: v in DEP.index for k, v in IDS.items()})
res = {}
for name in ("k562_gwps", "orion_hct116"):
    E, _ = energies(name)
    C = chron(name, list(E))
    common = [t for t in E if t in C and E[t] > 0]
    le = np.log([E[t] for t in common]); ch = np.array([C[t] for t in common])
    for lab, m in [("all", np.ones(len(ch), bool)), ("non-essential (chronos > -0.5)", ch > -0.5),
                   ("chronos in [-0.5, 0]", (ch > -0.5) & (ch < 0))]:
        r = np.corrcoef(ch[m], le[m])[0, 1]
        print(f"{name:13s} within-line: log E vs chronos, {lab:30s} n {m.sum():5d}  r {r:+.3f}")
    res[name] = (E, C)
# cross-line pair
Ek, Ck = res["k562_gwps"]; Eh, Ch = res["orion_hct116"]
common = [t for t in Ek if t in Eh and t in Ck and t in Ch and Ek[t] > 0 and Eh[t] > 0]
lr = np.log(np.array([Eh[t] / Ek[t] for t in common]))
dc = np.array([Ch[t] - Ck[t] for t in common])
mean_c = np.array([(Ch[t] + Ck[t]) / 2 for t in common])
print(f"\nHCT116 vs K562: {len(common)} shared targets with dependency in both")
for lab, m in [("all", np.ones(len(dc), bool)), ("both non-essential (> -0.5)", (np.array([Ch[t] for t in common]) > -0.5) & (np.array([Ck[t] for t in common]) > -0.5)),
               ("|dchronos| > 0.3", np.abs(dc) > 0.3)]:
    r = np.corrcoef(dc[m], lr[m])[0, 1]
    b = np.polyfit(dc[m], lr[m], 1)[0]
    print(f"  log(E_hct/E_k562) vs (chr_hct - chr_k562), {lab:28s} n {m.sum():5d}  r {r:+.3f}  slope {b:+.2f} per unit chronos")
# panel view: A/C targets with divergent dependency
pdp = pd.read_csv(ROOT / "experiments/360-depmap-three-lines/work/panel_dependency.csv")
for c in ("A", "C"):
    dd = pdp[f"chronos_{c}"] - pdp["donor_mean"]
    print(f"panel, context {c}: targets with chronos_{c} - donor_mean < -0.3: {(dd < -0.3).sum()}; > +0.3: {(dd > 0.3).sum()}")
