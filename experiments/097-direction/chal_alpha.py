"""Does the per-GENE key have the same shape on the real 300-target panel as on the harnesses?

The harness model carries TWO donors; production carries THREE, and
`MultiDonorModel.donor_agreement` requires sign unanimity across every carrying donor. With three
donors the unanimity bar is stricter AND the gene axis is the intersection of their coverages, so
the fraction of genes with a nonzero agreement value -- which is what the midrank alpha grades --
should be smaller on the challenge than on the harnesses. If it is much smaller, the alpha
distribution differs and the harness-measured multiplier is being read off a different estimator.
"""
from __future__ import annotations
import sys
import numpy as np, pandas as pd
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
from vcc2026.donor_transfer import Donor
PROD = ("k562_gwps", "orion_hct116", "orion_hek293t")
targets = sorted(set(pd.read_csv(f"{ROOT}/data/challenge/pert_counts.csv")["target_gene"].astype(str)))
genes = np.asarray(pd.read_csv(f"{ROOT}/data/challenge/gene_names.csv").iloc[:, 0].astype(str))
Ms, hv = {}, {}
for nm in PROD:
    dn = Donor.load(f"{ROOT}/data/donors/{nm}"); m, s = dn.decompose()
    gi = {g: i for i, g in enumerate(dn.genes)}
    cols = np.array([gi.get(g, -1) for g in genes]); ok = cols >= 0
    pi = {p: i for i, p in enumerate(dn.perts)}
    M = np.zeros((len(targets), genes.size), np.float32); h = np.zeros(len(targets), bool)
    for r, t in enumerate(targets):
        if t in pi: M[r, ok] = s[pi[t]][cols[ok]]; h[r] = True
    Ms[nm] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
    hv[nm] = h; del dn, s
frac, nd = [], []
for r in range(len(targets)):
    have = [n for n in PROD if hv[n][r]]
    nd.append(len(have))
    if len(have) < 2:
        frac.append(np.nan); continue
    st = np.vstack([Ms[n][r] for n in have])
    sg = np.sign(st)
    unan = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
    frac.append(float(unan.mean()))
frac = np.array(frac); nd = np.array(nd)
print(f"# challenge panel, fraction of the {genes.size} genes with a NONZERO donor-agreement value")
for k in (2, 3):
    m = nd == k
    if m.sum():
        print(f"  targets with {k} donors ({int(m.sum())}):  median {100*np.nanmedian(frac[m]):.1f} %"
              f"   p5 {100*np.nanpercentile(frac[m],5):.1f} %  p95 {100*np.nanpercentile(frac[m],95):.1f} %")
print(f"  ALL scoreable targets ({int(np.isfinite(frac).sum())}): median "
      f"{100*np.nanmedian(frac):.1f} %")
print(f"  harness comparison (2 donors, measured in tiefix.py): 20.8 % (h1-test), 20.6 % (orion)")
print(f"  -> the midrank alpha therefore grades a {100*np.nanmedian(frac):.0f} % minority of genes "
      f"on the challenge and damps the rest to one shared value; on the harnesses that minority "
      f"is 21 %.")
