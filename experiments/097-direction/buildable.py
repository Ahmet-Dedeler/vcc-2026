"""Two things the recommendation has to survive: can it be COMPUTED on the challenge panel, and
what does it do to the DE side?

(1) `donor_cos` on the real 300.  Production carries THREE donor slots
    (k562_gwps, orion_hct116, orion_hek293t; `make_submission_pm_cd4.donor_weights`), so the key
    generalises to the mean PAIRWISE cosine of the donors' specific-effect vectors.  A target
    carried by one donor has no pairwise cosine at all; it must be given the floor, and this
    prints how many such targets there are.

(2) A DE-side proxy.  Per-gene damping never flips a sign (alpha >= 0) so `fid`'s precision can
    only move through MEMBERSHIP: which genes end up in the top of the call order.  So take the
    top-K genes by |dh| per perturbation and measure the fraction whose sign matches the
    reference.  That is not `fid` -- it is the ordering-driven part of it, which is the part
    rule 1 says transfers.
"""
from __future__ import annotations
import sys
import numpy as np
import pandas as pd
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from vcc2026.donor_transfer import Donor
from decompose import terms
from msekit import group_stats

PROD = ("k562_gwps", "orion_hct116", "orion_hek293t")
targets = sorted(set(pd.read_csv(f"{ROOT}/data/challenge/pert_counts.csv")["target_gene"].astype(str)))
genes = np.asarray(pd.read_csv(f"{ROOT}/data/challenge/gene_names.csv").iloc[:, 0].astype(str))
print(f"### (1) challenge panel: {len(targets)} targets, {genes.size} genes, donors {PROD}")
Ms, hv = {}, {}
for nm in PROD:
    dn = Donor.load(f"{ROOT}/data/donors/{nm}")
    m, s = dn.decompose()
    gi = {g: i for i, g in enumerate(dn.genes)}
    cols = np.array([gi.get(g, -1) for g in genes]); ok = cols >= 0
    pi = {p: i for i, p in enumerate(dn.perts)}
    M = np.zeros((len(targets), genes.size), np.float32); h = np.zeros(len(targets), bool)
    for r, t in enumerate(targets):
        if t in pi:
            M[r, ok] = s[pi[t]][cols[ok]]; h[r] = True
    Ms[nm] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
    hv[nm] = h
    print(f"  {nm:<16} carries {int(h.sum()):>3}/{len(targets)} targets, "
          f"{int(ok.sum()):>6}/{genes.size} genes on the panel axis")
    del dn, s
ndon = sum(hv[n].astype(int) for n in PROD)
print(f"  donors per target: " + ", ".join(f"{k}: {int((ndon==k).sum())}" for k in (0, 1, 2, 3)))
pairs, cnt = np.zeros(len(targets)), np.zeros(len(targets))
for i in range(3):
    for j in range(i + 1, 3):
        a, b = Ms[PROD[i]], Ms[PROD[j]]
        m = hv[PROD[i]] & hv[PROD[j]]
        na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
        c = np.einsum("ij,ij->i", a, b) / np.maximum(na * nb, 1e-12)
        pairs += np.where(m, c, 0.0); cnt += m
dcos = np.where(cnt > 0, pairs / np.maximum(cnt, 1), np.nan)
ok = np.isfinite(dcos)
print(f"  donor_cos computable on {int(ok.sum())}/{len(targets)} targets "
      f"({int((~ok).sum())} get the floor).  range [{np.nanmin(dcos):+.4f}, "
      f"{np.nanmax(dcos):+.4f}] median {np.nanmedian(dcos):+.4f}")
q = np.nanpercentile(dcos, [5, 25, 50, 75, 95])
print(f"  quantiles p5/25/50/75/95: " + " / ".join(f"{v:+.4f}" for v in q))
z = np.zeros(len(targets))
z[ok] = (np.argsort(np.argsort(dcos[ok])).astype(float) + 0.5) / ok.sum()
print(f"  weight z=rank01^2: mean {np.mean(z**2):.4f}, so the global alpha would rise by "
      f"~{1/np.sqrt(np.mean(z**2)):.2f}x to hold the cliff")
np.savez(f"{ROOT}/experiments/097-direction/work/challenge_donor_cos.npz",
         targets=np.array(targets), dcos=dcos, ndon=ndon)

print("\n### (2) DE-side proxy: sign accuracy of the top-K genes by |dh|, base vs arms")
for nm, real, pred in [
    ("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")]:
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    dh, d = t["dh"], t["d"]
    tg = np.array([str(x) for x in pp[t["keep"]]])
    P, G = dh.shape
    Ds, hs = {}, {}
    for nmd in ("k562_gwps", "orion_hek293t"):
        dn = Donor.load(f"{ROOT}/data/donors/{nmd}")
        m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in rg]); okc = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        M = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, x in enumerate(tg):
            if x in pi:
                M[r, okc] = s[pi[x]][cols[okc]]; h[r] = True
        Ds[nmd] = M / (float(np.median(np.linalg.norm(s[:, cols[okc]], axis=1))) or 1.0)
        hs[nmd] = h
        del dn, s
    A_, B_ = Ds["k562_gwps"], Ds["orion_hek293t"]
    both = hs["k562_gwps"] & hs["orion_hek293t"]
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    dc = np.where(both, np.einsum("ij,ij->i", A_, B_) / np.maximum(na * nb, 1e-12), 0.0)
    zz = (np.argsort(np.argsort(dc)).astype(float) + 0.5) / P
    ag = np.sqrt(np.abs(A_)*np.abs(B_)) * ((np.sign(A_)==np.sign(B_)) & (A_!=0) & (B_!=0))
    al = np.empty_like(ag)
    for i in range(P):
        al[i] = (np.argsort(np.argsort(ag[i])).astype(float) + 0.5) / G
    print(f"\n  {nm}: {P} perturbations")
    print(f"  {'arm':<16}" + "".join(f"{'K='+str(K):>10}" for K in (20, 50, 100, 300, 1000)))
    for lbl, Wm in (("base", np.ones_like(dh)), ("pert^2", (zz**2)[:, None]*np.ones_like(dh)),
                    ("gene^1", al), ("stack p2g1", (zz**2)[:, None]*al),
                    ("stack p2g2", (zz**2)[:, None]*(al**2))):
        DH = Wm * dh
        row = []
        for K in (20, 50, 100, 300, 1000):
            acc = []
            for i in range(P):
                idx = np.argpartition(-np.abs(DH[i]), K)[:K]
                m = d[i][idx] != 0
                if m.sum():
                    acc.append(float(np.mean(np.sign(DH[i][idx][m]) == np.sign(d[i][idx][m]))))
            row.append(np.mean(acc))
        print(f"  {lbl:<16}" + "".join(f"{v:>10.4f}" for v in row))
