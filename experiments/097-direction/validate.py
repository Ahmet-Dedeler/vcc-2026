"""Validate the stacked donor key: permutation null, cross-surface exponent choice, bootstrap.

Nothing here is allowed to pick an exponent on the surface it is scored on.  The reported
headline for every arm is the CROSS-SURFACE number -- exponent fitted on the other harness --
and per RESEARCH-METHOD rule 1 the transferable estimate is the smaller of the two.
"""
from __future__ import annotations
import sys
import numpy as np
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/src")
from decompose import terms
from msekit import group_stats
from vcc2026.donor_transfer import Donor

HARNESS = [
    ("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    ("orion", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
     f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
]
DONORS = ("k562_gwps", "orion_hek293t")
GZ = [0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0]
GE = [0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]


def prep(real, pred):
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0; t["d"][k, c] = 0.0
    tg = np.array([str(x) for x in pp[t["keep"]]]); genes = rg
    P, G = t["dh"].shape
    Ms, hv = {}, {}
    for nmd in DONORS:
        dn = Donor.load(f"{ROOT}/data/donors/{nmd}")
        m, s = dn.decompose()
        gi = {g: i for i, g in enumerate(dn.genes)}
        cols = np.array([gi.get(str(g), -1) for g in genes]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(dn.perts)}
        M = np.zeros((P, G)); h = np.zeros(P, bool)
        for r, x in enumerate(tg):
            if x in pi:
                M[r, ok] = s[pi[x]][cols[ok]]; h[r] = True
        Ms[nmd] = M / (float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0)
        hv[nmd] = h
        del dn, s
    A_, B_ = Ms[DONORS[0]], Ms[DONORS[1]]
    both = hv[DONORS[0]] & hv[DONORS[1]]
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    dcos = np.where(both, np.einsum("ij,ij->i", A_, B_) / np.maximum(na * nb, 1e-12), 0.0)
    agree = np.sqrt(np.abs(A_) * np.abs(B_)) * ((np.sign(A_) == np.sign(B_)) & (A_ != 0) & (B_ != 0))
    al = np.empty_like(agree)
    for i in range(P):
        al[i] = (np.argsort(np.argsort(agree[i])).astype(float) + 0.5) / G
    return dict(dh=t["dh"], d=t["d"], w=t["w"], dcos=dcos, al=al, n=P,
                base=None, X=t["X"], A=t["A"], sumD=t["sumD"])


def W(DH, d, w):
    n = float((np.einsum("ij,ij->i", DH, d) * w).sum())
    q = float((np.einsum("ij,ij->i", DH, DH) * w).sum())
    return (n * n / q) if q > 0 else 0.0


def rank01(s):
    return (np.argsort(np.argsort(s)).astype(float) + 0.5) / len(s)


H = {}
for nm, r, p in HARNESS:
    h = prep(r, p)
    h["base"] = W(h["dh"], h["d"], h["w"])
    h["z"] = rank01(h["dcos"])
    H[nm] = h

print("### grid of sqrt(W/W_base) over (gamma on the per-PERT key, eta on the per-GENE key)")
grids = {}
for nm in H:
    h = H[nm]
    Gd = np.zeros((len(GZ), len(GE)))
    for i, gz in enumerate(GZ):
        zz = (h["z"] ** gz)[:, None]
        for j, ge in enumerate(GE):
            Gd[i, j] = np.sqrt(W(zz * (h["al"] ** ge) * h["dh"], h["d"], h["w"]) / h["base"])
    grids[nm] = Gd
    print(f"\n  {nm}   rows gamma={GZ}  cols eta={GE}")
    print("        " + "".join(f"{e:>8.2f}" for e in GE))
    for i, gz in enumerate(GZ):
        print(f"  g={gz:<4.1f}" + "".join(f"{v:>8.3f}" for v in Gd[i]))

print("\n### CROSS-SURFACE: exponents chosen on the OTHER harness")
a, b = "h1-test", "orion"
ia, ja = np.unravel_index(np.argmax(grids[a]), grids[a].shape)
ib, jb = np.unravel_index(np.argmax(grids[b]), grids[b].shape)
print(f"  argmax on {a}: gamma={GZ[ia]} eta={GE[ja]} -> {grids[a][ia,ja]:.3f}x in-sample, "
      f"{grids[b][ia,ja]:.3f}x on {b}")
print(f"  argmax on {b}: gamma={GZ[ib]} eta={GE[jb]} -> {grids[b][ib,jb]:.3f}x in-sample, "
      f"{grids[a][ib,jb]:.3f}x on {a}")
# the single setting that maximises the MINIMUM of the two -- the rule-1 estimator
mn = np.minimum(grids[a], grids[b])
im, jm = np.unravel_index(np.argmax(mn), mn.shape)
print(f"  max-min setting: gamma={GZ[im]} eta={GE[jm]} -> {grids[a][im,jm]:.3f}x / "
      f"{grids[b][im,jm]:.3f}x   transferable {mn[im,jm]:.3f}x")

print("\n### permutation null on the FULL procedure (both keys shuffled, argmax over the grid)")
NP_ = 400
for nm in H:
    h = H[nm]
    obs = grids[nm].max()
    rng = np.random.default_rng(2026)
    null = np.empty(NP_)
    for s in range(NP_):
        pz = rng.permutation(h["z"])
        pa = h["al"][rng.permutation(h["n"])]        # shuffle the per-gene key across rows
        best = 0.0
        for gz in GZ:
            zz = (pz ** gz)[:, None]
            for ge in GE:
                best = max(best, np.sqrt(W(zz * (pa ** ge) * h["dh"], h["d"], h["w"]) / h["base"]))
        null[s] = best
    print(f"  {nm:<10} observed {obs:.3f}   null p50 {np.percentile(null,50):.3f}  "
          f"p95 {np.percentile(null,95):.3f}  p99 {np.percentile(null,99):.3f}  "
          f"p-value {float((null>=obs).mean()):.4f}")

print("\n### bootstrap over perturbations at the FIXED max-min setting "
      f"(gamma={GZ[im]}, eta={GE[jm]})")
for nm in H:
    h = H[nm]
    z = (h["z"] ** GZ[im])[:, None]; A2 = h["al"] ** GE[jm]
    n = h["n"]; rng = np.random.default_rng(11)
    bs = np.empty(1500)
    for i in range(1500):
        k = rng.integers(0, n, n)
        bs[i] = np.sqrt(W((z * A2 * h["dh"])[k], h["d"][k], h["w"][k])
                        / max(W(h["dh"][k], h["d"][k], h["w"][k]), 1e-30))
    pt = np.sqrt(W(z * A2 * h["dh"], h["d"], h["w"]) / h["base"])
    print(f"  {nm:<10} point {pt:.3f}   boot p5 {np.percentile(bs,5):.3f}  "
          f"p50 {np.percentile(bs,50):.3f}  p95 {np.percentile(bs,95):.3f}  "
          f"P(>1) {np.mean(bs>1):.3f}")
    np.save(f"{ROOT}/experiments/097-direction/work/{nm}-boot.npy", bs)
