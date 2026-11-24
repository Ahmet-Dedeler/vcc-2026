"""Do the per-PERTURBATION key and the per-GENE key stack, or are they one signal seen twice?

096 §5 measured a per-GENE donor-agreement shrinkage at 1.381x / 1.172x on sqrt(W).
097 measures a per-PERTURBATION donor-cosine weight at 1.241x / 1.206x.
Both are built from the same two donor screens, so the natural worry is that they are the same
signal counted twice and the product is not available.  Measured here, on both harnesses:

    dh'_pg = z_p * (alpha_pg ** eta) * dh_pg ,   z_p = rank(donor_cos_p)^gamma

reported as sqrt(W/W_base) for each factor alone and for the product.  Also carries:
 * the noise-corrected common mode (route c), because the plug-in Sw*||dbar||^2 contains the
   reference's own sampling noise K/P + Kc and orion has only 150 cells per perturbation;
 * the ON-TARGET check (route d): the metric drops each row's own target coordinate, so the
   knockdown that transfers at rho=0.647 is worth EXACTLY zero here.  Verified numerically.
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


def load(real, pred):
    rp, rs, rc, rj, rg = group_stats(real)
    pp, ps, pc, pj, pg = group_stats(pred)
    t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
    t["dh_raw"] = t["dh"].copy(); t["d_raw"] = t["d"].copy()
    for k, c in t["cols"].items():
        t["dh"][k, c] = 0.0
        t["d"][k, c] = 0.0
    t["genes"] = rg
    t["targets"] = np.array([str(x) for x in pp[t["keep"]]])
    return t


def donor_mats(genes, targets):
    out = {}
    for name in DONORS:
        d = Donor.load(f"{ROOT}/data/donors/{name}")
        m, s = d.decompose()
        gi = {g: i for i, g in enumerate(d.genes)}
        cols = np.array([gi.get(str(g), -1) for g in genes]); ok = cols >= 0
        pi = {p: i for i, p in enumerate(d.perts)}
        M = np.zeros((len(targets), genes.size))
        hv = np.zeros(len(targets), bool)
        for r, tg in enumerate(targets):
            if tg in pi:
                M[r, ok] = s[pi[tg]][cols[ok]]; hv[r] = True
        scale = float(np.median(np.linalg.norm(s[:, cols[ok]], axis=1))) or 1.0
        out[name] = (M / scale, hv)
        del d, s
    return out


def W_of(DH, d, w):
    n = float((np.einsum("ij,ij->i", DH, d) * w).sum())
    q = float((np.einsum("ij,ij->i", DH, DH) * w).sum())
    return (n * n / q) if q > 0 else 0.0


def rank01(s):
    return (np.argsort(np.argsort(s)).astype(float) + 0.5) / len(s)


for name, real, pred in HARNESS:
    t = load(real, pred)
    dh, d, w, genes, tg = t["dh"], t["d"], t["w"], t["genes"], t["targets"]
    P, G = dh.shape
    base = W_of(dh, d, w)
    Dm = donor_mats(genes, tg)
    A_, hA = Dm[DONORS[0]]; B_, hB = Dm[DONORS[1]]
    both = hA & hB
    na, nb = np.linalg.norm(A_, axis=1), np.linalg.norm(B_, axis=1)
    dcos = np.where(both, np.einsum("ij,ij->i", A_, B_) / np.maximum(na * nb, 1e-12), 0.0)
    gate = (np.sign(A_) == np.sign(B_)) & (A_ != 0) & (B_ != 0)
    agree = np.sqrt(np.abs(A_) * np.abs(B_)) * gate           # 096 §5's per-gene alpha source
    # per-gene alpha in [0,1]: rank within each perturbation row
    al = np.zeros_like(agree)
    for i in range(P):
        r = agree[i]
        al[i] = (np.argsort(np.argsort(r)).astype(float) + 0.5) / G
    print(f"\n########## {name}  n={P}  base W={base:.6g} ##########")
    print(f"  mechanism check: donor_cos range [{dcos.min():+.3f}, {dcos.max():+.3f}] "
          f"median {np.median(dcos):+.3f}; rows with 2 donors {int(both.sum())}/{P}; "
          f"per-gene alpha mean {al.mean():.3f}")
    zc = rank01(dcos)
    print(f"\n  {'arm':<44}{'sqrt(W/W_base)':>16}")
    for lbl, DH in (
        ("per-PERT only   z=rank(donor_cos)^2", zc[:, None] ** 2 * dh),
        ("per-PERT only   z=rank(donor_cos)^3", zc[:, None] ** 3 * dh),
        ("per-GENE only   alpha^0.5", al ** 0.5 * dh),
        ("per-GENE only   alpha^1", al * dh),
        ("per-GENE only   alpha^2", al ** 2 * dh),
        ("STACK  z^2 * alpha^0.5", (zc[:, None] ** 2) * (al ** 0.5) * dh),
        ("STACK  z^2 * alpha^1", (zc[:, None] ** 2) * al * dh),
        ("STACK  z^3 * alpha^1", (zc[:, None] ** 3) * al * dh),
    ):
        print(f"  {lbl:<44}{np.sqrt(W_of(DH, d, w)/base):>15.3f}x")

    # ---------- route (c): noise-corrected common mode -----------------------------------
    Sw = float(w.sum())
    dbar = (w[:, None] * d).sum(0) / Sw
    plug = Sw * float(dbar @ dbar)
    K, Kc = t["K"], t["Kc"]
    noise = K / P + Kc                     # trVar(dbar) in the same w-weighted units
    print(f"\n  route (c) common mode: plug-in Sw||dbar||^2 = {plug:.6g}  "
          f"noise floor K/P+Kc = {noise:.6g}  ->  corrected {plug-noise:.6g} "
          f"= {(plug-noise)/t['sumD']:+.4f} raw mse   (plug-in said {plug/t['sumD']:+.4f})")

    # ---------- route (d): the on-target coordinate --------------------------------------
    dhr, dr = t["dh_raw"], t["d_raw"]
    ot_x = sum(w[k] * dhr[k, c] * dr[k, c] for k, c in t["cols"].items())
    ot_a = sum(w[k] * dhr[k, c] ** 2 for k, c in t["cols"].items())
    ot_r = sum(w[k] * dr[k, c] ** 2 for k, c in t["cols"].items())
    print(f"  route (d) on-target: rows whose own target is on the axis "
          f"{len(t['cols'])}/{P}.  If the metric COUNTED it, X would gain {ot_x:+.6g} "
          f"({100*ot_x/t['X']:+.1f} % of X), A {ot_a:+.6g}, R {ot_r:+.6g}. "
          f"The scorer excludes it, so the realised contribution is 0.")
