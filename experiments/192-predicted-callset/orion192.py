"""192 — the CROSS-LINE fit the mandate asks for: fit the membership predictor on the Orion
surfaces and score it on h1-test, and vice versa.

The h1-test LOTO fit is leave-one-TARGET-out but not leave-one-LINE-out: its three
coefficients are still fitted on H1's own reference DE calls.  On the challenge we have no
reference at all, so the coefficients must come from a DIFFERENT cell line.  This is the
falsifier for that step.

Caveats stated up front (AGENTS.md harness warning): the Orion THINNED surfaces carry a
target-independent technical blob (26.4% of HCT116's effect energy, 44.1% of HEK293T's) and
their absolute DE members are not comparable to the board.  Nothing absolute is quoted from
them here -- only whether COEFFICIENTS FITTED THERE reproduce the h1-test enrichment.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from debench import Surface                                            # noqa: E402
from enrich192 import NS, curve                                        # noqa: E402
from frontier192 import A15, INV, PowerSurface, order_key, topk_call   # noqa: E402
from keys192 import gene_stats                                         # noqa: E402

SHRINK = 0.4277


def orion_ctrl_stats(path, cache):
    """Per-gene control CPM / detection rate / cv2 from a harness's own control cells."""
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        return {k: z[k] for k in z.files}
    import anndata as ad
    from scipy.sparse import csr_matrix
    A = ad.read_h5ad(path, backed="r")
    genes = np.asarray([str(g) for g in A.var_names])
    idx = np.flatnonzero(np.asarray(A.obs["target_gene"].astype(str).values) == "non-targeting")
    G = genes.size
    n = idx.size
    s1 = np.zeros(G)
    s2 = np.zeros(G)
    nzc = np.zeros(G)
    for a in range(0, n, 2000):
        blk = A[idx[a:a + 2000]].to_memory().X
        blk = csr_matrix(blk) if not hasattr(blk, "toarray") else blk
        d = np.asarray(blk.todense(), dtype=np.float64)
        lib = d.sum(1, keepdims=True)
        lib[lib == 0] = 1.0
        c = d * (1e6 / lib)
        s1 += c.sum(0)
        s2 += (c ** 2).sum(0)
        nzc += (c > 0).sum(0)
        del d, c
    A.file.close()
    mu = s1 / n
    var = np.maximum(s2 / n - mu ** 2, 0.0)
    out = dict(genes=genes, cpm=mu, nz=nzc / n, cv2=var / np.maximum(mu ** 2, 1e-30))
    np.savez_compressed(cache, **out)
    return out


def feats(cpm, nz, cv2, shape):
    return np.stack([np.broadcast_to(np.log(cpm + 1)[None, :], shape),
                     np.broadcast_to(np.log(nz + 1e-3)[None, :], shape),
                     np.broadcast_to(np.log(cv2 + 1e-3)[None, :], shape)], -1)


def main() -> None:
    from sklearn.linear_model import LogisticRegression
    # ---- Orion HCT116 (thinned): its own reference DE calls and its own control statistics
    ost = orion_ctrl_stats(os.path.join(REPO, "data/harness/orion-hct116-thinned.h5ad"),
                           os.path.join(WORK, "orion_ctrl.npz"))
    oz = np.load(os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz"),
                 allow_pickle=True)
    og = np.asarray([str(g) for g in oz["genes"]])
    p = {str(g): i for i, g in enumerate(ost["genes"])}
    sel = np.array([p[g] for g in og])
    osig = (oz["real_padj"] < 0.05) & np.isfinite(oz["real_padj"])
    Xo = feats(ost["cpm"][sel], ost["nz"][sel], ost["cv2"][sel], osig.shape)
    mo, so = Xo.reshape(-1, 3).mean(0), Xo.reshape(-1, 3).std(0)
    rng = np.random.default_rng(11)
    Xf = ((Xo - mo) / so).reshape(-1, 3)
    yf = osig.reshape(-1)
    idx = rng.choice(Xf.shape[0], size=min(400_000, Xf.shape[0]), replace=False)
    mod = LogisticRegression(max_iter=400).fit(Xf[idx], yf[idx])
    print(f"Orion HCT116-thinned: {osig.shape[0]} perturbations x {osig.shape[1]} genes, "
          f"n_real median {np.median(osig.sum(1)):.0f}")
    print(f"  coefficients (standardised): log cpm {mod.coef_[0][0]:+.3f}, log nz "
          f"{mod.coef_[0][1]:+.3f}, log cv2 {mod.coef_[0][2]:+.3f}")

    # ---- score those coefficients on h1-test, using H1's OWN control statistics
    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    pos = {str(g): i for i, g in enumerate(E["genes"])}
    tpos = {str(t): i for i, t in enumerate(E["targets"])}
    s = Surface(A15)
    gi = np.array([pos[g] for g in s.genes])
    ti = np.array([tpos[t] for t in s.perts])
    cpm = E["cpm"][gi].astype(np.float64)
    absH = np.abs(E["H"][np.ix_(ti, gi)].astype(np.float64))
    gs = gene_stats(os.path.join(E185, "work/ctrl_cs.npy"), os.path.join(WORK, "genestats.npz"))
    sig = s.real_sig & np.isfinite(s.real_lfc)
    nreal = sig.sum(1)
    Xh = feats(cpm, gs["nz"][:], gs["cv2"][:], absH.shape)
    # the transfer is of the COEFFICIENTS; each line standardises with its own moments
    mh, sh = Xh.reshape(-1, 3).mean(0), Xh.reshape(-1, 3).std(0)
    key_x = ((Xh - mh) / sh) @ mod.coef_[0]

    # the h1-test-native LOTO key, for comparison
    Xz = (Xh - mh) / sh
    rngf = np.random.default_rng(7)
    key_n = np.zeros(absH.shape)
    for q in range(absH.shape[0]):
        tr = np.setdiff1d(np.arange(absH.shape[0]), [q])
        Xt = Xz[tr].reshape(-1, 3)
        yt = sig[tr].reshape(-1)
        ii = rngf.choice(Xt.shape[0], size=min(300_000, Xt.shape[0]), replace=False)
        key_n[q] = LogisticRegression(max_iter=400).fit(Xt[ii], yt[ii]).decision_function(Xz[q])

    print("\nFOLD-ENRICHMENT on h1-test — same-line LOTO fit vs CROSS-LINE (Orion) coefficients")
    print(f"  {'key':<34}" + "".join(f"{n:>9}" for n in NS) + f"{'n_real_p':>10}")
    keys = {"FIT LOTO (h1-test's own calls)": key_n,
            "FIT transferred from Orion HCT116": key_x,
            "cpm alone (no fit at all)": np.broadcast_to(cpm[None, :], absH.shape).copy()}
    res = {}
    for lab, kv in keys.items():
        c = curve(kv, sig, nreal)
        res[lab] = c
        print(f"  {lab:<34}" + "".join(f"{c[n]['fold']:>9.2f}" for n in NS)
              + f"{c['nreal']['fold']:>10.2f}")

    print("\nPRICED, at the retrodicted band f = 0.4 / 0.5, against a15's challenge point")
    for f in (0.4, 0.5):
        sf = PowerSurface(A15, f)
        a15c = topk_call(-np.log(np.maximum(np.load(A15, allow_pickle=True)["pred_p"], 1e-300)),
                         np.full(absH.shape[0], 4474))
        a0 = dict(fid=sf.fid(call=a15c), jac=sf.jac(call=a15c),
                  reach=sf.reach(key=order_key(a15c, absH)))
        for lab, kv in keys.items():
            for K in (2200, 2800):
                call = topk_call(kv, np.full(absH.shape[0], K))
                v = dict(fid=sf.fid(call=call), jac=sf.jac(call=call),
                         reach=sf.reach(key=order_key(call, absH)))
                d = {m: (v[m] - a0[m]) * INV[m] * (SHRINK if m == "reach" else 1) / 6 for m in v}
                res[f"{f}|{lab}|{K}"] = dict(raw=v, overall=d, sum3=sum(d.values()))
                print(f"  f={f} K={K:<5} {lab:<34} fid {v['fid']:.5f} jac {v['jac']:.5f} "
                      f"SUM(3) {sum(d.values()):+.5f}")

    with open(os.path.join(WORK, "orion192.json"), "w") as fh:
        json.dump(dict(orion_coef=mod.coef_[0].tolist(), rows=res), fh, indent=1, default=float)
    print(f"\nwrote {WORK}/orion192.json")


if __name__ == "__main__":
    main()
