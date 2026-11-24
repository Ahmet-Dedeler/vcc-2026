"""209 Phase 1 — estimate the subspace the two Orion donors share BEYOND what each shares
with K562, and decide whether it is anything other than the long-gene technical axis.

THE ESTIMATOR, AND EVERY CHOICE IN IT, DECLARED BEFORE THE ENDPOINT IS SEEN
--------------------------------------------------------------------------
frame        `E062.prepare_main_source`'s output -- the object `build_baseline` averages:
             the donor's delta on the 18,533-gene official axis, centred on that bank's own
             gene mean over ALL its rows, scaled by EFFECT_SCALE / raw_median_norm, 300
             panel columns zeroed.  This is the frame a deflation arm would act in.
rows         the **2,000 NON-PANEL** targets carried by all three donors (204's seeded draw,
             seed 20260903).  The 259 panel targets are EXCLUDED, so no fold's scored rows
             enter the estimate and the subspace is leak-free with respect to the panel.
columns      `C` = the genes all three donors measure on those rows (the K562 regression
             below is only defined where K562 can see).  A `wide` variant repeats the
             construction on the union of the two Orion donors' columns, where the K562
             residual is the raw row off C; both are reported and both are scored.
residual     the excess is what survives removing the K562-explained component TARGET BY
             TARGET:   eps_i(t) = d_i(t) - beta_i(t) d_k(t),  beta = <d_i,d_k>/||d_k||^2.
             K562 is the instrument for the shared factor; a per-target slope is the
             strongest available form and needs no global attenuation constant (204 proved
             none is identified).
statistic    S = Sym( eps_h^T eps_e / n ).  A CROSS-covariance between the two donors'
             residuals, so independent per-donor noise contributes only at O(1/sqrt(n))
             and only the SHARED part of the residual accumulates.  Symmetrised because a
             deflation is a projection and must be direction-, not order-, defined.
DIMENSION    **r = #{ i : lambda_i > lambda* }, where lambda* is the 95th percentile of the
             LARGEST eigenvalue of S over B = 200 permutations of the target correspondence
             between eps_h and eps_e.**  A family-wise (Westfall-Young) criterion: it
             controls the chance that ANY null direction is admitted, and it is a property
             of the estimate alone -- it never touches `pds`, `cos`, or any score.
basis        V = the top-r eigenvectors of S (orthonormal by construction), embedded in the
             official axis with zeros off C.

Computed exactly through the row space (rank(eps) <= n << p): stack A = [eps_h; eps_e],
take its economy SVD A = U Sig W^T, and eigendecompose the 2n x 2n matrix
Sym(B_h^T B_e)/n with B = U Sig; the eigenvectors of S are W Z.  The permutation null needs
no re-factorisation at all -- permuting the rows of eps_e permutes the rows of B_e.

FALSIFIER BUILT IN (mandate Phase 1, and the arm is VOID if it fires): if the subspace is
the long-gene / nuclear-fraction technical axis at higher rank, production's whitener
already removes it (0.11% of shipped energy) and there is nothing to take.  Three
independent checks: (a) the overlap of V with each Orion donor's own top PC -- the axis of
the standing harness warning, cos 0.872 between the two donors; (b) the energy share of V
in the equal-weight consensus BEFORE and AFTER the production operator; (c) re-estimating
the whole subspace after deflating that axis first, and asking whether anything survives.
"""
from __future__ import annotations

import json
import os
import resource
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

import prep204  # noqa: E402

prep204.WORK = WORK
from prep204 import prepared, shared_targets  # noqa: E402
from plolo import OFFICIAL, PANEL_SET  # noqa: E402

DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
SEED = 20260903
N_PERM = 200
GAMMA, LAM = 0.65, 0.0


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1 << 30)


# --------------------------------------------------------------------------- estimator
def excess_subspace(Xk, Xh, Xe, n_perm=N_PERM, seed=SEED, kmax=400):
    """S = Sym(eps_h^T eps_e / n); returns (eigvals, eigvecs, permutation threshold)."""
    n, p = Xh.shape
    nk = (Xk * Xk).sum(1)
    live = nk > 0
    bh = np.where(live, (Xh * Xk).sum(1) / np.maximum(nk, 1e-300), 0.0)
    be = np.where(live, (Xe * Xk).sum(1) / np.maximum(nk, 1e-300), 0.0)
    Eh = Xh - bh[:, None] * Xk
    Ee = Xe - be[:, None] * Xk

    A = np.vstack([Eh, Ee])
    U, sig, Wt = np.linalg.svd(A, full_matrices=False)
    keep = sig > sig[0] * 1e-10
    U, sig, Wt = U[:, keep], sig[keep], Wt[keep]
    B = U * sig                                     # A = B @ Wt, rows of B in the W basis
    Bh, Be = B[:n], B[n:]
    K = (Bh.T @ Be) / n
    K = 0.5 * (K + K.T)
    lam, Z = np.linalg.eigh(K)
    order = np.argsort(lam)[::-1]
    lam, Z = lam[order], Z[:, order]

    # ---- permutation null on the target correspondence, in the reduced basis ----------
    rng = np.random.default_rng(seed)
    m = K.shape[0]
    tops = []
    for _ in range(n_perm):
        pi = rng.permutation(n)
        Bep = Be[pi]

        def mv(v):
            return (Bh.T @ (Bep @ v) + Bep.T @ (Bh @ v)) / (2.0 * n)

        # power iteration on the shifted operator to get the largest ALGEBRAIC eigenvalue
        v = rng.standard_normal(m)
        v /= np.linalg.norm(v)
        shift = float(np.abs(lam).max())
        for _ in range(60):
            w = mv(v) + shift * v
            nv = np.linalg.norm(w)
            if nv == 0:
                break
            v = w / nv
        tops.append(float(v @ mv(v)))
    tops = np.asarray(tops)
    thresh = float(np.percentile(tops, 95))
    r = int((lam > thresh).sum())

    Vt = (Wt.T @ Z[:, :max(min(r, kmax), 1)]).T      # (r, p), orthonormal rows
    return dict(eigvals=lam, basis=Vt, r=r, thresh=thresh, null_tops=tops,
                eps=(Eh, Ee), beta=(bh, be))


def energy_share(X, V):
    """Fraction of ||X||_F^2 that lies in span(V) (V has orthonormal rows)."""
    e = float((X * X).sum())
    return float(((X @ V.T) ** 2).sum() / e) if e > 0 else 0.0


def top_pcs(X, k=3):
    Xc = X - X.mean(0, keepdims=True)
    _, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    return Vt[:k], (s ** 2 / (s ** 2).sum())[:k]


def main() -> None:
    tg = shared_targets()
    panel_rows = np.asarray([i for i, t in enumerate(tg) if t in PANEL_SET])
    est_rows = np.asarray([i for i, t in enumerate(tg) if t not in PANEL_SET])
    print(f"targets: {len(tg)} total, {len(panel_rows)} panel (EXCLUDED from the estimate),"
          f" {len(est_rows)} non-panel used for estimation", flush=True)

    P = {}
    for d in DONORS:
        p = prepared(d, tg, "shared")
        P[d] = np.asarray(p["effects"], dtype=np.float64)
    support = {d: (np.abs(P[d]).sum(0) > 0) for d in DONORS}
    C = np.flatnonzero(support["k562_gwps"] & support["orion_hct116"]
                       & support["orion_hek293t"])
    WIDE = np.flatnonzero(support["orion_hct116"] | support["orion_hek293t"])
    print(f"columns: C (all three measure) = {len(C)}; WIDE (either Orion) = {len(WIDE)}; "
          f"rss {rss():.1f}G", flush=True)

    out: dict = {"n_targets_total": len(tg), "n_panel": int(len(panel_rows)),
                 "n_estimation": int(len(est_rows)), "n_cols_C": int(len(C)),
                 "n_cols_WIDE": int(len(WIDE)), "seed": SEED, "n_perm": N_PERM}

    variants = {}
    for tag, cols in (("C", C), ("WIDE", WIDE)):
        Xk = P["k562_gwps"][np.ix_(est_rows, cols)]
        Xh = P["orion_hct116"][np.ix_(est_rows, cols)]
        Xe = P["orion_hek293t"][np.ix_(est_rows, cols)]
        res = excess_subspace(Xk, Xh, Xe)
        lam, V, r = res["eigvals"], res["basis"], res["r"]
        print(f"\n[{tag}] r = {r}   threshold {res['thresh']:.6g}   "
              f"null top: median {np.median(res['null_tops']):.6g} "
              f"max {res['null_tops'].max():.6g}", flush=True)
        print(f"[{tag}] top-12 eigenvalues: "
              + " ".join(f"{v:.4g}" for v in lam[:12]), flush=True)
        print(f"[{tag}] most negative: {lam[-1]:.4g}   sum(+)/sum|.| "
              f"{lam[lam > 0].sum() / np.abs(lam).sum():.3f}", flush=True)

        Vfull = np.zeros((V.shape[0], len(OFFICIAL)))
        Vfull[:, cols] = V
        # energy shares, on the estimation rows and on the PANEL rows (the population
        # production predicts, and which the estimate never saw)
        share = {}
        for d in DONORS:
            share[d] = dict(
                estimation=energy_share(P[d][np.ix_(est_rows, cols)], V),
                panel=energy_share(P[d][np.ix_(panel_rows, cols)], V),
                panel_full_axis=energy_share(P[d][panel_rows], Vfull),
                residual_estimation=None)
        # the same share, but of the RESIDUAL (the object the subspace was fitted to)
        Eh, Ee = res["eps"]
        share["orion_hct116"]["residual_estimation"] = energy_share(Eh, V)
        share["orion_hek293t"]["residual_estimation"] = energy_share(Ee, V)
        for d in DONORS:
            print(f"[{tag}] energy share in {d:15s} est {share[d]['estimation']:.4f}  "
                  f"panel {share[d]['panel']:.4f}  panel(full axis) "
                  f"{share[d]['panel_full_axis']:.4f}", flush=True)

        variants[tag] = dict(r=r, thresh=res["thresh"],
                             eigvals_top50=lam[:50].tolist(),
                             eigvals_bottom5=lam[-5:].tolist(),
                             null_top_median=float(np.median(res["null_tops"])),
                             null_top_max=float(res["null_tops"].max()),
                             energy_share=share)
        np.savez_compressed(WORK / f"subspace-{tag}.npz",
                            basis=V.astype(np.float64), cols=cols,
                            eigvals=lam[:V.shape[0]], r=r, thresh=res["thresh"],
                            null_tops=res["null_tops"])
        del Xk, Xh, Xe, res, V, Vfull
    out["variants"] = variants

    # ---------------------------------------------------------------- the long-gene axis
    print("\n=== the long-gene / nuclear-fraction technical axis ===", flush=True)
    pcs, var = {}, {}
    for d in ("orion_hct116", "orion_hek293t", "k562_gwps"):
        cols = WIDE if d != "k562_gwps" else np.flatnonzero(support["k562_gwps"])
        v, s = top_pcs(P[d][np.ix_(est_rows, cols)], k=3)
        full = np.zeros((3, len(OFFICIAL)))
        full[:, cols] = v
        pcs[d], var[d] = full, s
        print(f"  {d:15s} top-3 variance share {s.round(4).tolist()}", flush=True)
    cos_hh = float(abs(pcs["orion_hct116"][0] @ pcs["orion_hek293t"][0]))
    print(f"  cos(top PC HCT116, top PC HEK293T) = {cos_hh:.4f}   "
          f"(donor-axis value in 175 / the harness warning: 0.872)", flush=True)
    for d in ("orion_hct116", "orion_hek293t"):
        idx = np.argsort(-np.abs(pcs[d][0]))[:25]
        print(f"  {d} top-25 loadings: "
              + ", ".join(str(OFFICIAL[i]) for i in idx), flush=True)
    out["long_gene_axis"] = dict(
        cos_top_pc_hct_hek=cos_hh,
        variance_share={d: var[d].tolist() for d in var},
        top_loadings={d: [str(OFFICIAL[i]) for i in np.argsort(-np.abs(pcs[d][0]))[:30]]
                      for d in ("orion_hct116", "orion_hek293t")})

    # the rank-2 "blob" basis: the two Orion top PCs, orthonormalised
    Bblob = np.vstack([pcs["orion_hct116"][0], pcs["orion_hek293t"][0]])
    Q, _ = np.linalg.qr(Bblob.T)
    Bblob = Q.T
    np.savez_compressed(WORK / "blob-axis.npz", basis=Bblob,
                        pc_hct=pcs["orion_hct116"], pc_hek=pcs["orion_hek293t"],
                        var_hct=var["orion_hct116"], var_hek=var["orion_hek293t"])

    # overlaps: how much of each blob direction lies inside the excess subspace
    ov = {}
    for tag in ("C", "WIDE"):
        z = np.load(WORK / f"subspace-{tag}.npz")
        V = np.zeros((z["basis"].shape[0], len(OFFICIAL)))
        V[:, z["cols"]] = z["basis"]
        ent = {}
        for name, v in (("pc1_hct116", pcs["orion_hct116"][0]),
                        ("pc1_hek293t", pcs["orion_hek293t"][0]),
                        ("blob_rank2", None)):
            if v is None:
                ent[name] = float(((Bblob @ V.T) ** 2).sum() / Bblob.shape[0])
            else:
                ent[name] = float(((V @ v) ** 2).sum())
        # principal angles between span(V) and the rank-2 blob
        s = np.linalg.svd(V @ Bblob.T, compute_uv=False)
        ent["principal_cosines_vs_blob"] = s.tolist()
        # how much of V's own energy the blob explains
        ent["blob_share_of_V"] = float((s ** 2).sum() / V.shape[0])
        ov[tag] = ent
        print(f"  [{tag}] |P_S pc1_hct|^2 {ent['pc1_hct116']:.4f}   "
              f"|P_S pc1_hek|^2 {ent['pc1_hek293t']:.4f}   "
              f"blob share of V {ent['blob_share_of_V']:.4f}", flush=True)
    out["blob_overlap"] = ov

    # ------------------------------------------- the target-independent mean vectors (204 §3.2)
    mv = {}
    for d in DONORS:
        m = P[d][est_rows].mean(0)
        nrm = float(np.linalg.norm(m))
        mv[d] = dict(norm=nrm,
                     share_of_table=float(nrm ** 2 * len(est_rows)
                                          / (P[d][est_rows] ** 2).sum()))
    out["blob_mean_vectors"] = mv
    for tag in ("C", "WIDE"):
        z = np.load(WORK / f"subspace-{tag}.npz")
        V = np.zeros((z["basis"].shape[0], len(OFFICIAL)))
        V[:, z["cols"]] = z["basis"]
        for d in DONORS:
            m = P[d][est_rows].mean(0)
            u = m / max(np.linalg.norm(m), 1e-300)
            mv[d][f"in_subspace_{tag}"] = float(((V @ u) ** 2).sum())
    print(f"  target-independent mean vectors: "
          + "  ".join(f"{d.split('_')[-1]} share {mv[d]['share_of_table']:.5f} "
                      f"inS(C) {mv[d]['in_subspace_C']:.4f}" for d in DONORS), flush=True)

    # ------------------------------------ does the production whitener already remove it?
    print("\n=== does production's operator already remove the subspace? ===", flush=True)
    fold = np.load(WORK / "plolo/fold-kolf2.npz", allow_pickle=True)
    wcols = np.asarray(fold["wcols"])
    pre = np.asarray(fold["pre"], dtype=np.float64)
    with np.load(ROOT / "data/donors/whitener_union_pruned.npz") as z:
        evals = np.asarray(z["evals"], dtype=np.float64)
        evecs = np.asarray(z["evecs"], dtype=np.float64)
    w = (evals + LAM) ** (-GAMMA / 2.0)

    def apply_op(M):
        out_ = np.asarray(M, dtype=np.float64).copy()
        before = float(np.median(np.linalg.norm(out_[:, wcols], axis=1))) or 1.0
        t = ((out_[:, wcols] @ evecs) * w) @ evecs.T
        after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
        out_[:, wcols] = t * (before / after)
        return out_

    whit = apply_op(pre)
    surv = {}
    for tag in ("C", "WIDE"):
        z = np.load(WORK / f"subspace-{tag}.npz")
        V = np.zeros((z["basis"].shape[0], len(OFFICIAL)))
        V[:, z["cols"]] = z["basis"]
        surv[tag] = dict(
            share_in_raw_consensus=energy_share(pre, V),
            share_in_whitened_consensus=energy_share(whit, V),
            n_dim=int(V.shape[0]))
        print(f"  [{tag}] r={V.shape[0]:4d}  energy share in the equal-weight consensus: "
              f"raw {surv[tag]['share_in_raw_consensus']:.4f} -> whitened "
              f"{surv[tag]['share_in_whitened_consensus']:.4f}", flush=True)
    Vb = np.zeros((2, len(OFFICIAL)))
    Vb[:, :] = Bblob
    surv["blob_rank2"] = dict(share_in_raw_consensus=energy_share(pre, Bblob),
                              share_in_whitened_consensus=energy_share(whit, Bblob),
                              n_dim=2)
    print(f"  [blob rank-2] raw {surv['blob_rank2']['share_in_raw_consensus']:.4f} -> "
          f"whitened {surv['blob_rank2']['share_in_whitened_consensus']:.4f}   "
          f"(AGENTS.md: 9.8% raw -> 0.11% whitened for the rank-1 axis)", flush=True)
    # the same for a random subspace of the same dimension, as the scale
    rng = np.random.default_rng(SEED)
    for tag in ("C", "WIDE"):
        r = surv[tag]["n_dim"]
        cols = np.load(WORK / f"subspace-{tag}.npz")["cols"]
        G = rng.standard_normal((min(r, 200), len(cols)))
        Q, _ = np.linalg.qr(G.T)
        R = np.zeros((Q.shape[1], len(OFFICIAL)))
        R[:, cols] = Q.T
        surv[tag]["random_same_dim_raw"] = energy_share(pre, R) * (r / Q.shape[1])
        surv[tag]["random_same_dim_whitened"] = energy_share(whit, R) * (r / Q.shape[1])
        print(f"  [{tag}] random subspace of the same dimension: raw "
              f"{surv[tag]['random_same_dim_raw']:.4f} -> whitened "
              f"{surv[tag]['random_same_dim_whitened']:.4f}", flush=True)
    out["whitener_survival"] = surv

    # --------------------------------- blob-first: is anything left after the axis goes?
    print("\n=== blob control: re-estimate after deflating the long-gene axis ===",
          flush=True)
    Pb = {d: P[d].copy() for d in ("orion_hct116", "orion_hek293t")}
    for d in Pb:
        Pb[d] -= (Pb[d] @ Bblob.T) @ Bblob
    Xk = P["k562_gwps"][np.ix_(est_rows, C)]
    res2 = excess_subspace(Xk, Pb["orion_hct116"][np.ix_(est_rows, C)],
                           Pb["orion_hek293t"][np.ix_(est_rows, C)])
    print(f"  after blob deflation: r = {res2['r']} (was {variants['C']['r']}), "
          f"top-8 eigenvalues " + " ".join(f"{v:.4g}" for v in res2["eigvals"][:8]),
          flush=True)
    lam1 = variants["C"]["eigvals_top50"]
    print(f"  before:                             top-8 eigenvalues "
          + " ".join(f"{v:.4g}" for v in lam1[:8]), flush=True)
    V2 = res2["basis"]
    np.savez_compressed(WORK / "subspace-CblobDeflated.npz", basis=V2, cols=C,
                        eigvals=res2["eigvals"][:V2.shape[0]], r=res2["r"],
                        thresh=res2["thresh"])
    zC = np.load(WORK / "subspace-C.npz")
    ang = np.linalg.svd(V2 @ zC["basis"].T, compute_uv=False)
    out["blob_first"] = dict(
        r=res2["r"], thresh=res2["thresh"],
        eigvals_top20=res2["eigvals"][:20].tolist(),
        eigenvalue_ratio_top1=float(res2["eigvals"][0] / lam1[0]),
        principal_cosines_vs_original=ang[:10].tolist(),
        mean_sq_cosine=float((ang ** 2).mean()))
    print(f"  top-1 eigenvalue ratio after/before = "
          f"{out['blob_first']['eigenvalue_ratio_top1']:.4f}; "
          f"principal cosines vs the original subspace: "
          + " ".join(f"{v:.3f}" for v in ang[:6]), flush=True)

    out["peak_rss_gib"] = rss()
    (WORK / "subspace209.json").write_text(json.dumps(out, indent=2))
    print(f"\npeak RSS {rss():.2f} GiB -> work/subspace209.json")


if __name__ == "__main__":
    main()
