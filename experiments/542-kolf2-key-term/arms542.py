"""542 (PLOLO, cancer folds only) — an EXTERNAL truth proxy in the LODO key: kolf2 (iPSC; not a donor of the composition,
not a held-out line of these folds). Its term t_p = <F_p, T_kolf2,p> / |F_p|^2 (F = the PP composition row the LODO key
already uses, panel columns zeroed; T = kolf2's cached fold truth, 273 targets) is rescaled to the donor terms' scale by
the ratio of medians over the targets both cover, then averaged with the donor terms: key+ = (sum_j v_j + t') / (n_j + 1).
Keys: mean (KA), mean+kolf2, kolf2 alone (information only). Clips (0.5, 3) and (0.5, 6), q 1.75, realised 16.

(541's docstring follows.)
541 (PLOLO) — a robust LODO key: aggregate the per-donor leave-one-out projections by min / median / mean-minus-spread
instead of the mean (415's V1), with the cap at 3 (446bj2), 6 (495) and uncapped. Per target p and carrying donor j:
v_j(p) = <s (b^{-j}_p + radd_p), d_{j,p}> / |F_p|^2 (lib418.lodo_key's terms, reimplemented to keep them per j; KA: their
mean == lodo_key == 418's V1__PP at (0.3, 3)). Keys: mean (V1), min_j, median_j, mean - sd_j/sqrt(n) (a one-se lower
bound), and mean x agreement (fraction of j with v_j > 0). Everything else is 540b's path at q 1.75, realised 16.

(540b's docstring follows.)
540b (PLOLO) — the cap extension of 540 (caps 4 ... uncapped; same code, CLIPS changed).

540 (PLOLO) — the LODO key's clip (floor, cap) on 446bj2's stack at q 1.75, one fold per process (lib420, 420's path).

The fold's LODO-V1 key vector kv1 is computed as comp420 does (lodo_bases + lodo_key on PP; KA: weights(kv1, 0.3, 3) ==
418's saved `V1__PP`), then each clip w = weights(kv1, carried>=2, lo, hi) keys the fold's unkeyed stack composition (fST),
which is sharpened at q 1.75 on each frame's gate (row mode), rendered (minnorm) and realised at KMAX 16 against PP's
unkeyed carriers. KNOWN ANSWER (abort unless |d| < 1e-12): clip (0.5, 3) == comp420's `stack|V1_5_3+sh1.75` (designed and
realised 16).

    uv run python experiments/542-kolf2-key-term/arms542.py <fold>   # cancer folds only        (~5.5 GB peak, ~20 s)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib420  # noqa: E402
from lib420 import (COMPS, E062, FOLDS, PROD_OPS, TS, Donors, Operator, betas_for, build_fold,  # noqa: E402
                    cd4_frozen_and_partial, compose_arm, compose_pp, de_profiles, deflate, frames, key_scale,
                    pre_normrel, profile_from_effect, realise_k, sharpen, source_residual, union_rows, weights)
import lodo415  # noqa: E402
from comp420 import cmp, score_rows  # noqa: E402

W420 = ROOT / "experiments/420-bulk-stack/work"
W418 = ROOT / "experiments/418-bulk-composition/work"
WORK = HERE / "work"
TOL = 1e-12
CLIPS = [(0.5, 3.0), (0.5, 4.0), (0.5, 6.0), (0.5, 8.0), (0.5, 10.0), (0.5, 15.0), (0.5, 1e9), (0.3, 6.0), (0.4, 6.0),
         (0.3, 1e9)]


def tag(lo, hi):
    return f"V1_{lo:g}_{hi:g}"


def run(held: str) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    fold = build_fold(held)
    fold["scored"] = np.asarray(fold["scored"])
    scored = fold["scored"]
    pp_pos = fold["panel_positions"]
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    truth = truth_full[:, scored]
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production"
                  else FOLDS[held]["op"])
    nonpanel = np.ones(truth_full.shape[1], bool)
    nonpanel[pp_pos] = False
    betas = betas_for(held, nonpanel)
    res = dict(fold=held, n_targets=len(targets), ka={}, arms={}, keys={})
    arrays = {}
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in lodo415.UNION])
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    fb = pre228.copy()
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)
    f415, _, _, radd415, _ = compose_pp(op, pre, cd4, partial, h1, wcols, pp_pos)
    whit = op.apply(pre * lib420.ALPHA, wcols, lib420.GAMMA, lib420.LAM)
    fPP, pPP = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    assert np.array_equal(fPP, f415) and np.array_equal(pPP["radd"], radd415), "KA-PP compose"
    mc, gh = COMPS["stack"]
    fST, _ = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, g_h1=gh)
    import lib418 as _l418
    d_rows = {s_: _l418.donor_rows(s_, targets)[0] for s_ in sources}
    A_ = {s_: dn.A[s_] for s_ in sources}
    carried_all = sum(A_[s_].astype(int) for s_ in sources) >= 2
    bases = _l418.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = np.nan_to_num(_l418.lodo_key(fPP, pPP, bases, d_rows, A_, carried_all, sources, pp_pos), nan=0.0)
    assert np.array_equal(weights(kv1, carried_all, 0.3, 3.0)[0], np.load(W418 / f"arrays418-{held}.npz")["V1__PP"]), "KA V1"
    # per-donor terms (lodo_key's arithmetic, cols=None)
    pp = np.einsum("ij,ij->i", fPP, fPP)
    V = np.full((len(sources), len(targets)), np.nan)
    for a, j in enumerate(sources):
        Pj = pPP["s"] * (bases[j] + pPP["radd"])
        Pj[:, pp_pos] = 0.0
        num = np.einsum("ij,ij->i", Pj, d_rows[j])
        V[a] = np.where(A_[j], num / pp, np.nan)
    V[:, ~carried_all] = np.nan
    with np.errstate(all="ignore"):
        kmean = np.nan_to_num(np.nanmean(V, axis=0), nan=0.0)
    dmax = float(np.max(np.abs(kmean[carried_all] - kv1[carried_all]) / np.maximum(np.abs(kv1[carried_all]), 1e-300)))
    log(f"KA per-donor mean == lodo_key: max rel |d| {dmax:.1e}")
    assert dmax < 1e-12, "KA per-donor terms"
    nj = np.sum(~np.isnan(V), axis=0)
    with np.errstate(all="ignore"):
        kmin = np.nan_to_num(np.nanmin(V, axis=0), nan=0.0)
        kmed = np.nan_to_num(np.nanmedian(V, axis=0), nan=0.0)
        ksd = np.nan_to_num(np.nanstd(V, axis=0, ddof=1), nan=0.0)
        kagree = np.nan_to_num(np.nanmean(V > 0, axis=0), nan=0.0)
    zk2 = np.load(ROOT / "experiments/420-bulk-stack/work/plolo/fold-kolf2.npz", allow_pickle=False)
    k2pos = {str(t): i for i, t in enumerate(zk2["targets"])}
    T2 = np.asarray(zk2["truth"], np.float64)
    Fz = fPP.copy()
    Fz[:, pp_pos] = 0.0
    tk = np.full(len(targets), np.nan)
    for i, t in enumerate(targets):
        if t in k2pos and carried_all[i] and pp[i] > 0:
            tk[i] = float(Fz[i] @ T2[k2pos[t]]) / pp[i]
    both = ~np.isnan(tk) & carried_all
    scale = float(np.median(kv1[both]) / np.median(tk[both]))
    tks = tk * scale
    nsum = np.nansum(V, axis=0)
    kplus = np.where(both, (nsum + np.nan_to_num(tks)) / (nj + 1), kv1)
    kalone = np.where(both, tks, kv1)
    res["kolf2_term"] = dict(n_covered=int(both.sum()), scale=scale,
                             spearman_with_lodo=float(np.corrcoef(np.argsort(np.argsort(kv1[both])), np.argsort(np.argsort(tk[both])))[0, 1]))
    log(f"kolf2 term: covers {int(both.sum())} carried targets; scale {scale:.4f}; Spearman with the LODO key {res['kolf2_term']['spearman_with_lodo']:+.3f}")
    KEYS = {"mean": kv1, "mean+kolf2": kplus, "kolf2 alone": kalone}
    res["n_donors_hist"] = np.bincount(nj[carried_all]).tolist()
    del bases, d_rows, whit, f415, pre, pre228, fb, dn
    P_DE = de_profiles(fPP)
    log("KA-PP and KA V1 (== 418's V1__PP) PASS")

    def rows_for(bulk):
        dr = {}
        for ctx, (b, gate) in frames().items():
            F = sharpen(bulk, 1.75, gate, mode="row")
            des, rea = np.empty_like(F), np.empty_like(F)
            K = np.empty(F.shape[0], int)
            for i in range(F.shape[0]):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p_bk) - b
                p_re, _, K[i], _, _ = realise_k(P_DE[ctx][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            dr[ctx] = dict(design=des, realised=rea, K=K)
        return dr

    c420 = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]
    b540 = json.loads((ROOT / f"experiments/540-lodo-clip-scan/work/arms540b-{held}.json").read_text())["arms"]
    for kname, kv in KEYS.items():
        for lo, hi in ((0.5, 3.0), (0.5, 6.0)):
            W = weights(kv, carried_all, lo, hi)[0]
            dr = rows_for(key_scale(fST, W)[0])
            sr, ar = score_rows({c: dr[c]["realised"] for c in "ABC"}, scored, truth, betas, nonpanel)
            K = np.concatenate([dr[c]["K"] for c in "ABC"])
            sr["K_mean_emul"] = float(K.mean())
            name = f"{kname}|{lo:g}_{hi:g}"
            res["arms"][name] = dict(realised={"16": sr})
            res["keys"][name] = dict(frac_floor=float(np.mean(carried_all & (W <= lo + 1e-12))),
                                     frac_cap=float(np.mean(carried_all & (W >= hi - 1e-12))), w_max=float(W.max()))
            arrays[f"{name}@R16"] = ar
            log(f"{name:18s} R16 cos {sr['cos_agg']:.6f} pds {sr['pds']:.6f} pds@.5 {sr['pds_lam0.5']:.6f} K {K.mean():.2f}")
            if kname == "mean":
                ref = c420["stack|V1_5_3+sh1.75"]["realised"]["16"] if hi == 3.0 else b540[f"V1_{lo:g}_{hi:g}"]["realised"]["16"]
                dd = cmp(sr, ref)
                log(f"KA mean|{lo:g}_{hi:g}: realised |d| {dd:.1e}")
                assert dd < TOL, f"KA FAILED {held} {name}: {dd}"
                res["ka"][name] = dd
    res["peak_rss_gib"] = E062.peak_rss_gib()
    (WORK / f"arms542-{held}.json").write_text(json.dumps(res, indent=1, default=float))
    np.savez_compressed(WORK / f"arrays542-{held}.npz", names=np.asarray(list(arrays)), d=np.asarray((truth * truth).sum(1)),
                        **{f"{nm}__{k}": v for nm, a in arrays.items() for k, v in a.items()})
    log(f"done; peak RSS {res['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    assert sys.argv[1] != "kolf2", "kolf2 is the proxy; it cannot be the held-out fold"
    run(sys.argv[1])
