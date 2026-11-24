"""418 — shared machinery: the SHIPPED bulk composition ("PP", 415's reconstruction) with every composition constant
exposed as a one-variable knob, the shipped per-target key, 417a's bulk sharpening, and 414's read-back instruments.

Composition (sens415.compose_pp, builder order):
    pre  = normrel-weighted donor mean (cap, exp -1.5) + PM50 correction at the uniform K562 share   (knob: normrel)
    W    = Operator(gamma 0.65)(ALPHA * pre)  (on-axis preserve c)
    add  = f_cd4 * CD4_deflated(k=5)  +  f_h1 * H1union  +  f_part * CD4partial(g 0.5)                 (knobs a-d)
           f_cd4 = CD4_F(0.34) * m_cd4 / per_unit(all rows);  f_h1, f_part solved on carrying rows
    F    = s * (W + add)  (one post-add median preserve s),  panel columns zeroed
Key (held fixed; the 410 formula is composition-free): F -> F * w_p * c  (c restores the non-panel Frobenius norm),
    w = clip(key/median, 0.3, 3) of `key410.shared_energy_key` on the fold's source donors (the builder's `cols=None`).
Sharpening (417a): sign(h)|h|^1.25 on each frame's CPM gate, the row's gated norm restored (169 `sharpen(mode="row")`
    on the panel-zeroed row == 417's `bulk_sharpen` with panel columns protected).
Read-back (414): plolo.deliver(minnorm) in the challenge frames A/B/C -> pds, cos_agg, render-space A_c/B_c, and the
    board-emulated pds (414's artifact beta at lam 0.5, rescaled exactly as board414/extra414 do).

Nothing here writes outside experiments/418-bulk-composition/work.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
for p in ("experiments/162-retrodicting-validator", "experiments/218-independent-donors",
          "experiments/228-theta-gain-plolo", "experiments/410-pertscale-revisit",
          "experiments/169-sharpening-six-member", "experiments/415-lodo-key", "experiments/414-bulk-reprice"):
    sys.path.insert(0, str(ROOT / p))

import plolo  # noqa: E402
import arms414  # noqa: E402  (imports 228/410/169; repoints plolo.WORK -- reset below)
import sens415  # noqa: E402  (repoints plolo.WORK to 415/work -- reset below; its own caches are read-only)
from arms414 import frames, key_scale  # noqa: E402
from lodo415 import ALPHA, CD4_F, GAMMA, LAM, UNION, UNION_G, whiten_split  # noqa: E402
from sens415 import DEFLATE_K, PARTIAL_G, Donors, cd4_frozen_and_partial, compose_pp, deflate  # noqa: E402
from compose218 import solve_fraction_carrying, union_rows  # noqa: E402
from resid218 import source_residual  # noqa: E402
from key410 import raw_whiten, shared_energy_key, weights  # noqa: E402
from size410 import donor_rows  # noqa: E402
from sharp import sharpen  # noqa: E402

plolo.WORK = WORK / "plolo"          # cloned (cp -c) fold caches + challenge b_ctrl; nothing else is read or written here
from plolo import E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, endpoints, solve_fraction  # noqa: E402

W414 = ROOT / "experiments/414-bulk-reprice/work"
W415 = ROOT / "experiments/415-lodo-key/work"
NORMREL_EXP = -1.5
SHARP_Q = 1.25
LAMS = (0.5, 1.0)


# ------------------------------------------------------------------------------------------------ donor base
def pre_normrel(dn: Donors, subset, mode: str | None, exp: float | None, fallback_rows=None, pm: str = "prod"):
    """`Donors.pre` with the normrel MODE and EXPONENT exposed (sens415 hard-codes cap at -1.5).

    mode None (or exp None) = no per-target weighting (the builder with VCC_NORMREL_EXP unset): the uniform mean, PM50
    correction at the uniform K562 share.  'cap': w = max(1, |row|/median)^exp.  'plain': w = (|row|/median)^exp.
    The weighted row is norm-matched to the uniform row (the builder's shape-only rule), then the PM correction is added.
    """
    n, G = next(iter(dn.F.values())).shape
    out = np.zeros((n, G))
    weighted = mode is not None and exp is not None
    for t in range(n):
        have = [s for s in subset if dn.A[s][t]]
        if not have:
            if fallback_rows is not None:
                out[t] = fallback_rows[t]
            continue
        rows = np.stack([dn.F[s][t] for s in have])
        u = rows.mean(0)
        if weighted and len(have) > 1:
            r = [float(np.linalg.norm(dn.F[s][t])) / dn.med[s] for s in have]
            if mode == "cap":
                r = [max(1.0, x) for x in r]
            w = np.asarray([x ** exp for x in r])
            cand = (w[:, None] * rows).sum(0) / w.sum()
            nu, nc = float(np.linalg.norm(u)), float(np.linalg.norm(cand))
            if nu > 0 and nc > 0:
                cand = cand * (nu / nc)
        else:
            cand = u
        if "k562_gwps" in have:
            cand = cand + dn.pmcorr[pm][t] / len(have)
        out[t] = cand
    return out


# ------------------------------------------------------------------------------------------------ composition
def compose_arm(whitened, cd4, partial, h1, wcols, panel_positions, *, m_cd4=1.0, g_h1=UNION_G, g_part=PARTIAL_G):
    """compose_pp's residual stage with the knobs exposed; bit-identical to it at the defaults (asserted by callers).

    m_cd4 multiplies CD4's target f (0.34), which multiplies the solved fraction exactly (the solve is linear).
    g_h1 / g_part = 0 means the component is not added at all (the builder with the env var unset).
    Returns (final with panel zeroed, parts) where parts carries W, the merged residual, s and the fractions.
    """
    f_cd4, _ = solve_fraction(whitened, cd4, CD4_F * m_cd4)
    merged = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0) * f_cd4
    support = np.asarray(cd4.gene_votes) > 0
    fr = {"cd4": float(f_cd4), "h1_union": 0.0, "cd4_partial": 0.0}
    for tag, res, g in (("h1_union", h1, g_h1), ("cd4_partial", partial, g_part)):
        if not g:
            continue
        f, _ = solve_fraction_carrying(whitened, res, g)
        v = np.asarray(res.gene_votes) > 0
        merged += np.where(v, np.asarray(res.mean_effect), 0.0) * f
        support |= v
        fr[tag] = float(f)
    added = E062.CD4Prepared(mode="merged", mean_effect=merged, gene_votes=support.astype(np.uint8),
                             target_condition_count=np.asarray(cd4.target_condition_count), condition_details=[])
    final, info = E062.add_post_whiten_residual(whitened, added, 1.0, wcols)
    final[:, panel_positions] = 0.0
    radd = np.where(support, merged, 0.0)
    return final, dict(radd=radd, s=float(info["post_add_global_preserve_scale"]), fractions=fr)


def key_on_base_only(whitened, parts, w, panel_positions, frob_ref):
    """f1: the per-target key scales the donor base W only; the CD4/H1/partial residual is added unscaled; one
    constant then restores the non-panel Frobenius norm of the unkeyed final row (apply_pert_scale's rule)."""
    K = parts["s"] * (w[:, None] * whitened + parts["radd"])
    K[:, panel_positions] = 0.0
    return K * (np.linalg.norm(frob_ref) / np.linalg.norm(K))


# ------------------------------------------------------------------------------------------------ LODO key (415 V1)
def lodo_bases(op, dn, sources, mode, exp, pre_full, wcols):
    """Per left-out donor j: the whitened base rebuilt without j, with the on-axis preserve c held at the full row's
    value (sens415.loo_rows before the residual add).  Depends only on the normrel setting, not on knobs a-d."""
    c, _ = whiten_split(op, pre_full * ALPHA, wcols)
    out = {}
    for j in sources:
        pj = pre_normrel(dn, [i for i in sources if i != j], mode, exp, fallback_rows=None)
        bj = pj * ALPHA
        bj[:, wcols] = c * raw_whiten(op, bj, wcols)[:, wcols]
        out[j] = bj
    return out


def lodo_key(final, parts, bases, d, A, carried, sources, panel_positions, cols=None):
    """key_p = mean_j <s (b^{-j}_p + radd_p), d_{j,p}> / |F_p|^2 over carrying donors (415 V1). cols=None is the
    builder's convention (every column); cols=scored reproduces sens415's PLOLO convention (its KA)."""
    sel = (lambda M: M) if cols is None else (lambda M: M[:, cols])
    Fc = sel(final)
    pp = np.einsum("ij,ij->i", Fc, Fc)
    num = {}
    for j in sources:
        Pj = parts["s"] * (bases[j] + parts["radd"])
        Pj[:, panel_positions] = 0.0
        num[j] = np.einsum("ij,ij->i", sel(Pj), sel(d[j]))
    key = np.full(final.shape[0], np.nan)
    for p in np.flatnonzero(carried):
        key[p] = float(np.mean([num[j][p] for j in sources if A[j][p]])) / pp[p]
    return key


# ------------------------------------------------------------------------------------------------ read-back
def betas_for(held, nonpanel):
    """414's artifact beta per frame, rescaled to the 228 base's PLOLO row energy exactly as board414/extra414."""
    ref = json.loads((W414 / f"arms414-{held}.json").read_text())["arms"]
    bz = np.load(W414 / "beta414.npz")
    dec = json.loads((W414 / "decomp414.json").read_text())
    out = {}
    for ctx in "ABC":
        bb = np.asarray(bz[f"beta_{ctx}"], float).copy()
        bb[~nonpanel] = 0.0
        out[ctx] = bb * np.sqrt(ref["base"]["AB"][ctx]["A"] / (dec["395"][ctx]["A"] * 18532 / 300))
    return out


def score(bulk, fold, truth, betas, sharp_q=None):
    """arms414.score's read-back (identical arithmetic) plus the board-emulated pds at lam 0.5 / 1 and per-target arrays."""
    scored = fold["scored"]
    nonpanel = np.ones(bulk.shape[1], dtype=bool)
    nonpanel[fold["panel_positions"]] = False
    out, X, Aa, pdsp, pdsl = {}, [], [], [], []
    AB = {}
    for ctx, (b, gate) in frames().items():
        F = bulk if sharp_q is None else sharpen(bulk, sharp_q, gate, mode="row")
        dl = deliver(F, b, gate, "minnorm")
        P = dl[:, scored]
        e = endpoints(P, truth)
        e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
        for lam in LAMS:
            e[f"pds_lam{lam:g}"] = endpoints((dl + lam * betas[ctx][None, :])[:, scored], truth)["pds"]
        out[ctx] = e
        X.append((P * truth).sum(1))
        Aa.append((P * P).sum(1))
        Tu = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-300)
        for Pm, acc in ((P, pdsp), ((dl + 0.5 * betas[ctx][None, :])[:, scored], pdsl)):
            Pu = Pm / np.maximum(np.linalg.norm(Pm, axis=1, keepdims=True), 1e-300)
            S = Pu @ Tu.T
            own = np.diag(S).copy()
            n = S.shape[0]
            acc.append(((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1))
        H = dl[:, nonpanel]
        a_ = float((H * H).sum(1).mean())
        Hc = H - H.mean(0, keepdims=True)
        b_ = float((Hc * Hc).sum(1).mean())
        AB[ctx] = dict(A=a_, B=b_, BoverA=b_ / a_)
    keys = ("pds", "cos", "sse", "cos_agg", "sgn_w", "pds_lam0.5", "pds_lam1")
    summ = {k: float(np.mean([out[c][k] for c in "ABC"])) for k in keys}
    summ["per_frame"] = {c: {k: out[c][k] for k in ("pds", "cos_agg", "pds_lam0.5")} for c in "ABC"}
    summ["AB"] = AB
    arrays = dict(x=np.asarray(X), a=np.asarray(Aa), pds=np.asarray(pdsp), pdsl=np.asarray(pdsl))
    return summ, arrays


def agg_parts(P, T):
    return (P * T).sum(1), (P * P).sum(1), float((T * T).sum())


def agg(w, pt, pp, tt):
    return float((w * pt).sum() / np.sqrt((w * w * pp).sum() * tt))
