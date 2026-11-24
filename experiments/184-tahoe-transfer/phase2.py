"""184 Phase 2 -- is a perturbation's response in a HELD-OUT cell line predictable beyond the
mean of the other lines, given only the held line's basal (DMSO) profile?

Leave-one-LINE-out over all 50 Tahoe lines.  For each held line L and each training-line set S:

  M            baseline: mean over S of the centred response to drug d.  Production's estimator,
               transplanted.
  gain-mean    M reweighted per gene by a gain that is the SAME for every recipient (no context).
               161's lesson: a flexible learner pointed at this objective finds a context-free
               per-entry reshaping, so that reshaping is the control the context arm must beat.
  gain-basal   M reweighted per gene by a gain that is a linear function of the recipient's basal
               log-expression of that gene.
  gain-ORACLE  the same diagonal family with the gain fitted on the held line's OWN truth.  Leaky
               by construction: it is the ceiling of everything diagonal, not an arm.
  bilinear     low-rank operator: in a rank-r response basis fitted on S, predict the held line's
               coefficients as W(u_L) c_bar(d), W(u_L) = W_0 + sum_k u_L[k] W_k, with u_L the
               recipient's basal embedding.  This is the only family that can ROTATE, and hence
               the only one that can move sign accuracy at all -- a positive diagonal gain leaves
               every sign untouched.
  bilinear-nc  the same operator with the context block removed (W_0 only).  Isolates the
               context terms.

Nulls, run on every context-carrying arm:
  perm-line    the recipient's basal is replaced by a random OTHER line's.
  perm-drug    the consensus fed to the operator is the consensus for a random OTHER drug.

Everything is fitted on training lines only: the response basis, the basal standardisation, the
basal embedding, the gains, the operator and the ridge parameter.

Writes work/phase2-*.json only.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

from common import SEED, WORK, load, unit

R2_EPS = 1e-30


# ----------------------------------------------------------------------------- endpoints
def top_index(truth, top=200):
    """Per-row indices of the row's own top-|truth| genes.  Depends only on the truth."""
    k = min(top, truth.shape[1])
    return np.argpartition(-np.abs(truth), k - 1, axis=1)[:, :k]


def sign_acc_idx(pred, truth, idx):
    st = np.take_along_axis(truth, idx, 1)
    sp = np.take_along_axis(pred, idx, 1)
    s = np.sign(st) * np.sign(sp)
    return ((s > 0).sum(1) + 0.5 * (s == 0).sum(1)) / idx.shape[1]


def cos_rows(A, B):
    return (unit(A) * unit(B)).sum(1)


def orth(A, B):
    """Component of each row of A orthogonal to the corresponding row of B."""
    u = unit(B)
    return A - (A * u).sum(1, keepdims=True) * u


def score(pred, truth, idx, base=None, weak=None):
    """sign / cos, plus two endpoints that a rescaling of the consensus cannot fake:

    cos_orth   161's endpoint -- cosine between the arm's and the truth's components ORTHOGONAL
               to the baseline consensus.  This is the line-specific residual, directly.
    sign_weak  174's framing -- sign accuracy restricted to the half of each row's significant
               genes where the consensus is WEAKEST.  That half is where a context model would
               have to earn its keep on the challenge; the strong half is already carried by the
               shared component.
    """
    out = {"sign": float(np.nanmean(sign_acc_idx(pred, truth, idx))),
           "cos": float(np.nanmean(cos_rows(pred, truth)))}
    if base is not None:
        out["cos_orth"] = float(np.nanmean(cos_rows(orth(pred, base), orth(truth, base))))
    if weak is not None:
        out["sign_weak"] = float(np.nanmean(sign_acc_idx(pred, truth, weak)))
    return out


# ----------------------------------------------------------------------------- gain family
def fit_gain(Rl, Mcons, lam_frac=1e-3):
    """Per-gene least-squares gain g[j] minimising sum_d (Rl[d,j] - g[j]*Mcons[d,j])^2."""
    num = (Rl * Mcons).sum(0)
    den = (Mcons * Mcons).sum(0)
    return num / (den + lam_frac * den.mean() + R2_EPS)


def gain_features(gbar, z):
    """[1, gbar, z, gbar*z] for one line: gbar (G,) context-free gain, z (G,) basal deviation."""
    return np.stack([np.ones_like(gbar), gbar, z, gbar * z], 1)


# ----------------------------------------------------------------------------- ridge helper
def ridge_solve(XtX, XtY, lam):
    n = XtX.shape[0]
    return np.linalg.solve(XtX + lam * np.eye(n) * np.trace(XtX) / n, XtY)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ks", default="3,5,10,20,40,49")
    ap.add_argument("--draws", type=int, default=4)
    ap.add_argument("--rank", type=int, default=32)
    ap.add_argument("--nctx", type=int, default=4)
    ap.add_argument("--top", type=int, default=200)
    ap.add_argument("--centring", default="plate")
    ap.add_argument("--lines", type=int, default=0, help="0 = all held lines")
    ap.add_argument("--compose", action="store_true", help="operator on the gain-corrected consensus")
    ap.add_argument("--out", default="work/phase2.json")
    a = ap.parse_args()
    ks = [int(x) for x in a.ks.split(",")]

    rng = np.random.default_rng(SEED)
    R, D, genes, lines, treats, basal, mask = load(a.centring)
    T, L, G = R.shape
    ok = np.isfinite(R[:, :, 0])
    R = np.nan_to_num(R, nan=0.0)
    logb = np.log10(basal + 1e-8)                       # (L, G)

    # per-line Gram, so any training set's response basis is a sum of precomputed blocks
    t0 = time.time()
    Gram = np.zeros((L, G, G), dtype=np.float32)
    for l in range(L):
        X = R[:, l, :][ok[:, l]]
        Gram[l] = X.T @ X
    print(f"per-line Grams in {time.time()-t0:.1f}s  ({Gram.nbytes/1e9:.2f} GB)")

    held_lines = list(range(L)) if not a.lines else list(range(a.lines))
    rows = []
    for held in held_lines:
        truth = R[:, held, :]
        idx = top_index(truth, a.top)
        good = ok[:, held]
        others = [l for l in range(L) if l != held]
        for k in ks:
            ndraw = 1 if k >= L - 1 else a.draws
            for draw in range(ndraw):
                S = np.asarray(others if k >= L - 1 else
                               rng.choice(others, size=k, replace=False))
                res = one_fold(R, ok, Gram, logb, S, held, truth, idx, good, a, rng)
                for arm, v in res.items():
                    rows.append({"held": int(held), "line": str(lines[held]), "k": int(k),
                                 "draw": int(draw), "arm": arm, **v})
        print(f"held {lines[held]:>10} ({held+1}/{len(held_lines)}) done  "
              f"{time.time()-t0:.0f}s")

    out = {"config": vars(a), "n_treatments": T, "n_lines": L, "n_genes": G, "rows": rows}
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), a.out), "w") as fh:
        json.dump(out, fh)
    print("wrote", a.out, len(rows), "rows")


def one_fold(R, ok, Gram, logb, S, held, truth, idx, good, a, rng):
    G = R.shape[2]
    k = len(S)
    res = {}

    # ---------- baseline: the equal-weight consensus of the training lines
    Mfull = R[:, S, :].mean(1)
    am = np.abs(np.take_along_axis(Mfull, idx, 1))
    half = idx.shape[1] // 2
    weak = np.take_along_axis(idx, np.argpartition(am, half - 1, axis=1)[:, :half], 1)
    res["M"] = score(Mfull, truth, idx, weak=weak)
    sc = lambda P: score(P, truth, idx, Mfull, weak)  # noqa: E731
    # rule 10: two rows per fold whose answer is known before the run
    res["KNOWN-truth"] = score(truth, truth, idx, weak=weak)
    res["KNOWN-zero"] = score(np.zeros_like(truth), truth, idx, weak=weak)

    # ---------- basal features, standardised on the TRAINING lines only
    mu = logb[S].mean(0)
    sd = logb[S].std(0) + 1e-6
    Z = (logb - mu) / sd                                    # (L, G)
    zh = Z[held]

    # ---------- diagonal gain family
    gh = []
    for i, l in enumerate(S):
        rest = np.delete(S, i)
        Mi = R[:, rest, :].mean(1)
        gh.append(fit_gain(R[:, l, :], Mi))
    gh = np.asarray(gh)                                     # (k, G) empirical per-line gains
    gbar = gh.mean(0)

    Xs, ys = [], []
    for i, l in enumerate(S):
        Xs.append(gain_features(gbar, Z[l]))
        ys.append(gh[i])
    X = np.concatenate(Xs)
    y = np.concatenate(ys)
    cf = np.linalg.lstsq(X[:, :2], y, rcond=None)[0]        # context-free: [1, gbar]
    cx = np.linalg.lstsq(X, y, rcond=None)[0]               # + context: [1, gbar, z, gbar*z]
    g_mean = gain_features(gbar, zh)[:, :2] @ cf
    g_basal = gain_features(gbar, zh) @ cx
    res["gain-mean"] = sc(Mfull * g_mean)
    res["gain-basal"] = sc(Mfull * g_basal)
    # a response can be attenuated or amplified by the recipient, not inverted: clipping the
    # gain at zero makes the family sign-preserving, which is the shippable version
    res["gain-basal-clip"] = sc(Mfull * np.maximum(g_basal, 0.0))
    # null: someone else's basal
    zn = Z[rng.choice([l for l in range(R.shape[1]) if l != held])]
    res["gain-basal/perm-line"] = sc(Mfull * (gain_features(gbar, zn) @ cx))
    # ORACLE ceiling of the diagonal family (leaky, not an arm)
    g_or = fit_gain(truth, Mfull)
    res["gain-ORACLE"] = sc(Mfull * g_or)

    # ---------- low-rank bilinear operator
    Gs = Gram[S].sum(0).astype(np.float64)
    w, V = np.linalg.eigh(Gs)
    V = V[:, ::-1][:, :a.rank]                              # (G, r), training-only basis
    r = V.shape[1]

    # basal embedding of the lines, PCA on the training lines' standardised log-basal
    Btr = (logb[S] - mu)
    Ub, Sb, Vb = np.linalg.svd(Btr - Btr.mean(0), full_matrices=False)
    Wb = Vb[:a.nctx].T                                      # (G, p)
    proj = (logb - mu - Btr.mean(0)) @ Wb                   # (L, p)
    pm, ps = proj[S].mean(0), proj[S].std(0) + 1e-9
    U = (proj - pm) / ps                                    # (L, p) standardised
    p = U.shape[1]

    # baseline for the operator: either the plain consensus, or the gain-corrected one
    Gtr = np.stack([gain_features(gbar, Z[l]) @ cx for l in S])      # (k, G) per-training-line
    Bl = Mfull * g_basal if a.compose else Mfull

    # Parameterised as a CORRECTION to the consensus: target = c(d, l') - cbar_{-l'}(d).
    # W = 0 reproduces the baseline exactly, so the ridge shrinks towards M and the operator can
    # only be credited for what it adds.
    Rs = R[:, S, :]
    CY = np.einsum("tlg,gr->tlr", Rs, V)                    # (T, k, r) each line's own response
    CB = np.empty_like(CY)                                  # (T, k, r) its leave-one-out input
    for i, l in enumerate(S):
        rest = [j for j in range(k) if j != i]
        base = Rs[:, rest, :].mean(1)
        if a.compose:
            base = base * Gtr[i]
        CB[:, i, :] = base @ V
    XtX = np.zeros((r * (1 + p), r * (1 + p)))
    XtY = np.zeros((r * (1 + p), r))
    blocks, feats = [], []
    for i, l in enumerate(S):
        u = np.concatenate([[1.0], U[l]])
        Xi = (CB[:, i, None, :] * u[None, :, None]).reshape(len(CB), -1)
        Yi = CY[:, i, :] - CB[:, i, :]
        m = good & ok[:, l]
        feats.append((Xi, Yi, m))
        blocks.append((Xi[m].T @ Xi[m], Xi[m].T @ Yi[m]))
        XtX += blocks[-1][0]
        XtY += blocks[-1][1]

    # ridge chosen by leave-one-TRAINING-line-out inside the training set (no held-line data)
    lams = [1e-4, 1e-3, 1e-2, 1e-1, 1.0]
    best, best_lam = -np.inf, lams[0]
    for lam in lams:
        sse = 0.0
        for i in range(k):
            Wm = ridge_solve(XtX - blocks[i][0], XtY - blocks[i][1], lam)
            Xi, Yi, m = feats[i]
            sse += float(((Xi[m] @ Wm - Yi[m]) ** 2).sum())
        if -sse > best:
            best, best_lam = -sse, lam
    Wop = ridge_solve(XtX, XtY, best_lam)
    res["ridge_lambda"] = {"sign": float(best_lam), "cos": float(best_lam)}

    cbar = Bl @ V
    other = rng.choice([l for l in range(R.shape[1]) if l != held])
    perm = rng.permutation(len(cbar))

    def bil(u_vec, cin):
        Xh = (cin[:, None, :] * np.concatenate([[1.0], u_vec])[None, :, None]).reshape(len(cin), -1)
        return Bl + (Xh @ Wop) @ V.T

    tag = "bilinear+gain" if a.compose else "bilinear"
    res[tag] = sc(bil(U[held], cbar))
    res[f"{tag}/perm-line"] = sc(bil(U[other], cbar))
    res[f"{tag}/perm-drug"] = sc(bil(U[held], cbar[perm]))

    # context-free operator: same fit with the context block dropped
    XtX0 = np.zeros((r, r))
    XtY0 = np.zeros((r, r))
    for i in range(k):
        m = feats[i][2]
        cb = CB[m, i, :]
        XtX0 += cb.T @ cb
        XtY0 += cb.T @ (CY[m, i, :] - cb)
    W0 = ridge_solve(XtX0, XtY0, best_lam)
    res[f"{tag}-nc"] = sc(Bl + (cbar @ W0) @ V.T)
    res["M/perm-drug"] = score(Mfull[perm], truth, idx, Mfull, weak)
    return res


if __name__ == "__main__":
    main()
