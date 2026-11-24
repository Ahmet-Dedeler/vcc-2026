"""224 Phase 1 — the ladder, cheapest first, leave-one-LINE-out.

Endpoint: **transfer cosine on the held LINE, macro over targets**, against the equal-weight
donor mean of the other lines.  `pds` (300-decoy) is a free secondary and is NOT the gate.

Every rung is fitted by inner leave-one-training-LINE-out inside the donor pool `D`; the held
line `A` contributes no row to any fit at any stage.  Targets are split once per surface into
`T_fit` / `T_eval`; hyperparameters are selected inside `T_fit` on a line-AND-target holdout
(see `core224.select_operator` for why the target half of that matters).  Both the held-LINE
and the held-PERTURBATION endpoints are read on `T_eval`.

    uv run python ladder224.py --surface S7
"""
from __future__ import annotations

import argparse
import itertools
import json
import time

import numpy as np

from core224 import (DIAG_SHRINK, WORK, basal_bins, decoy_cos, diag_fit,  # noqa: E402
                     donor_mean, inner_pairs, load_surface, macro_cos,
                     per_target_cos, pds, ridge_svd, select_operator, unit)
from features224 import pca, target_features  # noqa: E402

SEED = 20260904
MAX_SUBSETS = 4
STACK_B = (0.0, 0.02, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6)
BASAL_BINS = 12


def run_config(C, basal, D, A, tsel, tval, tfit, teval, Xfeat, rng, perm=None):
    """One (held line A, donor subset D) configuration.  `perm` deranges the fit labels'
    target correspondence (Phase-2 falsifier 2) and is applied to fit rows only."""
    G = C.shape[2]
    ps = None if perm is None else np.argsort(np.argsort(perm[: len(tsel)]))
    pv = None if perm is None else np.argsort(np.argsort(perm[: len(tval)]))
    pf = None if perm is None else np.argsort(np.argsort(perm[: len(tfit)]))
    Xp, Yv, Lv = inner_pairs(C, D, tfit, pf)

    Pheld = donor_mean(C, D, teval)
    Theld = C[A][teval].astype(np.float64)
    Pp, Yp = [], []
    for B in D:
        Pp.append(donor_mean(C, [x for x in D if x != B], teval))
        Yp.append(C[B][teval].astype(np.float64))
    Ppert, Ypert = np.concatenate(Pp), np.concatenate(Yp)

    out, pt = {}, {}

    def record(tag, Pl, Pp_, **kw):
        v = per_target_cos(Pl, Theld)
        dec = decoy_cos(Pl, Theld)
        pt[tag] = v.astype(np.float32)
        out[tag] = dict(line=float(v.mean()), pert=macro_cos(Pp_, Ypert),
                        pds=pds(Pl, Theld), decoy=dec, excess=float(v.mean()) - dec,
                        cos_to_R0=macro_cos(Pl, Pheld), **kw)

    record("R0", Pheld, Ppert)

    # ---- R1: ridge on target features -> response, REPLACING the donor mean (132's family)
    if Xfeat is not None:
        Z = pca(Xfeat, np.asarray(tfit), 50)
        Zi = np.tile(Z[tfit], (len(D), 1))
        Gm = Zi.T @ Zi + 10.0 * np.eye(Zi.shape[1])
        W = np.linalg.solve(Gm, Zi.T @ Yv)
        Fl, Fpv = Z[teval] @ W, np.tile(Z[teval], (len(D), 1)) @ W
        record("R1_feat_replace", Fl, Fpv)
        Ffit = Zi @ W
        b = max(STACK_B, key=lambda b_: macro_cos(Xp + b_ * unit(Ffit), Yv))
        record("R2_feat_stack", Pheld + b * unit(Fl), Ppert + b * unit(Fpv), b=float(b))
        Wr = np.linalg.solve(Gm, Zi.T @ (Yv - Xp))
        Rfit = Zi @ Wr
        br = max(STACK_B, key=lambda b_: macro_cos(Xp + b_ * Rfit, Yv))
        record("R2_feat_resid", Pheld + br * (Z[teval] @ Wr),
               Ppert + br * (np.tile(Z[teval], (len(D), 1)) @ Wr), b=float(br))

    # ---- R3: learned gene-space linear operator P -> P M
    if len(D) >= 3:
        hold = int(rng.integers(len(D)))
        (rank, lam), inner, _ = select_operator(C, D, tsel, tval, hold, ps, pv)
        apply, sv = ridge_svd(Xp, Yv)
        record("R3_operator", apply(Pheld, rank, lam), apply(Ppert, rank, lam),
               rank=int(rank), lam=float(lam), inner_cos=inner, n_sv=int(len(sv)),
               hold=hold)
        # a shrunk blend with the donor mean, weight chosen on the same honest holdout
        Pv_, Yv_, Fv_ = inner_pairs(C, D, tval, pv)
        Ps_, Ys_, Fs_ = inner_pairs(C, D, tsel, ps)
        ap2, _ = ridge_svd(Ps_[Fs_ != hold], Ys_[Fs_ != hold])
        mv = Fv_ == hold
        cand = ap2(Pv_[mv], rank, lam)
        w = max((0.0, 0.25, 0.5, 0.75, 1.0),
                key=lambda w_: macro_cos(w_ * unit(cand) + (1 - w_) * unit(Pv_[mv]), Yv_[mv]))
        record("R3_blend", w * unit(apply(Pheld, rank, lam)) + (1 - w) * Pheld,
               w * unit(apply(Ppert, rank, lam)) + (1 - w) * Ppert, w=float(w))

    # ---- R4n: recipient-INDEPENDENT per-gene diagonal (the null rung 4 must beat)
    gene_bins = np.arange(G)
    Ps_, Ys_, Fs_ = inner_pairs(C, D, tsel, ps)
    Pv_, Yv_, Fv_ = inner_pairs(C, D, tval, pv)
    hold4 = int(rng.integers(len(D)))
    ms, mv = Fs_ != hold4, Fv_ == hold4
    sh = max(DIAG_SHRINK, key=lambda s: macro_cos(
        Pv_[mv] * diag_fit(Ps_[ms], Ys_[ms], gene_bins, G, s), Yv_[mv]))
    g = diag_fit(Xp, Yv, gene_bins, G, sh)
    record("R4_diag_pergene", Pheld * g, Ppert * g, shrink=float(sh),
           gain_p10=float(np.quantile(g, 0.1)), gain_p90=float(np.quantile(g, 0.9)),
           gain_corr_R0norm=float(np.corrcoef(
               g, np.abs(Pheld).mean(0))[0, 1]))
    # MOMENT-MATCHING control: the per-gene gain that makes the prediction's per-gene sd
    # equal the labels', computed WITHOUT ever using the (target -> target) correspondence.
    # If this matches R4_diag_pergene then the "learned" gain is variance matching and no
    # cross-line pairing was needed to find it; if it does not, the pairing is load-bearing.
    pt["_gain"] = g.astype(np.float32)
    gmm = Yv.std(0) / np.maximum(Xp.std(0), 1e-300)
    record("R4_diag_moment_match", Pheld * gmm, Ppert * gmm,
           corr_to_fitted=float(np.corrcoef(g, gmm)[0, 1]))
    # 004's own null: the SAME gain vector with its gene assignment permuted.  A gain that
    # still helps here is not reading gene identity and is therefore not a transfer model.
    record("R4_diag_geneshuffle", Pheld * g[rng.permutation(G)],
           Ppert * g[rng.permutation(G)], shrink=float(sh))
    # 004's estimator, verbatim: fitted on SINGLE-donor -> line pairs, applied to the
    # CONSENSUS.  If this loses where the consensus-fitted one wins, the mechanism that
    # explains 004's negative is "the weights were fitted for a different predictor".
    Ps1, Ys1 = [], []
    for B in D:
        for x in D:
            if x != B:
                Ps1.append(unit(C[x][tfit].astype(np.float64)))
                Ys1.append(unit(C[B][tfit].astype(np.float64)))
    g1 = diag_fit(np.concatenate(Ps1), np.concatenate(Ys1), gene_bins, G, sh)
    record("R4_diag_from_single", Pheld * g1, Ppert * g1, shrink=float(sh),
           cos_g_to_gconsensus=float(np.corrcoef(g, g1)[0, 1]))
    # ORACLE: the diagonal fitted on the HELD line itself, on T_fit.  Leaky by construction
    # and labelled as such — it is the ceiling of the diagonal family, not a candidate.
    go = diag_fit(donor_mean(C, D, tfit), unit(C[A][tfit].astype(np.float64)),
                  gene_bins, G, sh)
    record("R4_diag_ORACLE_leaky", Pheld * go, Ppert, shrink=float(sh),
           cos_g_to_oracle=float(np.corrcoef(g, go)[0, 1]))

    # ---- R4a: recipient-BASAL-conditioned per-gene gain (184's shape, on CRISPRi)
    nb = BASAL_BINS
    num = np.zeros(nb)
    den = np.zeros(nb)
    for j, B in enumerate(D):
        rest = [x for x in D if x != B]
        bj = basal_bins(basal[B].astype(np.float64) - basal[rest].mean(0).astype(np.float64), nb)
        m = Lv == j
        num += np.bincount(bj, weights=(Xp[m] * Yv[m]).sum(0), minlength=nb)
        den += np.bincount(bj, weights=(Xp[m] * Xp[m]).sum(0), minlength=nb)
    gv = num / np.maximum(den, 1e-300)
    bA = basal_bins(basal[A].astype(np.float64) - basal[list(D)].mean(0).astype(np.float64), nb)
    gp = np.concatenate([np.broadcast_to(
        gv[basal_bins(basal[B].astype(np.float64)
                      - basal[[x for x in D if x != B]].mean(0).astype(np.float64), nb)],
        (len(teval), G)) for B in D])
    record("R4_basal_gain", Pheld * gv[bA], Ppert * gp, bins=nb,
           gain_spread=float(gv.max() / max(gv.min(), 1e-9)))

    # ---- R4b: basal-similarity donor weights (010's family)
    def wts(ref, pool, beta):
        s = np.asarray([float(np.corrcoef(basal[ref].astype(np.float64),
                                          basal[x].astype(np.float64))[0, 1]) for x in pool])
        w = np.exp(beta * (s - s.max()))
        return w / w.sum()
    bb, bv = 0.0, -np.inf
    for beta in (0.0, 1.0, 2.0, 4.0, 8.0, 16.0):
        tot = 0.0
        for B in D:
            rest = [x for x in D if x != B]
            w = wts(B, rest, beta)
            tot += macro_cos(np.sum([w[i] * unit(C[x][tfit].astype(np.float64))
                                     for i, x in enumerate(rest)], axis=0),
                             C[B][tfit].astype(np.float64))
        if tot > bv:
            bv, bb = tot, beta
    w = wts(A, list(D), bb)
    Pw = np.sum([w[i] * unit(C[x][teval].astype(np.float64)) for i, x in enumerate(D)], axis=0)
    Pwp = []
    for B in D:
        rest = [x for x in D if x != B]
        wj = wts(B, rest, bb)
        Pwp.append(np.sum([wj[i] * unit(C[x][teval].astype(np.float64))
                           for i, x in enumerate(rest)], axis=0))
    record("R4_basal_weight", Pw, np.concatenate(Pwp), beta=float(bb))
    return out, pt


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--surface", default="S7")
    ap.add_argument("--tag", default=None)
    ap.add_argument("--kmin", type=int, default=2)
    ap.add_argument("--permute-labels", action="store_true",
                    help="Phase-2 falsifier 2: derange the fit labels' target correspondence")
    ap.add_argument("--max-subsets", type=int, default=MAX_SUBSETS)
    args = ap.parse_args()

    rng = np.random.default_rng(SEED)
    C, basal, members, groups, targets, genes, dropped = load_surface(args.surface)
    L, T = len(members), C.shape[1]
    print(f"surface {args.surface}: {L} LINES {members}  ({T} targets x {C.shape[2]} genes)"
          f"   dropped-as-same-line: {dropped}")

    perm = rng.permutation(T)
    tfit = np.sort(perm[: T // 2])
    teval = np.sort(perm[T // 2:])
    q = rng.permutation(len(tfit))
    tsel = np.sort(tfit[q[: int(0.6 * len(tfit))]])
    tval = np.sort(tfit[q[int(0.6 * len(tfit)):]])
    Xfeat, have, fmeta = target_features(targets)
    lab_perm = rng.permutation(max(len(tfit), len(tsel), len(tval))) \
        if args.permute_labels else None
    print(f"features: {fmeta}\nT_fit {len(tfit)} (sel {len(tsel)} / val {len(tval)})  "
          f"T_eval {len(teval)}   permute_labels={args.permute_labels}")

    rows, ptall, t0 = [], {}, time.time()
    for A in range(L):
        others = [x for x in range(L) if x != A]
        for k in range(args.kmin, len(others) + 1):
            combos = list(itertools.combinations(others, k))
            if len(combos) > args.max_subsets:
                combos = [combos[i] for i in
                          rng.choice(len(combos), args.max_subsets, replace=False)]
            for D in combos:
                res, pt = run_config(C, basal, list(D), A, tsel, tval, tfit, teval, Xfeat,
                                     rng, lab_perm)
                cid = len(rows)
                rows.append(dict(cid=cid, held=members[A], k=k,
                                 donors=[members[x] for x in D], **res))
                for tag, v in pt.items():
                    ptall[f"{cid}|{tag}"] = v
                print(f"  [{time.time() - t0:6.0f}s] {members[A]:14s} k={k} "
                      f"R0 {res['R0']['line']:.5f} | " +
                      " ".join(f"{m.split('_')[0]}:{res[m]['line']:+.5f}"
                               for m in res if m != "R0"), flush=True)

    out = dict(surface=args.surface, members=members, groups=groups, seed=SEED,
               n_targets=int(T), n_genes=int(C.shape[2]), t_fit=len(tfit),
               t_sel=len(tsel), t_val=len(tval), t_eval=len(teval), features=fmeta,
               permute_labels=bool(args.permute_labels), dropped_same_line=dropped, rows=rows)
    p = WORK / f"ladder224-{args.tag or args.surface}.json"
    p.write_text(json.dumps(out, indent=1))
    np.savez_compressed(WORK / f"ptcos-{args.tag or args.surface}.npz",
                        teval=teval, **ptall)
    print(f"\nwrote {p}   ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
