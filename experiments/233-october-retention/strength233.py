"""233 — THE SELECTION CONFOUND, measured: 199's synthetic panels are RANDOM draws; Arc's
panel is a SELECTED one.

AGENTS.md §2.1: the 300 targets were chosen from an H1 pilot screen, binned by DE-gene count,
sampled for response-mode diversity, and filtered for clean guide assignment and high basal
expression of the target -- and for "adequate representation in existing public datasets".
199's two synthetic bundles draw 300 genes at random from donor coverage.  They therefore
model the DISJOINTNESS and the COVERAGE of October's panel, but NOT its selection.  A random
300 contains perturbations whose true effect is near zero; those are unretrievable by
construction, and `pds` is a retrieval metric.

So a raw retention ratio confounds three things:

  (a) the decoy count            -- handled by `pds_at`, matched n
  (b) donor/residual coverage    -- handled by `thin233.py`
  (c) THE STRENGTH OF THE TARGETS THEMSELVES -- handled here

This measures (c) directly.  For each (panel, fold) it takes the held library's own effect
magnitude ||truth|| per target -- the honest, model-free measure of "how much did this
perturbation do" -- and then reports the members on STRENGTH-MATCHED, SIZE-MATCHED subsets:
common quantile bins over the pooled distribution, the same number of targets drawn from each
panel per bin, averaged over draws.  Anything left after that matching is the panel effect
that is NOT explained by the synthetic draws being unselected.

    uv run python strength233.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import panels233 as P  # noqa: E402
import score233 as S  # noqa: E402
import fold233  # noqa: E402
import plolo  # noqa: E402
from plolo import Operator, deliver  # noqa: E402

NBINS = 5
DRAWS = 64
SEED = 2331


def pds_cos_on(sub_P, sub_T):
    A = sub_P / np.maximum(np.linalg.norm(sub_P, axis=1, keepdims=True), 1e-300)
    B = sub_T / np.maximum(np.linalg.norm(sub_T, axis=1, keepdims=True), 1e-300)
    Smat = A @ B.T
    n = Smat.shape[0]
    own = np.diag(Smat).copy()
    beaten = (Smat < own[:, None]).sum(1)
    ties = (Smat == own[:, None]).sum(1) - 1
    pds = float(np.mean((beaten + 0.5 * ties) / (n - 1)))
    a = float((sub_P * sub_P).sum()); x = float((sub_P * sub_T).sum())
    d = float((sub_T * sub_T).sum())
    cos_agg = float(np.sqrt(max(0.0, x * x / (a * d)))) if a > 0 and d > 0 else 0.0
    return pds, cos_agg, float(np.mean(own))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(P.CANCER_FOLDS))
    ap.add_argument("--panels", nargs="*", default=["abc", "real", "pess"])
    args = ap.parse_args()
    resid218, compose218 = S._resid_module()

    bctx = {c: plolo.challenge_bctrl(c) for c in S.CONTEXTS}
    gate = {}
    for c in S.CONTEXTS:
        e = np.expm1(bctx[c])
        gate[c] = (1e6 * e / float(e.sum()) > plolo.MIN_CPM) & (bctx[c] > 0)

    store: dict[tuple, dict] = {}
    for panel in args.panels:
        S.bind_panel(panel, resid218)
        for held in args.folds:
            if not P.op_path(panel, held).exists():
                print(f"[{panel}/{held}] operator not fitted — skipped")
                continue
            fold = fold233.build(panel, held)
            sc = np.asarray(fold["scored"])
            wcols = np.asarray(fold["wcols"])
            Ts = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
            targets = fold["targets"]
            cd4 = S.cd4_prepared(panel, targets)
            h1 = compose218.union_rows(
                [resid218.source_residual(n, targets, "gene_norm") for n in S.UNION])
            op = Operator(fold["op_path"])
            pre, _ = fold233.consensus(fold, pm=True, only=None)
            final = S.build_final(fold, op, pre, cd4, h1, wcols, S.CD4_F, S.UNION_G,
                                  (resid218, compose218))
            dl = deliver(final, bctx["A"], gate["A"], "minnorm")[:, sc]
            strength = np.linalg.norm(Ts, axis=1)

            # ---- THE TRANSFERABILITY CEILING, model-free.
            # ||truth|| cannot separate a strong response from a noisy one.  The quantity
            # that actually bounds what any cross-line method can achieve on a target is how
            # much the SAME perturbation agrees between two OTHER lines: cos(source_i,
            # source_j) on that target's row.  It uses no prediction, so it is a property of
            # the PANEL, not of us.  AGENTS.md §2.1's selection criteria ("clean guide
            # assignment", "adequate representation in existing public datasets", binned by
            # DE-gene count) are exactly the criteria that would raise it.
            srcs = list(fold["sources"])
            E = {n: np.asarray(fold[f"eff__{n}"], dtype=np.float64)[:, sc] for n in srcs}
            A0 = E[srcs[0]] / np.maximum(np.linalg.norm(E[srcs[0]], axis=1, keepdims=True), 1e-300)
            A1 = E[srcs[1]] / np.maximum(np.linalg.norm(E[srcs[1]], axis=1, keepdims=True), 1e-300)
            both = (np.asarray(fold[f"avail__{srcs[0]}"], bool)
                    & np.asarray(fold[f"avail__{srcs[1]}"], bool))
            ceiling = np.einsum("ij,ij->i", A0, A1)
            ceiling[~both] = np.nan
            U = Ts / np.maximum(np.linalg.norm(Ts, axis=1, keepdims=True), 1e-300)
            truth_src = {n: np.einsum("ij,ij->i", U,
                                      E[n] / np.maximum(np.linalg.norm(E[n], axis=1,
                                                                       keepdims=True), 1e-300))
                         for n in srcs}
            for n in srcs:
                truth_src[n][~np.asarray(fold[f"avail__{n}"], bool)] = np.nan

            store[(panel, held)] = dict(P=dl, T=Ts, s=strength, ceiling=ceiling,
                                        truth_src=truth_src, srcs=srcs)
            print(f"[{panel}/{held}] transferability ceiling cos(src0, src1) median "
                  f"{np.nanmedian(ceiling):+.5f} mean {np.nanmean(ceiling):+.5f} "
                  f"(n={int(np.isfinite(ceiling).sum())}) | cos(truth, src) medians "
                  + " ".join(f"{n}={np.nanmedian(v):+.5f}" for n, v in truth_src.items()),
                  flush=True)
            q = np.quantile(strength, [0.1, 0.25, 0.5, 0.75, 0.9])
            # ---- PANEL GEOMETRY: how confusable is this decoy set, before any model?
            # AGENTS.md §2.1: today's 300 were sampled to MAXIMISE response-mode diversity.
            # A random draw is not, so its 300 truths sit closer together and retrieval is
            # harder BY CONSTRUCTION.  Same statistics AGENTS.md uses for the Orion blob.
            U = Ts / np.maximum(np.linalg.norm(Ts, axis=1, keepdims=True), 1e-300)
            G = U @ U.T
            off = ~np.eye(G.shape[0], dtype=bool)
            mean_pair = float(G[off].mean())
            mu = Ts.mean(0)
            blob = float((Ts @ mu) ** 2 @ np.ones(1) / 1) if False else float(
                ((Ts @ mu) ** 2).sum() / (float((mu * mu).sum()) * float((Ts * Ts).sum()))) \
                if float((mu * mu).sum()) > 0 else 0.0
            store[(panel, held)].update(mean_pairwise_cos_truth=mean_pair,
                                        truth_energy_in_mean=blob)
            print(f"[{panel}/{held}] n {len(strength)}  ||truth|| deciles "
                  f"p10 {q[0]:.3f} p25 {q[1]:.3f} MED {q[2]:.3f} p75 {q[3]:.3f} "
                  f"p90 {q[4]:.3f}  mean {strength.mean():.3f}", flush=True)
            # PREDICTION-side hubness (§4.2): if our own 300 predictions sit closer
            # together on one panel than the other, retrieval degrades even at equal cosine.
            V = dl / np.maximum(np.linalg.norm(dl, axis=1, keepdims=True), 1e-300)
            Gp = V @ V.T
            rr = float(Gp[off].mean())
            pmu = dl.mean(0)
            pblob = (float(((dl @ pmu) ** 2).sum()
                           / (float((pmu * pmu).sum()) * float((dl * dl).sum())))
                     if float((pmu * pmu).sum()) > 0 else 0.0)
            own = np.einsum("ij,ij->i", V, U)
            margin = float(np.median(own - np.sort(Gp * 0 + (V @ U.T), axis=1)[:, -2]))
            # per-target retrieval, to test the ceiling as the CAUSE within a panel
            Sm = V @ U.T
            n_ = Sm.shape[0]
            ownd = np.diag(Sm).copy()
            beaten = (Sm < ownd[:, None]).sum(1)
            ties = (Sm == ownd[:, None]).sum(1) - 1
            pds_t = (beaten + 0.5 * ties) / (n_ - 1)
            ok = np.isfinite(ceiling)
            r_ceiling = float(np.corrcoef(ceiling[ok], pds_t[ok])[0, 1])
            r_strength = float(np.corrcoef(strength, pds_t)[0, 1])
            store[(panel, held)].update(rr_pred=rr, pred_energy_in_mean=pblob,
                                        own_cos_median=float(np.median(own)),
                                        margin_median=margin, pds_t=pds_t,
                                        r_pds_ceiling=r_ceiling, r_pds_strength=r_strength)
            print(f"{'':17}decoy-set geometry: mean pairwise cos(truth) {mean_pair:+.5f}"
                  f"   truth energy in the mean vector {100 * blob:.2f}%", flush=True)
            print(f"{'':17}prediction hubness: mean pairwise cos(pred) {rr:+.5f}"
                  f"   pred energy in the mean vector {100 * pblob:.2f}%"
                  f"   own-cos median {np.median(own):+.5f}"
                  f"   retrieval margin (own - best decoy) {margin:+.5f}", flush=True)
            print(f"{'':17}within-panel r(per-target pds, ceiling) = {r_ceiling:+.4f}"
                  f"   r(per-target pds, ||truth||) = {r_strength:+.4f}", flush=True)
            del op, final

    out: dict[str, object] = {}
    print("\nSTRENGTH-MATCHED, SIZE-MATCHED retention (common bins over the pooled "
          "||truth||, equal draws per bin)\n")
    print(f"{'fold':15}{'panel':7}{'n_match':>9}{'pds':>10}{'cos_agg':>10}"
          f"{'pds credit':>12}{'retention':>11}")
    rng = np.random.default_rng(SEED)
    for held in args.folds:
        avail = [p for p in args.panels if (p, held) in store]
        if "abc" not in avail or len(avail) < 2:
            continue
        pooled = np.concatenate([store[(p, held)]["s"] for p in avail])
        edges = np.quantile(pooled, np.linspace(0, 1, NBINS + 1))
        edges[0], edges[-1] = -np.inf, np.inf
        bins = {p: [np.flatnonzero((store[(p, held)]["s"] >= edges[b])
                                   & (store[(p, held)]["s"] < edges[b + 1]))
                    for b in range(NBINS)] for p in avail}
        take = [min(len(bins[p][b]) for p in avail) for b in range(NBINS)]
        n_match = int(sum(take))
        base = None
        for p in avail:
            vals = []
            for _ in range(DRAWS):
                idx = np.concatenate([rng.choice(bins[p][b], size=take[b], replace=False)
                                      for b in range(NBINS) if take[b] > 0])
                vals.append(pds_cos_on(store[(p, held)]["P"][idx],
                                       store[(p, held)]["T"][idx]))
            pds = float(np.mean([v[0] for v in vals]))
            cos = float(np.mean([v[1] for v in vals]))
            cred = pds - 0.5
            if p == "abc":
                base = cred
            ret = 100 * cred / base if base else float("nan")
            out[f"{held}|{p}"] = dict(n_match=n_match, pds=pds, cos_agg=cos,
                                      pds_credit=cred, retention_pct=ret,
                                      bin_counts=[int(t) for t in take])
            print(f"{held:15}{p:7}{n_match:9d}{pds:10.5f}{cos:10.5f}{cred:12.5f}"
                  f"{ret:10.1f}%")
        print()
    print("\nCEILING-MATCHED, SIZE-MATCHED retention (bins over the pooled "
          "cos(src0, src1) transferability ceiling)\n")
    print(f"{'fold':15}{'panel':7}{'n_match':>9}{'pds':>10}{'cos_agg':>10}"
          f"{'pds credit':>12}{'retention':>11}")
    for held in args.folds:
        avail = [p for p in args.panels if (p, held) in store]
        if "abc" not in avail or len(avail) < 2:
            continue
        pooled = np.concatenate([store[(p, held)]["ceiling"] for p in avail])
        pooled = pooled[np.isfinite(pooled)]
        edges = np.quantile(pooled, np.linspace(0, 1, NBINS + 1))
        edges[0], edges[-1] = -np.inf, np.inf
        bins = {p: [np.flatnonzero(np.isfinite(store[(p, held)]["ceiling"])
                                   & (store[(p, held)]["ceiling"] >= edges[b])
                                   & (store[(p, held)]["ceiling"] < edges[b + 1]))
                    for b in range(NBINS)] for p in avail}
        take = [min(len(bins[p][b]) for p in avail) for b in range(NBINS)]
        n_match = int(sum(take))
        base = None
        for p in avail:
            vals = []
            for _ in range(DRAWS):
                idx = np.concatenate([rng.choice(bins[p][b], size=take[b], replace=False)
                                      for b in range(NBINS) if take[b] > 0])
                vals.append(pds_cos_on(store[(p, held)]["P"][idx],
                                       store[(p, held)]["T"][idx]))
            pds = float(np.mean([v[0] for v in vals]))
            cos = float(np.mean([v[1] for v in vals]))
            cred = pds - 0.5
            if p == "abc":
                base = cred
            out[f"ceilmatch|{held}|{p}"] = dict(n_match=n_match, pds=pds, cos_agg=cos,
                                                pds_credit=cred,
                                                retention_pct=100 * cred / base if base else None,
                                                bin_counts=[int(t) for t in take])
            print(f"{held:15}{p:7}{n_match:9d}{pds:10.5f}{cos:10.5f}{cred:12.5f}"
                  f"{(100 * cred / base if base else float('nan')):10.1f}%")
        print()

    (P.WORK / "strength233.json").write_text(json.dumps(
        {k: v for k, v in out.items()}
        | {f"dist|{p}|{h}": dict(
            n=int(len(store[(p, h)]["s"])),
            mean_pairwise_cos_truth=store[(p, h)]["mean_pairwise_cos_truth"],
            truth_energy_in_mean=store[(p, h)]["truth_energy_in_mean"],
            ceiling_median=float(np.nanmedian(store[(p, h)]["ceiling"])),
            ceiling_mean=float(np.nanmean(store[(p, h)]["ceiling"])),
            truth_src_median={n: float(np.nanmedian(v))
                              for n, v in store[(p, h)]["truth_src"].items()},
            rr_pred=store[(p, h)]["rr_pred"],
            pred_energy_in_mean=store[(p, h)]["pred_energy_in_mean"],
            own_cos_median=store[(p, h)]["own_cos_median"],
            margin_median=store[(p, h)]["margin_median"],
            r_pds_ceiling=store[(p, h)]["r_pds_ceiling"],
            r_pds_strength=store[(p, h)]["r_pds_strength"],
            quantiles=[float(x) for x in np.quantile(store[(p, h)]["s"],
                                                     [0, .1, .25, .5, .75, .9, 1])],
            mean=float(store[(p, h)]["s"].mean()))
           for (p, h) in store}, indent=2))
    print("-> work/strength233.json")


if __name__ == "__main__":
    main()
