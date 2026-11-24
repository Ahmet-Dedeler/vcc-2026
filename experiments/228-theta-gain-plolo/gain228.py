"""228 Phase 1 — fit 224's per-gene diagonal gain ON THE OFFICIAL GENE AXIS, per PLOLO fold.

WHAT IS DIFFERENT FROM 224
--------------------------
224 fitted the gain on a surface of its own making (an intersection gene axis of 5,070-7,312
columns, its own mean-removal, its own re-fitted whitener).  This script fits the SAME
estimator on the objects PRODUCTION actually builds:

  * every fit line is mapped onto the official 18,533-gene axis by `E062.prepare_main_source`
    -- production's own loader, including its per-line `EFFECT_SCALE / raw_median` rescale,
    its own mean-removal, and its zeroing of the 300 panel coordinates;
  * the consensus is `E062.build_baseline` -- production's own vote-weighted mean;
  * the whitening is the FOLD'S OWN operator at production's `gamma = 0.65`, `lam = 0`,
    applied through `plolo.Operator` -- not a re-fit of the family.

So the gain that comes out is directly multiplicable against the production effect matrix, and
no symbol remapping happens between fit and deploy.

LEAKAGE
-------
The fit pool excludes the PLOLO held line by cell-line GROUP.  Fit targets are the NON-PANEL
intersection of the fit pool's perturbations; the PLOLO evaluation targets are panel genes.
`fit and eval` is asserted empty.  `h1_2025` is never in the fit pool.

FOUR NULLS, all produced by the SAME code path with one input changed:
  shuf      g permuted over genes                                  (224: -7.2%)
  mm        sd(Y)/sd(P) per gene, no target correspondence         (224: -3.5%)
  derange   g fitted with the target labels deranged               (224: -25.7%)
  basal12   224's R4a: a 12-bin gain conditioned on the RECIPIENT's basal profile, fitted
            with the RIGHT recipient and with a WRONG one          (224: +0.41%, wrong>right)

    uv run python gain228.py [--folds ...]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, OFFICIAL, PANEL, PANEL_SET, PROD_OPS,  # noqa: E402
                   Operator, build_fold)

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
SEED = 20260904
DIAG_SHRINK = (0.0, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)
CHUNK = 1024

# name -> (path, cell-line GROUP).  Groups are what a LOLO split may not straddle.
POOL: dict[str, tuple[str, str]] = {
    "k562_gwps":     ("data/donors/k562_gwps", "K562"),
    "orion_hct116":  ("data/donors/orion_hct116", "HCT116"),
    "orion_hek293t": ("data/donors/orion_hek293t", "HEK293T"),
    "rpe1":          ("data/donors/rpe1", "RPE1"),
    "nadig_jurkat":  ("data/donors/nadig_jurkat", "JURKAT"),
    "nadig_hepg2":   ("experiments/166-donor-acquisition/donors/nadig_hepg2", "HEPG2"),
}
HELD_GROUP = {"k562_gwps": "K562", "orion_hct116": "HCT116",
              "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}


def unit32(A: np.ndarray) -> np.ndarray:
    n = np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)
    return (A / n).astype(np.float32)


def accum(P, Y, rows=None, yperm=None):
    """float64 (num, den) accumulation over float32 storage.

    `yperm` supplies the row indices of Y (the target-label derangement)."""
    idx = np.arange(P.shape[0]) if rows is None else np.flatnonzero(rows)
    num = np.zeros(P.shape[1]); den = np.zeros(P.shape[1])
    for s in range(0, len(idx), CHUNK):
        j = idx[s:s + CHUNK]
        p = P[j].astype(np.float64)
        y = Y[j if yperm is None else yperm[j]].astype(np.float64)
        num += (p * y).sum(0); den += (p * p).sum(0)
    return num, den


def accum_masked(P, Y, F, cov, n_folds, yperm=None):
    """(num, den) with each inner fold's contribution restricted to the genes that fold's
    LABEL LINE actually observes.

    `Y_ij == 0` for a gene the inner-held line does not carry means UNOBSERVED, not
    "the response is zero".  Summing `P_ij^2` over those rows biases `g_j` down by the
    gene's coverage fraction, which is a coverage weight masquerading as a reliability
    weight.  224 never met this because its gene axis was the INTERSECTION of its lines.
    """
    num = np.zeros(P.shape[1]); den = np.zeros(P.shape[1])
    for j in range(n_folds):
        n_j, d_j = accum(P, Y, F == j, yperm)
        num += n_j * cov[j]; den += d_j * cov[j]
    return num, den


def shrunk(num, den, shrink):
    """224 `core224.diag_fit`: g = <P,Y> / <P,P> shrunk toward 1."""
    if shrink <= 0:
        return num / np.maximum(den, 1e-300)
    pri = np.median(den[den > 0]) if (den > 0).any() else 1.0
    return (num + shrink * pri) / np.maximum(den + shrink * pri, 1e-300)


def shrunk_pooled(num, den, shrink):
    """The same family, shrunk toward the POOLED coefficient instead of toward 1.

    `shrunk`'s prior target is 1, but the natural scale of `g` here is the transfer cosine
    (~0.02 in the whitened challenge-panel frame, ~0.2 on 224's raw surface), so a prior at
    1 is 50x away from the data and the grid collapses the per-gene spread almost at once.
    Shrinking toward `sum(num)/sum(den)` is the same interpolation with a prior on scale.
    """
    live = den > 0
    c0 = float(num[live].sum() / max(den[live].sum(), 1e-300)) if live.any() else 1.0
    if shrink <= 0:
        return num / np.maximum(den, 1e-300)
    pri = np.median(den[live]) if live.any() else 1.0
    return (num + shrink * pri * c0) / np.maximum(den + shrink * pri, 1e-300)


def neutral(g, live):
    """Genes the fit pool cannot adjudicate get the covered median, i.e. no reweighting."""
    out = np.asarray(g, dtype=np.float64).copy()
    med = float(np.median(out[live])) if live.any() else 1.0
    out[~live] = med
    return out


def binned(num, den, bins, n_bins, shrink):
    return shrunk(np.bincount(bins, weights=num, minlength=n_bins),
                  np.bincount(bins, weights=den, minlength=n_bins), shrink)[bins]


def macro_cos_g(P, Y, g, rows=None):
    """mean_i cos(P_i * g, Y_i); P, Y already row-unit."""
    idx = np.arange(P.shape[0]) if rows is None else np.flatnonzero(rows)
    tot = 0.0
    for s in range(0, len(idx), CHUNK):
        j = idx[s:s + CHUNK]
        p = P[j].astype(np.float64) * g
        n = np.maximum(np.linalg.norm(p, axis=1), 1e-300)
        tot += float(((p / n[:, None]) * Y[j].astype(np.float64)).sum())
    return tot / max(len(idx), 1)


def basal_bins(b: np.ndarray, n_bins: int) -> np.ndarray:
    q = np.quantile(b, np.linspace(0, 1, n_bins + 1)[1:-1])
    return np.clip(np.searchsorted(q, b), 0, n_bins - 1).astype(np.int64)


def read_perts(rel: str) -> set[str]:
    return {x.strip() for x in (ROOT / rel / "perts.txt").read_text().splitlines() if x.strip()}


def fit_one(held: str, args) -> dict:
    t0 = time.time()
    rng = np.random.default_rng(args.seed)
    fold = build_fold(held)
    wcols = np.asarray(fold["wcols"])
    op_path = (PROD_OPS / "whitener_union_pruned.npz"
               if FOLDS[held]["op_kind"] == "production" else Path(FOLDS[held]["op"]))

    fit_lines = [n for n, (_, g) in POOL.items() if g != HELD_GROUP[held]]
    fit_targets = sorted(set.intersection(*[read_perts(POOL[n][0]) for n in fit_lines])
                         - PANEL_SET)
    eval_targets = set(fold["targets"])
    assert not (set(fit_targets) & eval_targets), "KA-4: fit/eval target overlap"
    assert not (set(fit_targets) & PANEL_SET), "a panel target reached the fit pool"
    print(f"[{held}] fit lines {fit_lines}\n"
          f"[{held}] fit targets {len(fit_targets)} NON-panel; eval {len(eval_targets)} "
          f"panel; overlap 0 (asserted)", flush=True)

    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos],
                                 dtype=np.int64)

    G = len(OFFICIAL)
    prepared, covered = {}, {}
    for n in fit_lines:
        bank = E062.Bank.load(n, ROOT / POOL[n][0])
        amap = E062.map_bank_to_axis(bank, OFFICIAL)
        c = np.zeros(G, dtype=bool)
        c[amap.unique_target_positions] = True
        c[panel_positions] = False
        covered[n] = c
        prepared[n] = E062.prepare_main_source(bank, fit_targets, OFFICIAL, panel_positions,
                                               chunk_rows=256, max_rss_gib=8.0)
        assert prepared[n].target_available.all(), f"{n} missing a fit target"
        print(f"  [{held}] prepared {n:15s} {prepared[n].mapped_gene_count} mapped genes, "
              f"{int(c.sum())} non-panel covered ({time.time() - t0:.0f}s)", flush=True)
        del bank

    T = len(fit_targets)
    op = Operator(op_path)
    P = np.empty((len(fit_lines) * T, G), dtype=np.float32)
    Y = np.empty_like(P)
    F = np.empty(len(fit_lines) * T, dtype=np.int32)
    for j, B in enumerate(fit_lines):
        base, votes = E062.build_baseline([prepared[n] for n in fit_lines if n != B])
        assert (votes > 0).all(), "a fit target has no donor vote"
        Pw = op.apply(base * ALPHA, wcols, GAMMA, LAM)
        P[j * T:(j + 1) * T] = unit32(Pw)
        Y[j * T:(j + 1) * T] = unit32(np.asarray(prepared[B].effects, dtype=np.float64))
        F[j * T:(j + 1) * T] = j
        del base, Pw
        print(f"  [{held}] inner fold {B:15s} done ({time.time() - t0:.0f}s)", flush=True)
    del op, prepared
    cov = [covered[n].astype(np.float64) for n in fit_lines]
    n_fold = len(fit_lines)
    cov_count = np.sum(cov, axis=0)
    full_cov = cov_count >= n_fold - 1e-9          # the 224-faithful intersection subset

    # ---- shrinkage selected on a held inner LINE **and** a disjoint target block (224 3.2)
    q = rng.permutation(T)
    tsel = np.zeros(T, dtype=bool)
    tsel[q[: int(0.6 * T)]] = True
    hold = int(rng.integers(n_fold))
    rowsel = np.tile(tsel, n_fold) & (F != hold)
    rowval = np.tile(~tsel, n_fold) & (F == hold)
    nsel = np.zeros(G); dsel = np.zeros(G)
    for j in range(n_fold):
        if j == hold:
            continue
        n_j, d_j = accum(P, Y, rowsel & (F == j))
        nsel += n_j * cov[j]; dsel += d_j * cov[j]
    scores = {s: macro_cos_g(P, Y, neutral(shrunk(nsel, dsel, s), dsel > 0), rowval)
              for s in DIAG_SHRINK}
    shrink = max(scores, key=scores.get)
    scores_p = {s: macro_cos_g(P, Y, neutral(shrunk_pooled(nsel, dsel, s), dsel > 0), rowval)
                for s in DIAG_SHRINK}
    shrink_p = max(scores_p, key=scores_p.get)
    print(f"  [{held}] shrink {shrink:g} (toward 1) / {shrink_p:g} (toward pooled) selected "
          f"on held line {fit_lines[hold]} + disjoint target block\n"
          f"        toward-1     {' '.join(f'{k:g}:{v:.5f}' for k, v in scores.items())}\n"
          f"        toward-pool  {' '.join(f'{k:g}:{v:.5f}' for k, v in scores_p.items())}",
          flush=True)

    # ---- the arm and the three gene-space nulls, same code path, one input changed
    num, den = accum_masked(P, Y, F, cov, n_fold)
    live = den > 0
    g = neutral(shrunk(num, den, shrink), live)
    g_reg = neutral(shrunk_pooled(num, den, shrink_p), live)
    g_shuf = g[rng.permutation(G)]

    # moment-matching control: the same estimator with the TARGET CORRESPONDENCE removed
    N = P.shape[0]
    m1 = np.zeros(G); m2 = np.zeros(G); y1 = np.zeros(G); y2 = np.zeros(G)
    cnt = np.zeros(G)
    for j in range(n_fold):
        idx = np.flatnonzero(F == j)
        for s in range(0, len(idx), CHUNK):
            r = idx[s:s + CHUNK]
            p = P[r].astype(np.float64); y = Y[r].astype(np.float64)
            m1 += p.sum(0) * cov[j]; m2 += (p * p).sum(0) * cov[j]
            y1 += y.sum(0) * cov[j]; y2 += (y * y).sum(0) * cov[j]
        cnt += len(idx) * cov[j]
    c = np.maximum(cnt, 1.0)
    sdP = np.sqrt(np.maximum(m2 / c - (m1 / c) ** 2, 0.0))
    sdY = np.sqrt(np.maximum(y2 / c - (y1 / c) ** 2, 0.0))
    g_mm = neutral(np.where(sdP > 0, sdY / np.maximum(sdP, 1e-300), 0.0), live & (sdP > 0))

    # target-label derangement, permuted WITHIN each inner fold so the coverage mask holds
    yperm = np.arange(N)
    for j in range(n_fold):
        idx = np.flatnonzero(F == j)
        yperm[idx] = idx[rng.permutation(len(idx))]
    numd, dend = accum_masked(P, Y, F, cov, n_fold, yperm)
    g_der = neutral(shrunk(numd, dend, shrink), dend > 0)

    # ---- null 4: 224's R4a shape -- a 12-bin gain conditioned on the RECIPIENT's basal
    b_right = np.asarray(fold["b_ctrl"], dtype=np.float64)
    wrong = [f for f in FOLDS if f != held][0]
    b_wrong = np.asarray(build_fold(wrong)["b_ctrl"], dtype=np.float64)
    g_b_right = neutral(binned(num, den, basal_bins(b_right, 12), 12, shrink), live)
    g_b_wrong = neutral(binned(num, den, basal_bins(b_wrong, 12), 12, shrink), live)
    # the 224-faithful restriction: reweight ONLY genes every fit line carries
    g_full = np.where(full_cov, g, float(np.median(g[live])))

    cands = dict(none=np.ones(G), gain=g, gain_reg=g_reg, gain_fullcov=g_full,
                 shuf=g_shuf, mm=g_mm, derange=g_der,
                 basal12=g_b_right, basal12_wrong=g_b_wrong)
    diag = {k: macro_cos_g(P, Y, v) for k, v in cands.items()}
    del P, Y

    out = dict(held=held, fit_lines=fit_lines, n_fit_targets=T,
               n_eval_targets=len(eval_targets), shrink=float(shrink),
               shrink_grid={f"{k:g}": float(v) for k, v in scores.items()},
               shrink_pooled=float(shrink_p),
               shrink_pooled_grid={f"{k:g}": float(v) for k, v in scores_p.items()},
               gain_reg_p10=float(np.quantile(g_reg[live], 0.10)),
               gain_reg_p50=float(np.quantile(g_reg[live], 0.50)),
               gain_reg_p90=float(np.quantile(g_reg[live], 0.90)),
               hold_line=fit_lines[hold], wrong_recipient=wrong,
               median_gain=float(np.median(g[live])),
               gain_p10=float(np.quantile(g[live], 0.10)),
               gain_p50=float(np.quantile(g[live], 0.50)),
               gain_p90=float(np.quantile(g[live], 0.90)),
               gain_spread_p90_p10=float(np.quantile(g[live], 0.90)
                                         / max(np.quantile(g[live], 0.10), 1e-300)),
               n_genes_live=int(live.sum()), n_genes_full_cov=int(full_cov.sum()),
               op=str(op_path.relative_to(ROOT)), in_sample_macro_cos=diag,
               in_sample_rel_pct={k: 100.0 * (v - diag["none"]) / abs(diag["none"])
                                  for k, v in diag.items()},
               seconds=time.time() - t0)
    np.savez_compressed(WORK / f"gain-{held}.npz", gain=g, gain_reg=g_reg,
                        gain_fullcov=g_full,
                        shuf=g_shuf, mm=g_mm, derange=g_der, basal12=g_b_right,
                        basal12_wrong=g_b_wrong, live=live, den=den,
                        full_cov=full_cov, cov_count=cov_count, meta=json.dumps(out))
    print(f"[{held}] gain p10/p50/p90 {out['gain_p10']:.4f}/{out['gain_p50']:.4f}/"
          f"{out['gain_p90']:.4f} over {out['n_genes_live']} live genes\n"
          f"[{held}] in-sample rel% " +
          "  ".join(f"{k}={v:+.2f}" for k, v in out["in_sample_rel_pct"].items()) +
          f"\n[{held}] {time.time() - t0:.0f}s", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--seed", type=int, default=SEED)
    args = ap.parse_args()
    out_path = WORK / "gain228.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else {}
    for held in args.folds:
        if (WORK / f"gain-{held}.npz").exists() and held in rows:
            print(f"[{held}] cached", flush=True)
            continue
        rows[held] = fit_one(held, args)
        out_path.write_text(json.dumps(rows, indent=2))
    out_path.write_text(json.dumps(rows, indent=2))
    print(f"\n-> {out_path}")


if __name__ == "__main__":
    main()
