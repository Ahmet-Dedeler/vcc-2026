"""027 — the per-gene emission variance channel, at count level.

Experiment 013 proposed it and could not test it; experiment 023 supplied the missing half of
the argument. A gene's Wilcoxon p-value depends on its mean effect **and** on its cell-to-cell
spread, and we set those independently. So:

  - genes we trust  -> emit at the variance floor -> significant, and first in `reach`'s order
  - genes we do not -> emit with realistic overdispersion -> weaker p-value, out of our set

**The mean profile is byte-identical either way**, so `pds`, `mse` and `nmae` should not move
while `fid` and `reach` do. That is the only lever found so far that improves one metric
through a degree of freedom the others cannot see.

Why a minimum-variance emitter calls essentially every gene, which is what this fixes: 400
near-identical values sit at the control distribution's *mean*, and for skewed counts the mean
is above the median, so the rank-sum against 18,400 real controls is significant even at zero
effect. Our `fid` of 0.4985 on the board is the sign accuracy of ~10,000 genes, i.e. chance.

The confidence signal is experiment 013's winner -- the geometric mean of the donors'
magnitudes, zeroed where they disagree on sign -- computed on the **un-whitened** effects,
because experiment 023 showed whitening makes `|effect|` a worse per-gene confidence signal
even as it improves the direction.


LODO scores `pds` exactly and everything else only by proxy. Whitening is a large `pds` gain
(+0.060, exp 019) but it reshapes the *relative magnitudes of genes inside each prediction*,
which is what the four DE metrics read through a Wilcoxon on emitted counts. The proxy in
experiment 018 says the cost is small; the proxy has no standing on `reach`, whose absolute
values it gets wrong by 4x.

So this runs the real `cell-eval2` stack on the `h1-fast` harness (40 perturbations, 2025 H1
hESC, its own measured baseline/replicate anchors). Two caveats carried from the log:

  - `pds` **saturates** on a 40-perturbation panel and must be ignored here.
  - `h1-fast` demonstrably lied about `jac` once (0.127 predicted, 0.0035 actual), and its
    `jac` seed noise is 6.7% of span against the competition's 0.8-2.2%. Read `jac` as
    directional at best.

`fid`, `reach` and `nmae` are what this run is for.

    uv run python experiments/027-variance-channel/run.py --gammas 0 0.65 --alpha 0.7 1.08 1.5
"""
from __future__ import annotations

import argparse
import itertools
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.donor_transfer import Donor, emit_counts, TS   # noqa: E402
from vcc2026.multi_donor import MultiDonorModel              # noqa: E402
from vcc2026.whiten import Whitener                          # noqa: E402

CONTROL, PERT_COL = "non-targeting", "target_gene"
DONORS = {"k562_gwps": "data/donors/k562_gwps",
          "orion_hct116": "data/donors/orion_hct116",
          "orion_hek293t": "data/donors/orion_hek293t"}


def confidence(model, perts, weights):
    """Per-(perturbation, gene) confidence: geometric mean of donor |effect|, zeroed on sign
    disagreement. Computed on un-whitened effects (exp 023)."""
    G = len(model.target_genes)
    out = np.zeros((len(perts), G))
    for i, p in enumerate(perts):
        have = [n for n in model.names if p in model.pert_index[n]]
        if not have:
            continue
        V = np.zeros((len(have), G))
        for j, n in enumerate(have):
            V[j, model.tcols[n]] = model.s[n][model.pert_index[n][p]][model.dcols[n]]
        agree = np.all(np.sign(V) == np.sign(V[0]), axis=0) & (V[0] != 0)
        out[i] = np.exp(np.log(np.abs(V) + 1e-12).mean(0)) * agree
    return out


def build(real_path, model, wcols, W, alpha, gamma, disp, mode, min_cpm, out, seed=0,
          tight_k=0, spread=0.0):
    a = ad.read_h5ad(real_path)
    X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
    labels = a.obs[PERT_COL].astype(str).to_numpy()
    genes = np.asarray(a.var_names)
    cm = labels == CONTROL
    Xc = X[cm].tocsr()
    Pc = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * Pc / Pc.sum())
    libs = np.asarray(Xc.sum(1)).ravel()
    perts = sorted(set(labels) - {CONTROL})
    n_per = int(pd.Series(labels[~cm]).value_counts().iloc[0])
    w = {n: 1.0 for n in model.names}

    E = np.vstack([model.effect(p, alpha, 0.0, weights=w) for p in perts]).astype(np.float64)
    if gamma > 0:
        before = float(np.median(np.linalg.norm(E, axis=1)))
        E[:, wcols] = W.apply(E[:, wcols], gamma, 0.0, renorm=False)
        E *= before / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)

    C = confidence(model, perts, w) if tight_k else None
    rng = np.random.default_rng(seed)
    data, indices, indptr, obs = [], [], [0], []
    for i, p in enumerate(perts):
        prof = model.profile_from_effect(E[i], b_ctrl, min_cpm=min_cpm, gate_unexpressed=True)
        pick = rng.choice(len(libs), n_per, replace=len(libs) < n_per)
        if tight_k:
            k = min(tight_k, int((C[i] > 0).sum()))
            tight = np.zeros(len(genes), bool)
            if k:
                tight[np.argpartition(-C[i], k - 1)[:k]] = True
            if spread > 0:
                spread_vec = np.where(tight, 0.0, spread)
                d, idx, ip = emit_counts(prof, libs[pick], disp, rng, mode="quasi",
                                         spread=spread_vec)
            else:
                d, idx, ip = emit_counts(prof, libs[pick], disp, rng, mode="mixed",
                                         tight_mask=tight)
        else:
            d, idx, ip = emit_counts(prof, libs[pick], disp, rng, mode=mode)
        data.append(d); indices.append(idx)
        indptr.extend(list(ip[1:] + indptr[-1]))
        obs += [p] * n_per
    data.append(Xc.data); indices.append(Xc.indices)
    indptr.extend(list(Xc.indptr[1:] + indptr[-1]))
    obs += [CONTROL] * Xc.shape[0]

    M = sp.csr_matrix((np.concatenate(data).astype(np.float32),
                       np.concatenate(indices).astype(np.int32),
                       np.asarray(indptr, dtype=np.int64)), shape=(len(obs), len(genes)))
    out.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=M, obs=pd.DataFrame({PERT_COL: obs},
                                     index=[f"c{i}" for i in range(len(obs))]),
               var=pd.DataFrame(index=genes)).write_h5ad(out, compression="lzf")
    return out


def score_one(t):
    tag, pred, real, bundle, wd = t
    run = Path(wd) / "run"
    r = subprocess.run(["uv", "run", "cell-eval2", "run", "-ap", pred, "-ar", real,
                        "-o", str(run), "--preset", "vcc2026", "--pert-col", PERT_COL,
                        "--control", CONTROL], capture_output=True, text=True)
    agg = run / "agg_results.csv"
    if r.returncode != 0 or not agg.exists():
        return tag, {"error": ("run: " + (r.stderr or r.stdout))[-800:]}
    out_csv = Path(wd) / "score.csv"
    sc = subprocess.run(["uv", "run", "cell-eval2", "score", "--user-agg", str(agg),
                         "--real-bundle", bundle, "-o", str(out_csv)],
                        capture_output=True, text=True)
    if sc.returncode != 0:
        return tag, {"error": ("score: " + (sc.stderr or sc.stdout))[-800:]}
    d = pd.read_csv(out_csv).set_index("metric")["from_replicate"]
    raw = pd.read_csv(agg).set_index("statistic")
    return tag, {"scaled": d.to_dict(),
                 "raw": {m: float(raw.loc["mean", m]) for m in raw.columns
                         if m != "statistic"}}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real", default="data/harness/h1-fast.h5ad")
    ap.add_argument("--bundle", default="experiments/002-calibration-sweep/fast/bundle")
    ap.add_argument("--outdir", default="experiments/027-variance-channel/fast")
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.0, 0.65])
    ap.add_argument("--alpha", type=float, nargs="+", default=[1.08])
    ap.add_argument("--dispersion", type=float, nargs="+", default=[0.5])
    ap.add_argument("--mode", nargs="+", default=["even"])
    ap.add_argument("--min-cpm", type=float, nargs="+", default=[5.0])
    ap.add_argument("--spread", type=float, nargs="+", default=[0.0],
                    help="fraction of Poisson noise given to LOW-confidence genes under "
                         "mode='quasi'. 0 reproduces the even emitter exactly.")
    ap.add_argument("--tight-k", type=int, nargs="+", default=[0],
                    help="genes per perturbation emitted at the variance floor; 0 = all "
                         "(the incumbent even emitter)")
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--seeds", type=int, nargs="+", default=[0],
                    help="emission seeds. `reach` responded non-monotonically to spread, so "
                         "any config kept has to survive a seed replicate before it counts.")
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    genes = np.asarray(ad.read_h5ad(args.real, backed="r").var_names)
    donors = {n: Donor.load(Path(p)) for n, p in DONORS.items()}
    model = MultiDonorModel(donors, genes, gene_norm=False, eb=False)
    shared = set(model.target_genes[model.tcols[model.names[0]]])
    for n in model.names[1:]:
        shared &= set(model.target_genes[model.tcols[n]])
    ti = {g: i for i, g in enumerate(model.target_genes)}
    wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
    print(f"whitening axis {len(wcols)} genes on the harness axis", flush=True)

    cache = outdir / "whitener_h1fast.npz"
    if cache.exists():
        W = Whitener.load(cache)
    else:
        rng = np.random.default_rng(0)
        rows = []
        for n in model.names:
            S = np.zeros((model.s[n].shape[0], len(model.target_genes)), np.float32)
            S[:, model.tcols[n]] = model.s[n][:, model.dcols[n]]
            if S.shape[0] > 6000:
                S = S[np.sort(rng.choice(S.shape[0], 6000, False))]
            rows.append(model._renorm[n] * S[:, wcols].astype(np.float64))
        B = np.vstack(rows)
        print(f"fitting whitener {B.shape}", flush=True)
        W = Whitener.fit(B); W.save(cache); del B

    jobs = []
    for g, al, di, mo, mc, tk, spr, sd in itertools.product(
            args.gammas, args.alpha, args.dispersion, args.mode, args.min_cpm,
            args.tight_k, args.spread, args.seeds):
        tag = f"g{g}_a{al}_d{di}_{mo}_cpm{mc:g}_k{tk}_s{spr:g}_seed{sd}"
        wd = outdir / "runs" / tag; wd.mkdir(parents=True, exist_ok=True)
        pred = wd / "pred.h5ad"
        if not pred.exists():
            print(f"  build {tag}", flush=True)
            build(Path(args.real), model, wcols, W, al, g, di, mo, mc, pred, seed=sd,
                  tight_k=tk, spread=spr)
        jobs.append((tag, str(pred), args.real, args.bundle, str(wd)))

    print(f"\nscoring {len(jobs)} configs, {args.jobs} workers\n", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(score_one, j): j[0] for j in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            tag, res = fut.result()
            res["tag"] = tag
            rows.append(res)
            if "error" in res:
                print(f"[{n}/{len(jobs)}] {tag}: {res['error'][:200]}", flush=True)
            else:
                r, s = res["raw"], res["scaled"]
                print(f"[{n}/{len(jobs)}] {tag}: avg={s.get('avg_score', float('nan')):.4f} "
                      f"fid={r['de_wilcoxon_direction_fidelity_yield_raw']:.4f} "
                      f"reach={r['de_wilcoxon_direction_reach_raw']:.4f} "
                      f"nmae={r['de_wilcoxon_lfc_nmae']:.4f} "
                      f"jac={r['de_wilcoxon_sig_jaccard']:.4f} "
                      f"mse={r['expr_mse_unbiased_capped_norm']:.4f}", flush=True)
            (outdir / "results.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
