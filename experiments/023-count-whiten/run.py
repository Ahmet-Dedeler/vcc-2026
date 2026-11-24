"""023 — count-level check of whitening, on the metrics LODO cannot see.

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

    uv run python experiments/023-count-whiten/run.py --gammas 0 0.65 --alpha 0.7 1.08 1.5
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


def build(real_path, model, wcols, W, alpha, gamma, disp, mode, min_cpm, out, seed=0):
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

    rng = np.random.default_rng(seed)
    data, indices, indptr, obs = [], [], [0], []
    for i, p in enumerate(perts):
        prof = model.profile_from_effect(E[i], b_ctrl, min_cpm=min_cpm, gate_unexpressed=True)
        pick = rng.choice(len(libs), n_per, replace=len(libs) < n_per)
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
    ap.add_argument("--outdir", default="experiments/023-count-whiten/fast")
    ap.add_argument("--gammas", type=float, nargs="+", default=[0.0, 0.65])
    ap.add_argument("--alpha", type=float, nargs="+", default=[1.08])
    ap.add_argument("--dispersion", type=float, nargs="+", default=[0.1])
    ap.add_argument("--mode", nargs="+", default=["even"])
    ap.add_argument("--min-cpm", type=float, nargs="+", default=[20.0])
    ap.add_argument("--jobs", type=int, default=6)
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
    for g, al, di, mo, mc in itertools.product(args.gammas, args.alpha, args.dispersion,
                                               args.mode, args.min_cpm):
        tag = f"g{g}_a{al}_d{di}_{mo}_cpm{mc:g}"
        wd = outdir / "runs" / tag; wd.mkdir(parents=True, exist_ok=True)
        pred = wd / "pred.h5ad"
        if not pred.exists():
            print(f"  build {tag}", flush=True)
            build(Path(args.real), model, wcols, W, al, g, di, mo, mc, pred)
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
