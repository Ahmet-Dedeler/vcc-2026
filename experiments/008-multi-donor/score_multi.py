"""Experiment 008 — score multi-donor transfer against the real scorer.

Offline cosine says 2:1 K562-gwps : Orion-HCT116 beats K562 alone by 14.1% (0.0816 vs
0.0715). Cosine is a good proxy but not the objective, and experiment 006 already produced
one case where an intervention that raised cosine damaged a scored metric. So this runs the
same configs through `cell-eval2`.

    uv run python experiments/008-multi-donor/score_multi.py --weights 1:0 2:1 1:1 --jobs 3
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, emit_counts, TS  # noqa: E402
from vcc2026.multi_donor import MultiDonorModel  # noqa: E402

CONTROL, PERT_COL = "non-targeting", "target_gene"
DONORS = {"k562_gwps": "data/donors/k562_gwps",
          "orion_hct116": "data/donors/orion_hct116"}
# gene_norm helps a dense donor axis and hurts a sparse one -- see MultiDonorModel.
GENE_NORM = {"k562_gwps": True, "orion_hct116": False}


def build(ref_path: Path, model: MultiDonorModel, weights: dict, alpha: float,
          beta: float, disp: float, out: Path, mode: str, min_cpm: float,
          seed: int = 0) -> Path:
    a = ad.read_h5ad(ref_path)
    X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
    labels = a.obs[PERT_COL].astype(str).to_numpy()
    genes = np.asarray(a.var_names)

    ctrl = labels == CONTROL
    Xc = X[ctrl].tocsr()
    P = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * P / P.sum())
    libs = np.asarray(Xc.sum(1)).ravel()

    rng = np.random.default_rng(seed)
    perts = sorted(set(labels) - {CONTROL})
    n_per = int(pd.Series(labels[~ctrl]).value_counts().iloc[0])

    data, indices, indptr, obs = [], [], [0], []
    for p in perts:
        prof = model.profile(p, b_ctrl, alpha, beta, weights=weights, min_cpm=min_cpm)
        pick = rng.choice(len(libs), n_per, replace=len(libs) < n_per)
        d, i, ip = emit_counts(prof, libs[pick], disp, rng, mode=mode)
        data.append(d); indices.append(i)
        indptr.extend(list(ip[1:] + indptr[-1]))
        obs += [p] * n_per

    data.append(Xc.data); indices.append(Xc.indices)
    indptr.extend(list(Xc.indptr[1:] + indptr[-1]))
    obs += [CONTROL] * Xc.shape[0]

    M = sp.csr_matrix((np.concatenate(data).astype(np.float32),
                       np.concatenate(indices).astype(np.int32),
                       np.asarray(indptr, dtype=np.int64)),
                      shape=(len(obs), len(genes)))
    out.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=M, obs=pd.DataFrame({PERT_COL: obs},
                                     index=[f"c{i}" for i in range(len(obs))]),
               var=pd.DataFrame(index=genes)).write_h5ad(out, compression="lzf")
    return out


def score_one(t):
    tag, pred, real, bundle, wd = t
    wd = Path(wd); run = wd / "run"
    r = subprocess.run(["uv", "run", "cell-eval2", "run", "-ap", pred, "-ar", real,
                        "-o", str(run), "--preset", "vcc2026",
                        "--pert-col", PERT_COL, "--control", CONTROL],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return tag, {"error": (r.stderr or r.stdout)[-1000:]}
    agg = run / "agg_results.csv"
    if not agg.exists():
        return tag, {"error": "no agg"}
    out_csv = wd / "score.csv"
    sc = subprocess.run(["uv", "run", "cell-eval2", "score", "--user-agg", str(agg),
                         "--real-bundle", bundle, "-o", str(out_csv)],
                        capture_output=True, text=True)
    if sc.returncode != 0:
        return tag, {"error": (sc.stderr or sc.stdout)[-1000:]}
    raw = pd.read_csv(agg).set_index("statistic")
    return tag, {"scaled": pd.read_csv(out_csv).set_index("metric")["from_replicate"].to_dict(),
                 "raw": {m: float(raw.loc["mean", m]) for m in raw.columns if m != "statistic"}}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real", default="data/harness/h1-fast.h5ad")
    ap.add_argument("--bundle", default="experiments/002-calibration-sweep/fast/bundle")
    ap.add_argument("--outdir", default="experiments/008-multi-donor/out")
    ap.add_argument("--weights", nargs="+", default=["1:0", "2:1"],
                    help="k562:orion weight ratios")
    ap.add_argument("--alpha", type=float, default=0.7)
    ap.add_argument("--beta", type=float, default=0.0)
    ap.add_argument("--dispersion", type=float, default=0.1)
    ap.add_argument("--mode", default="even")
    ap.add_argument("--min-cpm", type=float, default=20.0)
    ap.add_argument("--jobs", type=int, default=3)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    genes = np.asarray(ad.read_h5ad(args.real, backed="r").var_names)
    donors = {n: Donor.load(p) for n, p in DONORS.items() if Path(p).exists()}
    print(f"donors: {list(donors)}", flush=True)
    model = MultiDonorModel(donors, genes, gene_norm=GENE_NORM)

    jobs = []
    for w in args.weights:
        k, o = (float(x) for x in w.split(":"))
        tag = f"w{w.replace(':','-')}_a{args.alpha}_{args.mode}"
        wd = outdir / "runs" / tag; wd.mkdir(parents=True, exist_ok=True)
        pred = wd / "pred.h5ad"
        if not pred.exists():
            print(f"  build {tag}", flush=True)
            build(Path(args.real), model, {"k562_gwps": k, "orion_hct116": o},
                  args.alpha, args.beta, args.dispersion, pred, args.mode, args.min_cpm)
        jobs.append((tag, str(pred), args.real, args.bundle, str(wd)))

    print(f"\nscoring {len(jobs)} configs\n", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(score_one, j): j[0] for j in jobs}
        for n, f in enumerate(as_completed(futs), 1):
            tag, res = f.result(); res["tag"] = tag; rows.append(res)
            if "error" in res:
                msg = res["error"][:150]
            else:
                msg = f"avg={res['scaled']['avg_score']:.4f}"
            print(f"[{n}/{len(jobs)}] {tag}: {msg}", flush=True)
            (outdir / "results.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
