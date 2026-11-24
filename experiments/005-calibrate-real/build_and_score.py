"""Experiment 005 — calibrate the real donor-transfer model against the real scorer.

Experiment 002 swept (alpha, dispersion) against a *synthetic* direction with a controlled
cosine. This does it with the actual K562-gwps donor transfer, whose measured cosine on this
target is 0.073 (exp 004) rather than the 0.15-0.5 that sweep explored. The optimum moves
with cos, so the previous grid's answer does not carry over.

Three knobs:
  alpha   scale on the gene-specific effect `s` -- pds-invariant, drives mse/fid/reach/nmae
  beta    scale on the donor's mean response `m` -- earns nothing directly, but shifts the
          profile toward where real perturbed cells sit, which can move mse and the DE calls
  disp    NB emission dispersion -- sets how many genes clear the Wilcoxon/BH bar

    uv run python experiments/005-calibrate-real/build_and_score.py --jobs 5 \
        --alpha 0.5 1.0 2.0 --beta 0 1 --dispersion 0.1 0.5
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

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from vcc2026.donor_transfer import Donor, TransferModel, emit_counts, TS  # noqa: E402

CONTROL, PERT_COL = "non-targeting", "target_gene"


def build_prediction(ref_path: Path, model: TransferModel, alpha: float, beta: float,
                     disp: float, out: Path, seed: int = 0,
                     gate: bool = True, space: str = "log1p",
                     mode: str = "nb", min_cpm: float = 20.0) -> Path:
    a = ad.read_h5ad(ref_path)
    X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
    labels = a.obs[PERT_COL].astype(str).to_numpy()
    genes = np.asarray(a.var_names)

    ctrl_mask = labels == CONTROL
    Xc = X[ctrl_mask]
    Pc = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * Pc / Pc.sum())
    libs = np.asarray(Xc.sum(1)).ravel()

    rng = np.random.default_rng(seed)
    perts = sorted(set(labels) - {CONTROL})
    n_per = int(pd.Series(labels[~ctrl_mask]).value_counts().iloc[0])

    data, indices, indptr, obs = [], [], [0], []
    for p in perts:
        prof = model.profile(p, b_ctrl, alpha, beta, space=space,
                             gate_unexpressed=gate, min_cpm=min_cpm)
        pick = rng.choice(len(libs), n_per, replace=len(libs) < n_per)
        d, i, ip = emit_counts(prof, libs[pick], disp, rng, mode=mode)
        data.append(d); indices.append(i)
        indptr.append(indptr[-1]); indptr[-1:] = list(ip[1:] + indptr[-1])
        obs += [p] * n_per

    # Carry the reference's control cells through: the scorer uses them as the origin on
    # both sides, and the vcc2026 contract expects the control label present.
    Xc = Xc.tocsr()
    data.append(Xc.data); indices.append(Xc.indices)
    indptr.extend(list(Xc.indptr[1:] + indptr[-1]))
    obs += [CONTROL] * Xc.shape[0]

    M = sp.csr_matrix((np.concatenate(data).astype(np.float32),
                       np.concatenate(indices).astype(np.int32),
                       np.asarray(indptr, dtype=np.int64)),
                      shape=(len(obs), len(genes)))
    out.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=M,
               obs=pd.DataFrame({PERT_COL: obs}, index=[f"c{i}" for i in range(len(obs))]),
               var=pd.DataFrame(index=genes)).write_h5ad(out, compression="lzf")
    return out


def score_one(args_tuple):
    tag, pred, real, bundle, wd = args_tuple
    wd = Path(wd); run = wd / "run"
    r = subprocess.run(["uv", "run", "cell-eval2", "run", "-ap", pred, "-ar", real,
                        "-o", str(run), "--preset", "vcc2026",
                        "--pert-col", PERT_COL, "--control", CONTROL],
                       capture_output=True, text=True)
    if r.returncode != 0:
        return tag, {"error": ("run: " + (r.stderr or r.stdout))[-1200:]}
    agg = run / "agg_results.csv"
    if not agg.exists():
        return tag, {"error": f"no agg in {run}"}
    out_csv = wd / "score.csv"
    sc = subprocess.run(["uv", "run", "cell-eval2", "score", "--user-agg", str(agg),
                         "--real-bundle", bundle, "-o", str(out_csv)],
                        capture_output=True, text=True)
    if sc.returncode != 0:
        return tag, {"error": ("score: " + (sc.stderr or sc.stdout))[-1200:]}
    d = pd.read_csv(out_csv).set_index("metric")["from_replicate"]
    raw = pd.read_csv(agg).set_index("statistic")
    return tag, {"scaled": d.to_dict(),
                 "raw": {m: float(raw.loc["mean", m]) for m in raw.columns if m != "statistic"}}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real", default="data/harness/h1-fast.h5ad")
    ap.add_argument("--bundle", default="experiments/002-calibration-sweep/fast/bundle")
    ap.add_argument("--donor", default="data/donors/k562_gwps")
    ap.add_argument("--outdir", default="experiments/005-calibrate-real/fast")
    ap.add_argument("--alpha", type=float, nargs="+", default=[1.0])
    ap.add_argument("--beta", type=float, nargs="+", default=[0.0])
    ap.add_argument("--dispersion", type=float, nargs="+", default=[0.1])
    ap.add_argument("--mode", nargs="+", default=["nb"], choices=["nb", "even"],
                    help="emission model; 'even' is minimum-variance allocation")
    ap.add_argument("--min-cpm", type=float, default=20.0)
    ap.add_argument("--jobs", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    outdir = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    genes = np.asarray(ad.read_h5ad(args.real, backed="r").var_names)
    model = TransferModel(Donor.load(Path(args.donor)), genes)
    print(f"donor gene coverage {model.coverage:.1%}", flush=True)

    jobs = []
    for al, be, di, mo in itertools.product(args.alpha, args.beta, args.dispersion, args.mode):
        tag = f"a{al}_b{be}_d{di}" + ("" if mo == "nb" else f"_{mo}")
        wd = outdir / "runs" / tag; wd.mkdir(parents=True, exist_ok=True)
        pred = wd / "pred.h5ad"
        if not pred.exists():
            print(f"  build {tag}", flush=True)
            build_prediction(Path(args.real), model, al, be, di, pred, seed=args.seed,
                             mode=mo, min_cpm=args.min_cpm)
        jobs.append((tag, str(pred), args.real, args.bundle, str(wd)))

    print(f"\nscoring {len(jobs)} configs, {args.jobs} workers\n", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(score_one, j): j[0] for j in jobs}
        for n, fut in enumerate(as_completed(futs), 1):
            tag, res = fut.result()
            res["tag"] = tag
            rows.append(res)
            avg = res.get("scaled", {}).get("avg_score")
            print(f"[{n}/{len(jobs)}] {tag}: "
                  f"{res.get('error','')[:150] or f'avg={avg:.4f}'}", flush=True)
            (outdir / "results.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {outdir/'results.json'}")


if __name__ == "__main__":
    main()
