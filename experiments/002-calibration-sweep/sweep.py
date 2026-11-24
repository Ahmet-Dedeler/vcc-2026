"""Experiment 002 — map directional accuracy to achievable score.

Experiment 001 showed analytically that predicted-effect magnitude `alpha` is free for
`pds_cosine` and decisive for `expr_mse`. This runs the real scorer to (a) confirm that
end to end, (b) find the joint optimum over `alpha` and emission `dispersion`, and
(c) produce the efficient frontier: **what overall score is reachable at each level of
directional accuracy?** That last one tells us how good the model actually has to be.

Method. Take the harness reference's true per-perturbation log-pseudobulk effect
`delta_p = b_p - b_ctrl`. Build a synthetic prediction direction with a *controlled*
cosine similarity to it:

    d_p = delta_p + t * noise_p,   noise orthogonal to delta,  t = sqrt(1/cos^2 - 1)

then emit counts from `b_ctrl + alpha * d_p`. Sweeping cos isolates directional accuracy
from every other modelling choice, so the resulting table is a property of the *metric
suite*, not of any particular model.

    uv run python experiments/002-calibration-sweep/sweep.py --smoke
    uv run python experiments/002-calibration-sweep/sweep.py --cos 0.3 0.5 0.7 --alpha 0.1 0.25 0.5 1.0
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

TS = 5e4
CONTROL = "non-targeting"
PERT_COL = "target_gene"
HERE = Path(__file__).parent


# ---------------------------------------------------------------- reference


def load_reference(path: Path):
    """Return (control counts matrix, b_ctrl, {pert: delta}, gene names, lib sizes)."""
    a = ad.read_h5ad(path)
    X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
    labels = a.obs[PERT_COL].to_numpy()

    is_ctrl = labels == CONTROL
    Pc = np.asarray(X[is_ctrl].sum(0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * Pc / Pc.sum())

    deltas = {}
    for g in pd.unique(labels):
        if g == CONTROL:
            continue
        P = np.asarray(X[labels == g].sum(0)).ravel().astype(np.float64)
        deltas[g] = np.log1p(TS * P / P.sum()) - b_ctrl

    libs = np.asarray(X[is_ctrl].sum(1)).ravel()
    return X[is_ctrl], b_ctrl, deltas, np.asarray(a.var_names), libs


def make_direction(delta: np.ndarray, cos_target: float, mask: np.ndarray,
                   rng: np.random.Generator) -> np.ndarray:
    """A vector whose cosine similarity with `delta` is exactly `cos_target`."""
    if cos_target >= 0.999:
        return delta.copy()
    noise = rng.normal(0, 1, delta.size) * mask
    noise -= noise.dot(delta) / delta.dot(delta) * delta
    n = np.linalg.norm(noise)
    if n == 0:
        return delta.copy()
    noise *= np.linalg.norm(delta) / n
    return delta + np.sqrt(1.0 / cos_target**2 - 1.0) * noise


# ---------------------------------------------------------------- emission


def emit(profile_p: np.ndarray, libs: np.ndarray, dispersion: float | None,
         rng: np.random.Generator) -> sp.csr_matrix:
    """Sample `len(libs)` count vectors from a relative-expression profile."""
    n, G = len(libs), profile_p.size
    data, indices, indptr = [], [], [0]
    for lib in libs:
        mu = profile_p * float(lib)
        if dispersion:
            mu = rng.gamma(shape=1.0 / dispersion, scale=mu * dispersion)
        c = rng.poisson(mu)
        idx = np.flatnonzero(c)
        data.append(c[idx].astype(np.float32))
        indices.append(idx.astype(np.int32))
        indptr.append(indptr[-1] + idx.size)
    return sp.csr_matrix(
        (np.concatenate(data), np.concatenate(indices), np.array(indptr)), shape=(n, G)
    )


def build_prediction(ref, cos_target: float, alpha: float, dispersion: float | None,
                     out: Path, cells_per_pert: int, seed: int) -> Path:
    """Write a prediction h5ad with the reference's control cells carried through."""
    Xc, b_ctrl, deltas, genes, libs = ref
    rng = np.random.default_rng(seed)
    mask = (b_ctrl > 0).astype(float)

    blocks, labels = [], []
    for g, delta in deltas.items():
        d = make_direction(delta, cos_target, mask, rng)
        b_hat = b_ctrl + alpha * d
        p = np.clip(np.expm1(b_hat) / TS, 0, None)
        p /= p.sum()
        pick = rng.choice(len(libs), cells_per_pert, replace=len(libs) < cells_per_pert)
        blocks.append(emit(p, libs[pick], dispersion, rng))
        labels += [g] * cells_per_pert

    # The scorer treats the reference's controls as the origin on both sides, and the
    # vcc2026 contract says the prediction's labels match the reference "including the
    # non-targeting control label". Carry them through unchanged.
    blocks.append(Xc)
    labels += [CONTROL] * Xc.shape[0]

    a = ad.AnnData(
        X=sp.vstack(blocks).tocsr(),
        obs=pd.DataFrame({PERT_COL: labels}, index=[f"c{i}" for i in range(len(labels))]),
        var=pd.DataFrame(index=genes),
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out, compression="lzf")  # lzf: ~10x faster than gzip, files are transient
    return out


# ---------------------------------------------------------------- scoring


def score(pred: Path, real: Path, bundle: Path, workdir: Path) -> dict:
    """One scoring job. cell-eval2 is single-threaded, so parallelism lives here,
    across configs, rather than inside a run."""
    run = workdir / "run"
    r = subprocess.run(
        ["uv", "run", "cell-eval2", "run", "-ap", str(pred), "-ar", str(real),
         "-o", str(run), "--preset", "vcc2026",
         "--pert-col", PERT_COL, "--control", CONTROL],
        capture_output=True, text=True,
    )
    if r.returncode != 0:
        return {"error": ("run: " + (r.stderr or r.stdout))[-1500:]}

    agg = next(iter(sorted(run.glob("*agg*.csv"))), None)
    if agg is None:
        return {"error": f"no agg csv in {run}; files={[p.name for p in run.iterdir()]}"}

    out_csv = workdir / "score.csv"
    sc = subprocess.run(
        ["uv", "run", "cell-eval2", "score", "--user-agg", str(agg),
         "--real-bundle", str(bundle), "-o", str(out_csv)],
        capture_output=True, text=True,
    )
    if sc.returncode != 0:
        return {"error": ("score: " + (sc.stderr or sc.stdout))[-1500:], "agg": str(agg)}

    res = {"agg": str(agg), "score_csv": str(out_csv)}
    if out_csv.exists():
        df = pd.read_csv(out_csv)
        res["score_table"] = df.to_dict("records")
    res["stdout"] = sc.stdout[-1500:]
    return res


def _job(args_tuple):
    """Top-level so ProcessPoolExecutor can pickle it."""
    tag, pred, real, bundle, wd = args_tuple
    out = score(Path(pred), Path(real), Path(bundle), Path(wd))
    return tag, out


# ---------------------------------------------------------------- driver


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--real", default="data/harness/h1-test.h5ad")
    ap.add_argument("--bundle", default="experiments/002-calibration-sweep/out/bundle")
    ap.add_argument("--outdir", default="experiments/002-calibration-sweep/out")
    ap.add_argument("--cos", type=float, nargs="+", default=[0.5])
    ap.add_argument("--alpha", type=float, nargs="+", default=[1.0])
    ap.add_argument("--dispersion", type=float, nargs="+", default=[0.5])
    ap.add_argument("--cells-per-pert", type=int, default=400)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=6, help="concurrent scoring processes")
    ap.add_argument("--smoke", action="store_true", help="one tiny config, keep artefacts")
    args = ap.parse_args()

    if args.smoke:
        args.cos, args.alpha, args.dispersion, args.cells_per_pert = [0.5], [0.25], [0.5], 100
        args.jobs = 1

    outdir = Path(args.outdir)
    print(f"loading reference {args.real} ...", flush=True)
    ref = load_reference(Path(args.real))
    print(f"  {len(ref[2])} perturbations, {ref[1].size} genes, "
          f"{ref[0].shape[0]} control cells", flush=True)

    grid = list(itertools.product(args.cos, args.alpha, args.dispersion))
    jobs, meta = [], {}
    for cos, alpha, disp in grid:
        tag = f"cos{cos}_a{alpha}_d{disp}"
        wd = outdir / "runs" / tag
        wd.mkdir(parents=True, exist_ok=True)
        pred = wd / "pred.h5ad"
        if not pred.exists():
            print(f"  build {tag}", flush=True)
            build_prediction(ref, cos, alpha, disp, pred, args.cells_per_pert, args.seed)
        jobs.append((tag, str(pred), args.real, args.bundle, str(wd)))
        meta[tag] = dict(cos=cos, alpha=alpha, dispersion=disp)

    print(f"\nscoring {len(jobs)} configs with {args.jobs} workers "
          f"(~7 min each, single-threaded)\n", flush=True)

    rows, done = [], 0
    with ProcessPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(_job, j): j[0] for j in jobs}
        for fut in as_completed(futs):
            tag, res = fut.result()
            res.update(meta[tag], tag=tag)
            rows.append(res)
            done += 1
            brief = res.get("error", "")[:200] or _brief(res)
            print(f"[{done}/{len(jobs)}] {tag}: {brief}", flush=True)
            (outdir / "sweep_results.json").write_text(json.dumps(rows, indent=2))

    (outdir / "sweep_results.json").write_text(json.dumps(rows, indent=2))
    print(f"\nwrote {outdir / 'sweep_results.json'}")


def _brief(res: dict) -> str:
    t = res.get("score_table")
    if not t:
        return "(no score table)"
    keys = [k for k in t[0] if "score" in k.lower() or "avg" in k.lower()]
    return " ".join(f"{k}={t[0][k]}" for k in keys[:4])


if __name__ == "__main__":
    sys.exit(main())
