#!/usr/bin/env python3
"""803 item 5 — the real-scorer cross-check: write the cells of one engine run as an h5ad and compare the real cell-eval2's
per-target DE members with the engine's f 1 numbers (evaluate323's own DE + cell-eval2's member code).

  write  --surface h1test --arm zs --seed 0 --real data/harness/h1-test.h5ad --out <pred.h5ad>
      Reads the run's dumped pre-pass blocks (run496 --dump; de803 leaves them in the scratch dir when no rnup stage asked
      for them), re-applies the production zs pass exactly as 630/rnup630.py does (vcc2026.rank_mag.postpass with
      dn_k250_zs), asserts each block's stats equal the arm's stats file bit for bit, and writes the 85 x 400 cells plus the
      real file's own non-targeting cells (the scorer's controls) on the real file's gene axis. No emitter code is changed.
  compare --surface h1test --arm zs --seed 0 --run <cell-eval2 run dir>
      results.csv per target vs members.npz ce_* (and lib521 f 1); writes work/out/xcheck/<arm>-s<seed>.json.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE / "work" / "out"
sys.path.insert(0, str(HERE))
import de803  # noqa: E402

FULL = {"nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
        "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard", "pds": "pds_cosine",
        "mse": "expr_mse_unbiased_capped_norm"}


def cmd_write(a):
    import anndata as ad
    import scipy.sparse as sp
    sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
    import run400 as H
    from vcc2026 import rank_mag as ZS
    R, E = H.R, H.E
    R.WORK = de803.SURF[a.surface]["refdir"]
    man = json.loads((OUT / a.surface / a.arm / f"s{a.seed}" / "manifest.json").read_text())
    stats = ROOT / man["stats"]
    rundir = stats.parent
    zst = np.load(stats, allow_pickle=False)
    pre = [m for m in rundir.glob("metrics-*.json") if "pre_tag" not in json.loads(m.read_text())]
    tag = json.loads(pre[0].read_text())["tag"]
    files = sorted((rundir / "work" / "blocks").glob(f"{tag}_*.npz"), key=lambda f: int(re.search(r"_(\d+)\.npz$", f.name).group(1)))
    r, c, _ = E.load_reference()
    assert len(files) == len(r["targets"]), len(files)
    spec = ZS.parse_spec("dn_k250_zs")
    is_zs = man["arm"]["stage"] == "zs"
    blocks, labels = [], []
    for i, f in enumerate(files):
        z = np.load(f, allow_pickle=False)
        X = sp.csr_matrix((z["data"], z["indices"], z["indptr"]), shape=tuple(z["shape"])).toarray().astype(np.int64)
        assert str(z["target"]) == str(r["targets"][i])
        s, j, p, q, l = E.block_stats(X, c)
        if is_zs:
            X, _ = ZS.postpass(X, int(z["n_car"]), c, l, spec)
            s, j, p, q, l = E.block_stats(X, c)
        for nm, v in (("sums", s), ("p", p), ("q", q), ("lfc", l)):
            if not np.array_equal(np.asarray(v), zst[nm][i], equal_nan=True):
                raise AssertionError(f"target {i}: re-built block {nm} differs from the arm's stats")
        blocks.append(sp.csr_matrix(X.astype(np.float32)))
        labels += [str(r["targets"][i])] * X.shape[0]
    real = ad.read_h5ad(a.real)
    genes = np.asarray(r["genes"]).astype(str)
    assert list(real.var_names.astype(str)) == list(genes), "real gene axis != reference axis"
    ctrl = real[real.obs["target_gene"].astype(str) == "non-targeting"]
    Xc = sp.csr_matrix(ctrl.X).astype(np.float32)
    X = sp.vstack(blocks + [Xc]).tocsr()
    X.eliminate_zeros()
    obs = {"target_gene": labels + ["non-targeting"] * Xc.shape[0]}
    import pandas as pd
    out = ad.AnnData(X=X, obs=pd.DataFrame(obs, index=[f"c{k}" for k in range(X.shape[0])]), var=pd.DataFrame(index=genes))
    out.write_h5ad(a.out)
    print(f"[cells803] wrote {a.out}: {X.shape}, {X.nnz} nonzeros, blocks re-verified against {stats.name}")


def cmd_compare(a):
    import polars as pl
    z = np.load(OUT / a.surface / a.arm / f"s{a.seed}" / "members.npz")
    targets = [str(t) for t in z["targets"]]
    res = pl.read_csv(Path(a.run) / "results.csv")
    cols = res.columns
    pc = [c for c in cols if c in ("perturbation", "target_gene", "pert", "target")][0]
    rec = dict(arm=a.arm, seed=a.seed, run=a.run, members={})
    agg = pl.read_csv(Path(a.run) / "agg_results.csv")
    for m, fn in FULL.items():
        if "metric" in cols:
            sub = res.filter(pl.col("metric") == fn)
            d = dict(zip(sub[pc].cast(str).to_list(), sub["value"].to_list()))
        elif fn in cols:
            d = dict(zip(res[pc].cast(str).to_list(), res[fn].to_list()))
        else:
            continue
        real_ce = np.array([np.nan if d.get(t) is None else float(d[t]) for t in targets])
        agg_v = float(agg[fn][0]) if fn in agg.columns else None
        if m == "mse":
            rec["members"][m] = dict(cell_eval2_agg=agg_v, engine=float(np.nansum(z["mse_num"]) / np.nansum(z["mse_den"])))
            continue
        eng = z[f"ce_{m}"]
        both = np.isfinite(real_ce) & np.isfinite(eng)
        r = dict(cell_eval2_mean=float(np.nanmean(real_ce)), cell_eval2_agg=agg_v, engine_ce_mean=float(np.nanmean(eng)),
                 n_both=int(both.sum()), max_abs_target_diff=float(np.max(np.abs(real_ce[both] - eng[both]), initial=0.0)),
                 n_target_diff_gt_1e6=int((np.abs(real_ce[both] - eng[both]) > 1e-6).sum()),
                 nan_mismatch=int((np.isfinite(real_ce) != np.isfinite(eng)).sum()))
        if m != "pds":
            r["engine_lib521_f1_mean"] = float(np.nanmean(z[f"{m}_f1"]))
        rec["members"][m] = r
    od = OUT / "xcheck"
    od.mkdir(parents=True, exist_ok=True)
    (od / f"{a.surface}-{a.arm}-s{a.seed}.json").write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    for n in ("write", "compare"):
        p = sub.add_parser(n)
        p.add_argument("--surface", default="h1test")
        p.add_argument("--arm", default="zs")
        p.add_argument("--seed", type=int, default=0)
        if n == "write":
            p.add_argument("--real", required=True)
            p.add_argument("--out", required=True)
        else:
            p.add_argument("--run", required=True)
    a = ap.parse_args()
    {"write": cmd_write, "compare": cmd_compare}[a.cmd](a)


if __name__ == "__main__":
    main()
