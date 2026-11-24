#!/usr/bin/env python3
"""630 — A/B/C capacity pre-check for the prereg's capacity gate (each context's mean ur_frac_got >= 0.127, ur_n_moved >=
500/target), BEFORE any build: run the production module (src_rank_up.apply, b_k100) on the final blocks of a real A/B/C
package that carries zs (446bj2z), with each context's own scorer-order control block.

usage: uv run python experiments/630-rank-neutral-up/capacity630.py <restored.h5ad> [n_targets_per_context] [variant]
The restored h5ad comes from scripts/restore_prediction_h5ad.py on the package. Blocks are read per target (400 rows); the
ballast rows are identified by library (> 200,000 counts; the builder appends them after the carriers, checked).
Writes work/capacity630_<variant>.json. ~5 GB peak; run through gate602.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "src"))
import src_rank_up as M  # noqa: E402
from vcc2026 import rank_mag as ZS  # noqa: E402

path = Path(sys.argv[1])
n_per = int(sys.argv[2]) if len(sys.argv) > 2 else 60
variant = sys.argv[3] if len(sys.argv) > 3 else "b_k100"
genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
gpos = {g: i for i, g in enumerate(genes)}
pred = ad.read_h5ad(path, backed="r")
assert list(pred.var_names.astype(str)) == genes, "gene order differs from gene_names.csv"
obs = pred.obs[["target_gene", "context"]].astype(str)
out = {"package_h5ad": str(path), "variant": variant, "contexts": {}}
rng = np.random.default_rng(630)
for ctx in sorted(obs["context"].unique()):
    t0 = time.monotonic()
    cands = [p for p in (ROOT / "data/challenge").glob("context_*.h5ad")]
    cmat = None
    for cp in cands:
        a = ad.read_h5ad(cp, backed="r")
        if str(a.obs["context"].astype(str).iloc[0]) == ctx:
            cmat = a.X[:]
            break
    assert cmat is not None, f"no control file for context {ctx}"
    cc = ZS.scorer_control_block(cmat)
    tie = M.scorer_tie(cc)
    lib_med = float(np.median(np.asarray(cmat.sum(1)).ravel()))
    kpos = np.cumsum(np.asarray(cc.keep, dtype=bool)) - 1
    del cmat
    rows_ctx = np.flatnonzero((obs["context"] == ctx).to_numpy())
    targets = sorted(obs.iloc[rows_ctx]["target_gene"].unique())
    pick = sorted(rng.choice(len(targets), size=min(n_per, len(targets)), replace=False).tolist())
    recs = []
    for ti in pick:
        t = targets[ti]
        rows = rows_ctx[(obs.iloc[rows_ctx]["target_gene"] == t).to_numpy()]
        X = np.asarray(pred.X[rows[0]:rows[-1] + 1].toarray() if rows[-1] - rows[0] + 1 == rows.size
                       else pred.X[rows].toarray(), dtype=np.int64)
        lib = X.sum(1)
        isb = lib > 200_000
        K = int(isb.sum())
        if K and not isb[-K:].all():                 # put ballast rows last if the writer did not keep them there
            X = np.vstack([X[~isb], X[isb]])
        n_car = X.shape[0] - K
        full = gpos.get(t, -1)
        own = int(kpos[full]) if full >= 0 and bool(np.asarray(cc.keep, dtype=bool)[full]) else -1
        Xu, d = M.apply(X, n_car, cc, variant, tie=tie, lib_med=lib_med, own_col=own)
        assert np.array_equal(Xu.sum(0), X.sum(0)) and np.array_equal(Xu.sum(1), X.sum(1))
        recs.append(dict(target=t, K=K, ballast_last=bool(not K or isb[-K:].all()),
                         **{k: v for k, v in d.items() if k in ("ur_frac_got", "ur_n_moved", "ur_n_sel", "ur_n_reverted",
                                                                 "ur_n_no_candidate", "ur_max_dAUC", "ur_guard_excluded",
                                                                 "ur_max_abs_dlog10p", "ur_n_neg_donor_rows",
                                                                 "ur_n_no_pos_donor", "ur_n_mean_revert", "ur_guard_failed")}))
    fr = np.array([r["ur_frac_got"] for r in recs])
    nm = np.array([r["ur_n_moved"] for r in recs])
    out["contexts"][ctx] = dict(n_targets=len(recs), median_control_library=lib_med, K_median=float(np.median([r["K"] for r in recs])),
                                frac_got_mean=float(fr.mean()), frac_got_q10=float(np.quantile(fr, .1)),
                                n_moved_mean=float(nm.mean()), n_sel_mean=float(np.mean([r["ur_n_sel"] for r in recs])),
                                gate_frac_ok=bool(fr.mean() >= 0.127), gate_moved_ok=bool(nm.mean() >= 500),
                                max_dAUC=float(max(r["ur_max_dAUC"] for r in recs)),
                                guard_excluded_mean=float(np.mean([r["ur_guard_excluded"] for r in recs])),
                                # 605 fix counters (absent from pre-fix runs -> 0)
                                neg_donor_rows_sum=float(sum(r.get("ur_n_neg_donor_rows", 0) for r in recs)),
                                no_pos_donor_sum=float(sum(r.get("ur_n_no_pos_donor", 0) for r in recs)),
                                mean_revert_sum=float(sum(r.get("ur_n_mean_revert", 0) for r in recs)),
                                guard_failed_sum=float(sum(r.get("ur_guard_failed", 0) for r in recs)),
                                seconds=time.monotonic() - t0, per_target=recs)
    c = out["contexts"][ctx]
    print(f"{ctx}: {c['n_targets']} targets, frac_got mean {c['frac_got_mean']:.3f} (q10 {c['frac_got_q10']:.3f}), "
          f"moved {c['n_moved_mean']:.0f}/{c['n_sel_mean']:.0f}, K {c['K_median']:.0f}, gate frac {c['gate_frac_ok']} "
          f"moved {c['gate_moved_ok']}, max dAUC {c['max_dAUC']:.1e}, {c['seconds']:.0f}s; 605 counters: neg-donor rows "
          f"{c['neg_donor_rows_sum']:.0f}, no-pos-donor {c['no_pos_donor_sum']:.0f}, mean reverts {c['mean_revert_sum']:.0f}, "
          f"guard failed {c['guard_failed_sum']:.0f}", flush=True)
(HERE / f"work/capacity630_{variant}.json").write_text(json.dumps(out, indent=1, default=float))
