#!/usr/bin/env python3
"""700 — blend PIE's predicted response directions into a production bulk dump (VCC_BULK_LOAD format), norm-preserving.

PIE's lfc_pred per (context, target, gene) is mapped into the scorer's pseudobulk effect space with the context's own
control pseudobulk: p_ctrl = expm1(b_ctrl) / TS (renormalised), p_pert = p_ctrl * 2**lfc (renormalised), effect =
log1p(TS p_pert) - log1p(TS p_ctrl), TS = 5e4. Genes PIE does not predict (no gene_text) get lfc 0. The shared PIE row is
the mean over the contexts given (A/B/C); the per-context rows are written too, for 672's VCC_BULK_LOAD_CTX.

Blend, per target, on non-panel genes (the panel columns, incl. the row's own knockdown, keep the base values, as the
builder's sharpening does): u = base / |base|, v = pie / |pie|, new = |base| * normalise((1 - w) u + w v). Every row keeps
its non-panel norm, so the LODO key's per-target magnitudes, ||bulk|| / ||carriers|| (the tercile cut variable) and the
energy split are unchanged; only directions move. `--carriers` applies the same blend to the carriers (the DE design),
else they are copied unchanged.

usage: .venv/bin/python make_pie_dump700.py --base <dump.npz> --pred <predictions.parquet> --dataset vcc26abc
         --ctx A=jurkat --ctx B=CVCL_0030 --ctx C=CVCL_1108 --bundle data/challenge --w 0.5 [--carriers] --out <npz>
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pyarrow.parquet as pq

TS = 5e4


def b_ctrl_of(bundle: Path, ctx: str, genes: np.ndarray) -> np.ndarray:
    a = ad.read_h5ad(bundle / f"context_{ctx}.h5ad")
    assert list(map(str, a.var_names)) == list(map(str, genes)), "bundle gene order differs from the dump"
    P = np.asarray(a.X.sum(0)).ravel().astype(np.float64)
    return np.log1p(TS * P / P.sum())


MIN_CPM = 5.0   # 700 amendment 2: PIE's lfc head is supervised only on tested (>= 5 CPM) DE genes; below that it is
                # unconstrained (|lfc| median 8.4, p95 26 on A/B/C) and blows the effect row up, so it is zeroed there.


def pie_effect(lfc: np.ndarray, bc: np.ndarray) -> np.ndarray:
    p = np.expm1(bc)
    p = p / p.sum()
    lfc = np.where(1e6 * p >= MIN_CPM, lfc, 0.0)
    q = p * np.exp2(lfc)
    q = q / q.sum()
    return np.log1p(TS * q) - np.log1p(TS * p)


def blend(base: np.ndarray, pie: np.ndarray, w: float, nonpanel: np.ndarray) -> np.ndarray:
    out = base.copy()
    for i in range(base.shape[0]):
        b, v = base[i, nonpanel], pie[i, nonpanel]
        nb, nv = np.linalg.norm(b), np.linalg.norm(v)
        if nb == 0 or nv == 0:
            continue
        m = (1.0 - w) * b / nb + w * v / nv
        nm = np.linalg.norm(m)
        if nm > 0:
            out[i, nonpanel] = nb * m / nm
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True)
    ap.add_argument("--pred", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--ctx", action="append", required=True, help="bundle_context=pie_context")
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--w", type=float, required=True)
    ap.add_argument("--carriers", action="store_true")
    ap.add_argument("--seen", default=None, help="JSON list of PIE train-split perturbations: blend only these targets")
    ap.add_argument("--center", action="store_true",
                    help="702c: remove each gene's mean over the panel's targets from PIE's rows before blending")
    ap.add_argument("--ctx-out", default=None,
                    help="also write a 672-format per-context file (perts, genes, contexts, bulk_<ctx>): each context's "
                         "base row blended with THAT context's PIE row (same w, seen mask, centring)")
    ap.add_argument("--keep-env", action="store_true",
                    help="write the base dump's env unchanged (def_day's dump hook requires perts/genes/env identical); "
                         "provenance goes to <out>.pie700.json instead")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    z = np.load(a.base, allow_pickle=False)
    base, carriers = np.asarray(z["bulk"], np.float64), np.asarray(z["carriers"], np.float64)
    perts, genes = z["perts"].astype(str), z["genes"].astype(str)
    env = json.loads(str(z["env"]))
    t = pq.read_table(a.pred)
    pg = json.loads(t.schema.metadata[b"pie"])["genes"][a.dataset]
    gpos = {g: i for i, g in enumerate(pg)}
    have = np.array([g in gpos for g in genes])
    gix = np.array([gpos.get(g, 0) for g in genes])
    df = t.to_pandas()
    df = df[df["dataset"] == a.dataset]
    panel = set(perts)
    nonpanel = np.array([g not in panel for g in genes])
    per_ctx = {}
    for spec in a.ctx:
        bctx, pctx = spec.split("=")
        bc = b_ctrl_of(Path(a.bundle), bctx, genes)
        sub = df[df["context"] == pctx].set_index("perturbation")
        E = np.zeros_like(base)
        for i, p in enumerate(perts):
            if p not in sub.index:
                continue
            lf = np.where(have, np.asarray(sub.loc[p, "lfc_pred"], np.float64)[gix], 0.0)
            E[i] = pie_effect(lf, bc)
        per_ctx[bctx] = E
        print(f"{bctx} ({pctx}): {int((np.abs(E).sum(1) > 0).sum())}/{len(perts)} targets, median |row| "
              f"{float(np.median(np.linalg.norm(E[:, nonpanel], axis=1))):.3f}", flush=True)
    shared = np.mean(np.stack(list(per_ctx.values())), axis=0)
    # the seen mask (700): targets with same-perturbation evidence in PIE's training split; a function of (panel, split)
    seen_mask = np.ones(len(perts), dtype=bool)
    if a.seen:
        seen = set(json.load(open(a.seen)))
        seen_mask = np.array([p in seen for p in perts])
    if a.center:
        shared = shared - shared.mean(axis=0, keepdims=True)
        per_ctx = {c: v - v.mean(axis=0, keepdims=True) for c, v in per_ctx.items()}
    shared_m = np.where(seen_mask[:, None], shared, 0.0)    # zero rows are skipped by blend(): the base row is kept
    new_bulk = blend(base, shared_m, a.w, nonpanel)
    new_car = blend(carriers, shared_m, a.w, nonpanel) if a.carriers else carriers
    # direction diagnostics: cosine of base vs PIE vs blend (non-panel)
    def rowcos(x, y):
        x, y = x[:, nonpanel], y[:, nonpanel]
        return (x * y).sum(1) / np.maximum(np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1), 1e-300)
    diag = {"w": a.w, "center": a.center, "carriers_blended": a.carriers, "seen_rule": a.seen, "n_seen": int(seen_mask.sum()), "cos_base_pie_median": float(np.median(rowcos(base, shared))),
            "cos_base_new_median": float(np.median(rowcos(base, new_bulk))),
            "norm_ratio_max_dev": float(np.max(np.abs(np.linalg.norm(new_bulk[:, nonpanel], axis=1) /
                                                       np.maximum(np.linalg.norm(base[:, nonpanel], axis=1), 1e-300) - 1)))}
    env2 = dict(env)
    env2["PIE700"] = json.dumps({"base": str(a.base), "base_sha256": hashlib.sha256(Path(a.base).read_bytes()).hexdigest(),
                                 "pred": str(a.pred), "pred_sha256": hashlib.sha256(Path(a.pred).read_bytes()).hexdigest(),
                                 "ctx": a.ctx, **diag})
    out = Path(a.out)
    if a.ctx_out:
        np.savez(a.ctx_out, perts=z["perts"], genes=z["genes"], contexts=np.array(list(per_ctx)),
                 **{f"bulk_{c}": blend(base, np.where(seen_mask[:, None], v, 0.0), a.w, nonpanel) for c, v in per_ctx.items()})
        print("ctx-out sha256", hashlib.sha256(Path(a.ctx_out).read_bytes()).hexdigest(), a.ctx_out, flush=True)
    if a.keep_env:
        np.savez(out, bulk=new_bulk, carriers=new_car, perts=z["perts"], genes=z["genes"], env=z["env"])
        Path(str(out) + ".pie700.json").write_text(env2["PIE700"])
    else:
        np.savez(out, bulk=new_bulk, carriers=new_car, perts=perts, genes=genes, env=np.array(json.dumps(env2)))
    np.savez(out.with_name(out.stem + "_pie_ctx.npz"), perts=perts, genes=genes, contexts=np.array(list(per_ctx)),
             **{f"bulk_{c}": v for c, v in per_ctx.items()}, shared=shared)
    print(json.dumps(diag), flush=True)
    print("sha256", hashlib.sha256(out.read_bytes()).hexdigest(), out, flush=True)


if __name__ == "__main__":
    main()
