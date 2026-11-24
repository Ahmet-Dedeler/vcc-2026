#!/usr/bin/env python3
"""700 — write PIE's per-context predictions as a VCC_SAMELINE artifact, so the existing same-line plumbing carries them.

The same-line mode (builder 451/517; `scripts/def_sameline.py` writes the real ones) does, per matched context and covered
target: design row AND bulk row -> unit(ours) + w * unit(bank row) over non-panel genes at our row's own norm (direction
moves, magnitude kept, panel columns untouched); with `pt`, the target's DE block prior takes the bank row's sign wherever
`signs` is non-zero; without `noprog`, the context's default block prior becomes the bank programme. Here the "bank" of
every context is PIE's prediction for that context (context-specific: PIE reads each line's control mean and Cellosaurus
text), so the knob gives a per-context PIE blend with no builder change.

Rows: PIE lfc_pred mapped into the scorer's effect space with the context's control pseudobulk (see make_pie_dump700.py).
Signs (`pt` only): sign(lfc_pred) where p_de >= --sign-pde (default 0.5, PIE's own DE threshold) and the gene is
predicted, else 0. Programme: required by the format as +-1; written as -1 everywhere with programme_ok False, and the arm
must run with `noprog` (CLEAN4 stays the default block prior).
usage: .venv/bin/python make_pie_sameline700.py --pred <predictions.parquet> --dataset vcc26abc --ctx A=jurkat
         --ctx B=CVCL_0030 --ctx C=CVCL_1108 --bundle data/challenge [--sign-pde 0.5] --out <npz>
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pyarrow.parquet as pq

TS = 5e4


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--ctx", action="append", required=True)
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--sign-pde", type=float, default=0.5)
    ap.add_argument("--seen", default=None, help="JSON list of PIE train-split perturbations: cover only these targets")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    seen = set(json.load(open(a.seen))) if a.seen else None
    bundle = Path(a.bundle)
    genes = [r["gene_name"] for r in csv.DictReader(open(bundle / "gene_names.csv"))]
    perts = [r["target_gene"] for r in csv.DictReader(open(bundle / "pert_counts.csv")) if r["target_gene"] != "non-targeting"]
    t = pq.read_table(a.pred)
    pg = json.loads(t.schema.metadata[b"pie"])["genes"][a.dataset]
    gpos = {g: i for i, g in enumerate(pg)}
    have = np.array([g in gpos for g in genes])
    gix = np.array([gpos.get(g, 0) for g in genes])
    df = t.to_pandas()
    df = df[df["dataset"] == a.dataset]
    G = len(genes)
    arrays = {"genes": np.array(genes), "perts": np.array(perts)}
    ctxs = []
    summary = {}
    for spec in a.ctx:
        bctx, pctx = spec.split("=")
        h = ad.read_h5ad(bundle / f"context_{bctx}.h5ad")
        assert list(map(str, h.var_names)) == genes
        P = np.asarray(h.X.sum(0)).ravel().astype(np.float64)
        p = P / P.sum()
        sub = df[df["context"] == pctx].set_index("perturbation")
        tg, rows, signs = [], [], []
        for i, x in enumerate(perts):
            if x not in sub.index or (seen is not None and x not in seen):
                continue
            lf = np.where(have, np.asarray(sub.loc[x, "lfc_pred"], np.float64)[gix], 0.0)
            pde = np.where(have, np.asarray(sub.loc[x, "p_de"], np.float64)[gix], 0.0)
            q = p * np.exp2(lf)
            q = q / q.sum()
            rows.append(np.log1p(TS * q) - np.log1p(TS * p))
            signs.append(np.where(have & (pde >= a.sign_pde) & (lf != 0), np.sign(lf), 0.0))
            tg.append(i)
        arrays[f"targets_{bctx}"] = np.array(tg, dtype=np.int64)
        arrays[f"rows_{bctx}"] = np.stack(rows)
        arrays[f"signs_{bctx}"] = np.stack(signs)
        arrays[f"programme_{bctx}"] = -np.ones(G)
        arrays[f"programme_ok_{bctx}"] = np.array(False)
        ctxs.append(bctx)
        S = arrays[f"signs_{bctx}"]
        summary[bctx] = {"pie_context": pctx, "targets": len(tg), "signed_per_target_median": float(np.median((S != 0).sum(1))),
                         "down_frac_of_signed": float((S < 0).sum() / max((S != 0).sum(), 1))}
    arrays["contexts"] = np.array(ctxs)
    out = Path(a.out)
    np.savez(out, **arrays)
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    meta = {"seen_rule": a.seen, "pred": a.pred, "pred_sha256": hashlib.sha256(Path(a.pred).read_bytes()).hexdigest(), "sign_pde": a.sign_pde,
            "summary": summary, "sha256": sha}
    out.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    print(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
