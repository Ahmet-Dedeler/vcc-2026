"""194 Phase 1 pre-pass — per-cell METADATA only, no expression.

Parquet is columnar, so reading `gene_target`/`guide_target`/`sample`/`total_counts`/
`pct_counts_mt`/`n_genes_by_counts`/`pass_guide_filter` skips the two huge list columns
entirely.  That gives us, cheaply and before any expression is touched:

  * how many distinct guide CONSTRUCTS each target actually has  (candidate 2 in the
    mandate: guide-level aggregation.  If a donor screen is one construct per gene the
    whole candidate is structurally unavailable and must be reported as such);
  * the per-(target, batch) cell counts needed to build the batch-composition matching
    weights (candidate 4) BEFORE the weighted accumulation pass;
  * the depth / mito / gene-count distributions the QC thresholds (candidate 5) are set
    from, and the depth deciles the depth-matched permutation falsifier needs;
  * the true `pass_guide_filter` rate, which the incumbent Orion builder relies on.

Writes one npz per line into work/prepass-<line>.npz.
"""
from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

CONTROL_PATTERNS = (r"^non.?targeting", r"^ntc", r"^control", r"^safe.?harbou?r", r"^neg")

SMALL_COLS = ["gene_target", "guide_target", "sample", "total_counts",
              "pct_counts_mt", "n_genes_by_counts", "pass_guide_filter"]


def is_control(name: str) -> bool:
    n = str(name).lower()
    return any(re.match(p, n) for p in CONTROL_PATTERNS)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--line", default="HCT116")
    ap.add_argument("--data", default=str(ROOT / "data/external/orion/data"))
    args = ap.parse_args()

    files = sorted(Path(args.data).glob(f"{args.line}_Batch*.parquet"),
                   key=lambda p: int(re.search(r"Batch(\d+)", p.name)[1]))
    if not files:
        raise SystemExit(f"no {args.line} parquet in {args.data}")

    tgt_all: list[np.ndarray] = []
    guide_all: list[np.ndarray] = []
    batch_all: list[np.ndarray] = []
    tot_all: list[np.ndarray] = []
    mt_all: list[np.ndarray] = []
    ng_all: list[np.ndarray] = []
    pass_all: list[np.ndarray] = []

    t0 = time.time()
    for i, f in enumerate(files, 1):
        t = pq.read_table(f, columns=SMALL_COLS).to_pandas()
        tgt_all.append(t["gene_target"].to_numpy().astype("U"))
        guide_all.append(t["guide_target"].to_numpy().astype("U"))
        batch_all.append(np.full(len(t), i - 1, dtype=np.int32))
        tot_all.append(t["total_counts"].to_numpy().astype(np.float32))
        mt_all.append(t["pct_counts_mt"].to_numpy().astype(np.float32))
        ng_all.append(t["n_genes_by_counts"].to_numpy().astype(np.int32))
        pass_all.append(t["pass_guide_filter"].to_numpy().astype(np.int8))
        if i % 25 == 0 or i == len(files):
            print(f"  {i}/{len(files)}  {time.time() - t0:.0f}s", flush=True)

    tgt = np.concatenate(tgt_all)
    guide = np.concatenate(guide_all)
    batch = np.concatenate(batch_all)
    tot = np.concatenate(tot_all)
    mt = np.concatenate(mt_all)
    ng = np.concatenate(ng_all)
    passf = np.concatenate(pass_all)
    del tgt_all, guide_all, batch_all, tot_all, mt_all, ng_all, pass_all

    n = len(tgt)
    ctrl_names = {x for x in np.unique(tgt) if is_control(x)}
    ctrl = np.isin(tgt, list(ctrl_names))

    # guides per target, on the cells the incumbent keeps
    keep = passf == 1
    uniq_t = np.unique(tgt[keep & ~ctrl])
    order = np.argsort(tgt[keep & ~ctrl], kind="stable")
    tt = tgt[keep & ~ctrl][order]
    gg = guide[keep & ~ctrl][order]
    bounds = np.searchsorted(tt, uniq_t, side="left"), np.searchsorted(tt, uniq_t, side="right")
    gpt = np.array([len(np.unique(gg[a:b])) for a, b in zip(*bounds)], dtype=np.int32)

    summary = dict(
        line=args.line, n_batches=len(files), n_cells=int(n),
        pass_guide_filter_rate=float(keep.mean()),
        n_control_cells_kept=int((keep & ctrl).sum()),
        control_labels=sorted(ctrl_names),
        n_targets_kept=int(len(uniq_t)),
        guides_per_target=dict(
            mean=float(gpt.mean()), median=float(np.median(gpt)),
            frac_single=float((gpt == 1).mean()), max=int(gpt.max()),
            hist={str(k): int(v) for k, v in zip(*np.unique(np.clip(gpt, 0, 6),
                                                            return_counts=True))}),
        total_counts=dict(p1=float(np.percentile(tot[keep], 1)),
                          p10=float(np.percentile(tot[keep], 10)),
                          median=float(np.median(tot[keep])),
                          p90=float(np.percentile(tot[keep], 90))),
        pct_mt=dict(median=float(np.median(mt[keep])), p90=float(np.percentile(mt[keep], 90)),
                    p99=float(np.percentile(mt[keep], 99))),
        n_genes=dict(p1=float(np.percentile(ng[keep], 1)), median=float(np.median(ng[keep]))),
        control_vs_pert_depth=dict(
            control_median=float(np.median(tot[keep & ctrl])),
            pert_median=float(np.median(tot[keep & ~ctrl]))),
        elapsed_s=round(time.time() - t0, 1),
    )
    print(json.dumps(summary, indent=2))

    out = WORK / f"prepass-{args.line}.npz"
    np.savez_compressed(
        out, target=tgt, guide=guide, batch=batch, total_counts=tot, pct_mt=mt,
        n_genes=ng, pass_filter=passf, is_control=ctrl,
        summary=json.dumps(summary))
    (WORK / f"prepass-{args.line}.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
