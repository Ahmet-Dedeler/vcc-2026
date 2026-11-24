"""318 step 1 — recover the SHIPPED arm's per-target delivered magnitude profile, on the
REAL 300-target panel, from the two board-scored builds on disk.

Everything the mandate quotes about the distortion ("row norms 0.806-2.155, whole matrix
+6.0 %") is a FIXTURE number: 316's `test316.py` measured it on 313's synthetic fixture
(20 targets x 60 genes, 20 of them "panel"). This script measures the same quantity where
it actually acted -- the two 12.6 GB predictions that produced board points
`Q519I2pLWRCW7PAEx5ul` (shipped, 0.2177189) and `2K5ppxkowMSj2ugoZqYb` (parity, 0.2169733).

The DELIVERED log-effect of a block is what every DE member of the scorer reads:

    p(t,c)   = sum of the 400 emitted count vectors of target t in context c
    prof     = log1p(TS * p / p.sum())
    effect   = prof - b_ctrl(c)                 # b_ctrl from the challenge's own controls

and the per-target magnitude ratio is  m(t,c) = ||effect_shipped|| / ||effect_parity||.

Peak RSS is bounded by one block's nonzeros (~2M) plus two gene-length accumulators.
"""

from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
TS = 50_000.0

BUILDS = {
    "shipped": REPO / "experiments/312-normrel-arm/work/build/the-cat-sat-on-the-modem.h5ad",
    "parity": REPO / "experiments/316-parity-arm/work/build/a-duck-reviews-the-diff.h5ad",
}
CONTROLS = {c: REPO / f"data/challenge/context_{c}.h5ad" for c in ("A", "B", "C")}


def rss() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3) if sys.platform == "darwin" else r / (1024 ** 2)


def labels(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        cat = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                          for x in o["categories"][:]])
        return cat[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def control_profile(path: Path, n_genes: int) -> np.ndarray:
    """log1p(TS * cpm) of the pooled control block of one challenge context."""
    with h5py.File(path, "r") as h:
        X = h["X"]
        shape = tuple(int(x) for x in X.attrs["shape"])
        indptr = X["indptr"][:]
        tot = np.zeros(n_genes, dtype=np.float64)
        step = 4000
        for s in range(0, shape[0], step):
            e = min(s + step, shape[0])
            lo, hi = int(indptr[s]), int(indptr[e])
            tot += np.bincount(X["indices"][lo:hi].astype(np.int64),
                               weights=X["data"][lo:hi].astype(np.float64),
                               minlength=n_genes)
    return np.log1p(TS * tot / tot.sum())


def block_effects(path: Path, b_ctrl: dict[str, np.ndarray]):
    """Per (context, target) delivered log-effect row norms and the rows themselves are too
    big to keep, so accumulate norms + cross terms streaming block by block."""
    with h5py.File(path, "r") as h:
        X = h["X"]
        n_cells, n_genes = (int(x) for x in X.attrs["shape"])
        indptr = X["indptr"][:]
        ctx = labels(h, "obs/context")
        tgt = labels(h, "obs/target_gene")
        key = np.char.add(np.char.add(ctx, "|"), tgt)
        bounds = np.flatnonzero(np.r_[True, key[1:] != key[:-1], True])
        out = {}
        data, indices = X["data"], X["indices"]
        for bi in range(bounds.size - 1):
            s, e = int(bounds[bi]), int(bounds[bi + 1])
            lo, hi = int(indptr[s]), int(indptr[e])
            p = np.bincount(indices[lo:hi].astype(np.int64),
                            weights=data[lo:hi].astype(np.float64), minlength=n_genes)
            eff = np.log1p(TS * p / p.sum()) - b_ctrl[ctx[s]]
            out[(ctx[s], tgt[s])] = eff.astype(np.float32)
    return out, n_genes


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    genes = None
    with h5py.File(BUILDS["parity"], "r") as h:
        genes = labels(h, "var/_index")
    n_genes = genes.size
    panel_targets = sorted({t for t in
                            np.loadtxt(REPO / "data/challenge/pert_counts.csv", dtype=str,
                                       delimiter=",", skiprows=1, usecols=0)})
    panel_pos = np.asarray([i for i, g in enumerate(genes) if g in set(panel_targets)],
                           dtype=np.int64)
    non_panel = np.setdiff1d(np.arange(n_genes), panel_pos)
    print(f"[318] {n_genes} genes; {len(panel_targets)} panel targets, "
          f"{panel_pos.size} of them on the gene axis", flush=True)

    b = {c: control_profile(p, n_genes) for c, p in CONTROLS.items()}
    print(f"[318] control profiles built ({time.time() - t0:.0f}s, RSS {rss():.2f} GiB)",
          flush=True)

    eff = {}
    for name, path in BUILDS.items():
        eff[name], _ = block_effects(path, b)
        print(f"[318] {name}: {len(eff[name])} blocks ({time.time() - t0:.0f}s, "
              f"RSS {rss():.2f} GiB)", flush=True)

    rows = []
    for kk in sorted(eff["parity"]):
        a = eff["shipped"][kk].astype(np.float64)
        p = eff["parity"][kk].astype(np.float64)
        na, np_ = float(np.linalg.norm(a)), float(np.linalg.norm(p))
        naw = float(np.linalg.norm(a[non_panel]))
        npw = float(np.linalg.norm(p[non_panel]))
        cos = float(a @ p / max(na * np_, 1e-300))
        rows.append({"context": kk[0], "target": kk[1],
                     "norm_shipped": na, "norm_parity": np_,
                     "m_full": na / np_, "m_nonpanel": naw / npw, "row_cosine": cos})

    m = np.array([r["m_full"] for r in rows])
    whole = float(np.sqrt(sum(r["norm_shipped"] ** 2 for r in rows)
                          / sum(r["norm_parity"] ** 2 for r in rows)))
    summary = {
        "n_blocks": len(rows),
        "m_full": {"median": float(np.median(m)), "mean": float(m.mean()),
                   "min": float(m.min()), "max": float(m.max()),
                   "sd_log": float(np.std(np.log(m))),
                   "q05": float(np.quantile(m, 0.05)), "q95": float(np.quantile(m, 0.95)),
                   "frac_gt1": float((m > 1).mean()),
                   "frac_within_1e6": float((np.abs(m - 1) < 1e-6).mean())},
        "whole_matrix_norm_ratio": whole,
        "row_cosine": {"median": float(np.median([r["row_cosine"] for r in rows])),
                       "min": float(min(r["row_cosine"] for r in rows))},
        "peak_rss_gib": rss(),
        "seconds": round(time.time() - t0, 1),
    }
    (OUT / "profile318.json").write_text(json.dumps(
        {"summary": summary, "rows": rows}, indent=1))
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
