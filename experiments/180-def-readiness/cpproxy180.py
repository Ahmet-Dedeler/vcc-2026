"""180 - reproduce 096's `C_p` NTC proxy from the bundle alone, and pin the selection rule.

096's MSE-MECHANISM states the proxy as "the delete-1 jackknife of a single ntc construct's
400 cells - same cell count, same 20k median depth, same construct-level structure as a
reference perturbation group", 16 constructs per context, reading
**14.544 / 14.772 / 16.372** on A / B / C (sd 0.39 / 0.53 / 0.73). The cached
`work/cp-proxy-ntc.npz` holds those 16 values per context but NOT the selection rule, and a
D/E/F builder has to recompute the proxy from a bundle it has never seen.

So: compute the delete-1 jackknife for EVERY one of the 46 ntc constructs in each context,
then ask which 16-subset rule reproduces the cached vector. The known-answer row is the
cached npz itself -- if no rule reproduces it the proxy is not portable and the knee cannot
be computed on a new bundle, which is a finding and not a detail.

    uv run python cpproxy180.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
TS = 50_000.0
CTX = ("A", "B", "C")


def per_construct_jackknife(context: str) -> tuple[list[str], np.ndarray]:
    """`C_p` for every ntc construct in one context, in the file's own category order."""
    import anndata as ad
    from scipy.sparse import csr_matrix
    from cell_eval2.moments import jackknife_correction

    handle = ad.read_h5ad(REPO / f"data/challenge/context_{context}.h5ad")
    matrix = handle.X.tocsr()
    ids = handle.obs["ntc_id"].astype(str).to_numpy()
    order = list(dict.fromkeys(ids))          # first-appearance order, before any sorting
    values = []
    for name in order:
        rows = np.flatnonzero(ids == name)
        sub = csr_matrix(matrix[rows])
        values.append(
            jackknife_correction(sub, np.zeros(rows.size, np.intp), 1, TS, chunk=64)[0]
        )
    del handle, matrix
    return order, np.asarray(values, dtype=np.float64)


def main() -> None:
    WORK.mkdir(exist_ok=True)
    cached = dict(np.load(REPO / "experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())
    out: dict[str, object] = {}
    for context in CTX:
        t0 = time.time()
        order, values = per_construct_jackknife(context)
        target = np.asarray(cached[context], dtype=np.float64)
        alpha = sorted(order)
        rules = {
            "appearance_first16": [order.index(n) for n in order[:16]],
            "sorted_first16": [order.index(n) for n in alpha[:16]],
            "natural_first16": [
                order.index(n)
                for n in sorted(order, key=lambda s: int(s.rsplit("-", 1)[1]))[:16]
            ],
            "rng0_choice16": list(
                np.sort(np.random.default_rng(0).choice(len(order), 16, replace=False))
            ),
        }
        hits = {}
        for name, idx in rules.items():
            drawn = np.sort(values[np.asarray(idx, dtype=int)])
            hits[name] = float(np.abs(drawn - np.sort(target)).max())
        print(
            f"[cp] {context}: 46-construct mean {values.mean():.6f} sd {values.std():.6f}; "
            f"cached-16 mean {target.mean():.6f} sd {target.std():.6f}; "
            f"rule max|diff| " + " ".join(f"{k}={v:.3e}" for k, v in hits.items())
            + f"  [{time.time() - t0:.0f}s]",
            flush=True,
        )
        # Does the cached 16 sit inside the 46 at all? (set membership, order-free)
        inside = sum(
            1 for t in target if np.min(np.abs(values - t)) < 1e-9
        )
        print(f"[cp] {context}: {inside}/16 cached values are exactly among the 46",
              flush=True)
        out[context] = {
            "order": order,
            "values": values.tolist(),
            "mean46": float(values.mean()),
            "sd46": float(values.std()),
            "cached_mean16": float(target.mean()),
            "cached_sd16": float(target.std()),
            "rule_max_abs_diff": hits,
            "cached_inside_46": inside,
        }
        (WORK / "cp_proxy_46.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({c: {k: v for k, v in out[c].items()
                          if k not in ("order", "values")} for c in CTX}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
