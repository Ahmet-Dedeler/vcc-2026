"""194 — build the FULL-row `debias` tables into `data/donors/194-debias-<donor>/`.

`debias` is the one candidate that needs no re-accumulation: it is a closed-form
correction applied to the two profiles the shipped table is already the difference of.

    delta_shipped = b_T - b_C ,   b = log1p(TS * p_hat)

`log1p` is concave and `p_hat` is a finite-count estimate, so
`E[b] = log1p(TS p) - TS^2 p / (2 (1+TS p)^2 N)` to second order, with `N` the total
counts in the pool.  The control pool carries 10^2-10^3 times more counts than any one
target, so the bias does NOT cancel in the difference; what survives is a systematic,
target-INDEPENDENT offset on low-expressed genes in every row -- the shape of a shared
blob, which is the metric's own zero (003) and the thing AGENTS.md 6.4 warns about.

Adding the correction back needs only `N_T` (total counts per target) and `N_C`, both of
which come from the per-cell depths already computed by `accumulate.py --stage obs`.  So
these tables are exact full-row rebuilds, byte-compatible with the incumbent: same gene
axis, same target list, same dtype, same file layout.

    uv run python fulltables.py --donor all
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

import accumulate as A
from tables import _debias

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
TS = 5e4


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(8 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def total_counts(donor: str) -> tuple[dict[str, float], float]:
    """Total counts per target and in the control pool, on the donor's own gene axis."""
    cells = A.cell_table(donor)
    depth = np.load(WORK / f"obs-{donor}.npz")["depth"].astype(np.float64)
    keep = cells["keep"]
    ctrl = cells["is_ctrl"] & keep
    pert = (~cells["is_ctrl"]) & keep
    names, inv = np.unique(cells["target"][pert], return_inverse=True)
    tot = np.bincount(inv, weights=depth[pert], minlength=len(names))
    return {str(n): float(t) for n, t in zip(names, tot)}, float(depth[ctrl].sum())


def build(donor: str, k: float = 1.0) -> dict:
    src = ROOT / "data/donors" / donor
    out = ROOT / "data/donors" / f"194-debias-{donor}"
    delta = np.load(src / "delta.npy").astype(np.float64)
    b_ctrl = np.load(src / "b_ctrl.npy").astype(np.float64)
    perts = (src / "perts.txt").read_text().split("\n")
    genes = (src / "genes.txt").read_text().split("\n")
    assert delta.shape == (len(perts), len(genes))

    per_target, n_ctrl = total_counts(donor)
    missing = [p for p in perts if p not in per_target]
    assert not missing, f"{donor}: {len(missing)} targets have no cell counts, e.g. {missing[:3]}"
    n_t = np.array([per_target[p] for p in perts])

    b_t = delta + b_ctrl                      # the target profile the table differences
    b_t_c = _debias(b_t, n_t[:, None], k)
    b_c_c = _debias(b_ctrl, np.array(n_ctrl), k)
    new = (b_t_c - b_c_c).astype(np.float32)

    out.mkdir(parents=True, exist_ok=True)
    np.save(out / "delta.npy", new)
    np.save(out / "b_ctrl.npy", b_c_c.astype(np.float32))
    (out / "perts.txt").write_text("\n".join(perts))
    (out / "genes.txt").write_text("\n".join(genes))
    old = delta.astype(np.float32)
    n0 = np.linalg.norm(old, axis=1)
    n1 = np.linalg.norm(new, axis=1)
    rel = np.linalg.norm(new - old, axis=1) / np.maximum(n0, 1e-30)
    meta = json.loads((src / "meta.json").read_text())
    meta.update(
        variant="debias", experiment="194-donor-reestimation", debias_k=k,
        derived_from=dict(path=str(src.relative_to(ROOT)),
                          delta_sha256=sha256(src / "delta.npy"),
                          b_ctrl_sha256=sha256(src / "b_ctrl.npy")),
        total_counts=dict(control=n_ctrl, target_median=float(np.median(n_t)),
                          ratio_control_to_median_target=float(n_ctrl / np.median(n_t))),
        change=dict(rel_change_median=float(np.median(rel)),
                    rel_change_p90=float(np.percentile(rel, 90)),
                    amplification_median=float(np.median(n1 / np.maximum(n0, 1e-30)))),
        effect_norm=dict(min=float(n1.min()), median=float(np.median(n1)),
                         max=float(n1.max())))
    (out / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"{donor}: N_ctrl/N_target = {n_ctrl / np.median(n_t):.1f}x, "
          f"median relative change {np.median(rel):.4f}, "
          f"amplification {np.median(n1 / np.maximum(n0, 1e-30)):.5f} -> {out}")
    return meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donor", default="all")
    ap.add_argument("--k", type=float, default=1.0)
    args = ap.parse_args()
    todo = (["k562_gwps", "orion_hct116", "orion_hek293t"]
            if args.donor == "all" else [args.donor])
    metas = {d: build(d, args.k) for d in todo}
    (WORK / "fulltables.json").write_text(json.dumps(metas, indent=2))


if __name__ == "__main__":
    main()
