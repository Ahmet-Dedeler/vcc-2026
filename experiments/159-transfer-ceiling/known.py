"""159 §1 -- reproduce the brief's known-answer row before doing anything else.

The brief reports, on the two thinned Orion harnesses, in the scorer's own profile space
(b = log1p(TS * P / sum P), TS = 5e4, effect = b_target - b_ctrl), over the targets the two
surfaces share, with all panel target-gene coordinates dropped exactly as `pds_cosine` drops
them:

    same-target cross-line cosine   +0.1881  (median 0.1965, sd 0.1031)
    diff-target cross-line cosine   +0.1790
    target-specific excess          +0.0091
    within-HCT116  diff-target      +0.2476
    within-HEK293T diff-target      +0.4175
    shared / line-specific energy   60.1% / 39.9%
    context-independent ceiling     cos 0.7700

If this script cannot reproduce those numbers, STOP: one of the two measurements has a bug.
It reads 156's frozen exact-space truth (`experiments/156-delivered-geometry/work/truth-*.npz`)
read-only; nothing here is rebuilt.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
T156 = os.path.join(REPO, "experiments/156-delivered-geometry/work/truth-%s.npz")


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def offdiag_mean(S):
    n = S.shape[0]
    return float(S[~np.eye(n, dtype=bool)].mean())


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    A = np.load(T156 % "hct116", allow_pickle=True)
    B = np.load(T156 % "hek293t", allow_pickle=True)
    ga = np.asarray([str(x) for x in A["genes"]])
    gb = np.asarray([str(x) for x in B["genes"]])
    assert np.array_equal(ga, gb), "gene axes differ"
    ta = np.asarray([str(x) for x in A["targets"]])
    tb = np.asarray([str(x) for x in B["targets"]])

    shared = np.asarray(sorted(set(ta) & set(tb)))
    ia = np.asarray([int(np.flatnonzero(ta == s)[0]) for s in shared])
    ib = np.asarray([int(np.flatnonzero(tb == s)[0]) for s in shared])
    EA = A["effect"].astype(np.float64)[ia]
    EB = B["effect"].astype(np.float64)[ib]

    out = {"n_shared_targets": int(len(shared)),
           "n_targets_hct116": int(len(ta)), "n_targets_hek293t": int(len(tb)),
           "cells_per_target_hct116": int(np.median(A["cells"])),
           "cells_per_target_hek293t": int(np.median(B["cells"])),
           "ctrl_cells_hct116": int(A["ctrl_cells"]), "ctrl_cells_hek293t": int(B["ctrl_cells"])}

    # `pds_cosine` drops EVERY panel target's own coordinate, not just the row's own.
    # Three candidate drop sets, so the reproduction cannot be ambiguous about which is meant.
    gpos = {g: i for i, g in enumerate(ga)}
    drops = {
        "shared120": {gpos[s] for s in shared if s in gpos},
        "union_panel": {gpos[s] for s in set(ta) | set(tb) if s in gpos},
        "none": set(),
    }

    for label, drop in drops.items():
        mask = np.asarray([i for i in range(len(ga)) if i not in drop])
        a, b = unit(EA[:, mask]), unit(EB[:, mask])
        S = a @ b.T
        same = np.diag(S)
        diff = offdiag_mean(S)
        wa = offdiag_mean(a @ a.T)
        wb = offdiag_mean(b @ b.T)

        # energy split over the same targets: shared = the two-line mean, line-specific = the
        # half-difference.  ||x||^2 = ||m||^2 + ||d||^2 with m = (A+B)/2, d = (A-B)/2.
        m = 0.5 * (EA[:, mask] + EB[:, mask])
        dd = 0.5 * (EA[:, mask] - EB[:, mask])
        em, ed = float((m * m).sum()), float((dd * dd).sum())
        # cosine of the two-line mean to each line, per target, macro-averaged
        cm = np.concatenate([(unit(m) * unit(EA[:, mask])).sum(1),
                             (unit(m) * unit(EB[:, mask])).sum(1)])

        out[label] = {
            "n_scored_genes": int(len(mask)), "n_dropped": int(len(drop)),
            "same_target_mean": float(same.mean()), "same_target_median": float(np.median(same)),
            "same_target_sd": float(same.std(ddof=1)),
            "diff_target_mean": diff, "excess": float(same.mean() - diff),
            "within_hct116_diff": wa, "within_hek293t_diff": wb,
            "energy_shared_frac": em / (em + ed), "energy_linespec_frac": ed / (em + ed),
            "ceiling_cos_mean_to_line": float(cm.mean()),
        }
        r = out[label]
        print(f"[{label:12}] genes {r['n_scored_genes']:6d} (dropped {r['n_dropped']:3d})  "
              f"same {r['same_target_mean']:+.4f} (med {r['same_target_median']:+.4f} "
              f"sd {r['same_target_sd']:.4f})  diff {r['diff_target_mean']:+.4f}  "
              f"excess {r['excess']:+.4f}")
        print(f"{'':14} within-HCT {r['within_hct116_diff']:+.4f}  "
              f"within-HEK {r['within_hek293t_diff']:+.4f}  "
              f"energy shared {r['energy_shared_frac']:.3f}/{r['energy_linespec_frac']:.3f}  "
              f"ceiling cos {r['ceiling_cos_mean_to_line']:.4f}")

    with open(os.path.join(WORK, "known.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(f"\nwrote {os.path.join(WORK, 'known.json')}")


if __name__ == "__main__":
    main()
