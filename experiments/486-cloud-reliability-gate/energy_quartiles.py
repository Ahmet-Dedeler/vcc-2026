"""Post-hoc: mean per-target cosine by quartile of truth energy.

Reads work/results/per_target.npz (written by evaluate.py) and prints the
quartile table. Not part of the pre-registered comparison.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

LINES = ("k562", "rpe1", "hepg2", "jurkat")
HERE = Path(__file__).resolve().parent
NPZ = HERE / "work" / "results" / "per_target.npz"
OUT = HERE / "work" / "results" / "energy_quartiles.csv"


def _corr(a: np.ndarray, b: np.ndarray) -> float:
    a = a - a.mean()
    b = b - b.mean()
    return float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))


def main() -> None:
    z = np.load(NPZ)
    rows = ["held_out,quartile,n,equal,fixed,learned,fixed_minus_equal,learned_minus_equal"]
    for line in LINES:
        energy = z[f"{line}__equal__energy"]
        equal = z[f"{line}__equal__cos"]
        fixed = z[f"{line}__fixed__cos"]
        learned = z[f"{line}__learned__cos"]
        edges = np.quantile(energy, [0.25, 0.5, 0.75])
        bin_id = np.digitize(energy, edges)
        for b, name in enumerate(("Q1", "Q2", "Q3", "Q4")):
            m = bin_id == b
            ge = float(fixed[m].mean() - equal[m].mean())
            le = float(learned[m].mean() - equal[m].mean())
            rows.append(
                f"{line},{name},{int(m.sum())},{equal[m].mean():.6f},"
                f"{fixed[m].mean():.6f},{learned[m].mean():.6f},{ge:.6f},{le:.6f}"
            )
        print(
            f"{line} corr(fixed-equal, energy)={_corr(fixed - equal, energy):+.4f} "
            f"corr(learned-equal, energy)={_corr(learned - equal, energy):+.4f}"
        )
    OUT.write_text("\n".join(rows) + "\n")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
