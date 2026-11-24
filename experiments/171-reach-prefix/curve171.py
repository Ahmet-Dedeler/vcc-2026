"""171 — draw the P(k) curve and the head composition.  The mandate asks for the CURVE."""
from __future__ import annotations

import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, f"{ROOT}/experiments/171-reach-prefix")
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from diag171 import Walk  # noqa: E402

W_DIR = f"{ROOT}/experiments/171-reach-prefix/work"
ARMS = [("a15", "a15 as shipped"), ("armC_a8_flat", "a8 flat 0.25/phi4"),
        ("kt20000", "keep_top >= gated"), ("flat100_phi4", "flat 1.0/phi4")]
KS = np.arange(1, 101)


def main():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for nm, lab in ARMS:
        try:
            W = Walk(f"{W_DIR}/de171.{nm}.npz")
        except FileNotFoundError:
            continue
        rows = W.rows()
        pk, tf = [], []
        for k in KS:
            v = [r["purity"][k - 1] for r in rows if r["n_adj"] >= k]
            pk.append(np.mean(v) if v else np.nan)
            t = [r["tail_walked"][k - 1] for r in rows if r["n_adj"] >= k]
            tf.append(np.mean(t) if t else np.nan)
        ax[0].plot(KS, pk, label=lab, lw=1.6)
        ax[1].plot(KS, tf, label=lab, lw=1.6)
    ax[0].axhline(0.9, color="k", ls="--", lw=1, label="$P_0=0.9$")
    ax[0].axhline(0.5, color="grey", ls=":", lw=1, label="chance")
    ax[0].set_xscale("log")
    ax[0].set_xlabel("depth $k$ (adjudicable genes walked)")
    ax[0].set_ylabel("mean cumulative purity $P(k)$")
    ax[0].set_title("h1-test: $P(k)$ dies before $k=2$")
    ax[0].set_ylim(0.3, 1.02)
    ax[0].legend(fontsize=7, loc="upper right")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("position in the walk")
    ax[1].set_ylabel("fraction from the zero-effect TAIL block")
    ax[1].set_title("who occupies the head of the walk")
    ax[1].set_ylim(0, 1.02)
    ax[1].legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    out = f"{ROOT}/experiments/171-reach-prefix/pk-curve.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
