"""Verify that pds_cosine is scale-invariant while expr_mse is quadratic in magnitude.

If true, the predicted-effect magnitude alpha is a free parameter: it moves the mse
metric without touching pds. See AGENTS.md section 4.1.

Run:  uv run python experiments/001-magnitude-calibration/calibration_math.py
"""

import anndata as ad
import numpy as np

TS = 5e4
# Reference anchors for context A/B/C, from the official vcc2026 metric spec.
MSE_B, MSE_R = 0.989, 0.036


def scaled_mse(raw: float) -> float:
    return min(1.0, max(0.0, (raw - MSE_B) / (MSE_R - MSE_B)))


def cosine_pds_rank(dhat, deltas, true_idx):
    """PDS_p as the normalised rank of the true match under cosine distance."""
    dn = dhat / np.linalg.norm(dhat)
    D = deltas / np.linalg.norm(deltas, axis=1, keepdims=True)
    d = 1.0 - D @ dn
    n = len(deltas)
    k = np.sum(d < d[true_idx]) + 0.5 * (np.sum(d == d[true_idx]) - 1)
    return 1.0 - k / (n - 1)


def main():
    a = ad.read_h5ad("data/challenge/context_A.h5ad")
    P = np.asarray(a.X.sum(0)).ravel().astype(float)
    b_ctrl = np.log1p(TS * P / P.sum())
    G = len(b_ctrl)
    expressed = b_ctrl > 0
    rng = np.random.default_rng(0)

    # A synthetic 300-perturbation panel.
    n_pert = 300
    deltas = rng.normal(0, 0.25, (n_pert, G)) * expressed

    print("== 1. Is pds_cosine invariant to rescaling the predicted effect? ==")
    delta = deltas[0]
    noise = rng.normal(0, 1, G) * expressed
    noise -= noise.dot(delta) / delta.dot(delta) * delta
    noise *= np.linalg.norm(delta) / np.linalg.norm(noise)
    dhat = delta + np.sqrt(1 / 0.5**2 - 1) * noise  # cos = 0.5

    for alpha in (0.01, 0.1, 0.25, 1.0, 10.0, 100.0):
        print(f"   alpha={alpha:7.2f}  PDS_p = {cosine_pds_rank(alpha * dhat, deltas, 0):.6f}")

    print("\n== 2. How does raw mse move with alpha? ==")
    print(f"   {'cos':>5} {'mse@a=1':>9} {'a*':>6} {'mse@a*':>8} {'1-cos^2':>8} {'scaled':>7}")
    for cos_target in (0.2, 0.3, 0.5, 0.7, 0.9):
        noise = rng.normal(0, 1, G) * expressed
        noise -= noise.dot(delta) / delta.dot(delta) * delta
        noise *= np.linalg.norm(delta) / np.linalg.norm(noise)
        dhat = delta + np.sqrt(1 / cos_target**2 - 1) * noise
        cos = dhat.dot(delta) / np.linalg.norm(dhat) / np.linalg.norm(delta)

        D = delta.dot(delta)
        a_star = dhat.dot(delta) / dhat.dot(dhat)
        m1 = np.sum((dhat - delta) ** 2) / D
        ms = np.sum((a_star * dhat - delta) ** 2) / D
        print(f"   {cos:5.2f} {m1:9.3f} {a_star:6.3f} {ms:8.4f} {1 - cos**2:8.4f} {scaled_mse(ms):7.3f}")

    print("\n== 3. Where the leaderboard sits (raw mse -> scaled, contribution to overall) ==")
    for team, raw in [("wombat", 1.562), ("IsoCell", 1.974), ("Vivai", 3.481),
                      ("Sidechain", 3.608), ("Gia Khang", 6.535), ("p19", 6.965),
                      ("control-copy", 1.003)]:
        s = scaled_mse(raw)
        print(f"   {team:14} raw {raw:6.3f} -> scaled {s:6.3f}  (overall {s / 6:+.4f})")

    print("\n   min over alpha is always <= mse(alpha=0) ~= 1.003, so any raw mse > 1")
    print("   is provably not alpha-optimised. All top-10 teams have raw mse > 1.56.")


if __name__ == "__main__":
    main()
