"""195 — known-answer checks for the reduced-coordinate fit.

Three things can silently be wrong here and all three would produce plausible numbers:
the reduced algebra could not equal the full-space cosine, the analytic gradient could not
equal the loss's derivative, and the M = 0 row could not equal the incumbent.  Each is
checked against an independent full-space computation.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fitlib import Adapter, BatchStats, apply_adapter, batch_forward  # noqa: E402


def full_cos(Z, T, Ek, M):
    Zp = apply_adapter(Z, Ek, M)
    Zn = Zp / np.linalg.norm(Zp, axis=1, keepdims=True)
    Tn = T / np.linalg.norm(T, axis=1, keepdims=True)
    return Zn @ Tn.T


def full_loss(Z, T, Ek, M, tau):
    cos = full_cos(Z, T, Ek, M)
    s = cos / tau
    s = s - s.max(1, keepdims=True)
    p = np.exp(s) / np.exp(s).sum(1, keepdims=True)
    return float(-np.log(np.diag(p)).mean())


def main() -> None:
    rng = np.random.default_rng(7)
    G, B, k, r, tau = 200, 24, 12, 3, 0.05
    Q = rng.normal(size=(B, G))
    T = rng.normal(size=(B, G))
    Ek = np.linalg.qr(rng.normal(size=(G, k)))[0]
    bs = BatchStats(Q, T, Ek)

    # 1. reduced cosine == full-space cosine, for a non-trivial M
    M = rng.normal(0, 0.3, (k, k))
    Gm = bs.Zp + bs.Zp @ M
    num = bs.R + Gm @ bs.Tb.T
    q = np.sqrt(bs.c + np.einsum("ij,ij->i", Gm, Gm))
    cos_red = num / (q[:, None] * bs.tn[None, :])
    err = np.abs(cos_red - full_cos(Q, T, Ek, M)).max()
    print(f"reduced vs full cosine   max |err| {err:.3e}")
    assert err < 1e-10, err

    # 2. analytic dM == central finite difference of the loss
    _, _, dM = batch_forward(bs, M, tau)
    h, num_g = 1e-6, np.zeros_like(M)
    for i in range(k):
        for j in range(k):
            Mp, Mm = M.copy(), M.copy()
            Mp[i, j] += h
            Mm[i, j] -= h
            num_g[i, j] = (full_loss(Q, T, Ek, Mp, tau) - full_loss(Q, T, Ek, Mm, tau)) / (2 * h)
    rel = np.abs(dM - num_g).max() / np.abs(num_g).max()
    print(f"analytic vs numeric dM   max rel err {rel:.3e}")
    assert rel < 1e-5, rel

    # 3. M = 0 is exactly the un-adapted operator, and Adam at init is a no-op on it
    ad = Adapter(k, r, diag=True, seed=0)
    assert np.abs(ad.M()).max() == 0.0 or ad.rank  # V starts at zero => U V^T == 0
    assert np.abs(ad.M()).max() == 0.0, ad.M().max()
    assert np.abs(apply_adapter(Q, Ek, ad.M()) - Q).max() == 0.0
    print("M=0 identity row         exact")

    # 4. in-batch pds at M=0 equals the direct rank computation
    _, pds, _ = batch_forward(bs, ad.M(), tau)
    cos = full_cos(Q, T, Ek, None)
    sd = np.diag(cos)[:, None]
    ref = float((((cos < sd).sum(1) + 0.5 * ((cos == sd).sum(1) - 1)) / (B - 1)).mean())
    print(f"in-batch pds             {pds:.12f} vs {ref:.12f}")
    assert abs(pds - ref) < 1e-12
    print("\nALL CHECKS PASS")


if __name__ == "__main__":
    main()
