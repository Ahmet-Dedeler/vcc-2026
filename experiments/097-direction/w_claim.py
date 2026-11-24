"""The one assumption the `sqrt(W)` framing leaves unchecked: does an arm move `claim`?

`MSE_min - 1 = (Kc - 2*sqrt(claim*W))/sumD`, so the gain is `sqrt(claim*W)` and an arm that
raised `W` while cutting `claim` would be worth less than its `W` ratio says.  096 (b) measured
`claim` flat under AMPLITUDE changes (0.01% over s in 0.15..0.85), but a mixture change moves
the profile itself, which is a different perturbation.

`claim` needs the prediction's own delete-1 jackknife `Chat_p`, which needs real emitted cells
-- the exact-total repair kills emission noise in the MEAN but `Chat_p` is precisely the
cell-to-cell variance, so it cannot be reconstructed.  At ~0.2 s per group this is cheap
anyway: emit, jackknife, done.
"""
from __future__ import annotations

import sys

import numpy as np
import scipy.sparse as sp
from cell_eval2.metrics.delta import PRED_TRACE_CAP_K
from cell_eval2.moments import jackknife_correction
from w_arms import key_normrel, soft
from wkit import MIN_CPM, N_EMIT, SPREAD, TS, Ctx, uniform

sys.path.insert(0, "~/Code/virtual-cell-challenge-2026/src")
from vcc2026.donor_transfer import emit_counts  # noqa: E402


def chat(c, E, seed: int = 0) -> np.ndarray:
    """`Chat_p` for every scored row, through the production emitter."""
    rng = np.random.default_rng(seed)
    P = c.profiles(E)
    data, idx, indptr, codes = [], [], [0], []
    for i in range(len(c.perts)):
        pick = rng.choice(len(c.libs), N_EMIT, replace=len(c.libs) < N_EMIT)
        d, j, ip = emit_counts(P[i], c.libs[pick], 0.0, rng, mode="quasi",
                               spread=np.full(c.genes.size, SPREAD))
        data.append(d)
        idx.append(j)
        indptr.extend(list(ip[1:] + indptr[-1]))
        codes += [i] * N_EMIT
    M = sp.csr_matrix((np.concatenate(data).astype(np.float32),
                       np.concatenate(idx).astype(np.int32),
                       np.asarray(indptr, np.int64)),
                      shape=(len(codes), c.genes.size))
    return jackknife_correction(M, np.asarray(codes), len(c.perts), TS, chunk=64)


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else "h1-test"
    c = Ctx(name)
    Cp = c.t["Cp"]                       # reference's own jackknife, per scored row
    w = c.w
    b = c.run(uniform)
    keys = key_normrel(c)
    nb = np.linalg.norm(b["E"], axis=1)

    Edir = c.effects(soft(c, keys, -1.5))
    Edir = Edir * (nb / np.maximum(np.linalg.norm(Edir, axis=1), 1e-30))[:, None]
    cases = {"base uniform": b["E"], "normrel^-1.5 dir-only": Edir}

    print(f"########## {name}: does the arm move `claim`?")
    print(f"  scorer's shipped claim (from the cached prediction) = {c.claim:.8g}")
    print(f"  {'arm':<26}{'mean Chat':>12}{'claim':>12}{'claim/base':>12}"
          f"{'W':>12}{'sqrt(cl*W) ratio':>18}")
    ref = None
    for lbl, E in cases.items():
        Ch = chat(c, E)
        cl = float((np.minimum(Ch, PRED_TRACE_CAP_K * Cp) * w).sum())
        W = c.XAW(c.dh_exact(E))[2]
        if ref is None:
            ref = (cl, W)
        print(f"  {lbl:<26}{Ch.mean():>12.5g}{cl:>12.7g}{cl/ref[0]:>12.5f}"
              f"{W:>12.6g}{np.sqrt(cl*W/(ref[0]*ref[1])):>18.4f}")


if __name__ == "__main__":
    main()
