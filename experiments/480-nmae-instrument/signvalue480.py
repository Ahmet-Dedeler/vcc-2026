"""480 step 2c: what is sign information worth, per gene class, in nmae units (H1, our shipped construction)?

For a class C, raise our |l|-weighted gate sign accuracy by +da: flip a random subset of our wrong-signed gate cells in
C to the right sign (iid within C, drawn until the |l|-weighted accuracy gain reaches da), keep our magnitudes, and
read nmae; then also re-fit a global MDE add (c on the base, in-sample), since better signs support larger magnitudes.
20 random draws per cell of the table; the spread across draws is reported. This is an oracle-flavoured price: it
assumes the new information is iid across the class (434f's caveat: magnitudes that follow confidence reach the same
nmae at a lower average accuracy).

    uv run python signvalue480.py    # -> work/signvalue480.json
"""
from __future__ import annotations

import json

import numpy as np

from common480 import HERE, nmae
from decomp480 import h1_inputs

DA = [0.02, 0.05, 0.10]
NDRAW = 20
C_GRID = np.linspace(0, 0.4, 41)


def main():
    h, T, cnt, genes, targets, Lc, qc, L0, cons = h1_inputs()
    G = np.isfinite(T)
    Tz = np.where(G, T, 0.0)
    s = np.sign(Lc)
    s = np.where(s == 0, np.sign(cons), s)
    s = np.where(s == 0, 1.0, s)
    add1 = (np.minimum(1 / np.sqrt(np.maximum(cnt, 1e-3)), 1.0) * (cnt >= 0.3))[None, :]
    base = L0
    ship = nmae(base + s * 0.1 * add1, T)
    print(f"shipped analytic {ship:.5f}")
    q1 = np.quantile(np.broadcast_to(cnt, T.shape)[G], 2 / 3)
    classes = {"all gate genes": G, "uncalled": G & (qc >= 0.05), "called": G & (qc < 0.05),
               "uncalled, control count < 2.6/cell": G & (qc >= 0.05) & (cnt[None, :] < q1)}
    rng = np.random.default_rng(0)
    out = {}
    for name, M in classes.items():
        w = np.abs(Tz)
        right = (np.sign(base + s * 0.1 * add1) == np.sign(Tz)) & M
        acc0 = (w * right).sum() / (w * M).sum()
        wrong_idx = np.flatnonzero((M & ~right).ravel())
        row = {"acc0": float(acc0), "mass_share": float((w * M).sum() / (w * G).sum())}
        for da in DA:
            vals, valc, cs = [], [], []
            for _ in range(NDRAW):
                perm = rng.permutation(wrong_idx)
                cum = np.cumsum(w.ravel()[perm]) / (w * M).sum()
                nflip = int(np.searchsorted(cum, da)) + 1
                sf = s.copy().ravel()
                bf = base.copy().ravel()
                idx = perm[:nflip]
                sf[idx] = -sf[idx]
                bf[idx] = -bf[idx]                      # flip our base's sign too (same magnitude)
                sf, bf = sf.reshape(s.shape), bf.reshape(base.shape)
                vals.append(nmae(bf + sf * 0.1 * add1, T) - ship)
                v = [nmae(bf + sf * c * add1, T) for c in C_GRID]
                valc.append(min(v) - ship)
                cs.append(float(C_GRID[int(np.argmin(v))]))
            row[f"+{da:.2f}"] = dict(d_nmae=float(np.mean(vals)), sd=float(np.std(vals)),
                                     d_nmae_refit_c=float(np.mean(valc)), c_refit=float(np.mean(cs)))
            print(f"{name:36s} acc {acc0:.3f} (+{da:.2f}) mass {row['mass_share']:.3f}: Δnmae {np.mean(vals):+.4f} "
                  f"(sd {np.std(vals):.4f}); with c refit {np.mean(valc):+.4f} (c {np.mean(cs):.2f})", flush=True)
        out[name] = row
    json.dump(out, open(HERE / "work/signvalue480.json", "w"), indent=1)


if __name__ == "__main__":
    main()
