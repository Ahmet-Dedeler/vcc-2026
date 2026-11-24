"""480 step 3b: class-based magnitude rules, fitted on training targets and tested on held-out targets.

H1 (A/B/C power), our real base. For each oracle-free class variable and rule shape (c-only: base + s*c_C*MDE; m+c:
m_C*base + s*c_C*MDE), the per-class parameters are chosen on 4/5 of the scored targets and applied to the held-out
fifth (5-fold over targets, 50 random partitions). Per target, the held-out Δ vs the shipped rule (m 1, c 0.1) is
averaged over partitions; the interval is a target bootstrap (4,000 draws, unit = scored target) of those per-target
means. Panel-like subset (no essential-screen coverage) reported beside it.

Per-(grid point, target, class) contributions are precomputed once, so every fit is an argmin over arrays.

    uv run python cv480.py      # -> work/cv480.json
"""
from __future__ import annotations

import json

import numpy as np

from common480 import HERE, ROOT, nmae
from decomp480 import h1_inputs
from rules480 import C_GRID, M_GRID, class_vars

R, K, NB = 50, 5, 4000


def contrib_table(T, base, s, add1, cls, n, ms, cs):
    """tab[i_m, i_c, target, class] = sum_{g in class} |m*base + s*c*add - l| / D_p (NaN rows for unscored targets)."""
    G = np.isfinite(T)
    Tz = np.where(G, T, 0.0)
    D = np.abs(Tz).sum(1)
    D[D == 0] = np.nan
    onehot = [(G & (cls == c)) for c in range(n)]
    tab = np.zeros((len(ms), len(cs), T.shape[0], n))
    for i, m in enumerate(ms):
        mb = m * base
        for j, c in enumerate(cs):
            E = np.abs(mb + s * c * add1 - Tz)
            for k in range(n):
                tab[i, j, :, k] = np.where(onehot[k], E, 0.0).sum(1) / D
    return tab


def cv_rows(tab, scored, i0, j0, rng):
    """Held-out per-target Δ vs shipped (grid index i0, j0), averaged over R random K-fold partitions."""
    idx = np.flatnonzero(scored)
    nm, nc, _, n = tab.shape
    ship = tab[i0, j0].sum(1)                                    # per-target shipped nmae
    acc = np.zeros(tab.shape[2])
    cnt = np.zeros(tab.shape[2])
    picks = []
    for _ in range(R):
        perm = rng.permutation(idx)
        folds = np.array_split(perm, K)
        for f in folds:
            train = np.setdiff1d(idx, f)
            m = tab[:, :, train, :].mean(2)                       # grid x grid x class
            flat = m.reshape(nm * nc, n)
            best = flat.argmin(0)                                 # per class
            bi, bj = np.unravel_index(best, (nm, nc))
            held = sum(tab[bi[k], bj[k], f, k] for k in range(n))
            acc[f] += held - ship[f]
            cnt[f] += 1
            picks.append((bi, bj))
    d = np.full(tab.shape[2], np.nan)
    d[idx] = acc[idx] / cnt[idx]
    return d, picks


def boot(d, sub, rng):
    x = d[sub & np.isfinite(d)]
    bs = x[rng.integers(0, x.size, size=(NB, x.size))].mean(1)
    return float(x.mean()), float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))


def main():
    h, T, cnt, genes, targets, Lc, qc, L0, cons = h1_inputs()
    s = np.sign(Lc)
    s = np.where(s == 0, np.sign(cons), s)
    s = np.where(s == 0, 1.0, s)
    add1 = (np.minimum(1 / np.sqrt(np.maximum(cnt, 1e-3)), 1.0) * (cnt >= 0.3))[None, :]
    scored = np.isfinite(T).any(1)
    ess = set()
    for d in ["k562_essential", "rpe1", "nadig_jurkat"]:
        ess |= {l.strip() for l in open(ROOT / f"data/donors/{d}/perts.txt")}
    panel_like = np.array([t not in ess for t in targets])
    i0 = int(np.argmin(np.abs(M_GRID - 1.0)))
    j0 = int(np.argmin(np.abs(C_GRID - 0.1)))
    print(f"shipped analytic {nmae(L0 + s * 0.1 * add1, T):.5f}; grid shipped point m {M_GRID[i0]} c {C_GRID[j0]}")
    rng = np.random.default_rng(0)
    out = {}
    for name, (cls, n) in class_vars(T, Lc, qc, cnt, cons, targets, genes, s).items():
        row = {}
        for shape, ms in [("c", np.array([1.0])), ("m+c", M_GRID)]:
            tab = contrib_table(T, L0, s, add1, cls, n, ms, C_GRID)
            ii0 = 0 if shape == "c" else i0
            ship_check = np.nanmean(tab[ii0, j0].sum(1)[scored])
            assert abs(ship_check - nmae(L0 + s * 0.1 * add1, T)) < 1e-9
            d, picks = cv_rows(tab, scored, ii0, j0, rng)
            b_all = boot(d, scored, rng)
            b_pl = boot(d, scored & panel_like, rng)
            insample = tab.reshape(-1, T.shape[0], n)[:, scored, :].mean(1).min(0).sum() - ship_check
            row[shape] = dict(insample=float(insample), heldout=b_all, heldout_panel_like=b_pl)
            print(f"{name:42s} {shape:4s} in-sample {insample:+.4f} | held-out {b_all[0]:+.4f} "
                  f"[{b_all[1]:+.4f}, {b_all[2]:+.4f}] | panel-like {b_pl[0]:+.4f} [{b_pl[1]:+.4f}, {b_pl[2]:+.4f}]",
                  flush=True)
        out[name] = row
    json.dump(out, open(HERE / "work/cv480.json", "w"), indent=1)


if __name__ == "__main__":
    main()
