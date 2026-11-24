"""188 Phase 1 — fit the magnitude map leave-one-line-out and report ITS SHAPE.

The mandate: "Report the fitted map's shape explicitly: is it convex, concave, S-shaped?
Where does it cross the identity? What fraction of the gated mass does it move up versus
down? The shape IS the finding, whatever the arm's fate."

WHAT IS FITTED, AND ON WHAT
---------------------------
Population: `scored` (the held donor's own columns, minus all 300 panel targets -- exactly
what `pds_cosine` scores) INTERSECTED with the production challenge gate (pb CPM > 10 in all
of A/B/C), because that is the coordinate set the transform is allowed to touch.

Three estimands, all binned-median + PAVA, all leave-one-line-out:
  `signed`  y = r * sign(p)   -- the L1 optimum, derived in `mapfit.py`
  `abs`     y = |r|           -- the mandate's literal phrasing; reported to show the gap
  `de`      y = r * sign(p) restricted to the top-2000 |r| genes of each row, which is the
            selection `nmae`'s own gate imposes (098 §4: the challenge reference calls
            ~2,000 genes per row)
in two spaces: the builder's `h` = log1p(TS*p) and log2 fold change.

KNOWN-ANSWER ROWS
-----------------
  * the identity map (`m(x) = x`) must leave the effect matrix bit-identical after row
    restoration;
  * the fitted map must preserve the `|effect|` ORDERING (monotone by construction), which is
    what makes `keep_top 8000` overlap 1.00000 -- 169 §1's check, re-run here;
  * `median |h|/|r|` on the PLOLO folds must land in the neighbourhood of 169 §6.2's
    0.1266 (HCT116) / 0.2825 (h1-test) -- different surfaces, so the gate is order-of-
    magnitude, and it is reported not asserted.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, HERE)
from mapfit import Map, apply_map, apply_map_lfc, fit_map, to_lfc  # noqa: E402

FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
N_DE = 2000                      # 098 §4's estimate of the challenge reference's per-row calls


def load(fold):
    z = np.load(os.path.join(WORK, f"pair-{fold}.npz"), allow_pickle=True)
    return dict(pred=z["pred"].astype(np.float64), truth=z["truth"].astype(np.float64),
                b=z["b_ctrl"], gate=z["gate"], scored=z["scored"],
                targets=np.asarray([str(t) for t in z["targets"]]))


def chal_gate():
    z = np.load(os.path.join(WORK, "chalframes.npz"))
    return z["gA"] & z["gB"] & z["gC"], {c: (z[f"b{c}"], z[f"g{c}"]) for c in "ABC"}


def pairs(F, cols, space):
    """(x, y_signed, y_abs, de_mask) over every (target, gene) pair on `cols`."""
    P = F["pred"][:, cols]
    R = F["truth"][:, cols]
    if space == "lfc":
        g = np.ones(cols.size, bool)
        P = to_lfc(P, F["b"][cols], g)
        R = to_lfc(R, F["b"][cols], g)
    s = np.sign(P)
    live = (P != 0.0) & (R != 0.0)
    de = np.zeros(R.shape, bool)
    k = min(N_DE, R.shape[1])
    idx = np.argpartition(-np.abs(R), k - 1, axis=1)[:, :k]
    np.put_along_axis(de, idx, True, axis=1)
    return (np.abs(P)[live], (R * s)[live], np.abs(R)[live], de[live],
            dict(n_pairs=int(live.sum()), n_rows=int(P.shape[0]), n_cols=int(cols.size)))


def main() -> None:
    cg, frames = chal_gate()
    print(f"challenge gate (pb CPM > 10 in all of A/B/C): {int(cg.sum())} genes\n")

    F = {f: load(f) for f in FOLDS}
    cols = {}
    for f in FOLDS:
        m = np.zeros(cg.size, bool)
        m[F[f]["scored"]] = True
        cols[f] = np.flatnonzero(m & cg)
        med = np.median(np.abs(F[f]["pred"][:, cols[f]])[
            (F[f]["pred"][:, cols[f]] != 0) & (F[f]["truth"][:, cols[f]] != 0)] /
            np.abs(F[f]["truth"][:, cols[f]])[
            (F[f]["pred"][:, cols[f]] != 0) & (F[f]["truth"][:, cols[f]] != 0)])
        print(f"  {f:<14} {len(F[f]['targets']):>3} targets  {cols[f].size:>5} scored&gated "
              f"cols   median |h|/|r| = {med:.4f}")

    data = {}
    for space in ("h", "lfc"):
        for f in FOLDS:
            data[(space, f)] = pairs(F[f], cols[f], space)

    out = {}
    print(f"\n{'=' * 100}\nLEAVE-ONE-LINE-OUT FITS — the map read on the line that never "
          f"contributed to it")
    for space in ("h", "lfc"):
        for kind in ("signed", "abs", "de"):
            for held in FOLDS:
                xs, ys, ya, de = [], [], [], []
                for f in FOLDS:
                    if f == held:
                        continue
                    x, y, a, d, _ = data[(space, f)]
                    xs.append(x); ys.append(y); ya.append(a); de.append(d)
                x = np.concatenate(xs); y = np.concatenate(ys)
                a = np.concatenate(ya); d = np.concatenate(de)
                if kind == "abs":
                    mm = fit_map(x, a, meta=dict(space=space, kind=kind, held=held))
                elif kind == "de":
                    mm = fit_map(x[d], y[d], meta=dict(space=space, kind=kind, held=held))
                else:
                    mm = fit_map(x, y, meta=dict(space=space, kind=kind, held=held))
                out[f"{space}|{kind}|{held}"] = mm.to_dict()
        # the pooled (all four folds) map, for shape reporting only -- never used to score
        for kind in ("signed", "abs", "de"):
            x = np.concatenate([data[(space, f)][0] for f in FOLDS])
            y = np.concatenate([data[(space, f)][1] for f in FOLDS])
            a = np.concatenate([data[(space, f)][2] for f in FOLDS])
            d = np.concatenate([data[(space, f)][3] for f in FOLDS])
            src = dict(signed=(x, y), abs=(x, a), de=(x[d], y[d]))[kind]
            out[f"{space}|{kind}|POOLED"] = fit_map(*src, meta=dict(space=space, kind=kind,
                                                                   held="POOLED")).to_dict()

    with open(os.path.join(WORK, "maps.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    # ------------------------------------------------------------------ shape report
    for space in ("h", "lfc"):
        print(f"\n{'-' * 100}\nSPACE = {space}   POOLED map, per quantile bin of |p| "
              f"(the shape; LOO maps are within the last column's spread)")
        M = {k: Map.from_dict(out[f"{space}|{k}|POOLED"]) for k in ("signed", "abs", "de")}
        L = {k: [Map.from_dict(out[f"{space}|{k}|{h}"]) for h in FOLDS]
             for k in ("signed", "abs", "de")}
        print(f"  {'pctile':>7} {'|p| (x)':>10} {'m_signed':>10} {'m/x':>8} "
              f"{'m_abs':>10} {'m_de':>10} {'LOO sd(m_signed)':>18}")
        xs = M["signed"].xs
        for i in range(0, len(xs), 3):
            xv = xs[i]
            sd = float(np.std([mm(xv) for mm in L["signed"]]))
            print(f"  {100.0*(i+0.5)/len(xs):>7.1f} {xv:>10.5f} {M['signed'](xv):>10.5f} "
                  f"{M['signed'](xv)/max(xv,1e-12):>8.3f} {M['abs'](xv):>10.5f} "
                  f"{M['de'](xv):>10.5f} {sd:>18.5f}")
        xv = xs[-1]
        sd = float(np.std([mm(xv) for mm in L["signed"]]))
        print(f"  {100.0*(len(xs)-0.5)/len(xs):>7.1f} {xv:>10.5f} {M['signed'](xv):>10.5f} "
              f"{M['signed'](xv)/max(xv,1e-12):>8.3f} {M['abs'](xv):>10.5f} "
              f"{M['de'](xv):>10.5f} {sd:>18.5f}")

        # convexity: sign of the second difference of m on the quantile grid
        for k in ("signed", "abs", "de"):
            m = np.asarray(out[f"{space}|{k}|POOLED"]["ms"])
            x = np.asarray(out[f"{space}|{k}|POOLED"]["xs"])
            d1 = np.diff(m) / np.maximum(np.diff(x), 1e-12)
            d2 = np.diff(d1)
            zero = float(m[0])
            print(f"  {k:<7}: m(x_min) = {zero:.5f}, slope rises on "
                  f"{100.0*np.mean(d2 > 0):.0f}% of adjacent bin pairs "
                  f"(CONVEX if ~100, CONCAVE if ~0), first slope {d1[0]:.4f}, "
                  f"last slope {d1[-1]:.4f}")

    # ------------------------------------------------------------------ delivered shape
    print(f"\n{'=' * 100}\nWHAT THE MAP DOES TO THE EFFECT MATRIX (row L2 restored) — held "
          f"fold, LOO map, challenge frame A")
    bA, gA = frames["A"]
    print(f"  {'fold':<14}{'space':>6}{'kind':>8}{'cross k=1 at pct':>18}{'L1 up':>9}"
          f"{'L1 down':>9}{'top-8000 overlap':>18}{'corr with q=1.15':>18}")
    sys.path.insert(0, os.path.join(os.path.dirname(HERE), "169-sharpening-six-member"))
    from sharp import sharpen, sharpen_lfc  # noqa: E402
    rows = []
    for held in FOLDS:
        P = F[held]["pred"]
        for space in ("h", "lfc"):
            for kind in ("signed", "abs", "de"):
                mm = Map.from_dict(out[f"{space}|{kind}|{held}"])
                if space == "h":
                    Q = apply_map(P, mm, gA, mode="row")
                    Pw = sharpen(P, 1.15, gA, mode="row")
                else:
                    Q = apply_map_lfc(P, mm, bA, gA, restore="h")
                    Pw = sharpen_lfc(P, 1.15, bA, gA, mode="row")
                G = np.where(gA[None, :], P, 0.0)
                ratio = np.abs(Q) / np.maximum(np.abs(G), 1e-300)
                live = np.abs(G) > 0
                # where does the delivered magnitude cross the input's? (identity crossing)
                up = float(np.abs(G)[live & (ratio > 1)].sum() / np.abs(G)[live].sum())
                # percentile of |p| at which ratio crosses 1, per row then median
                pct = []
                for i in range(P.shape[0]):
                    li = live[i]
                    if li.sum() < 10:
                        continue
                    o = np.argsort(np.abs(G[i][li]))
                    rr = ratio[i][li][o]
                    c = np.flatnonzero(rr > 1.0)
                    pct.append(100.0 * c[0] / rr.size if c.size else 100.0)
                # keep_top ordering preservation
                k = min(8000, int(gA.sum()))
                s0 = set(np.argsort(-np.abs(G[0]))[:k].tolist())
                s1 = set(np.argsort(-np.abs(Q[0]))[:k].tolist())
                ov = len(s0 & s1) / k
                cc = float(np.corrcoef(Q[live].ravel(), Pw[live].ravel())[0, 1])
                rows.append((held, space, kind, float(np.median(pct)), up, 1 - up, ov, cc))
                print(f"  {held:<14}{space:>6}{kind:>8}{np.median(pct):>18.1f}{up:>9.3f}"
                      f"{1-up:>9.3f}{ov:>18.5f}{cc:>18.4f}")

    with open(os.path.join(WORK, "shape188.json"), "w") as fh:
        json.dump(dict(rows=rows), fh, indent=1)
    print(f"\nwrote {WORK}/maps.json and shape188.json")


if __name__ == "__main__":
    sys.exit(main())
