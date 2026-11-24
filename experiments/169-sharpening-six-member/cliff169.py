"""169 — the `nmae` cliff, decomposed.  Is the optimal `q` derivable or only fitted?

168 recorded a 9x jump in `d nmae` for a 0.10 change in `q` and its own first explanation --
"we start overshooting the reference" -- was REFUTED there: the overshoot fraction rises only
1.12% -> 3.40% while the MEDIAN |h|/|r| falls and the mean |h| moves 3.4x TOWARD the mean |r|.
That mandate asked for `Sum|h-r|` split by sign-correct/sign-wrong and overshoot/undershoot.

THE DECOMPOSITION IS AN IDENTITY, NOT A FIT
-------------------------------------------
For reals `h`, `r`:  `|h - r| = |h| + |r| - 2m`,  `m = min(|h|,|r|) * 1[sign h == sign r]`.
Summing over the reference's DE set for one perturbation and dividing by `Sum|r|`:

        nmae_p  =  1  +  ( Sum_g |h_g|  -  2 Sum_g m_g )  /  Sum_g |r_g|              (I)

So the metric charges the L1 BUDGET we spend and refunds TWICE the SIGN-CORRECT OVERLAP,
where the overlap is capped at `|r|` per gene.  There are no other terms.  (I) is checked
against `debench.Surface.nmae()` below and must agree to machine precision; if it does not,
everything after it is void.

WHAT (I) IMMEDIATELY GIVES
--------------------------
Moving one unit of L1 onto a gene changes `nmae_p * Sum|r|` by `1 - 2*(dm/d|h|)`, and
`dm/d|h|` is exactly `1[sign correct AND |h| < |r|]`.  So:

        an L1 reallocation helps `nmae` iff it lands where the SIGN-CORRECT-AND-UNSATURATED
        rate exceeds 1/2, and hurts otherwise.

That is a threshold with a known value, not a fitted knee.  Sharpening moves L1 from the
small entries to the large ones, so the question becomes empirical but *pointed*: the
marginal efficiency `e(decile) = P(correct & |h| < |r|)` is measured per decile of `|h|`,
and the predicted optimal `q` is the one where the L1 taken from the tail (efficiency
`e_tail`) equals the L1 given to the head (efficiency `e_head`) at the 1/2 line.

Both the exact per-q decomposition and the derived optimum are printed, and the derived
optimum is compared against the measured argmin on both surfaces.

Deterministic.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, HERE)
from debench import Surface                              # noqa: E402
from massfix import mass_fix, production_effect          # noqa: E402
from sharp import sharpen                                # noqa: E402

EPS = 1e-9
QS = (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.40, 1.50)
NDEC = 10
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")


# ------------------------------------------------------------------ surfaces
def surface_hct116():
    """Delivered lfc on the Orion HCT116 harness, against 095's cached reference DE table."""
    D, T, b, b_ref, gate, mask, ptar = production_effect("hct116")
    s = Surface(os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz"))
    z = np.load(FROZEN % "hct116", allow_pickle=True)
    ctrl = np.load(os.path.join(REPO, "experiments/156-delivered-geometry/work",
                                "truth-hct116.npz"), allow_pickle=True)["ctrl_counts"]
    ctrl_cpm = 1e6 * ctrl.astype(np.float64) / ctrl.astype(np.float64).sum()
    gpos = {g: i for i, g in enumerate(s.genes)}
    cols = np.asarray([gpos.get(str(g), -1) for g in z["genes"]])
    prow = {pp: i for i, pp in enumerate(s.perts)}
    rows = np.asarray([prow.get(pp, -1) for pp in ptar])

    def lfc(q):
        eff = mass_fix(sharpen(D, q, gate, mode="row"), b, gate)
        lin = np.clip(np.expm1(b[None, :] + eff), 0.0, None)
        p = lin / lin.sum(1, keepdims=True)
        L = np.log2((1e6 * p + EPS) / (ctrl_cpm[None, :] + EPS))
        ok, keep = cols >= 0, rows >= 0
        M = np.full((len(s.perts), len(s.genes)), 0.0)
        M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep),
                                                   np.flatnonzero(ok))]
        return M

    return "hct116", s, lfc


def surface_h1test():
    """Delivered lfc on h1-test, from the effect `probe_effect.py` recovered, against the
    reference DE table `realize169.py` computed from the harness's own perturbed cells."""
    path = os.path.join(HERE, "de169.q1.00_mf.npz")
    if not os.path.exists(path):
        return None
    E = np.load(os.path.join(WORK, "h1test-effect.npz"), allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    gate = E["gate"]
    egenes = np.asarray([str(g) for g in E["genes"]])
    etar = np.asarray([str(t) for t in E["targets"]])
    ctrl = E["ctrl_totals"].astype(np.float64)
    s = Surface(path)
    # `real_lfc` in the de169 tables is against the harness's MEAN-PER-CELL cpm, so the
    # predicted side must use the same denominator or the two are not comparable.
    import anndata as ad
    cc = np.load(os.path.join(WORK, "h1test-cellcpm.npy")) if os.path.exists(
        os.path.join(WORK, "h1test-cellcpm.npy")) else None
    if cc is None:
        raise SystemExit("run `nmae_h1_prep.py` first (caches the mean-per-cell control CPM)")
    del ad
    gpos = {g: i for i, g in enumerate(s.genes)}
    cols = np.asarray([gpos.get(g, -1) for g in egenes])
    prow = {p: i for i, p in enumerate(s.perts)}
    rows = np.asarray([prow.get(p, -1) for p in etar])
    ok, keep = cols >= 0, rows >= 0

    def lfc(q):
        eff = np.stack([sharpen(H0[i], q, gate, mode="row") for i in range(H0.shape[0])])
        lin = np.clip(np.expm1(b[None, :] + eff), 0.0, None)
        p = lin / lin.sum(1, keepdims=True)
        L = np.log2((1e6 * p + EPS) / (cc[None, :] + EPS))
        M = np.full((len(s.perts), len(s.genes)), 0.0)
        M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep),
                                                   np.flatnonzero(ok))]
        return M

    return "h1-test", s, lfc


# ------------------------------------------------------------------ the identity
def decompose(s, M):
    """(I) per perturbation, plus the per-decile marginal-efficiency table."""
    num, tot_b, tot_m, tot_r = [], 0.0, 0.0, 0.0
    dec_h = np.zeros(NDEC)
    dec_m = np.zeros(NDEC)
    dec_e = np.zeros(NDEC)        # sum of 1[correct & |h|<|r|] weighted by |h|
    dec_n = np.zeros(NDEC)
    over_correct = over_wrong = under_correct = under_wrong = 0.0
    for p in range(len(s.perts)):
        i = np.flatnonzero(s.real_sig[p] & np.isfinite(s.real_lfc[p]))
        if i.size < 10:
            continue
        r = s.real_lfc[p, i]
        h = np.where(np.isfinite(M[p, i]), M[p, i], 0.0)
        den = np.abs(r).sum()
        if not np.isfinite(den) or den == 0:
            continue
        same = np.sign(h) == np.sign(r)
        m = np.where(same, np.minimum(np.abs(h), np.abs(r)), 0.0)
        num.append((np.abs(h).sum() - 2 * m.sum()) / den)
        tot_b += np.abs(h).sum()
        tot_m += m.sum()
        tot_r += den
        sat = np.abs(h) >= np.abs(r)
        over_correct += np.abs(h)[same & sat].sum()
        over_wrong += np.abs(h)[~same & sat].sum()
        under_correct += np.abs(h)[same & ~sat].sum()
        under_wrong += np.abs(h)[~same & ~sat].sum()
        a = np.abs(h)
        if a.max() > 0:
            k = np.clip((np.searchsorted(np.sort(a), a, side="right") * NDEC) // a.size,
                        0, NDEC - 1)
            np.add.at(dec_h, k, a)
            np.add.at(dec_m, k, m)
            np.add.at(dec_e, k, a * (same & ~sat))
            np.add.at(dec_n, k, 1.0)
    return dict(nmae=1.0 + float(np.mean(num)), budget=tot_b, overlap=tot_m, denom=tot_r,
                over_correct=over_correct, over_wrong=over_wrong,
                under_correct=under_correct, under_wrong=under_wrong,
                dec_h=dec_h, dec_m=dec_m, dec_e=dec_e, dec_n=dec_n)


def main() -> None:
    out = {}
    for maker in (surface_hct116, surface_h1test):
        got = maker()
        if got is None:
            print("h1-test DE table not built yet -- skipping that surface")
            continue
        name, s, lfc = got
        print(f"\n{'=' * 116}\n{name}: {len(s.perts)} perturbations, {s.genes.size} genes")

        # ---- known-answer gate: (I) must reproduce debench's own nmae ----
        M0 = lfc(1.0)
        d0 = decompose(s, M0)
        ref = s.nmae(pred_lfc=M0)
        print(f"  IDENTITY GATE: debench nmae {ref:.10f} vs decomposition "
              f"{d0['nmae']:.10f}  |d| {abs(ref - d0['nmae']):.2e}")
        if abs(ref - d0["nmae"]) > 1e-9:
            sys.exit("ABORT: the decomposition is not the metric")

        rows = []
        print(f"\n  {'q':>5} {'nmae':>9} {'d nmae':>9} {'L1 budget':>11} {'d bud %':>9} "
              f"{'overlap':>10} {'d ovl %':>9} {'2.ovl/bud':>10} "
              f"{'|h| correct-sat':>16} {'|h| wrong':>11}")
        for q in QS:
            d = decompose(s, lfc(q))
            rows.append(dict(q=q, **{k: (v.tolist() if isinstance(v, np.ndarray) else v)
                                     for k, v in d.items()}))
            print(f"  {q:>5.2f} {d['nmae']:9.5f} {d['nmae'] - d0['nmae']:+9.5f} "
                  f"{d['budget']:11.2f} {100 * (d['budget'] / d0['budget'] - 1):+8.2f}% "
                  f"{d['overlap']:10.2f} {100 * (d['overlap'] / d0['overlap'] - 1):+8.2f}% "
                  f"{2 * d['overlap'] / d['budget']:10.4f} "
                  f"{d['over_correct'] + d['over_wrong']:16.2f} "
                  f"{d['over_wrong'] + d['under_wrong']:11.2f}")
        out[name] = rows

        print(f"\n  MARGINAL EFFICIENCY by decile of |h| at q=1 -- `e` is the share of that "
              f"decile's L1\n  landing on a gene that is sign-correct AND still below |r|. "
              f"L1 moved onto a decile\n  with e > 0.5 IMPROVES nmae; e < 0.5 makes it worse. "
              f"The 1/2 line is the whole mechanism.")
        print(f"  {'decile':>7} {'share of L1':>12} {'e':>8} {'2e-1':>8}")
        for k in range(NDEC):
            e = d0["dec_e"][k] / max(d0["dec_h"][k], 1e-30)
            print(f"  {k + 1:>7} {d0['dec_h'][k] / d0['dec_h'].sum():12.4f} {e:8.4f} "
                  f"{2 * e - 1:+8.4f}")

        print(f"\n  WHERE THE L1 ACTUALLY MOVES (q=1 -> each q), by decile, and the price "
              f"the 1/2 line puts on it")
        print(f"  {'q':>5} " + " ".join(f"{'d' + str(k + 1):>9}" for k in range(NDEC))
              + f" {'predicted d excess':>19} {'actual d excess':>16}")
        for row in rows:
            dh = np.asarray(row["dec_h"]) - d0["dec_h"]
            pred = float(np.sum(dh * (1.0 - 2.0 * (d0["dec_e"] / np.maximum(d0["dec_h"], 1e-30)))))
            act = (row["budget"] - 2 * row["overlap"]) - (d0["budget"] - 2 * d0["overlap"])
            print(f"  {row['q']:>5.2f} " + " ".join(f"{x:9.1f}" for x in dh)
                  + f" {pred:19.2f} {act:16.2f}")

    with open(os.path.join(WORK, "cliff169.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'cliff169.json')}")


if __name__ == "__main__":
    sys.exit(main())
