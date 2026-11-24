"""169 — `nmae` with the REFERENCE GATE DEPTH controlled, and in both candidate spaces.

Experiment 098 §4 is the reason this script exists, and it was nearly missed:

  "every `nmae` number this project has ever taken was read on a reference gate 2-20x
   shallower than the board's, and the SIGN of the amplitude response flips across that
   range ... the challenge's reference calls ~2,000 genes per row; our harnesses call
   89 (orion) and 1,120 (h1) at BH 0.05."

That is a rule-9-class artefact on `nmae`: the gate is reference-only, its SIZE is nothing
like the challenge's on any harness we hold, and 098 measured arm rankings inverting across
it.  Everything in `cliff169.py` / `nmaechal169.py` was read at native BH 0.05, so it has to
be re-read with the depth swept -- exactly as 157 §3b had to sweep `reach`'s pool.

098 also already screened a sharpening-shaped arm (its route 5, "sharpen |lfc|^+0.25 at held
A"), and read `-0.0027 / -0.0036` at BH 0.05 and `+0.0005 / -0.0028` at the deep gates --
i.e. NEGATIVE (an improvement), where `cliff169.py` reads `+0.0144` at q = 1.25.  One of the
two must be wrong, and the candidate explanation is the SPACE the power is applied in:

  `h`   -- the builder's log1p-pseudobulk effect, which is what `profile_from_effect` takes
           and therefore the only thing a builder flag can transform.  This is 169's arm.
  `lfc` -- the delivered log2 fold change, which is what 098 transformed.  For a
           high-expressed gene `h ~ ln2 * lfc`, but for a low-expressed one `h` is nearly
           linear in the profile while `lfc` is still logarithmic, so the two reweight the
           gate's low-expression half completely differently.

Both are measured here, at every depth, on both surfaces.  If they disagree the disagreement
is the finding; if the ranking is stable the arm survives its own strongest prior.

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
from sharp import sharpen, sharpen_lfc                   # noqa: E402

EPS = 1e-9
QS = (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30)
BHS = (0.05, 0.15, 0.35, 0.60, 0.85)
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")


def nmae_at(real_padj, real_lfc, M, bh):
    """`nmae` and the saturated-L1 share on the gate `real_padj < bh`."""
    out, satL1, allL1 = [], 0.0, 0.0
    sizes = []
    for p in range(real_lfc.shape[0]):
        s = np.flatnonzero((real_padj[p] < bh) & np.isfinite(real_padj[p])
                           & np.isfinite(real_lfc[p]) & (real_lfc[p] != 0.0))
        if s.size < 10:
            continue
        r = real_lfc[p, s]
        h = np.where(np.isfinite(M[p, s]), M[p, s], 0.0)
        den = np.abs(r).sum()
        if not np.isfinite(den) or den == 0:
            continue
        out.append(np.abs(h - r).sum() / den)
        sizes.append(s.size)
        sat = np.abs(h) > np.abs(r)
        satL1 += np.abs(h)[sat].sum()
        allL1 += np.abs(h).sum()
    return float(np.mean(out)), float(np.mean(sizes)), satL1 / max(allL1, 1e-30)


def hct116_lfcs():
    D, T, b, b_ref, gate, mask, ptar = production_effect("hct116")
    s = Surface(os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz"))
    z = np.load(FROZEN % "hct116", allow_pickle=True)
    ctrl = np.load(os.path.join(REPO, "experiments/156-delivered-geometry/work",
                                "truth-hct116.npz"), allow_pickle=True)["ctrl_counts"]
    ccpm = 1e6 * ctrl.astype(np.float64) / ctrl.astype(np.float64).sum()
    gpos = {g: i for i, g in enumerate(s.genes)}
    cols = np.asarray([gpos.get(str(g), -1) for g in z["genes"]])
    prow = {pp: i for i, pp in enumerate(s.perts)}
    rows = np.asarray([prow.get(pp, -1) for pp in ptar])
    ok, keep = cols >= 0, rows >= 0

    def place(L):
        M = np.zeros((len(s.perts), len(s.genes)))
        M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep), np.flatnonzero(ok))]
        return M

    def deliver(eff):
        lin = np.clip(np.expm1(b[None, :] + eff), 0.0, None)
        p = lin / lin.sum(1, keepdims=True)
        return np.log2((1e6 * p + EPS) / (ccpm[None, :] + EPS))

    def in_h(q):
        return place(deliver(mass_fix(sharpen(D, q, gate, mode="row"), b, gate)))

    def in_lfc(q):
        # the buildable lfc-space arm: power the log2 fold change, then restore the
        # BUILDER-space row norm so the arm carries no amplitude change (`sharpen_lfc`).
        return place(deliver(mass_fix(
            sharpen_lfc(D, q, b, gate, mode="row", restore="h"), b, gate)))

    return "hct116", s.real_padj, s.real_lfc, in_h, in_lfc


def h1test_lfcs():
    path = os.path.join(HERE, "de169.q1.00_mf.npz")
    if not os.path.exists(path):
        return None
    E = np.load(os.path.join(WORK, "h1test-effect.npz"), allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    gate = E["gate"]
    egenes = np.asarray([str(g) for g in E["genes"]])
    etar = np.asarray([str(t) for t in E["targets"]])
    cc = np.load(os.path.join(WORK, "h1test-cellcpm.npy"))
    s = Surface(path)
    gpos = {g: i for i, g in enumerate(s.genes)}
    cols = np.asarray([gpos.get(g, -1) for g in egenes])
    prow = {p: i for i, p in enumerate(s.perts)}
    rows = np.asarray([prow.get(p, -1) for p in etar])
    ok, keep = cols >= 0, rows >= 0

    def place(L):
        M = np.zeros((len(s.perts), len(s.genes)))
        M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep), np.flatnonzero(ok))]
        return M

    def deliver(eff):
        lin = np.clip(np.expm1(b[None, :] + eff), 0.0, None)
        p = lin / lin.sum(1, keepdims=True)
        return np.log2((1e6 * p + EPS) / (cc[None, :] + EPS))

    def in_h(q):
        return place(deliver(sharpen(H0, q, gate, mode="row")))

    def in_lfc(q):
        return place(deliver(sharpen_lfc(H0, q, b, gate, mode="row", restore="h")))

    return "h1-test", s.real_padj, s.real_lfc, in_h, in_lfc


def main() -> None:
    out = {}
    for maker in (hct116_lfcs, h1test_lfcs):
        got = maker()
        if got is None:
            print("h1-test tables absent -- skipping")
            continue
        name, padj, rlfc, in_h, in_lfc = got
        print(f"\n{'=' * 118}\n{name}")
        for space, fn in (("h  (log1p effect -- 169's BUILDABLE arm)", in_h),
                          ("lfc (log2 fold change, builder-norm restored -- the ALTERNATIVE arm)", in_lfc)):
            print(f"\n  space: {space}")
            print(f"  {'BH':>6} {'mean |S_p|':>11} {'nmae q=1':>10} {'sat L1':>8} "
                  + " ".join(f"{'d@' + f'{q:.2f}':>9}" for q in QS[1:]))
            Ms = {q: fn(q) for q in QS}
            for bh in BHS:
                n0, sz, sh = nmae_at(padj, rlfc, Ms[1.00], bh)
                cells = []
                for q in QS[1:]:
                    nq, _, _ = nmae_at(padj, rlfc, Ms[q], bh)
                    cells.append(f"{nq - n0:+9.5f}")
                    out[f"{name}|{space[:3].strip()}|{bh}|{q}"] = nq - n0
                print(f"  {bh:>6.2f} {sz:>11.0f} {n0:>10.5f} {sh:>8.4f} " + " ".join(cells))

    print(f"\n{'=' * 118}\nWHAT 098 MEASURED, for comparison (its route 5, "
          f"`sharpen |lfc|^+0.25 at held A`, on a12):")
    print(f"  h1 @BH0.05 -0.0027 | h1 @BH0.35 (|S| 2600) +0.0005 | "
          f"orion @BH0.05 -0.0036 | orion @BH0.60 (|S| 2204) -0.0028")
    print(f"\nThe challenge's own gate size is NOT directly measurable. Two live estimates:")
    print(f"  ~240 per row  (164 §3.3, from a21's call set)   and")
    print(f"  ~2,000 per row (098 (d2), back-solved from v14/v15's `fid` coverage);")
    print(f"  AGENTS.md rule 9 flags the 093 `n_real ~ 2,670` inference as CIRCULAR.")
    print(f"So the honest reading is the whole ladder, not one row of it.")

    with open(os.path.join(WORK, "nmaegate169.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)


if __name__ == "__main__":
    sys.exit(main())
