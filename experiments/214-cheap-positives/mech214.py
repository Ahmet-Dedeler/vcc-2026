"""214 — the mechanism behind arm A: WHY does PerturbMatch cost `reach` and `nmae`?

`reach` is a sign-purity walk down our own ordering of the reference's DE genes, and under
`rc_callset` that ordering IS `|H|` (the delivered set is the prescribed set, `rc_emit.py`
:983).  So a `reach` move is a statement about the KEY.  `nmae` reads the delivered log2FC
magnitudes against the reference's, and `docs/BOARD-LEDGER.md` §5 / `docs/LOG.md`'s
`nmae(k,q) = q|k-1| + (1-q)(k+1)` identity says the only thing that matters is the ratio
`k` = our magnitude / theirs, with the optimum at `k = 1` and a slope of +1 above it.

Three statistics, per perturbation and averaged, on the h1-test reference's own BH set:

  1. cumulative sign purity `P(1) / P(3) / P(10)` down the `|H|` order -- the exact object
     BOARD-LEDGER §5 fits, and the one `reach` integrates;
  2. Spearman rho between `|H|` and the reference's `|lfc|` -- how well the key ranks;
  3. the magnitude ratio `k` = median(|our lfc|) / median(|their lfc|) over that pool.

    uv run python mech214.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
from debench import ALPHA, Surface                                       # noqa: E402

ARMS = ("reb_pm50_floor", "reb_pm00_floor", "reb_pm50_ham", "reb_pm00_ham",
        "art_floor", "art_ham")
LABEL = {"reb_pm50_floor": "incumbent geometry (PM 0.50)",
         "reb_pm00_floor": "A  drop PerturbMatch (PM 0)",
         "reb_pm50_ham": "B  Hamilton (PM 0.50)",
         "reb_pm00_ham": "A∘B",
         "art_floor": "artifact (095 geometry, PM 0)",
         "art_ham": "artifact + Hamilton"}


def purity(z, s):
    """P(k) down our own |H| order, over the reference's BH set, per perturbation."""
    absH = z["absH"]
    out = {1: [], 3: [], 10: []}
    ks, ratios, rhos = [], [], []
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size < 10:
            continue
        order = pool[np.argsort(-absH[p, pool], kind="stable")]
        agree = (np.sign(s.pred_lfc[p, order]) == np.sign(s.real_lfc[p, order]))
        cum = np.cumsum(agree) / np.arange(1, agree.size + 1)
        for k in out:
            if cum.size >= k:
                out[k].append(float(cum[k - 1]))
        r = spearmanr(absH[p, pool], np.abs(s.real_lfc[p, pool])).statistic
        if np.isfinite(r):
            rhos.append(float(r))
        num = float(np.median(np.abs(s.pred_lfc[p, pool])))
        den = float(np.median(np.abs(s.real_lfc[p, pool])))
        if den > 0:
            ks.append(num / den)
        ratios.append(float(agree.mean()))
    return dict(P1=float(np.mean(out[1])), P3=float(np.mean(out[3])),
                P10=float(np.mean(out[10])), sign_acc=float(np.mean(ratios)),
                rho_absH_lfc=float(np.mean(rhos)), k=float(np.median(ks)),
                n=len(ratios))


def main() -> None:
    rows = {}
    print("h1-test, reference BH pool, ordered by our own |H| (which IS `rc_callset`'s "
          "delivered order)\n")
    print(f"  {'arm':<32}{'P(1)':>8}{'P(3)':>8}{'P(10)':>8}{'sign acc':>10}"
          f"{'rho(|H|,|lfc|)':>16}{'k = ours/theirs':>17}{'n':>5}")
    for a in ARMS:
        p = WORK / f"de214.{a}.npz"
        if not p.exists():
            continue
        s = Surface(p)
        z = np.load(p, allow_pickle=True)
        r = purity(z, s)
        rows[a] = r
        print(f"  {LABEL[a]:<32}{r['P1']:>8.4f}{r['P3']:>8.4f}{r['P10']:>8.4f}"
              f"{r['sign_acc']:>10.4f}{r['rho_absH_lfc']:>16.4f}{r['k']:>17.4f}{r['n']:>5}")

    a, b = rows["reb_pm00_floor"], rows["reb_pm50_floor"]
    print(f"\n  ARM A moves, on the incumbent's own emitter and at a FIXED median row norm "
          f"(h1parts200's\n  amplitude convention, so this is a pure SHAPE comparison):")
    for k, lab in (("P1", "head sign purity P(1)"), ("P3", "P(3)"), ("P10", "P(10)"),
                   ("sign_acc", "sign accuracy over the whole pool"),
                   ("rho_absH_lfc", "rho(|H|, reference |lfc|)"),
                   ("k", "magnitude ratio k (nmae optimum is 1)")):
        print(f"    {lab:<44}{b[k]:>10.4f} -> {a[k]:>8.4f}   {a[k]-b[k]:>+9.4f}")
    q = 0.537
    dk = a["k"] - b["k"]
    pred = (1 - 2 * q) * dk                      # docs/LOG.md: slope -0.074 below k = 1
    print(f"\n  `nmae(k, q) = q|k-1| + (1-q)(k+1)` (docs/LOG.md 2026-09-03).  Both arms sit "
          f"far BELOW k = 1,\n  where the slope is 1 - 2q = {1-2*q:+.3f}, so the aggregate "
          f"identity predicts d nmae = {pred:+.5f}\n  for arm A's dk = {dk:+.5f}.  The "
          f"MEASURED move is -0.01148 -- 17x larger and of the opposite\n  sign, so arm A's "
          f"`nmae` gain is NOT the aggregate magnitude ratio.  That is 187's finding "
          f"exactly:\n  the uniform-k identity is an artefact (the |t|-weighted per-gene k "
          f"has median 0.12, p99 1.11),\n  so this is a PER-GENE shape effect and must be "
          f"read off the measurement, not the identity.")
    json.dump(rows, open(WORK / "mech214.json", "w"), indent=1, default=float)
    print(f"\nwrote {WORK/'mech214.json'}")


if __name__ == "__main__":
    main()
