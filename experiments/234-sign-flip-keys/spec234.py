"""234 — the SPECIFICATION: how good would a flip key have to be, in a currency any future
head can be measured in?  (220 §7's instrument, moved onto the sign-flip channel.)

THE EXACT IDENTITY.  Let `o` be a per-(target, gene) bit that equals `1[our sign is right]`
with probability `a` and is flipped otherwise, and flip every coordinate with `o = 0`.  Then

    E[Dq_w] = -(1-a) q_w + a (1 - q_w) = a - q_w

so a binary flip key pays IFF its per-pair accuracy exceeds our own delivered sign accuracy,
and rule 13's gate needs `a = q_w + 0.00476`.  The identity is verified numerically below
against 20 seeds per rung -- it is not assumed.

THE OTHER HALF, which 220 never had to face.  `q_w` is measured on the REFERENCE's DE pool,
which is not knowable at prediction time, so a real key flips over the whole gene axis.  The
flips that land outside the pool earn nothing and cost cosine.  This script therefore also
reports the per-key CALIBRATION CURVE -- `P(correct | key decile)` -- because a flip rule on
a key can only ever pay in bins where that curve dips BELOW 0.5, and it reports each key's
AUC and correlation with correctness on PLOLO so the two surfaces can be compared.

    uv run python spec234.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common234 import CANCER, KA220_AUC, KA220_CORR, QW_NEEDED, RATE_PER_QW, SEED, WORK  # noqa: E402
from screen234 import Fold  # noqa: E402

ALL_FOLDS = CANCER + ("kolf2",)
A_GRID = (0.50, 0.52, 0.55, 0.57, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00)
NSEED = 20
DECILES = 10


def wauc(sig, correct, w):
    """|t|-weighted Mann-Whitney AUC of `sig` for predicting `correct`.

    `sum_{i in correct} sum_{j in wrong} w_i w_j [1(s_i > s_j) + 0.5 1(s_i == s_j)]`
    over `W1 * W0`.  Ties get half credit; only NEGATIVE weight is accumulated below a
    positive, which is the part the naive rank-sum form gets wrong when the weights differ
    between the two classes.
    """
    o = np.argsort(sig, kind="stable")
    s, c, ww = sig[o], correct[o].astype(np.float64), w[o]
    wneg = ww * (1.0 - c)
    cumneg = np.concatenate([[0.0], np.cumsum(wneg)])
    below = np.empty(s.size)
    i = 0
    while i < s.size:
        j = i
        while j + 1 < s.size and s[j + 1] == s[i]:
            j += 1
        below[i:j + 1] = cumneg[i] + 0.5 * (cumneg[j + 1] - cumneg[i])
        i = j + 1
    W1 = float((ww * c).sum())
    W0 = float(wneg.sum())
    if W1 <= 0 or W0 <= 0:
        return float("nan")
    return float((ww * c * below).sum() / (W1 * W0))


def main() -> None:
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    folds = {h: Fold(h) for h in ALL_FOLDS}
    C = [folds[h] for h in CANCER]
    q_cancer = float(np.mean([f.q_w for f in C]))
    audit = dict(q_w_cancer=q_cancer,
                 a_required={k: q_cancer + v for k, v in QW_NEEDED.items()})

    # ------------------------------------------------- the noisy-oracle identity, verified
    print(f"\nTHE NOISY-ORACLE SPECIFICATION.  Our delivered `q_w` = {q_cancer:.5f}; the gate "
          f"needs d q_w >= +{QW_NEEDED['k1.384']:.5f}\n"
          f"  => a binary flip key must be RIGHT about 'is our sign wrong' with per-pair "
          f"accuracy a >= {q_cancer + QW_NEEDED['k1.384']:.5f}\n")
    print(f"{'a':>6}{'measured d q_w':>17}{'a - q_w (theory)':>19}{'corr(o, correct)':>19}"
          f"{'AUC':>8}{'overall':>10}")
    print("-" * 79)
    rows = []
    for a in A_GRID:
        dq, cr, au = [], [], []
        for sd in range(NSEED if a < 1.0 else 1):
            r = np.random.default_rng(1000 + sd)
            per_fold = []
            for f in C:
                flip_bit = r.random(f.truth.shape) > a
                o = np.where(flip_bit, ~f.correct, f.correct)          # the noisy key
                sel = ~o & f.pool                                       # flip where o says wrong
                per_fold.append(float(np.mean((f.g * sel).sum(1)[f.live])))
                if sd == 0:
                    m = f.pool
                    cr.append(float(np.corrcoef(o[m].astype(float),
                                                f.correct[m].astype(float))[0, 1]))
                    au.append(wauc(o[m].astype(np.float64), f.correct[m], f.w[m]))
            dq.append(float(np.mean(per_fold)))
        d = float(np.mean(dq))
        rows.append(dict(a=a, dqw=d, theory=a - q_cancer, corr=float(np.mean(cr)),
                         auc=float(np.mean(au)), sd=float(np.std(dq, ddof=1)) if len(dq) > 1
                         else 0.0, overall=d * RATE_PER_QW["k1.384"]))
        print(f"{a:>6.2f}{d:>+17.5f}{a-q_cancer:>+19.5f}{np.mean(cr):>+19.4f}"
              f"{np.mean(au):>8.4f}{d*RATE_PER_QW['k1.384']:>+10.5f}")
    err = max(abs(r["dqw"] - r["theory"]) for r in rows)
    print(f"\n  identity `E[d q_w] = a - q_w` holds to {err:.2e} over {len(rows)} rungs "
          f"x {NSEED} seeds  ->  REQUIRED a = {q_cancer + QW_NEEDED['k1.384']:.5f} "
          f"(corr ~ {2*(q_cancer + QW_NEEDED['k1.384'])-1:+.4f}, AUC ~ "
          f"{q_cancer + QW_NEEDED['k1.384']:.4f})")
    audit["noisy_oracle"] = rows
    audit["identity_max_err"] = err

    # -------------------------------------- where every key we hold actually sits, on PLOLO
    print(f"\n{'='*104}\nEVERY KEY WE HOLD, MEASURED ON PLOLO (|t|-weighted, inside the "
          f"`wald89` pool, 3 cancer folds)\n{'='*104}")
    hdr = (f"{'key':<12}{'corr':>9}{'AUC':>8}{'min decile acc':>16}"
           f"{'mass below 0.5':>16}{'best decile Dq_w':>18}{'220 h1-test corr':>18}")
    print(hdr)
    print("-" * len(hdr))
    kt = []
    for kn in sorted(C[0].keys):
        cr, au, mn, mb, bq = [], [], [], [], []
        for f in C:
            m = f.pool
            k = f.keys[kn][m].astype(np.float64)
            c = f.correct[m]
            w = f.w[m]
            if k.std() > 0:
                cr.append(float(np.corrcoef(k, c.astype(float))[0, 1]))
                au.append(wauc(k, c, w))
            e = np.unique(np.quantile(k, np.linspace(0, 1, DECILES + 1)))
            e[0], e[-1] = -np.inf, np.inf
            b = np.digitize(k, e[1:-1])
            num = np.bincount(b, weights=w * c, minlength=len(e) - 1)
            den = np.bincount(b, weights=w, minlength=len(e) - 1)
            acc = np.where(den > 0, num / np.maximum(den, 1e-300), 1.0)
            mn.append(float(acc.min()))
            mb.append(float(den[acc < 0.5].sum() / max(den.sum(), 1e-300)))
            # the best SINGLE decile as a flip set, evaluated in-fold (optimistic)
            best = 0.0
            for bi in range(len(acc)):
                sel = np.zeros(f.truth.shape, dtype=bool)
                sel[m] = (b == bi)
                best = max(best, float(np.mean((f.g * sel).sum(1)[f.live])))
            bq.append(best)
        r = dict(key=kn, corr=float(np.mean(cr)) if cr else float("nan"),
                 auc=float(np.mean(au)) if au else float("nan"),
                 min_decile_acc=float(np.mean(mn)),
                 mass_below_half=float(np.mean(mb)),
                 best_decile_dqw=float(np.mean(bq)))
        kt.append(r)
        ref = KA220_CORR.get(kn)
        print(f"{kn:<12}{r['corr']:>+9.4f}{r['auc']:>8.4f}{r['min_decile_acc']:>16.4f}"
              f"{r['mass_below_half']:>16.4f}{r['best_decile_dqw']:>+18.5f}"
              + (f"{ref:>+18.4f}" if ref is not None else f"{'-':>18}"))
    audit["keys_on_plolo"] = kt
    print(f"\n  REQUIRED (to clear the gate at full dose): corr "
          f"{2*(q_cancer + QW_NEEDED['k1.384'])-1:+.4f}   AUC "
          f"{q_cancer + QW_NEEDED['k1.384']:.4f}")
    print(f"  220 §7's requirement on the ORDERING channel, for scale: corr +0.30 / "
          f"AUC 0.639  ({KA220_CORR['absH']:+.4f} measured for |H| there)")

    audit["seconds"] = time.time() - t0
    (WORK / "spec234.json").write_text(json.dumps(audit, indent=1))
    print(f"\n-> {WORK/'spec234.json'}  ({audit['seconds']:.0f}s)")


if __name__ == "__main__":
    main()
