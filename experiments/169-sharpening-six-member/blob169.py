"""169 — adjudicating the two inboxes: is the sharpening `pds` gain an Orion blob artefact?

`INBOX-from-168.md` UPDATE 3 closes the family as a "TERMINAL NO-GO" on the grounds that
h1-test reads GLOBAL q = 1.15 at `+0.00112` against the Orion pair's `+0.01145` — a 10.2x
overstatement — and lays down "never size an arm on Orion; use h1-test".
`INBOX-from-161.md` asks whether 161's `+0.0206` mean raw `pds` for `|C|^2` across five lines
(K562 / RPE1 / Jurkat / HCT116 / HEK293T) is the same overstatement re-measured, and proposes
the discriminator: each surface's **target-independent energy share** (h1-test 3.6 %,
orion-hct116 26.4 %, orion-hek293t 44.1 %).

Both are answerable from what PLOLO already computed, and the answer is not the one either
inbox assumed:

  * The `pds` a `pds` arm must be sized on is `pds_cosine` in the SCORER's space — the delivered
    log-pseudobulk effect over all 18,533 genes minus the 300 panel coordinates, ranked against
    the whole panel's decoys. `debench.Surface.pds()` on an h1-test DE table is a DIFFERENT
    statistic: log2 fold change, a CPM-5 gate, and **85** decoys against the challenge's 300.
    It has never been validated against a board `pds` point.
  * PLOLO has: 6/6 board signs, LOO RMS 0.001494 scaled (162). Its four folds are
    `k562_gwps`, `orion_hct116`, `orion_hek293t` and `kolf2` — **two of them are not Orion**,
    and one (`kolf2`) is the leak-free production fold.

So: report each PLOLO fold's target-independent energy share next to its own Δ`pds` at
q = 1.15. If the gain tracked the blob, the two non-Orion folds would show it.

Deterministic, no emitter, no cell-eval2, no build, no submission. Nothing outside this
directory is written.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
P162 = os.path.join(REPO, "experiments/162-retrodicting-validator")
sys.path.insert(0, P162)
import plolo                                        # noqa: E402
from plolo import FOLDS, build_fold                 # noqa: E402

PDS_SPAN = 0.95228 - 0.50020
REF = {"h1-test": 0.036, "orion-hct116": 0.264, "orion-hek293t": 0.441}


def main() -> None:
    P = json.loads(open(os.path.join(WORK, "plolo169.json")).read())
    R = {(r["fold"], r["arm"], r["q"], r["space"]): r for r in P}

    def chal(f, q):
        return float(np.mean([R[(f, "a23", q, f"chal{c}/delivered")]["pds"] for c in "ABC"]))

    print("published shares (157 §1 / 159 §2): "
          + ", ".join(f"{k} {v * 100:.1f}%" for k, v in REF.items()))
    print(f"\n{'fold':<16}{'targets':>9}{'scored':>9}{'blob share':>12}"
          f"{'d pds q=1.15':>14}{'scaled':>10}")
    rows = []
    for held in FOLDS:
        fold = build_fold(held)
        T = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        mu = T.mean(0)
        blob = float(T.shape[0] * (mu @ mu) / (T * T).sum())
        d = chal(held, 1.15) - chal(held, 1.00)
        rows.append((held, blob, d))
        print(f"{held:<16}{T.shape[0]:>9}{T.shape[1]:>9}{blob * 100:>11.1f}%"
              f"{d:>+14.6f}{d / PDS_SPAN:>+10.6f}")

    b = np.array([r[1] for r in rows])
    d = np.array([r[2] for r in rows])
    c = float(np.corrcoef(b, d)[0, 1])
    print(f"\ncorrelation(blob share, d pds) over the four folds: {c:+.3f}  (n = 4)")
    lo = [r for r in rows if r[1] < 0.15]
    print(f"folds with a blob share BELOW 15%: "
          + ", ".join(f"{r[0]} ({r[1]*100:.1f}%, d pds {r[2]:+.6f})" for r in lo))
    print("\nIf the gain were blob alignment, the low-blob folds would show little or none.")


if __name__ == "__main__":
    sys.exit(main())
