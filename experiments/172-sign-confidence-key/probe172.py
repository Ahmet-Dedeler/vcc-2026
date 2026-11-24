"""172 — probe: is the UN-WHITENED consensus a usable stand-in for the delivered `|H|`?

This decides the whole design.  The learned key is fitted on LODO folds, where no leak-free
whitening operator exists, so its inputs are un-whitened.  If the un-whitened consensus is a
much weaker key than the whitened `|H|` that the brief benchmarks against, then every learned
key would be handicapped by a variable that has nothing to do with learning.

Also measures how often whitening FLIPS a sign, because the label on h1-test is the
correctness of the *emitted* (whitened) sign while the label in a LODO fold can only be the
correctness of the un-whitened consensus sign.  If whitening rarely flips, the two labels are
the same problem; if it flips often, they are not, and that has to be said out loud.
"""
from __future__ import annotations

import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from feats172 import Recipient, build            # noqa: E402
from harness172 import H1Test, check             # noqa: E402

SOURCES = ("k562_gwps", "orion_hek293t")         # 095's h1-test artifact, INBOX-4 §"Leak check"


def main() -> None:
    T = H1Test()
    print("known-answer rows:")
    check(T)
    live = T.gate.copy()
    print(f"\nrecipient h1-test: {len(T.perts)} targets x {len(T.genes)} DE-axis genes, "
          f"{live.sum()} gated (min_cpm 5.0)")
    rec = Recipient("h1test", T.genes, T.b_ctrl, T.perts, live)
    B = build(rec, SOURCES)
    cons, geo, agree, F = B["cons"], B["geo"], B["agree"], B["F"]

    absH = np.abs(T.H)
    print(f"\ntargets carried by at least one source donor: {int(B['carried_any'].sum())}/85")

    # ---- 1. is |cons| a stand-in for |H| as an ORDERING KEY? ---------------------------
    print(f"\n{'key':<44}{'reach':>10}{'d1 acc':>9}   deciles")
    rows = [("production (our emitted p-values)", None),
            ("|H|  (whitened, delivered -- the benchmark)", absH),
            ("|cons| (un-whitened consensus)", np.abs(cons)),
            ("geo x concord (013/095's key)", agree),
            ("geo (magnitude only)", geo),
            ("|cons| / donor-panel sd", np.abs(cons) / np.maximum(
                np.asarray(F["log_sd"])[0] * 0 + 1.0, 1.0)),
            ("ORACLE (our emitted signs)", T.key_oracle())]
    for name, k in rows:
        r = T.reach(k)
        d = T.deciles(k if k is not None else _prodkey(T))
        print(f"{name:<44}{r:>10.5f}{d[0]:>9.3f}   " + " ".join(f"{x:.3f}" for x in d))

    # ---- 2. how close are the two magnitudes, and do the SIGNS agree? ------------------
    pool = T.pool
    rs, ss = [], []
    for p in range(len(T.perts)):
        i = np.flatnonzero(pool[p])
        if i.size < 5:
            continue
        rs.append(spearmanr(absH[p, i], np.abs(cons[p, i])).statistic)
        ss.append(float((np.sign(T.H[p, i]) == np.sign(cons[p, i])).mean()))
    print(f"\nwithin-target Spearman(|H|, |cons|) over the DE pool: median {np.median(rs):.4f}")
    print(f"sign(H) == sign(cons) over the DE pool:               median {np.median(ss):.4f}")
    m = pool & (cons != 0)
    print(f"sign agreement over pool pairs with a non-zero consensus: "
          f"{float((np.sign(T.H[m]) == np.sign(cons[m])).mean()):.4f}  ({int(m.sum())} pairs)")

    # ---- 3. sign accuracy of each sign source against the reference --------------------
    real = T.S.real_lfc
    for nm, v in (("emitted pred_lfc (what reach scores)", T.S.pred_lfc),
                  ("recovered H", T.H), ("un-whitened consensus", cons)):
        acc = float((np.sign(v[pool]) == np.sign(real[pool])).mean())
        print(f"pooled sign accuracy over the DE pool, {nm:<38} {acc:.4f}")

    # ---- 4. per-feature univariate reach, so the model has a reference point -----------
    print(f"\n{'single feature as the key':<20}{'reach':>10}{'d1':>8}{'rho(sign-correct)':>20}")
    ok = T.correct
    for k in F:
        v = np.asarray(F[k], dtype=np.float64)
        r = T.reach(v)
        d = T.deciles(v)[0]
        rho = spearmanr(v[pool], ok[pool].astype(float)).statistic
        print(f"{k:<20}{r:>10.5f}{d:>8.3f}{rho:>20.4f}")
    np.savez_compressed(os.path.join(HERE, "work", "h1test-feats.npz"),
                        cons=cons, geo=geo, agree=agree, live=live,
                        **{f"F_{k}": np.asarray(F[k], dtype=np.float32) for k in F})
    print(f"\nwrote work/h1test-feats.npz")


def _prodkey(T):
    q = np.nan_to_num(T.S.pred_padj, nan=1.0)
    return np.where(T.S.pred_sig, 1e12 - q, -q)


if __name__ == "__main__":
    main()
