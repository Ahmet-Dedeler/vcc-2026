"""220 — HOW GOOD WOULD A CONFIDENCE SIGNAL HAVE TO BE?

The transfer curve (`curve220.py`) says `reach` responds to front sign accuracy at
`d reach / d acc@10 = 0.143`, and rule 13's profile-preserving gate needs **+0.00783 raw
h1-test `reach`**.  So the family is not capped by the 0.90 bar; it is capped by the quality
of the confidence signal available.  This turns that into a specification.

A NOISY ORACLE: a per-(target, gene) bit `o` that equals `1[our sign is right]` with
probability `a` and is flipped otherwise.  `a = 0.5` is pure noise, `a = 1` is the oracle.
The key is `|H| * exp(c * o)`.  Sweeping `a` prices exactly what a confidence head would have
to deliver, in a currency (`a`) that any future head can be measured in.  20 seeds per cell.

Reported alongside: the Spearman of `o` against the truth, so the requirement can also be
read as a rank-correlation target.

    uv run python noisy220.py
"""
from __future__ import annotations

import json

import numpy as np

import common220 as C

GATE_RAW = 0.0019 * 6 / (1.290 * 1.12927)          # +0.0019 overall in raw h1-test reach


def main() -> None:
    T = C.Instrument(C.DE_RC, "rc_callset")
    old = C.Instrument(C.DE_A15, "a15")
    C.check_old(old)
    base = float(T.reach_vec(T.key_absH()).mean())
    truth = T.correct.astype(np.float64)
    pool = T.pool
    print(f"\n  incumbent `reach(|H|)` {base:.5f}; oracle {T.reach(T.key_oracle()):.5f}")
    print(f"  rule 13's profile-preserving gate +0.0019 overall = "
          f"**+{GATE_RAW:.5f} raw h1-test `reach`** "
          f"(x 1.290 board factor, /span 1.12927, /6)\n")

    print("=" * 100)
    print("A NOISY CONFIDENCE ORACLE: `o = 1[correct]` seen with per-pair accuracy `a`")
    print("=" * 100)
    print(f"{'a':>6}{'rho(o, truth) in pool':>24}" +
          "".join(f"{'c='+str(c):>12}" for c in (0.2, 0.5, 1.0, 3.0)))
    rows = []
    for a in (0.50, 0.55, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00):
        vals = {c: [] for c in (0.2, 0.5, 1.0, 3.0)}
        rhos = []
        for sd in range(20 if a < 1.0 else 1):
            rng = np.random.default_rng(1000 + sd)
            flip = rng.random(truth.shape) > a
            o = np.where(flip, 1.0 - truth, truth)
            rhos.append(float(np.corrcoef(o[pool], truth[pool])[0, 1]))
            for c in vals:
                k = np.abs(T.H) * np.exp(c * o)
                vals[c].append(float(T.reach_vec(k).mean()) - base)
        r = dict(a=a, rho=float(np.mean(rhos)),
                 **{f"c{c}": float(np.mean(vals[c])) for c in vals},
                 **{f"sd{c}": float(np.std(vals[c], ddof=1)) if len(vals[c]) > 1 else 0.0
                    for c in vals})
        rows.append(r)
        print(f"{a:>6.2f}{np.mean(rhos):>24.4f}" +
              "".join(f"{r[f'c{c}']:>+12.5f}" for c in (0.2, 0.5, 1.0, 3.0)))

    best = [(r["a"], max(r[f"c{c}"] for c in (0.2, 0.5, 1.0, 3.0))) for r in rows]
    cross = None
    for (a0, v0), (a1, v1) in zip(best, best[1:]):
        if v0 < GATE_RAW <= v1:
            cross = a0 + (GATE_RAW - v0) * (a1 - a0) / (v1 - v0)
            break
    print(f"\n  The gate (+{GATE_RAW:.5f} raw) is cleared at a per-pair confidence accuracy "
          f"of a ~ {cross:.3f}" if cross else
          f"\n  The gate (+{GATE_RAW:.5f} raw) is never cleared on this grid")
    if cross:
        # what rank correlation with the truth does that correspond to?
        j = min(range(len(rows)), key=lambda i: abs(rows[i]["a"] - cross))
        print(f"  i.e. a per-(target, gene) confidence whose correlation with "
              f"`1[our sign is right]` is ~ {rows[j]['rho']:.3f} inside the reference DE "
              f"pool.")
    print("\n  FOR SCALE, the same statistic for every key ever measured on this surface:")
    ref = {"|H| itself (172)": None, "189 theta (189 section 6.1)": "theta",
           "control expression (189 F2)": "expr"}
    absH = np.abs(T.H)
    r_absH = float(np.corrcoef(absH[pool], truth[pool])[0, 1])
    th = np.broadcast_to(T.theta[None, :], T.H.shape)
    r_th = float(np.corrcoef(th[pool], truth[pool])[0, 1])
    ex = np.broadcast_to(T.b_ctrl[None, :], T.H.shape)
    r_ex = float(np.corrcoef(ex[pool], truth[pool])[0, 1])
    def auc(sig):
        """Mann-Whitney AUC of `sig` for predicting `1[our sign is right]` in the pool."""
        s = np.asarray(sig)[pool]
        y = truth[pool] > 0.5
        r = np.empty(s.size)
        o = np.argsort(s, kind="stable")
        sv = s[o]
        i = 0
        while i < s.size:
            j = i
            while j + 1 < s.size and sv[j + 1] == sv[i]:
                j += 1
            r[o[i:j + 1]] = 0.5 * (i + 1 + j + 1)
            i = j + 1
        n1, n0 = int(y.sum()), int((~y).sum())
        return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

    a_absH, a_th, a_ex = auc(absH), auc(th), auc(ex)
    print(f"    |H|         corr {r_absH:+.4f}   AUC {a_absH:.4f}")
    print(f"    189's theta corr {r_th:+.4f}   AUC {a_th:.4f}")
    print(f"    expression  corr {r_ex:+.4f}   AUC {a_ex:.4f}")
    if cross:
        print(f"    REQUIRED    corr {rows[j]['rho']:+.4f}   AUC {cross:.4f}   "
              "(a symmetric binary bit of accuracy `a` has AUC = a)")
    print("  (|H| is the incumbent key and it is already the strongest of the three; the "
          "requirement above is on a NEW signal, orthogonal to |H|.)")

    json.dump(dict(gate_raw=GATE_RAW, base=base, rows=rows, cross=cross,
                   corr=dict(absH=r_absH, theta=r_th, expr=r_ex),
                   auc=dict(absH=a_absH, theta=a_th, expr=a_ex)),
              open(C.WORK / "noisy220.json", "w"), indent=1, default=float)
    print(f"\nwrote {C.WORK/'noisy220.json'}")


if __name__ == "__main__":
    main()
