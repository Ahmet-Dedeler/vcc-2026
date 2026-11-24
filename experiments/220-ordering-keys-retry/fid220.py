"""220 — the one loose end 190 left: its `fid` row, null-tested.

190 section 5 reported `ct-prop-k25`'s precision at `n_pred = 50` as 0.5332 against `|H|`'s
0.5254, a `+0.00555` Delta`fid` at one grid cell, called it one sigma, and said explicitly
that *"`fid` was not itself null-tested -- stated as a gap, not smuggled in as a finding"*.
189 section 4 tested the same statistic for its own key and found the gain reproduced by a
SHUFFLED prior (z = -0.04).  This closes 190's gap the same way: 40 shuffles of the network's
signs, best cell re-selected inside the null.

Precision at fixed `n` is the transferable quantity (rule 9: absolute `fid` is unreadable on
this surface because `n_real` is 222 here against ~90-110 on the board).

    uv run python fid220.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

import common220 as C

sys.path.insert(0, str(C.ROOT / "experiments/190-regulon-ordering"))
NS = (10, 25, 50, 100, 200, 400)


def precision_at(T: C.Instrument, key, ns=NS) -> dict:
    """Top-`n` sign accuracy over adjudicable genes, per target, then averaged (macro)."""
    S = T.S
    acc = {n: [] for n in ns}
    for p in range(len(T.perts)):
        idx = np.flatnonzero(S.adj[p])
        if idx.size == 0:
            continue
        k = np.asarray(key[p], dtype=np.float64)[idx]
        ok = (np.sign(S.pred_lfc[p, idx]) == np.sign(S.real_lfc[p, idx])).astype(float)
        ok = ok[np.lexsort((T.genes[idx], -k))]
        cum = np.cumsum(ok)
        for n in ns:
            m = min(n, ok.size)
            acc[n].append(float(cum[m - 1] / m))
    return {n: np.asarray(acc[n]) for n in ns}


def main() -> None:
    T = C.Instrument(C.DE_RC, "rc_callset")
    old = C.Instrument(C.DE_A15, "a15")
    C.check_old(old)
    import nets190

    base = precision_at(T, T.key_absH())
    N, _ = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25)
    R = C.agreement_R(N, T.H)
    keys = {"189 lam 0.50": C.key189(T.H, T.theta, 0.5),
            "190 lam -0.25": np.abs(T.H) * np.exp(np.clip(-0.25 * R, -50, 50))}

    print("\nPRECISION AT FIXED n (macro over targets), rc_callset instrument")
    print(f"{'key':<18}" + "".join(f"{'n='+str(n):>11}" for n in NS))
    print(f"{'|H| (incumbent)':<18}" +
          "".join(f"{base[n].mean():>11.4f}" for n in NS))
    real = {}
    for nm, k in keys.items():
        pr = precision_at(T, k)
        real[nm] = {n: float(pr[n].mean() - base[n].mean()) for n in NS}
        print(f"{nm:<18}" + "".join(f"{pr[n].mean():>11.4f}" for n in NS))
    print(f"\n{'delta vs |H|':<18}" + "".join(f"{'n='+str(n):>11}" for n in NS))
    for nm in keys:
        print(f"{nm:<18}" + "".join(f"{real[nm][n]:>+11.5f}" for n in NS))

    print("\nNULL (40 draws).  189: theta permuted across the gene axis.  "
          "190: the network's signs permuted, topology held.")
    out = {"base": {str(n): float(base[n].mean()) for n in NS}, "real": {}, "null": {}}
    for nm in keys:
        draws = {n: [] for n in NS}
        for sd in range(40):
            rng = np.random.default_rng(9000 + sd)
            if nm.startswith("189"):
                k = C.key189(T.H, rng.permutation(T.theta), 0.5)
            else:
                Ns, _ = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25,
                                        rng=rng, shuffle="signs")
                k = np.abs(T.H) * np.exp(np.clip(-0.25 * C.agreement_R(Ns, T.H), -50, 50))
            pr = precision_at(T, k)
            for n in NS:
                draws[n].append(float(pr[n].mean() - base[n].mean()))
        print(f"  {nm}")
        print(f"{'':<6}{'n':>6}{'real':>11}{'null mean':>12}{'null sd':>10}{'z':>8}"
              f"{'P(null>=real)':>15}")
        for n in NS:
            d = np.asarray(draws[n])
            z = (real[nm][n] - d.mean()) / max(d.std(ddof=1), 1e-12)
            print(f"{'':<6}{n:>6}{real[nm][n]:>+11.5f}{d.mean():>+12.5f}"
                  f"{d.std(ddof=1):>10.5f}{z:>8.2f}{(d >= real[nm][n]).mean():>15.3f}")
            out["null"].setdefault(nm, {})[str(n)] = dict(
                real=real[nm][n], null=float(d.mean()), sd=float(d.std(ddof=1)),
                z=float(z), p=float((d >= real[nm][n]).mean()))
        # best cell re-selected inside the null (190's own discipline)
        best_real = max(real[nm].values())
        best_null = np.asarray([max(draws[n][s] for n in NS) for s in range(40)])
        print(f"       best cell over the {len(NS)} sizes, re-selected inside the null: "
              f"real {best_real:+.5f}   null {best_null.mean():+.5f} +- "
              f"{best_null.std(ddof=1):.5f}   P(null >= real) "
              f"{(best_null >= best_real).mean():.3f}")
        out["real"][nm] = {str(n): real[nm][n] for n in NS}
        out.setdefault("best", {})[nm] = dict(real=float(best_real),
                                              null=float(best_null.mean()),
                                              sd=float(best_null.std(ddof=1)),
                                              p=float((best_null >= best_real).mean()))

    json.dump(out, open(C.WORK / "fid220.json", "w"), indent=1, default=float)
    print(f"\nwrote {C.WORK/'fid220.json'}")


if __name__ == "__main__":
    main()
