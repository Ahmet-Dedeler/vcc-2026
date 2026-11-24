"""521 — does the bottom tail of the strength ranking stay pure off H1? The rank policy that survives the selection null
(bottom few % by strength -> ~10 calls) lives entirely on the purity of that tail.

(1) Tail purity per surface: of the bottom q % of targets by strength, the fraction whose truth n_real sits in the truth's
    own bottom 35 % (H1's '< 20' base rate is 41 %; ranks with ties split), and the fraction in the truth's TOP half (the
    mistakes that cost fid's coverage). H1 (f 1) and the three cancer folds (480's gates, sources exclude the held line).
(2) Break-even purity on H1: per-target utilities of shrinking to K (f 1.0 and cmp, variant b and a), split by whether
    the target is in H1's bottom 35 % of n_real; pi* = the negligible fraction at which a shrink breaks even.
(3) Transfer: the value of 'bottom q % -> K' at a surface's measured purity, using H1's conditional utilities.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/purity521.py
"""
from __future__ import annotations

import json
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import lib521 as L  # noqa: E402
from policy521 import Table, per_target_contrib  # noqa: E402
from pred521 import BANKS, cd4_feats, fold_surface, meanrank, strength  # noqa: E402

LOWQ = 0.35


def trank(n):
    """truth percentile rank in [0, 1] with ties split."""
    from scipy.stats import rankdata
    return (rankdata(n, method="average") - 1) / (n.size - 1)


def tail(z, n, qs=(5, 10, 15, 20)):
    ok = np.isfinite(z)
    z, n = z[ok], n[ok]
    tr = trank(n.astype(float))
    out = {}
    for q in qs:
        k = max(1, int(round(q / 100 * z.size)))
        sel = np.argsort(z, kind="stable")[:k]
        out[q] = dict(k=k, low=float(np.mean(tr[sel] <= LOWQ)), tophalf=float(np.mean(tr[sel] > 0.5)),
                      med_n=float(np.median(n[sel])), med_pct=float(np.median(tr[sel])))
    return out


def main():
    res = {"tail": {}, "breakeven": {}, "transfer": {}}
    P = np.load(L.WORK / "pred521_h1.npz")
    n1 = P["nreal"]
    for fn, key in (("S4", "S4__S3_p_vip_thr4"), ("S4avg", "S4_thr-avg"), ("K562", "k562_t4")):
        res["tail"][f"H1|{fn}"] = tail(P[key], n1)
    for fold in ("k562", "hct116", "hek293t"):
        Tf, nf = fold_surface(fold)
        G = {b: strength(BANKS[b], Tf, 4.0) for b in ("k562", "hct116", "hek293t", "vip", "kolf2")}
        G["cd4rob"], G["cd4nde"] = cd4_feats(Tf)
        prod = [b for b in ("k562", "hct116", "hek293t", "vip") if b != fold and not (fold == "k562" and b == "vip")]
        allsrc = prod + ["kolf2", "cd4rob", "cd4nde"]
        res["tail"][f"{fold}|prod"] = tail(meanrank([G[b] for b in prod]), nf)
        res["tail"][f"{fold}|all"] = tail(meanrank([G[b] for b in allsrc]), nf)
    print(f"(1) tail purity: fraction of the bottom q % by strength in the truth's bottom {int(LOWQ * 100)} % "
          f"(base rate {LOWQ:.2f}) | in the truth's top half | median n_real")
    for k, v in res["tail"].items():
        print(f"  {k:16s} " + "  ".join(f"q{q}: {d['low']:.2f}|{d['tophalf']:.2f}|{d['med_n']:.0f} (k {d['k']})" for q, d in v.items()))
    # (2) break-even purity from H1's per-target utilities
    T = Table(0)
    low = trank(n1.astype(float)) <= LOWQ
    print(f"\n(2) H1 per-target utility of shrinking one target to K (board overall units, per target), by class "
          f"(negligible = H1 bottom {int(LOWQ * 100)} % of n_real, {int(low.sum())} targets)")
    for reg in ("1.0", "cmp"):
        for var in ("b", "a"):
            R0 = T.rows[f"{reg}|base"]
            for K in (0, 10, 25, 50, 100, 400):
                u = per_target_contrib(T.assemble(reg, "headcpm", var, [K] * 85), R0)
                ul, uh = float(u[low].mean()), float(u[~low].mean())
                pi = uh / (uh - ul) if ul > 0 > uh else (0.0 if ul > 0 and uh >= 0 else 1.0)
                res["breakeven"][f"{reg}|{var}|K{K}"] = dict(u_low=ul, u_high=uh, pi_star=pi,
                                                             u_low_sd=float(u[low].std(ddof=1)), u_high_sd=float(u[~low].std(ddof=1)))
                print(f"  {reg:4s} {var} K {K:4d}: E[u|negligible] {ul:+.6f} (sd {u[low].std(ddof=1):.6f})  "
                      f"E[u|not] {uh:+.6f} (sd {u[~low].std(ddof=1):.6f})  break-even purity {pi:.2f}")
    # (3) transfer: value of 'bottom q % -> K' on a 300-target context at each surface's measured purity
    print("\n(3) expected board-context value of 'bottom q % -> K' = q x [pi E(u|negl) + (1 - pi) E(u|not)] x 85 "
          "(H1 conditional utilities, per-target units scaled to the member means)")
    for reg in ("1.0", "cmp"):
        for var in ("b",):
            for K in (10, 25):
                be = res["breakeven"][f"{reg}|{var}|K{K}"]
                for surf, v in res["tail"].items():
                    for q in (5, 10):
                        pi = v[q]["low"]
                        val = (q / 100) * 85 * (pi * be["u_low"] + (1 - pi) * be["u_high"])
                        res["transfer"][f"{reg}|{var}|K{K}|{surf}|q{q}"] = dict(pi=pi, value=val)
                print(f"  {reg} {var} K {K}: " + "  ".join(
                    f"{s}:q5 {res['transfer'][f'{reg}|{var}|K{K}|{s}|q5']['value']:+.4f}/q10 {res['transfer'][f'{reg}|{var}|K{K}|{s}|q10']['value']:+.4f}"
                    for s in res["tail"]))
    (L.WORK / "purity521.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
