"""185 — where the reach gain lives: per-perturbation paired deltas by n_real stratum
(171 §7a(i) / 181 §15's kill-shot test applied to this arm), the full P(k) curve, and the
twin spreads that set the resolution of every number in this experiment."""
from __future__ import annotations

import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, f"{ROOT}/experiments/153-de-bench")
sys.path.insert(0, HERE)
from debench import Surface  # noqa: E402
from six185 import PowerSurface, walk  # noqa: E402

P = os.path.join(WORK, "de185.%s.npz")


def per_pert_reach(s, key=None):
    out = np.full(len(s.perts), np.nan)
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        if key is None:
            order = np.lexsort((s.genes[pool], -np.abs(s.pred_lfc[p, pool]), s.pred_p[p, pool],
                               s.pred_padj[p, pool], ~s.pred_sig[p, pool]))
        else:
            order = np.lexsort((s.genes[pool], -key[p][pool]))
        pool = pool[order]
        pool = pool[s.adj[p, pool]]
        if pool.size == 0:
            out[p] = 0.0
            continue
        ok = (np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool]))
        pur = np.cumsum(ok) / np.arange(1, pool.size + 1)
        deep = np.flatnonzero(pur >= 0.9)
        out[p] = (int(deep[-1] + 1) if deep.size else 0) / int(s.real_sig[p].sum())
    return out


def main():
    arms = [a for a in ("a15", "a15_b", "rc_hist_a", "rc_hist_a_b", "rc_hist_s", "rc_hybrid",
                        "rc_lfc", "rc_perm", "rc_k2000", "a21_like") if os.path.exists(P % a)]
    S = {a: Surface(P % a) for a in arms}
    absH = np.load(P % "a15", allow_pickle=True)["absH"].astype(np.float64)
    nr = np.array([int(S["a15"].real_sig[p].sum()) for p in range(len(S["a15"].perts))])
    base = 0.5 * (per_pert_reach(S["a15"]) + per_pert_reach(S["a15_b"]))
    cf = per_pert_reach(S["a15"], key=absH)
    strata = [("n_real <= 13", nr <= 13), ("14 - 800", (nr > 13) & (nr <= 800)),
              ("> 800", nr > 800), ("n_real >= 20", nr >= 20), ("n_real >= 100", nr >= 100),
              ("ALL", nr > 0)]
    print("PAIRED delta reach vs a15 (twin mean) by n_real stratum; 'share' = stratum's share of the panel sum")
    print(f"{'arm':<12}" + "".join(f"{lab:>16}" for lab, _ in strata))
    rows = {"|H| counterfactual": cf}
    for a in arms:
        if a in ("a15", "a15_b"):
            continue
        rows[a] = per_pert_reach(S[a])
    for name, r in rows.items():
        d = r - base
        tot = np.nansum(d)
        cells = []
        for lab, m in strata:
            dm = d[m & np.isfinite(d)]
            share = 100 * dm.sum() / tot if tot != 0 else 0
            cells.append(f"{dm.mean():+.4f} ({share:4.0f}%)" if lab != "ALL" else f"{dm.mean():+.4f} n={m.sum()}")
        print(f"{name:<12}" + "".join(f"{c:>16}" for c in cells))
    # paired t and the most influential perturbation
    for a in ("rc_hist_a", "rc_hybrid"):
        if a not in rows:
            continue
        d = rows[a] - base
        d = d[np.isfinite(d)]
        t = d.mean() / (d.std(ddof=1) / np.sqrt(d.size))
        j = np.argmax(np.abs(d))
        print(f"  {a}: paired t {t:+.2f} ({(d>0).sum()} up / {(d<0).sum()} down / {(d==0).sum()} same); "
              f"top-5 |d| share of the sum {100*np.sort(np.abs(d))[-5:].sum()/np.abs(d).sum():.0f}%; "
              f"most influential n_real {nr[np.isfinite(rows[a]-base)][j]} d {d[j]:+.3f}")

    print("\nP(k) curve (mean over perturbations with >= k adjudicable pool genes)")
    ks = (1, 2, 3, 5, 10, 20, 50, 100)
    print(f"{'arm':<12}" + "".join(f"{'P('+str(k)+')':>8}" for k in ks) + f"{'reach':>9}")
    for a in arms:
        r, pk, _ = walk(S[a])
        print(f"{a:<12}" + "".join(f"{pk[k]:>8.4f}" for k in ks) + f"{r:>9.5f}")
    r, pk, _ = walk(S["a15"], key=absH)
    print(f"{'|H| key':<12}" + "".join(f"{pk[k]:>8.4f}" for k in ks) + f"{r:>9.5f}")

    print("\nTWIN SPREADS (same policy, re-seeded), the resolution floor of this experiment")
    for a, b in (("a15", "a15_b"), ("rc_hist_a", "rc_hist_a_b")):
        if a in S and b in S:
            ra, rb = S[a].reach(), S[b].reach()
            print(f"  {a:<10} {ra:.5f}  {b:<12} {rb:.5f}   |diff| {abs(ra-rb):.5f}")

    print("\nREACH AT PROPORTIONAL POWER f (rc_hist_a - a15, both single seed 90000)")
    for f in (0.5, 1.0, 1.5, 2.0, 3.5, 5.0):
        ra = PowerSurface(P % "a15", f).reach()
        rc = PowerSurface(P % "rc_hist_a", f).reach()
        rh = PowerSurface(P % "rc_hybrid", f).reach() if "rc_hybrid" in S else float("nan")
        print(f"  f={f:<4} a15 {ra:.5f}  rc_hist_a {rc:.5f} ({rc-ra:+.5f})  rc_hybrid {rh:.5f} ({rh-ra:+.5f})")


if __name__ == "__main__":
    main()
