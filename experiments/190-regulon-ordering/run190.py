"""190 — the measurement.  Coverage, ceilings, key sweep, all four falsifiers, `fid`.

Known-answer rows, checked before anything prints (rule 10 + its corollary):
  * 172's four reach rows and two decile rows, via harness172.check
  * reach(midrank(|H|)) == reach(|H|)                        (the rank transform is order-safe)
  * 174's co-expression rows: |r| median 0.0256 at pc0,
    sign accuracy 0.5158 for -sign(r) and 0.4842 for +sign(r) on the reference pool
    -- this row already caught a flipped sign convention in the first draft of this experiment.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.join(ROOT, "experiments", "190-regulon-ordering")
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(ROOT, "experiments", "172-sign-confidence-key"))
sys.path.insert(0, HERE)

import nets190  # noqa: E402
from key190 import agreement, key_mult, key_rank, midrank_within  # noqa: E402

BIG = 1e9
LAM = (-8.0, -4.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
# board anchors, 094 §1 (debench.ANCHOR): overall = raw / span / 6
SPAN = {"reach": 0.96414 - 0.07902, "fid": 0.80964 - 0.51239, "jac": 0.39727 - 0.03038}


def to_overall(member: str, d_raw: float) -> float:
    return d_raw / SPAN[member] / 6.0


def build_all(T, H, rng=None, shuffle=None):
    """The network library, in one place so the falsifiers rebuild exactly the same objects."""
    tg, gn = T.perts, T.genes
    out = {}
    out["ct-direct"] = nets190.ct_direct(tg, gn, H, rng=rng, shuffle=shuffle)[0]
    out["ct-prop-self"] = nets190.ct_prop(tg, gn, H, rng=rng, shuffle=shuffle,
                                          self_only=True)[0]
    for k in (5, 25, None):
        out[f"ct-prop-k{k}"] = nets190.ct_prop(tg, gn, H, topk=k, rng=rng, shuffle=shuffle)[0]
    out["coexpr-pc0"] = nets190.coexpr(tg, gn, H, pc=0, rng=rng, shuffle=shuffle)[0]
    out["coexpr-pc2"] = nets190.coexpr(tg, gn, H, pc=2, rng=rng, shuffle=shuffle)[0]
    out["coexpr-pc0-top500"] = nets190.coexpr(tg, gn, H, pc=0, topk=500, rng=rng,
                                              shuffle=shuffle)[0]
    return out


def main() -> None:
    import harness172

    T = harness172.H1Test()
    print("=" * 100)
    print("PHASE 0 — KNOWN-ANSWER ROWS.  Any failure aborts (AGENTS.md §11.1 rule 10).")
    print("=" * 100)
    harness172.check(T)
    H = T.H.astype(np.float64)
    pool = T.pool
    sT = np.sign(T.S.real_lfc)
    sH = np.sign(H)
    rH = midrank_within(np.abs(H))
    base = T.reach(np.abs(H))
    prod = T.reach()
    assert abs(base - 0.23533) < 5e-6 and abs(prod - 0.18283) < 5e-6
    assert abs(T.reach(rH) - base) < 1e-9
    print(f"  OK  reach(midrank(|H|)) == reach(|H|)     {T.reach(rH):.5f}")
    Cz = np.load(os.path.join(WORK, "coex190.npz"), allow_pickle=True)
    C0 = Cz["C_pc0"].astype(np.float64)
    kar = {"coex_median_abs_r_pc0": float(np.median(np.abs(C0))),
           "coex_acc_minus_sign_r": float(np.mean(-np.sign(C0)[pool] == sT[pool])),
           "coex_acc_plus_sign_r": float(np.mean(np.sign(C0)[pool] == sT[pool]))}
    exp174 = {"coex_median_abs_r_pc0": 0.0256, "coex_acc_minus_sign_r": 0.5158,
              "coex_acc_plus_sign_r": 0.4842}
    for k, v in kar.items():
        ok = abs(v - exp174[k]) <= 0.0005
        print(f"  {'OK ' if ok else 'FAIL'} 174 {k:<28} {v:.4f}   expected {exp174[k]:.4f}")
        if not ok:
            raise SystemExit(f"KNOWN-ANSWER ROW FAILED: {k}")
    print(f"\n  surface: 85 targets x {len(T.genes):,} genes | pool {int(pool.sum()):,} pairs | "
          f"|H| sign accuracy on the pool {float(np.mean(sH[pool] == sT[pool])):.4f}")
    out: dict = {"base": {"absH": base, "production": prod, "oracle": T.reach(T.key_oracle()),
                          "H_oracle_172": 0.61862}}

    nets = build_all(T, H)

    # ---------------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("PHASE 1/3b — COVERAGE, and whether the network's sign beats ours where it speaks")
    print("=" * 100)
    print(f"{'network':22} {'pool pairs':>11} {'% pool':>8} {'targets':>8} "
          f"{'net sign':>9} {'|H| there':>10} {'net right':>10}")
    print(f"{'':22} {'':>11} {'':>8} {'covered':>8} {'':>9} {'':>10} {'& H wrong':>10}")
    cov = {}
    for name, N in nets.items():
        m = pool & (N != 0) & (H != 0)
        n = int(m.sum())
        ntg = int(((N != 0).any(1)).sum())
        if n == 0:
            print(f"{name:22} {0:11,}")
            cov[name] = dict(n=0)
            continue
        a_net = float(np.mean(np.sign(N[m]) == sT[m]))
        a_H = float(np.mean(sH[m] == sT[m]))
        rescue = float(np.mean((np.sign(N[m]) == sT[m]) & (sH[m] != sT[m])))
        cov[name] = dict(n=n, frac_pool=n / int(pool.sum()), targets=ntg,
                         net_acc=a_net, H_acc_there=a_H, net_right_H_wrong=rescue)
        print(f"{name:22} {n:11,} {100*n/pool.sum():7.2f}% {ntg:8d} {a_net:9.4f} "
              f"{a_H:10.4f} {rescue:10.4f}")
    out["coverage"] = cov

    # ---------------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("PHASE 2a — CEILING: an oracle confined to the pairs each network covers.")
    print("           No prior built on that network can beat this number.")
    print("=" * 100)
    correct = (np.sign(T.S.pred_lfc) == np.sign(T.S.real_lfc))
    ceil = {}
    for name, N in nets.items():
        c = N != 0
        k = np.where(c & correct, BIG + np.abs(H),
                     np.where(c & ~correct, -BIG + np.abs(H), np.abs(H)))
        r = T.reach(k)
        ceil[name] = r
        print(f"   {name:22} ceiling {r:.5f}  raw {r-base:+.5f}  "
              f"overall {to_overall('reach', r-base):+.5f}")
    out["ceiling"] = ceil

    # ---------------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("PHASE 2b — THE KEY SWEEP.  reach(|H| * exp(lambda * R)); lambda=0 IS |H| = 0.23533")
    print("=" * 100)
    print(f"   {'network':22}" + "".join(f"{l:>9.2f}" for l in LAM) + f"{'best':>9}{'raw d':>9}")
    sweeps, best = {}, {}
    for name, N in nets.items():
        R = agreement(N, H)
        row = [T.reach(key_mult(H, R, l)) for l in LAM]
        rowr = [T.reach(key_rank(rH, R, l)) for l in LAM]
        sweeps[name] = {"mult": row, "rank": rowr}
        j = int(np.argmax(row))
        best[name] = {"lambda": LAM[j], "reach": row[j], "d_raw": row[j] - base,
                      "d_overall": to_overall("reach", row[j] - base)}
        print(f"   {name:22}" + "".join(f"{v:9.5f}" for v in row)
              + f"{row[j]:9.5f}{row[j]-base:+9.5f}")
    print(f"\n   rank-additive: reach(midrank(|H|) + lambda * R)")
    print(f"   {'network':22}" + "".join(f"{l:>9.2f}" for l in LAM))
    for name in nets:
        print(f"   {name:22}" + "".join(f"{v:9.5f}" for v in sweeps[name]["rank"]))
    out["sweep"] = {"lambdas": list(LAM), "reach": sweeps, "best": best}

    # ---------------------------------------------------------------------------------
    print("\n" + "=" * 100)
    print("PHASE 3a — FALSIFIER 1: the SHUFFLED-REGULON control tydymy did not report.")
    print("           'targets' redraws each TF's target list (degree + sign multiset kept)")
    print("           'signs'   permutes mor globally (topology kept)")
    print("           'tfs'     reassigns whole regulons to other targets")
    print("           Each: 20 seeds, best-lambda reach recomputed per seed.")
    print("=" * 100)
    NSEED = 20
    SHUF = ("targets", "signs", "tfs")
    raw = {name: {sh: [] for sh in SHUF} for name in nets}
    for sh in SHUF:
        for s in range(NSEED):
            rng = np.random.default_rng(1000 + s)
            ns = build_all(T, H, rng=rng, shuffle=sh)
            for name, Ns in ns.items():
                Rs = agreement(Ns, H)
                raw[name][sh].append(max(T.reach(key_mult(H, Rs, l)) for l in LAM))
        print(f"   ... shuffle={sh} done", flush=True)
    nulls = {}
    print(f"\n   {'network':22} {'real best':>10}  " +
          "".join(f"{'shuffle=' + s:>26}" for s in SHUF))
    for name in nets:
        line = f"   {name:22} {best[name]['reach']:10.5f}  "
        nulls[name] = {}
        for sh in SHUF:
            v = np.array(raw[name][sh])
            nulls[name][sh] = {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
                               "max": float(v.max()), "min": float(v.min()),
                               "p_ge_real": float(np.mean(v >= best[name]["reach"]))}
            line += f"  {v.mean():.5f}+-{v.std(ddof=1):.5f} p={nulls[name][sh]['p_ge_real']:.2f}"
        print(line)
    out["null"] = nulls

    with open(os.path.join(WORK, "run190.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nwrote {WORK}/run190.json")


if __name__ == "__main__":
    main()
