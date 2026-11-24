"""220 — THE LOOPHOLE IN 189's CAP ARGUMENT, tested properly.

189's cap is about the MEDIAN perturbation.  `reach` averages `k*/n_real` over perturbations
and 172 measured that 94.4% of h1-test's `reach` lives in ~32 of 85 targets.  So a key that
helps only where the front is already nearly pure could pay even if the median is untouched.
Nobody had tested that.  Three questions, in order:

  1  Is the gain in fact concentrated where `reach` lives?  (per-target decomposition)
  2  What is the CEILING of any gating rule -- i.e. apply the key only on the perturbations
     where it helps?  That is an oracle over targets and it bounds the whole family.
  3  Is any such gate LEARNABLE from something we can see at prediction time?  Every
     stratifier here is computable from our own `H` and our own call set alone -- it never
     touches the reference.  Selected on half the targets, scored on the other half.

Rule: no stratifier may read `real_sig`, `real_lfc`, or anything else on the reference side.

    uv run python loop220.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

import common220 as C

sys.path.insert(0, str(C.ROOT / "experiments/190-regulon-ordering"))


def observables(T: C.Instrument) -> dict:
    """Per-target statistics computable from OUR side only (no reference is touched)."""
    sc = T.scored_mask()
    A = np.abs(T.H)
    srt = -np.sort(-A, axis=1)
    n_pred = (T.S.pred_sig & T.S.adj).sum(1).astype(float)
    o = {
        "|H| row norm": np.linalg.norm(T.H, axis=1),
        "|H| top10 mean": srt[:, :10].mean(1),
        "|H| head concentration": srt[:, :10].sum(1) / np.maximum(srt.sum(1), 1e-300),
        "|H| top10 / top100": (srt[:, :10].mean(1)
                               / np.maximum(srt[:, :100].mean(1), 1e-300)),
        "our n_pred": n_pred,
        "theta of our top10": np.zeros(len(T.perts)),
        "mean |H| gap at 10": srt[:, 9] - srt[:, 19],
    }
    order = np.argsort(-A, axis=1)
    for i in range(len(T.perts)):
        o["theta of our top10"][i] = T.theta[order[i, :10]].mean()
    return {k: v[sc] for k, v in o.items()}


def main() -> None:
    out = {}
    T = C.Instrument(C.DE_RC, "rc_callset")
    old = C.Instrument(C.DE_A15, "a15")
    C.check_old(old)
    base = T.reach_vec(T.key_absH())
    orc = T.reach_vec(T.key_oracle())
    n = base.size

    import nets190
    N, _ = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25)
    R = C.agreement_R(N, T.H)
    KEYS = {"189 lam 0.50": C.key189(T.H, T.theta, 0.5),
            "189 lam 0.25": C.key189(T.H, T.theta, 0.25),
            "190 lam -0.25": np.abs(T.H) * np.exp(np.clip(-0.25 * R, -50, 50))}

    print("=" * 100)
    print("1. IS THE GAIN WHERE `reach` LIVES?  (per-target, rc_callset counterfactual)")
    print("=" * 100)
    share = base / base.sum()
    o = np.argsort(-base)
    cum = np.cumsum(share[o])
    n944 = int(np.searchsorted(cum, 0.944) + 1)
    print(f"  {n} scored targets; 94.4% of `reach` lives in the top {n944} of them "
          f"(172 measured 32 of 85 on a15)")
    D = {}
    for nm, k in KEYS.items():
        d = T.reach_vec(k) - base
        D[nm] = d
        in_top = d[o[:n944]].sum() / n
        print(f"  {nm:<16} delta {d.mean():+.5f}   of which {in_top:+.5f} "
              f"({100*in_top/max(d.mean(),1e-12):.0f}%) is inside those {n944} targets")
    out["concentration"] = dict(n_targets_944=n944,
                                delta={k: float(v.mean()) for k, v in D.items()})

    print("\n" + "=" * 100)
    print("2. THE CEILING OF ANY GATING RULE — apply the key only where it helps (target "
          "oracle)")
    print("=" * 100)
    print(f"{'key':<16}{'plain':>10}{'gate=oracle':>14}{'up':>6}{'down':>6}"
          f"{'gate=oracle board*':>20}")
    for nm, d in D.items():
        ceil = float(np.maximum(d, 0).mean())
        board = ceil * 1.290 * 1.12927 / 6
        print(f"{nm:<16}{d.mean():>+10.5f}{ceil:>+14.5f}{int((d>1e-12).sum()):>6}"
              f"{int((d<-1e-12).sum()):>6}{board:>+20.6f}")
        out.setdefault("gate_ceiling", {})[nm] = dict(plain=float(d.mean()), ceiling=ceil,
                                                      board=float(board))
    print("  *board = raw x 1.290 (185's measured h1-test -> board factor for an "
          "order-changing arm) x 1/span 1.12927 / 6")
    print("  The rule-13 gate for a profile-PRESERVING emission arm is +0.0019 overall, "
          "i.e. +0.00783 raw h1-test `reach`.")

    print("\n" + "=" * 100)
    print("3. IS THE GATE LEARNABLE?  Stratifiers from OUR side only; split-half honest.")
    print("=" * 100)
    O = observables(T)
    print(f"{'stratifier':<26}" + "".join(f"{nm:>18}" for nm in D))
    corr = {}
    for s, v in O.items():
        cs = []
        for nm, d in D.items():
            c = float(np.corrcoef(v, d)[0, 1]) if np.std(v) > 0 else 0.0
            cs.append(c)
            corr.setdefault(s, {})[nm] = c
        print(f"{s:<26}" + "".join(f"{c:>+18.3f}" for c in cs))
    out["stratifier_corr"] = corr

    print("\n  Held-out value of the best (stratifier, direction, threshold) gate, 400 "
          "random halves:")
    rng = np.random.default_rng(20260904)
    stat = list(O.items())
    res = {}
    for nm, d in D.items():
        cands = []
        for s, v in stat:
            for sign in (+1, -1):
                for q in (0.25, 0.50, 0.75):
                    thr = np.quantile(sign * v, 1 - q)
                    cands.append((f"{s} {'hi' if sign > 0 else 'lo'} top{int(q*100)}",
                                  (sign * v) >= thr))
        cands.append(("ALL targets (no gate)", np.ones(n, bool)))
        held = []
        for _ in range(400):
            idx = rng.permutation(n)
            a, b = idx[: n // 2], idx[n // 2:]
            j = int(np.argmax([float(np.where(m, d, 0.0)[a].mean()) for _, m in cands]))
            held.append(float(np.where(cands[j][1], d, 0.0)[b].mean()))
        held = np.asarray(held)
        best_in = max((float(np.where(m, d, 0.0).mean()), lbl) for lbl, m in cands)
        print(f"  {nm:<16} in-sample best gate {best_in[0]:+.5f} ({best_in[1]})"
              f"   |   HELD-OUT {held.mean():+.5f} +- {held.std(ddof=1)/20:.5f}"
              f"   P(>0) {(held > 0).mean():.3f}")
        res[nm] = dict(in_sample=best_in[0], label=best_in[1], held=float(held.mean()),
                       held_se=float(held.std(ddof=1) / 20),
                       p=float((held > 0).mean()))
    out["gate_learnable"] = res

    print("\n" + "=" * 100)
    print("4. WHERE THE HEADROOM ACTUALLY IS — the oracle's own per-target profile")
    print("=" * 100)
    hr = orc - base
    print(f"  oracle reach {orc.mean():.5f}, incumbent {base.mean():.5f}, headroom "
          f"{hr.mean():+.5f}")
    print(f"  headroom in the top-{n944} reach targets: {hr[o[:n944]].sum()/n:+.5f} "
          f"({100*hr[o[:n944]].sum()/n/hr.mean():.0f}%)")
    q = np.argsort(-hr)
    print(f"  94.4% of the HEADROOM lives in "
          f"{int(np.searchsorted(np.cumsum(hr[q]/hr.sum()), 0.944) + 1)} targets")
    print(f"  correlation(headroom, incumbent reach) = "
          f"{np.corrcoef(hr, base)[0,1]:+.3f}   -- the targets with the most to gain are "
          f"{'the same as' if np.corrcoef(hr, base)[0,1] > 0.3 else 'NOT the same as'} "
          f"the ones already scoring")
    out["headroom"] = dict(mean=float(hr.mean()),
                           in_top=float(hr[o[:n944]].sum() / n),
                           corr=float(np.corrcoef(hr, base)[0, 1]))

    print("\n" + "=" * 100)
    print("5. WHY theta CANNOT MOVE THE FRONT — its dynamic range collapses exactly there")
    print("=" * 100)
    import sys as _sys
    _sys.path.insert(0, str(C.ROOT / "experiments/185-rank-controlled-emission"))
    import harness185 as H185
    Hd = H185.load_H()
    H0 = Hd["H"]
    g = np.asarray(Hd["genes"]).astype(str)
    th = C.theta_on(g)
    A = np.abs(H0)
    inhead = np.zeros(len(g), bool)
    for i in range(A.shape[0]):
        inhead[np.argsort(-A[i], kind="stable")[:400]] = True
    print(f"  genes that enter SOME target's top-400 |H| head: {int(inhead.sum())} / {len(g)}")
    print(f"  sd(theta):  overall {th.std():.4f}   INSIDE those genes {th[inhead].std():.4f}"
          f"   outside {th[~inhead].std():.4f}   ratio {th[~inhead].std()/th[inhead].std():.2f}x")
    print(f"  mean(theta): overall {th.mean():+.4f}   inside {th[inhead].mean():+.4f}")
    print("  A gene-level tilt whose dynamic range collapses on the front cannot reorder the")
    print("  front -- which is why the REAL theta perturbs the delivered call-set head LESS")
    print("  (Jaccard 0.804) than its own permutation does (0.728-0.754).")
    out["theta_range"] = dict(n_head_genes=int(inhead.sum()), n_genes=int(len(g)),
                              sd_all=float(th.std()), sd_in=float(th[inhead].std()),
                              sd_out=float(th[~inhead].std()),
                              mean_all=float(th.mean()), mean_in=float(th[inhead].mean()))

    json.dump(out, open(C.WORK / "loop220.json", "w"), indent=1, default=float)
    print(f"\nwrote {C.WORK/'loop220.json'}")


if __name__ == "__main__":
    main()
