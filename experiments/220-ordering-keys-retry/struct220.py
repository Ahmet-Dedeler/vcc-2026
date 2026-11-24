"""220 Phase 1b — is the ordering-key family NOISE-LIMITED or STRUCTURALLY CAPPED?

Everything here is deterministic and runs on both instruments, so it settles the structural
question independently of the emission noise the mandate is about.

  1  the pure-prefix anatomy under `rc_callset` -- 189 quoted a MEDIAN k* of 3 from 172
     section 7.2, measured on a15's `quasi` emitter.  `rc_callset` changed both the call set
     and the delivered order, so the number is re-measured rather than inherited.
  2  the 0.90 bar in genes: how many of the first d must be right, and how many actually are.
  3  the concentration of `reach` (172: 94.4% in 32 of 85 targets), re-measured.
  4  the LOOPHOLE the mandate names: `reach` averages `k*/n_real` over perturbations, so a key
     that pays only where the front is already nearly pure could pay even if the median is
     untouched.  Per-target stratification of both keys' gains.
  5  HONEST SELECTION.  189 chose lambda on its evaluation surface and said so; 190 took a max
     over 13 lambdas x 8 networks.  Here lambda is chosen on half the targets and scored on
     the other half, over 400 random splits.
  6  the CEILING of the whole gene-level-key family, by cross-fold gene oracle -- the `reach`
     analogue of 189 section 6.3's 22.3x trap.
  7  permutation nulls at 400 seeds with the selection INSIDE the null.

    uv run python struct220.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

import common220 as C

sys.path.insert(0, str(C.ROOT / "experiments/190-regulon-ordering"))

LAM189 = (-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0)
LAM190 = (-4.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0)


def kstar_table(T: C.Instrument, label: str, key=None) -> dict:
    ks = T.kstar_vec(key)
    rv = T.reach_vec(key)
    nr = T.n_real[T.scored_mask()]
    share = rv / max(rv.sum(), 1e-300)
    order = np.argsort(-rv)
    cum = np.cumsum(share[order])
    n32 = int(np.searchsorted(cum, 0.944) + 1)
    return dict(label=label, reach=float(rv.mean()),
                kstar_med=float(np.median(ks)), kstar_mean=float(ks.mean()),
                kstar_p25=float(np.percentile(ks, 25)),
                kstar_p75=float(np.percentile(ks, 75)),
                kstar_p90=float(np.percentile(ks, 90)),
                frac_zero=float((ks == 0).mean()), n_scored=int(ks.size),
                n_real_med=float(np.median(nr)),
                n_targets_for_944=n32,
                top32_share=float(cum[min(31, cum.size - 1)]))


def main() -> None:
    out = {}
    old = C.Instrument(C.DE_A15, "a15")
    new = C.Instrument(C.DE_RC, "rc_callset")
    C.check_old(old)
    assert abs(new.reach() - 0.223341) < 5e-6, "the rc_callset table moved"
    print()

    # ================================================================= 1. prefix anatomy
    print("=" * 100)
    print("1. THE PURE-PREFIX ANATOMY.  189's cap argument rests on `median k* = 3` "
          "(172 section 7.2, a15).")
    print("=" * 100)
    rows = [
        kstar_table(old, "a15  DELIVERED (production order)"),
        kstar_table(old, "a15  counterfactual |H|", old.key_absH()),
        kstar_table(new, "rc_callset  DELIVERED (production order)"),
        kstar_table(new, "rc_callset  counterfactual |H|", new.key_absH()),
        kstar_table(new, "rc_callset  ORACLE", new.key_oracle()),
    ]
    print(f"{'surface':<40}{'reach':>9}{'med k*':>8}{'mean':>8}{'p75':>7}{'p90':>7}"
          f"{'k*=0':>8}{'n_real':>8}{'#tgts for 94.4%':>17}")
    for r in rows:
        print(f"{r['label']:<40}{r['reach']:>9.5f}{r['kstar_med']:>8.1f}"
              f"{r['kstar_mean']:>8.2f}{r['kstar_p75']:>7.0f}{r['kstar_p90']:>7.0f}"
              f"{r['frac_zero']:>8.3f}{r['n_real_med']:>8.0f}{r['n_targets_for_944']:>17d}")
    out["prefix"] = rows

    # ================================================================= 2. the 0.90 bar
    print("\n" + "=" * 100)
    print("2. THE 0.90 BAR IN GENES.  To have k* >= d you need ceil(0.9 d) of the first d "
          "right.")
    print("=" * 100)
    depths = (5, 10, 20, 50)
    fa_d = new.front_acc(None, depths=depths)          # delivered order
    fa_h = new.front_acc(new.key_absH(), depths=depths)
    print(f"{'d':>4}{'need':>6}{'mean acc |H|':>14}{'mean acc delivered':>20}"
          f"{'P(acc >= 0.9) |H|':>19}{'median shortfall (genes)':>26}")
    bar = {}
    for d in depths:
        need = int(np.ceil(0.9 * d))
        a = fa_h[d][np.isfinite(fa_h[d])]
        b = fa_d[d][np.isfinite(fa_d[d])]
        short = need - np.round(a * d)
        bar[d] = dict(need=need, acc_absH=float(a.mean()), acc_delivered=float(b.mean()),
                      p_clear=float((a >= 0.9 - 1e-9).mean()),
                      median_shortfall=float(np.median(short)),
                      mean_shortfall=float(short.mean()))
        print(f"{d:>4}{need:>6}{a.mean():>14.4f}{b.mean():>20.4f}"
              f"{(a >= 0.9 - 1e-9).mean():>19.4f}{np.median(short):>26.1f}")
    print("  ('shortfall' = extra genes that would have to be flipped from wrong to right "
          "in the first d)")
    out["bar"] = bar

    # required accuracy for the family: what front accuracy crosses 0.90 at each depth
    print("\n  189's own framing: front accuracy 0.587 -> 0.669 inside |H| decile 1 does not "
          "cross 0.90.")
    print("  The bar is not an average -- it is PER PERTURBATION and it is a ceiling "
          "function:")
    for d in depths:
        print(f"    d = {d:>2}: need {int(np.ceil(0.9*d))}/{d} = "
              f"{np.ceil(0.9*d)/d:.3f} exactly, i.e. at most "
              f"{d - int(np.ceil(0.9*d))} error(s) in the first {d}")

    # ================================================================= 3/4. keys, per-target
    print("\n" + "=" * 100)
    print("3. THE TWO KEYS ON BOTH INSTRUMENTS, COUNTERFACTUALLY, WITH PER-TARGET "
          "STRATIFICATION")
    print("=" * 100)
    import nets190
    keys = {}
    for T, tag in ((old, "a15"), (new, "rc")):
        th = T.theta
        keys[(tag, "189")] = C.key189(T.H, th, 0.5)
        N, _ = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25)
        R = C.agreement_R(N, T.H)
        keys[(tag, "190")] = np.abs(T.H) * np.exp(np.clip(-0.25 * R, -50, 50))
        keys[(tag, "R")] = R

    strat = {}
    for T, tag in ((old, "a15"), (new, "rc")):
        base = T.reach_vec(T.key_absH())
        orc = T.reach_vec(T.key_oracle())
        nr = T.n_real[T.scored_mask()]
        print(f"\n  --- {tag} ---   base reach {base.mean():.5f}   oracle {orc.mean():.5f}")
        for nm in ("189", "190"):
            v = T.reach_vec(keys[(tag, nm)])
            d = v - base
            bb = C.paired_bootstrap(v, base)
            nz = np.flatnonzero(np.abs(d) > 1e-12)
            # concentration: how much of the total gain comes from targets already in the
            # top tercile of the incumbent's own reach
            hi = base >= np.percentile(base, 66.7)
            lo = base <= np.percentile(base, 33.3)
            print(f"  {nm}: delta {bb['mean']:+.5f}  target-se {bb['se']:.5f}  "
                  f"P(better) {bb['p_better']:.3f}  |  targets moved {nz.size}/{d.size}  "
                  f"up {int((d > 1e-12).sum())} down {int((d < -1e-12).sum())}")
            print(f"       gain in the HIGH-reach tercile {d[hi].sum()/d.size:+.5f}  "
                  f"MID {d[~hi & ~lo].sum()/d.size:+.5f}  LOW {d[lo].sum()/d.size:+.5f}"
                  f"   (they sum to the delta)")
            print(f"       largest single-target contributions: " +
                  " ".join(f"{x:+.4f}" for x in np.sort(d)[::-1][:5]) + "  ... " +
                  " ".join(f"{x:+.4f}" for x in np.sort(d)[:3]))
            strat[f"{tag}|{nm}"] = dict(delta=bb["mean"], se=bb["se"],
                                        p=bb["p_better"], moved=int(nz.size),
                                        hi=float(d[hi].sum() / d.size),
                                        mid=float(d[~hi & ~lo].sum() / d.size),
                                        lo=float(d[lo].sum() / d.size),
                                        top5=[float(x) for x in np.sort(d)[::-1][:5]],
                                        bot3=[float(x) for x in np.sort(d)[:3]],
                                        corr_with_base=float(
                                            np.corrcoef(d, base)[0, 1]),
                                        corr_with_nreal=float(np.corrcoef(d, nr)[0, 1]))
    out["strat"] = strat

    # ============================================== 5. honest held-out lambda selection
    print("\n" + "=" * 100)
    print("5. HONEST SELECTION.  lambda picked on half the targets, scored on the other "
          "half, 400 splits.")
    print("=" * 100)
    rng = np.random.default_rng(20260904)
    sel = {}
    for T, tag in ((old, "a15"), (new, "rc")):
        base = T.reach_vec(T.key_absH())
        # 189's family
        V189 = {lam: T.reach_vec(C.key189(T.H, T.theta, lam)) - base for lam in LAM189}
        # 190's family: 13 lambdas on the one network 190 selected
        R = keys[(tag, "R")]
        V190 = {lam: T.reach_vec(np.abs(T.H) * np.exp(np.clip(lam * R, -50, 50))) - base
                for lam in LAM190}
        for nm, V, grid in (("189", V189, LAM189), ("190", V190, LAM190)):
            M = np.vstack([V[l] for l in grid])
            insample = M.mean(1)
            best_in = float(insample.max())
            lam_in = grid[int(np.argmax(insample))]
            held = []
            for _ in range(400):
                idx = rng.permutation(M.shape[1])
                a, b = idx[: M.shape[1] // 2], idx[M.shape[1] // 2:]
                j = int(np.argmax(M[:, a].mean(1)))
                held.append(float(M[j, b].mean()))
            held = np.asarray(held)
            print(f"  {tag:<4} {nm}: in-sample max {best_in:+.5f} at lambda {lam_in:+.2f}"
                  f"   |   HELD-OUT {held.mean():+.5f} +- {held.std(ddof=1)/20:.5f}"
                  f"   (P(held-out > 0) = {(held > 0).mean():.3f})")
            print(f"        the whole grid, in-sample: " +
                  " ".join(f"{l:+.2f}:{v:+.5f}" for l, v in zip(grid, insample)))
            sel[f"{tag}|{nm}"] = dict(in_sample=best_in, lam=lam_in,
                                      held=float(held.mean()),
                                      held_se=float(held.std(ddof=1) / 20),
                                      p_held=float((held > 0).mean()),
                                      grid={str(l): float(v)
                                            for l, v in zip(grid, insample)})
    out["selection"] = sel

    # ============================================== 6. the gene-level-key family ceiling
    print("\n" + "=" * 100)
    print("6. THE CEILING OF THE WHOLE GENE-LEVEL-KEY FAMILY (the `reach` analogue of "
          "189 section 6.3)")
    print("=" * 100)
    print("  A gene-level key is `|H| x u_g`.  Its BEST possible u is the gene's own sign")
    print("  accuracy.  Estimated IN FOLD it is an oracle; estimated CROSS FOLD it is the")
    print("  honest ceiling.  12 random halves.")
    ceil = {}
    for T, tag in ((old, "a15"), (new, "rc")):
        base = T.reach_vec(T.key_absH())
        ok = T.correct.astype(np.float64)
        inpool = T.pool.astype(np.float64)
        rng2 = np.random.default_rng(7)
        TILTS = (1.0, 2.0, 4.0, 8.0)
        infold = {t: [] for t in TILTS}
        cross = {t: [] for t in TILTS}
        rhos = []
        sc = T.scored_mask()
        for _ in range(12):
            idx = rng2.permutation(len(T.perts))
            A, B = idx[: len(idx) // 2], idx[len(idx) // 2:]
            nA, nB = inpool[A].sum(0), inpool[B].sum(0)
            accA = np.where(nA > 0, (ok[A] * inpool[A]).sum(0) / np.maximum(nA, 1), 0.5)
            accB = np.where(nB > 0, (ok[B] * inpool[B]).sum(0) / np.maximum(nB, 1), 0.5)
            m = (nA > 0) & (nB > 0)
            rhos.append(float(np.corrcoef(accA[m], accB[m])[0, 1]))
            for tilt in TILTS:
                for u_fit, ev, bucket in ((accA, B, cross), (accB, A, cross),
                                          (accA, A, infold), (accB, B, infold)):
                    k = np.abs(T.H) * np.exp(tilt * (u_fit - 0.5))[None, :]
                    mask = np.zeros(len(T.perts), bool)
                    mask[ev] = True
                    bucket[tilt].append(
                        float((T.reach_vec(k) - base)[mask[sc]].mean()))
        bi = max(TILTS, key=lambda t: np.mean(infold[t]))
        bc = max(TILTS, key=lambda t: np.mean(cross[t]))
        print(f"  {tag:<4} split-half rho(per-gene accuracy) {np.mean(rhos):+.4f}   "
              f"IN-FOLD gene oracle {np.mean(infold[bi]):+.5f} (tilt {bi})   "
              f"CROSS-FOLD (honest) {np.mean(cross[bc]):+.5f} +- "
              f"{np.std(cross[bc], ddof=1)/np.sqrt(len(cross[bc])):.5f} (tilt {bc})")
        print("        per tilt, cross-fold: " +
              "  ".join(f"{t:g}:{np.mean(cross[t]):+.5f}" for t in TILTS))
        ceil[tag] = dict(rho=float(np.mean(rhos)), infold=float(np.mean(infold[bi])),
                         infold_tilt=bi, cross=float(np.mean(cross[bc])), cross_tilt=bc,
                         cross_se=float(np.std(cross[bc], ddof=1)
                                        / np.sqrt(len(cross[bc]))),
                         cross_all={str(t): float(np.mean(cross[t])) for t in TILTS})
    out["family_ceiling"] = ceil

    # ============================================== 7. permutation nulls, 400 seeds
    print("\n" + "=" * 100)
    print("7. PERMUTATION NULLS AT 400 SEEDS, WITH THE SELECTION INSIDE THE NULL "
          "(190's discipline)")
    print("=" * 100)
    nulls = {}
    for T, tag in ((old, "a15"), (new, "rc")):
        base = T.reach_vec(T.key_absH())
        real189 = float((T.reach_vec(C.key189(T.H, T.theta, 0.5)) - base).mean())
        rng3 = np.random.default_rng(4242)
        vals = []
        for _ in range(400):
            th = rng3.permutation(T.theta)
            vals.append(float((T.reach_vec(C.key189(T.H, th, 0.5)) - base).mean()))
        v = np.asarray(vals)
        z = (real189 - v.mean()) / v.std(ddof=1)
        print(f"  {tag:<4} 189 (theta permuted across the gene axis, lambda 0.5 fixed): "
              f"real {real189:+.5f}   null {v.mean():+.5f} +- {v.std(ddof=1):.5f}   "
              f"z {z:+.2f}   P(null >= real) {(v >= real189).mean():.4f}")
        # with lambda re-selected per null draw, exactly as for the real key
        rng4 = np.random.default_rng(555)
        sel_vals = []
        for _ in range(120):
            th = rng4.permutation(T.theta)
            sel_vals.append(max(float((T.reach_vec(C.key189(T.H, th, l)) - base).mean())
                                for l in LAM189))
        sv = np.asarray(sel_vals)
        real_sel = max(float((T.reach_vec(C.key189(T.H, T.theta, l)) - base).mean())
                       for l in LAM189)
        print(f"       with lambda RE-SELECTED inside the null (120 draws): real "
              f"{real_sel:+.5f}   null {sv.mean():+.5f} +- {sv.std(ddof=1):.5f}   "
              f"P(null >= real) {(sv >= real_sel).mean():.4f}")
        nulls[tag] = dict(real=real189, null_mean=float(v.mean()),
                          null_sd=float(v.std(ddof=1)), z=float(z),
                          p=float((v >= real189).mean()),
                          real_sel=real_sel, null_sel_mean=float(sv.mean()),
                          null_sel_sd=float(sv.std(ddof=1)),
                          p_sel=float((sv >= real_sel).mean()))
    out["nulls"] = nulls

    json.dump(out, open(C.WORK / "struct220.json", "w"), indent=1, default=float)
    print(f"\nwrote {C.WORK/'struct220.json'}")


if __name__ == "__main__":
    main()
