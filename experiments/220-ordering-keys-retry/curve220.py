"""220 — the FRONT-ACCURACY -> `reach` TRANSFER CURVE, and what it costs to cross 0.90.

189's cap argument is verbal: *"lifting front accuracy from 0.587 to 0.669 does not cross
0.90, so `k*` does not move."*  This measures the function it is asserting.  An oracle-tilted
key `|H| x exp(c . 1[correct])` interpolates continuously from `|H|` (c = 0) to the oracle
(c -> inf).  Sweeping `c` traces (front accuracy, reach) pairs on the real surface, so the
question "what front accuracy would be needed" has a measured answer instead of an argument.

Also here:
  * the sign-identity check between a15's and `rc_callset`'s DE tables -- the reason every
    counterfactual key row in `struct220.py` came out numerically identical on both;
  * what 189's and 190's keys actually do to front accuracy, at the depths that matter.

    uv run python curve220.py
"""
from __future__ import annotations

import json
import sys

import numpy as np

import common220 as C

sys.path.insert(0, str(C.ROOT / "experiments/190-regulon-ordering"))


def main() -> None:
    out = {}
    old = C.Instrument(C.DE_A15, "a15")
    new = C.Instrument(C.DE_RC, "rc_callset")
    C.check_old(old)

    # ---------------------------------------------------------------- sign identity
    print("\n" + "=" * 100)
    print("0. WHY EVERY COUNTERFACTUAL ROW IS IDENTICAL ON BOTH TABLES")
    print("=" * 100)
    sa = np.sign(old.S.pred_lfc)
    sb = np.sign(new.S.pred_lfc)
    pool = old.pool
    agree_all = float((sa == sb).mean())
    agree_pool = float((sa == sb)[pool].mean())
    same_real = bool(np.array_equal(np.nan_to_num(old.S.real_lfc),
                                    np.nan_to_num(new.S.real_lfc)))
    same_sig = bool(np.array_equal(old.S.real_sig, new.S.real_sig))
    print(f"  reference side identical (real_lfc / real_sig): {same_real} / {same_sig}")
    print(f"  sign(pred_lfc) agreement:  whole table {agree_all:.6f}   "
          f"inside the reference DE pool {agree_pool:.6f}")
    print(f"  pool pairs whose sign DIFFERS: {int((sa != sb)[pool].sum()):,} of "
          f"{int(pool.sum()):,}")
    print("  `reach(key)` reads ONLY sign(pred_lfc) and the reference -- it never reads "
          "pred_p.")
    print("  So a counterfactual key is measured on the SAME instrument under both emitters,")
    print("  and 214's 77x seed tightening cannot touch it.  What the new emitter buys is")
    print("  DELIVERY (reach with key=None), which is measured in `six220.py`.")
    out["sign_identity"] = dict(all=agree_all, pool=agree_pool,
                                differing=int((sa != sb)[pool].sum()),
                                pool_size=int(pool.sum()),
                                real_identical=same_real and same_sig)

    # ---------------------------------------------------------------- the curve
    print("\n" + "=" * 100)
    print("1. THE FRONT-ACCURACY -> reach TRANSFER CURVE (oracle-tilted key on the real "
          "surface)")
    print("=" * 100)
    depths = (5, 10, 20, 50)
    rows = []
    for T, tag in ((old, "a15"), (new, "rc_callset")):
        base = T.reach_vec(T.key_absH())
        print(f"\n  --- {tag} ---")
        print(f"{'c':>7}{'reach':>10}{'delta':>10}{'acc@5':>9}{'acc@10':>9}{'acc@20':>9}"
              f"{'acc@50':>9}{'med k*':>9}{'P(k*>=10)':>11}")
        for c in (0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 6.0):
            k = np.abs(T.H) * np.exp(c * T.correct.astype(np.float64))
            rv = T.reach_vec(k)
            ks = T.kstar_vec(k)
            fa = T.front_acc(k, depths=depths)
            r = dict(tag=tag, c=c, reach=float(rv.mean()),
                     delta=float(rv.mean() - base.mean()),
                     kstar_med=float(np.median(ks)),
                     p_k10=float((ks >= 10).mean()),
                     **{f"acc{d}": float(np.nanmean(fa[d])) for d in depths})
            rows.append(r)
            print(f"{c:>7.2f}{r['reach']:>10.5f}{r['delta']:>+10.5f}"
                  f"{r['acc5']:>9.4f}{r['acc10']:>9.4f}{r['acc20']:>9.4f}"
                  f"{r['acc50']:>9.4f}{r['kstar_med']:>9.1f}{r['p_k10']:>11.4f}")
    out["curve"] = rows

    # local slope of reach in front accuracy, at the incumbent
    r0 = [r for r in rows if r["tag"] == "rc_callset"]
    a = np.asarray([r["acc10"] for r in r0])
    y = np.asarray([r["reach"] for r in r0])
    j = 3
    slope = (y[j] - y[0]) / (a[j] - a[0])
    print(f"\n  LOCAL SLOPE at the incumbent (c = 0 -> {r0[j]['c']}): "
          f"d reach / d acc@10 = {slope:.4f}")
    print(f"  189 measured a front-accuracy gain of +0.0820 inside |H| decile 1 "
          f"(0.5873 -> 0.6694).")
    print(f"  Decile 1 is 10% of a median-89-gene pool, i.e. ~9 genes, so the whole-front "
          f"gain it could deliver is at most that.")
    print(f"  On this measured curve +0.0820 of acc@10 is worth "
          f"{slope * 0.0820:+.5f} raw `reach` -- against 189's structural estimate of "
          f"<= +0.0005.")
    out["slope"] = dict(slope=float(slope), implied=float(slope * 0.0820))

    # ---------------------------------------------------------------- what the keys do
    print("\n" + "=" * 100)
    print("2. WHAT 189's AND 190's KEYS ACTUALLY DO TO THE FRONT")
    print("=" * 100)
    import nets190
    T = new
    N, _ = nets190.ct_prop(list(T.perts), list(T.genes), T.H, topk=25)
    R = C.agreement_R(N, T.H)
    ks = {"|H| (incumbent)": T.key_absH(),
          "189 lam 0.50": C.key189(T.H, T.theta, 0.5),
          "189 lam 0.25": C.key189(T.H, T.theta, 0.25),
          "190 ct-prop-k25 lam -0.25": np.abs(T.H) * np.exp(np.clip(-0.25 * R, -50, 50))}
    print(f"{'key':<28}{'reach':>10}{'acc@5':>9}{'acc@10':>9}{'acc@20':>9}{'acc@50':>9}"
          f"{'med k*':>9}{'P(k*>=10)':>11}")
    keyrows = {}
    for nm, k in ks.items():
        rv = T.reach_vec(k)
        kk = T.kstar_vec(k)
        fa = T.front_acc(k, depths=depths)
        keyrows[nm] = dict(reach=float(rv.mean()), kstar_med=float(np.median(kk)),
                           p_k10=float((kk >= 10).mean()),
                           **{f"acc{d}": float(np.nanmean(fa[d])) for d in depths})
        r = keyrows[nm]
        print(f"{nm:<28}{r['reach']:>10.5f}{r['acc5']:>9.4f}{r['acc10']:>9.4f}"
              f"{r['acc20']:>9.4f}{r['acc50']:>9.4f}{r['kstar_med']:>9.1f}"
              f"{r['p_k10']:>11.4f}")
    out["keys"] = keyrows

    json.dump(out, open(C.WORK / "curve220.json", "w"), indent=1, default=float)
    print(f"\nwrote {C.WORK/'curve220.json'}")


if __name__ == "__main__":
    main()
