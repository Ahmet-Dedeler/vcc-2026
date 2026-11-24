"""190 Phase 2c/3c/4 — the comparator, the `fid` curve, error bars, and the direction blend.

Sections
  F  FALSIFIER 3 -- does the regulon beat `|H| x expression`, the cheap baseline 189 set?
  G  the `fid` curve, against an ASSUMED n_real distribution (rule 9)
  H  paired-bootstrap error bars on every headline reach delta
  I  PHASE 4 -- the DIRECTION blend (108's use, tydymy's use), all five h1-test members
  J  FALSIFIER 4 -- retrodicting tydymy's published board arm
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
from key190 import agreement, key_mult  # noqa: E402
from run190 import LAM, SPAN, build_all, to_overall  # noqa: E402


def reach_per_target(S, key) -> np.ndarray:
    """debench's own reach rule, but returned per perturbation so deltas can be paired."""
    out = []
    for p in range(len(S.perts)):
        pool = np.flatnonzero(S.real_sig[p])
        if pool.size == 0:
            out.append(np.nan)
            continue
        k = np.asarray(key[p], dtype=float)[pool]
        order = np.lexsort((S.genes[pool], -k))
        sel = pool[order]
        sel = sel[S.adj[p, sel]]
        if sel.size == 0:
            out.append(0.0)
            continue
        ok = np.sign(S.pred_lfc[p, sel]) == np.sign(S.real_lfc[p, sel])
        purity = np.cumsum(ok) / np.arange(1, sel.size + 1)
        deep = np.flatnonzero(purity >= 0.9)
        out.append((int(deep[-1] + 1) if deep.size else 0) / int(S.real_sig[p].sum()))
    return np.asarray(out, dtype=float)


def paired_boot(a: np.ndarray, b: np.ndarray, n=20000, seed=0):
    """Paired bootstrap over targets of mean(a) - mean(b).  172 §4.1's estimator."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, a.size, size=(n, a.size))
    d = a[idx].mean(1) - b[idx].mean(1)
    return dict(delta=float(a.mean() - b.mean()), se=float(d.std(ddof=1)),
                lo=float(np.quantile(d, 0.025)), hi=float(np.quantile(d, 0.975)),
                p_better=float(np.mean(d > 0)))


def main() -> None:
    import harness172

    T = harness172.H1Test()
    print("KNOWN-ANSWER ROWS")
    harness172.check(T)
    S = T.S
    H = T.H.astype(np.float64)
    absH = np.abs(H)
    pool = T.pool
    sT = np.sign(S.real_lfc)
    base = T.reach(absH)
    assert abs(base - 0.23533) < 5e-6
    out: dict = {}
    nets = build_all(T, H)

    # =================================================================== F
    print("\n" + "=" * 100)
    print("F. FALSIFIER 3 — the cheap comparator: control EXPRESSION.")
    print("   189 (`work/pilot189.json`): control expression predicts per-gene cross-line")
    print("   sign reliability with a held-out decile spread of +0.0505, against +0.0174")
    print("   for the residual after expression is removed.  A regulon prior has to beat")
    print("   THAT, not beat random.")
    print("=" * 100)
    cpm = np.maximum(T.cpm, 1e-6)
    lcpm = np.log10(cpm)
    z = (lcpm - lcpm.mean()) / lcpm.std()
    Z = np.broadcast_to(z, H.shape)
    rows = {"|H| (benchmark)": base,
            "expression alone (log cpm)": T.reach(Z.copy())}
    for mu in (-1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 2.0):
        rows[f"|H| * exp({mu:+.2f} * z_expr)"] = T.reach(absH * np.exp(mu * Z))
    for name, v in rows.items():
        print(f"   {name:34} reach {v:.5f}   raw {v - base:+.5f}   "
              f"overall {to_overall('reach', v - base):+.5f}")
    best_expr = max(v for k, v in rows.items() if k.startswith("|H| *"))
    print(f"\n   best expression-modulated key: {best_expr:.5f} "
          f"(raw {best_expr - base:+.5f})")
    # h1-test's own version of 189's statistic, for scale
    dec_e = T.deciles(Z.copy())
    dec_h = T.deciles(absH)
    print(f"   h1-test decile sign-accuracy spread (d1 - d10): "
          f"expression {dec_e[0] - dec_e[-1]:+.4f}   |H| {dec_h[0] - dec_h[-1]:+.4f}")
    out["falsifier3_expression"] = {"rows": rows, "best_expr_key": best_expr,
                                    "decile_spread_expr": dec_e[0] - dec_e[-1],
                                    "decile_spread_absH": dec_h[0] - dec_h[-1]}

    # network best, for the comparison
    netbest = {}
    for name, N in nets.items():
        R = agreement(N, H)
        vals = [T.reach(key_mult(H, R, l)) for l in LAM]
        j = int(np.argmax(vals))
        netbest[name] = (LAM[j], vals[j])
    print("\n   REGULON vs EXPRESSION, both selected the same way (max over the same grid):")
    print(f"   {'expression-modulated |H|':34} {best_expr:.5f}  raw {best_expr - base:+.5f}")
    for name, (l, v) in netbest.items():
        print(f"   {name + f' (lam={l:g})':34} {v:.5f}  raw {v - base:+.5f}")
    out["net_best"] = {k: {"lambda": v[0], "reach": v[1]} for k, v in netbest.items()}

    # =================================================================== G
    print("\n" + "=" * 100)
    print("G. THE `fid` CURVE.  F_p = k / max(n_pred, n_real).  `n_real` is a DISTRIBUTION")
    print("   (audit 181: board median 90-110, heavy tail), so it is swept, not assumed.")
    print("=" * 100)
    adj = S.adj
    keys = {"|H|": absH, "expression": Z.copy()}
    for name, (l, _) in netbest.items():
        keys[f"{name}(lam={l:g})"] = key_mult(H, agreement(nets[name], H), l)
    NPRED = (25, 50, 90, 150, 300)
    NREAL = (70, 90, 110, 150, 250)
    fid_out = {}
    print(f"   {'key':30} {'n_pred':>7} {'precision':>10} " +
          "".join(f"{'fid@nr=' + str(nr):>12}" for nr in NREAL))
    for kname, K in keys.items():
        fid_out[kname] = {}
        for m in NPRED:
            hits, tot = [], []
            for p in range(len(S.perts)):
                cand = np.flatnonzero(adj[p])
                if cand.size == 0:
                    continue
                kk = np.asarray(K[p], dtype=float)[cand]
                sel = cand[np.lexsort((S.genes[cand], -kk))][:m]
                hits.append(int((np.sign(S.pred_lfc[p, sel]) ==
                                 np.sign(S.real_lfc[p, sel])).sum()))
                tot.append(sel.size)
            hits = np.array(hits, float)
            tot = np.array(tot, float)
            prec = float(hits.sum() / tot.sum())
            fids = {nr: float(np.mean(hits / np.maximum(tot, nr))) for nr in NREAL}
            fid_out[kname][m] = {"precision": prec, **{str(k): v for k, v in fids.items()}}
            print(f"   {kname if m == NPRED[0] else '':30} {m:7d} {prec:10.4f} " +
                  "".join(f"{fids[nr]:12.4f}" for nr in NREAL))
    out["fid"] = fid_out
    dfid = {}
    for kname in keys:
        if kname == "|H|":
            continue
        dfid[kname] = {m: {nr: fid_out[kname][m][str(nr)] - fid_out["|H|"][m][str(nr)]
                           for nr in NREAL} for m in NPRED}
    worst = max(abs(v) for k in dfid for m in dfid[k] for v in dfid[k][m].values())
    print(f"\n   largest |Δfid| of ANY regulon/expression key vs |H|, over the whole grid: "
          f"{worst:.5f}  ->  overall {to_overall('fid', worst):+.5f}")
    out["fid_delta_max"] = worst

    # =================================================================== H
    print("\n" + "=" * 100)
    print("H. ERROR BARS — paired bootstrap over the 85 targets, 20,000 resamples (172 §4.1)")
    print("=" * 100)
    pt_base = reach_per_target(S, absH)
    assert abs(np.nanmean(pt_base) - base) < 1e-9
    print(f"   {'key':34} {'reach':>9} {'Δ vs |H|':>10} {'SE':>9} "
          f"{'95% CI':>20} {'P(better)':>10}")
    errs = {}
    cand = {"expression-modulated |H| (best)": None}
    for name, (l, _) in netbest.items():
        cand[f"{name} (lam={l:g})"] = key_mult(H, agreement(nets[name], H), l)
    for mu in (-1.0, -0.5, -0.25, 0.25, 0.5, 1.0, 2.0):
        if abs(T.reach(absH * np.exp(mu * Z)) - best_expr) < 1e-12:
            cand["expression-modulated |H| (best)"] = absH * np.exp(mu * Z)
    for name, K in cand.items():
        pt = reach_per_target(S, K)
        e = paired_boot(pt, pt_base)
        errs[name] = e
        print(f"   {name:34} {np.nanmean(pt):9.5f} {e['delta']:+10.5f} {e['se']:9.5f} "
              f"  [{e['lo']:+.5f},{e['hi']:+.5f}] {e['p_better']:10.3f}")
    out["bootstrap"] = errs

    # =================================================================== I / J
    print("\n" + "=" * 100)
    print("I/J. PHASE 4 — the DIRECTION blend (108's use, and tydymy's).")
    print("     pred_lfc' = pred_lfc + c * unit(correction), per target, c = fraction of")
    print("     the backbone's own L2.  tydymy published c = 0.25.")
    print("     All five h1-test members are exact functions of pred_lfc' given a fixed")
    print("     call set, so the ORDER channel and the DIRECTION channel are separated:")
    print("       key held at |H|            -> the direction channel alone")
    print("       key = |pred_lfc'|          -> what a magnitude correction actually")
    print("                                     delivers to the scorer's ordering")
    print("=" * 100)
    P0 = S.pred_lfc.astype(np.float64)
    n0 = np.linalg.norm(P0, axis=1, keepdims=True)
    m0 = {"pds": S.pds(P0), "nmae": S.nmae(P0), "reach|H|": T.reach(absH),
          "reach|pred|": T.reach(np.abs(P0)), "fid": S.fid(), "jac": S.jac(),
          "sign_acc": float(np.mean(np.sign(P0)[pool] == sT[pool]))}
    print(f"   {'correction':22} {'c':>6} {'pds':>9} {'nmae':>9} {'reach|H|':>10} "
          f"{'reach|pred|':>12} {'fid':>8} {'sgnacc':>8} {'||err||rel':>11}")
    print(f"   {'(shipped a15)':22} {'-':>6} {m0['pds']:9.5f} {m0['nmae']:9.5f} "
          f"{m0['reach|H|']:10.5f} {m0['reach|pred|']:12.5f} {m0['fid']:8.5f} "
          f"{m0['sign_acc']:8.4f} {1.0:11.4f}")
    blends = {}
    err0 = float(np.abs(P0 - S.real_lfc)[np.isfinite(S.real_lfc)].sum())
    for cname in ("ct-direct", "ct-prop-k25", "ct-prop-kNone", "coexpr-pc0-top500"):
        N = nets[cname]
        u = N / np.maximum(np.linalg.norm(N, axis=1, keepdims=True), 1e-30)
        blends[cname] = {}
        for c in (0.05, 0.10, 0.25, 0.50, 1.00):
            P = P0 + c * n0 * u
            Sp = S.pred_lfc
            S.pred_lfc = P                      # every member below reads the blended lfc
            try:
                r = {"pds": S.pds(), "nmae": S.nmae(), "fid": S.fid(), "jac": S.jac(),
                     "sign_acc": float(np.mean(np.sign(P)[pool] == sT[pool])),
                     "reach|H|": S.reach(key=absH), "reach|pred|": S.reach(key=np.abs(P))}
            finally:
                S.pred_lfc = Sp
            r["err_rel"] = float(np.abs(P - S.real_lfc)[np.isfinite(S.real_lfc)].sum()) / err0
            blends[cname][c] = r
            print(f"   {cname if c == 0.05 else '':22} {c:6.2f} {r['pds']:9.5f} "
                  f"{r['nmae']:9.5f} {r['reach|H|']:10.5f} {r['reach|pred|']:12.5f} "
                  f"{r['fid']:8.5f} {r['sign_acc']:8.4f} {r['err_rel']:11.4f}")
    out["direction_blend"] = {"base": m0, "blends": blends}

    print("\n   tydymy's published board move, for the retrodiction:")
    print("     reach 0.2196 -> 0.2549 (+0.0353) | fid 0.4789 -> 0.4922 (+0.0133)")
    print("     jac +0.0030 | pds -0.0013 | nmae 0.9792 -> 1.0213 (WORSE by 0.0421)")
    print("     log mse +0.108 (both ends clamped, so free) | net overall +0.0021")

    with open(os.path.join(WORK, "comp190.json"), "w") as f:
        json.dump(out, f, indent=2, default=float)
    print(f"\nwrote {WORK}/comp190.json")


def T2reach(T, S, key):
    """reach with S.pred_lfc temporarily swapped (T.S IS S, so debench sees the new lfc)."""
    return S.reach(key=key)


if __name__ == "__main__":
    main()
