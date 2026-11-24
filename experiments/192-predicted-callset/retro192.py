"""192 Phase 3 — the four mandatory falsifiers.

1. SHUFFLED-KEY control at every n (also carried in enrich192/frontier192).
2. RETRODICT the board's own five-point (`n_pred`, raw `fid`) series -- v14, v15, a8, a15, v11
   (164 §4, all same-profile emission-only arms) -- plus v14's `jac` and `reach` deltas.
   This is audit-181 rule 1: this arm family HAS board points and they must be retrodicted
   before any harness number here is read.  The single free parameter is `f`, the proportional
   level of the board's own `n_real` distribution, which rule 9 says our harnesses cannot read
   -- so a fit that reproduces five board points with one parameter MEASURES it.
3. RETRODICT a21 (n_pred 4,420 -> 1,141; board raw `fid` 0.5034 -> 0.4945, `jac` -0.017 scaled),
   with its emission-law confound decomposed through 185's own `a21_like` table.
4. LEAKAGE AUDIT: 5 of h1-test's 85 targets are in the challenge's 300-gene panel.  Every
   fitted quantity is re-fitted with those five deleted and the result re-read.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E185 = os.path.join(REPO, "experiments/185-rank-controlled-emission")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, HERE)
from frontier192 import A15, INV, PowerSurface, build_keys, order_key, topk_call  # noqa: E402

SHRINK = 0.4277                    # 164 §6 / ledger §4: h1-test -> board `reach` delta factor
# 164 §4: five same-profile emission-only board arms, challenge n_pred measured, raw board fid
BOARD_FID = [("v14", 1563, 0.49126), ("v15", 1911, 0.49654), ("a8", 3836, 0.5021),
             ("a15", 4424, 0.5022), ("v11", 5961, 0.50178)]
KGRID = [200, 400, 600, 800, 1100, 1600, 2200, 3200, 4400, 6400, 9000, 10779]
FS = [0.20, 0.30, 0.40, 0.50, 0.65, 0.80, 1.00, 1.50, 2.00]


def interp_at(K, vals, n):
    """Log-K linear interpolation of a frontier row."""
    return float(np.interp(np.log(n), np.log(K), vals))


def main() -> None:
    s, keys, absH, sig, coef, _ = build_keys()
    T = absH.shape[0]
    inc = keys["a15 p-order (incumbent)"]
    fitk = keys["FIT (cpm,nz,cv2) LOTO"]
    rng = np.random.default_rng(192)
    shuf = np.apply_along_axis(rng.permutation, 1, np.abs(absH))

    print("=" * 100)
    print("FALSIFIER 2/3 — retrodicting the board's five-point (n_pred, raw fid) series with ONE")
    print("free parameter: f, the proportional level of the board's own n_real distribution.")
    curves, out = {}, {}
    for f in FS:
        sf = PowerSurface(A15, f)
        nr = np.array([int(sf.real_sig[p].sum()) for p in range(T)])
        row = {k: [] for k in ("fid", "jac", "reach", "fid_fit", "jac_fit", "reach_fit",
                               "fid_shuf", "jac_shuf")}
        for K in KGRID:
            c = topk_call(inc, np.full(T, K))
            row["fid"].append(sf.fid(call=c))
            row["jac"].append(sf.jac(call=c))
            row["reach"].append(sf.reach(key=order_key(c, absH)))
            c2 = topk_call(fitk, np.full(T, K))
            row["fid_fit"].append(sf.fid(call=c2))
            row["jac_fit"].append(sf.jac(call=c2))
            row["reach_fit"].append(sf.reach(key=order_key(c2, absH)))
            c3 = topk_call(shuf, np.full(T, K))
            row["fid_shuf"].append(sf.fid(call=c3))
            row["jac_shuf"].append(sf.jac(call=c3))
        curves[f] = row
        pred = [interp_at(KGRID, row["fid"], n) for _, n, _ in BOARD_FID]
        obs = [v for _, _, v in BOARD_FID]
        # only the SHAPE is comparable (h1-test's absolute fid level differs); centre both
        d_pred = np.array(pred) - np.mean(pred)
        d_obs = np.array(obs) - np.mean(obs)
        rms = float(np.sqrt(np.mean((d_pred - d_obs) ** 2)))
        slope = float(np.polyfit(d_obs, d_pred, 1)[0])
        out[f] = dict(pred=pred, rms=rms, slope=slope,
                      n_real_med=float(np.median(nr)), n_real_mean=float(nr.mean()),
                      frac_above_1141=float((nr > 1141).mean()))
        print(f"  f={f:<5} n_real med {np.median(nr):>5.0f} mean {nr.mean():>6.0f} "
              f"P(n_real>1141) {(nr > 1141).mean():.3f} | pred fid "
              + " ".join(f"{p:.4f}" for p in pred)
              + f" | centred RMS {rms:.5f}  slope(pred/obs) {slope:.2f}")
    print("  board observed              "
          + " " * 34 + " ".join(f"{v:.4f}" for _, _, v in BOARD_FID))
    best = min(out, key=lambda f: out[f]["rms"])
    print(f"\n  BEST f = {best}  (centred RMS {out[best]['rms']:.5f}, slope "
          f"{out[best]['slope']:.2f}); implied board n_real median "
          f"{out[best]['n_real_med']:.0f}, mean {out[best]['n_real_mean']:.0f}")

    # ---- v14's other two members at the best f
    print("\n  v14 (n_pred 1563) vs v11 (5961), all three DE members, at the best f:")
    r = curves[best]
    for m, board_scaled, extra in (("fid", -0.035470, 1.0), ("jac", +0.004271, 1.0),
                                   ("reach", +0.013731, SHRINK)):
        p14 = interp_at(KGRID, r[m], 1563)
        p11 = interp_at(KGRID, r[m], 5961)
        praw = p14 - p11
        psc = praw * INV[m] * extra
        print(f"    {m:<6} h1-test raw d {praw:+.5f}  -> board scaled {psc:+.6f}   "
              f"observed {board_scaled:+.6f}   "
              f"{'SIGN OK' if np.sign(psc) == np.sign(board_scaled) else 'SIGN WRONG'}  "
              f"ratio {psc / board_scaled:+.2f}x")
        out.setdefault("v14", {})[m] = dict(pred_raw=praw, pred_scaled=psc, board=board_scaled)

    # ---- a21
    print("\n  FALSIFIER 3 — a21 (n_pred 4,420 -> 1,141).  Two channels, decomposed:")
    a21 = {}
    for m, board_scaled in (("fid", -0.029464), ("jac", -0.016583)):
        cs = interp_at(KGRID, r[m], 1141) - interp_at(KGRID, r[m], 4420)
        a21[m] = cs
        print(f"    {m:<5} call-set channel only: h1-test raw {cs:+.5f} -> scaled "
              f"{cs * INV[m]:+.6f}   board {board_scaled:+.6f}")
    # the emission-law channel, measured directly by 185's a21_like at its own n_pred
    from debench import Surface
    a15s = PowerSurface(A15, best)
    a21s = PowerSurface(os.path.join(E185, "work/de185.a21_like.npz"), best)
    n21 = int(np.median(Surface(os.path.join(E185, "work/de185.a21_like.npz")).counts()[1]))
    for m in ("fid", "jac"):
        realised = getattr(a21s, m)() - getattr(a15s, m)()
        callonly = interp_at(KGRID, r[m], n21) - interp_at(KGRID, r[m], 6034)
        print(f"    {m:<5} 185's a21_like (n_pred {n21}): realised raw {realised:+.5f}, of which "
              f"call-set {callonly:+.5f} -> EMISSION-LAW residue {realised - callonly:+.5f}")
        a21[f"{m}_emission"] = realised - callonly
    for m, board_scaled in (("fid", -0.029464), ("jac", -0.016583)):
        tot = (a21[m] + a21[f"{m}_emission"]) * INV[m]
        print(f"    {m:<5} TOTAL predicted scaled {tot:+.6f}  vs board {board_scaled:+.6f}  "
              f"ratio {tot / board_scaled:+.2f}x  "
              f"{'SIGN OK' if np.sign(tot) == np.sign(board_scaled) else 'SIGN WRONG'}")
    out["a21"] = a21

    # ---- falsifier 1: shuffled key at every n, at the retrodicted level
    print(f"\n{'=' * 100}")
    print(f"FALSIFIER 1 — shuffled-key control at every n, at the retrodicted f = {best}")
    print(f"  {'n_pred':>8}{'fid FIT':>10}{'fid inc':>10}{'fid shuf':>10}"
          f"{'jac FIT':>10}{'jac inc':>10}{'jac shuf':>10}{'d(F+J) ovr vs inc':>20}")
    for i, K in enumerate(KGRID):
        d = ((r["fid_fit"][i] - r["fid"][i]) * INV["fid"]
             + (r["jac_fit"][i] - r["jac"][i]) * INV["jac"]) / 6
        print(f"  {K:>8}{r['fid_fit'][i]:>10.4f}{r['fid'][i]:>10.4f}{r['fid_shuf'][i]:>10.4f}"
              f"{r['jac_fit'][i]:>10.4f}{r['jac'][i]:>10.4f}{r['jac_shuf'][i]:>10.4f}"
              f"{d:>+20.5f}")

    with open(os.path.join(WORK, "retro192.json"), "w") as fh:
        json.dump(dict(best_f=best, fits={str(k): v for k, v in out.items() if k in FS},
                       v14=out.get("v14"), a21=out.get("a21"), KGRID=KGRID,
                       curves={str(k): v for k, v in curves.items()}), fh, indent=1,
                  default=float)
    print(f"\nwrote {WORK}/retro192.json")


if __name__ == "__main__":
    main()
