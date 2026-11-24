"""521 — separate the two channels a per-target call-set policy mixes:
  (1) GLOBAL size: every target at the same K (CallSet(400, K - 400, 1.9); K < 400: the head's top K) — a constant, no
      predictor;
  (2) PER-TARGET information: the strength key's increment over the best global action, measured as
        Bayes(real key)  -  Bayes(permuted key)            (the permuted key's Bayes action = the best global action it
                                                             can learn out of sample; 40 permutations)
      and rank policies layered on the best global K (bottom q % by strength -> K_low, the rest at K*).
Every row on 3 emission seeds (435's s0/s1/s2), both emitter variants, every regime; target-bootstrap CI on seed 0.
Board-calibrated reading (docs/H1-BOARD-CALIBRATION.md; 192 §3.2 and the rc_callset board entry): fid and jac at 0.72x;
reach in the call-set class transferred at ~1x once (rc_callset 09-03) and inverted once (v14, board = -0.76 x H1): reported
as three columns (reach x 1, x 0, x -0.76); nmae x 1 (variant a only; b pins it).

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/global521.py
"""
from __future__ import annotations

import json
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import lib521 as L  # noqa: E402
from policy521 import SGN, Table, bayes_policy, boot, means, per_target_contrib, rank_policy, top3_share, value  # noqa: E402

KGLOB = (400, 800, 1200, 1600, 2000, 2200)
REG = ("1.0", "0.5", "0.35", "cmp")


def board_cols(dm):
    """board-overall predictions from raw H1 member deltas dm = (fid, jac, reach, nmae)."""
    base = 0.72 * (SGN[0] * dm[0] + SGN[1] * dm[1]) + SGN[3] * dm[3]
    return base + SGN[2] * dm[2], base, base - 0.76 * SGN[2] * dm[2]


def feats():
    P = np.load(L.WORK / "pred521_h1.npz")
    return {"S4": P["S4__S3_p_vip_thr4"], "S3": P["S3__k562phct116phek293t_nreal521_thr4"], "S4avg": P["S4_thr-avg"]}


def main():
    T = {s: Table(s) for s in (0, 1, 2)}
    F = feats()
    out = {}
    for reg in REG:
        n = T[0].nreal[reg]
        print(f"\n=== regime {reg}: n_real median {np.median(n):.0f}, mean {n.mean():.0f}, <20 {int((n < 20).sum())}, "
              f"zero {int((n == 0).sum())}, >2600 {int((n > 2600).sum())}, >1600 {int((n > 1600).sum())}")
        for var in ("b", "a"):
            print(f"  -- variant {var}: GLOBAL K (all 85 targets); value per seed s0/s1/s2, seed-0 CI, members (s0), "
                  f"board cols [reach x1 | x0 | x-0.76]")
            gv = {}
            for K in KGLOB:
                vals = []
                for s in (0, 1, 2):
                    R0 = T[s].rows[f"{reg}|base"]
                    R = T[s].assemble(reg, "headcpm", var, [K] * 85)
                    vals.append(value(R, R0))
                R0 = T[0].rows[f"{reg}|base"]
                R = T[0].assemble(reg, "headcpm", var, [K] * 85)
                lo, hi, p0 = boot(R, R0)
                dm = means(R) - means(R0)
                bc = board_cols(dm)
                gv[K] = float(np.mean(vals))
                out[f"{reg}|{var}|GLOBAL K{K}"] = dict(seeds=vals, ci=[lo, hi], p_le0=p0, d=dm.tolist(), board=list(bc))
                print(f"    K {K:5d}: {vals[0]:+.5f} {vals[1]:+.5f} {vals[2]:+.5f}  CI [{lo:+.5f}, {hi:+.5f}] P<=0 {p0:.3f}  "
                      f"fid {dm[0]:+.5f} jac {dm[1]:+.5f} reach {dm[2]:+.5f} nmae {dm[3]:+.5f}  board [{bc[0]:+.5f} | "
                      f"{bc[1]:+.5f} | {bc[2]:+.5f}]")
            Kstar = max(gv, key=gv.get)
            print(f"    best global K* = {Kstar} ({gv[Kstar]:+.5f} mean over seeds)")
            # per-target increment: Bayes(real key) - Bayes(permuted key)
            for fname, z in F.items():
                vb, vn = [], []
                for s in (0, 1, 2):
                    R0 = T[s].rows[f"{reg}|base"]
                    vb.append(value(T[s].assemble(reg, "headcpm", var, bayes_policy(T[s], reg, "headcpm", var, z)), R0))
                rng = np.random.default_rng(521)
                for d in range(40):
                    zp = z[rng.permutation(85)]
                    R0 = T[0].rows[f"{reg}|base"]
                    vn.append(value(T[0].assemble(reg, "headcpm", var, bayes_policy(T[0], reg, "headcpm", var, zp)), R0))
                vn = np.asarray(vn)
                inc = vb[0] - vn.mean()
                out[f"{reg}|{var}|{fname} BAYES vs permuted"] = dict(real=vb, perm_mean=float(vn.mean()),
                                                                    perm_q=[float(np.quantile(vn, 0.05)), float(np.quantile(vn, 0.95))],
                                                                    p_perm_ge_real=float(np.mean(vn >= vb[0])), increment=inc)
                print(f"    {fname:5s} BAYES real {vb[0]:+.5f} {vb[1]:+.5f} {vb[2]:+.5f} | permuted-key BAYES mean {vn.mean():+.5f} "
                      f"[5-95 % {np.quantile(vn, 0.05):+.5f}, {np.quantile(vn, 0.95):+.5f}]  P(perm >= real) {np.mean(vn >= vb[0]):.2f}"
                      f"  -> key increment {inc:+.5f}")
            # rank policies on top of the best global K*: bottom q % -> K_low, the rest at K*
            z = F["S4"]
            for q, Klow in ((10, 10), (10, 25), (10, 100), (20, 100), (20, 400), (30, 400)):
                if Klow >= Kstar:
                    continue
                vals = []
                for s in (0, 1, 2):
                    Kv = [Kstar] * 85
                    for i in np.argsort(z, kind="stable")[: int(round(q / 100 * 85))]:
                        Kv[i] = Klow
                    R0 = T[s].rows[f"{reg}|base"]
                    Rg = T[s].assemble(reg, "headcpm", var, [Kstar] * 85)
                    vals.append(value(T[s].assemble(reg, "headcpm", var, Kv), R0) - value(Rg, R0))
                Kv = [Kstar] * 85
                for i in np.argsort(z, kind="stable")[: int(round(q / 100 * 85))]:
                    Kv[i] = Klow
                Rg = T[0].assemble(reg, "headcpm", var, [Kstar] * 85)
                lo, hi, p0 = boot(T[0].assemble(reg, "headcpm", var, Kv), Rg)
                c = per_target_contrib(T[0].assemble(reg, "headcpm", var, Kv), Rg)
                out[f"{reg}|{var}|S4 on K*{Kstar} q{q} K{Klow}"] = dict(seeds=vals, ci=[lo, hi], p_le0=p0, top3=top3_share(c))
                print(f"    S4 bottom {q:2d}% -> {Klow:4d} on top of K*={Kstar}: {vals[0]:+.5f} {vals[1]:+.5f} {vals[2]:+.5f} "
                      f"CI [{lo:+.5f}, {hi:+.5f}] P<=0 {p0:.3f} top3 {top3_share(c):.2f}")
    (L.WORK / "global521.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
