"""Calibrate CLAIM-R7-DRAFT's board-like rule and HARM veto against 553's rehearsal.

Reads only tracked json (draws553, band553_h1r7, the five h1r7 accesses, the five-draw mean).
Does not open a truth h5ad. Prints the numbers RESULTS.md quotes.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
W553 = ROOT / "experiments/553-r7b-harness/work"
DRAWS = ("h1r7", "h1r7s2", "h1r7s3", "h1r7s4", "h1r7s5")


def Phi(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def p_le(t, mu, sd):
    if sd <= 0:
        return 1.0 if mu <= t else 0.0
    return Phi((t - mu) / sd)


def main() -> None:
    band = json.loads((W553 / "band553_h1r7.json").read_text())
    draws = json.loads((W553 / "draws553.json").read_text())
    reh = json.loads((W553 / "access553_rehearsal5/access553.json").read_text())
    print("band553_h1r7 between_draw_sd is null (run without --draws):",
          all(v.get("between_draw_sd") is None for v in band.values()))
    print("rehearsal5 CP raw jac", reh["arms"]["lcap-kept"]["effect"]["jac"]["mean"],
          "raw fid", reh["arms"]["lcap-kept"]["effect"]["fid"]["mean"])
    print("rehearsal5 has no retrodiction key:", "retrodiction" not in reh)
    print()
    print("CP raw jac / fid on each draw (board-like := jac>0 and fid<0):")
    n_like = 0
    for name in DRAWS:
        arm = json.loads((W553 / f"access553_{name}/access553.json").read_text())["arms"]["lcap-kept"]
        jac, fid = arm["effect"]["jac"]["mean"], arm["effect"]["fid"]["mean"]
        ok = jac > 0 and fid < 0
        n_like += int(ok)
        print(f"  {name}: jac {jac:+.6f} fid {fid:+.6f} board-like {ok}")
    print(f"  {n_like}/5 draws board-like")
    print()
    print(f"{'effect':<42} {'band_tot':>8} {'T=-2.5b':>9} {'mean':>8} {'sd':>7} "
          f"{'emp':>5} {'P0_b':>7} {'P0_sd':>7} {'Pmean':>7} {'ratio':>6} {'r_ctrl':>7}")
    family = []
    for k, row in band.items():
        d = draws[k]
        tot = row["band_total"]
        T = -2.5 * tot
        vals = d["values"]
        emp = sum(v <= T for v in vals)
        mu, sd = d["mean"], d["between_sd"]
        # claim's "1.3x (0.83-3.0)" is ratio_ctrl of the cells+controls band, not between/band_total
        print(f"{k:<42} {tot:8.5f} {T:9.5f} {mu:8.5f} {sd:7.5f} {emp:2d}/5 "
              f"{p_le(T, 0, tot):7.3f} {p_le(T, 0, sd):7.3f} {p_le(T, mu, sd):7.3f} "
              f"{sd / tot:6.2f} {d['ratio_ctrl']:7.2f}")
        family.append((k, vals, p_le(T, 0, sd)))
    print()
    print("Phi(-2.5) =", round(Phi(-2.5), 4), "(null tail if noise == band_total and the true effect is 0)")
    # K-on-CPH power
    k = "lcap-kept-n1600+c4h vs lcap-kept+c4h"
    tot = band[k]["band_total"]
    T = -2.5 * tot
    sd_user = 0.0023
    mu_harm = -0.0017
    prior = draws[k]["mean"]
    print()
    print("K-on-CPH power, noise = between-draw sd 0.0023 (the user's K-size harm):")
    print(f"  T = -2.5 * h1r7 band_total = {T:.5f}")
    print(f"  P(X<=T | mu=-0.0017, sd=0.0023) = {p_le(T, mu_harm, sd_user):.3f}")
    print(f"  P(X<=T | mu=rehearsal mean {prior:.5f}, sd=0.0023) = {p_le(T, prior, sd_user):.3f}")
    print(f"  P(X<=T | mu=0, sd=0.0023) = {p_le(T, 0, sd_user):.3f}")
    T_b = prior - 2.5 * tot
    T_s = prior - 2.5 * sd_user
    print(f"  proposed T = prior - 2.5*band_total = {T_b:.5f}; "
          f"P(veto | mu=-0.0017) = {p_le(T_b, mu_harm, sd_user):.3f}; "
          f"P(veto | mu=prior) = {p_le(T_b, prior, sd_user):.3f}")
    print(f"  proposed T = prior - 2.5*between_sd = {T_s:.5f}; "
          f"P(veto | mu=-0.0017) = {p_le(T_s, mu_harm, sd_user):.3f}")
    # family, independent, noise = between-draw
    groups = {
        "CP+hook+CPH": ["lcap-kept vs cpm", "cpm+c4h vs cpm", "lcap-kept+c4h vs cpm"],
        "hook + hook-on-CP": ["cpm+c4h vs cpm", "lcap-kept+c4h vs lcap-kept"],
        "K2000 + K-on-CPH": ["cpm-n1600 vs cpm", "lcap-kept-n1600+c4h vs lcap-kept+c4h"],
    }
    print()
    print("family-wise P(any veto | true 0, independent, noise=between-draw sd):")
    by = {k: p for k, _, p in family}
    vals_by = {k: v for k, v, _ in family}
    for name, keys in groups.items():
        p_any = 1.0 - math.prod(1.0 - by[k] for k in keys)
        emp = sum(any(vals_by[k][i] <= -2.5 * band[k]["band_total"] for k in keys) for i in range(5))
        print(f"  {name}: parametric {p_any:.3f}; empirical draws {emp}/5")


if __name__ == "__main__":
    main()
