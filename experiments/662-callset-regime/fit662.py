#!/usr/bin/env python3
"""662 — fit the board's regime for CALL-SET changes from its two call-set pairs (cp = 446bj2zcp − 446bj2z, composition;
K20 = 446bj2zk20 − 446bj2z, size), check it on the delivery pairs and on 618, and price K16 (10-06) and CPH (10-07) with it.

Instrument I-662 (defined here, before either read): H1's real-emitter prices (537's runs, 528's rows_at) at the PowerSurface
level f* where ONE fid multiplier fits both call-set pairs (board/H1 fid of cp = that of K20), with that multiplier for fid
and cp's board/H1 ratio for jac at f*. Band: every grid f where the two pairs' fid multipliers agree within ×1.5, plus the
per-seed roots. nmae: H1 at f* × T, T = 496a's board Δnmae / 496a's own H1 Δnmae at f* (after 10-05), and before 496a reads
the delivery pairs' board / H1 nmae range at f* on both panels (h1-test, 618). reach: unpriced (611 / 619).
Beside it, two reference instruments the reads also grade: "f 1" (the coordinator's reading: H1 at f 1, the K20-anchored fid
multiplier for size, cp's for composition) and 611's f 0.35 card (a995f0e2, already FAILED on K20).
Reads work/regime662.json, delivery662.json, tv662.json, t496_662.json and the board records; writes work/fit662.json.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/662-callset-regime/fit662.py"""
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = HERE / "work"
INV = dict(fid=3.32482, jac=2.72081, reach=1.12927, nmae=1.64787)      # 192's board INV (lib521.INV)
OV = {m: v / 6.0 for m, v in INV.items()}                               # overall per raw unit (nmae: lower is better)
RAW = dict(fid="de_wilcoxon_direction_fidelity_yield_raw", jac="de_wilcoxon_sig_jaccard",
           reach="de_wilcoxon_direction_reach_raw", nmae="de_wilcoxon_lfc_nmae")
SC = dict(fid="score_fid", jac="score_jac", reach="score_reach", nmae="score_nmae", pds="score_pds", mse="score_mse")
E_BASE, E_CP, E_K20 = "be90smZUDcIV7j5MgQ03", "K1Q0y2OFYjTNRqh9jKfk", "kJfe2dP9LEJawDEf0pMm"
LAB = dict(cp="REF cp (528, kept) on 446bj2z", K20="P2 K 2,000 alone (446bj2zk20's construction)",
           K16="P3 K 1,600 alone (446bj2zk16's construction)", P1="P1 CLEAN4 hook on the cp base",
           P8="P8 CLEAN4 hook r 0.17 on the cp base", P9="P9 sign-shuffled CLEAN4 hook on the cp base")
# 611's board values for the two delivery pairs (desk611b.py; 433a − 417e and 446bj2z − 446bj2)
BOARD_DELIV = {"zs (446bj2z-446bj2)": dict(nmae=-0.002179), "coarse MDE (433a-417e)": dict(nmae=0.919618 - 0.936405)}


def ent(i):
    return json.load(open(ROOT / f"data/board-entries/{i}.json"))["entry"]


def main():
    R = json.load(open(W / "regime662.json"))
    D = json.load(open(W / "delivery662.json"))
    T6 = json.load(open(W / "tv662.json"))
    H496 = json.load(open(W / "t496_662.json"))
    C = R["comparisons"]
    grid = [float(f) for f in R["FGRID"]]
    key = dict(zip(grid, R["FGRID"]))
    b0, cp, k20 = ent(E_BASE), ent(E_CP), ent(E_K20)
    B = {p: {m: e[RAW[m]] - b0[RAW[m]] for m in RAW} for p, e in (("cp", cp), ("K20", k20))}
    BOV = {p: {m: (e[SC[m]] - b0[SC[m]]) / 6.0 for m in SC} for p, e in (("cp", cp), ("K20", k20))}
    BOV_tot = {p: e["score_avg"] - b0["score_avg"] for p, e in (("cp", cp), ("K20", k20))}

    def h(p, m, f, seed=None):
        """H1 Δ of pair p, member m, at level f (linear in f between grid points)."""
        def at(g):
            r = C[LAB[p]][key[g]]
            return r["per_seed"][seed][m] if seed is not None else r["d_mean"][m]
        if f in key:
            return at(f)
        i = np.searchsorted(grid, f)
        a, b = grid[i - 1], grid[i]
        t = (f - a) / (b - a)
        return (1 - t) * at(a) + t * at(b)

    # ---- 1. the common-fid-multiplier roots: board_cp / H1_cp(f) == board_K20 / H1_K20(f)
    def roots(seed=None):
        out = []
        for a, b in zip(grid[:-1], grid[1:]):
            ca, cb, ka, kb = h("cp", "fid", a, seed), h("cp", "fid", b, seed), h("K20", "fid", a, seed), h("K20", "fid", b, seed)
            den = B["cp"]["fid"] * (kb - ka) - B["K20"]["fid"] * (cb - ca)
            if den == 0:
                continue
            t = (B["K20"]["fid"] * ca - B["cp"]["fid"] * ka) / den
            if 0 <= t <= 1:
                f = a + t * (b - a)
                c, k = h("cp", "fid", f, seed), h("K20", "fid", f, seed)
                if c < 0 and k < 0:
                    out.append(dict(f=f, m_fid=B["cp"]["fid"] / c))
        return out

    root_mean = roots()
    root_seed = {s: roots(s) for s in range(3)}
    assert root_mean, "no common-multiplier level"
    fstar = root_mean[0]["f"] if len(root_mean) == 1 else float(np.median([r["f"] for r in root_mean]))
    mstar = B["cp"]["fid"] / h("cp", "fid", fstar)

    # ---- 2. per-level table and the consistency band (both fid ratios > 0 and within x1.5 of each other)
    table, band = [], []
    for f in grid:
        row = dict(f=f)
        for p in ("cp", "K20"):
            for m in ("fid", "jac", "nmae"):
                v = h(p, m, f)
                row[f"{p}_{m}"] = v
                row[f"{p}_{m}_ratio"] = B[p][m] / v if abs(v) > 1e-12 else None
        for name, bd in BOARD_DELIV.items():
            row[f"h1_{name}_nmae_ratio"] = bd["nmae"] / D[name][key[f]]["nmae"]
            row[f"618_{name}_nmae_ratio"] = bd["nmae"] / T6[[k for k in T6 if k.startswith(name.split(" ")[0])][0]]["regimes"][key[f]]["nmae"]
        cpn = [k for k in T6 if k.startswith("cp")][0]
        row["618_cp_fid_ratio"] = B["cp"]["fid"] / T6[cpn]["regimes"][key[f]]["fid"]
        row["618_cp_jac_ratio"] = B["cp"]["jac"] / T6[cpn]["regimes"][key[f]]["jac"]
        # signs of the priced cells (|board| >= 1e-4 raw): cp fid / jac, K20 fid / jac, zs nmae, MDE nmae (h1-test)
        cells = [(B["cp"]["fid"], row["cp_fid"]), (B["cp"]["jac"], row["cp_jac"]), (B["K20"]["fid"], row["K20_fid"]),
                 (B["K20"]["jac"], row["K20_jac"])] + [(bd["nmae"], D[n][key[f]]["nmae"]) for n, bd in BOARD_DELIV.items()]
        row["sign_misses"] = int(sum(1 for b_, v in cells if np.sign(b_) != np.sign(v)))
        rc, rk = row["cp_fid_ratio"], row["K20_fid_ratio"]
        if rc and rk and rc > 0 and rk > 0 and abs(math.log(rc / rk)) <= math.log(1.5):
            band.append(f)
        table.append(row)
    fband = sorted(set(band) | {r["f"] for rs in root_seed.values() for r in rs} | {fstar})
    flo, fhi = min(fband), max(fband)

    # ---- 3. delivery check at f*
    deliv = {f"h1-test {n}": BOARD_DELIV[n]["nmae"] / (lambda g: np.interp(fstar, grid, [D[n][key[x]]["nmae"] for x in grid]))(0)
             for n in BOARD_DELIV}
    for n in BOARD_DELIV:
        k6 = [k for k in T6 if k.startswith(n.split(" ")[0])][0]
        deliv[f"618 {n}"] = BOARD_DELIV[n]["nmae"] / np.interp(fstar, grid, [T6[k6]["regimes"][key[x]]["nmae"] for x in grid])
    t_band = (min(deliv.values()), max(deliv.values()))
    t_mid = math.sqrt(t_band[0] * t_band[1])
    h496 = lambda f: float(np.interp(f, grid, [H496[key[x]]["d_mean"]["nmae"] for x in grid]))  # noqa: E731

    # ---- 4. K16 (10-06) under I-662, the f 1 reading, and 611's f 0.35 card
    def k16_at(f, fid_mult, jac_mult):
        d = dict(fid=fid_mult * h("K16", "fid", f), jac=jac_mult * h("K16", "jac", f), nmae=h("K16", "nmae", f))
        d["fidjac_ov"] = d["fid"] * OV["fid"] + d["jac"] * OV["jac"]
        return d
    jstar = B["cp"]["jac"] / h("cp", "jac", fstar)
    k16 = dict(point=k16_at(fstar, mstar, jstar))
    sweep = [k16_at(f, B["K20"]["fid"] / h("K20", "fid", f), B["cp"]["jac"] / h("cp", "jac", f)) for f in fband]
    k16["band"] = {m: [min(s[m] for s in sweep), max(s[m] for s in sweep)] for m in ("fid", "jac", "nmae", "fidjac_ov")}
    k16["f1"] = k16_at(1.0, B["K20"]["fid"] / h("K20", "fid", 1.0), B["cp"]["jac"] / h("cp", "jac", 1.0))
    k16["card611"] = dict(fid=0.0069, jac=0.0010, fidjac_ov=0.0043)

    # ---- 5. CPH (10-07): board(CPH − 446bj2zcp) = the hook-on-block increment; board(CPH − 446bj2z) = cp's read + increment
    def cph_at(f, src, m_fid, m_jac, T):
        inc = {m: h(src, m, f) for m in ("fid", "jac", "nmae", "reach")}
        fid, jac, nmae = inc["fid"] * m_fid, inc["jac"] * m_jac, inc["nmae"] * T
        prim = fid * OV["fid"] - nmae * OV["nmae"]
        ship = BOV["cp"]["fid"] + BOV["cp"]["jac"] + fid * OV["fid"] + jac * OV["jac"]
        tot = BOV_tot["cp"] + fid * OV["fid"] + jac * OV["jac"] - nmae * OV["nmae"]
        return dict(fid=fid, jac=jac, nmae=nmae, h1_reach=inc["reach"], primary_vs_cp=prim, fidjac_vs_446bj2z=ship,
                    overall_vs_446bj2z_ex_inc_reach=tot)

    def cph_instrument(fpt, fset, mf, mj, tband, tmid):
        pt = {src: cph_at(fpt, src, mf(fpt), mj(fpt), tmid) for src in ("P8", "P1", "P9")}
        lo_hi = [cph_at(f, src, mf(f), mj(f), T) for f in fset for src in ("P1", "P8", "P9") for T in tband]
        rng = {k: [min(x[k] for x in lo_hi), max(x[k] for x in lo_hi)] for k in ("primary_vs_cp", "fidjac_vs_446bj2z",
                                                                                    "overall_vs_446bj2z_ex_inc_reach")}
        rng8 = [cph_at(f, "P8", mf(f), mj(f), T) for f in fset for T in tband]
        rng["P8_primary_vs_cp"] = [min(x["primary_vs_cp"] for x in rng8), max(x["primary_vs_cp"] for x in rng8)]
        return dict(f=fpt, T_band=list(tband), T_mid=tmid, point=pt, range=rng)

    cpfid = lambda f: B["cp"]["fid"] / h("cp", "fid", f)  # noqa: E731
    cpjac = lambda f: B["cp"]["jac"] / h("cp", "jac", f)  # noqa: E731
    deliv1 = {n: BOARD_DELIV[n]["nmae"] / D[n]["1.0"]["nmae"] for n in BOARD_DELIV}
    for n in BOARD_DELIV:
        k6 = [k for k in T6 if k.startswith(n.split(" ")[0])][0]
        deliv1["618 " + n] = BOARD_DELIV[n]["nmae"] / T6[k6]["regimes"]["1.0"]["nmae"]
    t1 = (min(deliv1.values()), max(deliv1.values()))
    cph = {"I-662": cph_instrument(fstar, fband, lambda f: mstar if f == fstar else cpfid(f), cpjac, t_band, t_mid),
           "f1": cph_instrument(1.0, [1.0], cpfid, cpjac, t1, math.sqrt(t1[0] * t1[1]))}
    cph["T_after_496a"] = dict(I662=f"T = (496a board dnmae raw) / ({h496(fstar):+.6f})", f1=f"T = (496a board dnmae raw) / ({h496(1.0):+.6f})",
                               h1_496a_fstar=h496(fstar), h1_496a_f1=h496(1.0))

    out = dict(board=B, board_overall=BOV, board_total=BOV_tot, roots_mean=root_mean, roots_seed=root_seed, fstar=fstar,
               m_fid_star=mstar, m_jac_star=jstar, band_f=fband, delivery_nmae_ratio_at_fstar=deliv, T_band=t_band, T_mid=t_mid,
               delivery_nmae_ratio_f1=deliv1, table=table, K16=k16, CPH=cph)
    (W / "fit662.json").write_text(json.dumps(out, indent=1))

    print(f"board raw  cp: " + " ".join(f"{m} {B['cp'][m]:+.6f}" for m in RAW) + " | K20: " + " ".join(f"{m} {B['K20'][m]:+.6f}" for m in RAW))
    print(f"common-fid-multiplier roots (mean): {[(round(r['f'], 4), round(r['m_fid'], 3)) for r in root_mean]}")
    for s, rs in root_seed.items():
        print(f"   seed {s}: {[(round(r['f'], 4), round(r['m_fid'], 3)) for r in rs]}")
    print(f"f* = {fstar:.4f}  m_fid* = {mstar:.3f}  m_jac* (cp) = {jstar:.3f}  K20 jac ratio at f* = {B['K20']['jac'] / h('K20', 'jac', fstar):.3f}"
          f"  K20 nmae H1 at f* = {h('K20', 'nmae', fstar):+.5f} (board {B['K20']['nmae']:+.6f})")
    print(f"consistency band (fid multipliers within x1.5, + seed roots): f {flo:.2f} .. {fhi:.2f}  ({fband})")
    print("delivery nmae board/H1 at f*: " + ", ".join(f"{k} {v:.2f}" for k, v in deliv.items()) + f" -> T band {t_band[0]:.2f}..{t_band[1]:.2f}")
    print("\n f    | cp fid  ratio | cp jac ratio | K20 fid ratio | K20 jac ratio | K20 nmae | zs/MDE nmae ratio h1 | 618 | sign misses")
    for r in table:
        if r["f"] in (0.35, 0.4, 0.45, 0.5, 0.55, 0.58, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.67, 0.7, 0.75, 0.8, 0.9, 1.0):
            n = list(BOARD_DELIV)
            print(f" {r['f']:.2f} | {r['cp_fid']:+.5f} {r['cp_fid_ratio']:5.2f} | {r['cp_jac_ratio']:5.2f} | {r['K20_fid']:+.5f} {r['K20_fid_ratio']:6.2f} | "
                  f"{r['K20_jac']:+.5f} {r['K20_jac_ratio'] if r['K20_jac_ratio'] is not None else float('nan'):6.2f} | {r['K20_nmae']:+.5f} | "
                  f"{r['h1_' + n[0] + '_nmae_ratio']:.2f}/{r['h1_' + n[1] + '_nmae_ratio']:.2f} | {r['618_' + n[0] + '_nmae_ratio']:.2f}/{r['618_' + n[1] + '_nmae_ratio']:.2f} "
                  f"(cp fid {r['618_cp_fid_ratio']:.2f} jac {r['618_cp_jac_ratio']:.2f}) | {r['sign_misses']}")
    print("\nK16 vs 446bj2z (raw; fid+jac in overall units):")
    for k in ("point", "f1", "card611"):
        print(f"  {k:8s} " + " ".join(f"{m} {v:+.5f}" for m, v in k16[k].items()))
    print("  band     " + " ".join(f"{m} [{v[0]:+.5f}, {v[1]:+.5f}]" for m, v in k16["band"].items()))
    print("\nCPH:")
    for name in ("I-662", "f1"):
        c = cph[name]
        print(f"  {name}: f {c['f']:.3f}, T band {c['T_band'][0]:.2f}..{c['T_band'][1]:.2f} (mid {c['T_mid']:.2f})")
        for src, v in c["point"].items():
            print(f"    {src}: " + " ".join(f"{k} {x:+.5f}" for k, x in v.items()))
        print("    range: " + " ".join(f"{k} [{v[0]:+.5f}, {v[1]:+.5f}]" for k, v in c["range"].items()))
    print("  " + str(cph["T_after_496a"]))


if __name__ == "__main__":
    main()
