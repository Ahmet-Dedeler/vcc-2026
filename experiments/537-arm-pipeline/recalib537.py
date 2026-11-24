#!/usr/bin/env python3
"""537 — re-price 446bj2zcph (the block + CLEAN4 hook) and 446bj2zmzcp from the board reads as they land (10-01 cp, 10-05 hook).

usage: uv run python experiments/537-arm-pipeline/recalib537.py [--cp <446bj2zcp entry>] [--h <496a entry>]
Reads data/board-entries/<id>.json and work/price537.json (the H1 prices, 3 seeds). Prints:
  - the block's board read vs 446bj2z by member, and its ratio to H1's block-alone price (f 0.35 and cmp, zs form);
  - the hook's board factor T = 496a's board d nmae / H1's hook-alone d nmae on 446bj2's construction (-0.00146, 496);
  - CPH's predicted board row vs 446bj2z = the block's read + the H1 hook-on-block increment scaled member by member:
    fid by the block's own fid ratio (bounded to [0.07, 0.67] when the block's fid read is noise-sized), jac by the block's
    jac ratio, nmae by T (or 0.35 ... 1.76 before 496a), reach reported at x0.
Nothing here is a bin; the preregistered bins stay in PREREGISTRATION-draft-CPH.md / -CPMZ.md.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
INV = dict(fid=3.32482, jac=2.72081, reach=1.12927, nmae=1.64787)
Z, BJ2 = "be90smZUDcIV7j5MgQ03", "YvVru3DAF89Jiyr34olR"
H1_HOOK_496 = -0.00146


def entry(i):
    d = json.loads((ROOT / f"data/board-entries/{i}.json").read_text())
    return d.get("entry", d)


def ov(m, v):
    return v * INV[m] / 6 * (-1 if m == "nmae" else 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cp")
    ap.add_argument("--h")
    a = ap.parse_args()
    P = json.loads((HERE / "work/price537.json").read_text())["comparisons"]
    blk = P["REF cp (528, kept) on 446bj2z"]
    inc = P["P1 CLEAN4 hook on the cp base"]
    regs = ("0.35", "cmp")
    out = {}
    if a.cp:
        e, z = entry(a.cp), entry(Z)
        d = {m: e[RAW[m]] - z[RAW[m]] for m in RAW}
        out["cp_read"] = d
        print(f"446bj2zcp - 446bj2z (board, raw): " + " ".join(f"{m} {d[m]:+.6f}" for m in RAW)
              + f" | overall {e['score_avg'] - z['score_avg']:+.7f}")
        # 535's preregistered bins (528's, unchanged): primary d(fid + jac) overall vs 446bj2z
        prim = d["fid"] * INV["fid"] / 6 + d["jac"] * INV["jac"] / 6
        b = ("POSITIVE" if prim >= 0.0012 else "PARTIAL" if prim >= 0.0004 else "NULL" if prim > -0.0004 else "NEGATIVE")
        rider = b == "PARTIAL" and d["jac"] >= 0.002 and d["nmae"] <= 0.001
        rescue = b in ("NULL", "NEGATIVE") and d["jac"] >= 0.002 and d["fid"] < 0
        fals = [f for f, c in (("pds/mse moved > 1e-6 raw", abs(d["pds"]) > 1e-6 or abs(d["mse"]) > 1e-6),
                               ("d jac <= 0 (535's falsifier if d fid is at or above its prediction: check by hand)", d["jac"] <= 0), ("|d nmae| > 0.003 raw", abs(d["nmae"]) > 0.003)) if c]
        out["cp_bin"] = dict(primary=prim, bin=b, partial_rider=rider, rescue=rescue, falsifiers=fals,
                             cph_10_07=b in ("POSITIVE", "PARTIAL") or rescue)
        print(f"   535 primary d(fid + jac) overall {prim:+.5f} -> {b}" + (" (October rider: yes)" if rider else "")
              + (" | RESCUE flag (537 §0 C): the composition transferred, the unsteered signs cost fid" if rescue else "")
              + (f" | FALSIFIER: {fals}" if fals else "") + f" | 10-07 s1 = CPH: {out['cp_bin']['cph_10_07']}")
        for reg in regs:
            h = blk[reg]["d_mean"]
            print(f"   ratio board / H1 block-alone at f {reg}: " + " ".join(
                f"{m} {d[m] / h[m]:+.2f}" if abs(h[m]) > 1e-9 else f"{m} n/a" for m in ("fid", "jac", "reach", "nmae")))
    T = None
    if a.h:
        e, b = entry(a.h), entry(BJ2)
        dn = e[RAW["nmae"]] - b[RAW["nmae"]]
        T = dn / H1_HOOK_496
        out["T"] = T
        print(f"496a - 446bj2: d nmae raw {dn:+.6f} -> hook board factor T = {T:.2f} (H1 {H1_HOOK_496})")
    for reg in regs:
        hi = inc[reg]["d_mean"]
        Ts = [T] if T is not None else [0.35, 1.76]
        if a.cp:
            h = blk[reg]["d_mean"]
            rj = out["cp_read"]["jac"] / h["jac"] if abs(h["jac"]) > 1e-9 else 0.72
            rf = out["cp_read"]["fid"] / h["fid"] if abs(h["fid"]) > 1e-9 else 0.3
            rf = min(max(rf, 0.07), 0.67) if abs(out["cp_read"]["fid"]) < 0.0003 else rf
        else:
            rj, rf = 0.72, None
        for t in Ts:
            fids = [hi["fid"] * f for f in ((0.07, 0.67) if rf is None else (rf,))]
            for fv in fids:
                inc_ov = ov("fid", fv) + ov("jac", hi["jac"] * rj) + ov("nmae", hi["nmae"] * t)
                base = sum(ov(m, out["cp_read"][m]) for m in ("fid", "jac", "reach", "nmae")) if a.cp else None
                print(f"f {reg}: T {t:.2f}, fid factor {fv / hi['fid'] if hi['fid'] else 0:.2f}: hook-on-block increment "
                      f"(fid + jac + nmae) {inc_ov:+.5f} overall" + (f"; CPH vs 446bj2z (DE, reach x1 as read) "
                                                                        f"{base + inc_ov:+.5f}" if base is not None else ""))
    (HERE / "work/recalib537.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
