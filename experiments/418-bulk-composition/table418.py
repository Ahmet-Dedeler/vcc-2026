"""418: markdown tables for RESULTS.md from work/analyse418.json (per fold, 3-fold mean, range, sign count).

    uv run python table418.py > work/table418.md
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

WORK = Path(__file__).resolve().parent / "work"
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
KK = "k1.41"
KNOB = {
    "a": ("VCC_CD4_TARGET_F", 0.34, lambda m: f"{0.34 * m:.3g}"),
    "b": ("VCC_H1_UNION_G", 0.28, lambda m: "unset" if m == 0 else f"{0.28 * m:.3g}"),
}


def label(arm):
    fam, rest = arm.split(":", 1)
    if fam == "a":
        m = float(rest.replace("cd4x", ""))
        return f"CD4 f ×{m:g} (`VCC_CD4_TARGET_F={0.34 * m:.3g}`)"
    if fam == "b":
        m = float(rest.replace("h1x", ""))
        return f"H1 g ×{m:g} (`VCC_H1_UNION_G={'unset' if m == 0 else f'{0.28 * m:.3g}'}`)"
    if fam == "c":
        k = rest.replace("k", "")
        return f"CD4 deflation k {k} (`VCC_CD4_DEFLATE_K={'unset' if k == '0' else k}`)"
    if fam == "d":
        g = rest.replace("g", "")
        return f"partial g {g} (`VCC_CD4_PARTIAL_G={'unset' if g == '0' else g}`)"
    if fam == "e":
        x = rest.replace("exp", "")
        return f"normrel cap exp {x} (`VCC_NORMREL_EXP={'unset' if x == '0' else x}`)"
    if fam == "f1":
        return "f1: key on the donor base only (no knob)"
    if fam == "f2":
        return "f2: normrel PLAIN −1.5 (`VCC_NORMREL_MODE=plain`)"
    return arm


def main():
    an = json.loads((WORK / "analyse418.json").read_text())
    for variant, title in (("K3", "vs 410b's bulk (PP + fixed key, no sharpening)"),
                           ("K3+sh", "vs 417a's bulk (PP + fixed key + sharpening q 1.25)")):
        rows = []
        for arm, row in an["rows"].items():
            if variant not in row:
                continue
            per = row[variant]
            v1 = [per[f]["abc_k1"]["ov"] for f in CANCER]
            vK = [per[f]["abc_" + KK]["ov"] for f in CANCER]
            m = lambda key, fld: float(np.mean([per[f][key][fld] for f in CANCER]))
            hist = float(np.mean([per[f]["abc_hist"]["ov"] for f in CANCER])) if "abc_hist" in per[CANCER[0]] else None
            rows.append(dict(arm=arm, v1=v1, vK=vK, mean1=float(np.mean(v1)), meanK=float(np.mean(vK)),
                             rel=m("abc_k1", "d_rel"), dpds=m("abc_k1", "dpds_board"), dmse1=m("abc_k1", "dmse"),
                             dmseK=m("abc_" + KK, "dmse"), plolo=m("k1", "ov"), match=m("abc_matched_k1", "ov"), hist=hist,
                             Ar=[m("abc_k1", "Ar")] if False else [float(np.mean([per[f]["abc_k1"]["Ar"][c] for f in CANCER])) for c in "ABC"],
                             kolf=(per["kolf2"]["abc_k1"]["ov"], per["kolf2"]["abc_" + KK]["ov"]) if "kolf2" in per else None))
        rows.sort(key=lambda r: -r["mean1"])
        print(f"\n### Ranked {title}\n")
        print("Δoverall = (2.2138 Δpds_bd − 1.051 Δmse)/6. A/B/C-energy cliff route; board-emulated pds. "
              "Per fold k562 / hct116 / hek293t.\n")
        print("| # | arm | Δrel cos | Δpds_bd | Δmse (k=1 / 1.407) | **k=1** per fold | **k=1** mean [range] +n | "
              "**k=1.407** per fold | **k=1.407** mean [range] +n | hist | amp-matched | PLOLO-energy k=1 | A/B/C render A ratio | kolf2 k=1 / 1.407 |")
        print("|---:|---|---:|---:|---|---|---|---|---|---:|---:|---:|---|---|")
        for i, r in enumerate(rows, 1):
            f1 = " / ".join(f"{x:+.4f}" for x in r["v1"])
            fK = " / ".join(f"{x:+.4f}" for x in r["vK"])
            print(f"| {i} | {label(r['arm'])} | {r['rel']:+.4f} | {r['dpds']:+.5f} | {r['dmse1']:+.4f} / {r['dmseK']:+.4f} | "
                  f"{f1} | **{r['mean1']:+.5f}** [{min(r['v1']):+.4f}, {max(r['v1']):+.4f}] {sum(x > 0 for x in r['v1'])}/3 | "
                  f"{fK} | **{r['meanK']:+.5f}** [{min(r['vK']):+.4f}, {max(r['vK']):+.4f}] {sum(x > 0 for x in r['vK'])}/3 | "
                  f"{('%+.5f' % r['hist']) if r['hist'] is not None else '–'} | {r['match']:+.5f} | {r['plolo']:+.5f} | "
                  f"{r['Ar'][0]:.4f}/{r['Ar'][1]:.4f}/{r['Ar'][2]:.4f} | "
                  + (f"{r['kolf'][0]:+.4f} / {r['kolf'][1]:+.4f}" if r["kolf"] else "–") + " |")
    # secondary: own LODO-V1 key
    print("\n### Secondary: each arm with its OWN LODO-V1 key recomputed on its composition (vs PP + V1, and + sharpening)\n")
    print("| arm | vs PP+V1: k=1 mean [range] +n | k=1.407 mean | vs PP+V1+sh: k=1 mean [range] +n | k=1.407 mean |")
    print("|---|---|---:|---|---:|")
    rows = []
    for arm, row in an["rows"].items():
        if "V1" not in row:
            continue
        a = [row["V1"][f]["abc_k1"]["ov"] for f in CANCER]
        b = [row["V1+sh"][f]["abc_k1"]["ov"] for f in CANCER]
        aK = float(np.mean([row["V1"][f]["abc_" + KK]["ov"] for f in CANCER]))
        bK = float(np.mean([row["V1+sh"][f]["abc_" + KK]["ov"] for f in CANCER]))
        rows.append((float(np.mean(b)), f"| {label(arm)} | {np.mean(a):+.5f} [{min(a):+.4f}, {max(a):+.4f}] {sum(x > 0 for x in a)}/3 | "
                                         f"{aK:+.5f} | {np.mean(b):+.5f} [{min(b):+.4f}, {max(b):+.4f}] {sum(x > 0 for x in b)}/3 | {bK:+.5f} |"))
    for _, line in sorted(rows, key=lambda t: -t[0]):
        print(line)
    # retro
    for key, title in (("retro_abc", "cancer folds, A/B/C-panel energies"), ("retro_plolo", "cancer folds, PLOLO energies"),
                       ("retro_kolf2", "kolf2 only (3 production donors), A/B/C-panel energies")):
        print(f"\n### Retrodiction of the board's historical composition steps ({title})\n")
        print("| step | board Δpds | PLOLO Δpds | board Δmse | pred Δmse k=1 | pred k=1.407 | k_fit (mse) | "
              "board Δ(pds+mse) overall | pred k=1 |")
        print("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
        for st, v in an[key].items():
            bo = (2.2138 * v["board_dpds"] - 1.051 * v["board_dmse"]) / 6
            pr = (2.2138 * v["plolo_dpds"] - 1.051 * v["pred_dmse_k1"]) / 6
            print(f"| {v['what']} | {v['board_dpds']:+.5f} | {v['plolo_dpds']:+.5f} | {v['board_dmse']:+.5f} | "
                  f"{v['pred_dmse_k1']:+.5f} | {v['pred_dmse_kkey']:+.5f} | {v['k_fit']:.2f} | {bo:+.5f} | {pr:+.5f} |")
        s = an[key + "_summary"]
        print(f"\nsign agreement: overall {s['sign_overall']}/{s['n']}, pds {s['sign_pds']}/{s['n']}, mse {s['sign_mse']}/{s['n']}; "
              f"corr overall {s['corr_overall']:+.3f}, pds {s['corr_pds']:+.3f}, mse {s['corr_mse']:+.3f}")
    print("\nclass constants (A/B/C energies): " + "; ".join(
        f"{k}: k_pds {v['k_pds']:.2f}, k_mse {v['k_mse']:.2f}" for k, v in an["hist_constants"]["abc"].items()))


if __name__ == "__main__":
    main()
