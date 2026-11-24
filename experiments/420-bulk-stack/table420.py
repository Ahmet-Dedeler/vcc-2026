"""420: render the RESULTS tables from work/analyse420.json, abc420.json and abc420_leak.json (no hand transcription).

    uv run python table420.py > work/table420.md
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
A = json.loads((WORK / "analyse420.json").read_text())
ABC = json.loads((WORK / "abc420.json").read_text())["designs"]
LEAK = json.loads((WORK / "abc420_leak.json").read_text())
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
REF = "PP|K3+sh1.25"


def esc(name):
    """Design names carry '|' (composition|key); escape it inside markdown table cells."""
    return name.replace("|", "\\|")


def ff(v, n=4):
    return " / ".join(f"{x:+.{n}f}" for x in v)


def row_stats(name, key):
    r = A["rows"][name][key]
    v = [r[f]["ov"] for f in CANCER]
    return v, float(np.mean(v)), sum(x > 0 for x in v)


def eratio(name):
    return [ABC[name.split(":")[0].replace("V1own", "V1") if not name.startswith("CTRL") else "PP|V1+sh1.25"]["realised"]["16"][c]["A"]
            / ABC[REF]["realised"]["16"][c]["A"] for c in "ABC"]


def flag(er):
    return "**HI**" if max(er) > 1.03 else ("**LO**" if min(er) < 0.90 else "")


PS = {"K3": "wshare_lo0.3_hi3", "K5_3": "wshare_lo0.5_hi3", "V1": "lodo415_lo0.3_hi3", "V1_5_3": "lodo415_lo0.5_hi3*",
      "V1own": "<own-composition LODO>*", "Kinf": "wshare_lo0.3_inf"}
CENV = {"PP": None, "h1x2": "H1_G .56", "cd4x1.5": "CD4_F .51", "stack": "H1_G .56, CD4_F .51", "stack_h1x1.5": "H1_G .42, CD4_F .51",
        "stack_cd4x1.25": "H1_G .56, CD4_F .425"}


def short_builder(n, spec=None):
    """The env change vs 417a (Q 1.25, PS wshare_lo0.3_hi3, PP composition). * = artifact not yet wired in the builder."""
    if n.startswith("CTRL"):
        return "-"
    comp, rest = n.split(":")[0].split("|")
    key, q = (rest.split("+sh") + [None])[:2]
    parts = []
    if q and q != "1.25":
        parts.append(f"Q={q}")
    if not q:
        parts.append("Q unset")
    if CENV[comp]:
        parts.append(f"LOAD[{CENV[comp]}; PS {PS[key]}]")
    elif key != "K3":
        parts.append(f"PS={PS[key]}")
    if spec:
        parts.append(f"CTX={spec}")
    return "`" + " · ".join(parts) + "`" if parts else "(417a)"


def ranked_row(name, label=None):
    k1v, k1m, k1n = row_stats(name, "r417a_k1")
    kkv, kkm, kkn = row_stats(name, "r417a_k1.41")
    m = lambda fld, key="r417a_k1": float(np.mean([A["rows"][name][key][f][fld] for f in CANCER]))
    d1, dk = row_stats(name, "d418_k1")[1], row_stats(name, "d418_k1.41")[1]
    kf = A["rows"][name]["r417a_k1"].get("kolf2", {}).get("ov")
    er = eratio(name)
    b = short_builder(name)
    return (f"| {esc(label or name)} | {m('d_rel'):+.4f} | {m('dpds_bd'):+.5f} | {m('dmse'):+.4f} / {m('dmse', 'r417a_k1.41'):+.4f} | "
            f"{ff(k1v)} | **{k1m:+.5f}** {k1n}/3 | {ff(kkv)} | **{kkm:+.5f}** {kkn}/3 | {d1:+.5f} / {dk:+.5f} | "
            f"{'/'.join(f'{x:.3f}' for x in er)} {flag(er)} | {kf:+.4f} | {b} |")


HDR = ("| design | Δrel cos | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | **k=1 mean** n+ | Δov k=1.407 per fold | "
       "**k=1.407 mean** n+ | 418 route k=1 / 1.407 | A/B/C energy vs 417a | kolf2 k=1 | builder change vs 417a |\n"
       "|---|---:|---:|---|---|---|---|---|---|---|---:|---|")

print("### T-val. Energy-route validation on the shipped artifacts (A/B/C)\n")
print("| ratio | measured | 418 design route (err) | realised K16 route (err) |")
print("|---|---|---|---|")
for lab, v in A["validation"].items():
    me, de, re_ = v["measured"], v["design"], v["realised"]
    print(f"| {lab} | {' / '.join(f'{x:.4f}' for x in me)} | {' / '.join(f'{x:.4f}' for x in de)} ({ff([d - m for d, m in zip(de, me)])}) | "
          f"{' / '.join(f'{x:.4f}' for x in re_)} ({ff([d - m for d, m in zip(re_, me)])}) |")
print()
print("### T-a. The stack and its halves (vs 417a; q 1.25, key K3)\n")
print(HDR)
for n in ("stack|K3+sh1.25", "h1x2|K3+sh1.25", "cd4x1.5|K3+sh1.25", "stack_h1x1.5|K3+sh1.25", "stack_cd4x1.25|K3+sh1.25"):
    print(ranked_row(n))
print("\n### T-b. Sharpening dose on PP|K3 (and PP|V1), vs 417a\n")
print(HDR)
for n in ("PP|K3+sh1.35", "PP|K3+sh1.5", "PP|K3+sh1.75", "PP|K3+sh2", "PP|V1+sh1.5", "PP|V1+sh1.75"):
    print(ranked_row(n))
print("\n### T-c. Key variants on the sharpened bulk (q 1.25), vs 417a\n")
print(HDR)
for n in ("PP|V1+sh1.25", "PP|K5_3+sh1.25", "PP|V1_5_3+sh1.25", "CTRL:PP|V1shuf+sh1.25"):
    print(ranked_row(n))
print("\n### T-d. Best combinations found (top 15 of the 150-design grid by k=1 mean), vs 417a (no re-levelling)\n")
print(HDR)
for n in [x for x in A["ranked"] if not x.startswith("CTRL")][:15]:
    print(ranked_row(n))
print("\n### T-d2. Best combinations restricted to q <= 1.5 (energy within [0.90, 1.03] of 417a)\n")
print(HDR)
sel = [x for x in A["ranked"] if not x.startswith("CTRL") and "+sh" in x and float(x.split("+sh")[1]) <= 1.5
       and min(eratio(x)) >= 0.90 and max(eratio(x)) <= 1.03][:10]
for n in sel:
    print(ranked_row(n))

print("\n### T-lev. Per-context re-levelling on top (VCC_BULK_CTX_SCALE, EV-optimal; pds measured on PLOLO for the scaled design)\n")
print("| design | factors (k=1) | render scale A/B/C | level-only Δov per fold | level-only mean | Δraw mse | Δpds_bd | design + level vs 417a, k=1 per fold | mean n+ | k=1.407 mean |")
print("|---|---|---|---|---:|---:|---:|---|---|---:|")
for k, v in A["level"].items():
    if not k.endswith("|k1"):
        continue
    d = k[:-3]
    p = v["per"]
    vK = A["level"].get(f"{d}|kK", {}).get("per", {})
    lo = [p[f]["ov_level_only"] for f in CANCER]
    tot = [p[f]["ov_total"] for f in CANCER]
    totK = float(np.mean([vK[f]["ov_total"] for f in CANCER])) if vK else float("nan")
    meas = all(p[f]["pds_measured"] for f in CANCER)
    print(f"| {esc(d)} | `{v['spec'].split('=', 1)[1]}` | {'/'.join(f'{v['render_scale'][c]:.4f}' for c in 'ABC')} | {ff(lo)} | "
          f"{np.mean(lo):+.5f} | {np.mean([p[f]['dmse_level'] for f in CANCER]):+.4f} | "
          f"{np.mean([p[f]['dpds_level'] for f in CANCER]):+.5f}{'' if meas else ' (not measured)'} | {ff(tot)} | "
          f"**{np.mean(tot):+.5f}** {sum(x > 0 for x in tot)}/3 | {totK:+.5f} |")

print("\n### T-kmax. Ballast realisability and the KMAX 24 / 32 emulation (r417a route, cancer means; DE leak in units of 410's full weld)\n")
print("| design | capped K16 / K24 / K32 | cos(realised, design) q1 (K16) | Δov vs 417a k=1: K16 / K24 / K32 | k=1.407: K16 / K24 / K32 | "
      "DE leak K16 / K24 / K32 | mean K (K16 / K32) |")
print("|---|---|---:|---|---|---|---|")
for n in ("PP|K3", "PP|K3+sh1.25", "PP|K3+sh1.35", "PP|K3+sh1.5", "PP|K3+sh1.75", "PP|K3+sh2", "PP|V1+sh1.25", "PP|V1+sh1.5",
          "stack|K3+sh1.25", "stack|V1+sh1.5"):
    b = ABC[n]["ballast"]
    m = lambda key: float(np.mean([A["rows"][n][key][f]["ov"] for f in CANCER]))
    lk = LEAK.get(n, {})
    leak = " / ".join(f"{lk[k]['leak_vs_410weld']:.4f}" for k in ("16", "24", "32")) if lk else "-"
    mk = f"{lk['16']['mean_K']:.2f} / {lk['32']['mean_K']:.2f}" if lk else "-"
    print(f"| {esc(n)} | {100 * b['16']['frac_capped']:.1f} / {100 * b['24']['frac_capped']:.1f} / {100 * b['32']['frac_capped']:.1f} % | "
          f"{b['16']['cos_realised_vs_design_q'][0]:.4f} | {m('r417a_k1'):+.5f} / {m('r417a_K24_k1'):+.5f} / {m('r417a_K32_k1'):+.5f} | "
          f"{m('r417a_k1.41'):+.5f} / {m('r417a_K24_k1.41'):+.5f} / {m('r417a_K32_k1.41'):+.5f} | {leak} | {mk} |")

print("\n### T-int. Interactions: each lever vs the same q without it (Δov k=1 / k=1.407, cancer means; all 3/3 unless marked)\n")
qs = ("1.25", "1.35", "1.5", "1.75", "2")
print("| lever | " + " | ".join(f"q {q}" for q in qs) + " |")
print("|---|" + "---|" * len(qs))
for lab, v in A["interactions"].items():
    cells = []
    for q in qs:
        k1, kK = v[q]["k1"], v[q]["kK"]
        n = sum(x > 0 for x in k1)
        cells.append(f"{np.mean(k1):+.5f} / {np.mean(kK):+.5f}" + ("" if n == 3 else f" ({n}/3)"))
    print(f"| {lab} | " + " | ".join(cells) + " |")

print("\n### T-seq. Sequences (each arm one variable from the previous)\n")
for sname, sq in A["sequences"].items():
    print(f"**{sname}** — cumulative vs 417a: k=1 {sq['cumulative']['k1']:+.5f}, k=1.407 {sq['cumulative']['kK']:+.5f}, "
          f"class-calibrated {sq['cumulative']['cls']:+.5f}\n")
    print("| step | Δrel | Δpds_bd | Δraw mse k=1 / 1.407 | Δov k=1 per fold | k=1 mean n+ | k=1.407 | class-cal. | kolf2 k=1 | builder |")
    print("|---|---:|---:|---|---|---|---:|---:|---:|---|")
    for r in sq["steps"]:
        print(f"| {r['step']} | {r['d_rel']:+.4f} | {r['dpds_bd']:+.5f} | {r['dmse_k1']:+.4f} / {r['dmse_kK']:+.4f} | "
              f"{ff(r['per_fold_k1'])} | **{r['k1']:+.5f}** {r['n_pos_k1']}/3 | {r['kK']:+.5f} | {r['cls_ov']:+.5f} | "
              f"{r['kolf2_k1']:+.4f} | {short_builder(r['arm'], r['arm'].split(':')[1] if ':' in r['arm'] else None)} |")
    print()

print("### T-read. Reading 417a's board row: implied sharpening transfer, and the q steps' worth at that transfer\n")
R = A["read417a"]
print(f"PLOLO (realised, cancer mean) 417a vs 410b: rel cos {R['r_sh']:+.4f}, board-emulated pds {R['dp_sh']:+.5f}; measured "
      f"render amplitude 417a/410b {R['amp']:.4f}. k_pds is tied to k_sh at the kolf2 ratio (x0.895).\n")
print("| k_sh | predicted 417a raw mse | predicted 417a raw pds | Δoverall 417a vs 410b (bulk members) | q 1.35 | q 1.5 | q 1.75 | q 2 "
      "| levelled q 1.35 | levelled q 1.5 | levelled q 1.75 | levelled q 2 |")
print("|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
for r, rl in zip(R["rows"], R.get("levelled", [{}] * len(R["rows"]))):
    lq = rl.get("q_steps", {})
    print(f"| {r['k_sh']:.3f} | {r['raw_mse']:.4f} | {r['raw_pds']:.4f} | {r['d_overall']:+.4f} | "
          + " | ".join(f"{r['q_steps'][q]:+.5f}" for q in ("1.35", "1.5", "1.75", "2")) + " | "
          + " | ".join(f"{lq[q]:+.5f}" if q in lq else "-" for q in ("1.35", "1.5", "1.75", "2")) + " |")

print("\n### T-one. One-variable steps from 417a with class-calibrated transfer and bootstrap\n")
print("| step | k=1 per fold | k=1 mean n+ | k=1.407 mean | class-cal. | Δpds_bd (emulated) / literal | ov with literal pds | "
      "bootstrap P(Δrel>0)/P(Δpds_bd>0) k562 · hct116 · hek293t | kolf2 k=1 |")
print("|---|---|---|---:|---:|---|---:|---|---:|")
for lab, v in A["one_step"].items():
    p = v["per"]
    k1 = [p[f]["k1"]["ov"] for f in CANCER]
    kK = float(np.mean([p[f]["kK"]["ov"] for f in CANCER]))
    cl = float(np.mean([p[f]["cls"]["ov"] for f in CANCER]))
    bd = float(np.mean([p[f]["k1"]["dpds_bd"] for f in CANCER]))
    li = float(np.mean([p[f]["k1"]["dpds_lit"] for f in CANCER]))
    ovl = float(np.mean([p[f]["k1"]["ov"] + 2.2138 * (p[f]["k1"]["dpds_lit"] - p[f]["k1"]["dpds_bd"]) / 6 for f in CANCER]))
    bt = " · ".join(f"{v['boot'][f]['p_rel']:.2f}/{v['boot'][f]['p_pds']:.2f}" for f in CANCER)
    kf = p["kolf2"]["k1"]["ov"] if "kolf2" in p else float("nan")
    extra = f" (H1-then-CD4 class path {np.mean([v['cls_path_h1_then_cd4'][f] for f in CANCER]):+.5f})" if "cls_path_h1_then_cd4" in v else ""
    print(f"| {lab} | {ff(k1)} | **{np.mean(k1):+.5f}** {sum(x > 0 for x in k1)}/3 | {kK:+.5f} | {cl:+.5f}{extra} | "
          f"{bd:+.5f} / {li:+.5f} | {ovl:+.5f} | {bt} | {kf:+.4f} |")
