"""250 — the tables: known-answer rows, the two-instrument decomposition on A/B/C, 233's
retention reproduction, the k ladder + plateau test per panel, every control, leave-top-3-out,
the pinned A/B/C basis, the `+partial` composition, the six-member price (246 §4's convention)
and the per-panel verdict.

    uv run python analyse250.py            -> work/analyse250.json, stdout (tee to analyse250.log)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
from price218 import BOARD_COS, RATE, SPAN_PDS, SIGMA_PDS, GATE_SHAPE, GATE_AMPL  # noqa: E402

PANELS = ("abc", "real", "pess", "efp", "efr")
DISJOINT = ("real", "pess", "efp", "efr")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
KS = (2, 3, 5, 7, 10, 15)
BASE = "BASE"
PLATEAU = 0.90          # k is "inside the plateau" if d_pds(k) >= 0.9 x max_k d_pds
R233 = ROOT / "experiments/233-october-retention/work"
R246 = ROOT / "experiments/246-direction-hypotheses/work"


def load():
    rows, pt, meta = {}, {}, {}
    for f in sorted(WORK.glob("sweep250-*-*.json")):
        d = json.loads(f.read_text())
        panel, held = d["meta"]["panel"], d["meta"]["fold"]
        meta[(panel, held)] = d["meta"]
        for r in d["rows"]:
            rows[(panel, held, r["arm"])] = r
        z = np.load(WORK / f"pertarget250-{panel}-{held}.npz")
        for k in z.files:
            pt[(panel, held, k)] = z[k]
    bridge = {}
    for f in sorted(WORK.glob("bridge250-*.json")):
        d = json.loads(f.read_text())
        for r in d["rows"]:
            bridge[(r["fold"], r["operator"], r["arm"])] = r
    basis = json.loads((WORK / "basis250.json").read_text()) if (WORK / "basis250.json").exists() else None
    return rows, pt, meta, bridge, basis


def price6(rel_macro: float, d_pds: float) -> dict:
    """246 §4: cosine legs from `cd4f034`'s rate card at k_cos = 1 on the MACRO cosine; `pds`
    from PLOLO at transfer 0.7x / 1.0x / 2.0x (raw / 0.452 / 6 x 1.011 x T)."""
    d_cos_abs = BOARD_COS * rel_macro
    legs = {m: RATE[m] * d_cos_abs / 0.01 for m in ("mse", "reach", "nmae", "fid", "jac")}
    cos_total = sum(legs.values())
    pds1 = d_pds / SPAN_PDS / 6.0 * 1.011
    out = dict(legs=legs, cos_legs_total=cos_total,
               pds_0p7=0.7 * pds1, pds_1p0=pds1, pds_2p0=2.0 * pds1,
               total_0p7=cos_total + 0.7 * pds1, total_1p0=cos_total + pds1,
               total_2p0=cos_total + 2.0 * pds1, pds_sigma=abs(d_pds) / SIGMA_PDS)
    return out


def summarise(rows, panel, arm, folds, base=BASE):
    sel = [(f, rows[(panel, f, arm)], rows[(panel, f, base)]) for f in folds
           if (panel, f, arm) in rows and (panel, f, base) in rows]
    if not sel:
        return None
    agg = [r["cos_agg"] / b["cos_agg"] - 1 for _, r, b in sel]
    mac = [r["a23_cos"] / b["a23_cos"] - 1 for _, r, b in sel]
    dp = [r["pds"] - b["pds"] for _, r, b in sel]
    dq = [r["a15_sgn_w"] - b["a15_sgn_w"] for _, r, b in sel]
    aud = sel[0][1].get("audit", {})
    s = dict(n=len(sel), rel_agg=float(np.mean(agg)), rel_macro=float(np.mean(mac)),
             d_pds=float(np.mean(dp)), d_qw=float(np.mean(dq)),
             pds_pos=int(sum(x > 0 for x in dp)), cos_pos=int(sum(x > 0 for x in mac)),
             energy=float(np.mean([r["audit"].get("energy_removed_frac", np.nan) for _, r, _ in sel])),
             frac=float(np.mean([r["audit"]["frac_cd4"] for _, r, _ in sel])),
             per_fold={f: dict(rel_agg=a, rel_macro=m, d_pds=p, d_qw=q, frac=r["audit"]["frac_cd4"])
                       for (f, r, _), a, m, p, q in zip(sel, agg, mac, dp, dq)},
             audit_first=aud)
    s["price"] = price6(s["rel_macro"], s["d_pds"])
    return s


def concentration(pt, panel, arm, folds, base=BASE):
    res = {}
    for k in (1, 3, 5):
        rp, rc = [], []
        for f in folds:
            if (panel, f, f"{arm}|pds") not in pt:
                continue
            dp = pt[(panel, f, f"{arm}|pds")] - pt[(panel, f, f"{base}|pds")]
            dc = pt[(panel, f, f"{arm}|cos")] - pt[(panel, f, f"{base}|cos")]
            op = np.argsort(-dp)[:k]; keep_p = np.ones(dp.size, bool); keep_p[op] = False
            oc = np.argsort(-dc)[:k]; keep_c = np.ones(dc.size, bool); keep_c[oc] = False
            rp.append((dp[keep_p].mean(), dp.mean()))
            rc.append((dc[keep_c].mean(), dc.mean()))
        if not rp:
            continue
        fp = float(np.mean([t for _, t in rp])); lp = float(np.mean([l for l, _ in rp]))
        fc = float(np.mean([t for _, t in rc])); lc = float(np.mean([l for l, _ in rc]))
        res[k] = dict(pds_full=fp, pds_loo=lp, pds_retained=(lp / fp) if fp else float("nan"),
                      cos_full=fc, cos_loo=lc, cos_retained=(lc / fc) if fc else float("nan"))
    return res


def noncarry_cos(pt, meta, panel, arm, folds):
    out = {}
    for f in folds:
        if (panel, f, f"{arm}|cos") not in pt:
            continue
        carried = np.asarray(meta[(panel, f)]["carried_cd4"], bool)
        c = pt[(panel, f, f"{arm}|cos")]; b = pt[(panel, f, f"{BASE}|cos")]
        nc = ~carried
        out[f] = dict(n_noncarry=int(nc.sum()),
                      rel_out=float(c[nc].mean() / b[nc].mean() - 1) if nc.any() else 0.0,
                      rel_in=float(c[carried].mean() / b[carried].mean() - 1))
    return out


def fmt(s, arm):
    p = s["price"]
    return (f"{arm:22s} {s['n']:2d} {100*s['rel_agg']:+8.3f} {100*s['rel_macro']:+8.3f} "
            f"{s['d_pds']:+9.5f} {s['d_qw']:+9.5f} {s['pds_pos']:>3d}/{s['n']} {s['energy']:6.3f} "
            f"{s['frac']:6.3f} {p['total_0p7']:+8.5f} {p['total_1p0']:+8.5f} {p['total_2p0']:+8.5f}")


HDR = (f"{'arm':22s} {'nf':>2s} {'agg%':>8s} {'macro%':>8s} {'d pds':>9s} {'d q_w':>9s} "
       f"{'+pds':>5s} {'E rm':>6s} {'frac':>6s} {'px0.7':>8s} {'px1.0':>8s} {'px2.0':>8s}")


def main() -> None:
    rows, pt, meta, bridge, basis = load()
    out: dict = {}
    panels = [p for p in PANELS if any(k[0] == p for k in rows)]
    print("panels present:", panels, " folds:", sorted({k[1] for k in rows}))

    # ================================================================= KNOWN ANSWERS
    print("\n" + "=" * 120 + "\nKNOWN-ANSWER ROWS\n" + "=" * 120)
    ka = {}
    for (panel, held, arm), r in sorted(rows.items()):
        if "ka" in r:
            ka[f"{panel}|{held}|{arm}"] = r["ka"]
            lab = ("NO-REF (233 dropped kolf2; bridge vs 246 covers it)" if r['ka'].get('missing')
                   else ('PASS' if r['ka'].get('ok') else 'FAIL'))
            print(f"  {panel:5} {held:14} {arm:8} vs {r['ka']['ref']:5} in {r['ka']['source']}: "
                  f"max|d| {r['ka'].get('max_abs_diff', float('nan')):.1e}  {lab}")
    for (held, op, arm), r in sorted(bridge.items()):
        if "ka" in r:
            ka[f"bridge|{held}|{op}|{arm}"] = r["ka"]
            print(f"  bridge {held:14} {op:10} {arm:26} vs 246: max|d| {r['ka']['max_abs_diff']:.1e}  "
                  f"{'PASS' if r['ka']['ok'] else 'FAIL'}")
    if basis:
        ka["basis"] = basis["ka"]
        print(f"  basis KA-P1 (production residual bitwise): {basis['ka']['KA-P1']}")
        print(f"  basis KA-P2 (shipped deflation audit): {basis['ka']['KA-P2']}")
    out["ka"] = ka
    n_fail = sum(1 for v in ka.values() if isinstance(v, dict) and v.get("ok") is False)
    print(f"  known-answer failures: {n_fail}")

    # ================================================================= TWO INSTRUMENTS on A/B/C
    print("\n" + "=" * 120 + "\nTHE TWO INSTRUMENTS ON A/B/C -- 246's B2_cd4_svd5_resolved through three constructions\n"
          "  op062+246fold = 246 verbatim (062 held operators, plolo float32 consensus)\n"
          "  refit233+246fold = the same objects through 233's essential-pruned re-fit operator\n"
          "  250 sweep abc = 233's instrument (re-fit operator, fold233 float64 consensus)\n" + "=" * 120)
    inst = {}
    for held in CANCER + ("kolf2",):
        for label, get in (("op062+246fold", lambda a: bridge.get((held, "op062" if held != "kolf2" else "production", a))),
                           ("refit233+246fold", lambda a: bridge.get((held, "refit233", a))),
                           ("250 sweep abc", lambda a: rows.get(("abc", held, {"STACK": "BASE", "B2_cd4_svd5_resolved": "svd5", "B2_cd4_svd5": "svd5_fixed"}[a])))):
            b, r = get("STACK"), get("B2_cd4_svd5_resolved")
            if not (b and r):
                continue
            inst[f"{held}|{label}"] = dict(
                base_cos_agg=b["cos_agg"], base_macro=b["a23_cos"], base_pds=b["pds"],
                rel_agg=r["cos_agg"] / b["cos_agg"] - 1, rel_macro=r["a23_cos"] / b["a23_cos"] - 1,
                d_pds=r["pds"] - b["pds"], d_qw=r["a15_sgn_w"] - b["a15_sgn_w"],
                frac_stack=b["audit"]["frac_cd4"], frac_arm=r["audit"]["frac_cd4"])
    for held in CANCER + ("kolf2",):
        for label in ("op062+246fold", "refit233+246fold", "250 sweep abc"):
            v = inst.get(f"{held}|{label}")
            if v:
                print(f"  {held:14} {label:18} base cos_agg {v['base_cos_agg']:.6f} pds {v['base_pds']:.6f} | "
                      f"svd5 resolved: agg {100*v['rel_agg']:+.3f}% macro {100*v['rel_macro']:+.3f}% "
                      f"d pds {v['d_pds']:+.6f} d q_w {v['d_qw']:+.5f} frac {v['frac_stack']:.4f}->{v['frac_arm']:.4f}")
    for label in ("op062+246fold", "refit233+246fold", "250 sweep abc"):
        vs = [inst[f"{h}|{label}"] for h in CANCER if f"{h}|{label}" in inst]
        if len(vs) == 3:
            m = dict(rel_agg=float(np.mean([v["rel_agg"] for v in vs])),
                     rel_macro=float(np.mean([v["rel_macro"] for v in vs])),
                     d_pds=float(np.mean([v["d_pds"] for v in vs])),
                     d_qw=float(np.mean([v["d_qw"] for v in vs])))
            inst[f"MEAN3|{label}"] = m
            print(f"  {'MEAN 3 cancer':14} {label:18} agg {100*m['rel_agg']:+.4f}% macro {100*m['rel_macro']:+.4f}% "
                  f"d pds {m['d_pds']:+.6f} d q_w {m['d_qw']:+.5f}")
    out["instruments"] = inst

    # ================================================================= 233 RETENTION (KA-2)
    print("\n" + "=" * 120 + "\n233 RETENTION ROWS re-derived from 250's BASE233/null233 (credit above the null, panel/abc)\n" + "=" * 120)
    ret233 = json.loads((R233 / "retain233.json").read_text())
    an233 = json.loads((R233 / "analyse233.json").read_text())
    ret = {}
    for panel in [p for p in panels if p != "abc"]:
        for key in ("cos_agg", "cos", "pds"):
            vals = []
            for f in CANCER:
                a, b = rows.get(("abc", f, "BASE233")), rows.get((panel, f, "BASE233"))
                na, nb = rows.get(("abc", f, "null233")), rows.get((panel, f, "null233"))
                if not (a and b and na and nb):
                    vals.append(float("nan")); continue
                ca, cb = a[key] - na[key], b[key] - nb[key]
                vals.append(cb / ca)
            ref = ret233.get(f"retention|{panel}", {}).get(key)
            ret[f"{panel}|{key}"] = dict(mine=vals, mine_mean=float(np.nanmean(vals)), ref233=ref)
            print(f"  {panel:5} {key:8} mine " + " ".join(f"{100*v:7.1f}%" for v in vals)
                  + f"  mean {100*np.nanmean(vals):6.1f}%   233: "
                  + (" ".join(f"{100*v:7.1f}%" for v in ref) if ref else "n/a")
                  + (f"   max|d| {max(abs(x - y) for x, y in zip(vals, ref[:3])):.1e}" if ref else ""))
        for f in ("orion_hct116", "orion_hek293t"):
            k = f"{f}|{panel}"
            if k in an233["thinning"] and (panel, f, "BASE233") in rows:
                a, b = rows[("abc", f, "BASE233")], rows[(panel, f, "BASE233")]
                na, nb = rows[("abc", f, "null233")], rows[(panel, f, "null233")]
                rc = 100 * (b["cos_agg"] - nb["cos_agg"]) / (a["cos_agg"] - na["cos_agg"])
                rp = 100 * (b["pds"] - nb["pds"]) / (a["pds"] - na["pds"])
                ref = an233["thinning"][k]
                ret[f"analyse233|{k}"] = dict(ret_cos=rc, ret_pds=rp, ref=ref)
                print(f"  analyse233 {k:22} ret_cos {rc:.4f} (233 {ref['ret_cos']:.4f}, d {rc-ref['ret_cos']:.1e})  "
                      f"ret_pds {rp:.4f} (233 {ref['ret_pds']:.4f}, d {rp-ref['ret_pds']:.1e})")
    out["retention233"] = ret

    # ================================================================= Phase 1: k LADDER
    ladder = {}
    for panel in panels:
        folds = CANCER
        print("\n" + "=" * 120 + f"\nPHASE 1 -- k LADDER on `{panel}` (3 cancer folds vs BASE = the panel's stack, "
              f"sizing {meta[(panel, folds[0])]['sizing']}; resolved)\n" + "=" * 120)
        print(HDR)
        arms = ["BASE233", "BASE_carrying", "svd5_carrying", "svd5_allrows"] + [f"svd{k}" for k in KS] + ["svd5_fixed"]
        summ = {}
        for arm in arms:
            s = summarise(rows, panel, arm, folds)
            if s:
                summ[arm] = s
                print(fmt(s, arm))
        # per fold for the ladder
        for k in KS:
            s = summ.get(f"svd{k}")
            if s:
                print(f"    svd{k:<3d} per fold: " + "  ".join(
                    f"{f[:8]} agg {100*v['rel_agg']:+.2f}% pds {v['d_pds']:+.5f}" for f, v in s["per_fold"].items()))
        dps = {k: summ[f"svd{k}"]["d_pds"] for k in KS if f"svd{k}" in summ}
        if dps:
            kmax = max(dps, key=dps.get)
            plateau = [k for k in KS if dps[k] >= PLATEAU * dps[kmax]]
            knee = min(plateau) if plateau else None
            k5_in = 5 in plateau
            prices = {k: summ[f"svd{k}"]["price"]["total_1p0"] for k in dps}
            kbest = max(prices, key=prices.get)
            ladder[panel] = dict(summary=summ, d_pds=dps, k_max=kmax, plateau=plateau, knee=knee,
                                 k5_in_plateau=k5_in, price_1p0=prices, k_best_price=kbest,
                                 pds_pos_at5=summ["svd5"]["pds_pos"] if "svd5" in summ else None)
            print(f"  -> d pds max at k={kmax} ({dps[kmax]:+.5f}); plateau (>= {PLATEAU:.0%} of max): {plateau}; "
                  f"knee k={knee}; k=5 inside plateau: {k5_in}; best price(1.0x) k={kbest}")
        # kolf2 apart
        if (panel, "kolf2", BASE) in rows:
            print("  kolf2 apart: " + "; ".join(
                f"svd{k} agg {100*s['rel_agg']:+.2f}% pds {s['d_pds']:+.5f}"
                for k in KS if (s := summarise(rows, panel, f"svd{k}", ("kolf2",)))))
            ladder[panel]["kolf2"] = {f"svd{k}": summarise(rows, panel, f"svd{k}", ("kolf2",)) for k in KS}
    out["ladder"] = ladder

    # ---- plateau ratio, price-optimal k (refit and pinned), leave-one-fold-out k selection
    print("\n" + "=" * 120 + "\nk SELECTION per panel: d pds(5)/max, price-optimal k for the re-fit and the pinned A/B/C basis, LOFO\n" + "=" * 120)
    ksel = {}
    for panel in panels:
        if panel not in ladder:
            continue
        dps = ladder[panel]["d_pds"]
        ratio5 = dps[5] / max(dps.values())
        pin_arms = {k: f"pinABC{k}" for k in (3, 5, 10)}
        pin_s = {k: summarise(rows, panel, a, CANCER) for k, a in pin_arms.items()}
        pin_price = {k: v["price"]["total_1p0"] for k, v in pin_s.items() if v}
        pin_dp = {k: v["d_pds"] for k, v in pin_s.items() if v}
        refit_price = ladder[panel]["price_1p0"]
        # LOFO: pick k on two folds by the 1.0x price, read on the third (refit ladder)
        lofo = []
        for held in CANCER:
            others = [f for f in CANCER if f != held]
            best, bestv = None, -1e9
            for k in KS:
                s2 = summarise(rows, panel, f"svd{k}", others)
                if s2 and s2["price"]["total_1p0"] > bestv:
                    best, bestv = k, s2["price"]["total_1p0"]
            s3 = summarise(rows, panel, f"svd{best}", (held,))
            lofo.append((held, best, s3["price"]["total_1p0"], s3["d_pds"]))
        lofo_mean = float(np.mean([v for _, _, v, _ in lofo]))
        ksel[panel] = dict(ratio5=ratio5, refit_price_1p0=refit_price, pinned_price_1p0=pin_price,
                           pinned_d_pds=pin_dp, k_best_refit=max(refit_price, key=refit_price.get),
                           k_best_pinned=(max(pin_price, key=pin_price.get) if pin_price else None),
                           lofo=lofo, lofo_mean=lofo_mean)
        print(f"  {panel:5} d pds(5)/max {ratio5:.3f} | refit price(1.0x) by k: "
              + " ".join(f"k{k} {v:+.4f}" for k, v in refit_price.items())
              + f" -> best k={ksel[panel]['k_best_refit']} | pinned ABC price(1.0x): "
              + " ".join(f"k{k} {v:+.4f}" for k, v in pin_price.items())
              + f" -> best k={ksel[panel]['k_best_pinned']} | LOFO(refit) {lofo_mean:+.4f} "
              + ", ".join(f"{h[:8]}->k{k} {v:+.4f}" for h, k, v, _ in lofo))
    out["kselect"] = ksel

    # ---- cross-panel pinned bases (cross250), if present
    cross = {}
    for f in sorted(WORK.glob("cross250-*-*.json")):
        d = json.loads(f.read_text())
        for r in d["rows"]:
            cross[(r["panel"], r["fold"], r["arm"])] = r
    if cross:
        print("\n" + "=" * 120 + "\nPHASE 3b -- CROSS-PANEL PINNED BASES: rows = the panel deflated, columns = whose production V_k is pinned "
              "(d pds / rel agg %, 3 cancer folds, resolved)\n" + "=" * 120)
        xt = {}
        for k in (5, 10):
            print(f"  k={k:<3d}" + "".join(f"{'pin ' + src:>22}" for src in PANELS) + f"{'own refit':>22}")
            for panel in panels:
                cells = []
                for src in PANELS:
                    sel = [(cross[(panel, f, f"pin{src}{k}")], cross[(panel, f, "BASE")]) for f in CANCER
                           if (panel, f, f"pin{src}{k}") in cross]
                    if len(sel) != 3:
                        cells.append(("", None)); continue
                    dp = float(np.mean([r["pds"] - b["pds"] for r, b in sel]))
                    ra = float(np.mean([r["cos_agg"] / b["cos_agg"] - 1 for r, b in sel]))
                    rm = float(np.mean([r["a23_cos"] / b["a23_cos"] - 1 for r, b in sel]))
                    pos = sum(r["pds"] > b["pds"] for r, b in sel)
                    xt[f"{panel}|pin_{src}|{k}"] = dict(d_pds=dp, rel_agg=ra, rel_macro=rm, pds_pos=pos,
                                                       price_1p0=price6(rm, dp)["total_1p0"],
                                                       captured=float(np.mean([r["audit"]["captured_frac_of_centred"] for r, _ in sel])))
                    cells.append((f"{dp:+.4f}/{100*ra:+.1f}%({pos})", dp))
                own = ladder[panel]["summary"].get(f"svd{k}")
                ownc = f"{own['d_pds']:+.4f}/{100*own['rel_agg']:+.1f}%({own['pds_pos']})" if own else ""
                print(f"  {panel:6}" + "".join(f"{c[0]:>22}" for c in cells) + f"{ownc:>22}")
        out["cross_pinned"] = xt

    # ================================================================= Phase 2: CONTROLS
    ctrl = {}
    for panel in panels:
        folds = CANCER
        print("\n" + "=" * 120 + f"\nPHASE 2 -- CONTROLS at k = 5 on `{panel}` (vs BASE; resolved unless _fixed)\n" + "=" * 120)
        print(HDR)
        arms = ["svd5", "svd5_fixed", "randE5", "randE5b", "randE5_fixed", "randE5b_fixed", "svd5_SHUFDEF",
                "svd3_splithalf", "svd3_splithalfb", "svd5_splithalf", "svd5_splithalfb",
                "svd10_splithalf", "svd10_splithalfb"]
        cs = {}
        for arm in arms:
            s = summarise(rows, panel, arm, folds)
            if s:
                cs[arm] = s
                extra = ""
                if "rand_dim" in s["audit_first"]:
                    extra = f"   dims {s['audit_first']['rand_dim']}"
                print(fmt(s, arm) + extra)
        conc = {arm: concentration(pt, panel, arm, folds) for arm in
                ("svd3", "svd5", "svd7", "svd10", "svd5_fixed", "pinABC5", "svd5+partial@0.34", "partial@0.34")
                if (panel, folds[0], f"{arm}|pds") in pt}
        print("  leave-top-k-out retention (pds and cosine on their own orderings):")
        for arm, c in conc.items():
            if 3 in c:
                c3 = c[3]
                print(f"    {arm:20s} d pds full {c3['pds_full']:+.6f} LT3O {c3['pds_loo']:+.6f} "
                      f"({100*c3['pds_retained']:6.1f}% ret; LT1O {100*c[1]['pds_retained']:6.1f}%) | "
                      f"d cos full {c3['cos_full']:+.6f} LT3O {c3['cos_loo']:+.6f} ({100*c3['cos_retained']:6.1f}% ret)")
        ncx = {f: meta[(panel, f)]["noncarrying"] for f in folds if (panel, f) in meta}
        ncc = noncarry_cos(pt, meta, panel, "svd5", folds)
        print("  non-carrying rows (rows CD4 does not carry): exact test vs BASE x preserve ratio, and per-target cosine change:")
        for f in folds:
            if f in ncx:
                e = ncx[f]["svd5"]
                print(f"    {f:14} n {e.get('n_noncarry', 0):3d} max|d| {e.get('max_abs_diff', 0.0):.1e} "
                      f"(rel {e.get('max_rel', 0.0):.1e}, preserve ratio {e.get('preserve_ratio', 1.0):.6f}) | "
                      f"mean cos rel_out {100*ncc[f]['rel_out']:+.4f}% rel_in {100*ncc[f]['rel_in']:+.3f}%")
        split = [cs[a]["d_pds"] for a in ("svd5_splithalf", "svd5_splithalfb") if a in cs]
        ctrl[panel] = dict(summary=cs, concentration=conc, noncarry_exact=ncx, noncarry_cos=ncc,
                           splithalf5_mean=float(np.mean(split)) if split else None,
                           splithalf5_retained=(float(np.mean(split)) / cs["svd5"]["d_pds"]) if split and "svd5" in cs and cs["svd5"]["d_pds"] else None,
                           rand5_max=max(cs[a]["d_pds"] for a in ("randE5", "randE5b") if a in cs) if "randE5" in cs else None,
                           shufdef5=cs["svd5_SHUFDEF"]["d_pds"] if "svd5_SHUFDEF" in cs else None)
    out["controls"] = ctrl

    # ================================================================= Phase 3: PINNED BASIS
    pin = {}
    for panel in panels:
        folds = CANCER
        print("\n" + "=" * 120 + f"\nPHASE 3 -- the A/B/C PRODUCTION basis (291 rows, U_5 sha f06e0262) pinned vs the panel's own re-fit, `{panel}`\n" + "=" * 120)
        print(HDR)
        ps = {}
        for arm in ("svd3", "pinABC3", "svd5", "pinABC5", "pinABC5_fixed", "svd10", "pinABC10"):
            s = summarise(rows, panel, arm, folds)
            if s:
                ps[arm] = s
                cap = s["audit_first"].get("captured_frac_of_centred")
                print(fmt(s, arm) + (f"   captured {cap:.4f}" if cap is not None else ""))
        for k in (3, 5, 10):
            if f"svd{k}" in ps and f"pinABC{k}" in ps:
                print(f"    k={k}: pinned - refit: d pds {ps[f'pinABC{k}']['d_pds'] - ps[f'svd{k}']['d_pds']:+.5f}, "
                      f"agg {100*(ps[f'pinABC{k}']['rel_agg'] - ps[f'svd{k}']['rel_agg']):+.3f} pp; per fold pinned d pds: "
                      + " ".join(f"{f[:8]} {v['d_pds']:+.5f}" for f, v in ps[f"pinABC{k}"]["per_fold"].items()))
        pin[panel] = ps
    out["pin"] = pin
    if basis:
        print("\n  principal angles (deg) between production V_5 bases and cross-captured energy (basis250.json):")
        for key, ang in basis["angles"].items():
            a, b, k = key.split("|")
            if k == "5":
                print(f"    {a:5}|{b:5} k=5 angles " + " ".join(f"{x:5.1f}" for x in ang))
        for p in PANELS:
            c_own = basis["cross"].get(f"{p}|by_{p}|5", {}).get("captured")
            c_abc = basis["cross"].get(f"{p}|by_abc|5", {}).get("captured")
            if c_own is not None:
                print(f"    {p:5}: top-5 energy own {c_own:.4f}  captured by abc's V_5 {c_abc:.4f} "
                      f"(ratio {c_abc/c_own:.3f}; random-5 floor {basis['random_floor'][p]['random5_captured']:.4f}); "
                      f"n_carried {basis['audits'][f'{p}|5']['n_carried']}; spectrum "
                      f"{[round(x, 4) for x in basis['audits'][f'{p}|5']['spectrum_frac'][:6]]}")
        print("  fold-level own V_5 vs production bases (angles, deg):")
        for key, v in basis["fold_angles"].items():
            print(f"    {key:24} vs own production " + " ".join(f"{x:5.1f}" for x in v["vs_own_production"])
                  + " | vs abc production " + " ".join(f"{x:5.1f}" for x in v["vs_abc_production"])
                  + f" | top5 {v['audit']['top_k_energy_frac_of_centred']:.4f}")

    # ================================================================= Phase 4: +partial
    comp = {}
    for panel in panels:
        folds = CANCER
        print("\n" + "=" * 120 + f"\nPHASE 4 -- `+partial` composition on `{panel}` (partial rows carrying-row sized at 0.34)\n" + "=" * 120)
        print(HDR)
        cs = {}
        for arm in ("partial@0.34", "partial_SHUF@0.34", "svd5", "svd5+partial@0.34"):
            s = summarise(rows, panel, arm, folds)
            if s:
                cs[arm] = s
                print(fmt(s, arm) + f"   n_partial {s['audit_first'].get('n_partial', '-')}")
        if "svd5+partial@0.34" in cs and "svd5" in cs:
            inc = summarise(rows, panel, "svd5+partial@0.34", folds, base="svd5")
            cs["increment_over_svd5"] = inc
            print(f"    increment of +partial over svd5: agg {100*inc['rel_agg']:+.3f}% macro {100*inc['rel_macro']:+.3f}% "
                  f"d pds {inc['d_pds']:+.5f} ({inc['pds_pos']}/{inc['n']} folds)")
        comp[panel] = cs
    out["composition"] = comp

    # ================================================================= VERDICTS
    print("\n" + "=" * 120 + "\nVERDICT per panel\n" + "=" * 120)
    verdict = {}
    for panel in panels:
        if panel not in ladder:
            continue
        L5 = ladder[panel]["summary"].get("svd5")
        c = ctrl[panel]
        conc5 = c["concentration"].get("svd5", {}).get(3, {}).get("pds_retained")
        mech = dict(
            pds_positive_3of3=(L5["pds_pos"] == 3),
            pds_gain_ge_3sigma=(L5["price"]["pds_sigma"] >= 3.0),
            random_controls_below_half_the_arm=(c["rand5_max"] is not None and c["rand5_max"] < 0.5 * L5["d_pds"]),
            shuffled_projector_below_0p001=(c["shufdef5"] is not None and c["shufdef5"] < 0.001),
            splithalf_retains_ge_50pct=(c["splithalf5_retained"] is not None and c["splithalf5_retained"] >= 0.5),
            lt3o_retains_ge_25pct=(conc5 is not None and conc5 >= 0.25))
        pays = dict(price_1p0_above_shape_gate_k5=(L5["price"]["total_1p0"] > GATE_SHAPE),
                    price_1p0_above_shape_gate_bestk=(max(ladder[panel]["price_1p0"].values()) > GATE_SHAPE),
                    pinned_k5_price_1p0=ksel[panel]["pinned_price_1p0"].get(5),
                    pinned_k5_above_gate=(ksel[panel]["pinned_price_1p0"].get(5, -1) > GATE_SHAPE))
        tests = dict(**mech, k5_in_plateau=ladder[panel]["k5_in_plateau"], **{k: v for k, v in pays.items() if isinstance(v, bool)})
        transfers = all(mech.values())
        verdict[panel] = dict(tests=tests, mechanism=mech, pays=pays, transfers=transfers, d_pds5=L5["d_pds"], rel_agg5=L5["rel_agg"],
                              rel_macro5=L5["rel_macro"], knee=ladder[panel]["knee"],
                              plateau=ladder[panel]["plateau"], k_best_price=ladder[panel]["k_best_price"],
                              price_1p0=L5["price"]["total_1p0"], price_0p7=L5["price"]["total_0p7"],
                              price_2p0=L5["price"]["total_2p0"])
        print(f"  {panel:5} mechanism {'TRANSFERS' if transfers else 'DOES NOT'}; pays at k=5 (1.0x) {pays['price_1p0_above_shape_gate_k5']}, "
              f"at best k {pays['price_1p0_above_shape_gate_bestk']}, pinned ABC k=5 {pays['pinned_k5_price_1p0']:+.4f} ({pays['pinned_k5_above_gate']}): "
              f"k=5 d pds {L5['d_pds']:+.5f} "
              f"({L5['pds_pos']}/3, {L5['price']['pds_sigma']:.1f} sigma), agg {100*L5['rel_agg']:+.2f}%, macro {100*L5['rel_macro']:+.2f}%, "
              f"price 0.7/1.0/2.0x {L5['price']['total_0p7']:+.4f}/{L5['price']['total_1p0']:+.4f}/{L5['price']['total_2p0']:+.4f}; "
              f"knee k={ladder[panel]['knee']}, plateau {ladder[panel]['plateau']}, k=5 in plateau {ladder[panel]['k5_in_plateau']}; "
              f"failed tests: {[t for t, ok in tests.items() if not ok]}")
    out["verdict"] = verdict
    out["gates"] = dict(shape=GATE_SHAPE, amplitude=GATE_AMPL, plateau=PLATEAU)
    (WORK / "analyse250.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"\n-> {WORK / 'analyse250.json'}")


if __name__ == "__main__":
    main()
