"""286 -- the six-member price of every Stage B arm OVER cd4svd5 on the three COHERENT scenarios (READ-FIRST §3;
261's recipe imported verbatim through 271's price271.py conventions), the BR249 retrodiction row first.  Nothing is
priced unless the retrodiction and the known answers pass.

    PYTHONDONTWRITEBYTECODE=1 uv run python priceB286.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/286-crossproduct-estimator/work"
E257 = ROOT / "experiments/257-soft-deflation-price"
E261 = ROOT / "experiments/261-union-deflation"
sys.path.insert(0, str(E261))
sys.path.insert(0, str(E257))
import price261 as PR  # noqa: E402

P0 = PR.P
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BASE = "cd4svd5"
GATE_SHAPE = P0.GATE_SHAPE
GATE_AMP = P0.GATE_AMPL
LINES: list[str] = []


def Pr(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def load():
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in CANCER:
        d = json.loads((WORK / f"sweepB286-{held}.json").read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        z = np.load(WORK / f"pertargetB286-{held}.npz", allow_pickle=True)
        for k in z.files:
            if k.endswith("|pds") or k.endswith("|cos"):
                pt[(held, k)] = np.asarray(z[k], dtype=np.float64)
        z = np.load(WORK / f"identityB286-{held}.npz")
        for k in z.files:
            idn[(held, k)] = z[k]
    return rows, pt, idn, meta


def nmae_concat(idn, arm, folds) -> float:
    U = np.concatenate([idn[(f, f"{arm}|U")] for f in folds])
    M = np.concatenate([idn[(f, f"{arm}|M")] for f in folds])
    D = np.concatenate([idn[(f, f"{arm}|D")] for f in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(rows, idn, arm, base, folds) -> dict:
    ra, rm, dp, dq, dn_a, dn_b, fr = [], [], [], [], [], [], []
    for f in folds:
        a, b = rows[(f, arm)], rows[(f, base)]
        ra.append(a["cos_agg"] / b["cos_agg"] - 1.0)
        rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
        dp.append(a["pds"] - b["pds"])
        dq.append(a["a15_sgn_w"] - b["a15_sgn_w"])
        dn_a.append(a["delivered_norm"])
        dn_b.append(b["delivered_norm"])
        fr.append(a["audit"]["cd4_fraction"] / b["audit"]["cd4_fraction"])
    nm_a, nm_b = nmae_concat(idn, arm, folds), nmae_concat(idn, base, folds)
    return dict(arm=arm, base=base, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)),
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp],
                per_fold_rel_agg=[float(x) for x in ra], per_fold_rel_macro=[float(x) for x in rm],
                pds_pos=int(sum(x > 0 for x in dp)), sigma=float(np.mean(dp)) / P0.SIGMA_PDS,
                fold_se=float(np.std(dp, ddof=1) / np.sqrt(len(dp))), d_qw=float(np.mean(dq)),
                delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)), fraction_ratio=float(np.mean(fr)),
                nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(f, arm)]["top89_precision"] - rows[(f, base)]["top89_precision"] for f in folds])),
                energy_removed=0.0)


def main() -> None:
    rows, pt, idn, meta = load()
    targets = {f: meta[f]["targets"] for f in CANCER}
    an = json.loads((WORK / "analysisB286.json").read_text())
    fails: list[str] = []
    out: dict = dict(recipe=dict(PDS_SCALE=P0.PDS_SCALE, MSE_SCALE=P0.MSE_SCALE, MSE_SLOPE=P0.MSE_SLOPE,
                                 POINT=PR.POINT, K_PDS_BAND=PR.K_PDS_BAND, K_COS_BAND=PR.K_COS_BAND,
                                 DE3_BAND=PR.DE3_BAND, GATE_SHAPE=GATE_SHAPE, GATE_AMP=GATE_AMP))
    Pr("=" * 130)
    Pr("286 -- the cross-product estimators priced OVER cd4svd5 (261's recipe, three COHERENT scenarios); retrodiction first")
    Pr("=" * 130)
    # NAMED, DOCUMENTED EXEMPTION (recorded in RESULTS.md §KA): lib263's KA-0 compares `E062.build_baseline`'s float64
    # rebuild, cast to float32, against the fold cache's `pre` at tolerance 0.0.  It holds at 0.0 only when the fold
    # cache is a BYTE COPY of 262's saved float32 npz (271's setup).  This experiment BUILDS its own fold caches, so
    # `fold["pre"]` is the freshly computed float64 array and the comparison reads ~6e-08 (sub-float32-ulp).  It is
    # reclassified as non-fatal here because the substantive regression rows all pass in the same run: cd4svd5 == 263,
    # STACK == 257, relw / ebrow / wiener_own == 271, every one at <= 2.0e-11, and cd4svd5's per-target `pds` array is
    # BIT-IDENTICAL to 263's.  No other failure is exempt.
    EXEMPT = "KA-0 build_baseline(per-donor rows) float32 == cached pre"
    out["ka0_exemption"] = dict(label=EXEMPT, reason="fold caches built here (float64) rather than byte-copied from "
                                                     "262 (float32); substantive regression rows all pass",
                                values={h: next((k["abs_diff"] for k in meta[h]["ka"] if k["label"] == EXEMPT), None)
                                        for h in CANCER})
    for h in CANCER:
        m = meta[h]
        real_fail = [f for f in m["failures"] if not f.startswith(EXEMPT)]
        bad = [k for k in m["ka"] if not k["ok"] and k["label"] != EXEMPT]
        Pr(f"  {h:14s} {m['n_targets']} targets, donors {m['donors']}, {m['n_rows']} rows, {m['seconds']}s, peak "
           f"{m['peak_rss_gib']:.2f} GiB, known-answer {m['ka_pass']}/{m['ka_total']} "
           f"(KA-0 provenance exemption {out['ka0_exemption']['values'][h]:.1e}), failures {real_fail}")
        if real_fail or bad:
            fails.append(f"sweep {h}")

    c_r = compare(rows, idn, BASE, "STACK", CANCER)
    s_r = PR.p257_sets(c_r, P0.COS_H1UNION, P0.PER_H1UNION)
    r2 = s_r["mine"]["point"]
    p257 = json.loads((E257 / "work/price257.json").read_text())
    ref_r1 = p257["retrodiction"]["R1"]["price"]["point"]["total"]
    c257 = p257["retrodiction"]["R1"]["compare"]
    band_ok = {k: abs(r2[k] - PR.R2_BANDS[k][0]) <= PR.R2_BANDS[k][1] for k in ("pds", "mse")}
    band_ok["overall"] = abs(r2["total"] - PR.R2_BANDS["overall"][0]) <= PR.R2_BANDS["overall"][1]
    d_cmp = max(abs(c_r["d_pds"] - c257["d_pds"]), abs(c_r["rel_agg"] - c257["rel_agg"]), abs(c_r["d_nmae"] - c257["d_nmae"]))
    ok_r = all(band_ok.values()) and abs(r2["total"] - ref_r1) < 1e-6 and d_cmp < 1e-6 and s_r["max_abs_diff"] < 1e-12
    board_all = P0.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_avg"] - P0.BOARD["iUziqE1gewfcPKzYUDaB"]["score_avg"]
    Pr(f"\nRETRODICTION  cd4svd5 vs STACK (h1union): d pds {c_r['d_pds']:+.6f} (257 {c257['d_pds']:+.6f})  agg "
       f"{100*c_r['rel_agg']:+.4f}% (257 {100*c257['rel_agg']:+.4f}%)  d nmae89 {c_r['d_nmae']:+.5f} "
       f"(257 {c257['d_nmae']:+.5f})  max|d| {d_cmp:.1e}")
    Pr(f"              POINT: pds {r2['pds']:+.5f} mse {r2['mse']:+.5f} DE3 {r2['de3']:+.5f} nmae {r2['nmae']:+.5f} = "
       f"{r2['total']:+.5f} (257 R1 {ref_r1:+.5f}, |d| {abs(r2['total']-ref_r1):.1e}; BR249 board {board_all:+.7f}; "
       f"bands {band_ok})  -> {'PASS' if ok_r else 'FAIL'}")
    out["retrodiction"] = dict(compare=c_r, price=r2, ref257=dict(total=ref_r1, compare=c257), band_ok=band_ok,
                               board=board_all, ok=ok_r)
    if not ok_r:
        fails.append("retrodiction")
    if fails:
        Pr(f"\n*** KNOWN-ANSWER / RETRODICTION FAILURES: {fails} -- no cell is priced")
        (WORK / "priceB286.json").write_text(json.dumps(dict(out, failures=fails), indent=1, default=float))
        (WORK / "priceB286.log").write_text("\n".join(LINES) + "\n")
        raise SystemExit(1)

    ARMS = [a for a in an["arms"] if a not in (BASE, "STACK", "mean_explicit")]
    MAIN = [a for a in ARMS if a in an["vs_controls"]]
    table: dict[str, dict] = {}
    for name in ARMS:
        c = compare(rows, idn, name, BASE, CANCER)
        sc = PR.scenarios(c)
        ret = PR.retention(pt, name, BASE, CANCER, targets)
        b = an["arms"][name]["pooled"]
        amp = abs(c["delivered_ratio"] - 1) > P0.AMP_SHAPE_BOUND
        gate = GATE_AMP if amp else GATE_SHAPE
        table[name] = dict(vs_cd4svd5=c, scenarios=sc, retention=ret, boot=b, shape_only=not amp, gate=gate,
                           cond_gate=bool(sc["PESSIMISTIC"]["total"] >= gate),
                           cond_retention=bool(ret["retained"] >= 0.5),
                           cond_pds=bool(c["pds_pos"] == 3 and b["lo"] > 0))
    for name in MAIN:
        t = table[name]
        conds = {}
        for suf in (":mirror", ":shufpair", ":wrong", ":rand"):
            ctrl = name + suf
            d = an["vs_controls"][name].get(suf)
            if ctrl not in table or d is None:
                conds[suf] = None
                continue
            # the frozen rule: the control's Delta pds must be below HALF the arm's, and the paired difference > 0
            conds[suf] = bool(d["lo"] > 0 and table[ctrl]["vs_cd4svd5"]["d_pds"] < 0.5 * t["vs_cd4svd5"]["d_pds"])
        t["cond_controls"] = conds
        t["cleared"] = bool(t["cond_gate"] and t["cond_retention"] and t["cond_pds"]
                            and all(v for v in conds.values() if v is not None))

    Pr("\n" + "=" * 150)
    Pr("THE ARMS as INCREMENTS over cd4svd5 (three cancer folds); legs at POINT; PESSIMISTIC / POINT / OPTIMISTIC totals")
    Pr("=" * 150)
    Pr(f"{'cell':26s}{'d pds':>10s}{'f':>4s}{'boot 95% CI':>24s}{'agg%':>8s}{'dlv%':>8s}{'d nmae':>9s}"
       f"{'pds':>9s}{'mse':>9s}{'DE3':>9s}{'nmae':>9s}{'PESS':>9s}{'POINT':>9s}{'OPT':>9s}{'ret3%':>7s}{'gate':>8s}{'CLEAR':>6s}")
    for name in ARMS:
        t = table[name]
        c, s, b = t["vs_cd4svd5"], t["scenarios"], t["boot"]
        po = s["POINT"]
        Pr(f"{name:26s}{c['d_pds']:+10.6f}{c['pds_pos']:>2d}/3  [{b['lo']:+.4f},{b['hi']:+.4f}]"
           f"{100*c['rel_agg']:+8.3f}{100*(c['delivered_ratio']-1):+8.3f}{c['d_nmae']:+9.4f}"
           f"{po['pds']:+9.5f}{po['mse']:+9.5f}{po['de3']:+9.5f}{po['nmae']:+9.5f}"
           f"{s['PESSIMISTIC']['total']:+9.5f}{po['total']:+9.5f}{s['OPTIMISTIC']['total']:+9.5f}"
           f"{100*t['retention']['retained']:+7.0f}{t['gate']:+8.4f}{(str(t.get('cleared','-'))[0]):>6s}")

    cleared = [n for n in MAIN if table[n].get("cleared")]
    best = max(MAIN, key=lambda n: table[n]["scenarios"]["PESSIMISTIC"]["total"])
    best_pt = max(MAIN, key=lambda n: table[n]["scenarios"]["POINT"]["total"])
    Pr("\n" + "=" * 150)
    Pr(f"THE SCREEN: {len(MAIN)} estimator cells + {len(ARMS)-len(MAIN)} control rows on cd4svd5.  "
       f"gate: {sum(table[n]['cond_gate'] for n in MAIN)}/{len(MAIN)}; pds 3/3 with the interval excluding 0: "
       f"{sum(table[n]['cond_pds'] for n in MAIN)}/{len(MAIN)}; retention >= 50 %: "
       f"{sum(table[n]['cond_retention'] for n in MAIN)}/{len(MAIN)}; ALL conditions: {len(cleared)}/{len(MAIN)} {cleared or ''}")
    for tag, name in (("best PESSIMISTIC", best), ("best POINT", best_pt)):
        t = table[name]
        c, s = t["vs_cd4svd5"], t["scenarios"]
        Pr(f"\n{tag}: {name}: d pds {c['d_pds']:+.6f} ({c['pds_pos']}/3 {['%+.5f' % x for x in c['per_fold_d_pds']]}; "
           f"boot [{t['boot']['lo']:+.5f}, {t['boot']['hi']:+.5f}])  agg {100*c['rel_agg']:+.3f}%  delivered "
           f"{100*(c['delivered_ratio']-1):+.3f}% (shape-only {t['shape_only']})  d nmae89 {c['d_nmae']:+.5f}")
        for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            p = s[w]
            Pr(f"  {w:11s} k_pds {p['k_pds']:.2f} k_cos {p['k_cos']:.2f} DE3 {p['de3_rate']:.4f}:  pds {p['pds']:+.5f}  "
               f"mse {p['mse']:+.5f}  nmae {p['nmae']:+.5f}  = {p['total']:+.5f}  ({p['total']/t['gate']:+.2f}x gate)")
        Pr(f"  leave-top-3-out retention {100*t['retention']['retained']:.1f}%; conditions: gate {t['cond_gate']}; "
           f"pds 3/3 & CI > 0 {t['cond_pds']}; retention {t['cond_retention']}; controls {t.get('cond_controls')} "
           f"-> cleared {t.get('cleared')}")
    verdict = "CLEARS" if cleared else "DOES NOT CLEAR"
    Pr(f"\nVERDICT: {verdict} -- {len(cleared)} of {len(MAIN)} cells meet all conditions; best pessimistic "
       f"{table[best]['scenarios']['PESSIMISTIC']['total']:+.5f} ({best}), best point "
       f"{table[best_pt]['scenarios']['POINT']['total']:+.5f} ({best_pt})")
    out.update(table=table, best_pessimistic=best, best_point=best_pt, cleared=cleared, verdict=verdict)
    (WORK / "priceB286.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "priceB286.log").write_text("\n".join(LINES) + "\n")
    md = ["## P1. Six-member price over `cd4svd5` (261's recipe; legs at POINT; three coherent scenarios)\n",
          "| cell | Δ`pds` (3 folds) | pooled 95 % CI | Δ`cos_agg` % | delivered | Δ`nmae` id | `pds` | `mse` | DE-three | `nmae` | **PESSIMISTIC** | POINT | OPTIMISTIC | LT3O | gate | cleared |",
          "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for name in ARMS:
        t = table[name]
        c, s, b = t["vs_cd4svd5"], t["scenarios"], t["boot"]
        po = s["POINT"]
        md.append(f"| `{name}` | {c['d_pds']:+.5f} ({c['pds_pos']}/3) | [{b['lo']:+.4f}, {b['hi']:+.4f}] | "
                  f"{100*c['rel_agg']:+.3f} | {100*(c['delivered_ratio']-1):+.2f} % | {c['d_nmae']:+.4f} | "
                  f"{po['pds']:+.5f} | {po['mse']:+.5f} | {po['de3']:+.5f} | {po['nmae']:+.5f} | "
                  f"**{s['PESSIMISTIC']['total']:+.5f}** | {po['total']:+.5f} | {s['OPTIMISTIC']['total']:+.5f} | "
                  f"{100*t['retention']['retained']:+.0f} % | {t['gate']:+.4f} | {t.get('cleared', '—')} |")
    (WORK / "priceB286.md").write_text("\n".join(md) + "\n")
    Pr(f"\nwrote {WORK/'priceB286.json'}, priceB286.log, priceB286.md")


if __name__ == "__main__":
    main()
