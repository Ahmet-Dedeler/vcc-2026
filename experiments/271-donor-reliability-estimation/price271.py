"""271 -- the six-member price of every estimator OVER cd4svd5 on the three COHERENT scenarios (READ-FIRST §3; 261's
recipe imported verbatim: price261.price / factors / scenarios / retention with 257's constants), the retrodiction row
first (cd4svd5 vs STACK through the recipe at h1union's terms must read BR249: +0.0043 ± 0.0006 overall, and equal
257's R1 to 1e-6).  Nothing is priced unless it passes.

Bands: k_pds {0.85, 1.22}, k_cos {1.36, 3.76}, DE-three {0.0425, 0.24}, identity ×1.5 on a loss / ×0 on a gain
(pessimistic), ×1.0 (point); per leg the band end that lowers the total for THAT cell's sign (PESSIMISTIC) or raises it
(OPTIMISTIC); POINT = the set that retrodicts BR249.  The class of a base RE-ESTIMATION is unmeasured on the board: the
re-dose / re-shape bands are used as the pessimistic reading for the sign of each leg, and said so.

CLEARED (the mandate) = PESSIMISTIC >= +0.0019 AND d pds > 0 on 3/3 folds with the pooled paired-target bootstrap
interval excluding 0 AND leave-top-3-out retention >= 50 % AND beats (d) a13:nm and its shuffled control by more than
their bootstrap width (arm - control paired interval excludes 0 AND d pds(arm) - d pds(control) > half the control's
95 % width).
    PYTHONDONTWRITEBYTECODE=1 uv run python price271.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/271-donor-reliability-estimation"
WORK = HERE / "work"
E257 = ROOT / "experiments/257-soft-deflation-price"
E261 = ROOT / "experiments/261-union-deflation"
sys.path.insert(0, str(E261))
sys.path.insert(0, str(E257))
import price261 as PR  # noqa: E402   (imports price257 as PR.P; read-only at import)

P = PR.P
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
BASE = "cd4svd5"
GATE = P.GATE_SHAPE
MAIN = ["relw", "relw_eb", "relw_glob", "relw_wht", "ebrow", "ebrow_eb", "wiener_own", "wiener_gene", "a13:nm", "a13:lit"]
CONTROLS = ["relw:shuf", "relw:rev", "ebrow:shuf", "wiener_own:shuf", "random:all:s0:nm", "oracle_wht:all:nm"]
CONTROL_OF = {"relw": "relw:shuf", "relw_eb": "relw:shuf", "relw_glob": "relw:shuf", "relw_wht": "relw:shuf", "ebrow": "ebrow:shuf",
              "ebrow_eb": "ebrow:shuf", "wiener_own": "wiener_own:shuf", "wiener_gene": "wiener_own:shuf", "a13:nm": "relw:shuf", "a13:lit": "relw:shuf"}
LINES: list[str] = []


def Pr(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def load():
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in CANCER:
        d = json.loads((WORK / f"sweep271-{held}.json").read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        z = np.load(WORK / f"pertarget271-{held}.npz", allow_pickle=True)
        for k in z.files:
            if k.endswith("|pds") or k.endswith("|cos"):
                pt[(held, k)] = np.asarray(z[k], dtype=np.float64)
        z = np.load(WORK / f"identity271-{held}.npz")
        for k in z.files:
            idn[(held, k)] = z[k]
    return rows, pt, idn, meta


def nmae_concat(idn, arm, folds) -> float:
    U = np.concatenate([idn[(f, f"{arm}|U")] for f in folds])
    M = np.concatenate([idn[(f, f"{arm}|M")] for f in folds])
    D = np.concatenate([idn[(f, f"{arm}|D")] for f in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(rows, idn, arm, base, folds) -> dict:
    """261's `compare` (mean over folds of per-fold ratios / deltas), the fraction ratio read from `cd4_fraction`."""
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
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp], per_fold_rel_agg=[float(x) for x in ra],
                per_fold_rel_macro=[float(x) for x in rm], pds_pos=int(sum(x > 0 for x in dp)),
                sigma=float(np.mean(dp)) / P.SIGMA_PDS, fold_se=float(np.std(dp, ddof=1) / np.sqrt(len(dp))) if len(dp) > 1 else float("nan"),
                d_qw=float(np.mean(dq)), delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)), fraction_ratio=float(np.mean(fr)),
                nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(f, arm)]["top89_precision"] - rows[(f, base)]["top89_precision"] for f in folds])),
                energy_removed=0.0)


def main() -> None:
    rows, pt, idn, meta = load()
    targets = {f: meta[f]["targets"] for f in CANCER}
    an = json.loads((WORK / "analysis271.json").read_text())
    fails: list[str] = []
    out: dict = dict(recipe=dict(PDS_SCALE=P.PDS_SCALE, MSE_SCALE=P.MSE_SCALE, MSE_SLOPE=P.MSE_SLOPE, COS_CD4SVD5=P.COS_CD4SVD5, COS_H1UNION=P.COS_H1UNION,
                                 POINT=PR.POINT, K_PDS_BAND=PR.K_PDS_BAND, K_COS_BAND=PR.K_COS_BAND, DE3_BAND=PR.DE3_BAND, DE3_SPLIT=P.DE3_SPLIT,
                                 NMAE_SPAN=P.NMAE_SPAN, SIGMA_PDS=P.SIGMA_PDS, GATE=GATE, AMP_LEG_THRESHOLD=P.AMP_LEG_THRESHOLD))
    Pr("=" * 120)
    Pr("271 -- the noise-aware donor estimators priced OVER cd4svd5 (261's recipe, three COHERENT scenarios); retrodiction first")
    Pr("=" * 120)
    for h in CANCER:
        m = meta[h]
        Pr(f"  {h:14s} {m['n_targets']} targets, donors {m['donors']}, {m['n_rows']} rows, {m['seconds']}s, peak {m['peak_rss_gib']:.2f} GiB, "
           f"known-answer {m['ka_pass']}/{m['ka_total']}, failures {m['failures']}")
        if m["failures"] or m["ka_pass"] != m["ka_total"]:
            fails.append(f"sweep {h}")

    # ---- the retrodiction row: cd4svd5 vs STACK through the recipe at h1union's terms == BR249 and == 257's R1
    c_r = compare(rows, idn, BASE, "STACK", CANCER)
    s_r = PR.p257_sets(c_r, P.COS_H1UNION, P.PER_H1UNION)
    r2 = s_r["mine"]["point"]
    p257 = json.loads((E257 / "work/price257.json").read_text())
    ref_r1 = p257["retrodiction"]["R1"]["price"]["point"]["total"]
    c257 = p257["retrodiction"]["R1"]["compare"]
    band_ok = {k: abs(r2[k] - PR.R2_BANDS[k][0]) <= PR.R2_BANDS[k][1] for k in ("pds", "mse")}
    band_ok["overall"] = abs(r2["total"] - PR.R2_BANDS["overall"][0]) <= PR.R2_BANDS["overall"][1]
    d_cmp = max(abs(c_r["d_pds"] - c257["d_pds"]), abs(c_r["rel_agg"] - c257["rel_agg"]), abs(c_r["d_nmae"] - c257["d_nmae"]))
    ok_r = all(band_ok.values()) and abs(r2["total"] - ref_r1) < 1e-6 and d_cmp < 1e-6 and s_r["max_abs_diff"] < 1e-12
    board_all = P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_avg"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_avg"]
    Pr(f"\nRETRODICTION  cd4svd5 vs STACK (h1union): d pds {c_r['d_pds']:+.6f} (257 {c257['d_pds']:+.6f})  agg {100*c_r['rel_agg']:+.4f}% (257 {100*c257['rel_agg']:+.4f}%)  "
       f"d nmae89 {c_r['d_nmae']:+.5f} (257 {c257['d_nmae']:+.5f})  max|d| {d_cmp:.1e}")
    Pr(f"              through the recipe at h1union's terms, POINT: pds {r2['pds']:+.5f} mse {r2['mse']:+.5f} DE3 {r2['de3']:+.5f} nmae {r2['nmae']:+.5f} = "
       f"{r2['total']:+.5f} (257 R1 {ref_r1:+.5f}, |d| {abs(r2['total']-ref_r1):.1e}; BR249 board {board_all:+.7f}; bands {band_ok}); "
       f"code-path KA vs price257.price max|d| {s_r['max_abs_diff']:.1e}  -> {'PASS' if ok_r else 'FAIL'}")
    out["retrodiction"] = dict(compare=c_r, price=r2, ref257=dict(total=ref_r1, compare=c257), bands=PR.R2_BANDS, band_ok=band_ok, board=board_all, ok=ok_r)
    if not ok_r:
        fails.append("retrodiction")
    if fails:
        Pr(f"\n*** KNOWN-ANSWER / RETRODICTION FAILURES: {fails} -- no cell is priced")
        (WORK / "price271.json").write_text(json.dumps(dict(out, failures=fails), indent=1, default=float))
        (WORK / "price271.log").write_text("\n".join(LINES) + "\n")
        raise SystemExit(1)

    # ---- the cells
    table: dict[str, dict] = {}
    for name in MAIN + CONTROLS:
        c = compare(rows, idn, name, BASE, CANCER)
        sc = PR.scenarios(c)
        ret = PR.retention(pt, name, BASE, CANCER, targets)
        b = an["arms"][name]["pooled"]
        table[name] = dict(vs_cd4svd5=c, scenarios=sc, retention=ret, boot=b,
                           shape_only=abs(c["delivered_ratio"] - 1) <= P.AMP_SHAPE_BOUND,
                           cond_gate=bool(sc["PESSIMISTIC"]["total"] >= GATE), cond_retention=bool(ret["retained"] >= 0.5),
                           cond_pds=bool(c["pds_pos"] == 3 and b["lo"] > 0))
    for name in MAIN:
        t = table[name]
        vs = an["vs_controls"].get(name, {})
        conds = {}
        for ctrl in (CONTROL_OF[name], "a13:nm"):
            if ctrl == name:
                conds[ctrl] = True
                continue
            bc = table[ctrl]["boot"] if ctrl in table else an["arms"][ctrl]["pooled"]
            half_width = 0.5 * (bc["hi"] - bc["lo"])
            diff = vs.get(ctrl)
            if diff is None:
                # a13:nm vs itself handled above; for a13 vs its control compute from the pooled points
                conds[ctrl] = bool(t["vs_cd4svd5"]["d_pds"] - bc["point"] > half_width)
            else:
                conds[ctrl] = bool(diff["lo"] > 0 and (t["vs_cd4svd5"]["d_pds"] - bc["point"]) > half_width)
        t["cond_controls"] = conds
        t["cleared"] = bool(t["cond_gate"] and t["cond_retention"] and t["cond_pds"] and all(conds.values()))

    Pr("\n" + "=" * 120)
    Pr("THE ESTIMATORS as INCREMENTS over cd4svd5 (three cancer folds); legs at POINT; PESSIMISTIC / POINT / OPTIMISTIC totals")
    Pr("=" * 120)
    Pr(f"{'cell':18s}{'d pds':>10s}{'folds':>6s}{'boot 95% CI':>22s}{'agg%':>8s}{'macro%':>8s}{'dlv%':>7s}{'d nmae':>8s}{'pds':>9s}{'mse':>9s}{'DE3':>9s}{'nmae':>9s}"
       f"{'PESS':>9s}{'POINT':>9s}{'OPT':>9s}{'ret3%':>7s}{'gate':>5s}{'pds':>4s}{'ret':>4s}{'ctl':>4s}{'CLEAR':>6s}")
    for name in MAIN + CONTROLS:
        t = table[name]
        c, s, b = t["vs_cd4svd5"], t["scenarios"], t["boot"]
        pe, po, op = s["PESSIMISTIC"], s["POINT"], s["OPTIMISTIC"]
        ctl = t.get("cond_controls")
        Pr(f"{name:18s}{c['d_pds']:+10.6f}{c['pds_pos']:>4d}/3  [{b['lo']:+.4f}, {b['hi']:+.4f}]{100*c['rel_agg']:+8.3f}{100*c['rel_macro']:+8.3f}"
           f"{100*(c['delivered_ratio']-1):+7.3f}{c['d_nmae']:+8.4f}{po['pds']:+9.5f}{po['mse']:+9.5f}{po['de3']:+9.5f}{po['nmae']:+9.5f}"
           f"{pe['total']:+9.5f}{po['total']:+9.5f}{op['total']:+9.5f}{100*t['retention']['retained']:+7.0f}"
           f"{str(t['cond_gate'])[0]:>5s}{str(t['cond_pds'])[0]:>4s}{str(t['cond_retention'])[0]:>4s}{(str(all(ctl.values()))[0] if ctl else '-'):>4s}{(str(t['cleared'])[0] if 'cleared' in t else '-'):>6s}")
    Pr("  (dlv% = PLOLO frame-A delivered-norm ratio; the amplitude leg is applied beyond ±0.2 %; ret3% = leave-top-3-out retention)")

    best = max(MAIN, key=lambda n: table[n]["scenarios"]["PESSIMISTIC"]["total"])
    best_pt = max(MAIN, key=lambda n: table[n]["scenarios"]["POINT"]["total"])
    cleared = [n for n in MAIN if table[n]["cleared"]]
    Pr("\n" + "=" * 120)
    Pr(f"THE SCREEN: {len(MAIN)} estimator cells + {len(CONTROLS)} control rows on cd4svd5.  PESSIMISTIC >= +{GATE}: {sum(table[n]['cond_gate'] for n in MAIN)}/{len(MAIN)}; "
       f"pds 3/3 with the interval excluding 0: {sum(table[n]['cond_pds'] for n in MAIN)}/{len(MAIN)}; retention >= 50 %: {sum(table[n]['cond_retention'] for n in MAIN)}/{len(MAIN)}; "
       f"ALL conditions: {len(cleared)}/{len(MAIN)} {cleared or ''}")
    for tag, name in (("best PESSIMISTIC", best), ("best POINT", best_pt)):
        t = table[name]
        c, s = t["vs_cd4svd5"], t["scenarios"]
        Pr(f"\n{tag}: {name}: d pds {c['d_pds']:+.6f} ({c['pds_pos']}/3 {['%+.5f' % x for x in c['per_fold_d_pds']]}; boot [{t['boot']['lo']:+.5f}, {t['boot']['hi']:+.5f}])  "
           f"agg {100*c['rel_agg']:+.3f}% macro {100*c['rel_macro']:+.3f}%  delivered {100*(c['delivered_ratio']-1):+.3f}% (shape-only {t['shape_only']})  "
           f"d q_w {c['d_qw']:+.5f}  d nmae89 {c['d_nmae']:+.5f}  d prec89 {c['d_prec89']:+.5f}")
        for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            p = s[w]
            Pr(f"  {w:11s} k_pds {p['k_pds']:.2f} k_cos {p['k_cos']:.2f} DE3 {p['de3_rate']:.4f} id x{p['x_loss']:.1f}/{p['x_gain']:.1f}:  pds {p['pds']:+.5f}  mse {p['mse']:+.5f} "
               f"(d cos_board {p['dcos_board']:+.6f})  reach {p['reach']:+.5f} fid {p['fid']:+.5f} jac {p['jac']:+.5f}  nmae {p['nmae']:+.5f}  = {p['total']:+.5f}  ({p['total']/GATE:+.2f}x gate)")
        Pr(f"  leave-top-3-out retention {100*t['retention']['retained']:.1f}% (per fold {['%.0f%%' % (100*x) for x in t['retention']['per_fold_retained']]}); "
           f"top-3 per fold: " + "; ".join(f"{h}: {', '.join(t['retention']['top_names'][h])}" for h in CANCER))
        Pr(f"  conditions: gate {t['cond_gate']}; pds 3/3 & CI > 0 {t['cond_pds']}; retention {t['cond_retention']}; beats controls {t['cond_controls']} -> cleared {t['cleared']}")
    verdict = "CLEARS" if cleared else "DOES NOT CLEAR"
    Pr(f"\nVERDICT: {verdict} -- {len(cleared)} of {len(MAIN)} cells meet all conditions; best pessimistic {table[best]['scenarios']['PESSIMISTIC']['total']:+.5f} ({best}), "
       f"best point {table[best_pt]['scenarios']['POINT']['total']:+.5f} ({best_pt})")
    out.update(table=table, best_pessimistic=best, best_point=best_pt, cleared=cleared, verdict=verdict)
    (WORK / "price271.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price271.log").write_text("\n".join(LINES) + "\n")
    # a markdown table for RESULTS
    md = ["## P1. Six-member price over `cd4svd5` (261's recipe; legs at POINT; three coherent scenarios)\n",
          "| cell | Δ`pds` (3 folds) | pooled 95 % CI | Δ`cos_agg` % | delivered | Δ`nmae` id | `pds` | `mse` | DE-three | `nmae` | **PESSIMISTIC** | POINT | OPTIMISTIC | LT3O | cleared |",
          "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for name in MAIN + CONTROLS:
        t = table[name]
        c, s, b = t["vs_cd4svd5"], t["scenarios"], t["boot"]
        po = s["POINT"]
        md.append(f"| `{name}` | {c['d_pds']:+.5f} ({c['pds_pos']}/3) | [{b['lo']:+.4f}, {b['hi']:+.4f}] | {100*c['rel_agg']:+.3f} | {100*(c['delivered_ratio']-1):+.2f} % | "
                  f"{c['d_nmae']:+.4f} | {po['pds']:+.5f} | {po['mse']:+.5f} | {po['de3']:+.5f} | {po['nmae']:+.5f} | **{s['PESSIMISTIC']['total']:+.5f}** | {po['total']:+.5f} | "
                  f"{s['OPTIMISTIC']['total']:+.5f} | {100*t['retention']['retained']:+.0f} % | {t.get('cleared', '—')} |")
    (WORK / "price271.md").write_text("\n".join(md) + "\n")
    Pr(f"\nwrote {WORK/'price271.json'}, price271.log, price271.md")


if __name__ == "__main__":
    main()
