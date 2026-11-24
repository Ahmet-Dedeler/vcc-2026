"""272 -- the six-member price of every noise-corrected-operator cell OVER cd4svd5 on three COHERENT scenarios (READ-FIRST
section 3; 261's recipe re-implemented with the mandate's bands, code-path KA against price257.price), the
known-answer / retrodiction row first (KA-3: cd4svd5 vs STACK through the recipe at h1union's terms must read BR249,
+0.0043 +- 0.0008 overall, and equal 257's R1 to 1e-9).  Nothing is priced unless it passes.

    pds      2.2138 x k_pds x d pds(3 cancer folds, vs cd4svd5) / 6
    mse      d cos_board = k_cos x rel cos_agg x 0.105470;  d raw = -1.996 x d cos_board (+ 243's amplitude leg through
             197's identity at cd4svd5's own terms when the delivered amplitude moves > 0.2 %);  overall = -1.051 x d raw / 6
    DE-three rate x d cos_board, split reach / fid / jac by 218's card
    nmae     216/188's exact identity, -d nmae_wald89 / 0.60 / 6 (three cancer folds concatenated) x multiplier
Bands (the mandate): k_pds {0.85, 1.22}; k_cos {0.94, 3.76} -- the NEW-SOURCE band end (0.94) is the pessimistic
reading of a cosine GAIN and the RE-SHAPE end (3.76) of a cosine LOSS (the operator class is unmeasured on the board
beyond the <= +0.0004 hyperparameter points); DE-three {0.0425, 0.24} coherent with k_cos; identity x1.5 on a loss /
x0 on a gain (pessimistic), x1.0 (point).  POINT = 257's point set (k_pds 1.22, k_cos 3.76, DE-three 0.24 x 0.776,
identity x1.0), the set that retrodicts BR249.
CLEARED = PESSIMISTIC >= +0.0019 AND d pds > 0 on 3/3 folds with the pooled paired-bootstrap interval above 0 AND
leave-top-3-out retention >= 50 % AND every control below half the arm's d pds.

    PYTHONDONTWRITEBYTECODE=1 uv run python price272.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E257 = ROOT / "experiments/257-soft-deflation-price"
sys.path.insert(0, str(E257))
import price257 as P  # noqa: E402  (constants, raw_for, scaled_terms, price, AMP249, BOARD)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
K_PDS_BAND, K_COS_BAND, DE3_BAND = (0.85, 1.22), (0.94, 3.76), (0.0425, 0.24)
POINT = dict(k_pds=1.22, k_cos=3.76, de3=P.DE4_PER_DCOS * P.DE3_SHARE, x_loss=1.0, x_gain=1.0)
GATE = P.GATE_SHAPE
R_BANDS = dict(overall=(0.0043, 0.0008), pds=(0.0055, 0.0006), mse=(-0.0008, 0.0003))
ARMS_A = [f"wiener:k{k:g}" for k in (0.25, 0.5, 1, 2, 4)] + [f"floor:c{c:g}" for c in (0.25, 0.5, 1, 2, 4)] + ["hard:0.5"]
CONTROL_OF = {**{f"wiener:k{k:g}": [f"wiener:k{k:g}:shuf", f"wiener:k{k:g}:mirror"] for k in (0.5, 1, 2)},
              "wiener:k0.25": ["wiener:k0.5:shuf", "wiener:k0.5:mirror"], "wiener:k4": ["wiener:k2:shuf", "wiener:k2:mirror"],
              **{f"floor:c{c:g}": ["floor:c1:shuf", "floor:c1:mirror"] for c in (0.25, 0.5, 1, 2, 4)}, "hard:0.5": ["hard:0.5:mirror"]}
CONTROLS_A = sorted({c for v in CONTROL_OF.values() for c in v})
ARMS_B = [("refit_full", "cd4svd5"), ("op_A", "cd4svd5"), ("op_B", "cd4svd5"), ("xfit_avg", "cd4svd5"), ("half:inc", "cd4svd5"),
          ("half:full", "half:inc"), ("half:avgop", "half:inc"), ("half:xfit", "half:inc"), ("half:infit", "half:inc")]
BASE = "cd4svd5"
LINES: list[str] = []
MD: list[str] = []


def Pr(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def load(prefix: str, folds):
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in folds:
        d = json.loads((WORK / f"{prefix}-{held}.json").read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        tag = "pertarget272" + ("-refit" if prefix.startswith("refit") else "")
        z = np.load(WORK / f"{tag}-{held}.npz")
        for k in z.files:
            pt[(held, k)] = z[k]
        tag = "identity272" + ("-refit" if prefix.startswith("refit") else "")
        z = np.load(WORK / f"{tag}-{held}.npz")
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
        fa, fb_ = a["audit"].get("cd4_fraction"), b["audit"].get("cd4_fraction")
        fr.append(fa / fb_ if (fa and fb_) else float("nan"))
    nm_a, nm_b = nmae_concat(idn, arm, folds), nmae_concat(idn, base, folds)
    return dict(arm=arm, base=base, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)), d_pds=float(np.mean(dp)),
                per_fold_d_pds=[float(x) for x in dp], per_fold_rel_agg=[float(x) for x in ra], per_fold_rel_macro=[float(x) for x in rm],
                pds_pos=int(sum(x > 0 for x in dp)), sigma=float(np.mean(dp)) / P.SIGMA_PDS,
                fold_se=float(np.std(dp, ddof=1) / np.sqrt(len(dp))) if len(dp) > 1 else float("nan"),
                d_qw=float(np.mean(dq)), delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)), fraction_ratio=float(np.nanmean(fr)),
                nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(f, arm)]["top89_precision"] - rows[(f, base)]["top89_precision"] for f in folds])),
                energy_kept=float(np.mean([rows[(f, arm)]["audit"].get("whitened_energy_kept", 1.0) for f in folds])))


def retention(pt, arm, base, folds, targets, n=3) -> dict:
    full, loo, top = [], [], {}
    for f in folds:
        d = pt[(f, f"{arm}|pds")] - pt[(f, f"{base}|pds")]
        order = np.argsort(-d)
        keep = np.ones(d.size, bool)
        keep[order[:n]] = False
        full.append(d.mean())
        loo.append(d[keep].mean())
        top[f] = [f"{targets[f][i]} {d[i]:+.3f}" for i in order[:n]]
    fm, lm = float(np.mean(full)), float(np.mean(loo))
    return dict(full=fm, leave_top_out=lm, retained=(lm / fm if fm else float("nan")), n=n,
                per_fold_retained=[float(l / g) if g else float("nan") for l, g in zip(loo, full)], top_names=top)


def price(c: dict, *, cos_base: float, per_base: dict, k_pds: float, k_cos: float, de3: float, x_loss: float, x_gain: float) -> dict:
    dcos = k_cos * c["rel_agg"] * cos_base
    draw_cos = P.MSE_SLOPE * dcos
    amp = c["delivered_ratio"] - 1.0
    if abs(amp) > P.AMP_LEG_THRESHOLD:
        draw_amp = P.raw_for(P.scaled_terms(per_base, 1.0 + amp), cos_base + dcos) - P.raw_for(per_base, cos_base + dcos)
    else:
        draw_amp = 0.0
    leg_pds = P.PDS_SCALE * k_pds * c["d_pds"] / 6.0
    leg_mse = -P.MSE_SCALE * (draw_cos + draw_amp) / 6.0
    leg_de3 = de3 * dcos
    legs_de3 = {m: leg_de3 * P.DE3_SPLIT[m] for m in P.DE3_SPLIT}
    nm_point = -c["d_nmae"] / P.NMAE_SPAN / 6.0
    leg_nmae = nm_point * (x_loss if nm_point < 0 else x_gain)
    total = leg_pds + leg_mse + leg_de3 + leg_nmae
    return dict(k_pds=k_pds, k_cos=k_cos, de3_rate=de3, x_loss=x_loss, x_gain=x_gain, dcos_board=dcos, draw_cos=draw_cos, draw_amp=draw_amp,
                amp_pct=100 * amp, amp_leg_applied=bool(abs(amp) > P.AMP_LEG_THRESHOLD), pds=leg_pds, mse=leg_mse, de3=leg_de3, **legs_de3,
                nmae=leg_nmae, nmae_identity_point=nm_point, total=total)


def factors(which: str, c: dict) -> dict:
    gain_cos, gain_pds = c["rel_agg"] > 0, c["d_pds"] > 0
    if which == "POINT":
        return dict(POINT)
    if which == "PESSIMISTIC":
        return dict(k_pds=K_PDS_BAND[0] if gain_pds else K_PDS_BAND[1], k_cos=K_COS_BAND[0] if gain_cos else K_COS_BAND[1],
                    de3=DE3_BAND[0] if gain_cos else DE3_BAND[1], x_loss=1.5, x_gain=0.0)
    if which == "OPTIMISTIC":
        return dict(k_pds=K_PDS_BAND[1] if gain_pds else K_PDS_BAND[0], k_cos=K_COS_BAND[1] if gain_cos else K_COS_BAND[0],
                    de3=DE3_BAND[1] if gain_cos else DE3_BAND[0], x_loss=1.0, x_gain=1.0)
    raise ValueError(which)


def scenarios(c: dict, cos_base=P.COS_CD4SVD5, per_base=P.PER_CD4SVD5) -> dict:
    return {w: price(c, cos_base=cos_base, per_base=per_base, **factors(w, c)) for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC")}


def p257_sets(c: dict, cos_base, per_base) -> dict:
    mine = dict(point=price(c, cos_base=cos_base, per_base=per_base, k_pds=1.22, k_cos=3.76, de3=P.DE4_PER_DCOS * P.DE3_SHARE, x_loss=1.0, x_gain=1.0),
                conservative=price(c, cos_base=cos_base, per_base=per_base, k_pds=0.85, k_cos=3.76, de3=P.DE4_PER_DCOS, x_loss=1.5, x_gain=0.0))
    theirs = dict(point=P.price(c, cos_base=cos_base, per_base=per_base, k_pds=1.22, k_cos=3.76, de3_full=False, nmae_mode="point"),
                  conservative=P.price(c, cos_base=cos_base, per_base=per_base, k_pds=0.85, k_cos=3.76, de3_full=True, nmae_mode="conservative"))
    dmax = max(abs(mine[s][k] - theirs[s][k]) for s in mine for k in ("pds", "mse", "de3", "nmae", "total"))
    return dict(mine=mine, p257=theirs, max_abs_diff=dmax)


def main() -> None:
    rows, pt, idn, meta = load("sweep272", CANCER)
    targets = {f: meta[f]["targets"] for f in CANCER}
    ana = json.loads((WORK / "analysis272.json").read_text())
    fails: list[str] = []
    out: dict = dict(recipe=dict(PDS_SCALE=P.PDS_SCALE, MSE_SCALE=P.MSE_SCALE, MSE_SLOPE=P.MSE_SLOPE, COS_CD4SVD5=P.COS_CD4SVD5, COS_H1UNION=P.COS_H1UNION,
                                 POINT=POINT, K_PDS_BAND=K_PDS_BAND, K_COS_BAND=K_COS_BAND, DE3_BAND=DE3_BAND, DE3_SPLIT=P.DE3_SPLIT, NMAE_SPAN=P.NMAE_SPAN,
                                 SIGMA_PDS=P.SIGMA_PDS, GATE=GATE, AMP_LEG_THRESHOLD=P.AMP_LEG_THRESHOLD))
    Pr("=" * 120)
    Pr("272 -- the noise-corrected operator over cd4svd5: six-member price on three COHERENT scenarios (261's recipe; the mandate's bands)")
    Pr("=" * 120)
    Pr(f"recipe: PDS_SCALE {P.PDS_SCALE}, MSE_SCALE {P.MSE_SCALE}, MSE_SLOPE {P.MSE_SLOPE}, cos_cd4svd5 {P.COS_CD4SVD5:.6f}; POINT k_pds {POINT['k_pds']} k_cos {POINT['k_cos']} "
       f"DE-three {POINT['de3']:.4f} identity x1.0; bands k_pds {K_PDS_BAND} k_cos {K_COS_BAND} (0.94 new-source for a GAIN, 3.76 re-shape for a LOSS) DE-three {DE3_BAND} identity x1.5 loss / x0 gain")
    for h in CANCER:
        m = meta[h]
        Pr(f"  {h:14s} {m['n_targets']} targets, {m['n_rows']} rows, {m['seconds']}s, peak {m['peak_rss_gib']:.2f} GiB, known-answer {m['ka_pass']}/{m['ka_total']}, failures {m['failures']}")
        if m["failures"] or m["ka_pass"] != m["ka_total"]:
            fails.append(f"sweep {h}")

    # ---- KA-3 / R1: the retrodiction of BR249 through the recipe ----
    Pr("\n" + "=" * 120)
    Pr("KNOWN-ANSWER AND RETRODICTION ROWS (must pass before any cell is priced)")
    Pr("=" * 120)
    c_r = compare(rows, idn, BASE, "STACK", CANCER)
    p246 = json.loads((ROOT / "experiments/246-direction-hypotheses/work/price246.json").read_text())["B2_cd4_svd5_resolved"]
    ref_agg = P.AMP249["rel_agg_plolo"]
    d_r3 = max(abs(c_r["d_pds"] - p246["d_pds"]), abs(c_r["rel_macro"] - p246["rel_macro"]), abs(c_r["rel_agg"] - ref_agg))
    Pr(f"R3  cd4svd5 vs STACK: d pds {c_r['d_pds']:+.9f} (246 {p246['d_pds']:+.9f})  agg {100*c_r['rel_agg']:+.6f}% (249 {100*ref_agg:+.6f}%)  macro {100*c_r['rel_macro']:+.6f}% "
       f"(246 {100*p246['rel_macro']:+.6f}%)  max|d| {d_r3:.1e}  {'PASS' if d_r3 < P.TOL else 'FAIL'}")
    if d_r3 >= P.TOL:
        fails.append("R3")
    s_r = p257_sets(c_r, P.COS_H1UNION, P.PER_H1UNION)
    r2 = s_r["mine"]["point"]
    band_ok = {k: abs(r2[k] - R_BANDS[k][0]) <= R_BANDS[k][1] for k in ("pds", "mse")}
    band_ok["overall"] = abs(r2["total"] - R_BANDS["overall"][0]) <= R_BANDS["overall"][1]
    board_pds = (P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_pds"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_pds"]) / 6
    board_mse = (P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_mse"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_mse"]) / 6
    board_all = P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_avg"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_avg"]
    p257 = json.loads((E257 / "work/price257.json").read_text())
    ref_r1 = p257["retrodiction"]["R1"]["price"]["point"]["total"]
    r_ok = all(band_ok.values()) and s_r["max_abs_diff"] < 1e-12 and abs(r2["total"] - ref_r1) < 1e-9
    Pr(f"KA-3 BR249 through the recipe (cd4svd5 vs h1union's STACK at h1union's terms, POINT): d cos_board {r2['dcos_board']:+.6f} ({100*r2['dcos_board']/P.COS_H1UNION:+.3f}%; board -1.735 %); "
       f"d raw mse {r2['draw_cos']+r2['draw_amp']:+.5f} (board +{P.RAW_MSE_CD4SVD5-P.RAW_MSE_H1UNION:.5f}); pds {r2['pds']:+.5f} (board {board_pds:+.5f}) mse {r2['mse']:+.5f} (board {board_mse:+.5f}) "
       f"DE-three {r2['de3']:+.5f} nmae {r2['nmae']:+.5f} = OVERALL {r2['total']:+.5f} (band +0.0043 +- 0.0008; board {board_all:+.7f}; 257 R1 {ref_r1:+.5f}, |d| {abs(r2['total']-ref_r1):.1e}); "
       f"this file's price() == price257.price(): max|d| {s_r['max_abs_diff']:.1e}  -> {'PASS' if r_ok else 'FAIL'} {band_ok}")
    if not r_ok:
        fails.append("KA-3")
    out["retrodiction"] = dict(compare=c_r, price=s_r["mine"], bands=R_BANDS, band_ok=band_ok, board=dict(pds=board_pds, mse=board_mse, overall=board_all),
                               codepath_max_abs_diff=s_r["max_abs_diff"], ref_257_R1=ref_r1, R3_max_abs_diff=d_r3, ok=r_ok)
    if fails:
        Pr(f"\n*** KNOWN-ANSWER / RETRODICTION FAILURES: {fails} -- no cell is priced")
        (WORK / "price272.json").write_text(json.dumps(dict(out, failures=fails), indent=1, default=float))
        (WORK / "price272.log").write_text("\n".join(LINES) + "\n")
        raise SystemExit(1)

    # ---- construction A ----
    table: dict[str, dict] = {}
    for name in ARMS_A + CONTROLS_A:
        c = compare(rows, idn, name, BASE, CANCER)
        sc = scenarios(c)
        ret = retention(pt, name, BASE, CANCER, targets)
        pooled = ana["arms"][name]["pooled"]
        table[name] = dict(vs_cd4svd5=c, scenarios=sc, retention=ret, pooled=pooled, shape_only=abs(c["delivered_ratio"] - 1) <= P.AMP_SHAPE_BOUND,
                           cond_gate=bool(sc["PESSIMISTIC"]["total"] >= GATE), cond_retention=bool(ret["retained"] >= 0.5),
                           cond_pds=bool(c["pds_pos"] == 3 and pooled["lo"] > 0))
    for name in ARMS_A:
        arm = table[name]["vs_cd4svd5"]["d_pds"]
        ctl = {n: table[n]["vs_cd4svd5"]["d_pds"] for n in CONTROL_OF[name]}
        table[name]["controls_d_pds"] = ctl
        table[name]["cond_controls"] = bool(arm > 0 and all(v < 0.5 * arm for v in ctl.values()))
        t = table[name]
        t["cleared"] = bool(t["cond_gate"] and t["cond_retention"] and t["cond_pds"] and t["cond_controls"])

    Pr("\n" + "=" * 120)
    Pr("CONSTRUCTION A -- every cell as an INCREMENT over cd4svd5 (three cancer folds); PESSIMISTIC / POINT / OPTIMISTIC totals (legs at POINT)")
    Pr("=" * 120)
    hdr = (f"{'cell':20s}{'kept':>6s}{'d pds':>10s}{'folds':>6s}{'CI lo':>9s}{'CI hi':>9s}{'agg%':>8s}{'macro%':>8s}{'dlv%':>7s}{'d nmae':>8s}"
           f"{'pds':>9s}{'mse':>9s}{'DE3':>9s}{'nmae':>9s}{'PESS':>9s}{'POINT':>9s}{'OPT':>9s}{'ret3%':>7s}{'gate':>5s}{'pds':>4s}{'ret':>4s}{'ctl':>4s}")
    Pr(hdr)
    MD.append("## P1. Six-member price over `cd4svd5` (three cancer folds; legs at POINT; PESSIMISTIC / POINT / OPTIMISTIC totals; CLEARED = pessimistic ≥ +0.0019 AND `pds` 3/3 with the pooled interval above 0 AND retention ≥ 50 % AND every control below half)\n")
    MD.append("| cell | whitened energy kept | Δ`pds` (3 folds) | pooled 95 % CI | Δ`cos_agg` % | delivered | Δ`nmae` id | `pds` | `mse` | DE-three | `nmae` | **PESSIMISTIC** | POINT | OPTIMISTIC | LT3O | gate / pds / ret / ctl | cleared |")
    MD.append("|---|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")

    def line(name, is_arm=True):
        t = table[name]
        c, s = t["vs_cd4svd5"], t["scenarios"]
        pe, po, op = s["PESSIMISTIC"], s["POINT"], s["OPTIMISTIC"]
        Pr(f"{name:20s}{100*c['energy_kept']:6.1f}{c['d_pds']:+10.6f}{c['pds_pos']:>4d}/3{t['pooled']['lo']:+9.4f}{t['pooled']['hi']:+9.4f}{100*c['rel_agg']:+8.3f}{100*c['rel_macro']:+8.3f}"
           f"{100*(c['delivered_ratio']-1):+7.3f}{c['d_nmae']:+8.4f}{po['pds']:+9.5f}{po['mse']:+9.5f}{po['de3']:+9.5f}{po['nmae']:+9.5f}{pe['total']:+9.5f}{po['total']:+9.5f}{op['total']:+9.5f}"
           f"{100*t['retention']['retained']:+7.0f}{str(t['cond_gate'])[0]:>5s}{str(t['cond_pds'])[0]:>4s}{str(t['cond_retention'])[0]:>4s}{str(t.get('cond_controls', '-'))[0]:>4s}")
        flags = f"{'T' if t['cond_gate'] else 'F'} / {'T' if t['cond_pds'] else 'F'} / {'T' if t['cond_retention'] else 'F'} / {('T' if t['cond_controls'] else 'F') if is_arm else '—'}"
        MD.append(f"| `{name}` | {100*c['energy_kept']:.0f} % | {c['d_pds']:+.5f} ({c['pds_pos']}/3) | [{t['pooled']['lo']:+.4f}, {t['pooled']['hi']:+.4f}] | {100*c['rel_agg']:+.2f} | {100*(c['delivered_ratio']-1):+.2f} % | "
                  f"{c['d_nmae']:+.4f} | {po['pds']:+.5f} | {po['mse']:+.5f} | {po['de3']:+.5f} | {po['nmae']:+.5f} | **{pe['total']:+.5f}** | {po['total']:+.5f} | {op['total']:+.5f} | "
                  f"{100*t['retention']['retained']:+.0f} % | {flags} | {'**no**' if is_arm and not t['cleared'] else ('YES' if is_arm else 'control')} |")

    for name in ARMS_A:
        line(name)
    Pr("  -- controls --")
    for name in CONTROLS_A:
        line(name, is_arm=False)
    n_gate = sum(table[n]["cond_gate"] for n in ARMS_A)
    n_point = sum(table[n]["scenarios"]["POINT"]["total"] >= GATE for n in ARMS_A)
    n_opt = sum(table[n]["scenarios"]["OPTIMISTIC"]["total"] >= GATE for n in ARMS_A)
    cleared = [n for n in ARMS_A if table[n]["cleared"]]
    best = max(ARMS_A, key=lambda n: table[n]["scenarios"]["PESSIMISTIC"]["total"])
    best_pt = max(ARMS_A, key=lambda n: table[n]["scenarios"]["POINT"]["total"])
    Pr(f"\nTHE SCREEN (A): {len(ARMS_A)} cells.  PESSIMISTIC >= +{GATE}: {n_gate}/{len(ARMS_A)};  POINT: {n_point};  OPTIMISTIC: {n_opt};  pds 3/3 with CI > 0: "
       f"{sum(table[n]['cond_pds'] for n in ARMS_A)};  ALL conditions: {len(cleared)} {cleared or ''}")
    Pr(f"best PESSIMISTIC {best} {table[best]['scenarios']['PESSIMISTIC']['total']:+.5f} (point {table[best]['scenarios']['POINT']['total']:+.5f}); best POINT {best_pt} {table[best_pt]['scenarios']['POINT']['total']:+.5f}")
    for tag, name in (("best PESSIMISTIC", best), ("best POINT", best_pt)):
        t = table[name]
        c, s = t["vs_cd4svd5"], t["scenarios"]
        Pr(f"\n{tag}: {name}: d pds {c['d_pds']:+.6f} ({c['pds_pos']}/3 {['%+.5f' % x for x in c['per_fold_d_pds']]}; pooled [{t['pooled']['lo']:+.4f}, {t['pooled']['hi']:+.4f}])  agg {100*c['rel_agg']:+.3f}% "
           f"(per fold {['%+.3f%%' % (100*x) for x in c['per_fold_rel_agg']]})  macro {100*c['rel_macro']:+.3f}%  delivered {100*(c['delivered_ratio']-1):+.3f}% (shape-only {t['shape_only']})  d nmae89 {c['d_nmae']:+.5f}  d prec89 {c['d_prec89']:+.5f}")
        for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            p = s[w]
            Pr(f"  {w:11s} k_pds {p['k_pds']:.2f} k_cos {p['k_cos']:.2f} DE3 {p['de3_rate']:.4f} id x{p['x_loss']:.1f}/{p['x_gain']:.1f}:  pds {p['pds']:+.5f}  mse {p['mse']:+.5f} (d cos_board {p['dcos_board']:+.6f}, "
               f"d raw {p['draw_cos']+p['draw_amp']:+.5f}, amp leg {'on' if p['amp_leg_applied'] else 'off'})  reach {p['reach']:+.5f} fid {p['fid']:+.5f} jac {p['jac']:+.5f}  nmae {p['nmae']:+.5f}  = {p['total']:+.5f}  ({p['total']/GATE:+.2f}x gate)")
        r = t["retention"]
        Pr(f"  leave-top-3-out retention {100*r['retained']:.1f}% (per fold {['%.0f%%' % (100*x) for x in r['per_fold_retained']]}); top-3: " + "; ".join(f"{h}: {', '.join(r['top_names'][h])}" for h in CANCER))
        Pr(f"  controls d pds: {', '.join(f'{n} {v:+.6f}' for n, v in t['controls_d_pds'].items())} (half the arm {0.5*c['d_pds']:+.6f}; below-half {t['cond_controls']})")
        Pr(f"  conditions: gate {t['cond_gate']}; pds 3/3 & CI > 0 {t['cond_pds']}; retention {t['cond_retention']}; controls {t['cond_controls']} -> cleared {t['cleared']}")
    out.update(A=table, cleared=cleared, best_pessimistic=best, best_point=best_pt, n_gate_pessimistic=n_gate, n_gate_point=n_point, n_gate_optimistic=n_opt)

    # ---- construction B (one fold; priced for the record, cannot meet the 3/3 rule) ----
    refit_folds = [f for f in CANCER if (WORK / f"refit272-{f}.json").exists()]
    if refit_folds:
        rrows, rpt, ridn, rmeta = load("refit272", refit_folds)
        rt = {f: rmeta[f]["targets"] for f in refit_folds}
        out["B"], out["B_meta"] = {}, {}
        for rf in refit_folds:
            Pr("\n" + "=" * 120)
            Pr(f"CONSTRUCTION B ({rf}): the refitted / cross-fitted operators, priced on the ONE fold (cannot meet the 3/3 rule; for the record)")
            Pr("=" * 120)
            MD.append(f"\n## P2 ({rf}). Construction B on `{rf}` (one fold; priced for the record, three-fold conditions not testable): legs at POINT\n")
            MD.append("| cell | vs | Δ`pds` | Δ`cos_agg` % | delivered | Δ`nmae` id | `pds` | `mse` | DE-three | `nmae` | **PESSIMISTIC** | POINT | OPTIMISTIC | LT3O |")
            MD.append("|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
            tb = {}
            for name, base in ARMS_B:
                c = compare(rrows, ridn, name, base, (rf,))
                sc = scenarios(c)
                ret = retention(rpt, name, base, (rf,), rt)
                tb[f"{name}|{base}"] = dict(vs=c, scenarios=sc, retention=ret)
                pe, po, op = sc["PESSIMISTIC"], sc["POINT"], sc["OPTIMISTIC"]
                Pr(f"  {name:12s} vs {base:10s} d pds {c['d_pds']:+.6f}  agg {100*c['rel_agg']:+.3f}%  macro {100*c['rel_macro']:+.3f}%  dlv {100*(c['delivered_ratio']-1):+.3f}%  d nmae {c['d_nmae']:+.5f}  "
                   f"PESS {pe['total']:+.5f}  POINT {po['total']:+.5f}  OPT {op['total']:+.5f}  ret3 {100*ret['retained']:+.0f}%")
                MD.append(f"| `{name}` | `{base}` | {c['d_pds']:+.5f} | {100*c['rel_agg']:+.2f} | {100*(c['delivered_ratio']-1):+.2f} % | {c['d_nmae']:+.4f} | {po['pds']:+.5f} | {po['mse']:+.5f} | {po['de3']:+.5f} | "
                          f"{po['nmae']:+.5f} | **{pe['total']:+.5f}** | {po['total']:+.5f} | {op['total']:+.5f} | {100*ret['retained']:+.0f} % |")
            out["B"][rf] = tb
            out["B_meta"][rf] = dict(ka_pass=rmeta[rf]["ka_pass"], ka_total=rmeta[rf]["ka_total"], failures=rmeta[rf]["failures"], fit=rmeta[rf]["fit"])

    verdict = "CLEARS" if cleared else "DOES NOT CLEAR"
    Pr(f"\nVERDICT: {verdict} -- {len(cleared)} of {len(ARMS_A)} construction-A cells meet all four conditions; best pessimistic {table[best]['scenarios']['PESSIMISTIC']['total']:+.5f} ({best}), "
       f"best point {table[best_pt]['scenarios']['POINT']['total']:+.5f} ({best_pt})")
    out.update(verdict=verdict)
    (WORK / "price272.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price272.log").write_text("\n".join(LINES) + "\n")
    (WORK / "price272.md").write_text("\n".join(MD) + "\n")
    Pr(f"\nwrote {WORK/'price272.json'}, price272.log, price272.md")


if __name__ == "__main__":
    main()
