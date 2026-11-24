"""257 -- the six-member price of every softer-deflation cell OVER `cd4svd5`, on the transfer
factors the 2026-09-05 board point measured (BR249 section 2).  The recipe is the mandate's, fixed
here before any grid number was read; three retrodiction rows must pass before a cell is priced.

    pds   overall = PDS_SCALE * k_pds * d_pds(3 cancer folds, vs cd4svd5) / 6,  k_pds 1.22 point / 0.85 conservative
    mse   d cos_board = k_cos * rel_cos_agg(PLOLO, vs cd4svd5) * cos_cd4svd5 (0.10547),  k_cos 3.76 (re-shape class;
          2.15 macro read and 1.36 re-dose reported as sensitivity);  d raw = MSE_SLOPE (-1.996) * d cos_board
          (+ 243's amplitude leg, 197's identity at cd4svd5's own terms, when the delivered amplitude moves > 0.2 %);
          overall = -MSE_SCALE (1.051) * d raw / 6
    DE-three (fid, reach, jac)  overall = 0.24 * d cos_board * share  (share = the three members' part of 218's
          DE-four rate card, 0.776; conservative: the whole 0.24)
    nmae  216/188's exact identity (wald89 pool, the three cancer folds concatenated), x1.0 point;
          conservative x1.5 on a loss, 0 credit on a gain (BR237: the identity is a LOWER bound on the loss)
Retrodictions: R1 cd4svd5 vs h1union reads +0.0043 +- 0.0006 overall, pds +0.0055 +- 0.0006, mse -0.0008 +- 0.0003;
R2 k5_s1 == 246 B2_cd4_svd5_resolved (d pds +0.012262757, agg -0.4617 %, macro -0.8087 %) to 1e-9;
R3 k10_s1 vs k5_s1 == 255's increment to 1e-9 and prices NEGATIVE at the point factors.

    uv run python price257.py
"""
from __future__ import annotations

import json
import resource
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E246 = ROOT / "experiments/246-direction-hypotheses"
E249 = ROOT / "experiments/249-cd4-deflate"
E255 = ROOT / "experiments/255-cd4-deflate-k10"

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)

# ----------------------------------------------------------------- the recipe (MANDATE.md), fixed
PDS_SCALE = 2.2138                       # board scaled `pds` per raw `pds`
K_PDS = {"point": 1.22, "conservative": 0.85}
MSE_SCALE = 1.051                        # board scaled `mse` per raw `mse` (scaled falls as raw rises)
MSE_SLOPE = -1.996                       # d raw mse / d cos at cd4svd5's amplitude terms (197's identity, 249 section 9)
K_COS = {"point": 3.76, "macro_read": 2.15, "re_dose": 1.36}
COS_H1UNION = 0.107332                   # 243: h1union's board cosine, MEASURED off its own A_c/B_c
COS_CD4SVD5 = COS_H1UNION * (1.0 - 0.01735)   # BR249: -1.735 % -> 0.105470
DE4_PER_DCOS = 0.24                      # 243 section 4.3: DE-four overall per unit board d cos (the union point)
RATE = dict(mse=0.003634, reach=0.002446, nmae=0.000980, fid=0.000847, jac=0.000100)   # 218 section 7.4, SHARES only
DE3_SHARE = (RATE["reach"] + RATE["fid"] + RATE["jac"]) / (RATE["reach"] + RATE["fid"] + RATE["jac"] + RATE["nmae"])
DE3_SPLIT = {m: RATE[m] / (RATE["reach"] + RATE["fid"] + RATE["jac"]) for m in ("reach", "fid", "jac")}
NMAE_SPAN = 0.60
NMAE_CONS = 1.5
SIGMA_PDS = 0.001494                     # 162: LOO RMS on the raw PLOLO `pds` delta
GATE_SHAPE, GATE_AMPL = 0.0019, 0.0089
AMP_LEG_THRESHOLD = 0.002                # add 243's amplitude leg beyond this delivered-amplitude move
AMP_SHAPE_BOUND = 0.005                  # shape-only class within this
SUMD = 0.22914103411323902               # 197's truth energy (lib243 SUMD)
TOL = 1e-9
R1_BANDS = dict(overall=(0.0043, 0.0006), pds=(0.0055, 0.0006), mse=(-0.0008, 0.0003))

AMP249 = json.loads((E249 / "work/amp249.json").read_text())
PER_CD4SVD5 = AMP249["per_arm"]          # cd4svd5's own (A_c, B_c, Kc_c, cap_c) off the package (249 section 9)
PER_H1UNION = AMP249["per_inc"]
BOARD = {eid: json.loads((ROOT / f"data/board-entries/{eid}.json").read_text())["entry"]
         for eid in ("iUziqE1gewfcPKzYUDaB", "Xv0ayvuoHpsKYWv2iu42")}
RAW_MSE_H1UNION = BOARD["iUziqE1gewfcPKzYUDaB"]["expr_mse_unbiased_capped_norm"]
RAW_MSE_CD4SVD5 = BOARD["Xv0ayvuoHpsKYWv2iu42"]["expr_mse_unbiased_capped_norm"]
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


# ----------------------------------------------------------------- 197's identity at an arm's own terms
def raw_for(per: dict, cos: float) -> float:
    num = float(np.mean([per[c]["A"] + per[c]["Kc"] - min(per[c]["cap"], per[c]["B"]) for c in "ABC"]))
    msa = float(np.mean([np.sqrt(per[c]["A"]) for c in "ABC"]))
    return 1.0 + (num - 2.0 * cos * np.sqrt(SUMD) * msa) / SUMD


def scaled_terms(per: dict, ratio: float) -> dict:
    """A global delivered-amplitude rescale by `ratio`: A_c and B_c scale as ratio^2 (249 section 9
    measured the two moving together, -0.409 % / -0.410 %); Kc and the cap are reference constants."""
    return {c: dict(per[c], A=per[c]["A"] * ratio**2, B=per[c]["B"] * ratio**2) for c in "ABC"}


def exact_slope(per: dict) -> float:
    return -2.0 * float(np.mean([np.sqrt(per[c]["A"]) for c in "ABC"])) / np.sqrt(SUMD)


# ----------------------------------------------------------------- load the sweeps
def load():
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in FOLDS:
        d = json.loads((WORK / f"sweep257-{held}.json").read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        z = np.load(WORK / f"pertarget257-{held}.npz")
        for k in z.files:
            pt[(held, k)] = z[k]
        z = np.load(WORK / f"identity257-{held}.npz")
        for k in z.files:
            idn[(held, k)] = z[k]
    return rows, pt, idn, meta


def nmae_concat(idn, arm, folds) -> float:
    U = np.concatenate([idn[(f, f"{arm}|U")] for f in folds])
    M = np.concatenate([idn[(f, f"{arm}|M")] for f in folds])
    D = np.concatenate([idn[(f, f"{arm}|D")] for f in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(rows, idn, arm, base, folds) -> dict:
    """Everything the recipe reads, for `arm` against `base`, mean of per-fold ratios / deltas."""
    ra, rm, dp, dq, dn_a, dn_b, fr = [], [], [], [], [], [], []
    for f in folds:
        a, b = rows[(f, arm)], rows[(f, base)]
        ra.append(a["cos_agg"] / b["cos_agg"] - 1.0)
        rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
        dp.append(a["pds"] - b["pds"])
        dq.append(a["a15_sgn_w"] - b["a15_sgn_w"])
        dn_a.append(a["delivered_norm"])
        dn_b.append(b["delivered_norm"])
        fa = a["audit"].get("fraction", a["audit"].get("frac_cd4"))
        fb = b["audit"].get("fraction", b["audit"].get("frac_cd4"))
        fr.append(fa / fb if (fa and fb) else float("nan"))
    nm_a, nm_b = nmae_concat(idn, arm, folds), nmae_concat(idn, base, folds)
    return dict(arm=arm, base=base, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)),
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp],
                per_fold_rel_agg=[float(x) for x in ra], per_fold_rel_macro=[float(x) for x in rm],
                pds_pos=int(sum(x > 0 for x in dp)), sigma=float(np.mean(dp)) / SIGMA_PDS,
                d_qw=float(np.mean(dq)), delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)),
                fraction_ratio=float(np.nanmean(fr)), nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(f, arm)]["top89_precision"] - rows[(f, base)]["top89_precision"]
                                        for f in folds])),
                energy_removed=float(np.mean([rows[(f, arm)]["audit"].get("energy_removed_frac", 0.0) for f in folds])))


# ----------------------------------------------------------------- the price
def price(c: dict, *, cos_base: float, per_base: dict, k_pds: float, k_cos: float,
          de3_full: bool, nmae_mode: str, mse_form: str = "mandate") -> dict:
    dcos = k_cos * c["rel_agg"] * cos_base
    if mse_form == "mandate":
        draw_cos = MSE_SLOPE * dcos
    else:
        draw_cos = raw_for(per_base, cos_base + dcos) - raw_for(per_base, cos_base)
    amp = c["delivered_ratio"] - 1.0
    if abs(amp) > AMP_LEG_THRESHOLD:
        draw_amp = raw_for(scaled_terms(per_base, c["delivered_ratio"]), cos_base + dcos) - raw_for(per_base, cos_base + dcos)
    else:
        draw_amp = 0.0
    leg_pds = PDS_SCALE * k_pds * c["d_pds"] / 6.0
    leg_mse = -MSE_SCALE * (draw_cos + draw_amp) / 6.0
    de3 = DE4_PER_DCOS * dcos * (1.0 if de3_full else DE3_SHARE)
    legs_de3 = {m: de3 * DE3_SPLIT[m] for m in DE3_SPLIT}
    nm_point = -c["d_nmae"] / NMAE_SPAN / 6.0
    if nmae_mode == "point":
        leg_nmae = nm_point
    else:
        leg_nmae = nm_point * NMAE_CONS if nm_point < 0 else 0.0
    total = leg_pds + leg_mse + de3 + leg_nmae
    return dict(k_pds=k_pds, k_cos=k_cos, dcos_board=dcos, draw_cos=draw_cos, draw_amp=draw_amp,
                amp_pct=100 * amp, amp_leg_applied=bool(abs(amp) > AMP_LEG_THRESHOLD),
                pds=leg_pds, mse=leg_mse, de3=de3, **legs_de3, nmae=leg_nmae, nmae_identity_point=nm_point,
                total=total, mse_form=mse_form, de3_full=de3_full, nmae_mode=nmae_mode)


def price_both(c, cos_base, per_base):
    pt = price(c, cos_base=cos_base, per_base=per_base, k_pds=K_PDS["point"], k_cos=K_COS["point"],
               de3_full=False, nmae_mode="point")
    cn = price(c, cos_base=cos_base, per_base=per_base, k_pds=K_PDS["conservative"], k_cos=K_COS["point"],
               de3_full=True, nmae_mode="conservative")
    s215 = price(c, cos_base=cos_base, per_base=per_base, k_pds=K_PDS["point"], k_cos=K_COS["macro_read"],
                 de3_full=False, nmae_mode="point")
    s136 = price(c, cos_base=cos_base, per_base=per_base, k_pds=K_PDS["point"], k_cos=K_COS["re_dose"],
                 de3_full=False, nmae_mode="point")
    ex = price(c, cos_base=cos_base, per_base=per_base, k_pds=K_PDS["point"], k_cos=K_COS["point"],
               de3_full=False, nmae_mode="point", mse_form="exact")
    return dict(point=pt, conservative=cn, k_cos_2_15=s215, k_cos_1_36=s136, exact_identity=ex)


def retention(pt, arm, base, folds, n=3) -> dict:
    """246 analyse246's leave-top-n-out retention of the per-target d pds (mean over folds of
    mean(d without the n largest) / mean(d)), each fold on its own ordering."""
    full, loo = [], []
    for f in folds:
        d = pt[(f, f"{arm}|pds")] - pt[(f, f"{base}|pds")]
        keep = np.ones(d.size, bool)
        keep[np.argsort(-d)[:n]] = False
        full.append(d.mean())
        loo.append(d[keep].mean())
    fm, lm = float(np.mean(full)), float(np.mean(loo))
    return dict(full=fm, leave_top_out=lm, retained=(lm / fm if fm else float("nan")), n=n,
                per_fold_retained=[float(l / f) if f else float("nan") for l, f in zip(loo, full)])


def main() -> None:
    rows, pt, idn, meta = load()
    out: dict = dict(recipe=dict(PDS_SCALE=PDS_SCALE, K_PDS=K_PDS, MSE_SCALE=MSE_SCALE, MSE_SLOPE=MSE_SLOPE, K_COS=K_COS,
                                 COS_H1UNION=COS_H1UNION, COS_CD4SVD5=COS_CD4SVD5, DE4_PER_DCOS=DE4_PER_DCOS,
                                 DE3_SHARE=DE3_SHARE, NMAE_SPAN=NMAE_SPAN, NMAE_CONS=NMAE_CONS, SIGMA_PDS=SIGMA_PDS,
                                 GATE_SHAPE=GATE_SHAPE, GATE_AMPL=GATE_AMPL, AMP_LEG_THRESHOLD=AMP_LEG_THRESHOLD,
                                 AMP_SHAPE_BOUND=AMP_SHAPE_BOUND, SUMD=SUMD,
                                 exact_slope_cd4svd5=exact_slope(PER_CD4SVD5), exact_slope_h1union=exact_slope(PER_H1UNION)))
    fails: list[str] = []
    P(f"257 price -- recipe: PDS_SCALE {PDS_SCALE}, k_pds {K_PDS}, MSE_SCALE {MSE_SCALE}, MSE_SLOPE {MSE_SLOPE} "
      f"(exact at cd4svd5's terms {exact_slope(PER_CD4SVD5):.4f}, at h1union's {exact_slope(PER_H1UNION):.4f}), "
      f"k_cos {K_COS}, cos_cd4svd5 {COS_CD4SVD5:.6f}, DE-four/d cos {DE4_PER_DCOS} (DE-three share {DE3_SHARE:.3f}), "
      f"nmae identity x1.0 / x{NMAE_CONS} on a loss, sigma_pds {SIGMA_PDS}")
    fold_desc = ["%s: %d rows, %.1fs, peak %.2f GiB, %d KA failures" % (
        h, meta[h]["n_rows"], meta[h]["seconds"], meta[h]["peak_rss_gib"], len(meta[h]["failures"])) for h in FOLDS]
    P("folds: " + "; ".join(fold_desc))
    for h in FOLDS:
        if meta[h]["failures"]:
            fails.append(f"sweep {h}: {meta[h]['failures']}")
    ka_sweep = sum(1 for (h, a), r in rows.items() if "ka" in r)
    ka_sweep_ok = sum(1 for (h, a), r in rows.items() if "ka" in r and r["ka"]["ok"])
    P(f"sweep known-answer rows: {ka_sweep_ok}/{ka_sweep} at < 1e-9")
    if ka_sweep_ok != ka_sweep:
        fails.append("sweep KA rows failed")

    # ---- six246's identity rows (4 folds concatenated, five known answers at 5e-6) -----------
    six = json.loads((E246 / "work/six246.json").read_text())["arms"]
    dn_a15 = float(np.mean([rows[(f, "a15")]["delivered_norm"] for f in FOLDS]))
    P("\nKNOWN-ANSWER -- 188's exact nmae identity / precision / delivered norm vs six246.json (4 folds):")
    ka_six = []
    for mine, theirs in (("a15", "a15 f=0.10 SHIPPED"), ("STACK", "STACK h1union g=0.28"), ("k5_s1", "B2_cd4_svd5_resolved")):
        got = dict(nmae_wald89=nmae_concat(idn, mine, FOLDS),
                   top89_precision=float(np.mean([rows[(f, mine)]["top89_precision"] for f in FOLDS])),
                   delivered_norm_ratio=float(np.mean([rows[(f, mine)]["delivered_norm"] for f in FOLDS])) / dn_a15)
        for key, g in got.items():
            w = six[theirs][key]
            ok = abs(g - w) <= 5e-6
            ka_six.append(dict(arm=mine, key=key, got=g, want=w, ok=bool(ok)))
            P(f"  {'PASS' if ok else 'FAIL'}  {mine:8s} {key:22s} got {g:.8f} want {w:.8f} |d| {abs(g-w):.1e}")
            if not ok:
                fails.append(f"six246 KA {mine}/{key}")
    out["ka_six246"] = ka_six

    # ================================================================= RETRODICTION ROWS
    P("\n" + "=" * 110)
    P("RETRODICTION ROWS (must pass before any cell is priced)")
    P("=" * 110)
    # R2 -- k5_s1 (cd4svd5) vs STACK (h1union) == 246 / 249 to the digit
    c_r2 = compare(rows, idn, "k5_s1", "STACK", CANCER)
    p246 = json.loads((E246 / "work/price246.json").read_text())["B2_cd4_svd5_resolved"]
    ref_agg = AMP249["rel_agg_plolo"]
    d_r2 = max(abs(c_r2["d_pds"] - p246["d_pds"]), abs(c_r2["rel_macro"] - p246["rel_macro"]), abs(c_r2["rel_agg"] - ref_agg))
    r2_ok = d_r2 < TOL
    P(f"R2  k5_s1 vs STACK: d pds {c_r2['d_pds']:+.9f} (246 {p246['d_pds']:+.9f})  agg {100*c_r2['rel_agg']:+.6f}% "
      f"(249 {100*ref_agg:+.6f}%)  macro {100*c_r2['rel_macro']:+.6f}% (246 {100*p246['rel_macro']:+.6f}%)  "
      f"max|d| {d_r2:.1e}  {'PASS' if r2_ok else 'FAIL'};  pds {c_r2['pds_pos']}/3 folds {['%+.5f' % x for x in c_r2['per_fold_d_pds']]}, "
      f"{c_r2['sigma']:.1f} sigma; energy removed {c_r2['energy_removed']:.4f}; fraction x{c_r2['fraction_ratio']:.4f}; "
      f"delivered ratio {100*(c_r2['delivered_ratio']-1):+.3f}% (package: {AMP249['amplitude_pct']:+.3f}%)")
    if not r2_ok:
        fails.append("R2")

    # R1 -- the same increment through the recipe reads the board point
    pr1 = price_both(c_r2, COS_H1UNION, PER_H1UNION)
    r1 = pr1["point"]
    band_ok = {k: abs(r1[k] - R1_BANDS[k][0]) <= R1_BANDS[k][1] for k in ("pds", "mse")}
    band_ok["overall"] = abs(r1["total"] - R1_BANDS["overall"][0]) <= R1_BANDS["overall"][1]
    r1_ok = all(band_ok.values())
    board_pds = (BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_pds"] - BOARD["iUziqE1gewfcPKzYUDaB"]["score_pds"]) / 6
    board_mse = (BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_mse"] - BOARD["iUziqE1gewfcPKzYUDaB"]["score_mse"]) / 6
    board_all = BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_avg"] - BOARD["iUziqE1gewfcPKzYUDaB"]["score_avg"]
    P(f"R1  cd4svd5 vs h1union through the recipe (base cos {COS_H1UNION}, h1union's own terms, k_cos 3.76, k_pds 1.22):")
    P(f"      d cos_board {r1['dcos_board']:+.6f} ({100*r1['dcos_board']/COS_H1UNION:+.3f}% of cos; board -1.735 %)  "
      f"d raw mse: cosine {r1['draw_cos']:+.5f} + amplitude {r1['draw_amp']:+.5f} (leg {'applied' if r1['amp_leg_applied'] else 'not applied'} at "
      f"{r1['amp_pct']:+.3f}%) = {r1['draw_cos']+r1['draw_amp']:+.5f}  (board +{RAW_MSE_CD4SVD5-RAW_MSE_H1UNION:.5f})")
    P(f"      pds {r1['pds']:+.5f} (band {R1_BANDS['pds'][0]:+.4f}+-{R1_BANDS['pds'][1]}; board {board_pds:+.5f})  "
      f"mse {r1['mse']:+.5f} (band {R1_BANDS['mse'][0]:+.4f}+-{R1_BANDS['mse'][1]}; board {board_mse:+.5f})  "
      f"DE-three {r1['de3']:+.5f} (reach {r1['reach']:+.5f} fid {r1['fid']:+.5f} jac {r1['jac']:+.5f}; board -0.00034)  "
      f"nmae {r1['nmae']:+.5f} (identity; board -0.00005)")
    P(f"      OVERALL {r1['total']:+.5f} (band {R1_BANDS['overall'][0]:+.4f}+-{R1_BANDS['overall'][1]}; board {board_all:+.7f})  "
      f"-> {'PASS' if r1_ok else 'FAIL'} {band_ok}")
    P(f"      conservative factors: {pr1['conservative']['total']:+.5f};  exact identity form: mse {pr1['exact_identity']['mse']:+.5f}, "
      f"overall {pr1['exact_identity']['total']:+.5f};  at k_cos 2.15 / 1.36: {pr1['k_cos_2_15']['total']:+.5f} / {pr1['k_cos_1_36']['total']:+.5f}")
    # the identity at the PACKAGE's measured amplitude (-0.409 %) reproduces the board raw mse (how 3.76 was read)
    cos_svd5 = COS_H1UNION * (1 + K_COS["point"] * c_r2["rel_agg"])
    raw_pred = raw_for(scaled_terms(PER_H1UNION, AMP249["amplitude_ratio"]), cos_svd5)
    raw_pred_own = raw_for(PER_CD4SVD5, cos_svd5)
    P(f"      check: 197's identity at k_cos 3.76 with the package's own terms -> raw mse {raw_pred_own:.5f} "
      f"(h1union's terms rescaled by the package's {AMP249['amplitude_pct']:+.3f}%: {raw_pred:.5f}); board {RAW_MSE_CD4SVD5:.5f}; "
      f"h1union's own terms reproduce its board raw {raw_for(PER_H1UNION, COS_H1UNION):.7f} vs {RAW_MSE_H1UNION:.7f}")
    if not r1_ok:
        fails.append("R1")

    # R3 -- k10_s1 vs k5_s1 == 255's increment, and prices negative at the point factors
    c_r3 = compare(rows, idn, "k10_s1", "k5_s1", CANCER)
    s255 = json.loads((E255 / "work/repro255-summary.json").read_text())["increment"]
    d_r3 = max(abs(c_r3["d_pds"] - s255["d_pds"]), abs(c_r3["rel_agg"] - s255["rel_agg"]), abs(c_r3["rel_macro"] - s255["rel_macro"]))
    pr3 = price_both(c_r3, COS_CD4SVD5, PER_CD4SVD5)
    r3_ok = d_r3 < TOL and pr3["point"]["total"] < 0
    P(f"R3  k10_s1 vs k5_s1: d pds {c_r3['d_pds']:+.9f} (255 {s255['d_pds']:+.9f}, {c_r3['pds_pos']}/3 folds, {c_r3['sigma']:.2f} sigma)  "
      f"agg {100*c_r3['rel_agg']:+.4f}% (255 {100*s255['rel_agg']:+.4f}%)  macro {100*c_r3['rel_macro']:+.4f}%  max|d| {d_r3:.1e}")
    P(f"      price at the point factors: pds {pr3['point']['pds']:+.5f} mse {pr3['point']['mse']:+.5f} DE-three {pr3['point']['de3']:+.5f} "
      f"nmae {pr3['point']['nmae']:+.5f} = OVERALL {pr3['point']['total']:+.5f}  (conservative {pr3['conservative']['total']:+.5f}; "
      f"BR249 section 3: ~ -0.0014 mse vs +0.0002 pds)  -> {'PASS' if r3_ok else 'FAIL'}")
    if not r3_ok:
        fails.append("R3")
    out["retrodiction"] = dict(R1=dict(compare=c_r2, price=pr1, bands=R1_BANDS, band_ok=band_ok, ok=r1_ok,
                                       board=dict(pds=board_pds, mse=board_mse, overall=board_all),
                                       raw_mse_identity_check=dict(own_terms=raw_pred_own, rescaled_h1union_terms=raw_pred,
                                                                   board=RAW_MSE_CD4SVD5)),
                               R2=dict(compare=c_r2, ref=dict(d_pds=p246["d_pds"], rel_macro=p246["rel_macro"], rel_agg=ref_agg),
                                       max_abs_diff=d_r2, ok=r2_ok),
                               R3=dict(compare=c_r3, ref255=s255, max_abs_diff=d_r3, price=pr3, ok=r3_ok))
    if fails:
        P(f"\n*** RETRODICTION / KNOWN-ANSWER FAILURES: {fails} -- no cell is priced")
        (WORK / "price257.json").write_text(json.dumps(dict(out, failures=fails), indent=1, default=float))
        (WORK / "price257.log").write_text("\n".join(LINES) + "\n")
        raise SystemExit(1)

    # ================================================================= THE GRID, over cd4svd5
    cells = [f"k{k}_s{s:g}" for k in (3, 5, 8) for s in (0.25, 0.5, 0.75, 1.0)] + [f"k5_dir{j}" for j in range(1, 6)]
    table = {}
    P("\n" + "=" * 110)
    P("THE GRID -- every cell as an INCREMENT over cd4svd5 (k5_s1), three cancer folds; price at point / conservative factors")
    P("=" * 110)
    hdr = (f"{'cell':10s}{'E rm':>6s}{'frac x':>7s}{'d pds':>10s}{'folds':>6s}{'sig':>6s}{'agg%':>8s}{'macro%':>8s}{'dlv%':>7s}"
           f"{'dcosB':>9s}{'pds1.22':>9s}{'pds.85':>8s}{'mse':>9s}{'DE3':>8s}{'nmae':>8s}{'TOT pt':>9s}{'TOT cons':>9s}"
           f"{'@2.15':>8s}{'@1.36':>8s}{'kolf2 dp':>9s}{'kolf2 agg%':>11s}")
    P(hdr)
    for cell in cells + ["k10_s1"]:
        c = compare(rows, idn, cell, "k5_s1", CANCER)
        ck = compare(rows, idn, cell, "k5_s1", ("kolf2",))
        pr = price_both(c, COS_CD4SVD5, PER_CD4SVD5)
        cs = compare(rows, idn, cell, "STACK", CANCER)          # the cell's own reading vs the stack, for the exchange rate
        ctl = compare(rows, idn, f"randE_{cell}", "STACK", CANCER) if (CANCER[0], f"randE_{cell}") in rows else None
        table[cell] = dict(vs_cd4svd5=c, vs_stack=cs, kolf2_vs_cd4svd5=ck, price=pr,
                           control_vs_stack=ctl,
                           audit={f: rows[(f, cell)]["audit"] for f in CANCER},
                           shape_only=abs(c["delivered_ratio"] - 1) <= AMP_SHAPE_BOUND,
                           clears_conservative=pr["conservative"]["total"] > GATE_SHAPE,
                           clears_point=pr["point"]["total"] > GATE_SHAPE)
        pt_, cn_ = pr["point"], pr["conservative"]
        P(f"{cell:10s}{c['energy_removed']:6.3f}{c['fraction_ratio']:7.4f}{c['d_pds']:+10.6f}{c['pds_pos']:>4d}/3{c['sigma']:+6.1f}"
          f"{100*c['rel_agg']:+8.3f}{100*c['rel_macro']:+8.3f}{100*(c['delivered_ratio']-1):+7.3f}{pt_['dcos_board']:+9.5f}"
          f"{pt_['pds']:+9.5f}{cn_['pds']:+8.5f}{pt_['mse']:+9.5f}{pt_['de3']:+8.5f}{pt_['nmae']:+8.5f}{pt_['total']:+9.5f}"
          f"{cn_['total']:+9.5f}{pr['k_cos_2_15']['total']:+8.5f}{pr['k_cos_1_36']['total']:+8.5f}{ck['d_pds']:+9.5f}{100*ck['rel_agg']:+11.3f}")
    out["grid"] = table

    # ---- the per-direction reading: which direction carries the pds, which the cosine ------------
    P("\nPER-DIRECTION (k = 5, one direction removed alone, s = 1) -- vs the STACK (what each direction is worth) and vs cd4svd5")
    P(f"{'cell':10s}{'E share':>8s}{'E rm':>7s}{'d pds vs STACK':>15s}{'folds':>6s}{'agg% vs STACK':>14s}{'macro%':>8s}"
      f"{'pds per -1% agg':>16s}{'d nmae89':>9s}{'d prec89':>9s}{'| vs cd4svd5: d pds':>20s}{'agg%':>8s}{'TOT cons':>9s}")
    dir_energy = np.mean([rows[(f, "k5_s1")]["audit"]["per_direction_energy_frac_of_centred"] for f in CANCER], axis=0)
    for j in range(1, 6):
        cell = f"k5_dir{j}"
        cs = table[cell]["vs_stack"]
        c = table[cell]["vs_cd4svd5"]
        xr = cs["d_pds"] / (-100 * cs["rel_agg"]) if cs["rel_agg"] < 0 else float("inf")
        P(f"{cell:10s}{dir_energy[j-1]:8.4f}{cs['energy_removed']:7.4f}{cs['d_pds']:+15.6f}{cs['pds_pos']:>4d}/3{100*cs['rel_agg']:+14.3f}"
          f"{100*cs['rel_macro']:+8.3f}{xr:+16.5f}{cs['d_nmae']:+9.5f}{cs['d_prec89']:+9.5f}{c['d_pds']:+20.6f}{100*c['rel_agg']:+8.3f}"
          f"{table[cell]['price']['conservative']['total']:+9.5f}")
    sum_dirs = float(sum(table[f"k5_dir{j}"]["vs_stack"]["d_pds"] for j in range(1, 6)))
    P(f"  sum of the five single-direction d pds vs STACK {sum_dirs:+.6f} vs the joint k5_s1 {table['k5_s1']['vs_stack']['d_pds']:+.6f} "
      f"(the directions are not additive on pds: ratio {table['k5_s1']['vs_stack']['d_pds']/sum_dirs:.2f})")
    P(f"  per-direction share of the centred energy (mean of 3 folds): {[round(float(x), 4) for x in dir_energy]}")

    # ---- exchange rate vs the stack, every cell (246 section 3.5's statistic) ---------------------
    P("\nEXCHANGE RATE vs the STACK (246 section 3.5): d pds per -1 % aggregate cosine; cd4svd5 = 0.0266")
    for cell in cells + ["k10_s1"]:
        cs = table[cell]["vs_stack"]
        ctl = table[cell]["control_vs_stack"]
        xr = cs["d_pds"] / (-100 * cs["rel_agg"]) if cs["rel_agg"] < 0 else float("inf")
        cst = (f"control (random, {int(np.mean([rows[(f, f'randE_{cell}')]['audit']['rand_dim'] for f in CANCER]))} dims, "
               f"E {ctl['energy_removed']:.3f}): d pds {ctl['d_pds']:+.6f}, agg {100*ctl['rel_agg']:+.3f}%") if ctl else ""
        P(f"  {cell:10s} d pds {cs['d_pds']:+.6f} ({cs['pds_pos']}/3)  agg {100*cs['rel_agg']:+.3f}%  macro {100*cs['rel_macro']:+.3f}%  "
          f"rate {xr:+.4f}   {cst}")

    # ---- best cell ----------------------------------------------------------------------------
    cand = [c for c in cells if c != "k5_s1"]
    best_c = max(cand, key=lambda c: table[c]["price"]["conservative"]["total"])
    best_p = max(cand, key=lambda c: table[c]["price"]["point"]["total"])
    any_clear = [c for c in cand if table[c]["clears_conservative"] and table[c]["shape_only"]]
    any_clear_point = [c for c in cand if table[c]["clears_point"]]
    P("\n" + "=" * 110)
    P(f"BEST CELL (conservative total): {best_c}  {table[best_c]['price']['conservative']['total']:+.5f} overall "
      f"(point {table[best_c]['price']['point']['total']:+.5f}); best at the point factors: {best_p} "
      f"{table[best_p]['price']['point']['total']:+.5f} (conservative {table[best_p]['price']['conservative']['total']:+.5f})")
    P(f"cells clearing +{GATE_SHAPE} at the CONSERVATIVE factors (and shape-only): {any_clear or 'NONE'};  at the point factors: {any_clear_point or 'NONE'}")
    for tag, cell in (("best (conservative)", best_c), ("best (point)", best_p)):
        c = table[cell]["vs_cd4svd5"]
        pr = table[cell]["price"]
        ret_inc = retention(pt, cell, "k5_s1", CANCER)
        ret_stack = retention(pt, cell, "STACK", CANCER)
        ctl = table[cell]["control_vs_stack"]
        ctl_inc = compare(rows, idn, f"randE_{cell}", "k5_s1", CANCER)
        table[cell]["retention_vs_cd4svd5"] = ret_inc
        table[cell]["retention_vs_stack"] = ret_stack
        table[cell]["control_vs_cd4svd5"] = ctl_inc
        P(f"\n{tag}: {cell}  audit {json.dumps({k: v for k, v in rows[(CANCER[0], cell)]['audit'].items() if k in ('k', 's', 'dirs', 'energy_removed_frac', 'top_k_energy_frac_of_centred', 'fraction')})}")
        P(f"  increment over cd4svd5: d pds {c['d_pds']:+.6f} ({c['pds_pos']}/3 folds {['%+.5f' % x for x in c['per_fold_d_pds']]}, {c['sigma']:+.2f} sigma)  "
          f"agg {100*c['rel_agg']:+.3f}% (per fold {['%+.3f%%' % (100*x) for x in c['per_fold_rel_agg']]})  macro {100*c['rel_macro']:+.3f}%  "
          f"delivered {100*(c['delivered_ratio']-1):+.3f}% (shape-only {table[cell]['shape_only']})  d q_w {c['d_qw']:+.5f}  d nmae89 {c['d_nmae']:+.5f}  d prec89 {c['d_prec89']:+.5f}")
        P(f"  six members, point:        pds {pr['point']['pds']:+.5f}  mse {pr['point']['mse']:+.5f} (d cos_board {pr['point']['dcos_board']:+.6f}, d raw {pr['point']['draw_cos']+pr['point']['draw_amp']:+.5f})  "
          f"reach {pr['point']['reach']:+.5f} fid {pr['point']['fid']:+.5f} jac {pr['point']['jac']:+.5f}  nmae {pr['point']['nmae']:+.5f}  = {pr['point']['total']:+.5f}")
        P(f"  six members, conservative: pds {pr['conservative']['pds']:+.5f}  mse {pr['conservative']['mse']:+.5f}  "
          f"reach {pr['conservative']['reach']:+.5f} fid {pr['conservative']['fid']:+.5f} jac {pr['conservative']['jac']:+.5f}  nmae {pr['conservative']['nmae']:+.5f}  = {pr['conservative']['total']:+.5f}  "
          f"vs gate +{GATE_SHAPE}: {pr['conservative']['total']/GATE_SHAPE:+.2f}x")
        P(f"  sensitivity: k_cos 2.15 -> {pr['k_cos_2_15']['total']:+.5f};  k_cos 1.36 -> {pr['k_cos_1_36']['total']:+.5f};  exact-identity mse form -> {pr['exact_identity']['total']:+.5f}")
        P(f"  leave-top-3-out retention (246's statistic) of the increment's per-target d pds: {100*ret_inc['retained']:.1f}% "
          f"(per fold {['%.1f%%' % (100*x) for x in ret_inc['per_fold_retained']]});  of the cell vs STACK: {100*ret_stack['retained']:.1f}%")
        P(f"  energy-matched random control (246's rand, seed {rows[(CANCER[0], f'randE_{cell}')]['audit']['seed']}, "
          f"{[rows[(f, f'randE_{cell}')]['audit']['rand_dim'] for f in CANCER]} dims, E {ctl['energy_removed']:.4f} vs the cell's {table[cell]['vs_stack']['energy_removed']:.4f}): "
          f"vs STACK d pds {ctl['d_pds']:+.6f} ({ctl['pds_pos']}/3), agg {100*ctl['rel_agg']:+.3f}%  [the cell vs STACK: d pds {table[cell]['vs_stack']['d_pds']:+.6f}, agg {100*table[cell]['vs_stack']['rel_agg']:+.3f}%];  "
          f"control vs cd4svd5: d pds {ctl_inc['d_pds']:+.6f}, agg {100*ctl_inc['rel_agg']:+.3f}%")
    verdict = "CLEARS" if any_clear else "DOES NOT CLEAR"
    P(f"\nVERDICT: {verdict} -- {'a' if any_clear else 'no'} cell clears +{GATE_SHAPE} overall over cd4svd5 at the conservative factors"
      f"{' (' + ', '.join(any_clear) + ')' if any_clear else ''}; best conservative total {table[best_c]['price']['conservative']['total']:+.5f} ({best_c}), "
      f"best point total {table[best_p]['price']['point']['total']:+.5f} ({best_p})")
    out.update(best_conservative=best_c, best_point=best_p, clears_conservative=any_clear, clears_point=any_clear_point,
               verdict=verdict, peak_rss_gib=peak_gib())
    (WORK / "price257.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price257.log").write_text("\n".join(LINES) + "\n")
    P(f"\nwrote {WORK/'price257.json'} and price257.log; peak RSS {peak_gib():.2f} GiB")


if __name__ == "__main__":
    main()
