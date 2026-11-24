"""261 -- the six-member price of every union-deflation cell OVER `cd4svd5`, on three COHERENT scenarios
(READ-FIRST v2.8 section 3), the known-answer / retrodiction rows first.  Nothing is priced unless they pass.

    uv run python price261.py

Recipe = 257 section 1.3 (`price257.py` imported for its constants and its `price()`; this file's `price()` is the
same arithmetic with the DE-three rate and the identity multipliers as free arguments, asserted equal to 257's on
257's own factor sets):
    pds      2.2138 x k_pds x d pds(3 cancer folds, vs cd4svd5) / 6
    mse      d cos_board = k_cos x rel cos_agg x 0.105470;  d raw = -1.996 x d cos_board  (+ 243's amplitude leg,
             197's identity at cd4svd5's own terms, when the delivered amplitude moves > 0.2 %);  overall = -1.051 x d raw / 6
    DE-three rate x d cos_board, split reach / fid / jac by 218's card
    nmae     216/188's exact identity, -d nmae_wald89 / 0.60 / 6 (three cancer folds concatenated) x multiplier
Bands (the mandate): k_pds {0.85, 1.22}, k_cos {1.36, 3.76}, DE-three {0.0425, 0.24}, identity x1.5 on a loss / x0 on
a gain (pessimistic), x1.0 (point).  POINT = 257's point set (k_pds 1.22, k_cos 3.76, DE-three 0.24 x 0.776 = 0.186,
identity x1.0) -- the set that retrodicts BR249.  PESSIMISTIC / OPTIMISTIC: per leg, the band end that lowers /
raises the total for THAT cell's sign on that leg (one k_cos feeds `mse` and DE-three coherently).
CLEARED = PESSIMISTIC >= +0.0019 AND leave-top-3-out >= 50 % AND d pds > 0 on 3/3 folds at >= 2 LOO-sigma
(0.002988) AND the random and shuffled controls below half the arm's d pds.
"""
from __future__ import annotations

import json
import resource
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E246 = ROOT / "experiments/246-direction-hypotheses"
E249 = ROOT / "experiments/249-cd4-deflate"
E257 = ROOT / "experiments/257-soft-deflation-price"
sys.path.insert(0, str(E257))
import price257 as P  # noqa: E402  (constants, raw_for, scaled_terms, exact_slope, price -- read-only at import)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
FOLDS = CANCER + ("kolf2",)
K_GRID = (1, 2, 3, 5, 8)
K_PDS_BAND, K_COS_BAND, DE3_BAND = (0.85, 1.22), (1.36, 3.76), (0.0425, 0.24)
POINT = dict(k_pds=1.22, k_cos=3.76, de3=P.DE4_PER_DCOS * P.DE3_SHARE, x_loss=1.0, x_gain=1.0)
GATE = P.GATE_SHAPE
SIG2 = 2.0 * P.SIGMA_PDS
R2_BANDS = dict(overall=(0.0043, 0.0008), pds=(0.0055, 0.0006), mse=(-0.0008, 0.0003))
REAL = ["u_c0_g0.28"] + [f"u_k{k}_g0.28" for k in K_GRID] + [f"u_k{k}_g{g:g}" for g in (0.20, 0.35) for k in (3, 5)]
DOSE_REF = [f"u_k0_g{g:g}" for g in (0.20, 0.35)]
FIXED = ["u_k3_g0.28_fixed", "u_k5_g0.28_fixed"]
CONTROLS = {k: [f"randE_k{k}", f"shufP_k{k}", f"split_k{k}", f"split_k{k}b"] for k in K_GRID}
BASE = "cd4svd5"
LINES: list[str] = []


def Pr(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


# ----------------------------------------------------------------------------------- load
def load():
    rows, pt, idn, meta = {}, {}, {}, {}
    for held in FOLDS:
        d = json.loads((WORK / f"sweep261-{held}.json").read_text())
        meta[held] = d["meta"]
        for r in d["rows"]:
            rows[(held, r["arm"])] = r
        z = np.load(WORK / f"pertarget261-{held}.npz")
        for k in z.files:
            pt[(held, k)] = z[k]
        z = np.load(WORK / f"identity261-{held}.npz")
        for k in z.files:
            idn[(held, k)] = z[k]
    return rows, pt, idn, meta


def nmae_concat(idn, arm, folds) -> float:
    U = np.concatenate([idn[(f, f"{arm}|U")] for f in folds])
    M = np.concatenate([idn[(f, f"{arm}|M")] for f in folds])
    D = np.concatenate([idn[(f, f"{arm}|D")] for f in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(rows, idn, arm, base, folds) -> dict:
    """Everything the recipe reads for `arm` against `base`: mean over folds of per-fold ratios / deltas (257's
    `compare`, with the fraction ratio read from the union's audit where the cell is a union cell)."""
    ra, rm, dp, dq, dn_a, dn_b, fr, en = [], [], [], [], [], [], [], []
    for f in folds:
        a, b = rows[(f, arm)], rows[(f, base)]
        ra.append(a["cos_agg"] / b["cos_agg"] - 1.0)
        rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
        dp.append(a["pds"] - b["pds"])
        dq.append(a["a15_sgn_w"] - b["a15_sgn_w"])
        dn_a.append(a["delivered_norm"])
        dn_b.append(b["delivered_norm"])
        fr.append(a["audit"].get("fraction_ratio", a["audit"].get("fraction", np.nan) / b["audit"].get("fraction", np.nan)))
        en.append(a["audit"].get("energy_removed_frac", 0.0))
    nm_a, nm_b = nmae_concat(idn, arm, folds), nmae_concat(idn, base, folds)
    return dict(arm=arm, base=base, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)),
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp], per_fold_rel_agg=[float(x) for x in ra],
                per_fold_rel_macro=[float(x) for x in rm], pds_pos=int(sum(x > 0 for x in dp)),
                sigma=float(np.mean(dp)) / P.SIGMA_PDS, fold_se=float(np.std(dp, ddof=1) / np.sqrt(len(dp))) if len(dp) > 1 else float("nan"),
                d_qw=float(np.mean(dq)), delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)),
                fraction_ratio=float(np.nanmean(fr)), nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(f, arm)]["top89_precision"] - rows[(f, base)]["top89_precision"] for f in folds])),
                energy_removed=float(np.mean(en)))


def retention(pt, arm, base, folds, targets, n=3) -> dict:
    """Leave-top-n-out retention of the per-target d pds (each fold on its own ordering) and the top-n share."""
    full, loo, share, top = [], [], [], {}
    for f in folds:
        d = pt[(f, f"{arm}|pds")] - pt[(f, f"{base}|pds")]
        order = np.argsort(-d)
        keep = np.ones(d.size, bool)
        keep[order[:n]] = False
        full.append(d.mean())
        loo.append(d[keep].mean())
        share.append(float(d[order[:n]].sum() / d.sum()) if d.sum() != 0 else float("nan"))
        top[f] = [f"{targets[f][i]} {d[i]:+.3f}" for i in order[:n]]
    fm, lm = float(np.mean(full)), float(np.mean(loo))
    return dict(full=fm, leave_top_out=lm, retained=(lm / fm if fm else float("nan")), n=n,
                per_fold_retained=[float(l / f) if f else float("nan") for l, f in zip(loo, full)],
                share=float(np.nanmean(share)), per_fold_share=share, top_names=top)


# ---------------------------------------------------------------------------------- price
def price(c: dict, *, cos_base: float, per_base: dict, k_pds: float, k_cos: float, de3: float, x_loss: float,
          x_gain: float, amp_pct=None) -> dict:
    dcos = k_cos * c["rel_agg"] * cos_base
    draw_cos = P.MSE_SLOPE * dcos
    amp = (c["delivered_ratio"] - 1.0) if amp_pct is None else amp_pct / 100.0
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
    return dict(k_pds=k_pds, k_cos=k_cos, de3_rate=de3, x_loss=x_loss, x_gain=x_gain, dcos_board=dcos, draw_cos=draw_cos,
                draw_amp=draw_amp, amp_pct=100 * amp, amp_leg_applied=bool(abs(amp) > P.AMP_LEG_THRESHOLD),
                pds=leg_pds, mse=leg_mse, de3=leg_de3, **legs_de3, nmae=leg_nmae, nmae_identity_point=nm_point, total=total)


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


def scenarios(c: dict, cos_base=P.COS_CD4SVD5, per_base=P.PER_CD4SVD5, amp_pct=None) -> dict:
    return {w: price(c, cos_base=cos_base, per_base=per_base, amp_pct=amp_pct, **factors(w, c))
            for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC")}


def p257_sets(c: dict, cos_base, per_base) -> dict:
    """257's own two factor sets through this file's price() AND through price257.price() (code-path KA)."""
    mine = dict(point=price(c, cos_base=cos_base, per_base=per_base, k_pds=1.22, k_cos=3.76, de3=P.DE4_PER_DCOS * P.DE3_SHARE, x_loss=1.0, x_gain=1.0),
                conservative=price(c, cos_base=cos_base, per_base=per_base, k_pds=0.85, k_cos=3.76, de3=P.DE4_PER_DCOS, x_loss=1.5, x_gain=0.0))
    theirs = dict(point=P.price(c, cos_base=cos_base, per_base=per_base, k_pds=1.22, k_cos=3.76, de3_full=False, nmae_mode="point"),
                  conservative=P.price(c, cos_base=cos_base, per_base=per_base, k_pds=0.85, k_cos=3.76, de3_full=True, nmae_mode="conservative"))
    dmax = max(abs(mine[s][k] - theirs[s][k]) for s in mine for k in ("pds", "mse", "de3", "nmae", "total"))
    return dict(mine=mine, p257=theirs, max_abs_diff=dmax)


def main() -> None:
    rows, pt, idn, meta = load()
    targets = {f: meta[f]["targets"] for f in FOLDS}
    prod = json.loads((WORK / "prod261.json").read_text()) if (WORK / "prod261.json").exists() else None
    fails: list[str] = []
    out: dict = dict(recipe=dict(PDS_SCALE=P.PDS_SCALE, MSE_SCALE=P.MSE_SCALE, MSE_SLOPE=P.MSE_SLOPE, COS_CD4SVD5=P.COS_CD4SVD5,
                                 COS_H1UNION=P.COS_H1UNION, POINT=POINT, K_PDS_BAND=K_PDS_BAND, K_COS_BAND=K_COS_BAND, DE3_BAND=DE3_BAND,
                                 DE3_SPLIT=P.DE3_SPLIT, NMAE_SPAN=P.NMAE_SPAN, SIGMA_PDS=P.SIGMA_PDS, GATE=GATE, SIG2=SIG2,
                                 AMP_LEG_THRESHOLD=P.AMP_LEG_THRESHOLD, exact_slope_cd4svd5=P.exact_slope(P.PER_CD4SVD5)))
    Pr("=" * 120)
    Pr("261 -- deflate the H1-UNION residual's target-shared directions before its carrying-row sizing, on cd4svd5")
    Pr("       (PLOLO 162, three cancer folds, kolf2 apart; 246/257's construction; three COHERENT scenarios)")
    Pr("=" * 120)
    Pr(f"recipe: PDS_SCALE {P.PDS_SCALE}, MSE_SCALE {P.MSE_SCALE}, MSE_SLOPE {P.MSE_SLOPE} (exact at cd4svd5's terms {P.exact_slope(P.PER_CD4SVD5):.4f}), "
       f"cos_cd4svd5 {P.COS_CD4SVD5:.6f}; POINT k_pds {POINT['k_pds']} k_cos {POINT['k_cos']} DE-three {POINT['de3']:.4f} identity x1.0; "
       f"bands k_pds {K_PDS_BAND} k_cos {K_COS_BAND} DE-three {DE3_BAND} identity x1.5 loss / x0 gain")
    for h in FOLDS:
        m = meta[h]
        Pr(f"  {h:14s} {m['n_targets']} targets, union carries {m['union_carrying']} rows, {m['n_rows']} rows, {m['seconds']}s, peak {m['peak_rss_gib']:.2f} GiB, "
           f"known-answer {m['ka_pass']}/{m['ka_total']}, row KAs {sum(1 for r in rows.values() if r['fold'] == h and 'ka' in r and r['ka']['ok'])}/"
           f"{sum(1 for r in rows.values() if r['fold'] == h and 'ka' in r)}, failures {m['failures']}")
        if m["failures"] or m["ka_pass"] != m["ka_total"]:
            fails.append(f"sweep {h}")
    if prod is not None:
        Pr(f"  production frame: known-answer {prod['ka_pass']}/{prod['ka_total']}, cross-path cells {sum(c['ok'] for c in prod['cross_path'])}/{len(prod['cross_path'])}, "
           f"builder {prod['builder_sha256'][:16]}, {prod['seconds']}s, peak {prod['peak_rss_gib']:.2f} GiB")
        if prod["ka_pass"] != prod["ka_total"]:
            fails.append("prod261 KA")
    else:
        Pr("  production frame: prod261.json ABSENT (the amplitude leg is priced on PLOLO's proxy only)")

    # ---- identity known-answer rows vs six246.json (4 folds, as 257) ---------------------------------
    six = json.loads((E246 / "work/six246.json").read_text())["arms"]
    dn_a15 = float(np.mean([rows[(f, "a15")]["delivered_norm"] for f in FOLDS]))
    Pr("\nKNOWN-ANSWER -- 188's exact nmae identity / precision / delivered norm vs six246.json (4 folds concatenated):")
    n6 = 0
    for mine, theirs in (("a15", "a15 f=0.10 SHIPPED"), ("STACK", "STACK h1union g=0.28"), (BASE, "B2_cd4_svd5_resolved")):
        got = dict(nmae_wald89=nmae_concat(idn, mine, FOLDS),
                   top89_precision=float(np.mean([rows[(f, mine)]["top89_precision"] for f in FOLDS])),
                   delivered_norm_ratio=float(np.mean([rows[(f, mine)]["delivered_norm"] for f in FOLDS])) / dn_a15)
        for key, g in got.items():
            w = six[theirs][key]
            ok = abs(g - w) <= 5e-6
            n6 += ok
            Pr(f"  {'PASS' if ok else 'FAIL'}  {mine:8s} {key:22s} got {g:.8f} want {w:.8f} |d| {abs(g-w):.1e}")
            if not ok:
                fails.append(f"six246 KA {mine}/{key}")
    out["ka_six246"] = n6

    # ================================================================== R1 / R2 / R3
    Pr("\n" + "=" * 120)
    Pr("KNOWN-ANSWER AND RETRODICTION ROWS (must pass before any cell is priced)")
    Pr("=" * 120)
    # R1: k = 0 reproduces cd4svd5 on every fold and on the production frame
    r1_fold = all(x["ok"] for h in FOLDS for x in meta[h]["ka_rows"])
    r1_k0 = all(rows[(h, "u_k0_g0.28")]["ka"]["ok"] and rows[(h, "u_k0_g0.28")]["ka"]["max_abs_diff"] == 0.0 for h in FOLDS)
    r1_prod = prod is not None and all(x["ok"] for x in prod["known_answer"] if x["label"].startswith("KA-P"))
    Pr(f"R1  k = 0 through the resolved path == cd4svd5: {'0.0e+00 on 4/4 folds' if r1_k0 else 'FAIL'} "
       f"(max|d| {max(rows[(h, 'u_k0_g0.28')]['ka']['max_abs_diff'] for h in FOLDS):.1e}); cd4svd5 == 257 k5_s1 / 258 / 259 rows at 0.0: "
       f"{'PASS' if r1_fold else 'FAIL'}; production frame (249 effects-svd5.npy a28d5be2 re-composed, both code paths): "
       f"{'PASS' if r1_prod else ('ABSENT' if prod is None else 'FAIL')}  -> {'PASS' if (r1_k0 and r1_fold and r1_prod) else 'FAIL'}")
    if not (r1_k0 and r1_fold and r1_prod):
        fails.append("R1")
    # R3: the same deflate code on CD4 at k = 5 == 246 B2_cd4_svd5_resolved
    c_r3 = compare(rows, idn, BASE, "STACK", CANCER)
    p246 = json.loads((E246 / "work/price246.json").read_text())["B2_cd4_svd5_resolved"]
    ref_agg = P.AMP249["rel_agg_plolo"]
    d_r3 = max(abs(c_r3["d_pds"] - p246["d_pds"]), abs(c_r3["rel_macro"] - p246["rel_macro"]), abs(c_r3["rel_agg"] - ref_agg))
    r3_ok = d_r3 < P.TOL
    Pr(f"R3  lib246.deflate on CD4 at k = 5 (cd4svd5 vs STACK): d pds {c_r3['d_pds']:+.9f} (246 {p246['d_pds']:+.9f})  agg {100*c_r3['rel_agg']:+.6f}% "
       f"(249 {100*ref_agg:+.6f}%)  macro {100*c_r3['rel_macro']:+.6f}% (246 {100*p246['rel_macro']:+.6f}%)  max|d| {d_r3:.1e}  {'PASS' if r3_ok else 'FAIL'}; "
       f"pds {c_r3['pds_pos']}/3 folds {['%+.5f' % x for x in c_r3['per_fold_d_pds']]}, {c_r3['sigma']:.1f} sigma; energy removed {c_r3['energy_removed']:.4f}")
    if not r3_ok:
        fails.append("R3")
    # R2: the recipe retrodicts BR249 and 257's best cell
    s_r2 = p257_sets(c_r3, P.COS_H1UNION, P.PER_H1UNION)
    r2 = s_r2["mine"]["point"]
    band_ok = {k: abs(r2[k] - R2_BANDS[k][0]) <= R2_BANDS[k][1] for k in ("pds", "mse")}
    band_ok["overall"] = abs(r2["total"] - R2_BANDS["overall"][0]) <= R2_BANDS["overall"][1]
    board_pds = (P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_pds"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_pds"]) / 6
    board_mse = (P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_mse"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_mse"]) / 6
    board_all = P.BOARD["Xv0ayvuoHpsKYWv2iu42"]["score_avg"] - P.BOARD["iUziqE1gewfcPKzYUDaB"]["score_avg"]
    c_75 = compare(rows, idn, "k5_s0.75", BASE, CANCER)
    s_75 = p257_sets(c_75, P.COS_CD4SVD5, P.PER_CD4SVD5)
    p257 = json.loads((E257 / "work/price257.json").read_text())
    ref75 = p257["grid"]["k5_s0.75"]["price"]
    d_75 = max(abs(s_75["mine"]["point"]["total"] - ref75["point"]["total"]), abs(s_75["mine"]["conservative"]["total"] - ref75["conservative"]["total"]))
    ref_r1 = p257["retrodiction"]["R1"]["price"]["point"]["total"]
    r2_ok = all(band_ok.values()) and d_75 < 1e-9 and s_r2["max_abs_diff"] < 1e-12 and s_75["max_abs_diff"] < 1e-12 and abs(r2["total"] - ref_r1) < 1e-9
    Pr(f"R2  (a) BR249 through the recipe (cd4svd5 vs h1union, base cos {P.COS_H1UNION}, h1union's own terms, POINT): d cos_board {r2['dcos_board']:+.6f} "
       f"({100*r2['dcos_board']/P.COS_H1UNION:+.3f}%; board -1.735 %); d raw mse cosine {r2['draw_cos']:+.5f} + amplitude {r2['draw_amp']:+.5f} "
       f"(leg {'applied' if r2['amp_leg_applied'] else 'not applied'} at {r2['amp_pct']:+.3f}%) = {r2['draw_cos']+r2['draw_amp']:+.5f} (board +{P.RAW_MSE_CD4SVD5-P.RAW_MSE_H1UNION:.5f})")
    Pr(f"        pds {r2['pds']:+.5f} (band {R2_BANDS['pds'][0]:+.4f}+-{R2_BANDS['pds'][1]}; board {board_pds:+.5f})  mse {r2['mse']:+.5f} (band {R2_BANDS['mse'][0]:+.4f}+-{R2_BANDS['mse'][1]}; "
       f"board {board_mse:+.5f})  DE-three {r2['de3']:+.5f} (board -0.00034)  nmae {r2['nmae']:+.5f} (board -0.00005)  OVERALL {r2['total']:+.5f} "
       f"(band {R2_BANDS['overall'][0]:+.4f}+-{R2_BANDS['overall'][1]}; board {board_all:+.7f}; 257 R1 {ref_r1:+.5f}, |d| {abs(r2['total']-ref_r1):.1e})  {band_ok}")
    Pr(f"    (b) 257's best cell k5_s0.75 vs cd4svd5 on 257's factor sets: point {s_75['mine']['point']['total']:+.8f} (257 {ref75['point']['total']:+.8f}), "
       f"conservative {s_75['mine']['conservative']['total']:+.8f} (257 {ref75['conservative']['total']:+.8f}), max|d| {d_75:.1e}")
    Pr(f"    (c) this file's price() == price257.price() on both rows: max|d| {max(s_r2['max_abs_diff'], s_75['max_abs_diff']):.1e}  -> R2 {'PASS' if r2_ok else 'FAIL'}")
    if not r2_ok:
        fails.append("R2")
    out["retrodiction"] = dict(R1=dict(k0_identity=r1_k0, fold_rows=r1_fold, prod=r1_prod),
                               R2=dict(compare=c_r3, price=s_r2["mine"], bands=R2_BANDS, band_ok=band_ok, board=dict(pds=board_pds, mse=board_mse, overall=board_all),
                                       best257=dict(compare=c_75, price=s_75["mine"], ref=ref75, max_abs_diff=d_75), codepath_max_abs_diff=max(s_r2["max_abs_diff"], s_75["max_abs_diff"]), ok=r2_ok),
                               R3=dict(compare=c_r3, ref=dict(d_pds=p246["d_pds"], rel_macro=p246["rel_macro"], rel_agg=ref_agg), max_abs_diff=d_r3, ok=r3_ok))
    if fails:
        Pr(f"\n*** KNOWN-ANSWER / RETRODICTION FAILURES: {fails} -- no cell is priced")
        (WORK / "price261.json").write_text(json.dumps(dict(out, failures=fails), indent=1, default=float))
        (WORK / "price261.log").write_text("\n".join(LINES) + "\n")
        raise SystemExit(1)

    # ================================================================== the union's spectrum
    Pr("\n" + "=" * 120)
    Pr("THE UNION'S TARGET-SHARED SPECTRUM (centred carried-row Gram, share of the centred energy) vs CD4's")
    Pr("=" * 120)
    for h in FOLDS:
        s = meta[h]["union_spectrum"]
        Pr(f"  {h:14s} union {s['n_carried']:3d} rows: centred/raw {s['centred_energy_frac_of_raw']:.4f}; top-8 {[round(x, 4) for x in s['top_shares_of_centred'][:8]]}; "
           f"cum k=1/2/3/5/8 {[round(s['cum_shares_of_centred'][i-1], 3) for i in (1, 2, 3, 5, 8)]}")
    cd4_dir = np.mean([rows[(h, "k5_s0.75")]["audit"]["per_direction_energy_frac_of_centred"] for h in CANCER], axis=0)
    Pr(f"  CD4 (257, 3 cancer folds, 264-286 rows): per-direction shares {[round(float(x), 4) for x in cd4_dir]}, cum k=5 {float(cd4_dir.sum()):.3f}; "
       f"the union's rank-1 mean carries {100*np.mean([1-meta[h]['union_spectrum']['centred_energy_frac_of_raw'] for h in CANCER]):.2f} % of its raw energy "
       f"(CD4's: {100*(1-np.mean([rows[(h, 'k5_s0.75')]['audit']['centred_energy_frac_of_raw'] for h in CANCER])):.2f} %)")
    if prod is not None:
        s = prod["union_spectrum"]
        Pr(f"  production panel: union {s['n_carried']} rows, centred/raw {s['centred_energy_frac_of_raw']:.4f}, top-8 {[round(x, 4) for x in s['top_shares_of_centred'][:8]]}")

    # ================================================================== the grid
    amp_prod = (prod or {}).get("amplitude", {})
    table: dict[str, dict] = {}

    def cell(name):
        c = compare(rows, idn, name, BASE, CANCER)
        ck = compare(rows, idn, name, BASE, ("kolf2",))
        cs = compare(rows, idn, name, "STACK", CANCER)
        sc = scenarios(c)
        ap = amp_prod.get(name)
        sc_prod = scenarios(c, amp_pct=ap["calibrated_pct"]) if ap else None
        ret = retention(pt, name, BASE, CANCER, targets)
        d_nmae4 = nmae_concat(idn, name, FOLDS) - nmae_concat(idn, BASE, FOLDS)
        table[name] = dict(vs_cd4svd5=c, vs_stack=cs, kolf2=ck, scenarios=sc, scenarios_prod_amplitude=sc_prod, amplitude_prod=ap,
                           retention=ret, d_nmae_4fold=d_nmae4,
                           audit={f: {k: v for k, v in rows[(f, name)]["audit"].items() if k != "per_direction_energy_frac_of_centred"} for f in CANCER},
                           shape_only=abs(c["delivered_ratio"] - 1) <= P.AMP_SHAPE_BOUND,
                           cond_gate=bool(sc["PESSIMISTIC"]["total"] >= GATE), cond_retention=bool(ret["retained"] >= 0.5),
                           cond_pds=bool(c["pds_pos"] == 3 and c["d_pds"] >= SIG2))
        return table[name]

    for name in REAL + DOSE_REF + FIXED + [x for k in K_GRID for x in CONTROLS[k]]:
        cell(name)
    for k in K_GRID:
        arm = table[f"u_k{k}_g0.28"]["vs_cd4svd5"]["d_pds"]
        ctl = {n: table[n]["vs_cd4svd5"]["d_pds"] for n in CONTROLS[k]}
        table[f"u_k{k}_g0.28"]["controls_d_pds"] = ctl
        table[f"u_k{k}_g0.28"]["cond_controls"] = bool(arm > 0 and ctl[f"randE_k{k}"] < 0.5 * arm and ctl[f"shufP_k{k}"] < 0.5 * arm)
    for name in REAL:
        t = table[name]
        t["cleared"] = bool(t["cond_gate"] and t["cond_retention"] and t["cond_pds"] and t.get("cond_controls", False))

    Pr("\n" + "=" * 120)
    Pr("THE k x g GRID -- every cell as an INCREMENT over cd4svd5 (three cancer folds); PESSIMISTIC / POINT / OPTIMISTIC totals")
    Pr("=" * 120)
    hdr = (f"{'cell':18s}{'E rm':>6s}{'frac x':>7s}{'d pds':>10s}{'folds':>6s}{'sig':>6s}{'agg%':>8s}{'macro%':>8s}{'dlv%':>7s}{'prod%':>7s}"
           f"{'d nmae':>8s}{'pds':>9s}{'mse':>9s}{'DE3':>9s}{'nmae':>9s}{'PESS':>9s}{'POINT':>9s}{'OPT':>9s}{'ret3%':>7s}{'top3%':>7s}{'kolf2 dp':>9s}{'k2 agg%':>8s}")
    Pr(hdr)

    def line(name, tag=""):
        t = table[name]
        c, s = t["vs_cd4svd5"], t["scenarios"]
        pe, po, op = s["PESSIMISTIC"], s["POINT"], s["OPTIMISTIC"]
        ap = t["amplitude_prod"]
        Pr(f"{(name + tag):18s}{c['energy_removed']:6.3f}{c['fraction_ratio']:7.4f}{c['d_pds']:+10.6f}{c['pds_pos']:>4d}/3{c['sigma']:+6.1f}"
           f"{100*c['rel_agg']:+8.3f}{100*c['rel_macro']:+8.3f}{100*(c['delivered_ratio']-1):+7.3f}{(ap['calibrated_pct'] if ap else float('nan')):+7.3f}"
           f"{c['d_nmae']:+8.4f}{po['pds']:+9.5f}{po['mse']:+9.5f}{po['de3']:+9.5f}{po['nmae']:+9.5f}{pe['total']:+9.5f}{po['total']:+9.5f}{op['total']:+9.5f}"
           f"{100*t['retention']['retained']:+7.0f}{100*t['retention']['share']:+7.0f}{t['kolf2']['d_pds']:+9.5f}{100*t['kolf2']['rel_agg']:+8.3f}")

    for name in REAL[:6]:
        line(name)
    Pr("  -- the dose axis (k = 0 is 258's family, for reference; k = 3 / 5 re-dosed) --")
    for name in ("u_k0_g0.2", "u_k3_g0.2", "u_k5_g0.2", "u_k0_g0.35", "u_k3_g0.35", "u_k5_g0.35"):
        line(name)
    Pr("  -- the projector at the SHIPPED union fraction (sizing not re-solved) --")
    for name in FIXED:
        line(name)
    Pr("  (legs shown at POINT; 'dlv%' = PLOLO frame-A delivered-norm proxy, 'prod%' = production-panel A_c amplitude, calibrated; "
       "'ret3%' = leave-top-3-out retention of the per-target d pds, 'top3%' = the top-3 targets' share)")

    # ---- the controls -------------------------------------------------------------------------------
    Pr("\n" + "=" * 120)
    Pr("CONTROLS at every k (vs cd4svd5): energy-matched random subspace, shuffled projector, split-half projector (two seeds)")
    Pr("=" * 120)
    Pr(f"{'k':>3s}{'arm d pds':>11s}{'arm agg%':>9s}{'arm PESS':>10s} | {'randE dims':>10s}{'E':>6s}{'d pds':>10s}{'agg%':>8s}{'PESS':>9s} | {'shufP E':>8s}{'d pds':>10s}{'agg%':>8s}{'PESS':>9s} | "
       f"{'split E':>8s}{'d pds':>10s}{'agg%':>8s}{'PESS':>9s} | {'splitb d pds':>13s}{'agg%':>8s}{'PESS':>9s} | {'ctl<half':>9s}")
    for k in K_GRID:
        a = table[f"u_k{k}_g0.28"]
        r, sp, s1, s2 = (table[n] for n in CONTROLS[k])
        Pr(f"{k:>3d}{a['vs_cd4svd5']['d_pds']:+11.6f}{100*a['vs_cd4svd5']['rel_agg']:+9.3f}{a['scenarios']['PESSIMISTIC']['total']:+10.5f} | "
           f"{int(rows[(CANCER[0], f'randE_k{k}')]['audit']['rand_dim']):>10d}{r['vs_cd4svd5']['energy_removed']:6.3f}{r['vs_cd4svd5']['d_pds']:+10.6f}{100*r['vs_cd4svd5']['rel_agg']:+8.3f}{r['scenarios']['PESSIMISTIC']['total']:+9.5f} | "
           f"{sp['vs_cd4svd5']['energy_removed']:8.3f}{sp['vs_cd4svd5']['d_pds']:+10.6f}{100*sp['vs_cd4svd5']['rel_agg']:+8.3f}{sp['scenarios']['PESSIMISTIC']['total']:+9.5f} | "
           f"{s1['vs_cd4svd5']['energy_removed']:8.3f}{s1['vs_cd4svd5']['d_pds']:+10.6f}{100*s1['vs_cd4svd5']['rel_agg']:+8.3f}{s1['scenarios']['PESSIMISTIC']['total']:+9.5f} | "
           f"{s2['vs_cd4svd5']['d_pds']:+13.6f}{100*s2['vs_cd4svd5']['rel_agg']:+8.3f}{s2['scenarios']['PESSIMISTIC']['total']:+9.5f} | {str(a['cond_controls']):>9s}")

    # ---- exchange rate vs the stack ------------------------------------------------------------------
    Pr("\nEXCHANGE RATE vs the STACK (246 section 3.5): d pds per -1 % aggregate cosine; cd4svd5 = "
       f"{c_r3['d_pds'] / (-100 * c_r3['rel_agg']):.4f}; and vs cd4svd5 (the increment)")
    for name in REAL:
        cs, c = table[name]["vs_stack"], table[name]["vs_cd4svd5"]
        xr = cs["d_pds"] / (-100 * cs["rel_agg"]) if cs["rel_agg"] < 0 else float("inf")
        xi = c["d_pds"] / (-100 * c["rel_agg"]) if c["rel_agg"] < 0 else float("inf")
        Pr(f"  {name:14s} vs STACK d pds {cs['d_pds']:+.6f} ({cs['pds_pos']}/3) agg {100*cs['rel_agg']:+.3f}% rate {xr:+.4f} | increment d pds {c['d_pds']:+.6f} agg {100*c['rel_agg']:+.3f}% rate {xi:+.4f}")

    # ---- best cell --------------------------------------------------------------------------------------
    best = max(REAL, key=lambda n: table[n]["scenarios"]["PESSIMISTIC"]["total"])
    best_pt = max(REAL, key=lambda n: table[n]["scenarios"]["POINT"]["total"])
    best_pds = max(REAL, key=lambda n: table[n]["vs_cd4svd5"]["d_pds"])
    cleared = [n for n in REAL if table[n]["cleared"]]
    n_gate = sum(table[n]["cond_gate"] for n in REAL)
    n_point = sum(table[n]["scenarios"]["POINT"]["total"] >= GATE for n in REAL)
    n_opt = sum(table[n]["scenarios"]["OPTIMISTIC"]["total"] >= GATE for n in REAL)
    Pr("\n" + "=" * 120)
    Pr(f"THE SCREEN: {len(REAL)} real cells on cd4svd5.  PESSIMISTIC >= +{GATE}: {n_gate}/{len(REAL)};  POINT: {n_point}/{len(REAL)};  OPTIMISTIC: {n_opt}/{len(REAL)};  "
       f"ALL FOUR conditions: {len(cleared)}/{len(REAL)} {cleared or ''}")
    Pr(f"best PESSIMISTIC: {best} {table[best]['scenarios']['PESSIMISTIC']['total']:+.5f};  best POINT: {best_pt} {table[best_pt]['scenarios']['POINT']['total']:+.5f};  "
       f"largest d pds: {best_pds} {table[best_pds]['vs_cd4svd5']['d_pds']:+.6f}")
    for tag, name in (("best PESSIMISTIC", best), ("best POINT", best_pt), ("largest d pds", best_pds)):
        t = table[name]
        c, s = t["vs_cd4svd5"], t["scenarios"]
        Pr(f"\n{tag}: {name}  audit {json.dumps({k: (round(v, 6) if isinstance(v, float) else v) for k, v in t['audit'][CANCER[0]].items() if k in ('mode', 'k', 'g', 'n_carried', 'energy_removed_frac', 'top_k_energy_frac_of_centred', 'fraction', 'fraction_ratio', 'sizing')})}")
        Pr(f"  increment over cd4svd5: d pds {c['d_pds']:+.6f} ({c['pds_pos']}/3 folds {['%+.5f' % x for x in c['per_fold_d_pds']]}, {c['sigma']:+.2f} sigma, fold se {c['fold_se']:.5f})  "
           f"agg {100*c['rel_agg']:+.3f}% (per fold {['%+.3f%%' % (100*x) for x in c['per_fold_rel_agg']]})  macro {100*c['rel_macro']:+.3f}%  union fraction x{c['fraction_ratio']:.4f}  "
           f"delivered {100*(c['delivered_ratio']-1):+.3f}% PLOLO / {(t['amplitude_prod']['calibrated_pct'] if t['amplitude_prod'] else float('nan')):+.3f}% production (shape-only {t['shape_only']})  "
           f"d q_w {c['d_qw']:+.5f}  d nmae89 {c['d_nmae']:+.5f} (4-fold {t['d_nmae_4fold']:+.5f})  d prec89 {c['d_prec89']:+.5f}")
        for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC"):
            p = s[w]
            Pr(f"  {w:11s} k_pds {p['k_pds']:.2f} k_cos {p['k_cos']:.2f} DE3 {p['de3_rate']:.4f} id x{p['x_loss']:.1f}/{p['x_gain']:.1f}:  pds {p['pds']:+.5f}  mse {p['mse']:+.5f} "
               f"(d cos_board {p['dcos_board']:+.6f}, d raw {p['draw_cos']+p['draw_amp']:+.5f}, amp leg {'on' if p['amp_leg_applied'] else 'off'})  reach {p['reach']:+.5f} fid {p['fid']:+.5f} jac {p['jac']:+.5f}  "
               f"nmae {p['nmae']:+.5f}  = {p['total']:+.5f}  ({p['total']/GATE:+.2f}x gate)")
        if t["scenarios_prod_amplitude"]:
            sp = t["scenarios_prod_amplitude"]
            Pr(f"  with the production-panel amplitude ({t['amplitude_prod']['calibrated_pct']:+.3f}%): PESS {sp['PESSIMISTIC']['total']:+.5f} POINT {sp['POINT']['total']:+.5f} OPT {sp['OPTIMISTIC']['total']:+.5f}")
        r = t["retention"]
        Pr(f"  leave-top-3-out retention {100*r['retained']:.1f}% (per fold {['%.0f%%' % (100*x) for x in r['per_fold_retained']]}); top-3 share {100*r['share']:.1f}%; top-3 per fold: "
           + "; ".join(f"{h}: {', '.join(r['top_names'][h])}" for h in CANCER))
        if "controls_d_pds" in t:
            Pr(f"  controls d pds vs cd4svd5: {', '.join(f'{n} {v:+.6f}' for n, v in t['controls_d_pds'].items())}  (half the arm: {0.5*c['d_pds']:+.6f}; below-half {t['cond_controls']})")
        Pr(f"  conditions: pessimistic >= gate {t['cond_gate']}; retention >= 50% {t['cond_retention']}; pds 3/3 at >= 2 sigma {t['cond_pds']}; controls below half {t.get('cond_controls', 'n/a')}  -> cleared {t['cleared']}")

    # ---- sanity: best cell under the alternative k_cos reads -------------------------------------------
    verdict = "CLEARS" if cleared else "DOES NOT CLEAR"
    Pr(f"\nVERDICT: {verdict} -- {len(cleared)} of {len(REAL)} real cells meet all four conditions; best pessimistic {table[best]['scenarios']['PESSIMISTIC']['total']:+.5f} ({best}), "
       f"best point {table[best_pt]['scenarios']['POINT']['total']:+.5f} ({best_pt})")
    out.update(grid=table, best_pessimistic=best, best_point=best_pt, largest_d_pds=best_pds, cleared=cleared, verdict=verdict,
               n_gate_pessimistic=n_gate, n_gate_point=n_point, n_gate_optimistic=n_opt, peak_rss_gib=peak_gib())
    (WORK / "price261.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / "price261.log").write_text("\n".join(LINES) + "\n")
    Pr(f"\nwrote {WORK/'price261.json'} and price261.log; peak RSS {peak_gib():.2f} GiB")


if __name__ == "__main__":
    main()
