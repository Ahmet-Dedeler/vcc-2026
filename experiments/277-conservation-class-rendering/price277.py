"""277 -- the six-member price of every class-conditional rendering OVER `cd4svd5`, three COHERENT scenarios.

The recipe is 257's, re-used through `price272` (which itself KAs its `price()` against `price257.price`):

    pds       2.2138 x k_pds x d pds(3 cancer folds, vs cd4svd5) / 6
    mse       d cos_board = k_cos x rel cos_agg x 0.105470;  d raw = -1.996 x d cos_board (+ 243's amplitude leg
              through 197's identity at cd4svd5's own terms when the delivered amplitude moves > 0.2 %);
              overall = -1.051 x d raw / 6
    DE-three  rate x d cos_board, split reach / fid / jac by 218's card
    nmae      216/188's exact identity, -d nmae_wald89 / 0.60 / 6 (three cancer folds CONCATENATED) x multiplier

Bands.  V1-V3 / V5 / the controls re-shape OCCUPIED CD4 rows -- the deflation's own class, so k_cos
{1.36, 3.76} (re-dose / re-shape, BR249 section 2) and k_pds {0.85, 1.22}.  V4 changes the OPERATOR, so it
carries 272's band k_cos {0.94, 3.76}.  DE-three {0.0425, 0.24} coherent with k_cos; the identity x1.5 on a
loss / x0 on a gain (pessimistic), x1.0 (point).  POINT = k_pds 1.22, k_cos 3.76, DE-three 0.24 x 0.776,
identity x1.0 -- the set that retrodicts BR249.

KA-4 (nothing is priced unless it passes): `cd4svd5` vs `h1union` through the recipe at h1union's terms must
read BR249, +0.0043 +- 0.0008 overall, and equal 257's R1 to 1e-9.

    PYTHONDONTWRITEBYTECODE=1 uv run python price277.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib277 as Z  # noqa: E402

sys.path.insert(0, str(Z.E272))
import price272 as Q  # noqa: E402  (imports price257 as Q.P: constants, raw_for, scaled_terms, price)

P257 = Q.P
CANCER = Z.CANCER
FOLDS = Z.FOLDS
BASE = "V0_cd4svd5"
GATE = P257.GATE_SHAPE
K_PDS_BAND = (0.85, 1.22)
K_COS_RESHAPE = (1.36, 3.76)          # the deflation's own class (BR249 section 2): re-dose / re-shape
K_COS_OPERATOR = (0.94, 3.76)         # 272's band for an operator change
DE3_BAND = (0.0425, 0.24)
POINT = dict(k_pds=1.22, k_cos=3.76, de3=P257.DE4_PER_DCOS * P257.DE3_SHARE, x_loss=1.0, x_gain=1.0)
R_BANDS = dict(overall=(0.0043, 0.0008), pds=(0.0055, 0.0006), mse=(-0.0008, 0.0003))
OPERATOR_ARMS = ("V4_floor",)
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def load():
    rows, aud, arr, meta = {}, {}, {}, {}
    for h in FOLDS:
        d = json.loads((Z.WORK / f"render277-{h}.json").read_text())
        meta[h] = d["meta"]
        for a, r in d["arms"].items():
            rows[(h, a)] = r
        for a, r in d["audits"].items():
            aud[(h, a)] = r
        z = np.load(Z.WORK / f"render277-{h}.npz", allow_pickle=True)
        for k in z.files:
            arr[(h, k)] = z[k]
    return rows, aud, arr, meta


def nmae_concat(arr, arm, folds) -> float:
    U = np.concatenate([arr[(f, f"{arm}|U")] for f in folds])
    M = np.concatenate([arr[(f, f"{arm}|M")] for f in folds])
    D = np.concatenate([arr[(f, f"{arm}|D")] for f in folds])
    return float(np.mean(1 + (U - 2 * M) / D))


def compare(rows, aud, arr, arm, base, folds) -> dict:
    ra, rm, dp, dn_a, dn_b, fr = [], [], [], [], [], []
    for f in folds:
        a, b = rows[(f, arm)], rows[(f, base)]
        ra.append(a["cos_agg"] / b["cos_agg"] - 1.0)
        rm.append(a["a23_cos"] / b["a23_cos"] - 1.0)
        dp.append(a["pds"] - b["pds"])
        dn_a.append(a["delivered_norm"])
        dn_b.append(b["delivered_norm"])
        fa, fb = aud[(f, arm)].get("fraction"), aud[(f, base)].get("fraction")
        fr.append(fa / fb if (fa and fb) else float("nan"))
    nm_a, nm_b = nmae_concat(arr, arm, folds), nmae_concat(arr, base, folds)
    return dict(arm=arm, base=base, n=len(folds), rel_agg=float(np.mean(ra)), rel_macro=float(np.mean(rm)),
                d_pds=float(np.mean(dp)), per_fold_d_pds=[float(x) for x in dp],
                per_fold_rel_agg=[float(x) for x in ra], pds_pos=int(sum(x > 0 for x in dp)),
                sigma=float(np.mean(dp)) / P257.SIGMA_PDS,
                delivered_ratio=float(np.mean(dn_a) / np.mean(dn_b)),
                fraction_ratio=float(np.nanmean(fr)), nmae_arm=nm_a, nmae_base=nm_b, d_nmae=nm_a - nm_b,
                d_prec89=float(np.mean([rows[(f, arm)]["top89_precision"] - rows[(f, base)]["top89_precision"]
                                        for f in folds])))


def factors(which: str, c: dict, kcos_band) -> dict:
    gain_cos, gain_pds = c["rel_agg"] > 0, c["d_pds"] > 0
    if which == "POINT":
        return dict(POINT)
    if which == "PESSIMISTIC":
        return dict(k_pds=K_PDS_BAND[0] if gain_pds else K_PDS_BAND[1],
                    k_cos=kcos_band[0] if gain_cos else kcos_band[1],
                    de3=DE3_BAND[0] if gain_cos else DE3_BAND[1], x_loss=1.5, x_gain=0.0)
    if which == "OPTIMISTIC":
        return dict(k_pds=K_PDS_BAND[1] if gain_pds else K_PDS_BAND[0],
                    k_cos=kcos_band[1] if gain_cos else kcos_band[0],
                    de3=DE3_BAND[1] if gain_cos else DE3_BAND[0], x_loss=1.0, x_gain=1.0)
    raise ValueError(which)


def scenarios(c, kcos_band, cos_base=P257.COS_CD4SVD5, per_base=P257.PER_CD4SVD5) -> dict:
    return {w: Q.price(c, cos_base=cos_base, per_base=per_base, **factors(w, c, kcos_band))
            for w in ("PESSIMISTIC", "POINT", "OPTIMISTIC")}


def main() -> None:
    rows, aud, arr, meta = load()
    fails: list[str] = []
    ka_all = [(h, r) for h in FOLDS for r in meta[h]["ka"]]
    ka_ok = sum(1 for _, r in ka_all if r["ok"])
    P("=" * 132)
    P("277 -- rendering the identifiable LOW-CONSERVATION rows differently: the six-member price over cd4svd5")
    P("       (PLOLO 162, three cancer folds; kolf2 apart; the class = the bottom 30 % of 276's C_lolo per fold)")
    P("=" * 132)
    P(f"per-fold known-answer rows: {ka_ok}/{len(ka_all)} PASS  " +
      "; ".join(f"{h} {meta[h]['seconds']}s peak {meta[h]['peak_rss_gib']:.2f} GiB" for h in FOLDS))
    if ka_ok != len(ka_all):
        fails.append("render277 KA rows")
    for h in FOLDS:
        P(f"  {h:14s} {meta[h]['n_targets']} targets, {meta[h]['n_eligible']} with C_lolo, CLASS30 n = "
          f"{meta[h]['n_class']} (cut {meta[h]['class_cut']:+.5f}), V0 CD4 fraction {meta[h]['v0_fraction']:.6f}")

    # ------------------------------------------------------------------ KA-4: BR249 retrodicted
    c_r = compare(rows, aud, arr, BASE, "h1union", CANCER)
    s_r = Q.price(c_r, cos_base=P257.COS_H1UNION, per_base=P257.PER_H1UNION, **POINT)
    ref_r1 = json.loads((Z.E257 / "work/price257.json").read_text())["retrodiction"]["R1"]["price"]["point"]["total"]
    board = json.loads((Z.ROOT / "data/board-entries/Xv0ayvuoHpsKYWv2iu42.json").read_text())["entry"]["score_avg"]
    board_prev = json.loads((Z.ROOT / "data/board-entries/iUziqE1gewfcPKzYUDaB.json").read_text())["entry"]["score_avg"]
    band_ok = dict(pds=abs(s_r["pds"] - R_BANDS["pds"][0]) <= R_BANDS["pds"][1],
                   mse=abs(s_r["mse"] - R_BANDS["mse"][0]) <= R_BANDS["mse"][1],
                   overall=abs(s_r["total"] - R_BANDS["overall"][0]) <= R_BANDS["overall"][1],
                   vs257=abs(s_r["total"] - ref_r1) < 1e-9)
    ok_r = all(band_ok.values())
    P()
    P(f"KA-4 RETRODICTION  cd4svd5 vs h1union at h1union's terms: d pds {c_r['d_pds']:+.6f} (3/3 folds "
      f"{c_r['pds_pos']}/3), rel cos_agg {100*c_r['rel_agg']:+.3f} %, d nmae {c_r['d_nmae']:+.6f}")
    P(f"     -> pds {s_r['pds']:+.5f}  mse {s_r['mse']:+.5f}  DE-three {s_r['de3']:+.5f}  nmae "
      f"{s_r['nmae']:+.5f}  =  OVERALL {s_r['total']:+.5f}   (band +0.0043 +- 0.0008; 257 R1 {ref_r1:+.5f}, "
      f"|d| {abs(s_r['total']-ref_r1):.1e}; BOARD realised {board - board_prev:+.7f})  -> "
      f"{'PASS' if ok_r else 'FAIL'} {band_ok}")
    if not ok_r:
        fails.append("KA-4 retrodiction")

    # ------------------------------------------------------------------ the arm table
    arms = [a for a in rows if a[0] == CANCER[0]]
    names = [a[1] for a in arms if a[1] not in (BASE, "h1union", "BASE")]
    cells: dict[str, dict] = {}
    rng = np.random.default_rng(Z.SEED)
    for nm in names:
        band = K_COS_OPERATOR if nm in OPERATOR_ARMS else K_COS_RESHAPE
        if any((f, nm) not in rows for f in CANCER):
            continue
        c = compare(rows, aud, arr, nm, BASE, CANCER)
        sc = scenarios(c, band)
        diffs = [arr[(f, f"{nm}|pds")] - arr[(f, f"{BASE}|pds")] for f in CANCER]
        bf = Z.boot_folds(diffs, np.random.default_rng(Z.SEED), Z.B_BOOT)
        per_fold_boot = {f: Z.boot(d, np.random.default_rng(Z.SEED), Z.B_BOOT)
                         for f, d in zip(CANCER, diffs)}
        ret = Z.lt3o(np.concatenate(diffs))
        ret_fold = [Z.lt3o(d)["retention"] for d in diffs]
        # by class and by 262 group
        by_class, by_grp = {}, {}
        for f, d in zip(CANCER, diffs):
            m = arr[(f, "mask|CLASS30")].astype(bool)
            by_class.setdefault("class", []).append(float(d[m].mean()))
            by_class.setdefault("rest", []).append(float(d[~m].mean()))
            for g in ("e (base.truth <= 0)", "marginal (0, 0.05]", "saturated > 0.05"):
                gm = arr[(f, f"grp|{g}")].astype(bool)
                by_grp.setdefault(g, []).append(float(d[gm].mean()))
        kolf = (rows[("kolf2", nm)]["pds"] - rows[("kolf2", BASE)]["pds"]) if ("kolf2", nm) in rows else None
        cells[nm] = dict(compare=c, scen=sc, boot=bf, per_fold_boot=per_fold_boot, lt3o=ret,
                         lt3o_per_fold=ret_fold, by_class={k: float(np.mean(v)) for k, v in by_class.items()},
                         by_group={k: float(np.mean(v)) for k, v in by_grp.items()}, kolf2_d_pds=kolf,
                         kcos_band=band,
                         clears=bool(sc["PESSIMISTIC"]["total"] >= GATE and c["pds_pos"] == 3
                                     and bf["lo"] > 0 and ret["retention"] >= 0.5))

    P()
    P("ARM TABLE (3 cancer folds; d pds vs cd4svd5, its pooled paired target bootstrap, the class split, the "
      "three coherent six-member totals)")
    P(f"  {'arm':13s} {'d pds':>9s} {'[95% CI]':>21s} {'3/3':>4s} {'class':>9s} {'rest':>9s} "
      f"{'cos_agg':>8s} {'macro':>8s} {'dnmae':>9s} {'amp%':>6s} {'LT3O':>6s} {'PESS':>9s} {'POINT':>9s} "
      f"{'OPT':>9s} {'kolf2':>9s}")
    for nm in names:
        if nm not in cells:
            continue
        e = cells[nm]
        c, s = e["compare"], e["scen"]
        P(f"  {nm:13s} {c['d_pds']:+9.5f} [{e['boot']['lo']:+.5f},{e['boot']['hi']:+.5f}] "
          f"{c['pds_pos']:>3d}/3 {e['by_class']['class']:+9.5f} {e['by_class']['rest']:+9.5f} "
          f"{100*c['rel_agg']:+7.2f}% {100*c['rel_macro']:+7.2f}% {c['d_nmae']:+9.5f} "
          f"{100*(c['delivered_ratio']-1):+5.2f}% {e['lt3o']['retention']:+6.2f} "
          f"{s['PESSIMISTIC']['total']:+9.5f} {s['POINT']['total']:+9.5f} {s['OPTIMISTIC']['total']:+9.5f} "
          f"{(e['kolf2_d_pds'] if e['kolf2_d_pds'] is not None else float('nan')):+9.5f}")

    P()
    P("PER-FOLD d pds (cancer folds) and the 262-group split of the same difference")
    P(f"  {'arm':13s} {'k562':>9s} {'hct116':>9s} {'hek293t':>9s} | {'e<=0':>9s} {'marginal':>9s} {'saturated':>9s}")
    for nm in names:
        if nm not in cells:
            continue
        e = cells[nm]
        g = e["by_group"]
        P(f"  {nm:13s} " + " ".join(f"{x:+9.5f}" for x in e["compare"]["per_fold_d_pds"]) +
          " | " + f"{g['e (base.truth <= 0)']:+9.5f} {g['marginal (0, 0.05]']:+9.5f} {g['saturated > 0.05']:+9.5f}")

    # ------------------------------------------------------------------ the six legs of the headline cells
    P()
    P("SIX-MEMBER LEGS (pessimistic / point) for V1, V2, V3, V4")
    for nm in ("V1_undefl", "V2_mirror", "V3_noresid", "V4_floor"):
        if nm not in cells:
            continue
        for w in ("PESSIMISTIC", "POINT"):
            s = cells[nm]["scen"][w]
            P(f"  {nm:13s} {w:12s} k_pds {s['k_pds']:.2f} k_cos {s['k_cos']:.2f} de3 {s['de3_rate']:.4f} | "
              f"pds {s['pds']:+.5f} mse {s['mse']:+.5f} (dcos {s['dcos_board']:+.5f}) reach {s['reach']:+.5f} "
              f"fid {s['fid']:+.5f} jac {s['jac']:+.5f} nmae {s['nmae']:+.5f} = {s['total']:+.5f}")

    # ------------------------------------------------------------------ the controls and the verdict
    P()
    v1 = cells.get("V1_undefl")
    if v1 is not None:
        half = 0.5 * v1["compare"]["d_pds"]
        P(f"CONTROLS for V1 (d pds {v1['compare']['d_pds']:+.5f}; the rule asks every control BELOW half, "
          f"{half:+.5f})")
        for nm in ("C_rand1", "C_rand2", "C_rand3", "C_agree30", "C_agreeN", "C_top30", "V2_mirror"):
            if nm in cells:
                d = cells[nm]["compare"]["d_pds"]
                P(f"  {nm:13s} d pds {d:+.5f}   {'below half' if d < half else 'NOT below half'}")
    P()
    cleared = [nm for nm in cells if cells[nm]["clears"]]
    P(f"cells screened: {len(cells)}; clearing the pessimistic gate AND 3/3 folds AND interval > 0 AND "
      f"retention >= 50 %: {len(cleared)} {cleared}")
    P(f"pessimistic >= +{GATE}: {[nm for nm in cells if cells[nm]['scen']['PESSIMISTIC']['total'] >= GATE]}")
    P(f"point >= +{GATE}: {[nm for nm in cells if cells[nm]['scen']['POINT']['total'] >= GATE]}")

    out = dict(recipe=dict(POINT=POINT, K_PDS_BAND=K_PDS_BAND, K_COS_RESHAPE=K_COS_RESHAPE,
                           K_COS_OPERATOR=K_COS_OPERATOR, DE3_BAND=DE3_BAND, GATE=GATE,
                           PDS_SCALE=P257.PDS_SCALE, MSE_SCALE=P257.MSE_SCALE, MSE_SLOPE=P257.MSE_SLOPE,
                           COS_CD4SVD5=P257.COS_CD4SVD5, COS_H1UNION=P257.COS_H1UNION,
                           NMAE_SPAN=P257.NMAE_SPAN, SIGMA_PDS=P257.SIGMA_PDS),
               retrodiction=dict(compare=c_r, price=s_r, band_ok=band_ok, ref_257_R1=ref_r1,
                                 board_realised=board - board_prev),
               ka=dict(n=len(ka_all), ok=ka_ok), cells=cells, cleared=cleared, failures=fails)
    (Z.WORK / "price277.json").write_text(json.dumps(out, indent=1, default=float))
    (Z.WORK / "price277.log").write_text("\n".join(LINES) + "\n")
    if fails:
        P(f"*** {fails}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
