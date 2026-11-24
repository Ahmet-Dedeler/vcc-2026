"""293 -- pool the three splits into the challenge-regime construct ceiling (PREREGISTRATION.md sec 3.3, 4).

Each target is measured inside its OWN split's residual space (its own panel centring, its own scored
columns) -- the space the 2025 challenge would have scored it in.  The pooled statistic is the median over
the 17 two-construct targets, paired target bootstrap, 20,000 draws.

    PYTHONDONTWRITEBYTECODE=1 uv run python agg293.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib293 as X  # noqa: E402
from lib293 import L  # noqa: E402

TAG = sys.argv[1] if len(sys.argv) > 1 else ""
P = L.Logger(X.WORK / f"agg293{TAG}.log")
SPLITS = ("train", "validation", "test")

# 285 (F14) and 273 reference values, for the comparison table
REF285 = {"K562-GWPS (748 perts, n=25/50)": 0.09, "Jurkat": 0.17, "HepG2": 0.41,
          "CD4 (GWCD4i)": 0.04, "Broad-16 Cas9": 0.785}
REF273 = dict(construct_median=0.013512255772155394, lane_median=0.24291977562969289, n=5)


def boot_ratio(b: np.ndarray, a: np.ndarray, n: int = L.NBOOT, seed: int = L.SEED):
    """Paired target bootstrap of median(b)/median(a)."""
    m = np.isfinite(a) & np.isfinite(b)
    a, b = a[m], b[m]
    rng = np.random.default_rng(seed)
    ii = rng.integers(0, a.size, size=(n, a.size))
    r = np.median(b[ii], 1) / np.median(a[ii], 1)
    return dict(point=float(np.median(b) / np.median(a)), lo=float(np.percentile(r, 2.5)),
                hi=float(np.percentile(r, 97.5)), n=int(a.size),
                median_of_ratios=float(np.median(b / np.where(np.abs(a) > 1e-9, a, np.nan))),
                p_ge1=float((r >= 1.0).mean()))


def main() -> None:
    t0 = time.time()
    P(f"293 agg293.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    R = {s: json.loads((X.WORK / f"run293_{s}{TAG}.json").read_text()) for s in SPLITS}
    ka = json.loads((X.WORK / "ka293.json").read_text())
    P(f"  KA rows: {[(k['label'].split()[0], 'PASS' if k['ok'] else 'FAIL') for k in ka['ka']]}")
    assert not ka["failures"], "KA failure -- ceiling not read"

    KEYS = ("a_cos", "b_cos", "d_cos", "bfix_cos", "a2_cos", "x_cos", "a_pds", "b_pds", "d_pds",
            "a2_pds", "x_pds", "a_own", "b_own", "d_own")
    out: dict = {"ka": [dict(label=k["label"], ok=k["ok"]) for k in ka["ka"]], "rungs": {}, "per_target": {},
                 "reference": dict(f14_285=REF285, r273=REF273), "tag": TAG,
                 "thin": R["test"].get("thin", 1.0),
                 "depth": dict(native_median_umi=54256, h1test_median_umi=19965, vcc2026_ctrl_median_umi=20340)}
    rows_csv: list[dict] = []

    # ------------------------------------------------------------------ per-rung pooled statistics
    for n in X.RUNGS:
        pool: dict[str, list] = {k: [] for k in KEYS}
        tg, sp_, extra = [], [], []
        for s in SPLITS:
            r = R[s]["rungs"].get(str(n)) or R[s]["rungs"].get(n)
            if r is None:
                continue
            for t in r["cohort"]:
                v = r["per_target"][t]
                tg.append(t); sp_.append(s)
                for k in KEYS:
                    pool[k].append(np.nan if v.get(k) is None else v[k])
                c1, c2 = R[s]["multi"][t]
                big, small = ((c1, c2) if R[s]["cells_per_construct"][c1] >= R[s]["cells_per_construct"][c2]
                              else (c2, c1))
                extra.append(dict(target=t, split=s, c_big=big, c_small=small,
                                  cells_big=R[s]["cells_per_construct"][big],
                                  cells_small=R[s]["cells_per_construct"][small],
                                  kd_big=R[s]["kd_full"][big]["lfc"], kd_small=R[s]["kd_full"][small]["lfc"],
                                  kd_1=R[s]["kd_full"][c1]["lfc"], kd_2=R[s]["kd_full"][c2]["lfc"],
                                  mag=R[s]["mag_full"][t],
                                  panel_median=r["panel"]["cos"]["median"],
                                  panel_rank=r["panel"]["multi_rank"].get(t)))
        A = {k: np.asarray(pool[k], dtype=np.float64) for k in KEYS}
        row: dict = dict(n=n, n_targets=len(tg), targets=tg, splits=sp_)
        for k in KEYS:
            row[k] = dict(**L.stats(A[k]), boot=L.boot_median(A[k]))
        row["a_minus_b"] = L.boot_median_diff(A["a_cos"], A["b_cos"])
        row["a_minus_d"] = L.boot_median_diff(A["a_cos"], A["d_cos"])
        row["b_minus_x"] = L.boot_median_diff(A["b_cos"], A["x_cos"])   # (b) against the unrelated-pair null
        row["rho_guide"] = boot_ratio(A["b_cos"], A["a_cos"])
        row["rho_guide_own"] = boot_ratio(A["b_own"], A["a_own"])
        # disattenuated: (b) / sqrt((a) x (a2)) -- corrects for the SECOND construct's own reliability, so a
        # construct with no measurable response is no longer scored as a construct that disagrees.
        m2 = np.isfinite(A["a_cos"]) & np.isfinite(A["a2_cos"]) & np.isfinite(A["b_cos"])
        if m2.sum() >= 4:
            aa, a2, bb = A["a_cos"][m2], A["a2_cos"][m2], A["b_cos"][m2]
            rg = np.random.default_rng(L.SEED)
            jj = rg.integers(0, aa.size, size=(L.NBOOT, aa.size))

            def dis(x, y, z):
                den = np.sqrt(max(np.median(x), 0.0) * max(np.median(y), 0.0))
                return float(np.median(z) / den) if den > 1e-9 else np.nan

            dd = np.asarray([dis(aa[jj[i]], a2[jj[i]], bb[jj[i]]) for i in range(2000)])
            row["rho_dis"] = dict(point=dis(aa, a2, bb), lo=float(np.nanpercentile(dd, 2.5)),
                                  hi=float(np.nanpercentile(dd, 97.5)), n=int(m2.sum()),
                                  median_a2=float(np.median(a2)))
        # Spearman-Brown: (a) extrapolated to 400/400, and the implied between-construct cosine there
        k400 = 400.0 / n
        sb_a = np.asarray([L.spearman_brown(v, k400) for v in A["a_cos"]])
        sb_b = np.asarray([L.spearman_brown(v, k400) for v in A["b_cos"]])
        row["sb400"] = dict(a=L.boot_median(sb_a), b_direct=L.boot_median(sb_b),
                            b_implied=float(row["rho_guide"]["point"] * float(np.nanmedian(sb_a))),
                            b_implied_lo=float(row["rho_guide"]["lo"] * float(np.nanmedian(sb_a))),
                            b_implied_hi=float(row["rho_guide"]["hi"] * float(np.nanmedian(sb_a))))
        # covariates (PREREGISTRATION sec 3.3)
        ratio = A["b_cos"] / np.where(np.abs(A["a_cos"]) > 1e-9, A["a_cos"], np.nan)
        cov = {}
        for name, v in (("kd_weaker", np.asarray([min(e["kd_1"], e["kd_2"], key=abs) for e in extra])),
                        ("kd_mean", np.asarray([(e["kd_1"] + e["kd_2"]) / 2 for e in extra])),
                        ("kd_absdiff", np.asarray([abs(e["kd_1"] - e["kd_2"]) for e in extra])),
                        ("magnitude", np.asarray([e["mag"] for e in extra])),
                        ("min_cells", np.asarray([min(e["cells_big"], e["cells_small"]) for e in extra]))):
            for tgt_name, y in (("b_cos", A["b_cos"]), ("rho_ratio", ratio)):
                m = np.isfinite(v) & np.isfinite(y)
                if m.sum() >= 4:
                    r_ = spearmanr(v[m], y[m])
                    cov[f"{name}|{tgt_name}"] = dict(rho=float(r_.statistic), p=float(r_.pvalue), n=int(m.sum()))
        row["covariates"] = cov
        row["extra"] = extra
        out["rungs"][n] = row
        P(f"\n  n={n:3d}/side, {len(tg)} targets pooled over the three splits")
        P(f"    (a) same construct, cell halves  {L.fmt(row['a_cos']['boot'])}  (<=0 {row['a_cos']['share_le0']:.3f})")
        P(f"    (b) construct 1 vs construct 2   {L.fmt(row['b_cos']['boot'])}  (<=0 {row['b_cos']['share_le0']:.3f})")
        P(f"    (d) lane halves, one construct   {L.fmt(row['d_cos']['boot'])}  (<=0 {row['d_cos']['share_le0']:.3f})")
        P(f"    (a) - (b) {L.fmt(row['a_minus_b'])}   (a) - (d) {L.fmt(row['a_minus_d'])}")
        P(f"    (x) CROSS-TARGET null (C_big of t vs C_small of another target) {L.fmt(row['x_cos']['boot'])}"
          f"  (<=0 {row['x_cos']['share_le0']:.3f});  (b) - (x) {L.fmt(row['b_minus_x'])}")
        if "rho_dis" in row:
            P(f"    (a2) the OTHER construct's own cell halves  {L.fmt(row['a2_cos']['boot'])}  "
              f"(<=0 {row['a2_cos']['share_le0']:.3f}); disattenuated rho = med(b)/sqrt(med(a)med(a2)) "
              f"{row['rho_dis']['point']:+.4f} [{row['rho_dis']['lo']:+.4f}, {row['rho_dis']['hi']:+.4f}]")
        P(f"    rho_guide = med(b)/med(a) {row['rho_guide']['point']:+.4f} "
          f"[{row['rho_guide']['lo']:+.4f}, {row['rho_guide']['hi']:+.4f}]  "
          f"(median of per-target ratios {row['rho_guide']['median_of_ratios']:+.4f})")
        P(f"    scorer pds: (a) {row['a_pds']['median']:+.4f} / (b) {row['b_pds']['median']:+.4f} / "
          f"(d) {row['d_pds']['median']:+.4f}  [chance 0.5]")
        P(f"    SB to 400/400: (a) {L.fmt(row['sb400']['a'])}; between-construct implied "
          f"{row['sb400']['b_implied']:+.4f} [{row['sb400']['b_implied_lo']:+.4f}, {row['sb400']['b_implied_hi']:+.4f}]"
          f"; (b) extrapolated directly {L.fmt(row['sb400']['b_direct'])}")
        for k, v in sorted(cov.items()):
            if k.endswith("rho_ratio") or k.startswith(("kd_weaker", "magnitude")):
                P(f"      covariate {k:24s} Spearman {v['rho']:+.3f} (p {v['p']:.3f}, n {v['n']})")
        for i, t in enumerate(tg):
            rows_csv.append(dict(n=n, split=sp_[i], target=t, **{k: A[k][i] for k in KEYS},
                                 rho_ratio=ratio[i], c_big=extra[i]["c_big"], c_small=extra[i]["c_small"],
                                 cells_big=extra[i]["cells_big"], cells_small=extra[i]["cells_small"],
                                 kd_big=extra[i]["kd_big"], kd_small=extra[i]["kd_small"],
                                 magnitude=extra[i]["mag"], panel_median=extra[i]["panel_median"],
                                 panel_rank=extra[i]["panel_rank"]))

    # ------------------------------------------------------------------ per-split panel baselines
    out["panel"] = {}
    for s in SPLITS:
        out["panel"][s] = {}
        for n in X.RUNGS:
            r = R[s]["rungs"].get(str(n)) or R[s]["rungs"].get(n)
            if r:
                out["panel"][s][n] = dict(n_targets=r["panel"]["n_targets"], median=r["panel"]["cos"]["median"],
                                          mean=r["panel"]["cos"]["mean"], p10=r["panel"]["cos"]["p10"],
                                          p90=r["panel"]["cos"]["p90"], le0=r["panel"]["cos"]["share_le0"],
                                          boot=r["panel"]["cos"]["boot"], pds=r["panel"]["pds"]["median"])
        P(f"\n  panel baseline {s}: " + "  ".join(
            f"n={n} {v['median']:+.4f} ({v['n_targets']}t)" for n, v in out["panel"][s].items()))

    # ------------------------------------------------------------------ library-wide knockdown
    kd_lib = {}
    allk = []
    for s_ in SPLITS:
        v = np.asarray(sorted(R[s_]["kd_all"].values()))
        allk += list(v)
        kd_lib[s_] = dict(n=int(v.size), median=float(np.median(v)),
                          frac_abs_lt_02=float((np.abs(v) < 0.2).mean()),
                          frac_abs_lt_05=float((np.abs(v) < 0.5).mean()),
                          frac_abs_lt_10=float((np.abs(v) < 1.0).mean()),
                          frac_positive=float((v > 0).mean()))
        P(f"  knockdown, all {v.size} targets of {s_}: median lfc {np.median(v):+.3f}; |lfc| < 0.2 "
          f"{kd_lib[s_]['frac_abs_lt_02']:.1%}, < 0.5 {kd_lib[s_]['frac_abs_lt_05']:.1%}, < 1.0 "
          f"{kd_lib[s_]['frac_abs_lt_10']:.1%}; lfc > 0 {kd_lib[s_]['frac_positive']:.1%}")
    av = np.asarray(allk)
    kd_lib["pooled"] = dict(n=int(av.size), median=float(np.median(av)),
                            frac_abs_lt_02=float((np.abs(av) < 0.2).mean()),
                            frac_abs_lt_05=float((np.abs(av) < 0.5).mean()),
                            frac_abs_lt_10=float((np.abs(av) < 1.0).mean()))
    out["kd_library"] = kd_lib
    P(f"  POOLED {av.size} targets: median {np.median(av):+.3f}; |lfc| < 0.5 {kd_lib['pooled']['frac_abs_lt_05']:.1%}, "
      f"< 1.0 {kd_lib['pooled']['frac_abs_lt_10']:.1%}")

    # ------------------------------------------------------------------ the headline
    r200 = out["rungs"][200]; r100 = out["rungs"][100]
    head = r100 if r100["n_targets"] >= r200["n_targets"] else r200
    out["headline"] = dict(
        rung=head["n"], n_targets=head["n_targets"], rho_guide=head["rho_guide"],
        a=head["a_cos"]["boot"], b=head["b_cos"]["boot"], d=head["d_cos"]["boot"],
        sb400_a=head["sb400"]["a"], sb400_b_implied=head["sb400"]["b_implied"],
        pds_a=head["a_pds"]["median"], pds_b=head["b_pds"]["median"],
        construct_specific_fraction=1.0 - head["rho_guide"]["point"])
    P(f"\n  HEADLINE at n={head['n']}/side, {head['n_targets']} targets: rho_guide "
      f"{head['rho_guide']['point']:+.4f} [{head['rho_guide']['lo']:+.4f}, {head['rho_guide']['hi']:+.4f}]; "
      f"1 - rho_guide = {1 - head['rho_guide']['point']:.3f} of a 400-cell truth's reproducible variance is "
      f"construct-specific")

    with open(X.WORK / f"construct293{TAG}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_csv[0]))
        w.writeheader()
        w.writerows(rows_csv)
    out["seconds"] = round(time.time() - t0, 1)
    L.Logger.dump(P, out, X.WORK / f"agg293{TAG}.json")
    P(f"\ndone {out['seconds']}s; wrote construct293.csv ({len(rows_csv)} rows)")
    P.save()


if __name__ == "__main__":
    main()
