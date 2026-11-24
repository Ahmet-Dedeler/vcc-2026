"""277 -- the finer reads on top of `render277` / `flex277` / `price277`.

  1. every arm split by ITS OWN row set (touched / untouched), and per touched row;
  2. the class x 262-group crosstab for V1 (where the un-deflated residual gains and where it loses);
  3. the 0.02 knee: how many touched rows cross it, in each direction, and the per-target cosine move;
  4. WHERE THE DEFLATION ITSELF EARNS ITS `pds` -- `cd4svd5` - `h1union` by conservation decile (the
     scientific read the negative result leaves behind);
  5. the Flex-400 re-read next to the PLOLO one.

    PYTHONDONTWRITEBYTECODE=1 uv run python analyse277.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib277 as Z  # noqa: E402

CANCER = Z.CANCER
FOLDS = Z.FOLDS
BASE = "V0_cd4svd5"
# each arm -> the row set it actually edits (V4 edits the whitener on CLASS30)
# V2 keeps the class deflated (as cd4svd5 already does) and un-deflates EVERYTHING ELSE, so the rows it
# edits are the complement of the class -- marked with a leading "~".
ARM_MASK = {"V1_undefl": "CLASS30", "V1f_undefl": "CLASS30", "V2_mirror": "~CLASS30",
            "V2f_mirror": "~CLASS30", "V3_noresid": "CLASS30", "V3f_noresid": "CLASS30",
            "V5_20": "CLASS20", "V5_40": "CLASS40", "C_top30": "TOPCONS30",
            "C_agree30": "AGREE30", "C_agreeN": "AGREEN", "C_rand1": "RAND30_s1",
            "C_rand2": "RAND30_s2", "C_rand3": "RAND30_s3", "V4_floor": "CLASS30"}
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def main() -> None:
    arr, meta = {}, {}

    def getmask(f, mk):
        if mk.startswith("~"):
            return ~arr[(f, f"mask|{mk[1:]}")].astype(bool)
        return arr[(f, f"mask|{mk}")].astype(bool)

    for h in FOLDS:
        meta[h] = json.loads((Z.WORK / f"render277-{h}.json").read_text())["meta"]
        z = np.load(Z.WORK / f"render277-{h}.npz", allow_pickle=True)
        for k in z.files:
            arr[(h, k)] = z[k]
    out: dict = {}

    # ---------------------------------------------------------------- 1. per arm, its OWN row set
    P("=" * 122)
    P("1. EVERY ARM ON ITS OWN ROW SET (3 cancer folds; d pds vs cd4svd5, fold-mean; 'per row' = the summed")
    P("   difference over the rows the arm edits divided by their count -- the size-free reading)")
    P("=" * 122)
    P(f"  {'arm':13s} {'rows':>5s} {'fold-mean d':>12s} {'on edited':>11s} {'per edited row':>15s} "
      f"{'untouched':>11s}")
    own: dict[str, dict] = {}
    for nm, mk in ARM_MASK.items():
        if (CANCER[0], f"{nm}|pds") not in arr:
            continue
        rows_n, dfull, don, doff, per = [], [], [], [], []
        for f in CANCER:
            d = arr[(f, f"{nm}|pds")] - arr[(f, f"{BASE}|pds")]
            m = getmask(f, mk)
            rows_n.append(int(m.sum()))
            dfull.append(float(d.mean()))
            don.append(float(d[m].mean()))
            doff.append(float(d[~m].mean()))
            per.append(float(d[m].sum() / max(m.sum(), 1)))
        own[nm] = dict(mask=mk, n_rows=float(np.mean(rows_n)), d_full=float(np.mean(dfull)),
                       d_on=float(np.mean(don)), d_off=float(np.mean(doff)), per_row=float(np.mean(per)))
        e = own[nm]
        P(f"  {nm:13s} {e['n_rows']:5.0f} {e['d_full']:+12.5f} {e['d_on']:+11.5f} {e['per_row']:+15.5f} "
          f"{e['d_off']:+11.5f}")
    out["own_rows"] = own

    # ---------------------------------------------------------------- 2. class x group crosstab
    P()
    P("=" * 122)
    P("2. WHERE V1 GAINS AND LOSES -- 262's groups, ON THE EDITED ROWS ONLY (mean d pds and row count,")
    P("   pooled over the three cancer folds; the same for the mirror V2 and for a random row set)")
    P("=" * 122)
    groups = ("e (base.truth <= 0)", "marginal (0, 0.05]", "saturated > 0.05")
    tab: dict[str, dict] = {}
    for nm in ("V1_undefl", "V2_mirror", "V3_noresid", "C_rand1", "C_rand2", "C_rand3", "C_top30",
               "C_agree30", "V4_floor"):
        if (CANCER[0], f"{nm}|pds") not in arr:
            continue
        mk = ARM_MASK[nm]
        row = {}
        for g in groups:
            ds, ns = [], []
            for f in CANCER:
                d = arr[(f, f"{nm}|pds")] - arr[(f, f"{BASE}|pds")]
                m = getmask(f, mk) & arr[(f, f"grp|{g}")].astype(bool)
                if m.any():
                    ds.append(d[m])
                ns.append(int(m.sum()))
            row[g] = dict(mean=float(np.concatenate(ds).mean()) if ds else float("nan"),
                          n=int(np.sum(ns)),
                          total=float(np.sum([x.sum() for x in ds])) if ds else 0.0)
        tab[nm] = row
        P(f"  {nm:13s} " + "  ".join(f"{g.split(' ')[0]:9s} n={row[g]['n']:4d} mean {row[g]['mean']:+.5f}"
                                     for g in groups))
    out["class_by_group"] = tab

    # ---------------------------------------------------------------- 3. the 0.02 knee
    P()
    P("=" * 122)
    P("3. THE 0.02 KNEE on the edited rows (per-target own cosine to the truth, a23, three frames):")
    P("   how many cross it upward / downward, and the mean cosine move")
    P("=" * 122)
    knee: dict[str, dict] = {}
    for nm in ("V1_undefl", "V2_mirror", "V3_noresid", "C_rand1", "C_top30", "V4_floor"):
        if (CANCER[0], f"{nm}|cos") not in arr:
            continue
        mk = ARM_MASK[nm]
        up = dn = same = 0
        moves, dpds_up, dpds_dn = [], [], []
        for f in CANCER:
            c0 = arr[(f, f"{BASE}|cos")]
            c1 = arr[(f, f"{nm}|cos")]
            d = arr[(f, f"{nm}|pds")] - arr[(f, f"{BASE}|pds")]
            m = getmask(f, mk)
            a, b = c0[m] > 0.02, c1[m] > 0.02
            up += int((~a & b).sum())
            dn += int((a & ~b).sum())
            same += int((a == b).sum())
            moves.append(c1[m] - c0[m])
            dpds_up.append(d[m][~a & b])
            dpds_dn.append(d[m][a & ~b])
        mv = np.concatenate(moves)
        knee[nm] = dict(up=up, down=dn, unchanged=same, mean_cos_move=float(mv.mean()),
                        median_cos_move=float(np.median(mv)),
                        mean_dpds_crossers_up=float(np.concatenate(dpds_up).mean()) if up else float("nan"),
                        mean_dpds_crossers_down=float(np.concatenate(dpds_dn).mean()) if dn else float("nan"))
        e = knee[nm]
        P(f"  {nm:13s} crossings up {up:4d}  down {dn:4d}  unchanged {same:4d} | mean cosine move "
          f"{e['mean_cos_move']:+.5f} (median {e['median_cos_move']:+.5f}) | d pds of the up-crossers "
          f"{e['mean_dpds_crossers_up']:+.4f}, of the down-crossers {e['mean_dpds_crossers_down']:+.4f}")
    out["knee"] = knee

    # ---------------------------------------------------------------- 4. where the deflation earns its pds
    P()
    P("=" * 122)
    P("4. WHERE THE k = 5 DEFLATION ITSELF EARNS ITS `pds` -- (cd4svd5 - h1union) by conservation band")
    P("   (the same rows the map ranks; pooled over the three cancer folds, mean per (target, fold) row)")
    P("=" * 122)
    bands = [(0.0, 0.1, "bottom decile"), (0.1, 0.3, "10-30 %"), (0.3, 0.5, "30-50 %"),
             (0.5, 0.7, "50-70 %"), (0.7, 0.9, "70-90 %"), (0.9, 1.0, "top decile")]
    defl: dict[str, dict] = {}
    v1band: dict[str, dict] = {}
    for lo, hi, label in bands:
        ds, vs, ns, cs = [], [], [], []
        for f in CANCER:
            cl = arr[(f, "c_lolo")]
            ok = np.isfinite(cl)
            r = np.full(cl.shape, np.nan)
            order = np.argsort(cl[ok], kind="stable")
            rr = np.empty(order.size)
            rr[order] = (np.arange(order.size) + 0.5) / order.size
            r[ok] = rr
            m = ok & (r >= lo) & (r < hi if hi < 1.0 else r <= 1.0)
            dd = arr[(f, f"{BASE}|pds")] - arr[(f, "h1union|pds")]
            v1 = arr[(f, "V1_undefl|pds")] - arr[(f, f"{BASE}|pds")]
            ds.append(dd[m])
            vs.append(v1[m])
            cs.append(arr[(f, f"{BASE}|pds")][m])
            ns.append(int(m.sum()))
        dd_all, v1_all, c_all = np.concatenate(ds), np.concatenate(vs), np.concatenate(cs)
        defl[label] = dict(n=int(np.sum(ns)), mean_d_pds=float(dd_all.mean()),
                           share_of_total=float(dd_all.sum()), mean_pds=float(c_all.mean()))
        v1band[label] = dict(n=int(np.sum(ns)), mean_d_pds=float(v1_all.mean()))
    tot = sum(v["share_of_total"] for v in defl.values())
    P(f"  {'band':16s} {'n rows':>7s} {'mean cd4svd5 pds':>17s} {'deflation d pds':>16s} "
      f"{'share of its gain':>18s} {'V1 d pds':>10s}")
    for _, _, label in bands:
        e, v = defl[label], v1band[label]
        P(f"  {label:16s} {e['n']:7d} {e['mean_pds']:+17.5f} {e['mean_d_pds']:+16.5f} "
          f"{100*e['share_of_total']/tot:17.1f}% {v['mean_d_pds']:+10.5f}")
    out["deflation_by_band"] = defl
    out["v1_by_band"] = v1band

    # ---------------------------------------------------------------- 5. the Flex-400 re-read
    P()
    P("=" * 122)
    P("5. FLEX-400 (h1-test, union OFF, 273's instrument): the same arms against the union-OFF incumbent")
    P("=" * 122)
    fx = json.loads((Z.WORK / "flex277.json").read_text())
    P(f"  85 targets, CLASS30 n = {fx['meta']['n_class']} (cut {fx['meta']['class_cut']:+.5f}), "
      f"KA {'PASS' if not fx['meta']['failures'] else 'FAIL'}; the incumbent reads pds "
      f"{fx['arms']['INCUMBENT']['pds_a23']:.6f}, cos_agg {fx['arms']['INCUMBENT']['a23_cos_agg']:.5f}; "
      f"CD4 off reads pds {fx['arms']['cd4h1:off']['pds_a23']:.6f}, cos_agg "
      f"{fx['arms']['cd4h1:off']['a23_cos_agg']:.5f} (+{100*(fx['arms']['cd4h1:off']['a23_cos_agg']/fx['arms']['INCUMBENT']['a23_cos_agg']-1):.1f} % vs the incumbent)")
    P(f"  {'arm':13s} {'d pds':>9s} {'[95% CI]':>21s} {'on class':>10s} {'off class':>10s} "
      f"{'cos_agg':>8s} {'LT3O':>7s} {'PLOLO d pds':>12s}")
    plolo_d = {nm: json.loads((Z.WORK / "price277.json").read_text())["cells"][nm]["compare"]["d_pds"]
               for nm in json.loads((Z.WORK / "price277.json").read_text())["cells"]}
    for nm, r in fx["table"].items():
        P(f"  {nm:13s} {r['d_pds']:+9.5f} [{r['boot']['lo']:+.5f},{r['boot']['hi']:+.5f}] "
          f"{r['d_pds_class']:+10.5f} {r['d_pds_rest']:+10.5f} {100*r['rel_cos_agg']:+7.2f}% "
          f"{r['lt3o']['retention']:+7.2f} "
          f"{plolo_d.get(nm, float('nan')):+12.5f}")
    out["flex"] = fx["table"]

    # ---------------------------------------------------------------- 6. the map as a predictor of the deflation's gain
    from scipy.stats import rankdata, spearmanr

    def pool(per_fold_vals) -> np.ndarray:
        """276's convention: within-fold rank, concatenated over the three cancer folds; NaN kept."""
        acc = []
        for v in per_fold_vals:
            v = np.asarray(v, dtype=np.float64)
            r = np.full(v.size, np.nan)
            ok = np.isfinite(v)
            r[ok] = rankdata(v[ok]) / ok.sum()
            acc.append(r)
        return np.concatenate(acc)

    def boot_rho(x, y, b=2000, seed=Z.SEED):
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = np.asarray(x)[ok], np.asarray(y)[ok]
        rng = np.random.default_rng(seed)
        r0 = float(spearmanr(x, y).statistic)
        d = np.empty(b)
        for i in range(b):
            s = rng.integers(0, x.size, x.size)
            d[i] = spearmanr(x[s], y[s]).statistic
        return dict(rho=r0, ci=[float(v) for v in np.percentile(d, [2.5, 97.5])], n=int(x.size))

    def paired_drho(x1, x2, y, b=2000, seed=Z.SEED):
        ok = np.isfinite(x1) & np.isfinite(x2) & np.isfinite(y)
        x1, x2, y = np.asarray(x1)[ok], np.asarray(x2)[ok], np.asarray(y)[ok]
        d0 = spearmanr(x1, y).statistic - spearmanr(x2, y).statistic
        rng = np.random.default_rng(seed)
        d = np.empty(b)
        for i in range(b):
            s = rng.integers(0, x1.size, x1.size)
            d[i] = spearmanr(x1[s], y[s]).statistic - spearmanr(x2[s], y[s]).statistic
        return dict(delta=float(d0), ci=[float(v) for v in np.percentile(d, [2.5, 97.5])],
                    p_le0=float((d <= 0).mean()), n=int(x1.size))

    P()
    P("=" * 122)
    P("6. IS THE MAP A PREDICTOR OF WHERE THE DEFLATION PAYS?  Spearman with the per-target (cd4svd5 -")
    P("   h1union) `pds` difference, pooled by within-fold rank over the three cancer folds (276's")
    P("   convention), against the pipeline's own donor-donor agreement and the incumbent's own `pds`")
    P("=" * 122)
    y_defl = pool([arr[(f, f"{BASE}|pds")] - arr[(f, "h1union|pds")] for f in CANCER])
    y_v1 = pool([arr[(f, "V1_undefl|pds")] - arr[(f, f"{BASE}|pds")] for f in CANCER])
    x_cl = pool([arr[(f, "c_lolo")] for f in CANCER])
    x_dd = pool([arr[(f, "dd_agree")] for f in CANCER])
    x_bc = pool([arr[(f, "base_cos_truth")] for f in CANCER])
    x_p0 = pool([arr[(f, f"{BASE}|pds")] for f in CANCER])
    preds = {"C_lolo (276's conservation map)": x_cl, "donor-donor agreement (262)": x_dd,
             "base.truth cosine (262's mechanism)": x_bc, "cd4svd5's own per-target pds": x_p0}
    sp: dict[str, dict] = {}
    for label, x in preds.items():
        sp[label] = dict(vs_deflation=boot_rho(x, y_defl), vs_V1=boot_rho(x, y_v1))
        a, b_ = sp[label]["vs_deflation"], sp[label]["vs_V1"]
        P(f"  {label:38s} vs deflation d pds rho {a['rho']:+.3f} [{a['ci'][0]:+.3f}, {a['ci'][1]:+.3f}] "
          f"(n {a['n']}) | vs V1 d pds rho {b_['rho']:+.3f} [{b_['ci'][0]:+.3f}, {b_['ci'][1]:+.3f}]")
    # 6b: is the decile gradient of section 4 just HEADROOM?  partial Spearman controlling for the
    #     un-deflated arm's own per-target `pds` (the room a row has left to gain).
    def partial_rho(x, y, ctrl, b=2000, seed=Z.SEED):
        ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(ctrl)
        x, y, c = np.asarray(x)[ok], np.asarray(y)[ok], np.asarray(ctrl)[ok]

        def pr(xx, yy, cc):
            rx, ry, rc = rankdata(xx), rankdata(yy), rankdata(cc)
            rc = (rc - rc.mean()) / max(rc.std(), 1e-12)
            ex = rx - rx.mean() - np.dot(rx - rx.mean(), rc) / rc.size * rc
            ey = ry - ry.mean() - np.dot(ry - ry.mean(), rc) / rc.size * rc
            return float(np.corrcoef(ex, ey)[0, 1])
        r0 = pr(x, y, c)
        rng = np.random.default_rng(seed)
        d = np.empty(b)
        for i in range(b):
            s = rng.integers(0, x.size, x.size)
            d[i] = pr(x[s], y[s], c[s])
        return dict(rho=r0, ci=[float(v) for v in np.percentile(d, [2.5, 97.5])], n=int(x.size))

    x_h1 = pool([arr[(f, "h1union|pds")] for f in CANCER])
    sp["C_lolo vs deflation | h1union's own pds (headroom)"] = partial_rho(x_cl, y_defl, x_h1)
    sp["C_lolo vs h1union's own pds"] = boot_rho(x_cl, x_h1)
    a = sp["C_lolo vs deflation | h1union's own pds (headroom)"]
    h = sp["C_lolo vs h1union's own pds"]
    P(f"  partial rho(C_lolo, deflation d pds | h1union's own pds) {a['rho']:+.3f} "
      f"[{a['ci'][0]:+.3f}, {a['ci'][1]:+.3f}] (n {a['n']}); rho(C_lolo, h1union's own pds) "
      f"{h['rho']:+.3f} [{h['ci'][0]:+.3f}, {h['ci'][1]:+.3f}] -- the decile gradient of section 4 is "
      f"HEADROOM, not a deflation-specific signal")
    sp["paired C_lolo - donor-donor (vs deflation)"] = paired_drho(x_cl, x_dd, y_defl)
    pd_ = sp["paired C_lolo - donor-donor (vs deflation)"]
    P(f"  paired dRho (C_lolo - donor-donor agreement) on the deflation's gain: {pd_['delta']:+.3f} "
      f"[{pd_['ci'][0]:+.3f}, {pd_['ci'][1]:+.3f}], P(<=0) {pd_['p_le0']:.3f}, n {pd_['n']}")
    out["spearman"] = sp

    # ---------------------------------------------------------------- 7. the random-row-set null
    P()
    P("=" * 122)
    P("7. THE PROPER NULL FOR THE SELECTOR (null277.py): V1's construction on 24 size-matched RANDOM row")
    P("   sets per fold; the named row sets read as z and as a percentile of that null")
    P("=" * 122)
    nl = {f: json.loads((Z.WORK / f"null277-{f}.json").read_text()) for f in CANCER}
    P(f"  {'row set':11s} " + " ".join(f"{f[:11]:>28s}" for f in CANCER) + f" {'pooled z':>10s}")
    zsum: dict[str, float] = {}
    for k in ("CLASS30", "AGREE30", "AGREEN", "TOPCONS30"):
        zs, obs_m, null_m, var = [], [], [], []
        line = f"  {k:11s} "
        for f in CANCER:
            o = nl[f]["observed"][k]
            mu, sd = nl[f]["null"]["d_pds"]["mean"], nl[f]["null"]["d_pds"]["sd"]
            z = (o["d_pds"] - mu) / sd
            zs.append(z)
            obs_m.append(o["d_pds"])
            null_m.append(mu)
            var.append(sd ** 2)
            line += f"  {o['d_pds']:+.5f} z{z:+5.2f} ({100*o['null_frac_at_or_below_d_pds']:3.0f}%)"
        zp_ = (float(np.mean(obs_m)) - float(np.mean(null_m))) / (np.sqrt(np.sum(var)) / len(CANCER))
        zsum[k] = float(zp_)
        P(line + f" {zp_:+10.2f}")
    out["null_summary"] = dict(pooled_z=zsum,
                               per_fold={f: nl[f]["observed"] for f in CANCER},
                               null={f: nl[f]["null"] for f in CANCER})

    (Z.WORK / "analyse277.json").write_text(json.dumps(out, indent=1, default=float))
    (Z.WORK / "analyse277.log").write_text("\n".join(LINES) + "\n")


if __name__ == "__main__":
    main()
