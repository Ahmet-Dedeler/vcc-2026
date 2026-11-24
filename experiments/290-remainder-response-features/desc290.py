"""290 step 4 — POST-HOC, labelled (PREREG §9): attribution, block ablations, the magnitude
partial, the best/worst-predicted genes, and whether 284's two tails differ in (d).

Everything here is written AFTER the preregistered numbers are read and cannot change a verdict.

    PYTHONDONTWRITEBYTECODE=1 uv run python desc290.py

Writes work/desc290.json, work/desc290.log.
"""
from __future__ import annotations

import json
import resource
import time

import numpy as np
import pandas as pd
from scipy import stats

import lib290 as L

T0 = time.time()
LOG: list[str] = []
N_PERM = 20

# R284 §6, the two described tails (20 genes each), copied verbatim from that document
TAIL_ABOVE = ["EIF3H", "ZC3H13", "OGFOD1", "NARS", "TANGO6", "RPAP1", "C14orf178", "TRRAP", "SMG7",
              "SUDS3", "SMG5", "RBMX2", "BOD1L1", "ZFC3H1", "PPRC1", "TCERG1", "CASC3", "CMTR1",
              "UPF2", "MED7"]
TAIL_BELOW = ["CDC37", "SAP30BP", "WDR74", "NCL", "TUBGCP6", "ATP5MF", "TBCC", "TAF1C", "STX5",
              "VCP", "DTYMK", "VPS37A", "NCAPG", "PSMB4", "NCAPD2", "BOP1", "POP5", "HAUS6",
              "GPSM2", "IRS2"]


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def rank_resid(a: np.ndarray, z: np.ndarray) -> np.ndarray:
    """Residual of rank(a) after linear regression on rank(z) (a rank-partial helper)."""
    m = np.isfinite(a) & np.isfinite(z)
    out = np.full(len(a), np.nan)
    ra, rz = stats.rankdata(a[m]), stats.rankdata(z[m])
    A = np.column_stack([np.ones(m.sum()), rz])
    beta, *_ = np.linalg.lstsq(A, ra, rcond=None)
    out[m] = ra - A @ beta
    return out


def main() -> None:
    h = L.check_prereg()
    P(f"prereg {h} MATCH")

    pred = pd.read_csv(L.WORK / "pred290.csv")
    feat = pd.read_csv(L.WORK / "feat290.csv")
    held = L.load_held()
    fitj = json.loads((L.WORK / "fit290.json").read_text())

    gidx = {g: i for i, g in enumerate(held["gene"].astype(str))}
    u = np.array([gidx[g] for g in pred["gene"].astype(str)])
    y = pred["remainder_half2"].to_numpy(float)
    foldu = pred["fold"].to_numpy(int)
    compu = pred["component"].to_numpy()
    cols = [c if c == "on_axis" else f"{c}_h1" for c in L.FEATS290]
    X = feat.loc[u, cols].to_numpy(float)
    X21 = held.loc[u, L.FEATS].to_numpy(float)
    out: dict = {"prereg_sha256": h}

    # ---- 1. fold-by-fold reading (the unbalanced-fold caveat, 284 §7.3) ------------------------
    per_fold = {}
    for k in range(5):
        m = foldu == k
        per_fold[int(k)] = dict(n=int(m.sum()), n_train=int((~m).sum()),
                                rho=L.spearman(y[m], pred["oof_response_ridge"].to_numpy()[m]),
                                rho_annot=L.spearman(y[m],
                                                     pred["oof_annotation_ridge"].to_numpy()[m]),
                                rho_r1=L.spearman(y[m], pred["remainder_half1"].to_numpy()[m]))
    P("per-fold out-of-fold Spearman (train size in brackets):")
    for k, v in per_fold.items():
        P(f"  fold {k}: n {v['n']:5d} (train {v['n_train']:5d})  response {v['rho']:+.4f}  "
          f"annotation {v['rho_annot']:+.4f}  remainder_half1 {v['rho_r1']:+.4f}")
    out["per_fold"] = per_fold

    # a balanced comparison: gene-level 5-fold that does NOT hold families (an upper bound)
    rng = np.random.default_rng(L.SEED + 101)
    rf = rng.integers(0, 5, size=len(y))
    p_rf, _ = L.oof_ridge(X, y, rf)
    out["random_folds_upper_bound"] = L.spearman(y, p_rf)
    P(f"POST-HOC family-blind random 5-fold (upper bound, families NOT held): "
      f"{out['random_folds_upper_bound']:+.4f}")

    # ---- 2. permutation importance on held genes ------------------------------------------------
    models = {}
    for k in range(5):
        te, tr = foldu == k, foldu != k
        st = L.Standardiser().fit(X[tr])
        from sklearn.linear_model import Ridge
        a = fitj["primary"]["ridge_info"][str(k)]["alpha"]
        models[k] = (st, Ridge(alpha=a).fit(st.transform(X[tr]), y[tr]))

    def predict_with(Xp: np.ndarray) -> np.ndarray:
        p = np.empty(len(y))
        for k in range(5):
            te = foldu == k
            st, m = models[k]
            p[te] = m.predict(st.transform(Xp[te]))
        return p

    base_rho = L.spearman(y, predict_with(X))
    rperm = np.random.default_rng(L.SEED + 111)
    imp = {}
    for j, c in enumerate(L.FEATS290):
        drops = []
        for _ in range(N_PERM):
            Xp = X.copy()
            for k in range(5):                       # permute WITHIN the held fold
                te = np.flatnonzero(foldu == k)
                Xp[te, j] = X[te[rperm.permutation(len(te))], j]
            drops.append(base_rho - L.spearman(y, predict_with(Xp)))
        imp[c] = dict(mean_drop=float(np.mean(drops)), sd=float(np.std(drops)))
    out["permutation_importance"] = imp
    out["base_rho_refit"] = base_rho
    P(f"\npermutation importance (drop in Spearman, {N_PERM} repeats; refit base {base_rho:+.4f}):")
    for c, v in sorted(imp.items(), key=lambda x: -x[1]["mean_drop"]):
        P(f"  {c:34s} {v['mean_drop']:+.4f} +- {v['sd']:.4f}")

    # ---- 3. block ablations ---------------------------------------------------------------------
    BLOCKS = {
        "(a) magnitude": ["logmag"],
        "(b) breadth": ["breadth"],
        "(c) genericness": ["gen_cos", "gen_abscos", "gen_pc5", "gen_pc20"],
        "(d) programmes": L.PE_COLS + ["pe_prolif_abs"],
        "(f) on-target": ["ontarget", "on_axis"],
    }
    abl = {}
    for name, members in BLOCKS.items():
        jj = [L.FEATS290.index(c) for c in members]
        only, _ = L.oof_ridge(X[:, jj], y, foldu)
        rest = [j for j in range(len(L.FEATS290)) if j not in jj]
        wo, _ = L.oof_ridge(X[:, rest], y, foldu)
        abl[name] = dict(alone=L.spearman(y, only), without=L.spearman(y, wo), n_cols=len(members))
        P(f"  {name:20s} alone {abl[name]['alone']:+.4f}   without {abl[name]['without']:+.4f}")
    out["block_ablation"] = abl

    # ---- 4. is it magnitude?  partial the magnitude and reliability out -------------------------
    logmag = feat.loc[u, "logmag_h1"].to_numpy(float)
    rbar = pred["r_bar"].to_numpy(float)
    p_resp = pred["oof_response_ridge"].to_numpy(float)
    r1 = pred["remainder_half1"].to_numpy(float)
    part = {}
    for zname, zv in (("logmag_h1", logmag), ("r_bar", rbar),
                      ("n_pairs_half2", pred["n_pairs_half2"].to_numpy(float))):
        yr, pr_ = rank_resid(y, zv), rank_resid(p_resp, zv)
        r1r = rank_resid(r1, zv)
        part[zname] = dict(response=L.spearman(yr, pr_), remainder_half1=L.spearman(yr, r1r))
        P(f"  partialling out {zname:16s}: response {part[zname]['response']:+.4f}   "
          f"remainder_half1 {part[zname]['remainder_half1']:+.4f}")
    # both at once
    def resid2(a, z1, z2):
        m = np.isfinite(a) & np.isfinite(z1) & np.isfinite(z2)
        o = np.full(len(a), np.nan)
        A = np.column_stack([np.ones(m.sum()), stats.rankdata(z1[m]), stats.rankdata(z2[m])])
        b, *_ = np.linalg.lstsq(A, stats.rankdata(a[m]), rcond=None)
        o[m] = stats.rankdata(a[m]) - A @ b
        return o
    part["logmag+r_bar"] = dict(
        response=L.spearman(resid2(y, logmag, rbar), resid2(p_resp, logmag, rbar)),
        remainder_half1=L.spearman(resid2(y, logmag, rbar), resid2(r1, logmag, rbar)))
    P(f"  partialling out logmag_h1 + r_bar : response {part['logmag+r_bar']['response']:+.4f}   "
      f"remainder_half1 {part['logmag+r_bar']['remainder_half1']:+.4f}")
    # does the response model carry the SAME signal as half 1's own remainder?
    C1 = pred["C_half1"].to_numpy(float)
    part["remainder_half1"] = dict(
        response=L.spearman(rank_resid(y, r1), rank_resid(p_resp, r1)),
        remainder_half1=float("nan"))
    part["response_prediction"] = dict(
        response=float("nan"),
        remainder_half1=L.spearman(rank_resid(y, p_resp), rank_resid(r1, p_resp)))
    P(f"  response model vs half1's own remainder: rho(pred, r1) "
      f"{L.spearman(p_resp, r1):+.4f}, rho(pred, C1) {L.spearman(p_resp, C1):+.4f}")
    P(f"  rho(y, pred | r1) {part['remainder_half1']['response']:+.4f}   "
      f"rho(y, r1 | pred) {part['response_prediction']['remainder_half1']:+.4f}")
    Xr1 = np.column_stack([X, r1])
    p_comb, _ = L.oof_ridge(Xr1, y, foldu)
    out["combined_features_plus_r1"] = L.spearman(y, p_comb)
    out["pred_vs_half1"] = dict(rho_pred_r1=L.spearman(p_resp, r1),
                                rho_pred_C1=L.spearman(p_resp, C1))
    P(f"  POST-HOC features + remainder_half1 as a 18th column: {out['combined_features_plus_r1']:+.4f}")
    out["partials"] = part

    # magnitude-matched strata
    q = np.nanquantile(logmag, [1 / 3, 2 / 3])
    strat = {}
    for lab, m in (("logmag low", logmag <= q[0]), ("logmag mid", (logmag > q[0]) & (logmag <= q[1])),
                   ("logmag high", logmag > q[1])):
        strat[lab] = dict(n=int(m.sum()), response=L.spearman(y[m], p_resp[m]),
                          remainder_half1=L.spearman(y[m], r1[m]))
        P(f"  {lab:12s} n {strat[lab]['n']:5d}  response {strat[lab]['response']:+.4f}  "
          f"remainder_half1 {strat[lab]['remainder_half1']:+.4f}")
    out["magnitude_strata"] = strat

    # ---- 5. best / worst predicted genes ---------------------------------------------------------
    err = y - p_resp
    o = np.argsort(err)
    genes = pred["gene"].to_numpy().astype(str)
    def rows(ix):
        return [dict(gene=genes[i], y=float(y[i]), pred=float(p_resp[i]),
                     remainder_half1=float(r1[i]), logmag_h1=float(logmag[i]),
                     pe_prolif_abs_h1=float(feat.loc[u[i], "pe_prolif_abs_h1"]),
                     gen_pc5_h1=float(feat.loc[u[i], "gen_pc5_h1"]),
                     n_pairs_half2=int(pred["n_pairs_half2"].to_numpy()[i])) for i in ix]
    out["most_overpredicted"] = rows(o[:20])
    out["most_underpredicted"] = rows(o[::-1][:20])
    ap = np.argsort(np.abs(err))
    out["best_predicted"] = rows(ap[:20])
    P("\nmost UNDER-predicted (actual remainder far above prediction): "
      + ", ".join(r["gene"] for r in out["most_underpredicted"]))
    P("most OVER-predicted: " + ", ".join(r["gene"] for r in out["most_overpredicted"]))
    P("best predicted (|error| smallest): " + ", ".join(r["gene"] for r in out["best_predicted"]))

    # ---- 6. 284's two tails vs the programme features ---------------------------------------------
    gi = {g: i for i, g in enumerate(genes)}
    ia = np.array([gi[g] for g in TAIL_ABOVE if g in gi])
    ib = np.array([gi[g] for g in TAIL_BELOW if g in gi])
    rest = np.setdiff1d(np.arange(len(genes)), np.concatenate([ia, ib]))
    P(f"\n284's tails present in the unit: above {len(ia)}/20, below {len(ib)}/20")
    tails = {}
    for c in (["logmag", "breadth", "gen_cos", "gen_abscos", "gen_pc5", "gen_pc20"]
              + L.PE_COLS + ["pe_prolif_abs", "pe_stress_abs"]):
        col = feat.loc[u, f"{c}_h1"].to_numpy(float)
        ma, mb, mr = np.nanmean(col[ia]), np.nanmean(col[ib]), np.nanmean(col[rest])
        sd = np.nanstd(col[rest])
        tails[c] = dict(above=float(ma), below=float(mb), rest=float(mr), sd_rest=float(sd),
                        z_above=float((ma - mr) / sd) if sd > 0 else float("nan"),
                        z_below=float((mb - mr) / sd) if sd > 0 else float("nan"),
                        mwu_p=float(stats.mannwhitneyu(col[ia][np.isfinite(col[ia])],
                                                       col[ib][np.isfinite(col[ib])]).pvalue))
        P(f"  {c:32s} above {ma:+.4f}  below {mb:+.4f}  rest {mr:+.4f}  "
          f"z {tails[c]['z_above']:+.2f}/{tails[c]['z_below']:+.2f}  MWU p {tails[c]['mwu_p']:.4f}")
    out["tails"] = tails
    out["tail_model_prediction"] = dict(
        above_pred=float(np.mean(p_resp[ia])), below_pred=float(np.mean(p_resp[ib])),
        rest_pred=float(np.mean(p_resp[rest])), above_y=float(np.mean(y[ia])),
        below_y=float(np.mean(y[ib])), rest_y=float(np.mean(y[rest])),
        rho_on_tails=L.spearman(np.concatenate([y[ia], y[ib]]),
                                np.concatenate([p_resp[ia], p_resp[ib]])))
    t = out["tail_model_prediction"]
    P(f"  model prediction on the tails: above {t['above_pred']:+.4f} (actual {t['above_y']:+.4f}), "
      f"below {t['below_pred']:+.4f} (actual {t['below_y']:+.4f}), rest {t['rest_pred']:+.4f} "
      f"(actual {t['rest_y']:+.4f}); rho over the 40 tail genes {t['rho_on_tails']:+.4f}")

    out["seconds"] = time.time() - T0
    out["peak_rss_gib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30
    (L.WORK / "desc290.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "desc290.log").write_text("\n".join(LOG))
    P(f"peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
