"""269 -- the utility router: small models predicting each action's per-target utility, decisions with a margin,
three leakage-safe split settings reported apart, strong baselines, shuffled controls, paired target bootstraps.

    PYTHONDONTWRITEBYTECODE=1 uv run python model269.py

Inputs: work/rows_choose269.csv, work/rows_defl269.csv, work/table269.npz (build269.py).  Outputs: work/model269.json,
work/model269.log, work/tables269.md, work/decisions269.npz.
"""
from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import Ridge
from sklearn.model_selection import GroupKFold

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
CANCER = ["K562", "HCT116", "HEK293T"]
ALL = CANCER + ["KOLF2"]
DONORS = ["K562", "HCT116", "HEK293T"]
FAM = "family_km"
TAUS = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3]
TAU_REPORT = 0.02
B_BOOT = 20_000
SEED = 20260905
N_SHUF = 3
t0 = time.time()
LINES: list[str] = []
OUT: dict = dict(seed=SEED, b_boot=B_BOOT, taus=TAUS, tau_report=TAU_REPORT, settings={}, ka=[])
rng0 = np.random.default_rng(SEED)


def P(s=""):
    print(s, flush=True)
    LINES.append(s)


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


# ----------------------------------------------------------------------------------------------- data
dfc = pd.read_csv(WORK / "rows_choose269.csv")
dfd = pd.read_csv(WORK / "rows_defl269.csv")
npz = np.load(WORK / "table269.npz", allow_pickle=True)
all_targets = npz["all_targets"].tolist()
dm_profile = npz["dm_profile"]
prof = np.where(np.isfinite(dm_profile), dm_profile, np.nanmean(dm_profile, 1, keepdims=True))
prof_c = prof - prof.mean(1, keepdims=True)
prof_c /= np.maximum(np.linalg.norm(prof_c, axis=1, keepdims=True), 1e-12)
b_lines = npz["b_ctrl_lines"][:, npz["sh_axis"]]
line_names = npz["b_ctrl_line_names"].tolist()
for d in DONORS:
    dfc[f"d_is_{d}"] = (dfc.donor == d).astype(float)
dfc["s_dm_in_held_missing"] = dfc.s_dm_in_held.isna().astype(float)
dfd["s_dm_in_held_missing"] = dfd.s_dm_in_held.isna().astype(float)
dfc["r_dm_in_donor_missing"] = dfc.r_dm_in_donor.isna().astype(float)

SH = ["s_dd_raw", "s_dd_wht", "s_n_whitened", "s_cd4_norm", "s_cd4_energy_removed", "s_union_norm", "s_base_norm", "s_resid_norm",
      "s_resid_cos_base", "s_resid_share", "s_cd4_cos_base", "s_cd4_deflated_cos_base", "s_pred_nn", "s_pred_mean", "s_pred_norm",
      "s_cd4_de_log", "s_target_basal_held", "s_dm_mean", "s_dm_sd", "s_dm_in_held", "s_dm_in_held_missing", "s_collectri_out",
      "s_collectri_in", "s_dorothea_in", "s_is_tf", "s_basal_cc_pc1", "s_basal_cc_pc2", "s_basal_cc_pc3"]
REL = ["r_norm_raw", "r_norm_wht", "r_cd4_cos", "r_cd4_cos_w", "r_cd4d_cos", "r_cd4d_cos_w", "r_selfcons", "r_selfcons_w", "r_basal_sim",
       "r_target_basal_donor", "r_dm_in_donor", "r_dm_in_donor_missing"] + [f"d_is_{d}" for d in DONORS]
DMPC = [f"dm_pc{i+1}" for i in range(5)]
BPC = ["basal_tr_pc1", "basal_tr_pc2"]
FEAT_D = SH + DMPC + BPC
FEAT_C = FEAT_D + REL
P(f"features: deflation model {len(FEAT_D)}, choose model {len(FEAT_C)} (<= 50)")
assert len(FEAT_C) <= 50

ctx_n = {c: int(dfd[dfd.context == c].t_idx.max() + 1) for c in ALL}
ctx_targets = {c: dfd[dfd.context == c].sort_values("t_idx").target.tolist() for c in ALL}


def dm_pcs(train_gidx):
    pca = PCA(n_components=5, random_state=0).fit(prof_c[np.unique(train_gidx)])
    return pca.transform(prof_c)                                   # (300, 5)


def basal_pcs(train_ctx):
    idx = [line_names.index(c) for c in train_ctx]
    B = b_lines[idx]
    mu = B.mean(0)
    pca = PCA(n_components=2, random_state=0).fit(B - mu)
    Z = pca.transform(b_lines - mu)
    return {c: Z[line_names.index(c)] for c in ALL}


def design(df, feats, train_mask, pcs_dm, pcs_b):
    X = df[[f for f in feats if f not in DMPC + BPC]].to_numpy(float)
    X = np.column_stack([X, pcs_dm[df.g_idx.to_numpy()], np.stack([pcs_b[c] for c in df.context])])
    cols = [f for f in feats if f not in DMPC + BPC] + DMPC + BPC
    med = np.nanmedian(X[train_mask], 0)
    X = np.where(np.isfinite(X), X, med[None, :])
    mu, sd = X[train_mask].mean(0), X[train_mask].std(0) + 1e-9
    return (X - mu) / sd, cols


def fit_ridge(X, y, groups):
    alphas = [0.3, 1, 3, 10, 30, 100, 300, 1000]
    gkf = GroupKFold(n_splits=5)
    best, best_mse = None, np.inf
    for a in alphas:
        mse = 0.0
        for tr, va in gkf.split(X, y, groups):
            m = Ridge(alpha=a).fit(X[tr], y[tr])
            mse += float(((m.predict(X[va]) - y[va]) ** 2).mean())
        if mse < best_mse:
            best, best_mse = a, mse
    return Ridge(alpha=best).fit(X, y), best


def fit_gbm(X, y):
    return HistGradientBoostingRegressor(max_depth=3, max_iter=300, learning_rate=0.05, min_samples_leaf=20, l2_regularization=1.0,
                                         random_state=0).fit(X, y)


# ----------------------------------------------------------------------------------------- splits
def splits(setting):
    """yield (name, train_c mask, test_c mask, train_d mask, test_d mask, train contexts).  Test rows: cancer contexts only."""
    if setting == "A":
        for c in CANCER:
            trc = (dfc.context != c).to_numpy(); tec = (dfc.context == c).to_numpy()
            trd = (dfd.context != c).to_numpy(); ted = (dfd.context == c).to_numpy()
            yield f"held {c}", trc, tec, trd, ted, [x for x in ALL if x != c]
    elif setting == "B":
        for f in sorted(dfd[FAM].unique()):
            trc = (dfc[FAM] != f).to_numpy(); tec = ((dfc[FAM] == f) & dfc.context.isin(CANCER)).to_numpy()
            trd = (dfd[FAM] != f).to_numpy(); ted = ((dfd[FAM] == f) & dfd.context.isin(CANCER)).to_numpy()
            yield f"held family {f}", trc, tec, trd, ted, ALL
    elif setting == "C":
        for c in CANCER:
            for f in sorted(dfd[FAM].unique()):
                trc = ((dfc[FAM] != f) & (dfc.context != c)).to_numpy(); tec = ((dfc[FAM] == f) & (dfc.context == c)).to_numpy()
                trd = ((dfd[FAM] != f) & (dfd.context != c)).to_numpy(); ted = ((dfd[FAM] == f) & (dfd.context == c)).to_numpy()
                yield f"held {c} x family {f}", trc, tec, trd, ted, [x for x in ALL if x != c]


# --------------------------------------------------------------------------------- predictions per setting
def run_setting(setting, *, shuffle_features=None, shuffle_labels=None, feats_c=FEAT_C, feats_d=FEAT_D, tag=""):
    """Returns per-row out-of-sample predictions (choose rows and defl rows, cancer contexts), for ridge and gbm."""
    pred = {m: dict(c=np.full(len(dfc), np.nan), d=np.full(len(dfd), np.nan)) for m in ("ridge", "gbm")}
    alphas = []
    yc_all = dfc.u.to_numpy(float).copy(); yd_all = dfd.u_undefl.to_numpy(float).copy()
    rng = np.random.default_rng(shuffle_features if shuffle_features is not None else (shuffle_labels if shuffle_labels is not None else 0))
    Xc_base, cols_c = None, None
    for name, trc, tec, trd, ted, train_ctx in splits(setting):
        pcs_dm = dm_pcs(dfc.g_idx.to_numpy()[trc])
        pcs_b = basal_pcs(train_ctx)
        Xc, cols_c = design(dfc, feats_c, trc, pcs_dm, pcs_b)
        Xd, cols_d = design(dfd, feats_d, trd, pcs_dm, pcs_b)
        yc, yd = yc_all.copy(), yd_all.copy()
        if shuffle_features is not None:          # permute feature ROWS within (context, donor) / context: features carry no target information
            for c in ALL:
                for d in DONORS:
                    m = np.flatnonzero((dfc.context == c) & (dfc.donor == d))
                    if m.size:
                        Xc[m] = Xc[rng.permutation(m)]
                m = np.flatnonzero(dfd.context == c)
                Xd[m] = Xd[rng.permutation(m)]
        if shuffle_labels is not None:            # permute the TRAINING labels within (context, donor) / context
            for c in ALL:
                for d in DONORS:
                    m = np.flatnonzero(trc & (dfc.context == c).to_numpy() & (dfc.donor == d).to_numpy())
                    if m.size:
                        yc[m] = yc[rng.permutation(m)]
                m = np.flatnonzero(trd & (dfd.context == c).to_numpy())
                yd[m] = yd[rng.permutation(m)]
        gc = dfc.target.to_numpy()[trc]; gd = dfd.target.to_numpy()[trd]
        rc, ac = fit_ridge(Xc[trc], yc[trc], gc)
        rd, ad = fit_ridge(Xd[trd], yd[trd], gd)
        alphas.append((ac, ad))
        pred["ridge"]["c"][tec] = rc.predict(Xc[tec]); pred["ridge"]["d"][ted] = rd.predict(Xd[ted])
        gbc = fit_gbm(Xc[trc], yc[trc]); gbd = fit_gbm(Xd[trd], yd[trd])
        pred["gbm"]["c"][tec] = gbc.predict(Xc[tec]); pred["gbm"]["d"][ted] = gbd.predict(Xd[ted])
    return pred, alphas, cols_c


# ------------------------------------------------------------------------------------ decisions and scoring
def candidates(c):
    """Per target of context c: list of (action, label) from the choose rows and the undeflate row."""
    n = ctx_n[c]
    cand = [[] for _ in range(n)]
    sub = dfc[dfc.context == c]
    for t, d, u in zip(sub.t_idx, sub.donor, sub.u):
        cand[t].append((f"choose {d}", float(u)))
    subd = dfd[dfd.context == c]
    for t, u in zip(subd.t_idx, subd.u_undefl):
        cand[t].append(("undeflate", float(u)))
    return cand


CAND = {c: candidates(c) for c in ALL}
LAB_ORACLE = {c: np.array([max([0.0] + [u for _, u in CAND[c][t]]) for t in range(ctx_n[c])]) for c in ALL}
BEST_ACTION = {c: [("keep" if max([u for _, u in CAND[c][t]] + [0.0]) <= 0 else max(CAND[c][t], key=lambda x: x[1])[0]) for t in range(ctx_n[c])] for c in ALL}
ORACLE263 = {c: dfd[dfd.context == c].sort_values("t_idx").oracle263.to_numpy() for c in ALL}
for c in ALL:
    OUT["ka"].append(dict(label=f"[{c}] 263 realised oracle mean", got=float(ORACLE263[c].mean())))
P("labelled oracle per context: " + ", ".join(f"{c} {LAB_ORACLE[c].mean():+.4f} (263 realised {ORACLE263[c].mean():+.4f})" for c in ALL))


def decide(c, pred_c, pred_d, tau, *, choose_only=False):
    """pred_c: dict (t, donor) -> predicted utility; pred_d: dict t -> predicted.  Returns realised (n,), acted (n,), action list."""
    n = ctx_n[c]
    realised = np.zeros(n); acted = np.zeros(n, bool); action = ["keep"] * n
    for t in range(n):
        opts = []
        for a, u in CAND[c][t]:
            if a == "undeflate":
                if choose_only or t not in pred_d or not np.isfinite(pred_d[t]):
                    continue
                opts.append((pred_d[t], a, u))
            else:
                key = (t, a.split()[1])
                if key in pred_c and np.isfinite(pred_c[key]):
                    opts.append((pred_c[key], a, u))
        if not opts:
            continue
        p, a, u = max(opts, key=lambda x: x[0])
        if p > tau:
            realised[t] = u; acted[t] = True; action[t] = a
    return realised, acted, action


def boot(r, rng):
    n = r.size
    idx = rng.integers(0, n, size=(B_BOOT, n))
    means = r[idx].mean(1)
    return dict(point=float(r.mean()), lo=float(np.percentile(means, 2.5)), hi=float(np.percentile(means, 97.5)),
                sd=float(means.std()), p_le0=float((means <= 0).mean()))


def boot_strat(rs, rng):
    means = np.zeros(B_BOOT)
    for r in rs:
        idx = rng.integers(0, r.size, size=(B_BOOT, r.size))
        means += r[idx].mean(1) / len(rs)
    point = float(np.mean([r.mean() for r in rs]))
    return dict(point=point, lo=float(np.percentile(means, 2.5)), hi=float(np.percentile(means, 97.5)), sd=float(means.std()),
                p_le0=float((means <= 0).mean()))


def score_policy(realised_by_ctx, acted_by_ctx, action_by_ctx, rng):
    out = {}
    rs = []
    for c in CANCER:
        r = realised_by_ctx[c]; a = acted_by_ctx[c]
        b = boot(r, rng)
        acc = float(np.mean([action_by_ctx[c][t] == BEST_ACTION[c][t] for t in np.flatnonzero(a)])) if a.any() else np.nan
        out[c] = dict(**b, n=int(r.size), acted=int(a.sum()), acted_frac=float(a.mean()),
                      precision=float((r[a] > 0).mean()) if a.any() else np.nan,
                      captured_lab=float(r.mean() / LAB_ORACLE[c].mean()), captured_263=float(r.mean() / ORACLE263[c].mean()),
                      pick_acc=acc, sum_gain=float(r[r > 0].sum()), sum_loss=float(r[r < 0].sum()))
        rs.append(r)
    b = boot_strat(rs, rng)
    out["mean3"] = dict(**b, acted_frac=float(np.mean([acted_by_ctx[c].mean() for c in CANCER])),
                        captured_lab=float(b["point"] / np.mean([LAB_ORACLE[c].mean() for c in CANCER])),
                        captured_263=float(b["point"] / np.mean([ORACLE263[c].mean() for c in CANCER])),
                        folds_pos=int(sum(realised_by_ctx[c].mean() > 0 for c in CANCER)))
    return out


def policy_from_pred(pred, model, tau, *, choose_only=False):
    R, A, ACT = {}, {}, {}
    for c in CANCER:
        m = (dfc.context == c).to_numpy()
        pc = {(t, d): p for t, d, p in zip(dfc.t_idx[m], dfc.donor[m], pred[model]["c"][m])}
        md = (dfd.context == c).to_numpy()
        pd_ = {t: p for t, p in zip(dfd.t_idx[md], pred[model]["d"][md])}
        R[c], A[c], ACT[c] = decide(c, pc, pd_, tau, choose_only=choose_only)
    return R, A, ACT


def fixed_policy(action, tau=0.0):
    R, A, ACT = {}, {}, {}
    for c in CANCER:
        n = ctx_n[c]; r = np.zeros(n); a = np.zeros(n, bool); act = ["keep"] * n
        for t in range(n):
            for aa, u in CAND[c][t]:
                if aa == action:
                    r[t] = u; a[t] = True; act[t] = aa
        R[c], A[c], ACT[c] = r, a, act
    return R, A, ACT


def keep_policy():
    return ({c: np.zeros(ctx_n[c]) for c in CANCER}, {c: np.zeros(ctx_n[c], bool) for c in CANCER}, {c: ["keep"] * ctx_n[c] for c in CANCER})


def oracle_policy(kind):
    R, A, ACT = {}, {}, {}
    for c in CANCER:
        if kind == "labelled":
            R[c] = LAB_ORACLE[c].copy(); A[c] = LAB_ORACLE[c] > 0; ACT[c] = BEST_ACTION[c]
        else:
            R[c] = ORACLE263[c].copy(); A[c] = np.ones(ctx_n[c], bool); ACT[c] = ["oracle263"] * ctx_n[c]
    return R, A, ACT


def random_policy(acted_frac_by_ctx, rng):
    R, A, ACT = {}, {}, {}
    for c in CANCER:
        n = ctx_n[c]; r = np.zeros(n); a = np.zeros(n, bool); act = ["keep"] * n
        k = int(round(acted_frac_by_ctx[c] * n))
        for t in rng.choice(n, size=k, replace=False):
            if CAND[c][t]:
                aa, u = CAND[c][t][rng.integers(len(CAND[c][t]))]
                r[t] = u; a[t] = True; act[t] = aa
        R[c], A[c], ACT[c] = r, a, act
    return R, A, ACT


def target_only_policy(tau):
    """Setting A's target-level lookup: for the held context, each target's action is the one with the best MEAN label over
    the training contexts where that target has the action (066's same-target lookup, in utility form)."""
    R, A, ACT = {}, {}, {}
    for c in CANCER:
        n = ctx_n[c]; r = np.zeros(n); a = np.zeros(n, bool); act = ["keep"] * n
        other_c = dfc[dfc.context != c].groupby(["target", "donor"]).u.mean()
        other_d = dfd[dfd.context != c].groupby("target").u_undefl.mean()
        for t in range(n):
            tg = ctx_targets[c][t]
            opts = []
            for aa, u in CAND[c][t]:
                if aa == "undeflate":
                    if tg in other_d.index:
                        opts.append((float(other_d[tg]), aa, u))
                else:
                    key = (tg, aa.split()[1])
                    if key in other_c.index:
                        opts.append((float(other_c[key]), aa, u))
            if opts:
                p, aa, u = max(opts, key=lambda x: x[0])
                if p > tau:
                    r[t] = u; a[t] = True; act[t] = aa
        R[c], A[c], ACT[c] = r, a, act
    return R, A, ACT


def recipient_only_policy(tau):
    """The fixed action whose MEAN utility is best on the training contexts, applied everywhere in the held context."""
    R, A, ACT = {}, {}, {}
    for c in CANCER:
        tr_c = dfc[dfc.context != c]
        means = {f"choose {d}": float(tr_c[tr_c.donor == d].u.mean()) for d in DONORS}
        means["undeflate"] = float(dfd[dfd.context != c].u_undefl.mean())
        best = max(means, key=means.get)
        if means[best] > tau:
            r, a, act = fixed_policy(best)
            R[c], A[c], ACT[c] = r[c], a[c], act[c]
        else:
            R[c], A[c], ACT[c] = np.zeros(ctx_n[c]), np.zeros(ctx_n[c], bool), ["keep"] * ctx_n[c]
    return R, A, ACT


def pred_quality(pred, model):
    """Reference only: out-of-sample Spearman(pred, label) and AUC(label > 0) on the cancer-context rows."""
    from scipy.stats import spearmanr
    from sklearn.metrics import roc_auc_score
    q = {}
    for kind, pk, df, key in (("choose", "c", dfc, "u"), ("undeflate", "d", dfd, "u_undefl")):
        m = df.context.isin(CANCER).to_numpy() & np.isfinite(pred[model][pk])
        p, y = pred[model][pk][m], df[key].to_numpy()[m]
        q[kind] = dict(n=int(m.sum()), spearman=float(spearmanr(p, y).correlation), auc_pos=float(roc_auc_score(y > 0, p)),
                       auc_pos05=float(roc_auc_score(y > 0.05, p)) if (y > 0.05).any() else np.nan)
    return q


# ------------------------------------------------------------------------------------------ main loop
rng = np.random.default_rng(SEED)
base_rows = {}
P("\n== baselines (setting-independent; realised labels on the three cancer contexts) ==")
for nm, pol in [("keep (incumbent)", keep_policy()), ("oracle, labelled (ceiling)", oracle_policy("labelled")), ("oracle, 263 realised arm", oracle_policy("263"))] + \
        [(f"always {a}", fixed_policy(a)) for a in [f"choose {d}" for d in DONORS] + ["undeflate"]]:
    s = score_policy(*pol, rng)
    base_rows[nm] = s
    P(f"  {nm:32s} 3-ctx {s['mean3']['point']:+.4f} [{s['mean3']['lo']:+.4f}, {s['mean3']['hi']:+.4f}]  per ctx " +
      " / ".join(f"{s[c]['point']:+.4f}" for c in CANCER) + f"  acted {s['mean3']['acted_frac']:.2f}")
OUT["baselines"] = base_rows

DEC = {}
for setting in ["A", "B", "C"]:
    P(f"\n== setting {setting} ==")
    S = dict(rows={}, tau_curve={}, quality={}, alphas=None)
    pred, alphas, cols = run_setting(setting)
    S["alphas"] = alphas
    S["feature_columns"] = cols
    for model in ("ridge", "gbm"):
        S["quality"][model] = pred_quality(pred, model)
        P(f"  {model} out-of-sample quality: " + "; ".join(f"{k} n {v['n']} Spearman {v['spearman']:+.3f} AUC(u>0) {v['auc_pos']:.3f} AUC(u>0.05) {v['auc_pos05']:.3f}" for k, v in S["quality"][model].items()))
        curve = []
        for tau in TAUS:
            for co in (False, True):
                pol = policy_from_pred(pred, model, tau, choose_only=co)
                s = score_policy(*pol, rng)
                curve.append(dict(tau=tau, choose_only=co, **{k: s["mean3"][k] for k in ("point", "lo", "hi", "p_le0", "acted_frac", "captured_lab", "captured_263", "folds_pos")},
                                  per_ctx={c: s[c] for c in CANCER}))
                if tau == TAU_REPORT:
                    S["rows"][f"{model}{' (choose only)' if co else ''}"] = s
                    DEC[(setting, model, co)] = pol
        S["tau_curve"][model] = curve
        P(f"  {model} tau curve (all actions): " + " | ".join(f"τ {r['tau']:.3f}: {r['point']:+.4f} act {r['acted_frac']:.2f} +/3 {r['folds_pos']}" for r in curve if not r["choose_only"]))
        P(f"  {model} tau curve (choose only): " + " | ".join(f"τ {r['tau']:.3f}: {r['point']:+.4f} act {r['acted_frac']:.2f}" for r in curve if r["choose_only"]))
    # controls at the reporting tau: shuffled features, shuffled labels (N_SHUF seeds each), random action at the model's rate
    for ctrl in ("shuffled features", "shuffled labels"):
        for model in ("ridge", "gbm"):
            pts = []
            for sd in range(N_SHUF):
                kw = dict(shuffle_features=SEED + sd) if ctrl == "shuffled features" else dict(shuffle_labels=SEED + 100 + sd)
                pr, _, _ = run_setting(setting, **kw)
                pol = policy_from_pred(pr, model, TAU_REPORT)
                s = score_policy(*pol, rng)
                pts.append(s)
            S["rows"][f"{model}, {ctrl}"] = dict(seeds=[p["mean3"] for p in pts], per_ctx=[{c: p[c]["point"] for c in CANCER} for p in pts],
                                                mean3_point_mean=float(np.mean([p["mean3"]["point"] for p in pts])),
                                                mean3_point_range=[float(min(p["mean3"]["point"] for p in pts)), float(max(p["mean3"]["point"] for p in pts))],
                                                acted_frac_mean=float(np.mean([p["mean3"]["acted_frac"] for p in pts])))
    for model in ("ridge", "gbm"):
        af = {c: S["rows"][model][c]["acted_frac"] for c in CANCER}
        pts = [score_policy(*random_policy(af, np.random.default_rng(SEED + 200 + k)), rng) for k in range(N_SHUF)]
        S["rows"][f"random action at {model}'s rate"] = dict(seeds=[p["mean3"] for p in pts], mean3_point_mean=float(np.mean([p["mean3"]["point"] for p in pts])),
                                                             mean3_point_range=[float(min(p["mean3"]["point"] for p in pts)), float(max(p["mean3"]["point"] for p in pts))])
    if setting == "A":
        S["rows"]["target-only lookup (066 form)"] = score_policy(*target_only_policy(TAU_REPORT), rng)
        S["rows"]["target-only lookup, tau 0"] = score_policy(*target_only_policy(0.0), rng)
    S["rows"]["recipient-only (best fixed action on training)"] = score_policy(*recipient_only_policy(TAU_REPORT), rng)
    for nm, s in S["rows"].items():
        if "mean3" in s:
            P(f"  {nm:48s} 3-ctx {s['mean3']['point']:+.4f} [{s['mean3']['lo']:+.4f}, {s['mean3']['hi']:+.4f}] P(<=0) {s['mean3']['p_le0']:.3f}  per ctx " +
              " / ".join(f"{s[c]['point']:+.4f}" for c in CANCER) + f"  acted {s['mean3']['acted_frac']:.2f}  captured {s['mean3']['captured_lab']*100:+.1f} % (lab) {s['mean3']['captured_263']*100:+.1f} % (263)  +/3 {s['mean3']['folds_pos']}")
        else:
            P(f"  {nm:48s} 3-ctx mean {s['mean3_point_mean']:+.4f} range [{s['mean3_point_range'][0]:+.4f}, {s['mean3_point_range'][1]:+.4f}]")
    OUT["settings"][setting] = S
    P(f"  ({time.time() - t0:.0f} s, peak {peak_gib():.2f} GiB)")

# family sensitivity: setting C with the average-linkage families
P("\n== setting C with the average-linkage DepMap families (sensitivity) ==")
FAM = "family_dm"
pred, _, _ = run_setting("C")
sens = {}
for model in ("ridge", "gbm"):
    s = score_policy(*policy_from_pred(pred, model, TAU_REPORT), rng)
    sens[model] = s
    P(f"  {model:8s} 3-ctx {s['mean3']['point']:+.4f} [{s['mean3']['lo']:+.4f}, {s['mean3']['hi']:+.4f}]  per ctx " + " / ".join(f"{s[c]['point']:+.4f}" for c in CANCER) + f"  acted {s['mean3']['acted_frac']:.2f}")
OUT["setting_C_family_dm"] = sens
FAM = "family_km"

# ------------------------------------------------------------------------------------- failure analysis
P("\n== failure analysis ==")
FA = {}
from sklearn.metrics import roc_auc_score
# (i) held rows where the model acted and lost, setting C and A, both models, at the reporting tau: which features drove it
for setting in ("C", "A"):
    pred, _, cols = run_setting(setting)
    for model in ("ridge", "gbm"):
        R, A, ACT = policy_from_pred(pred, model, TAU_REPORT)
        recs = []
        for c in CANCER:
            for t in np.flatnonzero(A[c]):
                a = ACT[c][t]
                if a == "undeflate":
                    row = dfd[(dfd.context == c) & (dfd.t_idx == t)].iloc[0]
                    p = float(pred[model]["d"][row.name])
                else:
                    row = dfc[(dfc.context == c) & (dfc.t_idx == t) & (dfc.donor == a.split()[1])].iloc[0]
                    p = float(pred[model]["c"][row.name])
                recs.append(dict(context=c, t=t, target=ctx_targets[c][t], action=a, pred=p, realised=float(R[c][t]), group_e=bool(row.group_e), marginal=bool(row.marginal),
                                 **{k: float(row[k]) for k in SH if k in row.index}, **{k: float(row[k]) for k in REL if k in row.index}))
        fa = pd.DataFrame(recs)
        if len(fa) == 0:
            FA[f"{setting}/{model}"] = dict(n_acted=0)
            P(f"  [{setting}/{model}] acted on 0 targets at tau {TAU_REPORT}")
            continue
        lost, won = fa[fa.realised < 0], fa[fa.realised > 0]
        feats = [k for k in SH + REL if k in fa.columns and fa[k].notna().any()]
        sep = []
        for k in feats:
            x = fa[k].to_numpy(float); y = (fa.realised < 0).to_numpy()
            ok = np.isfinite(x)
            if y[ok].sum() >= 3 and (~y[ok]).sum() >= 3 and np.nanstd(x[ok]) > 0:
                sep.append((k, float(roc_auc_score(y[ok], x[ok])), float(np.nanmedian(lost[k])), float(np.nanmedian(won[k]))))
        sep.sort(key=lambda z: -abs(z[1] - 0.5))
        by_action = fa.groupby("action").realised.agg(["size", "mean", lambda v: float((v < 0).mean())]).rename(columns={"<lambda_0>": "loss_frac"})
        FA[f"{setting}/{model}"] = dict(n_acted=int(len(fa)), n_lost=int(len(lost)), n_won=int(len(won)), realised_mean=float(fa.realised.mean()),
                                        loss_sum=float(lost.realised.sum()), gain_sum=float(won.realised.sum()),
                                        lost_in_group_e=float(lost.group_e.mean()) if len(lost) else np.nan, won_in_group_e=float(won.group_e.mean()) if len(won) else np.nan,
                                        pred_lost_mean=float(lost.pred.mean()) if len(lost) else np.nan, pred_won_mean=float(won.pred.mean()) if len(won) else np.nan,
                                        by_action=by_action.reset_index().to_dict("records"),
                                        feature_auc_lost_vs_won=[dict(feature=k, auc_lost_higher=a, median_lost=ml, median_won=mw) for k, a, ml, mw in sep[:12]],
                                        spearman_pred_realised=float(fa[["pred", "realised"]].corr(method="spearman").iloc[0, 1]) if len(fa) > 3 else np.nan)
        P(f"  [{setting}/{model}] acted {len(fa)}: won {len(won)} (sum {won.realised.sum():+.2f}), lost {len(lost)} (sum {lost.realised.sum():+.2f}); "
          f"group (e) share lost {FA[f'{setting}/{model}']['lost_in_group_e']:.2f} vs won {FA[f'{setting}/{model}']['won_in_group_e']:.2f}; "
          f"Spearman(pred, realised) {FA[f'{setting}/{model}']['spearman_pred_realised']:+.3f}; by action " +
          ", ".join(f"{r['action']} n {r['size']} mean {r['mean']:+.3f} loss-frac {r['loss_frac']:.2f}" for r in FA[f"{setting}/{model}"]["by_action"]))
        P("     features separating lost from won (AUC lost-higher): " + ", ".join(f"{k} {a:.2f}" for k, a, _, _ in sep[:8]))

# (ii) the targets the oracle gains most on: does ANY feature separate them?  per-feature AUC vs a permutation reference
recs = []
for c in CANCER:
    top = LAB_ORACLE[c] >= np.quantile(LAB_ORACLE[c], 0.9)
    sub = dfd[dfd.context == c].set_index("t_idx")
    for t in range(ctx_n[c]):
        row = sub.loc[t]
        a = BEST_ACTION[c][t]
        rec = dict(context=c, t=t, top=bool(top[t]), gain=float(LAB_ORACLE[c][t]), group_e=bool(row.group_e), best=a, **{k: float(row[k]) for k in SH if k in row.index})
        if a.startswith("choose"):
            rc = dfc[(dfc.context == c) & (dfc.t_idx == t) & (dfc.donor == a.split()[1])]
            if len(rc):
                rec.update({k: float(rc.iloc[0][k]) for k in REL if k in rc.columns and not k.startswith("d_is")})
        recs.append(rec)
fo = pd.DataFrame(recs)
feats = [k for k in SH + REL if k in fo.columns and fo[k].notna().sum() > 30 and not k.startswith("d_is")]
aucs = {}
for k in feats:
    x = fo[k].to_numpy(float); y = fo.top.to_numpy(); ok = np.isfinite(x)
    if np.nanstd(x[ok]) > 0:
        aucs[k] = float(roc_auc_score(y[ok], x[ok]))
# permutation reference: the max |AUC - 0.5| over the same features with the top flag permuted within context
rngp = np.random.default_rng(SEED + 7)
maxdev = []
for _ in range(300):
    yp = fo.top.to_numpy().copy()
    for c in CANCER:
        m = np.flatnonzero((fo.context == c).to_numpy())
        yp[m] = yp[rngp.permutation(m)]
    dev = 0.0
    for k in aucs:
        x = fo[k].to_numpy(float); ok = np.isfinite(x)
        dev = max(dev, abs(roc_auc_score(yp[ok], x[ok]) - 0.5))
    maxdev.append(dev)
thr = float(np.percentile(maxdev, 95))
ranked = sorted(aucs.items(), key=lambda z: -abs(z[1] - 0.5))
FA["oracle_top_decile_separation"] = dict(n_top=int(fo.top.sum()), n=int(len(fo)), top_in_group_e=float(fo[fo.top].group_e.mean()),
                                          top_best_action=fo[fo.top].best.value_counts().to_dict(),
                                          feature_auc=[dict(feature=k, auc_top_higher=a) for k, a in ranked],
                                          perm_max_abs_dev_p95=thr, n_features_beyond_perm=int(sum(abs(a - 0.5) > thr for _, a in ranked)))
P(f"  oracle top decile (n {int(fo.top.sum())} of {len(fo)}; {fo[fo.top].group_e.mean()*100:.0f} % in group (e); best action {fo[fo.top].best.value_counts().to_dict()}):")
P("     per-feature AUC (top-decile higher): " + ", ".join(f"{k} {a:.2f}" for k, a in ranked[:12]))
P(f"     permutation reference for the max |AUC-0.5| over {len(aucs)} features: 95th pct {thr:.3f} -> {FA['oracle_top_decile_separation']['n_features_beyond_perm']} features beyond it")
# same for group (e) membership itself (262's anti-aligned group) as a second definition of "the targets that matter"
aucs_e = {}
for k in feats:
    x = fo[k].to_numpy(float); y = fo.group_e.to_numpy(); ok = np.isfinite(x)
    if np.nanstd(x[ok]) > 0:
        aucs_e[k] = float(roc_auc_score(y[ok], x[ok]))
ranked_e = sorted(aucs_e.items(), key=lambda z: -abs(z[1] - 0.5))
FA["group_e_separation"] = [dict(feature=k, auc_e_higher=a) for k, a in ranked_e]
P("     group (e) per-feature AUC ((e) higher): " + ", ".join(f"{k} {a:.2f}" for k, a in ranked_e[:12]))
OUT["failure_analysis"] = FA
OUT["seconds"] = time.time() - t0
OUT["peak_rss_gib"] = peak_gib()
(WORK / "model269.json").write_text(json.dumps(OUT, indent=1, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o)))
(WORK / "model269.log").write_text("\n".join(LINES) + "\n")
np.savez_compressed(WORK / "decisions269.npz", **{f"{s}|{m}|{'co' if co else 'all'}|{c}|realised": pol[0][c] for (s, m, co), pol in DEC.items() for c in CANCER},
                    **{f"{s}|{m}|{'co' if co else 'all'}|{c}|acted": pol[1][c] for (s, m, co), pol in DEC.items() for c in CANCER})
P(f"done in {time.time() - t0:.0f} s, peak RSS {peak_gib():.2f} GiB")
