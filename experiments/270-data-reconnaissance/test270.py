"""270 — frozen tests of the GPSAdb sign-stability feature against 262's per-target arrays.
Known-answer rows first (must reproduce 262), then the feature tests if work/feature270.csv exists.
Reads only experiments/262-where-pds-is-lost/work/*; writes only into experiments/270-data-reconnaissance/work/.
"""
import json, os, sys, numpy as np, pandas as pd
from scipy.stats import spearmanr
R = "~/Code/virtual-cell-challenge-2026/experiments/262-where-pds-is-lost/work"
W = "~/Code/virtual-cell-challenge-2026/experiments/270-data-reconnaissance/work"
CANCER = ["k562_gwps", "orion_hct116", "orion_hek293t"]; FOLDS = CANCER + ["kolf2"]
LINE = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}
rng = np.random.default_rng(20260905)
D = {f: np.load(f"{R}/diag262-{f}.npz", allow_pickle=True) for f in FOLDS}
M = {f: json.load(open(f"{R}/diag262-{f}.json"))["meta"] for f in FOLDS}
CACHE = {f: np.load(f"{R}/plolo/fold-{f}.npz", allow_pickle=True) for f in FOLDS}
out = {"ka": {}, "tests": {}}
# ---------------- KA rows
for f in FOLDS:
    z = D[f]; t = [str(x) for x in z["targets"]]
    assert list(t) == [str(x) for x in CACHE[f]["targets"]], "target order mismatch between diag and cache"
    rho = spearmanr(z["cov|base_cos_truth"], z["cd4svd5|pds"]).correlation
    out["ka"][f] = {"n": len(t), "mean_pds_cd4svd5": float(z["cd4svd5|pds"].mean()), "rho_base_cos_pds": float(rho)}
    print(f"KA {f}: n {len(t)} mean cd4svd5 pds {z['cd4svd5|pds'].mean():.6f} (262 KA-1: 0.815384/0.801659/0.780978/0.640294); "
          f"rho(base·truth, pds) {rho:+.3f} (262 §4.1: +0.913/+0.885/+0.855)")
# held-vs-donor cross-line truth agreement, recomputed exactly as analyse262.py (pooled by within-fold rank -> +0.630)
def unit_rows(A): return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)
truths = {f: np.asarray(CACHE[f]["truth"], dtype=np.float64) for f in FOLDS}
scored = {f: np.asarray(CACHE[f]["scored"]) for f in FOLDS}
tlist = {f: [str(x) for x in CACHE[f]["targets"]] for f in FOLDS}
tpos = {f: {t: i for i, t in enumerate(tlist[f])} for f in FOLDS}
pair_cos = {}
for i, a in enumerate(FOLDS):
    for b in FOLDS[i + 1:]:
        cols = np.intersect1d(scored[a], scored[b]); common = [t for t in tlist[a] if t in tpos[b]]
        A = truths[a][[tpos[a][t] for t in common]][:, cols]; B = truths[b][[tpos[b][t] for t in common]][:, cols]
        c = (unit_rows(A) * unit_rows(B)).sum(1); pair_cos[(a, b)] = dict(zip(common, c.tolist())); pair_cos[(b, a)] = pair_cos[(a, b)]
def ranks(x):
    x = np.asarray(x, float); r = np.empty_like(x); r[np.argsort(x)] = np.arange(len(x)); return pd.Series(x).rank().values
pooled_x, pooled_y = [], []
held = {}
for f in CANCER:
    donors = M[f]["sources"]; t = tlist[f]
    hv = np.array([np.mean([pair_cos[(f, d)].get(x, np.nan) for d in donors]) for x in t]); held[f] = hv
    ok = ~np.isnan(hv); pooled_x.append(ranks(hv[ok])); pooled_y.append(ranks(D[f]["cd4svd5|pds"][ok]))
rho_pool = spearmanr(np.concatenate(pooled_x), np.concatenate(pooled_y)).correlation
out["ka"]["rho_pooled_held_vs_donor_agreement_pds"] = float(rho_pool)
print(f"KA pooled rho(held-vs-donor truth agreement, pds) = {rho_pool:+.3f} (262 §4.1: +0.630)")
# ---------------- per-target reference table (cancer folds), kolf2 apart
rows = {}
for f in CANCER:
    z = D[f]
    for i, t in enumerate(tlist[f]):
        r = rows.setdefault(t, {"pds": [], "base": [], "held": [], "h1union_pds": []})
        r["pds"].append(float(z["cd4svd5|pds"][i])); r["base"].append(float(z["cov|base_cos_truth"][i]))
        r["held"].append(float(held[f][i])); r["h1union_pds"].append(float(z["h1union|pds"][i]))
ref = pd.DataFrame([{"target": t, "n_folds": len(r["pds"]), "pds_mean": np.mean(r["pds"]), "base_cos_mean": np.mean(r["base"]),
                     "frac_folds_anti": np.mean([b <= 0 for b in r["base"]]), "held_agree_mean": np.mean(r["held"]),
                     "anti_aligned": int(np.mean([b <= 0 for b in r["base"]]) >= 0.5)} for t, r in rows.items()])
zk = D["kolf2"]; kolf = pd.DataFrame({"target": tlist["kolf2"], "pds_kolf2": zk["cd4svd5|pds"], "base_cos_kolf2": zk["cov|base_cos_truth"]})
ref = ref.merge(kolf, on="target", how="left")
ref.to_csv(f"{W}/ref262_per_target.csv", index=False)
print(f"reference table: {len(ref)} targets on >=1 cancer fold; anti-aligned (majority rule) {int(ref.anti_aligned.sum())} "
      f"({ref.anti_aligned.mean():.1%}); mean pds {ref.pds_mean.mean():.4f}")
# ---------------- frozen tests
def boot_rho(x, y, B=2000):
    n = len(x); vals = []
    for _ in range(B):
        idx = rng.integers(0, n, n); vals.append(spearmanr(x[idx], y[idx]).correlation)
    return np.nanpercentile(vals, [2.5, 97.5])
def perm_p(x, y, B=2000):
    obs = spearmanr(x, y).correlation; cnt = 0
    for _ in range(B):
        cnt += abs(spearmanr(rng.permutation(x), y).correlation) >= abs(obs)
    return (cnt + 1) / (B + 1)
def auc(score, pos):  # P(score(pos) > score(neg)), ties 0.5
    s = np.asarray(score, float); p = np.asarray(pos, bool); a, b = s[p], s[~p]
    if len(a) == 0 or len(b) == 0: return np.nan
    return float(((a[:, None] > b[None, :]).mean() + 0.5 * (a[:, None] == b[None, :]).mean()))
def boot_auc(score, pos, B=2000):
    n = len(score); vals = []
    for _ in range(B):
        idx = rng.integers(0, n, n); vals.append(auc(score[idx], pos[idx]))
    return np.nanpercentile(vals, [2.5, 97.5])
def run_tests(fp, label, control_cols, primary):
    if not os.path.exists(fp): print(f"\n[{label}] no feature file at {fp}"); return
    F = pd.read_csv(fp); T = ref.merge(F, on="target", how="inner")
    print(f"\n[{label}] feature table: {len(F)} panel targets; matched to 262 reference: {len(T)}")
    cols = [c for c in F.columns if c.startswith("SS_") or c.startswith("RHO_")]
    for c in cols:
        x = T[c].values.astype(float); ok = ~np.isnan(x)
        if ok.sum() < 10: print(f"  {c}: only {ok.sum()} defined, skipped"); continue
        x = x[ok]; y = T.pds_mean.values[ok]; anti = T.anti_aligned.values[ok].astype(bool)
        rho = spearmanr(x, y).correlation; lo, hi = boot_rho(x, y); pp = perm_p(x, y)
        a = auc(x, ~anti); alo, ahi = boot_auc(x, ~anti)
        km = ok & ~T.pds_kolf2.isna().values
        rk = spearmanr(T[c].values[km], T.pds_kolf2.values[km]).correlation if km.sum() >= 10 else np.nan
        hm = ok & ~np.isnan(T.held_agree_mean.values)
        rh = spearmanr(T[c].values[hm], T.held_agree_mean.values[hm]).correlation if hm.sum() >= 10 else np.nan
        tag = "PRIMARY" if c == primary else "secondary"
        out["tests"][f"{label}:{c}"] = dict(n=int(ok.sum()), rho_pds=float(rho), rho_ci=[float(lo), float(hi)], perm_p=float(pp), auc_rest_vs_anti=float(a),
                               auc_ci=[float(alo), float(ahi)], n_anti=int(anti.sum()), rho_pds_kolf2=float(rk), rho_held_agreement=float(rh), tag=tag,
                               feature_mean=float(np.mean(x)), feature_median=float(np.median(x)))
        print(f"  [{tag}] {c}: n {ok.sum()} (anti {anti.sum()}) mean {np.mean(x):.3f} | rho(feature, mean cancer-fold pds) {rho:+.3f} CI[{lo:+.3f},{hi:+.3f}] perm p {pp:.3f} | "
              f"AUC rest>anti {a:.3f} CI[{alo:.3f},{ahi:.3f}] | rho kolf2 pds {rk:+.3f} | rho held-vs-donor agreement {rh:+.3f}")
    for c in control_cols:
        if c in T:
            x = T[c].values.astype(float); ok = ~np.isnan(x)
            if ok.sum() < 10: continue
            rho = spearmanr(x[ok], T.pds_mean.values[ok]).correlation; lo, hi = boot_rho(x[ok], T.pds_mean.values[ok])
            out["tests"][f"{label}:control_{c}"] = dict(rho=float(rho), ci=[float(lo), float(hi)], n=int(ok.sum()))
            print(f"  [control] {c}: rho with pds {rho:+.3f} CI[{lo:+.3f},{hi:+.3f}] n {ok.sum()}")
    T.to_csv(f"{W}/test270_{label}_merged.csv", index=False)
run_tests(f"{W}/feature270.csv", "gpsadb", ["n_groups", "n_studies", "n_lines"], "SS_cross")
run_tests(f"{W}/feature270_cmap.csv", "cmap", ["n_lines_sh", "n_lines_xpr", "n_pairs_sh"], "SS_cmap_sh")
json.dump(out, open(f"{W}/test270.json", "w"), indent=1)
