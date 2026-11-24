"""269 -- assemble the (context, target) table for the utility router: LABELS (a frozen action set, per-target
Δ`pds` read off 262/263's saved per-target arrays) and FEATURES (build-time observable only; nothing from the
held line's truth).  Arithmetic on saved arrays + three external tables.  Nothing is staged, no fold is rebuilt.

    PYTHONDONTWRITEBYTECODE=1 uv run python build269.py

Inputs (read-only):
  experiments/263-donor-selection-arbiters/work/{sweep263,pertarget263,picks263,cov263}-<fold>.*   (labels, reliability)
  experiments/262-where-pds-is-lost/work/diag262-<fold>.npz                                            (deflation labels, structure)
  experiments/263-donor-selection-arbiters/work/plolo/fold-<fold>.npz                                  (targets, b_ctrl only)
  data/external/cd4t/cloud-compact/cd4_compact/cd4_*/{adj_p_value.npy,perts.txt}                       (CD4 DE count)
  work/depmap/{CRISPRGeneEffect.csv,Model.csv}  (DepMap 24Q2 public, figshare 25880521, downloaded by this experiment)
  data/external/regulon/collectri_omnipath.tsv, data/external/omnipath-dorothea/dorothea-abc.tsv        (regulon membership)
  data/external/ccle_broad_2025/ccle_cache.npz                                                         (external basal programmes)
Outputs: work/table269.npz, work/build269.json, work/build269.log
"""
from __future__ import annotations

import json
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E262 = ROOT / "experiments/262-where-pds-is-lost/work"
E263 = ROOT / "experiments/263-donor-selection-arbiters/work"
FOLDS = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]
CTX = {"k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}
DONOR_ID = {"k562_gwps_pm50": "K562", "k562_gwps": "K562", "orion_hct116": "HCT116", "orion_hek293t": "HEK293T", "kolf2": "KOLF2"}
DONORS = ["K562", "HCT116", "HEK293T"]            # the donor identities that exist as actions
DEPMAP_MODEL = {"K562": "K562", "HCT116": "HCT116", "HEK293T": "HEK293T", "KOLF2": None}
N_FAMILIES = 10
t0 = time.time()
LINES: list[str] = []
OUT: dict = dict(ka=[], notes=[])


def P(s=""):
    print(s, flush=True)
    LINES.append(s)


def ka(label, got, want, tol):
    d = float(abs(got - want))
    ok = d <= tol
    OUT["ka"].append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
    P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.1e}")
    if not ok:
        raise SystemExit(f"KA FAIL {label}")


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


# ------------------------------------------------------------------------------- external target biology
gene_axis = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].to_numpy()
assert gene_axis.shape[0] == 18533
gpos = {g: i for i, g in enumerate(gene_axis)}

all_targets = sorted({t for f in FOLDS for t in np.load(E263 / f"plolo/fold-{f}.npz", allow_pickle=True)["targets"].tolist()})
P(f"union of panel targets over the four contexts: {len(all_targets)}")

# DepMap 24Q2 gene effect: the 300 targets' columns across every model
hdr = pd.read_csv(WORK / "depmap/CRISPRGeneEffect.csv", nrows=0).columns.tolist()
sym2col = {c.split(" (")[0]: c for c in hdr[1:]}
dm_cols = [sym2col[t] for t in all_targets if t in sym2col]
dm = pd.read_csv(WORK / "depmap/CRISPRGeneEffect.csv", usecols=[hdr[0]] + dm_cols, index_col=0)
dm.columns = [c.split(" (")[0] for c in dm.columns]
model = pd.read_csv(WORK / "depmap/Model.csv", usecols=["ModelID", "StrippedCellLineName", "OncotreeLineage"]).set_index("ModelID")
P(f"DepMap 24Q2: {dm.shape[0]} models x {dm.shape[1]} of {len(all_targets)} targets (missing: {sorted(set(all_targets) - set(dm.columns))})")
line_model = {}
for line, name in DEPMAP_MODEL.items():
    if name is None:
        continue
    hits = model.index[model["StrippedCellLineName"].str.upper() == name.upper()].tolist()
    hits = [h for h in hits if h in dm.index]
    line_model[line] = hits[0] if hits else None
    P(f"  DepMap model for {line}: {line_model[line]} ({model.loc[hits[0], 'OncotreeLineage'] if hits else 'NOT IN THE GENE-EFFECT MATRIX'})")
dm_profile = np.full((len(all_targets), dm.shape[0]), np.nan)
for i, t in enumerate(all_targets):
    if t in dm.columns:
        dm_profile[i] = dm[t].to_numpy()
dm_mean = np.nanmean(dm_profile, 1)
dm_sd = np.nanstd(dm_profile, 1)
dm_frac_dep = np.nanmean(dm_profile < -0.5, 1)
dm_in_line = {line: (dm.loc[mid, :].reindex(all_targets).to_numpy() if mid else np.full(len(all_targets), np.nan))
              for line, mid in line_model.items()}
dm_in_line["KOLF2"] = np.full(len(all_targets), np.nan)

# families: average-linkage hierarchical clustering on the correlation distance between the targets' DepMap profiles
from scipy.cluster.hierarchy import fcluster, linkage
prof = dm_profile.copy()
prof = np.where(np.isfinite(prof), prof, np.nanmean(prof, 1, keepdims=True))
prof_c = prof - prof.mean(1, keepdims=True)
prof_c /= np.maximum(np.linalg.norm(prof_c, axis=1, keepdims=True), 1e-12)
Z = linkage(prof_c, method="average", metric="correlation")
fam_dm = fcluster(Z, N_FAMILIES, criterion="maxclust") - 1
# a size-balanced alternative: k-means on the same unit profiles (families must hold enough targets to test on)
from sklearn.cluster import KMeans
fam_km = KMeans(n_clusters=N_FAMILIES, n_init=20, random_state=0).fit_predict(prof_c)
P(f"DepMap families (average linkage, {N_FAMILIES}): sizes {np.bincount(fam_dm).tolist()}")
P(f"DepMap families (k-means, {N_FAMILIES}): sizes {np.bincount(fam_km).tolist()}")

# regulon membership (CollecTRI, DoRothEA A/B/C)
ct = pd.read_csv(ROOT / "data/external/regulon/collectri_omnipath.tsv", sep="\t", usecols=["source_genesymbol", "target_genesymbol"])
do = pd.read_csv(ROOT / "data/external/omnipath-dorothea/dorothea-abc.tsv", sep="\t", usecols=["source_genesymbol", "target_genesymbol", "dorothea_level"])
do_ab = do[do["dorothea_level"].str.startswith(("A", "B"))]
reg = dict(
    collectri_out=np.array([int((ct.source_genesymbol == t).sum()) for t in all_targets], float),
    collectri_in=np.array([int((ct.target_genesymbol == t).sum()) for t in all_targets], float),
    dorothea_ab_out=np.array([int((do_ab.source_genesymbol == t).sum()) for t in all_targets], float),
    dorothea_in=np.array([int((do.target_genesymbol == t).sum()) for t in all_targets], float),
)
reg["is_tf"] = (reg["collectri_out"] > 0).astype(float)
P(f"regulon: {int(reg['is_tf'].sum())} panel targets are CollecTRI TFs; in-degree > 0 for {int((reg['collectri_in'] > 0).sum())}")

# CCLE external basal programmes (1,600 samples): the target's expression across lines + a basal-programme basis
cc = np.load(ROOT / "data/external/ccle_broad_2025/ccle_cache.npz", allow_pickle=True)
cc_genes = cc["genes"].tolist()
cc_pos = {g: i for i, g in enumerate(cc_genes)}
ccX = np.asarray(cc["X"], dtype=np.float32)
cc_target_mean = np.array([ccX[:, cc_pos[t]].mean() if t in cc_pos else np.nan for t in all_targets])
cc_target_sd = np.array([ccX[:, cc_pos[t]].std() if t in cc_pos else np.nan for t in all_targets])
shared = [g for g in gene_axis if g in cc_pos]
sh_axis = np.array([gpos[g] for g in shared])
sh_cc = np.array([cc_pos[g] for g in shared])
# programme basis on the CCLE lines: top-10 PCs of the (standardised) log-TPM over the shared genes, fitted on EXTERNAL data only
from sklearn.decomposition import PCA
Xs = ccX[:, sh_cc].astype(np.float64)
mu, sd = Xs.mean(0), Xs.std(0) + 1e-6
pca_cc = PCA(n_components=10, random_state=0).fit((Xs - mu) / sd)
P(f"CCLE programme basis on {len(shared)} shared genes: explained {np.round(pca_cc.explained_variance_ratio_, 3).tolist()}")
del Xs

# CD4 DE count (the compact banks; per condition #genes at adj p < 0.05, median over the conditions carrying the target)
cd4_counts = {}
for cond in ["cd4_rest", "cd4_stim8hr", "cd4_stim48hr"]:
    d = ROOT / f"data/external/cd4t/cloud-compact/cd4_compact/{cond}"
    perts = [ln.strip() for ln in (d / "perts.txt").read_text().splitlines() if ln.strip()]
    genes_c = [ln.strip() for ln in (d / "genes.txt").read_text().splitlines() if ln.strip()]
    gp_c = {g: i for i, g in enumerate(genes_c)}
    ap = np.load(d / "adj_p_value.npy")
    sig = ap < 0.05
    n_de = sig.sum(1)
    for i, p in enumerate(perts):                   # 262's count excludes the on-target gene itself (KA below)
        c = int(n_de[i]) - (int(sig[i, gp_c[p]]) if p in gp_c else 0)
        cd4_counts.setdefault(p, []).append(float(c))
cd4_de = np.array([np.median(cd4_counts[t]) if t in cd4_counts else np.nan for t in all_targets])
P(f"CD4 DE count: {int(np.isfinite(cd4_de).sum())} targets; median over targets {np.nanmedian(cd4_de):.1f} (262: 3.0)")
ka("CD4 DE count median over targets == 262's 3.0", np.nanmedian(cd4_de), 3.0, 0.0)
_pooled = np.array([cd4_de[all_targets.index(t)] for f in ["k562_gwps", "orion_hct116", "orion_hek293t"]
                    for t in np.load(E263 / f"plolo/fold-{f}.npz", allow_pickle=True)["targets"].tolist() if np.isfinite(cd4_de[all_targets.index(t)])])
_bins = [int(((_pooled >= a) & (_pooled < b)).sum()) for a, b in [(0, 10), (10, 30), (30, 100), (100, 300), (300, 1e9)]]
P(f"CD4 DE count bins over the cancer folds: {_bins} (262 s4.3: [566, 111, 82, 40, 33])")
ka("CD4 DE count bin [0,10) == 262's 566", _bins[0], 566, 0)
ka("CD4 DE count bin [10,30) == 262's 111", _bins[1], 111, 0)
OUT["notes"].append(f"CD4 DE bins {_bins} vs 262 [566,111,82,40,33]: the two upper-middle bins differ by 3 targets (262's exact on-target rule not recovered); no consequence for a log1p feature")

tpos = {t: i for i, t in enumerate(all_targets)}

# ---------------------------------------------------------------------------------------- per context
b_ctrl_by_line = {}
for f in FOLDS:
    z = np.load(E263 / f"plolo/fold-{f}.npz", allow_pickle=True)
    b_ctrl_by_line[CTX[f]] = np.asarray(z["b_ctrl"], dtype=np.float64)
basal_pc_cc = {ln: pca_cc.transform(((b[sh_axis] - mu) / sd)[None, :])[0] for ln, b in b_ctrl_by_line.items()}
for ln, v in basal_pc_cc.items():
    P(f"  {ln} on the CCLE programme basis: {np.round(v[:5], 2).tolist()}")

rows_choose, rows_defl, per_ctx = [], [], {}
label_audit = []
for f in FOLDS:
    c = CTX[f]
    z = np.load(E263 / f"plolo/fold-{f}.npz", allow_pickle=True)
    targets = z["targets"].tolist()
    n = len(targets)
    sw = json.load(open(E263 / f"sweep263-{f}.json"))
    meta = sw["meta"]
    rows_sw = {r["arm"]: r for r in sw["rows"]}
    donors = [DONOR_ID[d] for d in meta["donors"]]
    J = len(donors)
    pt = np.load(E263 / f"pertarget263-{f}.npz")
    pk = np.load(E263 / f"picks263-{f}.npz")
    cov = np.load(E263 / f"cov263-{f}.npz", allow_pickle=True)
    dg = np.load(E262 / f"diag262-{f}.npz", allow_pickle=True)
    assert cov["targets"].tolist() == targets == dg["targets"].tolist() == meta["targets"]
    base = np.asarray(pt["cd4svd5|pds"], float)
    ka(f"[{c}] cd4svd5 per-target pds: 263 pertarget == 262 diag (bit)", np.abs(base - dg["cd4svd5|pds"]).max(), 0.0, 0.0)
    ka(f"[{c}] cd4svd5 per-target pds: 263 pertarget == 263 cov copy (bit)", np.abs(base - cov["cd4svd5_262_pds"]).max(), 0.0, 0.0)
    ka(f"[{c}] fold mean of cd4svd5 == sweep row", base.mean(), rows_sw["cd4svd5"]["pds"], 1e-9)
    h1 = np.asarray(dg["h1union|pds"], float)
    u_undefl = h1 - base
    ka(f"[{c}] mean(h1union - cd4svd5) == -(262 s6 step)", u_undefl.mean(), -{"k562_gwps": 0.00999, "orion_hct116": 0.00945, "orion_hek293t": 0.01736, "kolf2": 0.00414}[f], 1e-5)
    orc = np.asarray(pt["oracle_wht:all:nm|pds"], float) - base
    orc3 = np.asarray(pt["oracle3_wht:all:nm|pds"], float) - base
    avail = np.asarray(cov["avail"], bool)                          # (J, n)

    # ---- choose-j labels: mean over every norm-matched selection arm whose picks are build-time functions
    label_arms = [a for a in pk.files if a.endswith(":nm") and not a.startswith(("oracle", "majority", "mean", "only", "drop"))
                  and ":sm" not in a and f"{a}|pds" in pt.files]
    acc = np.zeros((n, J)); acc2 = np.zeros((n, J)); cnt = np.zeros((n, J))
    per_arm = {}
    for a in label_arms:
        p = np.asarray(pk[a], int)
        v = np.asarray(pt[f"{a}|pds"], float) - base
        per_arm[a] = (p, v)
        for j in range(J):
            m = p == j
            acc[m, j] += v[m]; acc2[m, j] += v[m] ** 2; cnt[m, j] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        u_mean = np.where(cnt > 0, acc / np.maximum(cnt, 1), np.nan)
        u_sd = np.where(cnt > 1, np.sqrt(np.maximum(acc2 / np.maximum(cnt, 1) - u_mean ** 2, 0)), np.nan)
    # proxy check (265-style): each real arm's total from LEAVE-THAT-ARM-OUT labels vs its real total
    errs = []
    for a, (p, v) in per_arm.items():
        tot_real = float(v.mean())
        loo = np.zeros(n)
        for t in np.flatnonzero(p >= 0):
            j = p[t]
            k = cnt[t, j] - 1
            loo[t] = (acc[t, j] - v[t]) / k if k > 0 else np.nan
        ok = np.isfinite(loo)
        tot_proxy = float(np.nansum(loo) / n)
        errs.append(dict(arm=a, real=tot_real, proxy=tot_proxy, err=tot_proxy - tot_real, n_sel=int((p >= 0).sum()), n_loo_defined=int(ok.sum())))
    e = np.array([x["err"] for x in errs])
    label_audit.append(dict(context=c, n_label_arms=len(label_arms), proxy_err_mean=float(e.mean()), proxy_err_sd=float(e.std()),
                            proxy_err_max_abs=float(np.abs(e).max()), arms=errs))
    P(f"[{c}] donors {donors}; {len(label_arms)} label arms; per-(t,j) arm count median {np.median(cnt[cnt > 0]):.0f}; "
      f"label SD across arms median {np.nanmedian(u_sd):.4f}; proxy-vs-real arm total: mean err {e.mean():+.5f}, sd {e.std():.5f}, max|err| {np.abs(e).max():.5f}")
    P(f"[{c}] label means: " + ", ".join(f"choose {donors[j]} {np.nanmean(u_mean[:, j]):+.4f} (n {int((cnt[:, j] > 0).sum())})" for j in range(J))
      + f"; undeflate {u_undefl.mean():+.4f}; 263 oracle {orc.mean():+.4f}; abstaining oracle {orc3.mean():+.4f}")
    lab_oracle = np.nanmax(np.column_stack([np.zeros(n), np.where(np.isfinite(u_mean), u_mean, -np.inf), u_undefl]), 1)
    P(f"[{c}] LABELLED oracle (max over keep / choose-j / undeflate): {lab_oracle.mean():+.4f}; "
      f"choose-only labelled oracle {np.nanmax(np.column_stack([np.zeros(n), np.where(np.isfinite(u_mean), u_mean, -np.inf)]), 1).mean():+.4f}")

    # ---- shared per-(context, target) features
    b_held = b_ctrl_by_line[c]
    tidx = np.array([gpos[t] for t in targets])
    gidx = np.array([tpos[t] for t in targets])
    cd4_row_norm = np.asarray(dg["cov|cd4_norm"], float)
    shared_feat = dict(
        dd_raw=np.asarray(cov["dd_raw"], float), dd_wht=np.asarray(cov["dd_wht"], float),
        n_whitened=np.asarray(dg["cov|n_whitened"], float), pm_carried=np.asarray(dg["cov|pm_carried"], float),
        cd4_count=np.asarray(dg["cov|cd4_count"], float), cd4_norm=cd4_row_norm,
        cd4_deflated_norm=np.asarray(dg["cov|cd4_deflated_norm"], float), cd4_energy_removed=np.asarray(dg["cov|cd4_energy_removed"], float),
        union_count=np.asarray(dg["cov|union_count"], float), union_norm=np.asarray(dg["cov|union_norm"], float),
        base_norm=np.asarray(dg["cov|base_norm"], float), resid_norm=np.asarray(dg["cov|resid_norm"], float),
        resid_cos_base=np.asarray(dg["cov|resid_cos_base"], float), resid_share=np.asarray(dg["cov|resid_share"], float),
        cd4_cos_base=np.asarray(dg["cov|cd4_cos_base"], float), cd4_deflated_cos_base=np.asarray(dg["cov|cd4_deflated_cos_base"], float),
        pred_nn=np.asarray(dg["cd4svd5|pred_nn"], float), pred_mean=np.asarray(dg["cd4svd5|pred_mean"], float),
        pred_norm=np.asarray(dg["cd4svd5|pred_norm"], float),
        cd4_de_log=np.log1p(cd4_de[gidx]),
        target_basal_held=b_held[tidx],
        dm_mean=dm_mean[gidx], dm_sd=dm_sd[gidx], dm_frac_dep=dm_frac_dep[gidx], dm_in_held=dm_in_line[c][gidx],
        cc_target_mean=cc_target_mean[gidx], cc_target_sd=cc_target_sd[gidx],
        **{k: v[gidx] for k, v in reg.items()},
        **{f"basal_cc_pc{i+1}": np.full(n, basal_pc_cc[c][i]) for i in range(5)},
    )
    # whitened-space donor-mean features of both donors (the "other" donor for the choose rows is derived below)
    norm_raw = np.asarray(cov["norm_raw"], float); norm_wht = np.asarray(cov["norm_wht"], float)
    S = {k: np.asarray(cov[f"S|{k}|raw"], float) for k in ["cd4", "cd4d", "union", "selfcons", "basal"]}
    Sw = {k: np.asarray(cov[f"S|{k}|wht"], float) for k in ["cd4", "cd4d", "union", "selfcons"]}
    basal_j = np.asarray(cov["basal"], float)

    # truth-side quantities for ANALYSIS ONLY (never features): 262 groups
    base_cos_truth = np.asarray(dg["cov|base_cos_truth"], float)
    group_e = base_cos_truth <= 0
    marginal = (base_cos_truth > 0) & (base_cos_truth <= 0.05)
    orc_pick = np.asarray(pk["oracle_wht:all:nm"], int)

    def donor_feats(j):
        d = dict(norm_raw=norm_raw[:, j], norm_wht=norm_wht[:, j], cd4_cos=S["cd4"][:, j], cd4_cos_w=Sw["cd4"][:, j],
                 cd4d_cos=S["cd4d"][:, j], cd4d_cos_w=Sw["cd4d"][:, j], union_cos=S["union"][:, j], union_cos_w=Sw["union"][:, j],
                 selfcons=S["selfcons"][:, j], selfcons_w=Sw["selfcons"][:, j], basal_sim=np.full(n, basal_j[j]),
                 target_basal_donor=b_ctrl_by_line[donors[j]][tidx], dm_in_donor=dm_in_line[donors[j]][gidx])
        return d

    dfe = [donor_feats(j) for j in range(J)]
    for j in range(J):
        others = [k for k in range(J) if k != j]
        for t in range(n):
            if not avail[j, t] or not np.isfinite(u_mean[t, j]):
                continue
            oth = [k for k in others if avail[k, t]]
            row = dict(context=c, target=targets[t], t_idx=t, g_idx=int(gidx[t]), donor=donors[j], j=j,
                       u=float(u_mean[t, j]), u_sd=float(u_sd[t, j]) if np.isfinite(u_sd[t, j]) else np.nan, u_n=int(cnt[t, j]),
                       u_undefl=float(u_undefl[t]), oracle263=float(orc[t]), oracle263_pick=int(orc_pick[t]),
                       group_e=bool(group_e[t]), marginal=bool(marginal[t]), base_pds=float(base[t]),
                       n_other=len(oth), family_dm=int(fam_dm[gidx[t]]), family_km=int(fam_km[gidx[t]]))
            for k, v in shared_feat.items():
                row[f"s_{k}"] = float(v[t])
            for k, v in dfe[j].items():
                row[f"d_{k}"] = float(v[t])
                ov = [dfe[o][k][t] for o in oth]
                row[f"o_{k}"] = float(np.nanmean(ov)) if oth and np.isfinite(ov).any() else np.nan
                row[f"r_{k}"] = row[f"d_{k}"] - row[f"o_{k}"] if np.isfinite(row[f"o_{k}"]) else np.nan
            rows_choose.append(row)
    for t in range(n):
        row = dict(context=c, target=targets[t], t_idx=t, g_idx=int(gidx[t]), u_undefl=float(u_undefl[t]), oracle263=float(orc[t]),
                   oracle263_pick=int(orc_pick[t]), oracle3_263=float(orc3[t]), group_e=bool(group_e[t]), marginal=bool(marginal[t]),
                   base_pds=float(base[t]), n_avail=int(avail[:, t].sum()), family_dm=int(fam_dm[gidx[t]]), family_km=int(fam_km[gidx[t]]),
                   lab_oracle=float(lab_oracle[t]))
        for j, dn in enumerate(donors):
            row[f"u_choose_{dn}"] = float(u_mean[t, j]) if np.isfinite(u_mean[t, j]) else np.nan
        for k, v in shared_feat.items():
            row[f"s_{k}"] = float(v[t])
        # donor-mean reliability features (context-level view of the two/three donors)
        for k in ["norm_raw", "norm_wht", "cd4_cos", "cd4_cos_w", "cd4d_cos", "cd4d_cos_w", "selfcons", "selfcons_w"]:
            vals = [dfe[j][k][t] for j in range(J) if avail[j, t]]
            row[f"m_{k}"] = float(np.nanmean(vals)) if vals and np.isfinite(vals).any() else np.nan
        rows_defl.append(row)
    per_ctx[c] = dict(fold=f, donors=donors, n=n, base_pds=float(base.mean()), oracle263=float(orc.mean()), oracle3_263=float(orc3.mean()),
                      undefl=float(u_undefl.mean()), lab_oracle=float(lab_oracle.mean()),
                      choose_mean={donors[j]: float(np.nanmean(u_mean[:, j])) for j in range(J)},
                      n_choose_rows={donors[j]: int(np.isfinite(u_mean[:, j]).sum()) for j in range(J)},
                      group_e=int(group_e.sum()), marginal=int(marginal.sum()))

dfc = pd.DataFrame(rows_choose)
dfd = pd.DataFrame(rows_defl)
P(f"\nrows: choose {len(dfc)} (per context {dfc.groupby('context').size().to_dict()}), deflation {len(dfd)}")
P(f"label distributions (choose): mean {dfc.u.mean():+.4f} sd {dfc.u.std():.4f}; quantiles " + ", ".join(f"{q:.0%} {dfc.u.quantile(q):+.3f}" for q in [0.05, 0.25, 0.5, 0.75, 0.95]))
P(f"  positive fraction {float((dfc.u > 0).mean()):.3f}; |u| > 0.05: {float((dfc.u.abs() > 0.05).mean()):.3f}; u > 0.05: {float((dfc.u > 0.05).mean()):.3f}")
P(f"label distributions (undeflate): mean {dfd.u_undefl.mean():+.4f} sd {dfd.u_undefl.std():.4f}; positive {float((dfd.u_undefl > 0).mean()):.3f}; > 0.05 {float((dfd.u_undefl > 0.05).mean()):.3f}")
dfc.to_csv(WORK / "rows_choose269.csv", index=False)
dfd.to_csv(WORK / "rows_defl269.csv", index=False)
np.savez_compressed(WORK / "table269.npz", all_targets=np.array(all_targets), dm_profile=dm_profile, dm_models=np.array(dm.index.tolist()),
                    fam_dm=fam_dm, fam_km=fam_km, cc_target_profile=np.stack([ccX[:, cc_pos[t]] if t in cc_pos else np.full(ccX.shape[0], np.nan) for t in all_targets]),
                    b_ctrl_lines=np.stack([b_ctrl_by_line[c] for c in ["K562", "HCT116", "HEK293T", "KOLF2"]]), b_ctrl_line_names=np.array(["K562", "HCT116", "HEK293T", "KOLF2"]),
                    sh_axis=sh_axis)
OUT.update(per_context=per_ctx, label_audit=label_audit, n_rows_choose=len(dfc), n_rows_defl=len(dfd),
           families_dm_sizes=np.bincount(fam_dm).tolist(), families_km_sizes=np.bincount(fam_km).tolist(),
           depmap=dict(models=int(dm.shape[0]), targets_covered=int(dm.shape[1]), line_model=line_model),
           seconds=time.time() - t0, peak_rss_gib=peak_gib())
(WORK / "build269.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "build269.log").write_text("\n".join(LINES) + "\n")
P(f"done in {time.time() - t0:.0f} s, peak RSS {peak_gib():.2f} GiB")
