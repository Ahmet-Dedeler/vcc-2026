#!/usr/bin/env python3
"""700 EXPLORATORY: does PIE predict per-target response SIZE better than our key? (518: a per-target magnitude oracle at our
directions is worth +0.018..+0.023 overall.) h1-test truth: ||truth effect|| (pseudobulk, panel genes removed) and the
least-squares scale of our key direction m_t = <key, truth>/|key|. Predictors: |key| (ours, the LODO-like key), |PIE
effect|, PIE sum p_de, PIE sum |delta_p|. Spearman over targets, all / seen."""
import json, sys
from pathlib import Path
import numpy as np, pyarrow.parquet as pq
from scipy.stats import spearmanr
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal")); import lib521 as L
S = L.load_state(0); targets, genes = list(S["targets"]), list(S["genes"]); T, G = len(targets), len(genes)
t = pq.read_table(sys.argv[1]); pg = json.loads(t.schema.metadata[b"pie"])["genes"]["h1test"]; df = t.to_pandas()
pos = {g: i for i, g in enumerate(pg)}; have = np.array([g in pos for g in genes]); ix = np.array([pos.get(g, 0) for g in genes])
import anndata as ad
a = ad.read_h5ad(ROOT / "data/harness/h1-test.h5ad"); vp = {g: i for i, g in enumerate(a.var_names)}; cols = np.array([vp[g] for g in genes])
tg = a.obs["target_gene"].to_numpy(); X = a.X.tocsr()
def pb(m):
    P = np.asarray(X[m].sum(0)).ravel().astype(float); return np.log1p(5e4 * P / P.sum())[cols]
bc = pb(tg == "non-targeting"); tb = np.stack([pb(tg == x) - bc for x in targets]); keep = ~np.isin(np.array(genes), np.array(targets))
pc = np.expm1(bc); pc /= pc.sum()
E = np.zeros((T, G)); PDE = np.zeros(T); DP = np.zeros(T)
for _, r in df.iterrows():
    i = targets.index(r["perturbation"]); lf = np.where(have, np.asarray(r["lfc_pred"])[ix], 0)
    q = pc * np.exp2(lf); q /= q.sum(); E[i] = np.log1p(5e4 * q) - np.log1p(5e4 * pc)
    PDE[i] = float(np.asarray(r["p_de"]).sum()); DP[i] = float(np.abs(np.asarray(r["delta_p_pred"])).sum())
key = S["key"]
tn = np.linalg.norm(tb[:, keep], axis=1); kn = np.linalg.norm(key[:, keep], axis=1); en = np.linalg.norm(E[:, keep], axis=1)
m_t = (key[:, keep] * tb[:, keep]).sum(1) / np.maximum(kn, 1e-12)      # truth's projection on our unit direction
seen = set(json.load(open(HERE / "work/pie_train_perts.json"))); sm = np.array([x in seen for x in targets])
out = {}
for lab, msk in (("all", np.ones(T, bool)), ("seen", sm)):
    row = {}
    for tgt_name, tv in (("|truth|", tn), ("proj_on_our_dir", m_t)):
        for pn, pv in (("|key|", kn), ("|PIE effect|", en), ("PIE sum p_de", PDE), ("PIE sum|dp|", DP)):
            row[f"{tgt_name} ~ {pn}"] = round(float(spearmanr(pv[msk], tv[msk]).correlation), 3)
    out[lab] = row
print(json.dumps(out, indent=1))
