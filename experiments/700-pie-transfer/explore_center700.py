#!/usr/bin/env python3
"""700 EXPLORATORY (not preregistered; recorded as a try): does PIE's gain come from a shared response component?
Blend variants on h1-test and K562, seen targets only, w 0.25: (u) PIE effect rows as is; (c) PIE rows centred across the
panel's targets (each gene's mean over targets removed), i.e. only PIE's target-specific part. Reports ew cos, minipds."""
import json, sys
from pathlib import Path
import numpy as np, pyarrow.parquet as pq
HERE = Path(__file__).resolve().parent; ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal")); import lib521 as L
seen = set(json.load(open(HERE / "work/pie_train_perts.json")))
def ublend(base, pie, w, only):
    out = base.copy()
    for i in range(len(base)):
        if not only[i]: continue
        nb, nv = np.linalg.norm(base[i]), np.linalg.norm(pie[i])
        if nb == 0 or nv == 0: continue
        m = (1 - w) * base[i] / nb + w * pie[i] / nv; out[i] = nb * m / np.linalg.norm(m)
    return out
def metrics(pred, tb):
    na, nb = np.linalg.norm(pred, axis=1), np.linalg.norm(tb, axis=1)
    c = (pred * tb).sum(1) / np.maximum(na * nb, 1e-300); w = nb ** 2
    an = pred / np.maximum(na[:, None], 1e-12); bn = tb / np.maximum(nb[:, None], 1e-12); C = an @ bn.T; n = len(pred)
    r = [(C[i] > C[i, i]).sum() + 0.5 * ((C[i] == C[i, i]).sum() - 1) for i in range(n)]
    return round(float((c * w).sum() / w.sum()), 4), round(float(1 - np.mean(r) / (n - 1)), 4)
# ---- H1
S = L.load_state(0); targets, genes = list(S["targets"]), list(S["genes"]); T, G = len(targets), len(genes)
t = pq.read_table(sys.argv[1]); pg = json.loads(t.schema.metadata[b"pie"])["genes"]["h1test"]; df = t.to_pandas()
pos = {g: i for i, g in enumerate(pg)}; have = np.array([g in pos for g in genes]); ix = np.array([pos.get(g, 0) for g in genes])
import anndata as ad
a = ad.read_h5ad(ROOT / "data/harness/h1-test.h5ad"); vpos = {g: i for i, g in enumerate(a.var_names)}; cols = np.array([vpos[g] for g in genes])
tg = a.obs["target_gene"].to_numpy(); X = a.X.tocsr()
def pb(m):
    P = np.asarray(X[m].sum(0)).ravel().astype(float); return np.log1p(5e4 * P / P.sum())[cols]
bc = pb(tg == "non-targeting"); tb = np.stack([pb(tg == x) - bc for x in targets])
keep = ~np.isin(np.array(genes), np.array(targets))
pc = np.expm1(bc); pc /= pc.sum(); E = np.zeros((T, G))
for _, r in df.iterrows():
    i = targets.index(r["perturbation"]); lf = np.where(have, np.asarray(r["lfc_pred"])[ix], 0)
    q = pc * np.exp2(np.where(1e6 * pc >= 5.0, lf, 0.0)); q /= q.sum(); E[i] = np.log1p(5e4 * q) - np.log1p(5e4 * pc)
key = S["key"]; sm = np.array([x in seen for x in targets])
Ec = E - E.mean(0, keepdims=True)
for nm, v in (("base", None), ("uncentred", E), ("centred", Ec)):
    p = key if v is None else ublend(key[:, keep], v[:, keep], 0.25, sm)
    p = p[:, keep] if v is None else p
    print("H1  ", nm, metrics(p, tb[:, keep]))
# ---- K562 (donor mean base)
D = ROOT / "data/donors"
def donor(n): return (D/n/"genes.txt").read_text().split(), (D/n/"perts.txt").read_text().split(), np.load(D/n/"delta.npy", mmap_mode="r")
t2 = pq.read_table(sys.argv[2]); pg2 = json.loads(t2.schema.metadata[b"pie"])["genes"]["replogle"]; d2 = t2.to_pandas()
kg, kp, kd = donor("k562_gwps"); hg, hp, hd = donor("orion_hct116"); eg, ep, ed = donor("orion_hek293t")
import csv
panel = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv")) if r["target_gene"] != "non-targeting"]
common = sorted(set(pg2) & set(kg) & set(hg) & set(eg) - set(panel))
ix2 = {n: {g: i for i, g in enumerate(gl)} for n, gl in (("p", pg2), ("k", kg), ("h", hg), ("e", eg))}
ci = {n: np.array([ix2[n][g] for g in common]) for n in ix2}
kpi, hpi, epi = ({p: i for i, p in enumerate(pl)} for pl in (kp, hp, ep))
P2 = {r["perturbation"]: np.asarray(r["lfc_pred"], float) for _, r in d2.iterrows()}
tg2 = [x for x in panel if x in P2 and x in kpi]
bck = np.load(D / "k562_gwps" / "b_ctrl.npy").astype(float)[ci["k"]]; pck = np.expm1(bck); pck /= pck.sum()
z = lambda v: (v - v.mean()) / v.std() if v.std() > 0 else v * 0
TB, DM, EP = [], [], []
for x in tg2:
    TB.append(np.asarray(kd[kpi[x]], float)[ci["k"]])
    h = np.asarray(hd[hpi[x]], float)[ci["h"]] if x in hpi else np.zeros(len(common))
    e = np.asarray(ed[epi[x]], float)[ci["e"]] if x in epi else np.zeros(len(common))
    DM.append((z(h) + z(e)) / max(1, (x in hpi) + (x in epi)))
    q = pck * np.exp2(np.where(1e6 * pck >= 5.0, P2[x][ci["p"]], 0.0)); q /= q.sum(); EP.append(np.log1p(5e4 * q) - np.log1p(5e4 * pck))
TB, DM, EP = map(np.stack, (TB, DM, EP)); sm2 = np.array([x in seen for x in tg2]); EPc = EP - EP.mean(0, keepdims=True)
for nm, v in (("base", None), ("uncentred", EP), ("centred", EPc)):
    p = DM if v is None else ublend(DM, v, 0.25, sm2)
    print("K562", nm, metrics(p, TB))
