#!/usr/bin/env python3
"""702 post-mortem: were the held-out gains distinguishable from noise, and did the instruments see where the change sat?

Target bootstrap (4,000 draws, seed 702) of the energy-weighted-cosine (ew) gain of 702's blend (seen targets, w 0.25,
norm-preserving, 5 CPM rule) over each surface's base: K562 GWPS (Orion donor-mean base, eval700_k562's construction) and
h1-test (435's key base, eval700_h1's construction). Gene subsets by the surface's own control CPM: all genes of its axis,
CPM >= 20, CPM < 20. Target sets: the seen targets (where the blend acts) and all targets (what amendment 1's prediction
translated). No new choices: the constructions are copied from eval700_*.
usage: uv run python overread702.py [--k562-only]   (--k562-only keeps the h1-test rows already in work/overread702.json)
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CLOUD = ROOT / "experiments/700-pie-transfer/work/cloud/20261014T225548Z/out"
seen = set(json.load(open(ROOT / "experiments/700-pie-transfer/work/pie_train_perts.json")))
rng = np.random.default_rng(702)


def ublend(base, pie, w, only):
    out = base.copy()
    for i in range(len(base)):
        if not only[i]:
            continue
        nb, nv = np.linalg.norm(base[i]), np.linalg.norm(pie[i])
        if nb == 0 or nv == 0:
            continue
        m = (1 - w) * base[i] / nb + w * pie[i] / nv
        out[i] = nb * m / np.linalg.norm(m)
    return out


def ewcos_rows(pred, tb):
    na, nb = np.linalg.norm(pred, axis=1), np.linalg.norm(tb, axis=1)
    return (pred * tb).sum(1) / np.maximum(na * nb, 1e-300), nb ** 2


def boot(c_new, c_base, w, mask):
    idx = np.flatnonzero(mask)
    d = []
    for _ in range(4000):
        s = rng.choice(idx, len(idx))
        d.append(((c_new[s] - c_base[s]) * w[s]).sum() / w[s].sum())
    point = ((c_new[idx] - c_base[idx]) * w[idx]).sum() / w[idx].sum()
    base_ew = (c_base[idx] * w[idx]).sum() / w[idx].sum()
    return {"gain": round(float(point), 5), "ci95": [round(float(np.quantile(d, .025)), 5), round(float(np.quantile(d, .975)), 5)],
            "base_ew": round(float(base_ew), 5), "n": int(len(idx))}


def score(out, name, base_all, pie_all, tb, cpm, gene_ok, sm):
    for lab, gm in (("all genes", gene_ok), ("CPM>=20", gene_ok & (cpm >= 20)), ("CPM<20", gene_ok & (cpm < 20))):
        base, new = base_all[:, gm], ublend(base_all[:, gm], pie_all[:, gm], 0.25, sm)
        cb, w = ewcos_rows(base, tb[:, gm])
        cn, _ = ewcos_rows(new, tb[:, gm])
        out[f"{name} seen | {lab}"] = boot(cn, cb, w, sm)
        out[f"{name} all targets | {lab}"] = boot(cn, cb, w, np.ones(len(sm), bool))


def h1_part(out):
    sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
    import anndata as ad
    import lib521 as L
    S = L.load_state(0)
    targets, genes = list(S["targets"]), list(S["genes"])
    T, G = len(targets), len(genes)
    t = pq.read_table(CLOUD / "pred_h1test.parquet")
    pg = json.loads(t.schema.metadata[b"pie"])["genes"]["h1test"]
    df = t.to_pandas()
    pos = {g: i for i, g in enumerate(pg)}
    have = np.array([g in pos for g in genes])
    ix = np.array([pos.get(g, 0) for g in genes])
    a = ad.read_h5ad(ROOT / "data/harness/h1-test.h5ad")
    vp = {g: i for i, g in enumerate(a.var_names)}
    cols = np.array([vp[g] for g in genes])
    tg = a.obs["target_gene"].to_numpy()
    X = a.X.tocsr()

    def pb(m):
        P = np.asarray(X[m].sum(0)).ravel().astype(float)
        return np.log1p(5e4 * P / P.sum())[cols], P[cols]
    bc, Pc = pb(tg == "non-targeting")
    tb = np.stack([pb(tg == x)[0] - bc for x in targets])
    cpm = 1e6 * Pc / Pc.sum()
    keep = ~np.isin(np.array(genes), np.array(targets))
    pc = np.expm1(bc)
    pc /= pc.sum()
    E = np.zeros((T, G))
    for _, r in df.iterrows():
        i = targets.index(r["perturbation"])
        lf = np.where(have, np.asarray(r["lfc_pred"])[ix], 0)
        q = pc * np.exp2(np.where(1e6 * pc >= 5.0, lf, 0.0))
        q /= q.sum()
        E[i] = np.log1p(5e4 * q) - np.log1p(5e4 * pc)
    sm = np.array([x in seen for x in targets])
    score(out, "h1test", S["key"], E, tb, cpm, keep, sm)


def k562_part(out):
    D = ROOT / "data/donors"

    def donor(nm):
        return ((D / nm / "genes.txt").read_text().split(), (D / nm / "perts.txt").read_text().split(),
                np.load(D / nm / "delta.npy", mmap_mode="r"))
    t2 = pq.read_table(CLOUD / "pred_k562.parquet")
    pg2 = json.loads(t2.schema.metadata[b"pie"])["genes"]["replogle"]
    d2 = t2.to_pandas()
    kg, kp, kd = donor("k562_gwps")
    hg, hp, hd = donor("orion_hct116")
    eg, ep, ed = donor("orion_hek293t")
    panel = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))
             if r["target_gene"] != "non-targeting"]
    common = sorted(set(pg2) & set(kg) & set(hg) & set(eg) - set(panel))
    ix2 = {n: {g: i for i, g in enumerate(gl)} for n, gl in (("p", pg2), ("k", kg), ("h", hg), ("e", eg))}
    ci = {n: np.array([ix2[n][g] for g in common]) for n in ix2}
    kpi, hpi, epi = ({p: i for i, p in enumerate(pl)} for pl in (kp, hp, ep))
    P2 = {r["perturbation"]: np.asarray(r["lfc_pred"], float) for _, r in d2.iterrows()}
    tg2 = [x for x in panel if x in P2 and x in kpi]
    bck = np.load(D / "k562_gwps" / "b_ctrl.npy").astype(float)[ci["k"]]
    pck = np.expm1(bck)
    pck /= pck.sum()
    cpmk = 1e6 * pck

    def z(v):
        return (v - v.mean()) / v.std() if v.std() > 0 else v * 0
    TB, DM, EP = [], [], []
    for x in tg2:
        TB.append(np.asarray(kd[kpi[x]], float)[ci["k"]])
        h = np.asarray(hd[hpi[x]], float)[ci["h"]] if x in hpi else np.zeros(len(common))
        e = np.asarray(ed[epi[x]], float)[ci["e"]] if x in epi else np.zeros(len(common))
        DM.append((z(h) + z(e)) / max(1, (x in hpi) + (x in epi)))
        q = pck * np.exp2(np.where(cpmk >= 5.0, P2[x][ci["p"]], 0.0))
        q /= q.sum()
        EP.append(np.log1p(5e4 * q) - np.log1p(5e4 * pck))
    TB, DM, EP = map(np.stack, (TB, DM, EP))
    sm2 = np.array([x in seen for x in tg2])
    score(out, "K562", DM, EP, TB, cpmk, np.ones(len(common), bool), sm2)


def main():
    path = HERE / "work/overread702.json"
    if "--k562-only" in sys.argv:
        old = json.loads(path.read_text())
        out = {k: v for k, v in old.items() if k.startswith("h1test")}
    else:
        out = {}
        h1_part(out)
    k562_part(out)
    print(json.dumps(out, indent=1))
    path.write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
