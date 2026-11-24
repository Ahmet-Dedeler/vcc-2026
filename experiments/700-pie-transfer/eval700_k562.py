#!/usr/bin/env python3
"""700 — bulk direction on a line PIE never trained on: K562 (Replogle GWPS truth), the A/B/C panel's covered targets.

Truth: data/donors/k562_gwps (pseudobulk effect log1p(5e4 p) - control, 9,787 perts x 8,248 genes).
PIE: predictions.parquet from infer700.py for dataset `replogle`, context `k562` (PIE's own Replogle-Nadig dir supplies
     K562's control mean; the xdataset checkpoint has no Replogle training rows and no Replogle evidence).
Baselines, the donor-transfer idea in its plainest form: Orion HCT116 and HEK293T rows of the same target (data/donors/
orion_*), each z-scored per row, and their mean. Both are inputs PIE also sees as evidence, so "PIE minus donor mean" is
PIE's added value over pooling those donors.
Metrics on genes common to every source, panel targets removed (as pds does): per-target cosine to the truth, the
energy-weighted mean cosine (weights |truth|^2, the mse-relevant average), and a mini-pds (rank of the right target among
all evaluated truths).
usage: uv run python experiments/700-pie-transfer/eval700_k562.py <predictions.parquet> [out.json]
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D = ROOT / "data/donors"


def donor(name):
    genes = (D / name / "genes.txt").read_text().split()
    perts = (D / name / "perts.txt").read_text().split()
    return genes, perts, np.load(D / name / "delta.npy", mmap_mode="r")


def zrow(v):
    s = v.std()
    return (v - v.mean()) / s if s > 0 else v * 0


def main():
    panel = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
    panel = [t for t in panel if t != "non-targeting"]
    t = pq.read_table(sys.argv[1])
    meta = json.loads(t.schema.metadata[b"pie"])
    df = t.to_pandas()
    df = df[(df["dataset"] == "replogle") & (df["context"] == "k562")]
    pg = meta["genes"]["replogle"]
    P = {r["perturbation"]: (np.asarray(r["lfc_pred"], float), np.asarray(r["delta_p_pred"], float)) for _, r in df.iterrows()}
    kg, kp, kd = donor("k562_gwps")
    hg, hp, hd = donor("orion_hct116")
    eg, ep, ed = donor("orion_hek293t")
    common = sorted(set(pg) & set(kg) & set(hg) & set(eg) - set(panel))
    ix = {n: {g: i for i, g in enumerate(gl)} for n, gl in (("pie", pg), ("k", kg), ("h", hg), ("e", eg))}
    ci = {n: np.array([ix[n][g] for g in common]) for n in ix}
    kpi, hpi, epi = ({p: i for i, p in enumerate(pl)} for pl in (kp, hp, ep))
    tg = [x for x in panel if x in P and x in kpi]
    rows = {"truth": [], "pie_lfc": [], "pie_dp": [], "pie_effect": [], "hct116": [], "hek293t": [], "donor_mean": [], "pie_dp+donor": []}
    bck = np.load(D / "k562_gwps" / "b_ctrl.npy").astype(float)[ci["k"]]
    pck = np.expm1(bck); pck = pck / pck.sum()
    for x in tg:
        tr = np.asarray(kd[kpi[x]], float)[ci["k"]]
        lf, dp = P[x][0][ci["pie"]], P[x][1][ci["pie"]]
        h = np.asarray(hd[hpi[x]], float)[ci["h"]] if x in hpi else np.zeros(len(common))
        e = np.asarray(ed[epi[x]], float)[ci["e"]] if x in epi else np.zeros(len(common))
        dm = (zrow(h) + zrow(e)) / max(1, (x in hpi) + (x in epi))
        qk = pck * np.exp2(np.where(1e6 * pck >= 5.0, lf, 0.0)); qk = qk / qk.sum()
        rows["pie_effect"].append(np.log1p(5e4 * qk) - np.log1p(5e4 * pck))
        rows["truth"].append(tr); rows["pie_lfc"].append(lf); rows["pie_dp"].append(dp)
        rows["hct116"].append(h); rows["hek293t"].append(e); rows["donor_mean"].append(dm)
        rows["pie_dp+donor"].append(zrow(dp) + dm)
    R = {k: np.stack(v) for k, v in rows.items()}
    seen_set = set(json.load(open(HERE / "work/pie_train_perts.json")))
    seen_t = np.array([x in seen_set for x in tg])

    def ublend(base, pie, w, only):
        out = base.copy()
        for i in range(len(base)):
            if only is not None and not only[i]:
                continue
            nb, nv = np.linalg.norm(base[i]), np.linalg.norm(pie[i])
            if nb == 0 or nv == 0:
                continue
            m = (1 - w) * base[i] / nb + w * pie[i] / nv
            out[i] = nb * m / np.linalg.norm(m)
        return out
    for w in (0.25, 0.5, 0.75):
        R[f"blend_all_w{w}"] = ublend(R["donor_mean"], R["pie_effect"], w, None)
        R[f"blend_seen_w{w}"] = ublend(R["donor_mean"], R["pie_effect"], w, seen_t)
    tb = R["truth"]
    w = (tb ** 2).sum(1)

    def cos(a, b):
        na, nb = np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1)
        return (a * b).sum(1) / np.where(na * nb > 0, na * nb, np.inf)

    def minipds(a):
        an = a / np.maximum(np.linalg.norm(a, axis=1, keepdims=True), 1e-12)
        bn = tb / np.maximum(np.linalg.norm(tb, axis=1, keepdims=True), 1e-12)
        C = an @ bn.T
        n = len(a)
        r = [(C[i] > C[i, i]).sum() + 0.5 * ((C[i] == C[i, i]).sum() - 1) for i in range(n)]
        return float(1 - np.mean(r) / (n - 1))

    res = {"n_targets": len(tg), "n_genes": len(common), "methods": {}}
    for k in R:
        if k == "truth":
            continue
        c = cos(R[k], tb)
        res["methods"][k] = dict(mean_cos=float(np.mean(c)), energy_weighted_cos=float((c * w).sum() / w.sum()),
                                 minipds=minipds(R[k]))
    seen = set(json.load(open(HERE / "work/pie_train_perts.json")))
    sm = np.array([x in seen for x in tg])
    res["n_seen"] = int(sm.sum())
    res["by_seen"] = {}
    for lab, msk in (("seen", sm), ("unseen", ~sm)):
        if msk.sum() < 3:
            continue
        res["by_seen"][lab] = {"n": int(msk.sum())}
        for k in R:
            if k == "truth":
                continue
            c = cos(R[k][msk], tb[msk])
            ww = w[msk]
            res["by_seen"][lab][k] = dict(mean_cos=float(np.mean(c)), energy_weighted_cos=float((c * ww).sum() / ww.sum()))
    out = sys.argv[2] if len(sys.argv) > 2 else str(HERE / "work/eval700_k562.json")
    Path(out).write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
