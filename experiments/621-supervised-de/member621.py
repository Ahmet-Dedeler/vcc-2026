#!/usr/bin/env python3
"""621b — does the learned membership head pick more truly-DE genes than production's call set?

h1 (the only surface with production's call classes), non-essential targets. For each target, rank tested genes by:
production's composition (head = top-400 |key| first, then the CPM block by ref_lin, then the rest), |key| alone, the
recipient's control log-CPM alone, |Zbal| alone, and the 621 HGB membership score. Report pooled AUROC and the mean
per-target precision (phi) and recall of the gate at K = 400 / 1,000 / 2,600 calls.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
from inst636 import ESS, Inst  # noqa: E402
import lib521 as L  # noqa: E402


def auroc(score, y):
    from scipy.stats import rankdata
    r = rankdata(score)
    n1 = y.sum(); n0 = len(y) - n1
    return float((r[y].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def main() -> None:
    I = Inst(); D = I.S["h1"]; T, G = D["T"], D["G"]
    z = np.load(HERE / "work/pred-h1.npz", allow_pickle=True)
    member, Zb = z["member"].astype(float), z["Zbal"].astype(float)
    S = L.load_state(0)
    ti = {t: i for i, t in enumerate(S["targets"])}; gi = {g: j for j, g in enumerate(S["genes"])}
    rt = np.array([ti[t] for t in T]); cg = np.array([gi.get(g, -1) for g in G]); ok = cg >= 0

    def remap(A, fill):
        out = np.full(D["lfc"].shape, fill, dtype=float)
        out[:, ok] = np.asarray(A, float)[rt][:, cg[ok]]
        return out
    absk = remap(S["absk"], 0.0)
    head = remap(S["head"], 0.0) > 0
    block = remap(S["block"], 0.0) > 0
    reflin = np.broadcast_to(remap(np.broadcast_to(S["ref_lin"], S["absk"].shape), -np.inf), D["lfc"].shape)
    lcpm = np.broadcast_to(np.log10(np.load(ROOT / "experiments/480-nmae-instrument/work/gate-h1_full_ka.npz")["ctrl_cpm"] + 1e-3), D["lfc"].shape)
    # production composition score: head first (by |key|), then block (by ref_lin), then the rest (by ref_lin)
    prod = np.where(head, 2e6 + absk, np.where(block, 1e6 + reflin, reflin))
    noness = np.array([t not in ESS for t in T])
    tested = D["tested"] & ok[None, :] & noness[:, None]
    y_all = D["gate"]
    rankers = {"production composition": prod, "|key|": absk, "control log-CPM": lcpm,
               "|Zbal|": np.nan_to_num(np.abs(Zb), nan=-1.0), "621 HGB membership": np.nan_to_num(member, nan=-1.0),
               "HGB x lcpm rank-mix": None}
    # rank-mix: average of per-target ranks of HGB membership and control log-CPM
    out = {}
    for name, sc in rankers.items():
        if sc is None:
            from scipy.stats import rankdata
            sc = np.zeros_like(prod)
            for i in range(len(T)):
                sc[i] = rankdata(np.nan_to_num(member[i], nan=-1)) + rankdata(lcpm[i])
        m = tested
        au = auroc(sc[m], y_all[m])
        res = {"auroc_pooled": au}
        for K in (400, 1000, 2600):
            prec, rec = [], []
            for i in np.flatnonzero(noness):
                cand = np.flatnonzero(tested[i])
                if len(cand) == 0 or y_all[i, cand].sum() == 0:
                    continue
                top = cand[np.argsort(-sc[i, cand], kind="stable")[:K]]
                hit = y_all[i, top].sum()
                prec.append(hit / len(top)); rec.append(hit / y_all[i, cand].sum())
            res[f"phi@{K}"] = float(np.mean(prec)); res[f"recall@{K}"] = float(np.mean(rec))
        out[name] = res
        print(f"{name:26s} AUROC {au:.4f} | " + "  ".join(f"phi@{K} {res[f'phi@{K}']:.4f} rec@{K} {res[f'recall@{K}']:.3f}" for K in (400, 1000, 2600)))
    (HERE / "work/member621.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
