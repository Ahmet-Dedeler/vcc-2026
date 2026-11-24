"""481: is the efr synthetic panel a fair stand-in for an Arc-built panel? A donor-side strength proxy, calibrated on the
two panels where the builder's own design strength is known (A/B/C relabelled: 453 dump; efr: 465 dump), then applied to
Arc's other 300-target panel (the 2025 VCC H1 panel, essential targets removed = the essential-free stand-in).

Proxy per target t: the number of context-D render-gated genes (pooled control CPM > 10, the builder's gate) where the
three production donors that carry t agree in sign and their mean raw delta has |mean| > tau. Raw bank rows (log1p
pseudobulk shifts at TS 5e4, `data/donors/<d>/delta.npy`), no whitening: a proxy, judged only by how well it ranks and
reproduces the builder's strength on the two known panels.

usage: uv run python proxy481.py -> work/proxy481.json (peak < 1.5 GB: bank rows are read by target through mmap)
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
PANELS = ["abc", "h1_2025", "efr", "efp", "real199", "pess199"]
TAUS = (0.1, 0.15, 0.2, 0.3)
ESS = {ln.strip() for ln in (ROOT / "data/donors/k562_essential/perts.txt").read_text().split("\n") if ln.strip()}
OFFICIAL = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
OPOS = {g: i for i, g in enumerate(OFFICIAL)}


def gate_D() -> np.ndarray:
    z = np.load(HERE / "work/ctrl_def_D.npz", allow_pickle=False)
    tot = z["totals"]
    cpm = 1e6 * tot / tot.sum()
    return (cpm > 10.0) & (np.log1p(5e4 * tot / tot.sum()) > 0)


def banks():
    out = {}
    for d in DONORS:
        dd = ROOT / "data/donors" / d
        delta = np.load(dd / "delta.npy", mmap_mode="r")
        perts = [ln.strip() for ln in (dd / "perts.txt").read_text().split("\n") if ln.strip()]
        genes = [ln.strip() for ln in (dd / "genes.txt").read_text().split("\n") if ln.strip()]
        assert delta.shape == (len(perts), len(genes)), d
        col = np.array([OPOS.get(g, -1) for g in genes])
        out[d] = (delta, {p: i for i, p in enumerate(perts)}, col)
    return out


def proxy(panel: list[str], B, gate) -> dict[str, np.ndarray]:
    pset = set(panel)
    g_ok = gate & np.array([g not in pset for g in OFFICIAL])
    res = {tau: np.zeros(len(panel)) for tau in TAUS}
    ncar = np.zeros(len(panel), dtype=int)
    for i, t in enumerate(panel):
        s = np.zeros(len(OFFICIAL)); n = np.zeros(len(OFFICIAL)); pos = np.zeros(len(OFFICIAL)); neg = np.zeros(len(OFFICIAL))
        for d, (delta, pidx, col) in B.items():
            if t not in pidx:
                continue
            ncar[i] += 1
            r = np.asarray(delta[pidx[t]], dtype=np.float64)
            m = col >= 0
            s[col[m]] += r[m]; n[col[m]] += 1
            pos[col[m]] += r[m] > 0; neg[col[m]] += r[m] < 0
        have = (n > 0) & g_ok
        mean = np.where(have, s / np.maximum(n, 1), 0.0)
        agree = have & ((pos == n) | (neg == n))
        for tau in TAUS:
            res[tau][i] = float((agree & (np.abs(mean) > tau)).sum())
    res["ncar"] = ncar
    return res


def dump_strength(dump: Path, panel: list[str], gate) -> np.ndarray:
    z = np.load(dump, allow_pickle=False)
    assert [str(p) for p in z["perts"]] == panel, "dump panel order"
    car = np.asarray(z["carriers"], dtype=np.float64)
    return (np.abs(np.where(gate[None, :], car, 0.0)) > 0.2).sum(1).astype(float)


def main() -> None:
    gate = gate_D()
    B = banks()
    P = {n: [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{n}.csv")["target_gene"]] for n in PANELS}
    prox = {n: proxy(P[n], B, gate) for n in PANELS}
    known = {"abc": dump_strength(ROOT / "experiments/453-abc-dress-rehearsal/work/run/dump_h1x2_cd4x1.npz", P["abc"], gate),
             "efr": dump_strength(ROOT / "experiments/465-both-switches-efr/work/run/dump_h1x2_cd4x1.npz", P["efr"], gate)}
    out = {"calibration": {}, "panels": {}}
    for tau in TAUS:
        rec = {}
        for n, s in known.items():
            rho = spearmanr(prox[n][tau], s).statistic
            rec[n] = {"spearman": round(float(rho), 3)}
        # quantile map fitted on A/B/C: proxy -> builder strength (monotone, by matching ranks)
        xa, ya = np.sort(prox["abc"][tau]), np.sort(known["abc"])
        def qmap(x):
            return np.interp(x, xa, ya, left=ya[0], right=ya[-1])
        pred_efr = np.round(qmap(prox["efr"][tau]))
        rec["efr_predicted_frac_le1_from_abc_map"] = round(float((pred_efr <= 1).mean()), 3)
        rec["efr_actual_frac_le1"] = round(float((known["efr"] <= 1).mean()), 3)
        out["calibration"][str(tau)] = rec
        for n in PANELS:
            x = prox[n][tau]
            ne = np.array([t not in ESS for t in P[n]])
            pred = np.round(qmap(x))
            out["panels"].setdefault(n, {})[str(tau)] = {
                "proxy_quartiles": np.quantile(x, [0, .25, .5, .75, 1]).tolist(),
                "proxy_nonessential_quartiles": np.quantile(x[ne], [0, .25, .5, .75, 1]).tolist() if ne.any() else None,
                "pred_frac_le1": round(float((pred <= 1).mean()), 3),
                "pred_frac_le1_nonessential": round(float((pred[ne] <= 1).mean()), 3) if ne.any() else None,
                "pred_q33_threshold_nonessential": float(np.quantile(pred[ne], 0.33)) if ne.any() else None,
                "n_nonessential": int(ne.sum())}
    for n, s in known.items():
        out["panels"][n]["builder_strength_quartiles"] = np.quantile(s, [0, .25, .5, .75, 1]).tolist()
        out["panels"][n]["builder_frac_le1"] = round(float((s <= 1).mean()), 3)
    (HERE / "work/proxy481.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out["calibration"], indent=1))
    for n in PANELS:
        r = out["panels"][n]
        print(n, {tau: (r[tau]["proxy_nonessential_quartiles"], r[tau]["pred_frac_le1_nonessential"],
                        r[tau]["pred_q33_threshold_nonessential"]) for tau in ("0.15", "0.2")},
              r.get("builder_frac_le1"))


if __name__ == "__main__":
    main()
