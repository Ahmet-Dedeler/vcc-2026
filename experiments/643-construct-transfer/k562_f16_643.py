"""643: F16's construct seed off-target signature, re-measured in K562 GWPS cells carrying the SAME protospacer pairs.

Preregistration: PREREGISTRATION.md (sha256 in work/PREREGISTRATION.sha256, committed f8f95ecf before any K562 value).
    uv run python experiments/643-construct-transfer/k562_f16_643.py      -> work/k562_f16_643.json, .csv
"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
H5 = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
ENR = ROOT / "experiments/309-c9-calibration-repair/work/enrichment309.csv"
SPACERS = ROOT / "experiments/301-metadata-acquisitions/work/vcc2025_spacers.csv"
N_CTRL = 10000
N_NULL = 10
SEED = 643


def rss():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30


def log(m):
    print(f"{time.strftime('%H:%M:%S')} [rss {rss():.2f} GiB] {m}", flush=True)


def read_rows_once(path: Path, idx: np.ndarray) -> dict[int, np.ndarray]:
    """Rows `idx` of the GWPS X (dense, gzip-chunked (3886, 33)) in ONE sorted pass: each row block is read as a
    contiguous slice, so every chunk is decompressed once, not once per scattered-row request."""
    import h5py
    idx = np.unique(idx)
    out = {}
    with h5py.File(path, "r") as f:
        X = f["X"]
        rb = X.chunks[0]
        blocks = np.unique(idx // rb)
        for k, b in enumerate(blocks):
            lo, hi = int(b * rb), int(min((b + 1) * rb, X.shape[0]))
            want = idx[(idx >= lo) & (idx < hi)]
            blk = X[lo:hi]
            for i in want:
                out[int(i)] = np.asarray(blk[i - lo], dtype=np.float64)
            if k % 50 == 0:
                log(f"  block {k + 1}/{len(blocks)} rows {len(out)}/{idx.size}")
    return out


def mean_log1p_cp50k(X: np.ndarray) -> np.ndarray:
    assert np.all(X[: min(50, len(X))] == np.floor(X[: min(50, len(X))])), "X is not integer counts"
    tot = X.sum(axis=1)
    return np.log1p(X * (5e4 / np.maximum(tot, 1))[:, None]).mean(axis=0)


def main():
    enr = [r for r in csv.DictReader(open(ENR)) if r["unit"] == "targeting"]
    log(f"309 constructs (targeting): {len(enr)}")
    sp = {r["guide_id"].split("|")[0].rsplit("_", 1)[0]: r for r in csv.DictReader(open(SPACERS))}
    a = ad.read_h5ad(H5, backed="r")
    obs = a.obs
    pert = obs["perturbation"].astype(str).to_numpy()
    gid = obs["guide_id"].astype(str).to_numpy()
    ctrl_labels = [x for x in pd.unique(pert) if "control" in x.lower() or "non-targeting" in x.lower()]
    log(f"GWPS cells {a.n_obs}, genes {a.n_vars}; control labels {ctrl_labels}")
    genes = a.var_names.astype(str).to_numpy()
    gpos = {g: i for i, g in enumerate(genes)}
    pairset = {}
    for g in pd.unique(gid):
        pairset.setdefault(frozenset(g.split("|")), []).append(g)
    rows, cons = [], {}
    for r in enr:
        c = r["construct"]
        s = sp.get(c)
        key = frozenset({s["sgid_A"], s["sgid_B"]}) if s and s["sgid_A"] and s["sgid_B"] and s["route"] != "UNRESOLVED" else None
        gw = pairset.get(key, []) if key else []
        cons[c] = dict(target=r["intended_target"], candidates=[x for x in r["candidates"].split(",") if x], gw=gw,
                       h1_delta_own=float(r["delta_own"]), h1_kd=float(r["kd_lfc"]), h1_cells=int(r["cells"]),
                       h1_delta_repaired=float(r["delta_repaired"]))
    have = {c: v for c, v in cons.items() if v["gw"]}
    log(f"constructs found in GWPS by sgID pair: {len(have)}/{len(cons)}")
    # cells
    rng = np.random.default_rng(SEED)
    ctrl_idx = np.flatnonzero(np.isin(pert, ctrl_labels))
    ctrl_idx = np.sort(rng.choice(ctrl_idx, size=min(N_CTRL, ctrl_idx.size), replace=False))
    cells_of = {c: np.flatnonzero(np.isin(gid, v["gw"])) for c, v in have.items()}
    rows = read_rows_once(H5, np.concatenate([ctrl_idx, *cells_of.values()]))
    Xc = np.vstack([rows[int(i)] for i in ctrl_idx])
    ctrl = mean_log1p_cp50k(Xc)
    ctrl_cp = None
    log(f"controls: {ctrl_idx.size} cells")
    for c, v in have.items():
        idx = cells_of[c]
        v["k562_cells"] = int(idx.size)
        v["delta"] = mean_log1p_cp50k(np.vstack([rows[int(i)] for i in idx])) - ctrl
    # control CP50k mean per gene (expression filter): from the same controls
    tot = Xc.sum(axis=1)
    ctrl_cp = (Xc * (5e4 / np.maximum(tot, 1))[:, None]).mean(axis=0)
    del Xc
    dec = np.digitize(np.log1p(ctrl_cp), np.quantile(np.log1p(ctrl_cp[ctrl_cp > 0]), np.linspace(0, 1, 11)[1:-1]))
    out = []
    for c, v in have.items():
        partners = [p for p, w in have.items() if w["target"] == v["target"] and p != c]
        if len(partners) != 1:
            continue
        p = partners[0]
        S = [gpos[g] for g in v["candidates"] if g in gpos and ctrl_cp[gpos[g]] >= 1.0]
        if not S:
            out.append(dict(construct=c, partner=p, target=v["target"], n_cand=0))
            continue
        tgt = gpos.get(v["target"])
        own, par = float(np.mean(v["delta"][S])), float(np.mean(have[p]["delta"][S]))
        # matched null: random gene sets of the same size, matched on control-expression decile, excluding the target
        nulls = []
        for _ in range(N_NULL):
            pick = []
            for s in S:
                pool = np.flatnonzero((dec == dec[s]) & (ctrl_cp >= 1.0))
                pool = pool[pool != (tgt if tgt is not None else -1)]
                pick.append(int(rng.choice(pool)))
            nulls.append(float(np.mean(v["delta"][pick]) - np.mean(have[p]["delta"][pick])))
        out.append(dict(construct=c, partner=p, target=v["target"], n_cand=len(S), candidates=[genes[s] for s in S],
                        k562_cells=v["k562_cells"], partner_cells=have[p]["k562_cells"], own=own, partner_val=par,
                        D=own - par, null_D_mean=float(np.mean(nulls)), null_D=nulls,
                        k562_kd_own=float(v["delta"][tgt]) if tgt is not None else None,
                        h1_delta_own=v["h1_delta_own"], h1_delta_repaired=v["h1_delta_repaired"], h1_kd=v["h1_kd"]))
    ok = [r for r in out if r.get("n_cand", 0) > 0]
    D = np.array([r["D"] for r in ok])
    N = np.array([r["null_D_mean"] for r in ok])
    res = dict(n_constructs_309=len(cons), n_in_gwps=len(have), n_tested=len(ok), median_D=float(np.median(D)) if D.size else None,
               mean_D=float(np.mean(D)) if D.size else None, frac_negative=float(np.mean(D < 0)) if D.size else None,
               wilcoxon_one_sided_p=float(wilcoxon(D, alternative="less").pvalue) if D.size >= 6 else None,
               null_median=float(np.median(N)) if N.size else None,
               null_wilcoxon_p=float(wilcoxon(N, alternative="less").pvalue) if N.size >= 6 else None)
    if D.size >= 6:
        h1 = np.array([r["h1_delta_repaired"] for r in ok])
        res["corr_D_k562_vs_h1_delta_repaired"] = float(np.corrcoef(D, h1)[0, 1])
        boot = [np.median(D[rng.integers(0, D.size, D.size)]) for _ in range(5000)]
        res["median_D_ci95"] = [float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))]
    verdict = ("UNRESOLVABLE (matched null median below -0.01)" if res["null_median"] is not None and res["null_median"] < -0.01
               else "REPLICATES in K562" if res["wilcoxon_one_sided_p"] is not None and res["median_D"] < 0 and res["wilcoxon_one_sided_p"] < 0.05
               else "DOES NOT REPLICATE")
    res["verdict"] = verdict
    log(json.dumps(res, indent=1))
    (HERE / "work/k562_f16_643.json").write_text(json.dumps(dict(result=res, rows=out), indent=1, default=str))
    pd.DataFrame([{k: v for k, v in r.items() if k not in ("null_D", "candidates")} for r in out]).to_csv(HERE / "work/k562_f16_643.csv", index=False)


if __name__ == "__main__":
    main()
