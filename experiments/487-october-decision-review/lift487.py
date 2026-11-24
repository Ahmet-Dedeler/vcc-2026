"""487 decision 1c: how much direction does a same-line bank add, in the October-relevant cases, on ONE statistic?

450's +0.0934 was priced on H1, where the same-line bank (Pan 2025 hESC) is NEW information: none of our donors is an
ESC line. The October lines with high bank coverage (K562: GWPS; HCT116: Orion; 451's coverage table: 300/300) are
production donors ALREADY, so the blend adds that bank a second time. 450 §6 measured that case only on the essential
stratum (+15 % centred cosine). The 2026 panel is essential-free.

Everything here uses one construction and one statistic:
  rows: each bank's delta rows on the scenario's common genes, each target's own gene zeroed, centred by the bank's mean
        over the evaluated targets;
  mix3 = unit-mean of the three production donors (K562-GWPS, Orion HCT116, Orion HEK293T), unwhitened;
  blend(w) = unit(mix3) + w * unit(same-line row);
  statistic: centred cosine to the truth row, mean over targets (also the truth's top-100 genes' sign agreement).
Scenarios:
  S1 H1 (calibration: the bank is new information): truth = Arc's 2025 H1 screen (h1_2025_full, all cells),
     same = Pan 2025 hESC rows, clean and with a random fraction p of measured genes sign-flipped (450f's weakening,
     p 0.2 / 0.3), whose emitter-level prices exist (+0.0934 / +0.0268 / +0.0007 at w 2).
  S2 K562 (bank already a donor): truth = VIPerturb K562 (Satija lab), same = K562-GWPS.
  S3 HEK293 (bank already a donor): truth = XuCao 2023 HEK293, same = Orion HEK293T.
Ratios cos(blend)/cos(mix3) are invariant to independent truth noise (it scales every predictor's cosine alike), so
they compare across truths of different depth. Target bootstrap (2,000) on the ratio. Peak RSS ~1-2 GB (mmap rows).
usage: uv run python experiments/487-october-decision-review/lift487.py -> work/lift487.json
"""
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
D = ROOT / "data/donors"
T225 = ROOT / "experiments/225-source-sweep/work/tables"
ESS = set((D / "k562_essential/perts.txt").read_text().split())
WS = (0.5, 1, 2, 4)


def bank(path):
    return (np.load(path / "delta.npy", mmap_mode="r"), [l.strip() for l in open(path / "perts.txt")],
            [l.strip() for l in open(path / "genes.txt")])


def _cat(ds):
    return ds["categories"][:].astype(str)[ds["codes"][:]] if isinstance(ds, h5py.Group) else ds[:].astype(str)


def pan_bank(targets):
    db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
    with h5py.File(db, "r") as f:
        ds, pg, vg = _cat(f["obs/dataset"]), _cat(f["obs/perturbed_gene_name"]), _cat(f["var/gene_name"])
        sel = np.sort(np.where((ds == "Pan2025_hESC") & np.isin(pg, list(targets)))[0])
        X = np.where(f["layers/measure_mask"][sel, :] > 0, f["X"][sel, :], np.nan)
        pt = pg[sel]
    ts = sorted(set(pt))
    R = np.full((len(ts), len(vg)), np.nan)
    for i, t in enumerate(ts):
        with np.errstate(invalid="ignore"):
            R[i] = np.nanmean(X[pt == t], axis=0)
    return R, ts, list(vg)


def pan_targets():
    db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
    with h5py.File(db, "r") as f:
        ds, pg = _cat(f["obs/dataset"]), _cat(f["obs/perturbed_gene_name"])
    return set(pg[ds == "Pan2025_hESC"])


def matrix(b, ts, genes):
    delta, perts, g = b
    pi = {p: i for i, p in enumerate(perts)}
    gi = {x: i for i, x in enumerate(g)}
    rows = np.array([pi[t] for t in ts])
    cols = np.array([gi[x] for x in genes])
    order = np.argsort(rows)
    X = np.empty((len(ts), len(cols)))
    X[order] = np.asarray(delta[rows[order]], dtype=np.float64)[:, cols]
    return X


def prep(M, ts, genes):
    own = {x: j for j, x in enumerate(genes)}
    for X in M.values():
        X[~np.isfinite(X)] = 0.0
        for r, t in enumerate(ts):
            if t in own:
                X[r, own[t]] = 0.0
    return {k: X - X.mean(0) for k, X in M.items()}


u = lambda A: A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-12)
cos = lambda A, B: np.sum(u(A) * u(B), 1)


def flip(X, p, seed):
    rng = np.random.default_rng(20260926 + seed)
    Y = X.copy()
    for r in range(Y.shape[0]):
        m = np.flatnonzero(Y[r] != 0)
        f = rng.random(m.size) < p
        Y[r, m[f]] = -Y[r, m[f]]
    return Y


def evaluate(T, mix, same, mask, rng, nb=2000):
    """Per-target cosines of mix and blends; ratio of means with a target bootstrap."""
    c0 = cos(mix, T)
    res = {"n": int(mask.sum()), "mix3": round(float(c0[mask].mean()), 4), "same": round(float(cos(same, T)[mask].mean()), 4)}
    idx = np.flatnonzero(mask)
    for w in WS:
        cw = cos(u(mix) + w * u(same), T)
        ratio = cw[idx].mean() / c0[idx].mean()
        bs = []
        for _ in range(nb):
            j = idx[rng.integers(0, len(idx), len(idx))]
            bs.append(cw[j].mean() / c0[j].mean())
        res[f"w{w:g}"] = {"cos": round(float(cw[idx].mean()), 4), "ratio": round(float(ratio), 3),
                          "ratio_ci95": [round(float(np.quantile(bs, .025)), 3), round(float(np.quantile(bs, .975)), 3)],
                          "dcos": round(float(cw[idx].mean() - c0[idx].mean()), 4)}
    return res


def strength_masks(T, ts):
    e = np.linalg.norm(T, axis=1)
    ess = np.array([t in ESS for t in ts])
    q50 = np.quantile(e[~ess], 0.5) if (~ess).sum() else np.inf
    q75 = np.quantile(e[~ess], 0.75) if (~ess).sum() else np.inf
    return {"all": np.ones(len(ts), bool), "essential": ess, "non-essential": ~ess,
            "non-essential, truth energy top half": (~ess) & (e >= q50), "non-essential, truth energy top quartile": (~ess) & (e >= q75)}


def main():
    rng = np.random.default_rng(487)
    out = {}
    donors = {"k562_gwps": bank(D / "k562_gwps"), "hct116": bank(D / "orion_hct116"), "hek293t": bank(D / "orion_hek293t")}
    # ---- S1 H1: truth Arc H1 (all cells), same = Pan (new information)
    h1 = bank(D / "h1_2025_full")
    pts = pan_targets()
    ts = set(h1[1]) & pts
    for b in donors.values():
        ts &= set(b[1])
    ts = sorted(ts)
    P, pts_sorted, pg = pan_bank(ts)
    assert pts_sorted == ts
    genes = set(h1[2]) & set(pg)
    for b in donors.values():
        genes &= set(b[2])
    genes = sorted(genes)
    pgi = {x: i for i, x in enumerate(pg)}
    Pm = P[:, [pgi[x] for x in genes]]
    Pm = np.where(np.isfinite(Pm), Pm, 0.0)
    raw = {"truth": matrix(h1, ts, genes), **{k: matrix(b, ts, genes) for k, b in donors.items()}}
    M = prep({k: v.copy() for k, v in raw.items()}, ts, genes)
    mix = u(u(M["k562_gwps"]) + u(M["hct116"]) + u(M["hek293t"]))
    masks = strength_masks(M["truth"], ts)
    s1 = {"targets": len(ts), "genes": len(genes)}
    for p in (0.0, 0.2, 0.3):
        # 450f's weakening: flip the sign of a random fraction p of each RAW row's measured genes, then centre
        same = prep({"same": Pm.copy() if p == 0 else flip(Pm, p, 0)}, ts, genes)["same"]
        s1[f"flip{p:g}"] = {k: evaluate(M["truth"], mix, same, m, rng) for k, m in masks.items() if m.sum() >= 10}
    out["S1_H1_pan_new_information"] = s1
    print(f"S1 H1: {len(ts)} targets, {len(genes)} genes")
    for p in ("flip0", "flip0.2", "flip0.3"):
        for k, v in s1[p].items():
            print(f"  {p:8s} {k:42s} n {v['n']:4d} mix3 {v['mix3']:.4f} same {v['same']:.4f} "
                  + " ".join(f"w{w:g} x{v[f'w{w:g}']['ratio']:.2f}" for w in WS), flush=True)
    # ---- S2 K562: truth VIPerturb, same = GWPS (already a donor)
    vip = bank(D / "viperturb_k562")
    ts = set(vip[1])
    for b in donors.values():
        ts &= set(b[1])
    ts = sorted(ts)
    ts = sorted(np.random.default_rng(20260925).choice(ts, size=min(3000, len(ts)), replace=False).tolist())
    genes = set(vip[2])
    for b in donors.values():
        genes &= set(b[2])
    genes = sorted(genes)
    M = prep({"truth": matrix(vip, ts, genes), **{k: matrix(b, ts, genes) for k, b in donors.items()}}, ts, genes)
    mix = u(u(M["k562_gwps"]) + u(M["hct116"]) + u(M["hek293t"]))
    xline = u(u(M["hct116"]) + u(M["hek293t"]))
    masks = strength_masks(M["truth"], ts)
    s2 = {"targets": len(ts), "genes": len(genes),
          "in_mix": {k: evaluate(M["truth"], mix, M["k562_gwps"], m, rng) for k, m in masks.items()},
          "new_information (xline + same)": {k: evaluate(M["truth"], xline, M["k562_gwps"], m, rng) for k, m in masks.items()}}
    out["S2_K562_gwps_in_mix_truth_viperturb"] = s2
    print(f"S2 K562: {len(ts)} targets, {len(genes)} genes")
    for case in ("in_mix", "new_information (xline + same)"):
        for k, v in s2[case].items():
            print(f"  {case[:12]:12s} {k:42s} n {v['n']:4d} base {v['mix3']:.4f} same {v['same']:.4f} "
                  + " ".join(f"w{w:g} x{v[f'w{w:g}']['ratio']:.2f}" for w in WS), flush=True)
    # ---- S3 HEK293: truth XuCao, same = Orion HEK293T
    xu = bank(T225 / "xucao_hek293")
    ts = set(xu[1])
    for b in donors.values():
        ts &= set(b[1])
    ts = sorted(ts)
    genes = set(xu[2])
    for b in donors.values():
        genes &= set(b[2])
    genes = sorted(genes)
    if len(ts) >= 20:
        M = prep({"truth": matrix(xu, ts, genes), **{k: matrix(b, ts, genes) for k, b in donors.items()}}, ts, genes)
        mix = u(u(M["k562_gwps"]) + u(M["hct116"]) + u(M["hek293t"]))
        xline = u(u(M["k562_gwps"]) + u(M["hct116"]))
        masks = strength_masks(M["truth"], ts)
        s3 = {"targets": len(ts), "genes": len(genes),
              "in_mix": {k: evaluate(M["truth"], mix, M["hek293t"], m, rng) for k, m in masks.items() if m.sum() >= 10},
              "new_information (xline + same)": {k: evaluate(M["truth"], xline, M["hek293t"], m, rng) for k, m in masks.items() if m.sum() >= 10}}
        out["S3_HEK293_orion_in_mix_truth_xucao"] = s3
        print(f"S3 HEK293: {len(ts)} targets, {len(genes)} genes")
        for case in ("in_mix", "new_information (xline + same)"):
            for k, v in s3[case].items():
                print(f"  {case[:12]:12s} {k:42s} n {v['n']:4d} base {v['mix3']:.4f} same {v['same']:.4f} "
                      + " ".join(f"w{w:g} x{v[f'w{w:g}']['ratio']:.2f}" for w in WS), flush=True)
    else:
        print(f"S3 HEK293: only {len(ts)} common targets, skipped")
    (HERE / "work/lift487.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
