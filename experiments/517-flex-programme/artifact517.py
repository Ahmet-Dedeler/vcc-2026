#!/usr/bin/env python3
"""517 — the hook-source artifact with the better programme ESTIMATOR, for the A/B/C panel (a `file:` source for 517's
builder patch; same format as 458d's `arch1d_programme_abc.npz`). Not a block prior (w 0, no rows).

  meanz_ne = mean over Arc's NON-ESSENTIAL rows (K562-essential list out) of delta_g / sd_g (sd over all the programme's
             rows), the A/B/C panel's own rows excluded (396's rule).
Why: variant517 (H1, analytic, validated path): the gated hook from meanz_ne gives nmae −0.0183 vs −0.0148 for P_down over
all rows (the production rule), reach 0 vs −0.003; progest517 (Xaira, cross-line within a lab, Arc-sized 275-row
programmes): the essential rows do not transfer to non-essential targets' uncalled genes, and mean-z over the non-essential
rows is the best small-programme estimator (+0.01 … +0.03 sign accuracy over the P_down mix).
Encoding: pdown_<c> = clip(0.5 − 0.05 · z / z_q), z_q chosen so that |pdown − 0.5| ≥ 0.05 covers exactly as many official-axis
genes as the production programme (P_down over all rows) does at tau 0.05; so `file:<this>:<sha>:0.05` has the
production hook's coverage. programme_<c> = its sign (−1 where z ≤ 0 or uncovered), recorded but not meant for the block.
KA: with --estimator pdown_all the script rebuilds 458d's pdown_A/B/C exactly.
usage: uv run python experiments/517-flex-programme/artifact517.py [--estimator meanz_ne|pdown_all] [--out work/<name>.npz]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BANK = ROOT / "data/donors/h1_2025_full"
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}
REF458D = ROOT / "experiments/458-arc-lab-programme/work/arch1d_programme_abc.npz"


def estimate(delta, perts, genes, official, panel, estimator):
    keep = np.array([p not in panel for p in perts])
    X = np.asarray(delta[np.flatnonzero(keep)], dtype=np.float64)
    ne = np.array([p not in ESS for p in np.asarray(perts)[keep]])
    nz = (X != 0).sum(0); dn = (X < 0).sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        pd_all = np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    pos = {g: j for j, g in enumerate(genes)}
    col = np.array([pos.get(g, -1) for g in official])
    on = lambda v: np.where(col >= 0, v[np.maximum(col, 0)], np.nan)   # noqa: E731
    P_all = on(pd_all)
    info = {"rows": int(keep.sum()), "rows_nonessential": int(ne.sum())}
    if estimator == "pdown_all":
        return P_all, info
    sd = X.std(0); sd[sd == 0] = 1.0
    z = on((X[ne] / sd).mean(0))
    z = np.where(np.isfinite(P_all), z, np.nan)                         # the same coverage as the production programme
    n_cov = int((np.isfinite(P_all) & (np.abs(P_all - 0.5) >= 0.05)).sum())
    a = np.abs(np.nan_to_num(z))
    z_q = float(np.sort(a)[::-1][n_cov - 1])
    P = np.where(np.isfinite(z), np.clip(0.5 - 0.05 * z / z_q, 0.001, 0.999), np.nan)
    info.update(z_q=z_q, n_cov_tau05=n_cov, n_cov_check=int((np.abs(np.nan_to_num(P, nan=0.5) - 0.5) >= 0.05).sum()))
    return P, info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--estimator", choices=["meanz_ne", "pdown_all"], default="meanz_ne")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    perts_panel = [r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
    delta = np.load(BANK / "delta.npy", mmap_mode="r")
    bp = [l.strip() for l in open(BANK / "perts.txt")]
    bg = [l.strip() for l in open(BANK / "genes.txt")]
    P, info = estimate(delta, bp, bg, official, set(perts_panel), a.estimator)
    if a.estimator == "pdown_all":
        with np.load(REF458D, allow_pickle=False) as z:
            ok = all(np.array_equal(P, z[f"pdown_{c}"], equal_nan=True) for c in "ABC")
        print(f"KA: pdown_all == 458d's pdown_A/B/C: {'PASS' if ok else 'FAIL'}", flush=True)
        if not ok:
            raise SystemExit(1)
        return
    sign = np.where(np.isfinite(P) & (P <= 0.5), 1.0, -1.0)
    out = {"genes": np.asarray(official), "perts": np.asarray(perts_panel), "contexts": np.asarray(list("ABC"))}
    for c in "ABC":
        out.update({f"programme_{c}": sign, f"programme_ok_{c}": np.asarray(True), f"w_{c}": np.asarray(0.0), f"pdown_{c}": P,
                    f"targets_{c}": np.zeros(0, np.int64), f"rows_{c}": np.zeros((0, len(official)), np.float32),
                    f"signs_{c}": np.zeros((0, len(official)), np.int8)})
    bank_sha = hashlib.sha256((BANK / "delta.npy").read_bytes()).hexdigest()
    out["provenance"] = np.asarray(json.dumps({"experiment": "517", "estimator": a.estimator, "bank": "h1_2025_full",
                                               "bank_delta_sha256": bank_sha, **info,
                                               "note": "hook source (VCC_DE_MDE_SIGN=file:...); not a block prior"}))
    path = Path(a.out) if a.out else HERE / f"work/arch1{a.estimator}_hook_abc.npz"
    np.savez_compressed(path, **out)
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    Path(str(path) + ".json").write_text(json.dumps({"sha256": sha, **info, "estimator": a.estimator}, indent=1))
    print(f"wrote {path} sha256 {sha}; {json.dumps(info)}", flush=True)
    # coverage on each context's gated genes at tau .05, vs the production programme's
    with np.load(REF458D, allow_pickle=False) as z:
        Pp = z["pdown_A"]
    for c in "ABC":
        m = np.load(ROOT / f"experiments/491-own-control-programme/work/ctrl_meancpm_{c}.npy") > 5
        n_new = int((np.abs(np.nan_to_num(P[m], nan=0.5) - 0.5) >= 0.05).sum())
        n_old = int((np.abs(np.nan_to_num(Pp[m], nan=0.5) - 0.5) >= 0.05).sum())
        agree = float(np.mean(np.sign(0.5 - P[m]) == np.sign(0.5 - Pp[m])))
        print(f"  {c}: gated genes {int(m.sum())}; tau .05 signs {n_new} (production programme {n_old}); sign agreement with 458d's "
              f"programme {agree:.3f}", flush=True)


if __name__ == "__main__":
    main()
