"""282 known-answer rows.  KA-1 (the across-donor ceiling reproduces 280) must pass before any 282 number.

    PYTHONDONTWRITEBYTECODE=1 uv run python ka282.py
"""
from __future__ import annotations

import json
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from lib282 import (CONDS, DORDER, PAIRINGS, ROOT, W280, Bank, boot_median, build_index,  # noqa: E402
                    cosines, panel_set, rows_for)

P = lambda *a: print(*a, flush=True)  # noqa: E731
t0 = time.time()
OUT: dict = {}
bank = Bank(W280)
SC, CC = bank.scored, bank.cols_cd4
P(f"bank {bank.X.shape}, mapped {len(CC)}, scored {SC.sum()}  [{time.time()-t0:.0f}s]")

# ===================================================================== KA-1
# 280 R280 §2.1, keep_for_DE, per-condition target sets, 3 disjoint pairings.  Same code path as q1_280.py.
# full precision, straight out of 280's own table (not R280's 4-d.p. print)
_h = pd.read_csv(f"{W280}/q1_halves_keepDE.csv")
_h["key"] = _h.pairing.str.replace(" ", "").str.replace("+", "").str.replace("(3pairings)", "", regex=False)
R280 = {c: {r.key: (r.median, r.lo, r.hi, r.n) for r in g.itertuples()} for c, g in _h.groupby("cond")}
R280_CI = {"Stim48hr": (R280["Stim48hr"]["POOLED"][1], R280["Stim48hr"]["POOLED"][2])}

o = bank.obs
tg = o[(o.gtype == "targeting") & o.keep_for_DE]
IDX280 = {(c, d, t): g.sort_values("guide") for (c, d, t), g in tg.groupby(["cond", "D", "pert"])}
ALLT = sorted({t for (_, _, t) in IDX280})
ka1 = []
for cond in CONDS:
    tt = [t for t in ALLT if all((cond, d, t) in IDX280 for d in DORDER)]
    pooled = []
    for (a, b) in PAIRINGS:
        A = np.zeros((len(tt), len(CC)))
        B = np.zeros((len(tt), len(CC)))
        for side, pair, M in ((0, a, A), (1, b, B)):
            bb = bank.ctrl_profile(cond, pair)[CC]
            rw = rows_for(IDX280, cond, pair, tt)
            for k, t in enumerate(tt):
                M[k] = bank.profile(bank.counts(rw[t]))[CC] - bb
        for M in (A, B):
            M -= M.mean(0, keepdims=True)
            M[:, ~SC] = 0.0
        cs = cosines(A[:, SC], B[:, SC])
        pooled.append(cs)
        key = f"{''.join(a)}|{''.join(b)}"
        m, lo, hi, n = boot_median(cs)
        ref = R280[cond][key][0]
        d = m - ref
        ka1.append(dict(cond=cond, pairing=key, median=m, r280=ref, delta=d, n=n))
        P(f"  KA-1 {cond:<9} {key}  n={n:3d}  median {m:.10f}  (280 {ref:.10f})  Δ {d:+.2e}")
    pl = np.concatenate(pooled)
    m, lo, hi, n = boot_median(pl)
    rm, rlo, rhi, rn = R280[cond]["POOLED"]
    d = m - rm
    ka1.append(dict(cond=cond, pairing="POOLED", median=m, r280=rm, delta=d, n=n, lo=lo, hi=hi,
                    d_lo=lo - rlo, d_hi=hi - rhi, d_n=int(n - rn)))
    P(f"  KA-1 {cond:<9} POOLED     n={n:3d}  median {m:.10f} [{lo:.10f},{hi:.10f}]  "
      f"(280 {rm:.10f})  Δ {d:+.2e}  Δci {lo-rlo:+.1e}/{hi-rhi:+.1e}  Δn {int(n-rn)}")
mx = max(abs(r["delta"]) for r in ka1)
# preregistered bar: <= 5e-5 (R280 prints 4 d.p.); measured here against 280's full-precision CSV
OUT["KA1"] = dict(rows=ka1, max_abs_delta=float(mx), bar=5e-5, pass_=bool(mx <= 5e-5))
P(f"  KA-1 max|Δ| = {mx:.2e}  (bar 5e-5; this is floating-point identity, ~1 ulp)  -> {'PASS' if mx <= 5e-5 else 'FAIL'}")
if mx > 5e-5:
    P("ABORT: the ceiling does not reproduce 280.")
    json.dump(OUT, open("ka282.json", "w"), indent=1)
    sys.exit(1)

# ===================================================================== KA-2 extraction integrity
z = np.load(f"{W280}/obs280.npz", allow_pickle=True)
sel = np.load(f"{W280}/bank280_rows.npy")
tot = z["total"][sel]
rng = np.random.default_rng(20280906)
sample = np.sort(rng.choice(len(sel), 4000, replace=False))
mx2 = 0.0
for lo_ in range(0, len(sample), 500):
    idx = sample[lo_:lo_ + 500]
    got = np.asarray(bank.X[idx], dtype=np.float64).sum(1)
    mx2 = max(mx2, float(np.abs(got - tot[idx]).max()))
OUT["KA2"] = dict(n=4000, max_abs_delta=mx2, pass_=bool(mx2 == 0.0))
P(f"  KA-2 extraction integrity on 4,000 rows: max|rowsum - obs.total| = {mx2:g}  "
  f"-> {'PASS' if mx2 == 0 else 'FAIL'}")

# ===================================================================== KA-3 conservation join + C_noCD4
c276 = pd.read_csv(f"{ROOT}/experiments/276-conservation-map/work/conservation276.csv").set_index("target")
n276 = np.load(f"{ROOT}/experiments/276-conservation-map/work/cons276.npz", allow_pickle=True)
keys = [str(x) for x in n276["keys"]]
tgs = [str(x) for x in n276["targets"]]
Cc, pc = n276["Ccross"], n276["pairs_cross"]
c_raw_re = np.nanmean(np.where(np.isfinite(Cc), Cc, np.nan), axis=0)
d3 = np.nanmax(np.abs(c_raw_re - c276.reindex(tgs)["C_raw"].values))
cd4i = [i for i, k in enumerate(keys) if k.startswith("cd4_")]
mask_cd4 = np.isin(pc[:, 0], cd4i) | np.isin(pc[:, 1], cd4i)
with np.errstate(invalid="ignore"):
    C_noCD4 = np.nanmean(np.where(np.isfinite(Cc[~mask_cd4]), Cc[~mask_cd4], np.nan), axis=0)
n_nocd4 = np.isfinite(Cc[~mask_cd4]).sum(0)
targets, IDXc = build_index(bank, "common")
have = np.isin(targets, tgs).all()
OUT["KA3"] = dict(max_abs_delta=float(d3), pass_=bool(d3 <= 1e-6 and have),
                  n_analysis_targets=len(targets), all_have_row=bool(have),
                  n_pairs_dropped=int(mask_cd4.sum()), n_pairs_total=int(len(pc)),
                  cd4_sources=[keys[i] for i in cd4i])
P(f"  KA-3 re-derived C_raw vs published: max|Δ| = {d3:.2e}; every analysis target has a row: {have}; "
  f"dropped {mask_cd4.sum()}/{len(pc)} CD4 pairs -> {'PASS' if d3 <= 1e-6 and have else 'FAIL'}")
np.savez("cnocd4_282.npz", targets=np.asarray(tgs), C_noCD4=C_noCD4, n_pairs_noCD4=n_nocd4,
         C_raw=c276.reindex(tgs)["C_raw"].values)

# ===================================================================== KA-5 design reconciliation
rows5 = []
for filt, lab in (("common", "282 primary (state-common guides)"), ("keepDE", "280 filter"), ("all", "all guides")):
    tt, IDX = build_index(bank, filt)
    nrow = sum(len(IDX[(c, d, t)]) for c in CONDS for d in DORDER for t in tt)
    P(f"       {lab:<36} rows over the analysis set: {nrow}")
    for cond in CONDS:
        cells = np.array([sum(float(IDX[(cond, d, t)].n_cells.sum()) for d in DORDER) for t in tt])
        rows5.append(dict(filt=filt, cond=cond, n_targets=len(tt), med_cells_4donor=float(np.median(cells))))
    P(f"  KA-5 {lab:<36} n={len(tt):3d}  median 4-donor cells/(target,state): "
      + " / ".join(f"{r['med_cells_4donor']:.0f}" for r in rows5[-3:]))
OUT["KA5"] = dict(rows=rows5, r280_all_guides=[500, 476, 481])
P("       R280 KA-5 (all guides, its own per-condition sets): 500 / 476 / 481")
OUT["elapsed_s"] = time.time() - t0
json.dump(OUT, open("ka282.json", "w"), indent=1, default=float)
P(f"[{time.time()-t0:.0f}s]")
