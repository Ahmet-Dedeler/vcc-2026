"""343 -- DESCRIPTIVE decomposition of the H1 NMAE deficit for the cap15 (319-composition) arm.

Reads only 323's saved outputs (stats-cap15.npz, reference_de.parquet, reference.npz axis). No emission,
no scoring library, no gates. NMAE_p = sum|lhat-l| / sum|l| over the reference gate S_p (q<.05, finite
nonzero lfc, own target removed). With m = sum|lhat|/sum|l| and sign-weighted accuracy p_w, an all-correct-
sign prediction gives 1-m and an all-wrong one 1+m; the observed value sits between them. Counterfactuals
recomputed exactly per target, not from that identity.
"""
from __future__ import annotations
import json, sys
import numpy as np, polars as pl
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]
W = REPO / "experiments/323-current-family-validation/work"
OUT = Path(__file__).resolve().parent / "work"
ARM = sys.argv[1] if len(sys.argv) > 1 else "cap15"

s = np.load(W / f"stats-{ARM}.npz", allow_pickle=False)
r = np.load(W / "reference.npz", allow_pickle=False)
genes_all = r["genes"].astype(str); keep = r["keep"].astype(bool); genes = genes_all[keep]
targets = s["targets"].astype(str)
assert (targets == r["targets"].astype(str)).all()
ref = pl.read_parquet(W / "reference_de.parquet")
cols = ref.columns
tcol = next(c for c in cols if c in ("target", "perturbation", "target_gene"))
gcol = next(c for c in cols if c in ("feature", "gene", "gene_name"))
qcol = next(c for c in cols if "adj" in c or c == "q")
lcol = next(c for c in cols if "fold" in c or c == "lfc")
print("reference columns:", cols, "->", tcol, gcol, qcol, lcol)
gi = {g: i for i, g in enumerate(genes)}
lhat_all = np.asarray(s["lfc"], dtype=np.float64); qhat_all = np.asarray(s["q"], dtype=np.float64)
assert lhat_all.shape == (len(targets), len(genes)), lhat_all.shape

rows = []
for ti, t in enumerate(targets):
    sub = ref.filter(pl.col(tcol) == t)
    g = sub[gcol].to_numpy().astype(str); q = sub[qcol].to_numpy().astype(float); l = sub[lcol].to_numpy().astype(float)
    gate = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if gate.sum() < 10: continue
    idx = np.array([gi[x] for x in g[gate]]); l = l[gate]
    lh = lhat_all[ti, idx]; qh = qhat_all[ti, idx]
    lh = np.where(np.isfinite(lh), lh, 0.)
    D = np.abs(l).sum()
    nmae = np.abs(lh - l).sum() / D
    m = np.abs(lh).sum() / D
    pw = (np.abs(l) * (np.sign(lh) == np.sign(l))).sum() / D           # |l|-weighted sign accuracy (zero-pred counts wrong)
    pw_nz = (np.abs(l) * (np.sign(lh) == np.sign(l))).sum() / max(np.abs(l[lh != 0]).sum(), 1e-12)  # among genes we give any magnitude
    zero_share = np.abs(l[np.abs(lh) < 1e-6]).sum() / D                 # gate mass where we predict exactly nothing
    called = qh < .05
    called_share = np.abs(l[called]).sum() / D                           # gate mass inside our own significant set
    # counterfactuals, exact
    cf_sign = np.abs(np.abs(lh) - np.abs(l)).sum() / D                   # every sign right, our magnitudes
    a = (np.abs(lh) * np.abs(l)).sum() / max((lh**2).sum(), 1e-12)        # least-squares scalar on |.| ignoring sign
    cf_scale = np.abs(a * lh - l).sum() / D                              # one scalar gain per target, our signs
    cf_sign_scale = np.abs(a * np.abs(lh) - np.abs(l)).sum() / D
    # magnitude ratio inside our called set vs outside
    m_called = np.abs(lh[called]).sum() / max(np.abs(l[called]).sum(), 1e-12) if called.any() else np.nan
    m_uncalled = np.abs(lh[~called]).sum() / max(np.abs(l[~called]).sum(), 1e-12) if (~called).any() else np.nan
    rows.append(dict(target=t, n_gate=int(gate.sum()), nmae=nmae, m=m, pw=pw, pw_nonzero=pw_nz, zero_share=zero_share,
                     called_share=called_share, n_called_in_gate=int(called.sum()), m_called=m_called, m_uncalled=m_uncalled,
                     cf_all_signs_right=cf_sign, cf_best_scalar=cf_scale, cf_signs_and_scalar=cf_sign_scale, scalar=a))
T = pl.DataFrame(rows)
T.write_csv(OUT / f"per-target-{ARM}.csv")
def mean(c): return float(T[c].mean())
summary = {k: mean(k) for k in ("nmae","m","pw","pw_nonzero","zero_share","called_share","m_called","m_uncalled",
                                "cf_all_signs_right","cf_best_scalar","cf_signs_and_scalar","scalar")}
summary["n_targets"] = T.height; summary["arm"] = ARM
summary["nmae_323_recorded"] = None
try:
    mm = json.load(open(W / f"metrics-{ARM}.json")); summary["nmae_323_recorded"] = mm["raw"].get("de_wilcoxon_lfc_nmae") if isinstance(mm.get("raw"), dict) else None
except Exception as e: summary["nmae_323_recorded_error"] = str(e)
(OUT / f"summary-{ARM}.json").write_text(json.dumps(summary, indent=1))
print(json.dumps(summary, indent=1))
print(T.sort("nmae").head(8)); print(T.sort("nmae").tail(5))
