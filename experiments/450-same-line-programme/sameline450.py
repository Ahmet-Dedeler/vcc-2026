"""450: does a SAME-LINE bank from a different study carry the recipient's shared response programme?

441 showed the DE sign gap to the leaders is the recipient context's own mean response (oracle 0.645 on gate-block genes
vs CLEAN4's 0.580), and that no CROSS-line bank recovers it. It never tested a same-line, cross-study bank. On H1 we have
one: Pan 2025's hESC CRISPRi screen (GSE295214 pseudobulk database, 2,354 perturbations; a different lab, library and
chemistry from the 2025 VCC H1 data that h1-test is cut from). If Pan's mean response predicts h1-test's gate signs near
the oracle, same-line public banks carry the line's programme across studies. That decides whether context A's
Nadig-Jurkat tilt (builder `clean4q80ja`, never on the board) is worth a slot, and whether October needs a same-line mode.

Truth: h1-test's own Wilcoxon table (`323/work/reference.npz`), gate = q < 0.05, finite non-zero lfc, target gene out.
Block = the 2,200 top control-CPM genes (`ref_lin`). Every predictor is one sign per gene, the same for every target.
Pan rows for h1-test's 85 targets are excluded from Pan's programme (no per-target information enters).
usage: uv run python experiments/450-same-line-programme/sameline450.py
"""
import json
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)
genes_all = ref["genes"].astype(str)
keep = ref["keep"]
G = genes_all[keep]                      # the 10,779 tested genes, in the table's column order
T = ref["targets"].astype(str)
lfc, q = ref["lfc"], ref["q"]
ref_lin = ref["ref_lin"]
gidx = {g: i for i, g in enumerate(G)}
block = np.zeros(len(G), bool)
block[np.argsort(-ref_lin)[:2200]] = True

gate = (q < 0.05) & np.isfinite(lfc) & (lfc != 0)
for i, t in enumerate(T):
    if t in gidx:
        gate[i, gidx[t]] = False
adj = np.isfinite(lfc) & (lfc != 0)       # adjudicable (fid-like reads every called gene against this sign)
for i, t in enumerate(T):
    if t in gidx:
        adj[i, gidx[t]] = False
truth = np.sign(lfc)


def on_axis(names, values):
    """Map a per-gene vector on another gene axis onto G (NaN where absent)."""
    out = np.full(len(G), np.nan)
    for n, v in zip(names, values):
        j = gidx.get(n)
        if j is not None:
            out[j] = v
    return out


# --- predictors (each: predicted sign per gene, +1 / -1, NaN = no call) ---------------------------------------------
P = {}
# CLEAN4 (the shipped cross-line prior): P_down over K562-GWPS, K562-ess, Nadig Jurkat, Nadig HepG2
c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
c4p = on_axis(c4["genes"].astype(str), c4["P_down"])
P["clean4"] = np.where(np.isnan(c4p), np.nan, np.where(c4p > 0.5, -1.0, 1.0))

# Pan 2025 hESC (same line, different study): P_down and mean delta over its rows, h1-test targets excluded
db = ROOT / "data/external/gse295214/pseudo_bulk_perturbation_database.h5ad"
with h5py.File(db, "r") as f:
    ds = f["obs/dataset"]
    cats = ds["categories"][:].astype(str) if isinstance(ds, h5py.Group) else None
    codes = ds["codes"][:] if isinstance(ds, h5py.Group) else None
    dataset = cats[codes] if cats is not None else ds[:].astype(str)
    pg = f["obs/perturbed_gene_name"]
    pgn = (pg["categories"][:].astype(str)[pg["codes"][:]] if isinstance(pg, h5py.Group) else pg[:].astype(str))
    vg = f["var/gene_name"]
    vnames = (vg["categories"][:].astype(str)[vg["codes"][:]] if isinstance(vg, h5py.Group) else vg[:].astype(str))
    rows = np.where((dataset == "Pan2025_hESC") & ~np.isin(pgn, T))[0]
    X = f["X"][np.sort(rows), :]
    M = f["layers/measure_mask"][np.sort(rows), :]
X = np.where(M > 0, X, np.nan)
n_rows = len(rows)
pan_pdown = np.nanmean(np.where(np.isnan(X), np.nan, (X < 0).astype(float)), axis=0)
pan_mean = np.nanmean(X, axis=0)
pp = on_axis(vnames, pan_pdown)
pm = on_axis(vnames, pan_mean)
P["pan_pdown"] = np.where(np.isnan(pp), np.nan, np.where(pp > 0.5, -1.0, 1.0))
P["pan_mean"] = np.where(np.isnan(pm), np.nan, np.sign(pm))


def z(v):
    return (v - np.nanmean(v)) / np.nanstd(v)


# combined: z(P_down) average of CLEAN4 and Pan (down-probability space), both available
comb = np.nanmean(np.vstack([z(c4p), z(pp)]), axis=0)
P["clean4+pan"] = np.where(np.isnan(comb), np.nan, np.where(comb > 0, -1.0, 1.0))


def score(pred_gene, loo=False):
    """Pooled and per-target sign accuracy on subsets; pred is per gene unless loo (per target, per gene)."""
    res = {}
    for name, sub in (("gate_block", gate & block), ("gate_nonblock", gate & ~block), ("gate_all", gate),
                      ("block_all_adj", adj & block)):
        if loo:
            pr = pred_gene
        else:
            pr = np.broadcast_to(pred_gene, lfc.shape)
        ok = sub & np.isfinite(pr)
        hit = (pr == truth) & ok
        w = np.abs(np.where(ok, lfc, 0))
        res[name] = {"acc": float(hit.sum() / max(ok.sum(), 1)),
                     "wacc": float((w * hit).sum() / max(w.sum(), 1e-12)),
                     "coverage": float(ok.sum() / max(sub.sum(), 1))}
    return res


out = {"n_pan_rows": int(n_rows), "n_targets": int(len(T)), "block": int(block.sum()), "gate_pairs": int(gate.sum())}
for k, v in P.items():
    out[k] = score(v)
# the LOO oracle: sign of the mean reference lfc over the OTHER targets (NaN-safe), per target
L = np.where(np.isfinite(lfc), lfc, 0.0)
S = L.sum(0)
loo = np.sign(-(L - S[None, :]))  # sign of sum over others = sign(S - L_i); note the minus of (L - S)
out["h1_loo_oracle"] = score(np.where(loo == 0, np.nan, loo), loo=True)
# agreement between predictors on the block
bl = block
out["agree_block"] = {
    "clean4_vs_pan": float(np.nanmean((P["clean4"] == P["pan_pdown"])[bl & np.isfinite(P["clean4"]) & np.isfinite(P["pan_pdown"])])),
    "pan_vs_oracle_sign_of_total": float(np.nanmean((P["pan_pdown"] == np.sign(S))[bl & np.isfinite(P["pan_pdown"]) & (S != 0)])),
    "clean4_vs_oracle_sign_of_total": float(np.nanmean((P["clean4"] == np.sign(S))[bl & np.isfinite(P["clean4"]) & (S != 0)])),
}
cor = lambda a, b: float(np.corrcoef(a[np.isfinite(a) & np.isfinite(b)], b[np.isfinite(a) & np.isfinite(b)])[0, 1])
out["corr_with_h1_mean_response"] = {"pan_mean": cor(pm, S / len(T)), "clean4_pdown(neg)": cor(-c4p, S / len(T)),
                                     "pan_pdown(neg)": cor(-pp, S / len(T))}
(OUT / "sameline450.json").write_text(json.dumps(out, indent=1))
print(f"Pan rows {n_rows} (h1-test targets excluded); gate pairs {gate.sum()}; block {block.sum()}")
print(f"{'predictor':16s} " + " ".join(f"{s:>22s}" for s in ("gate_block acc/w/cov", "gate_nonblock", "gate_all", "block_all_adj")))
for k in list(P) + ["h1_loo_oracle"]:
    r = out[k]
    print(f"{k:16s} " + " ".join(f"{r[s]['acc']:.3f}/{r[s]['wacc']:.3f}/{r[s]['coverage']:.2f}".rjust(22) for s in
                                 ("gate_block", "gate_nonblock", "gate_all", "block_all_adj")))
print("agreement on block:", out["agree_block"])
print("corr with H1 mean response:", out["corr_with_h1_mean_response"])
