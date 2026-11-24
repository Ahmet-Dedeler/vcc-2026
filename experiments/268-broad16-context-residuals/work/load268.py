"""Shared loader for experiment 268: half-pseudobulk responses on a common gene axis, working lines only."""
import json, os, numpy as np, pandas as pd
from pathlib import Path
W = Path("~/Code/virtual-cell-challenge-2026/experiments/268-broad16-context-residuals/work")
RES = json.load(open(W / "reserve.json")); WORKING = RES["working"]; RESERVED = RES["reserved"]
if os.environ.get("LINES268"): WORKING = [l for l in os.environ["LINES268"].split(",")]; assert not set(WORKING) & set(RESERVED); print("SMOKE TEST on", WORKING)
INV = json.load(open(W / "inventory.json")); CONTROLS = INV["controls"]; ALL_TARGETS = INV["targets"]; COMPLETE = INV["complete"]
META = pd.read_csv("~/Code/virtual-cell-challenge-2026/data/external/depmap26q1/extracted/metadata/cell_line_metadata.csv").set_index("cell_line")
TS = 5e4  # the scorer's pseudobulk target sum (log1p(TS * p / sum p))

def load_halves(lines=None):
    lines = list(lines or WORKING)
    H = {}
    for l in lines:
        assert l not in RESERVED, l
        z = np.load(W / "halves" / f"{l}.npz", allow_pickle=True)
        H[l] = {k: z[k] for k in ("full", "A", "B", "v1", "v2", "ncell", "libmean")}; H[l]["targets"] = list(z["targets"]); genes = z["genes"]
    return H, genes

def logprof(counts, mask):
    """log1p(TS * p / sum p) over the gene mask; counts: (..., G) float."""
    c = counts[..., mask].astype(np.float64); tot = c.sum(-1, keepdims=True)
    return np.log1p(TS * c / np.maximum(tot, 1.0))

def build_responses(H, genes, lines, targets, keys=("full", "A", "B", "v1", "v2"), min_cpm=None, filter_lines=None):
    """Returns dict key -> array (n_targets, n_lines, G_kept) of responses (perturbed log-profile minus the same-half
    control log-profile, controls = the 10 OR knockouts summed), the kept gene mask, and per-(t,l) availability.
    Gene axis: genes measured (non-NaN) in every line of `lines`; the 100 Broad target genes removed; optional mean-control
    CPM filter fitted on `filter_lines` (default: `lines`)."""
    G = len(genes); ok = np.ones(G, bool)
    for l in lines:
        ok &= ~np.isnan(H[l]["full"]).any(0)
    ok &= ~np.isin(genes, ALL_TARGETS)
    if min_cpm is not None:
        fl = filter_lines or lines; cpm = []
        for l in fl:
            ci = [H[l]["targets"].index(c) for c in CONTROLS]; c = np.nansum(H[l]["full"][ci], 0)
            cpm.append(1e6 * c / c.sum())
        ok &= np.nanmean(np.array(cpm), 0) >= min_cpm
    out = {}; avail = np.zeros((len(targets), len(lines)), bool); avail_key = {}
    for k in keys:
        R = np.zeros((len(targets), len(lines), int(ok.sum())), np.float32); ak = np.zeros((len(targets), len(lines)), bool)
        for j, l in enumerate(lines):
            tl = H[l]["targets"]; ci = [tl.index(c) for c in CONTROLS]
            ctrl = logprof(np.nansum(H[l][k][ci], 0), ok)
            for i, t in enumerate(targets):
                if t not in tl: continue
                row = H[l][k][tl.index(t)]
                if np.isnan(row[ok]).any() or row[ok].sum() <= 0: continue
                R[i, j] = logprof(row, ok) - ctrl; ak[i, j] = True
        out[k] = R; avail_key[k] = ak
        if k == "full": avail = ak
    out["avail"] = avail_key
    return out, ok, avail

def control_profiles(H, genes, lines, ok, key="full"):
    """Basal (control) log-profiles per line on the kept gene axis."""
    C = np.zeros((len(lines), int(ok.sum())), np.float64)
    for j, l in enumerate(lines):
        ci = [H[l]["targets"].index(c) for c in CONTROLS]; C[j] = logprof(np.nansum(H[l][key][ci], 0), ok)
    return C

def usable_targets(H, lines, min_cells=6, targets=None):
    """Complete targets with >= min_cells assigned cells in EVERY line of `lines` (so both random halves have >= 3)."""
    targets = list(targets or COMPLETE); keep = []
    for t in targets:
        if all(H[l]["ncell"][H[l]["targets"].index(t), 0] >= min_cells for l in lines): keep.append(t)
    return keep
