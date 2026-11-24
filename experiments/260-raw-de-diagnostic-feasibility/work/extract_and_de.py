"""260 — DE-power measurement on K562-GWPS raw counts for the challenge panel targets.

One streaming pass over the dense gzip X (row blocks of 3,886, the file's chunk height),
keeping up to CAP cells per panel target present in GWPS and N_CTRL control cells drawn
uniformly (seed 0, so the draw is reproducible). Then, per target, the scorer's DE recipe
(docs/reference/vcc2026-metrics-BRIEF.md §3): per-cell CPM to 1e6, genes with control mean
CPM > 5, two-sided Mann-Whitney U (normal approximation, tie-corrected AND uncorrected),
BH per target, alpha 0.05; n_real excludes the target gene. NOT a harness: nothing is
written outside experiments/260-*/work, no h5ad, no scorer.
"""
import json, sys, time, resource
import h5py, numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import rankdata, norm

P = "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
OUT = Path("experiments/260-raw-de-diagnostic-feasibility/work")
CAP, N_CTRL, SEED, BLOCK = 400, 18_400, 0, 3886
MIN_CPM, ALPHA, GATE = 5.0, 0.05, 10

def rss_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 3

t_all = time.time()
panel = pd.read_csv("data/challenge/pert_counts.csv")["target_gene"].astype(str).tolist()
rng = np.random.default_rng(SEED)
with h5py.File(P, "r") as f:
    cats = f["obs/perturbation/categories"][:].astype(str); codes = f["obs/perturbation/codes"][:]
    genes = f["var/gene_name"][:].astype(str)
    idx = {c: i for i, c in enumerate(cats)}
    ctrl_code = idx["control"]
    keep = {}
    for t in panel:
        if t not in idx: continue
        rows = np.flatnonzero(codes == idx[t])
        if len(rows) > CAP: rows = rng.choice(rows, CAP, replace=False)
        keep[t] = np.sort(rows)
    ctrl_rows = np.sort(rng.choice(np.flatnonzero(codes == ctrl_code), N_CTRL, replace=False))
    sel = np.concatenate([ctrl_rows] + [keep[t] for t in keep])
    order = np.argsort(sel); sel_sorted = sel[order]
    print(f"selected {len(sel):,} rows: {N_CTRL:,} control + {len(sel) - N_CTRL:,} from {len(keep)} targets", flush=True)
    X = f["X"]; n = X.shape[0]
    M = np.zeros((len(sel), X.shape[1]), dtype=np.float32)
    t0 = time.time(); pos = 0
    for start in range(0, n, BLOCK):
        stop = min(start + BLOCK, n)
        lo = np.searchsorted(sel_sorted, start); hi = np.searchsorted(sel_sorted, stop)
        if hi > lo:
            blk = X[start:stop, :]
            M[order[lo:hi]] = blk[sel_sorted[lo:hi] - start]
        if (start // BLOCK) % 64 == 0:
            print(f"  block {start // BLOCK}/{n // BLOCK}  {time.time() - t0:.0f}s  rss {rss_gib():.2f} GiB", flush=True)
    t_read = time.time() - t0
print(f"pass done in {t_read / 60:.1f} min; rss {rss_gib():.2f} GiB", flush=True)
np.savez(OUT / "extract_rows.npz", ctrl_rows=ctrl_rows, **{f"t_{t}": r for t, r in keep.items()})

# ---- DE
tot = M.sum(1, keepdims=True)
ctrl = M[:N_CTRL]; ctrl_cpm = ctrl / tot[:N_CTRL] * 1e6
ctrl_mean = ctrl_cpm.mean(0)
gate = ctrl_mean > MIN_CPM
print(f"genes passing control mean CPM > {MIN_CPM}: {gate.sum()} / {len(genes)}", flush=True)
G = np.flatnonzero(gate); ctrl_cpm_g = ctrl_cpm[:, G]
ctrl_med_umi = float(np.median(tot[:N_CTRL]))

def bh(p):
    m = len(p); o = np.argsort(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)

def mwu_both(A, B, chunk=800):
    """A: pert (n1,G) B: ctrl (n2,G). Returns two-sided p (tie-corrected), p (uncorrected)."""
    n1, n2 = A.shape[0], B.shape[0]; N = n1 + n2
    p_tc = np.empty(A.shape[1]); p_uc = np.empty(A.shape[1])
    mu = n1 * (N + 1) / 2.0
    for s in range(0, A.shape[1], chunk):
        e = min(s + chunk, A.shape[1])
        Z = np.concatenate([A[:, s:e], B[:, s:e]], axis=0).astype(np.float64)
        R = rankdata(Z, axis=0)
        R1 = R[:n1].sum(0)
        # tie term per gene
        Zs = np.sort(Z, axis=0)
        T = np.zeros(e - s)
        for j in range(e - s):
            _, c = np.unique(Zs[:, j], return_counts=True)
            T[j] = (c ** 3 - c).sum()
        var_uc = n1 * n2 * (N + 1) / 12.0
        var_tc = n1 * n2 / 12.0 * ((N + 1) - T / (N * (N - 1)))
        z_uc = (R1 - mu) / np.sqrt(var_uc)
        z_tc = np.where(var_tc > 0, (R1 - mu) / np.sqrt(np.maximum(var_tc, 1e-300)), 0.0)
        p_uc[s:e] = 2 * norm.sf(np.abs(z_uc)); p_tc[s:e] = 2 * norm.sf(np.abs(z_tc))
    return p_tc, p_uc

rows = []
pos = N_CTRL
gene_pos = {g: i for i, g in enumerate(genes)}
t0 = time.time()
targets = list(keep)
for k, t in enumerate(targets):
    n_p = len(keep[t]); A = M[pos:pos + n_p]; pos += n_p
    A_cpm = A / tot[pos - n_p:pos] * 1e6
    Ag = A_cpm[:, G]
    p_tc, p_uc = mwu_both(Ag, ctrl_cpm_g)
    lfc = np.log2((Ag.mean(0) + 1e-9) / (ctrl_mean[G] + 1e-9))
    self_mask = np.array([genes[g] != t for g in G])
    q_tc = bh(p_tc); q_uc = bh(p_uc)
    sig_tc = (q_tc < ALPHA) & self_mask; sig_uc = (q_uc < ALPHA) & self_mask
    fin = np.isfinite(lfc) & (lfc != 0)
    rows.append(dict(target=t, cells=n_p, median_umi=float(np.median(tot[pos - n_p:pos])),
                     n_real_tc=int(sig_tc.sum()), n_real_uc=int(sig_uc.sum()),
                     gate_tc=int((sig_tc & fin).sum()), gate_uc=int((sig_uc & fin).sum()),
                     frac_neg_tc=float((lfc[sig_tc] < 0).mean()) if sig_tc.any() else np.nan,
                     target_in_axis=t in gene_pos, target_knockdown_lfc=float(lfc[G == gene_pos[t]][0]) if (t in gene_pos and gene_pos[t] in set(G.tolist())) else np.nan))
    if k % 20 == 0:
        print(f"  {k}/{len(targets)} {t}: cells {n_p} n_real tc/uc {rows[-1]['n_real_tc']}/{rows[-1]['n_real_uc']}  {time.time() - t0:.0f}s rss {rss_gib():.2f}", flush=True)
df = pd.DataFrame(rows); df.to_csv(OUT / "de_power_panel_targets.csv", index=False)
print(f"\nDE done in {(time.time() - t0) / 60:.1f} min; peak rss {rss_gib():.2f} GiB")
print("control median UMI", ctrl_med_umi)
for lo, hi in ((0, 100), (100, 150), (150, 200), (200, 300), (300, 401)):
    d = df[(df.cells >= lo) & (df.cells < hi)]
    if len(d):
        print(f"cells [{lo},{hi}): n={len(d):3d}  n_real(tc) median {d.n_real_tc.median():6.1f} IQR [{d.n_real_tc.quantile(.25):.0f},{d.n_real_tc.quantile(.75):.0f}]  "
              f"gate>=10: {(d.gate_tc >= GATE).sum():3d}  | uc median {d.n_real_uc.median():6.1f} gate>=10 {(d.gate_uc >= GATE).sum():3d}")
print("all:", len(df), "gate>=10 tc", (df.gate_tc >= GATE).sum(), "uc", (df.gate_uc >= GATE).sum(),
      "frac_neg median", df.frac_neg_tc.median())
# subsampling curve on the deepest targets: n_real vs cells
big = df.sort_values("cells", ascending=False).head(3).target.tolist()
curve = []
pos_map = {}; p = N_CTRL
for t in targets: pos_map[t] = (p, p + len(keep[t])); p += len(keep[t])
for t in big:
    a, b = pos_map[t]; A_all = M[a:b]; T_all = tot[a:b]
    for nc in (100, 150, 200, 300, 400):
        if nc > (b - a): continue
        sub = np.random.default_rng(1).choice(b - a, nc, replace=False)
        Ag = (A_all[sub] / T_all[sub] * 1e6)[:, G]
        p_tc, _ = mwu_both(Ag, ctrl_cpm_g); q = bh(p_tc)
        self_mask = np.array([genes[g] != t for g in G])
        curve.append(dict(target=t, cells=nc, n_real_tc=int(((q < ALPHA) & self_mask).sum())))
        print("  curve", t, nc, curve[-1]["n_real_tc"], flush=True)
pd.DataFrame(curve).to_csv(OUT / "de_power_curve.csv", index=False)
json.dump(dict(n_selected=int(len(sel)), n_targets=len(targets), pass_min=t_read / 60, genes_gated=int(gate.sum()),
               ctrl_median_umi=ctrl_med_umi, peak_rss_gib=rss_gib(), total_min=(time.time() - t_all) / 60),
          open(OUT / "de_power_summary.json", "w"), indent=1)
