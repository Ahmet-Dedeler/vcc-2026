"""266 — one streaming pass over K562-GWPS X: cache ALL control cells (75,328) plus the 50 candidate
targets' cached seed-0 rows from 260 (experiments/260-*/work/extract_rows.npz) as float32 raw counts in
work/counts266.npy (memmap), with per-row gemgroup / UMI / total / zero-fraction in work/rows266.npz.
Nothing written outside experiments/266-*/work. No h5ad, no scorer."""
import json, time, resource
import h5py, numpy as np, pandas as pd
from pathlib import Path

P = "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
W = Path("experiments/266-k562-de-matched-controls/work")
R260 = Path("experiments/260-raw-de-diagnostic-feasibility/work")
BLOCK = 3886

def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 3

t_all = time.time()
Z = np.load(R260 / "extract_rows.npz")
cand = pd.read_csv(R260 / "candidates_de_power.csv")
targets = cand.target.tolist(); assert len(targets) == 50
with h5py.File(P, "r") as f:
    cats = f["obs/perturbation/categories"][:].astype(str); codes = f["obs/perturbation/codes"][:]
    batch = f["obs/batch"][:]; ncounts = f["obs/ncounts"][:]
    genes = f["var/gene_name"][:].astype(str)
    idx = {c: i for i, c in enumerate(cats)}
    ctrl_all = np.flatnonzero(codes == idx["control"])
    trows = {t: np.asarray(Z[f"t_{t}"]) for t in targets}
    for t in targets: assert (codes[trows[t]] == idx[t]).all(), t
    assert np.isin(Z["ctrl_rows"], ctrl_all).all()
    sel = np.concatenate([ctrl_all] + [trows[t] for t in targets])
    assert len(np.unique(sel)) == len(sel)
    label = np.concatenate([np.zeros(len(ctrl_all), int)] + [np.full(len(trows[t]), i + 1) for i, t in enumerate(targets)])
    order = np.argsort(sel); sel_sorted = sel[order]
    print(f"selected {len(sel):,} rows: {len(ctrl_all):,} controls + {len(sel)-len(ctrl_all):,} cells of 50 targets", flush=True)
    X = f["X"]; n, G = X.shape
    M = np.lib.format.open_memmap(W / "counts266.npy", mode="w+", dtype=np.float32, shape=(len(sel), G))
    tot = np.zeros(len(sel)); zf = np.zeros(len(sel))
    t0 = time.time()
    for start in range(0, n, BLOCK):
        stop = min(start + BLOCK, n)
        lo = np.searchsorted(sel_sorted, start); hi = np.searchsorted(sel_sorted, stop)
        if hi > lo:
            blk = X[start:stop, :]
            sub = blk[sel_sorted[lo:hi] - start]
            M[order[lo:hi]] = sub
            tot[order[lo:hi]] = sub.sum(1, dtype=np.float64)
            zf[order[lo:hi]] = (sub == 0).mean(1)
        if (start // BLOCK) % 64 == 0:
            print(f"  block {start // BLOCK}/{n // BLOCK}  {time.time() - t0:.0f}s  rss {rss():.2f} GiB", flush=True)
    M.flush(); del M
    t_read = time.time() - t0
assert np.allclose(tot, ncounts[sel]), "row sums must equal obs.ncounts"
np.savez(W / "rows266.npz", sel=sel, label=label, gemgroup=batch[sel], umi=ncounts[sel], tot=tot, zero_frac=zf,
         targets=np.array(targets), genes=genes, ctrl_rows_260=np.asarray(Z["ctrl_rows"]))
json.dump(dict(n_rows=int(len(sel)), n_ctrl=int(len(ctrl_all)), pass_min=t_read / 60, peak_rss_gib=rss(),
               total_min=(time.time() - t_all) / 60), open(W / "extract266.json", "w"), indent=1)
print(f"pass {t_read/60:.1f} min; peak rss {rss():.2f} GiB; rowsum==ncounts OK", flush=True)
