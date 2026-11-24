"""266 — the scorer's Wilcoxon recipe (260's, verbatim arithmetic) on 260's 50 candidate targets under three
control-selection schemes, three seeds each:
  pooled     : 18,400 controls uniform over all 75,328 (seed 0 = 260's exact draw -> known-answer row; seeds 1,2 fresh)
  gem        : 18,400 controls drawn from the target's own gemgroups in proportion to its cells' gemgroup counts
  gemdepth   : as gem, and within each gemgroup the controls are the nearest-UMI matches to the target cells'
               within-gemgroup UMI distribution (jittered quantile grid, greedy nearest unused control)
CPM is computed exactly as 260 did (float32 counts / float32 row total * 1e6, float32), so the pooled seed-0 row
must reproduce 260's n_real_tc to the digit. Ranks are exact (midranks) and the tie term exact, so z/p/BH are
bitwise 260's formulas. Peak RSS: memmap'd counts (file-backed, shared) + ~0.4 GiB per worker.
Writes only into experiments/266-*/work."""
import os, sys, time, json, resource
import numpy as np, pandas as pd
from pathlib import Path
from scipy.stats import norm
import multiprocessing as mp

W = Path("experiments/266-k562-de-matched-controls/work")
N_CTRL_DRAW, ALPHA, MIN_CPM, GATE, CH = 18_400, 0.05, 5.0, 10, 512
SEEDS = (0, 1, 2)
N_WORKERS = int(os.environ.get("N_WORKERS", "6"))

rows = np.load(W / "rows266.npz", allow_pickle=False)
sel, label, gemgroup, umi, tot, zf = rows["sel"], rows["label"], rows["gemgroup"], rows["umi"], rows["tot"], rows["zero_frac"]
targets = rows["targets"].astype(str).tolist(); genes = rows["genes"].astype(str)
n_ctrl = int((label == 0).sum()); G = len(genes)
assert (label[:n_ctrl] == 0).all() and (label[n_ctrl:] > 0).all()
tot32 = tot.astype(np.float32)
gene_pos = {g: i for i, g in enumerate(genes)}

def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 3

# ---------------------------------------------------------------- control draws (obs-only)
ctrl_gem = gemgroup[:n_ctrl]; ctrl_umi = umi[:n_ctrl]
by_gem = {}
for g in np.unique(ctrl_gem):
    p = np.flatnonzero(ctrl_gem == g); p = p[np.argsort(ctrl_umi[p], kind="stable")]
    by_gem[int(g)] = p                      # control positions in gemgroup g, sorted by UMI

def apportion(weights, total):
    """Hamilton largest-remainder apportionment of `total` integer units proportional to weights."""
    w = np.asarray(weights, float); q = total * w / w.sum(); m = np.floor(q).astype(int)
    r = total - m.sum()
    if r > 0:
        m[np.argsort(-(q - m), kind="stable")[:r]] += 1
    return m

def allocate(gems, counts, total=N_CTRL_DRAW):
    """controls per gemgroup ∝ target cell counts, capped by availability, shortfall redistributed."""
    avail = np.array([len(by_gem.get(int(g), ())) for g in gems]); m = np.zeros(len(gems), int)
    active = avail > 0; remaining = total
    for _ in range(50):
        if remaining == 0 or not active.any(): break
        add = np.zeros(len(gems), int); add[active] = apportion(counts[active], remaining)
        cap = np.minimum(add, avail - m); m += cap
        active &= (m < avail)
        remaining = total - m.sum()
    return m

def draw_pooled(seed):
    if seed == 0:
        pos = np.searchsorted(sel[:n_ctrl], rows["ctrl_rows_260"]); assert (sel[pos] == rows["ctrl_rows_260"]).all()
        return np.sort(pos)
    return np.sort(np.random.default_rng([266, 0, seed]).choice(n_ctrl, N_CTRL_DRAW, replace=False))

def draw_gem(ti, seed, depth):
    rng = np.random.default_rng([266, 2 if depth else 1, seed, ti])
    t_pos = np.flatnonzero(label == ti + 1); tg = gemgroup[t_pos]; tu = umi[t_pos]
    gems, counts = np.unique(tg, return_counts=True)
    m = allocate(gems, counts)
    out = []
    for g, c, mg in zip(gems, counts, m):
        if mg == 0: continue
        cp = by_gem[int(g)]
        if not depth:
            out.append(rng.choice(cp, mg, replace=False)); continue
        cu = ctrl_umi[cp].astype(float); tug = tu[tg == g]
        q = (np.arange(mg) + rng.random(mg)) / mg
        u_star = np.quantile(tug, q)
        used = np.zeros(len(cp), bool)
        for j in rng.permutation(mg):
            d = np.abs(cu - u_star[j]); d[used] = np.inf; used[np.argmin(d)] = True
        out.append(cp[used])
    pos = np.sort(np.concatenate(out)); assert len(np.unique(pos)) == len(pos)
    return pos

# ---------------------------------------------------------------- the test (260's arithmetic, exact ranks)
def bh(p):
    m = len(p); o = np.argsort(p); q = np.empty(m); q[o] = np.minimum.accumulate((p[o] * m / np.arange(1, m + 1))[::-1])[::-1]
    return np.minimum(q, 1.0)

MM = None
def _init():
    global MM
    MM = np.load(W / "counts266.npy", mmap_mode="r")

def run_de(job):
    scheme, seed, ti, ctrl_pos = job
    t0 = time.time()
    t_pos = np.flatnonzero(label == ti + 1); t = targets[ti]
    n1, n2 = len(t_pos), len(ctrl_pos); N = n1 + n2
    A_tot = tot32[t_pos][:, None]; B_tot = tot32[ctrl_pos][:, None]
    R1 = np.empty(G); T = np.empty(G); mean_p = np.empty(G); mean_c = np.empty(G)
    for s in range(0, G, CH):
        e = min(s + CH, G); k = e - s
        Ac = np.array(MM[t_pos, s:e]); Ac /= A_tot; Ac *= 1e6        # == 260's `A / tot * 1e6` in float32
        Bc = np.array(MM[ctrl_pos, s:e]); Bc /= B_tot; Bc *= 1e6
        assert Ac.dtype == np.float32 and Bc.dtype == np.float32
        mean_p[s:e] = Ac.mean(0, dtype=np.float64); mean_c[s:e] = Bc.mean(0, dtype=np.float64)
        Zt = np.empty((k, N), np.float32); Zt[:, :n1] = Ac.T; Zt[:, n1:] = Bc.T
        idx = np.argsort(Zt, axis=1)
        Zs = np.take_along_axis(Zt, idx, axis=1)
        pert_sorted = idx < n1
        new = np.ones((k, N), bool); new[:, 1:] = Zs[:, 1:] != Zs[:, :-1]
        starts = np.flatnonzero(new.ravel()); lengths = np.diff(starts, append=k * N)
        midrank = (starts % N) + (lengths + 1) / 2.0                  # 1-based midrank of each tie run
        ranks = np.repeat(midrank, lengths).reshape(k, N)
        R1[s:e] = (ranks * pert_sorted).sum(1)
        T[s:e] = np.bincount(starts // N, weights=lengths.astype(np.float64) ** 3 - lengths, minlength=k)
    mu = n1 * (N + 1) / 2.0
    var_uc = n1 * n2 * (N + 1) / 12.0
    var_tc = n1 * n2 / 12.0 * ((N + 1) - T / (N * (N - 1)))
    z_uc = (R1 - mu) / np.sqrt(var_uc)
    z_tc = np.where(var_tc > 0, (R1 - mu) / np.sqrt(np.maximum(var_tc, 1e-300)), 0.0)
    p_uc = 2 * norm.sf(np.abs(z_uc)); p_tc = 2 * norm.sf(np.abs(z_tc))
    gate = mean_c > MIN_CPM                                           # the scorer's gate on THIS control pool
    Gi = np.flatnonzero(gate)
    lfc = np.log2((mean_p + 1e-9) / (mean_c + 1e-9))
    self_mask = genes != t
    q_tc = np.full(G, np.nan); q_uc = np.full(G, np.nan)
    q_tc[Gi] = bh(p_tc[Gi]); q_uc[Gi] = bh(p_uc[Gi])
    sig_tc = (q_tc < ALPHA) & self_mask; sig_uc = (q_uc < ALPHA) & self_mask
    fin = np.isfinite(lfc) & (lfc != 0)
    c_umi = umi[ctrl_pos]; t_umi = umi[t_pos]
    tg = gemgroup[t_pos]; cg = gemgroup[ctrl_pos]
    gu = np.union1d(tg, cg)
    tv = 0.5 * np.abs(np.bincount(np.searchsorted(gu, tg), minlength=len(gu)) / n1
                      - np.bincount(np.searchsorted(gu, cg), minlength=len(gu)) / n2).sum()
    rec = dict(scheme=scheme, seed=seed, target=t, cells=n1, n_ctrl=n2, genes_gated=int(gate.sum()),
               n_real_tc=int(sig_tc.sum()), n_real_uc=int(sig_uc.sum()),
               gate_tc=int((sig_tc & fin).sum()), gate_uc=int((sig_uc & fin).sum()),
               frac_neg_tc=float((lfc[sig_tc] < 0).mean()) if sig_tc.any() else np.nan,
               target_knockdown_lfc=float(lfc[gene_pos[t]]) if t in gene_pos else np.nan,
               median_umi_t=float(np.median(t_umi)), median_umi_c=float(np.median(c_umi)),
               depth_ratio=float(np.median(t_umi) / np.median(c_umi)),
               ctrl_umi_cv=float(c_umi.std() / c_umi.mean()), ctrl_umi_mean=float(c_umi.mean()),
               zero_gap=float(zf[t_pos].mean() - zf[ctrl_pos].mean()),
               n_gem_t=int(len(np.unique(tg))), n_gem_c=int(len(np.unique(cg))), gem_tv=float(tv),
               min_ctrl_mean_cpm=float(mean_c.min()), secs=time.time() - t0, rss_gib=rss())
    return rec, p_tc, lfc.astype(np.float32)

if __name__ == "__main__":
    t_all = time.time()
    jobs = []
    draws_meta = []
    for seed in SEEDS:
        pooled = draw_pooled(seed)
        for ti in range(len(targets)):
            jobs.append(("pooled", seed, ti, pooled))
            jobs.append(("gem", seed, ti, draw_gem(ti, seed, depth=False)))
            jobs.append(("gemdepth", seed, ti, draw_gem(ti, seed, depth=True)))
    print(f"{len(jobs)} jobs drawn in {time.time() - t_all:.0f}s; workers {N_WORKERS}", flush=True)
    recs = []; tables = {}
    with mp.get_context("fork").Pool(N_WORKERS, initializer=_init) as pool:
        for i, (rec, p_tc, lfc) in enumerate(pool.imap_unordered(run_de, jobs, chunksize=2)):
            recs.append(rec); tables[f"{rec['scheme']}_s{rec['seed']}_{rec['target']}_p"] = p_tc
            tables[f"{rec['scheme']}_s{rec['seed']}_{rec['target']}_lfc"] = lfc
            if i % 25 == 0:
                print(f"  {i}/{len(jobs)} {rec['scheme']} s{rec['seed']} {rec['target']}: n_real {rec['n_real_tc']} "
                      f"{rec['secs']:.1f}s worker-rss {rec['rss_gib']:.2f} GiB  {(time.time() - t_all)/60:.1f} min", flush=True)
    df = pd.DataFrame(recs).sort_values(["scheme", "seed", "target"])
    df.to_csv(W / "runs266.csv", index=False)
    np.savez_compressed(W / "de_tables266.npz", genes=genes, **tables)
    json.dump(dict(n_jobs=len(jobs), workers=N_WORKERS, total_min=(time.time() - t_all) / 60,
                   parent_rss_gib=rss(), worker_rss_gib_max=float(df.rss_gib.max())), open(W / "pilot266.json", "w"), indent=1)
    print(f"done {len(jobs)} runs in {(time.time() - t_all)/60:.1f} min; parent rss {rss():.2f} GiB; max worker rss {df.rss_gib.max():.2f} GiB", flush=True)
