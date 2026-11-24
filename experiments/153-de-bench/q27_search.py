"""Q27 — search (keep_top, keep_spread, tail_spread, phi) for a law landing n_pred ~1,000.

a20 measured -14.3% (4420 -> 3789), which leaves coverage saturated under BOTH H_A and H_B, so
it cannot separate them.  What separates them is an arm near n_pred ~1,000:

    H_B (n_real 240):   coverage 1        -> fid = precision       (~flat)
    H_A (n_real 2,500): coverage 0.40     -> fid ~0.50 -> ~0.20    (factor-2 collapse)

Design: a policy only assigns each gene to one (spread, phi) LEVEL, so P(significant) is
computed per gene per level ONCE and any policy composes by selection.  5 draws per gene per
level to damp the single-draw noise that made q26's delta unreliable.

CARRIES KNOWN-ANSWER ROWS: a15 and a20 are in the table with their measured artifact values,
so both calibrations are visible rather than assumed.
"""
import anndata as ad, numpy as np, scipy.sparse as sp, csv

TS, N, NDRAW = 5e4, 400, 5
rng = np.random.default_rng(0)

def lr(w, T):
    if T <= 0: return np.zeros(w.size)
    q = w / w.sum() * T; fl = np.floor(q); r = int(round(T - fl.sum()))
    if r > 0: fl[np.argsort(-(q - fl))[:r]] += 1
    return fl

A = ad.read_h5ad("../../data/challenge/context_A.h5ad")
X = A.X.tocsr().astype(np.float64)
lib = np.asarray(X.sum(1)).ravel(); lib[lib == 0] = 1.0
L = int(np.median(lib)); D = (sp.diags(1e6 / lib) @ X).tocsc(); m = D.shape[0]
gated = np.flatnonzero(np.asarray(D.mean(0)).ravel() > 5.0)
colsum = np.asarray(X.sum(0)).ravel(); p0 = colsum / colsum.sum()
fac = np.sqrt(12.0 * N * m / (m + N + 1.0)); Ls = np.full(N, float(L))
canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
z104 = np.load("../104-context-heldout-credit-adapter/work/hct116/frozen.npz", allow_pickle=True)
g104 = {str(g): i for i, g in enumerate(z104["genes"])}
eff = z104["baseline"][0][[g104[g] for g in canon]].astype(np.float64)
f = np.exp(eff); p_eff = p0 * f; p_eff /= p_eff.sum()

LEVELS = [(s, ph) for s in (0.25, 0.40, 0.55, 0.70, 1.00) for ph in (4.0, 8.0, 16.0)]
psig = np.zeros((gated.size, len(LEVELS)))
for k, g in enumerate(gated):
    T = int(np.floor(N * L * p_eff[g]))
    nz = np.sort(D.data[D.indptr[g]:D.indptr[g+1]]); n_zero = m - nz.size
    for j, (s, ph) in enumerate(LEVELS):
        if T <= 0: continue
        t_r = int(np.floor(T * s)); t_e = T - t_r
        base = lr(Ls, t_e); hits = 0
        for _ in range(NDRAW):
            c = base.copy()
            if t_r > 0:
                c = c + rng.poisson(rng.gamma(1.0/ph, ph * np.full(N, t_r / N)))
            vals, cts = np.unique(c * 1e6 / L, return_counts=True); tot = 0.0
            for v, ct in zip(vals, cts):
                if v == 0.0: below, eq = 0, n_zero
                else:
                    lo = np.searchsorted(nz, v, "left"); hi = np.searchsorted(nz, v, "right")
                    below, eq = n_zero + lo, hi - lo
                tot += ct * (below + 0.5 * eq) / m
            hits += abs((tot / N - 0.5) * fac) > 2.8
        psig[k, j] = hits / NDRAW

rank = np.argsort(-np.abs(eff[gated]))
def npred(keep_top, keep_lv, tail_lv):
    lvl = np.full(gated.size, tail_lv); lvl[rank[:keep_top]] = keep_lv
    return float(psig[np.arange(gated.size), lvl].sum())

idx = {lv: j for j, lv in enumerate(LEVELS)}
a15 = npred(8000, idx[(0.25, 4.0)], idx[(1.00, 4.0)])
a20 = npred(14000, idx[(0.25, 4.0)], idx[(0.25, 4.0)])
print(f"ANCHORS   model     artifact    ratio")
print(f"a15     {a15:7.0f}      4420.0     {4420/max(a15,1):.3f}")
print(f"a20     {a20:7.0f}      3789.0     {3789/max(a20,1):.3f}")
lvl_cal = (4420/max(a15,1) + 3789/max(a20,1)) / 2
print(f"level calibration (mean of the two anchors): x{lvl_cal:.3f}\n")
cands = []
for kt in (1000, 2000, 4000, 8000, 14000):
    for kl in LEVELS:
        for tl in LEVELS:
            v = npred(kt, idx[kl], idx[tl])
            cands.append((abs(v * lvl_cal - 1000), v, v * lvl_cal, kt, kl, tl))
cands.sort()
print(f"{'keep_top':>9} {'keep(spread,phi)':>18} {'tail(spread,phi)':>18} {'model':>7} {'x cal ->':>9}")
seen = set()
for _, v, vc, kt, kl, tl in cands:
    key = (kt, kl, tl)
    if kl == tl and kt != 14000: continue
    if len(seen) >= 8: break
    seen.add(key)
    print(f"{kt:>9} {str(kl):>18} {str(tl):>18} {v:7.0f} {vc:9.0f}")

# ---- global reachable range, and the mixed-law recipe -----------------------------------
allv = [npred(kt, idx[kl], idx[tl]) for kt in (500,1000,2000,4000,8000,14000)
        for kl in LEVELS for tl in LEVELS]
lo, hi = min(allv), max(allv)
print(f"\nPARAMETRIC REACHABLE RANGE over the whole grid: model {lo:.0f}-{hi:.0f}"
      f"  -> calibrated {lo*lvl_cal:.0f}-{hi*lvl_cal:.0f}")
print(f"target 1000 actual  = {1000/lvl_cal:.0f} model  -> "
      f"{'REACHABLE' if lo <= 1000/lvl_cal <= hi else 'NOT REACHABLE by any parametric policy'}")

# resampling sig-rate measured in q24: 1.2% of gated genes
r_rate = 0.012
best = psig[:, idx[(0.25, 4.0)]]                 # a15's keep level, per-gene P(sig)
p_rate = float(best.mean())
print(f"\nper-gene P(significant): parametric a15-keep level {p_rate:.4f}   resampling {r_rate:.4f}")
tgt_model = 1000 / lvl_cal
x = (tgt_model / gated.size - r_rate) / (p_rate - r_rate)
print(f"MIXED LAW: put x of genes on the parametric emitter, 1-x on resampling.")
print(f"   x = {x:.3f}  ->  {x*gated.size:.0f} genes parametric, {(1-x)*gated.size:.0f} resampled")
print(f"   predicted model n_pred {tgt_model:.0f}  -> calibrated ~{tgt_model*lvl_cal:.0f} actual")
for tgt in (700, 1000, 1500, 2000):
    tm = tgt / lvl_cal
    xx = (tm / gated.size - r_rate) / (p_rate - r_rate)
    print(f"   target {tgt:5d} actual -> x = {xx:.3f} ({xx*gated.size:5.0f} genes parametric)")
