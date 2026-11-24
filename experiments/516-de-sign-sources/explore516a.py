"""516 exploration: production's sign by block (top-2,200 H1 control CPM) vs rest, agreement with CLEAN4, and the
oracle per-target choice among same-target sources (an upper bound for any per-target gate)."""
import numpy as np
from common516 import PANEL_DONORS, ROOT, cd4_rows, donor_rows, known_answers, load_h1

d = load_h1(); got, ok, _ = known_answers(d); assert ok; print("KA PASS", got)
T, G = d["targets"], d["genes"]; pool, band, truth, prod = d["pool"], d["band"], d["truth"], d["prod"]
ctrl = d["ctrl"]
blk = np.zeros(len(G), bool); blk[np.argsort(-ctrl)[:2200]] = True
z = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
pidx = {g: i for i, g in enumerate(z["genes"].astype(str))}
Pd = np.array([z["P_down"][pidx[g]] if g in pidx else np.nan for g in G])
c4 = np.where(np.isfinite(Pd), np.sign(0.5 - Pd), 0)
for nm, m in (("block", blk), ("rest", ~blk)):
    M = pool & m[None, :] & (prod != 0)
    acc = (prod[M] == truth[M]).mean()
    agree_c4 = ((prod == c4[None, :]) & M & (c4[None, :] != 0)).sum() / (M & (c4[None, :] != 0)).sum()
    print(f"{nm}: pool cells {M.sum()} ({M.sum()/pool.sum():.3f}), production acc {acc:.4f}, "
          f"prod agrees with CLEAN4 sign on {agree_c4:.3f}; truth down share {(truth[M] < 0).mean():.3f}, "
          f"prod down share {(prod[M] < 0).mean():.3f}")
    for b in range(3):
        Mb = M & (band == b)[:, None]
        if Mb.sum():
            print(f"    band {b}: n {Mb.sum()}, prod acc {(prod[Mb] == truth[Mb]).mean():.4f}")
rows = {dn: donor_rows(dn, T, G)[0] for dn in PANEL_DONORS}
cd4 = cd4_rows(T, G, layers=("log_fc",))
for c, r in cd4.items():
    rows["cd4_" + c] = r["log_fc"]
pan = np.stack([rows[k] for k in PANEL_DONORS])
rows["cons3"] = np.nanmean(pan, 0)
for nm, m in (("block", blk), ("rest", ~blk)):
    print(f"-- {nm}: covered-cell accuracy, production vs source")
    for k, R in rows.items():
        S = np.sign(np.nan_to_num(R)); M = pool & m[None, :] & (prod != 0) & (S != 0)
        print(f"   {k:14s} n {M.sum():6d}  prod {(prod[M] == truth[M]).mean():.4f}  src {(S[M] == truth[M]).mean():.4f}")
# oracle per-target choice among production and sources (hybrid with production fallback)
cands = {"prod": prod}
for k, R in rows.items():
    S = np.sign(np.nan_to_num(R)); cands[k] = np.where(S != 0, S, prod)
M = pool & (prod != 0)
hits = {k: ((C == truth) & M).sum(1) for k, C in cands.items()}
n = M.sum(1)
best = np.max(np.stack(list(hits.values())), 0)
print(f"\noracle per-target best source: pooled acc {best.sum()/n.sum():.4f} vs production {hits['prod'].sum()/n.sum():.4f}")
for b in range(3):
    g = band == b
    print(f"   band {b}: oracle {best[g].sum()/n[g].sum():.4f} vs prod {hits['prod'][g].sum()/n[g].sum():.4f}")
# per-target accuracy table of sources for the >=100 band
print("\nper-target accuracy, targets with n_real >= 100 (sorted by n_real):")
idx = np.argsort(-d["n_real"])
print(f"{'target':10s} {'n_real':>6s} " + " ".join(f"{k[:9]:>9s}" for k in cands))
for k in idx:
    if band[k] != 2: continue
    print(f"{T[k]:10s} {d['n_real'][k]:6d} " + " ".join(f"{hits[c][k]/max(n[k],1):9.3f}" for c in cands))
