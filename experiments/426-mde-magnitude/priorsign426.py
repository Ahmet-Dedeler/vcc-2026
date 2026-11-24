"""426b analytic: which sign should the MDE add use on genes we do not call? Our delivered LFC sign (the bulk profile's)
vs the CLEAN4 per-gene prior (371's consensus P_down over K562-GWPS, K562-ess, Nadig Jurkat, Nadig HepG2), on H1's gates.
Base = 426's base_s0_bps delivered stats (the 410b construction). MDE add = min(0.1/sqrt(count), 0.1), 0 below 0.3."""
import sys
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
W323 = ROOT / "experiments/323-current-family-validation/work"
r = np.load(W323 / "reference.npz", allow_pickle=True)
keep = r["keep"].astype(bool)
genes_full = r["genes"].astype(str); genes = genes_full[keep]
s = np.load(Path(__file__).parent / "work/h1/stats-base_s0_bps.npz")
targets = s["targets"].astype(str)
L = np.nan_to_num(s["lfc"]); Q = s["q"]
count = r["ref_lin"] * np.median(r["libraries"]) / 1e6
add = np.where(count >= 0.3, np.minimum(0.1 / np.sqrt(np.maximum(count, 1e-3)), 0.1), 0.0)
cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float)
P = PB[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
col = np.array([pos.get(g, -1) for g in genes])
Pd = np.where(col >= 0, P[np.maximum(col, 0)], np.nan)          # P_down per gene on the DE axis
prior = np.where(np.isfinite(Pd), np.where(Pd > 0.5, -1.0, 1.0), 0.0)
conf = np.where(np.isfinite(Pd), np.abs(Pd - 0.5), 0.0)
ref = pl.read_parquet(W323 / "reference_de.parquet"); gi = {g: i for i, g in enumerate(genes)}
gates = []
for k, t in enumerate(targets):
    m = ref.filter(pl.col("target") == t)
    g = m["feature"].to_numpy().astype(str); q = m["p_adj"].to_numpy().astype(float); l = m["log2_fold_change"].to_numpy().astype(float)
    ok = (q < .05) & np.isfinite(l) & (l != 0) & (g != t)
    if ok.sum() >= 10:
        gates.append((k, np.array([gi[x] for x in g[ok]]), l[ok]))
def nmae(Lh):
    return float(np.mean([np.abs(Lh[k, idx] - l).sum() / np.abs(l).sum() for k, idx, l in gates]))
K = np.concatenate([np.full(i.size, k) for k, i, _ in gates]); I = np.concatenate([i for _, i, _ in gates]); LL = np.concatenate([l for _, _, l in gates])
called = Q[K, I] < 0.05
ours = np.sign(L[K, I]); pr = prior[I]
for lab, msk in (("called", called), ("uncalled", ~called)):
    cov = pr[msk] != 0
    print(f"{lab:8s}: n {msk.sum():6d}  sign acc ours {np.mean(ours[msk] == np.sign(LL[msk])):.4f}  prior {np.mean(pr[msk][cov] == np.sign(LL[msk][cov])):.4f} (covered {cov.mean():.3f})"
          f"  agree(ours, prior) {np.mean(ours[msk][cov] == pr[msk][cov]):.3f}  acc where agree {np.mean(ours[msk][cov][ours[msk][cov] == pr[msk][cov]] == np.sign(LL[msk][cov][ours[msk][cov] == pr[msk][cov]])):.4f}")
base = nmae(L)
sg = np.sign(L)
uncalled_full = Q >= 0.05
def arm(sign_uncalled):
    S = np.where(uncalled_full, sign_uncalled, sg)
    return nmae(L + S * add[None, :]) - base
print(f"\nbase nmae {base:.5f}")
print(f"MDE, our signs everywhere:                 {arm(sg):+.5f}")
pfull = np.broadcast_to(prior[None, :], L.shape)
print(f"MDE, prior sign on uncalled (covered):      {arm(np.where(pfull != 0, pfull, sg)):+.5f}")
for t in (0.05, 0.1, 0.15, 0.2):
    cf = np.broadcast_to(conf[None, :] > t, L.shape)
    print(f"MDE, prior on uncalled where |P-.5|>{t:.2f}:   {arm(np.where(cf & (pfull != 0), pfull, sg)):+.5f}  (share {cf[0].mean():.3f})")
print(f"MDE only where ours == prior (uncalled):    {nmae(L + np.where(uncalled_full & (sg != pfull), 0.0, sg) * add[None, :]) - base:+.5f}")

# 426e: drop the add (do not flip) where our uncalled sign disagrees with a covered prior; keep it where they agree or the
# prior is silent. No flips -> no extra class mass.
keep_add = ~(uncalled_full & (pfull != 0) & (sg != pfull))
print(f"MDE, add dropped where uncalled ours != covered prior: {nmae(L + sg * add[None, :] * keep_add) - base:+.5f}  "
      f"(drops {(~keep_add).mean():.3f} of gene-target cells)")
for t in (0.02, 0.05, 0.1):
    cf = np.broadcast_to(conf[None, :] > t, L.shape)
    ka = ~(uncalled_full & (pfull != 0) & cf & (sg != pfull))
    print(f"   ... only where the prior is confident |P-.5|>{t}: {nmae(L + sg * add[None, :] * ka) - base:+.5f}")
# and for called genes too (block genes carry the prior already, head genes carry the design)
ka_all = ~((pfull != 0) & (sg != pfull))
print(f"MDE, add dropped wherever ours != covered prior (called too): {nmae(L + sg * add[None, :] * ka_all) - base:+.5f}")
