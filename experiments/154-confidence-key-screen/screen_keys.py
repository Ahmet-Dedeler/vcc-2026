"""Does ANY available signal isolate a subset of our predictions with >=90% sign accuracy?

`reach` is the deepest prefix of the reference's DE pool, in OUR confidence order, whose
sign purity clears 0.90. The pool and each gene's correctness are fixed by the model; only
the order is free. So reach is unlocked iff some computable key puts a genuinely pure head
first. Ordering by |lfc| tops out near 0.736 in its best decile (peer session, exp 153).

This screens candidate keys DETERMINISTICALLY -- no emitter, no seeds, so the ~0.012 reach
seed-noise that invalidates single-seed phi arms cannot reach these numbers.

Leakage: the evaluation surface is the orion-hct116 harness, so donor `orion_hct116` is
EXCLUDED. Only donors the LODO prediction could legitimately see are used.

Caveat recorded up front: this harness's reference DE pool is low-power (n_real ~68-304)
versus the challenge (~2,670 by board inversion). Higher power admits weaker genes, where
our signs are worse -- so purity here is OPTIMISTIC. A negative result therefore closes the
channel robustly; a positive result needs re-checking in a higher-power regime.
"""
from __future__ import annotations
import sys
import numpy as np

ROOT = "~/Code/virtual-cell-challenge-2026"
ALPHA = 0.05
HARNESS_DONOR = "orion_hct116"          # excluded: it IS the evaluation surface
DONORS = ["k562_gwps", "orion_hek293t", "h1_2025", "k562_essential", "nadig_jurkat", "rpe1"]

d = np.load(f"{ROOT}/experiments/095-clean-harness/de_thinned.npz", allow_pickle=True)
genes = d["genes"].astype(str); perts = d["perts"].astype(str)
real_padj, real_lfc = d["real_padj"], d["real_lfc"]
pred_padj, pred_lfc = d["pred_padj"], d["pred_lfc"]
P, G = real_lfc.shape
print(f"{P} perturbations x {G} genes")

# ---- donor delta matrices projected onto the DE-table gene axis ----
donor_delta = {}
for name in DONORS:
    root = f"{ROOT}/data/donors/{name}"
    try:
        dg = np.array(open(f"{root}/genes.txt").read().split())
        dp = np.array(open(f"{root}/perts.txt").read().split())
        delta = np.load(f"{root}/delta.npy", mmap_mode="r")
    except FileNotFoundError:
        print(f"  {name}: missing, skipped"); continue
    gi = {g: i for i, g in enumerate(dg)}; cols = np.array([gi.get(g, -1) for g in genes])
    pi = {p: i for i, p in enumerate(dp)}
    M = np.zeros((P, G), dtype=np.float64); have = np.zeros(P, dtype=bool)
    ok = cols >= 0
    for j, p in enumerate(perts):
        r = pi.get(p)
        if r is None: continue
        row = np.asarray(delta[r], dtype=np.float64)
        M[j, ok] = row[cols[ok]]; have[j] = True
    donor_delta[name] = (M, have)
    print(f"  {name}: {have.sum()}/{P} perturbations, {ok.sum()}/{G} genes")

# ---- correctness on the reference's DE pool ----
pool = (real_padj < ALPHA) & np.isfinite(real_lfc) & (real_lfc != 0)
correct = np.sign(pred_lfc) == np.sign(real_lfc)
print(f"\npool: {int(pool.sum())} (pert,gene) pairs; overall sign accuracy {correct[pool].mean():.4f}")

# ---- candidate keys (higher = more confident) ----
stack = np.stack([donor_delta[n][0] for n in donor_delta])          # (D,P,G)
mask = np.stack([np.broadcast_to(donor_delta[n][1][:, None], (P, G)) for n in donor_delta])
sgn = np.sign(stack) * mask
n_avail = mask.sum(0).astype(np.float64)
net = np.abs(sgn.sum(0))
concord = np.divide(net, n_avail, out=np.zeros_like(net), where=n_avail > 0)
absmean = np.divide((np.abs(stack) * mask).sum(0), n_avail,
                    out=np.zeros_like(net), where=n_avail > 0)
geo = np.sqrt(np.maximum(absmean, 0.0))
agree_pred = np.divide(((np.sign(stack) == np.sign(pred_lfc)[None]) * mask).sum(0), n_avail,
                       out=np.zeros_like(net), where=n_avail > 0)

KEYS = {
    "|pred_lfc| (baseline)":            np.abs(pred_lfc),
    "-pred_padj":                       -pred_padj,
    "donor sign concordance":           concord,
    "concordance x |pred_lfc|":         concord * np.abs(pred_lfc),
    "geo-magnitude x concordance":      geo * concord,
    "agreement-with-pred":              agree_pred,
    "agreement x |pred_lfc|":           agree_pred * np.abs(pred_lfc),
    "n donors available":               n_avail,
    "concordance^2 x |pred_lfc|":       (concord ** 2) * np.abs(pred_lfc),
}

def segment_purity(key):
    """Pool all (pert,gene) pool members, rank by key desc, report purity of top fractions."""
    k = key[pool]; c = correct[pool]
    order = np.argsort(-k, kind="stable"); c = c[order]
    n = c.size; out = {}
    for frac in (0.01, 0.05, 0.10, 0.25, 0.50, 1.00):
        m = max(1, int(round(n * frac)))
        out[frac] = c[:m].mean()
    cum = np.cumsum(c) / np.arange(1, n + 1)
    # deepest prefix clearing 0.90, requiring at least 20 genes to avoid tiny-prefix flukes
    ok = np.flatnonzero(cum[19:] >= 0.90)
    out["depth90"] = int(ok[-1] + 20) if ok.size else 0
    out["n"] = n
    return out

print(f"\n{'key':<30} {'top1%':>7} {'top5%':>7} {'top10%':>7} {'top25%':>7} {'all':>7}  depth>=0.90")
print("-" * 82)
rows = []
for name, key in KEYS.items():
    r = segment_purity(np.nan_to_num(key, nan=-np.inf))
    rows.append((name, r))
    print(f"{name:<30} {r[0.01]:>7.4f} {r[0.05]:>7.4f} {r[0.10]:>7.4f} "
          f"{r[0.25]:>7.4f} {r[1.00]:>7.4f}  {r['depth90']:>6d} / {r['n']}")

best = max(rows, key=lambda x: x[1][0.01])
print(f"\nbest top-1% purity: {best[0]} at {best[1][0.01]:.4f}")
print("VERDICT:", "SOME KEY CLEARS 0.90" if best[1][0.01] >= 0.90 else
      "no key reaches 0.90 in even its top 1% -- reach is model-limited, not order-limited")


# ============================================================================
# PER-PERTURBATION reach -- this is what the scorer actually computes.
# Pooled purity above is an easier problem (it may cherry-pick across targets).
# ============================================================================
def reach_per_pert(key, min_pool=10):
    depths, sizes, ks = [], [], []
    for j in range(P):
        idx = np.flatnonzero(pool[j])
        if idx.size < min_pool: continue
        k = np.nan_to_num(key[j, idx], nan=-np.inf)
        order = idx[np.argsort(-k, kind="stable")]
        c = correct[j, order]
        cum = np.cumsum(c) / np.arange(1, c.size + 1)
        ok = np.flatnonzero(cum >= 0.90)
        kstar = int(ok[-1] + 1) if ok.size else 0
        depths.append(kstar / c.size); sizes.append(c.size); ks.append(kstar)
    return (float(np.mean(depths)), int(np.median(ks)),
            int(np.sum(np.array(ks) == 0)), len(depths), float(np.median(sizes)))

# oracle: correct signs first
oracle_key = correct.astype(np.float64)

print("\n" + "=" * 82)
print("PER-PERTURBATION reach (the scorer's definition)")
print("=" * 82)
print(f"{'key':<30} {'reach':>8} {'median k*':>10} {'k*=0':>7} {'perts':>6}")
print("-" * 66)
for name, key in list(KEYS.items()) + [("ORACLE (correct first)", oracle_key)]:
    r, medk, zeros, n, medpool = reach_per_pert(key)
    print(f"{name:<30} {r:>8.5f} {medk:>10d} {zeros:>7d} {n:>6d}")
base = reach_per_pert(KEYS["|pred_lfc| (baseline)"])[0]
best_name, best_r = None, -1
for name, key in KEYS.items():
    r = reach_per_pert(key)[0]
    if r > best_r: best_name, best_r = name, r
orc = reach_per_pert(oracle_key)[0]
print(f"\nmedian pool size: {reach_per_pert(oracle_key)[4]:.0f}")
print(f"baseline |pred_lfc| reach {base:.5f}")
print(f"best key: {best_name} at {best_r:.5f}  (delta {best_r-base:+.5f})")
print(f"oracle {orc:.5f}  -> best key captures {(best_r-base)/max(orc-base,1e-9)*100:.1f}% of the headroom")
