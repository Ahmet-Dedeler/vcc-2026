"""Q32 — the cross-context transfer CEILING, with zero estimation noise on either side.

~40 experiments have tried to improve per-gene direction by adding or reweighting donor sources,
and all failed.  Nobody has measured the thing that bounds them all:

    If you had a donor line's TRUE effect for a target -- no estimation noise, no whitening, no
    blending -- how well would its per-gene SIGN predict the same target's sign in a DIFFERENT
    cell line?

That is the ceiling for the entire donor-transfer paradigm on one donor.  We hold both sides
measured: `de_thinned.npz` is the HCT116 reference and `de_hek293t_ref.npz` the HEK293T
reference, on the SAME 132 VCC panel targets.  Both are real measured perturbation responses,
so this is truth-to-truth.

If the ceiling is ~0.60, our 0.52 is near it and better estimation cannot close the gap -- the
paradigm is capped and the remaining score needs a different mechanism.
If it is ~0.85, we are losing ~0.33 to estimation and denoising/more data is the whole game.
"""
import numpy as np
from debench import Surface

hct = Surface("../095-clean-harness/de_thinned.npz")
hek = Surface("de_hek293t_ref.npz")
gh = {g: i for i, g in enumerate(hek.genes)}
th = {t: i for i, t in enumerate(hek.perts)}
gcol = np.array([gh.get(g, -1) for g in hct.genes])
prow = np.array([th.get(t, -1) for t in hct.perts])
gok, pok = gcol >= 0, prow >= 0
print(f"shared: {pok.sum()}/{hct.perts.size} targets, {gok.sum()}/{hct.genes.size} genes\n")

num = den = 0                      # truth-to-truth, on HCT's DE pool
num_w = den_w = 0.0                # |lfc|-weighted
by_dec = [[0, 0] for _ in range(5)]
prod_num = prod_den = 0
for p in range(hct.perts.size):
    if not pok[p]:
        continue
    pool = np.flatnonzero(hct.real_sig[p] & hct.adj[p] & gok)
    if pool.size < 10:
        continue
    a = np.sign(hct.real_lfc[p, pool])                       # truth in HCT116
    b = np.sign(hek.real_lfc[prow[p], gcol[pool]])           # truth in HEK293T
    ok = (a == b) & (b != 0)
    num += int(ok.sum()); den += int((b != 0).sum())
    w = np.abs(hct.real_lfc[p, pool])
    num_w += float(w[ok].sum()); den_w += float(w[b != 0].sum())
    mag = np.abs(hek.real_lfc[prow[p], gcol[pool]])          # stratify by DONOR magnitude
    qs = np.percentile(mag, [20, 40, 60, 80])
    idx = np.digitize(mag, qs)
    for d in range(5):
        m = idx == d
        by_dec[d][0] += int(((a == b) & (b != 0) & m).sum()); by_dec[d][1] += int(((b != 0) & m).sum())
    ps = np.sign(hct.pred_lfc[p, pool])                      # our shipped prediction
    prod_num += int((ps == a).sum()); prod_den += pool.size

print(f"TRUTH-TO-TRUTH cross-context sign agreement (HEK293T truth -> HCT116 truth)")
print(f"   unweighted            {num/den:.4f}   (n={den:,})")
print(f"   weighted by |lfc_HCT| {num_w/den_w:.4f}")
print(f"   by DONOR |lfc| quintile (5 = donor's strongest effects):")
for d in range(5):
    n_, dn = by_dec[d]
    print(f"      Q{d+1}  {n_/max(dn,1):.4f}   (n={dn:,})")
print(f"\nOUR SHIPPED PREDICTION on the same pool: {prod_num/prod_den:.4f}   (n={prod_den:,})")
print(f"replicate anchor for reference: `reach` needs 90% purity; `fid` replicate is 0.795-0.832")
