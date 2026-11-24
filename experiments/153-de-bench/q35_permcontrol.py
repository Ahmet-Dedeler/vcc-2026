"""Q35 — is the 0.7038 real cross-context transfer, or a shared per-gene artefact?

Both harnesses come from Orion via the same builder with the same control-thinning, so a shared
systematic (composition drift, batch, depth) would make their per-gene lfc signs agree for
NON-biological reasons.  `docs/HARNESS-DEPTH-CONFOUND.md` documents exactly this class of error
costing the project six experiments.

Control: predict HCT116 target X's signs from HEK293T target **Y != X**.  Target-specific
transfer should collapse to chance; a shared per-gene artefact should not move at all.

This is the check that decides whether Q34's "0.13 lost in the donor table" is a finding or my
own version of the depth confound.
"""
import numpy as np
from debench import Surface

hct = Surface("../095-clean-harness/de_thinned.npz")
hek = Surface("de_hek293t_ref.npz")
gh = {g: i for i, g in enumerate(hek.genes)}; th = {t: i for i, t in enumerate(hek.perts)}
tgt = [p for p, t in enumerate(hct.perts) if t in th]
cg = [i for i, g in enumerate(hct.genes) if g in gh]
rng = np.random.default_rng(0)

def run(perm):
    num = den = 0
    for k, p in enumerate(tgt):
        pool = np.array([g for g in cg if hct.real_sig[p, g] and hct.adj[p, g]])
        if pool.size < 10: continue
        q = th[hct.perts[tgt[perm[k]]]]
        b = np.sign(hek.real_lfc[q, [gh[hct.genes[g]] for g in pool]])
        a = np.sign(hct.real_lfc[p, pool])
        num += int(((a == b) & (b != 0)).sum()); den += int((b != 0).sum())
    return num / max(den, 1)

ident = np.arange(len(tgt))
real = run(ident)
nulls = []
for s in range(20):
    perm = rng.permutation(len(tgt))
    while np.any(perm == ident):                      # no target may keep its own partner
        perm = rng.permutation(len(tgt))
    nulls.append(run(perm))
nulls = np.array(nulls)
print(f"matched target      {real:.4f}")
print(f"MISMATCHED target   {nulls.mean():.4f}  (20 permutations, sd {nulls.std():.4f}, "
      f"range {nulls.min():.4f}-{nulls.max():.4f})")
print(f"target-specific excess: {real - nulls.mean():+.4f}")
print(f"\nchance is 0.50.  A shared per-gene artefact would keep the mismatched value HIGH.")
