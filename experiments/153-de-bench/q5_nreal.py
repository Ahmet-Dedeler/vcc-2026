"""Q5 — is n_real on the challenge ~240 or ~2,000?  The two board constraints disagree.

Constraint A (094 A8): solving a8-vs-v11's jac CHANGE gives n_real ~ 2,000.
Constraint B (this file): the board's jac LEVEL is 0.0301.  If our call set on the challenge
behaved like it does on the repaired harness, n_real ~ 2,000 would put board jac at ~0.175.

They can only both hold if our challenge recall is far worse than our harness recall.  The
discriminator is whether recall is stable in n_real -- measurable here, for free.
"""
import numpy as np
from debench import Surface

print(f"{'regime':>12} {'n_real':>8} {'n_pred':>8} {'|R n R^|':>9} {'recall':>8} {'jac':>8} {'enrich vs random':>17}")
for tgt in (None, 250, 500, 1000, 2000, 3000):
    s = Surface("../095-clean-harness/de_thinned.npz", n_real_target=tgt)
    inter = rec = jc = rnd = 0.0
    n = 0
    for p in range(len(s.perts)):
        a, b = s.real_sig[p], s.pred_sig[p] & s.adj[p]
        nr, npd = int(a.sum()), int(b.sum())
        if nr == 0 or npd == 0:
            continue
        i = int((a & b).sum())
        inter += i; rec += i / nr; rnd += nr * npd / s.genes.size; n += 1
    s2 = Surface("../095-clean-harness/de_thinned.npz", n_real_target=tgt)
    lab = "native" if tgt is None else f"top-{tgt}"
    print(f"{lab:>12} {np.median(s.counts()[0]):8.0f} {np.median(s.counts()[1]):8.0f} "
          f"{inter/n:9.1f} {rec/n:8.4f} {s2.jac():8.5f} {inter/rnd:17.3f}")

print("\nBoard: n_pred ~5,961 (v11), jac 0.0321, ~9,900 gated genes.")
print(f"{'assumed n_real':>15} {'implied |inter|':>15} {'implied recall':>15} {'random |inter|':>15} {'enrichment':>11}")
for N in (240, 500, 1000, 2000, 4000):
    inter = 0.0321 * (5961 + N) / 1.0321
    print(f"{N:15d} {inter:15.1f} {inter/N:15.3f} {N*5961/9900:15.1f} {inter/(N*5961/9900):11.3f}")
