"""Map every arm's harness-measured (M, delta-eps) onto the challenge panel, exactly.

The challenge's own shape constants come from a shipped submission plus the board-anchored solve
in MSE-MECHANISM.md: A, B, Kc are computed from our own file, X and sumD were solved from two
board points.  So `eps` and `kappa` have CHALLENGE levels; only their CHANGE is transferred from
the harnesses, and the change is taken pessimistically (the worse of the two surfaces on each).
"""
import numpy as np
A, B, X = 0.282948, 0.280503, 0.018846
Kc, sD, CLAIM = 0.005410, 0.228572, 0.24653          # claim at the #247 cap = a12's policy
eps0 = (A - B) / B
kap0 = X / np.sqrt(B * sD)
BA, SPAN = 0.987863, 0.950163
A11_SC, A11_OVER = 0.082041, 0.166115


def raw(M, deps):
    return 1 + (Kc + CLAIM * (eps0 + deps) - 2 * np.sqrt(CLAIM * sD) * kap0 * M) / sD


def sc(r):
    return max(0.0, min(1.0, (BA - r) / SPAN))


print(f"challenge shape: eps0={eps0:.6f}  kappa0={kap0:.6f}  claim/B={CLAIM/B:.4f}")
print(f"a12 base (M=1, deps=0): raw={raw(1,0):.4f}  (MSE-MECHANISM quotes 0.8785)\n")
# (label, M_h1, M_ori, eps_h1, eps_ori)  eps values are the ARM's, base subtracted below
E0H, E0O = 0.01424, 0.01208
ARMS = [
    ("a12 as built (global alpha only)", 1.000, 1.000, E0H, E0O),
    ("(b) per-pert donor_cos^1", 1.208, 1.165, 0.01385, 0.01125),
    ("(b) per-pert donor_cos^2 (no floor)", 1.238, 1.202, 0.01339, 0.01035),
    ("(b) per-pert donor_cos^2, floor 0.3", 1.261, 1.220, 0.01344, 0.01042),
    ("(b) per-pert donor_cos^3", 1.240, 1.205, 0.01310, 0.00978),
    ("(096 #5) per-gene agreement^1 -> ZERO [midrank]", 1.101, 1.525, 0.01317, 0.00960),
    ("(096 #5) per-gene agreement^1 -> MEAN", 1.178, 1.329, 0.03072, 0.02245),
    ("(b+5) STACK pert^2 x gene^1, no floor", 1.467, 1.560, 0.01273, 0.00877),
    ("(b+5) STACK pert^2 fl0.3 x gene^1 [SHIP]", 1.478, 1.610, 0.01278, 0.00887),
    ("(b+5) STACK pert^3 fl0.5 x gene^1 (edge)", 1.500, 1.658, 0.01274, 0.00879),


    ("ORACLE per-perturbation z=x/a", 2.012, 1.967, 0.01257, 0.00931),
    ("ORACLE per-gene slope s_g (unbuildable)", 4.098, 7.236, 0.08303, 0.23598),
]
base_sc = sc(raw(1, 0))
print(f"{'arm':<46}{'M use':>7}{'d-eps':>9}{'raw mse':>9}{'score_mse':>10}"
      f"{'overall':>9}{'d vs a12':>10}{'d vs a11':>10}")
for lbl, m1, m2, e1, e2 in ARMS:
    M = min(m1, m2)
    de = 0.0 if not np.isfinite(e1) else max(e1 - E0H, e2 - E0O)
    r = raw(M, de); s = sc(r)
    print(f"{lbl:<46}{M:>7.3f}{de:>+9.5f}{r:>9.4f}{s:>10.4f}{s/6:>9.4f}"
          f"{(s-base_sc)/6:>+10.4f}{(s-A11_SC)/6:>+10.4f}")
print(f"\nboard context: a11 = 0.166115 (rank 7).  leader Aginglab.com 0.21122 with raw mse "
      f"0.82698.\nthe leader's mse member is beaten at M >= "
      f"{(Kc + CLAIM*eps0 - (0.82698-1)*sD) / (2*np.sqrt(CLAIM*sD)*kap0):.3f}")
