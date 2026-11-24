"""Q29 — what base rate does `fid` precision actually run against?

The concern: a21 calls 85% UP against a reference pool that 153 section 2 measured as 59.7% DOWN,
so precision should suffer.  But section 2's 59.7% was over the reference's **DE pool**.  `fid`'s
k counts our called genes against the reference's lfc SIGN whether or not the reference calls
them significant (BRIEF section 4: the call set is our significant genes that are *adjudicable*,
i.e. the reference assigns a defined non-zero lfc).  Since n_pred >> n_real, most of our calls
are reference-NON-significant genes -- and 094 measured those drifting UP.

So the base rate that matters is the reference lfc sign over the genes WE call, not over the
reference's DE set.  Those are different populations and they point opposite ways.
"""
import numpy as np
from debench import Surface

for name, path in (("HCT116", "../095-clean-harness/de_thinned.npz"),
                   ("HEK293T", "de_hek293t_ref.npz")):
    s = Surface(path)
    up_all = up_de = up_nonde = n_all = n_de = n_non = 0
    for p in range(len(s.perts)):
        adj = s.adj[p] & np.isfinite(s.real_lfc[p])
        de = adj & s.real_sig[p]
        non = adj & ~s.real_sig[p]
        up_all += int((s.real_lfc[p, adj] > 0).sum()); n_all += int(adj.sum())
        up_de += int((s.real_lfc[p, de] > 0).sum()); n_de += int(de.sum())
        up_nonde += int((s.real_lfc[p, non] > 0).sum()); n_non += int(non.sum())
    q_all, q_de, q_non = up_all/n_all, up_de/max(n_de,1), up_nonde/max(n_non,1)
    print(f"{name}: reference lfc UP-fraction")
    print(f"   over ALL adjudicable gated genes   {q_all:.4f}   (n={n_all:,})")
    print(f"   over the reference's DE set         {q_de:.4f}   (n={n_de:,})   <- 153 sec2's population")
    print(f"   over reference NON-significant      {q_non:.4f}   (n={n_non:,})  <- where our calls mostly live")
    print(f"   expected fid precision at our up-rate u, no skill:  u*q + (1-u)*(1-q), q={q_all:.4f}")
    for u in (0.62, 0.70, 0.85):
        print(f"       u={u:.2f} -> {u*q_all + (1-u)*(1-q_all):.4f}")
    print()
