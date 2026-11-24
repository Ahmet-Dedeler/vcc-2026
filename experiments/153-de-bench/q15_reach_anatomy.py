"""Q15 — is the lam arm's reach gain REAL, or the consensus-direction artefact?

The 155 subagent found that raw `reach` is not monotone in direction quality on this surface:
an arm with WORSE head purity (0.6093 vs 0.6646) posted nearly 2x the reach, via a heavy tail
of long >=90%-pure prefixes on large-n_real perturbations -- mean k* 14.8 vs 5.0 but MEDIAN
k* 1 vs 3, and 36% of perturbations at zero.  A target-independent consensus direction can buy
reach that way.

If tonight's ship candidate gains reach through that channel, the gain is not real quality and
should not be shipped.  Signature to look for: mean k* up but median k* flat or down, mass
moving into the large-n_real tail, and k*=0 count NOT improving.
"""
import numpy as np, csv
from debench import Surface

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
wcols = np.load("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
                "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")
wgenes = [canon[i] for i in wcols]
W = np.load("../095-clean-harness/out/whitener_union.npz"); ev, V = W["evals"], W["evecs"]

def anatomy(de_npz, froz, label):
    z = np.load(froz, allow_pickle=True)
    g104 = {str(g): i for i, g in enumerate(z["genes"])}
    base = z["baseline"][:, [g104[g] for g in wgenes]].astype(np.float32)
    t104 = [str(t) for t in z["targets"]]
    s = Surface(de_npz)
    gpos = {g: i for i, g in enumerate(s.genes)}; tpos = {t: i for i, t in enumerate(t104)}
    gsel = np.array([i for i, g in enumerate(wgenes) if g in gpos])
    gdst = np.array([gpos[wgenes[i]] for i in gsel])
    prow = np.array([tpos.get(p, -1) for p in s.perts]); pok = prow >= 0
    kp = np.zeros_like(s.real_sig); kp[np.ix_(pok, gdst)] = True
    s.real_sig &= kp; s.adj &= kp
    U_q = (base @ V) * (ev ** (0.65 / 2.0)); med0 = [None]

    def ks(gamma, lam):
        P = ((U_q * ((ev + lam) ** (-gamma / 2.0))).astype(np.float32) @ V.T)
        m = float(np.median(np.linalg.norm(P, axis=1)))
        if med0[0] is None: med0[0] = m
        P = P * (med0[0] / m)
        F = np.zeros_like(s.pred_lfc); F[np.ix_(pok, gdst)] = P[np.ix_(prow[pok], gsel)]
        out, nre = [], []
        for p in range(len(s.perts)):
            pool = np.flatnonzero(s.real_sig[p])
            if pool.size == 0: continue
            o = np.lexsort((s.genes[pool], -np.abs(F[p, pool])))
            q = pool[o]; q = q[s.adj[p, q]]
            if q.size == 0: out.append(0); nre.append(int(s.real_sig[p].sum())); continue
            ok = np.sign(F[p, q]) == np.sign(s.real_lfc[p, q])
            pur = np.cumsum(ok) / np.arange(1, q.size + 1)
            d = np.flatnonzero(pur >= 0.9)
            out.append(int(d[-1] + 1) if d.size else 0); nre.append(int(s.real_sig[p].sum()))
        return np.array(out), np.array(nre)

    print(f"\n=== {label} ===")
    print(f"{'arm':22} {'reach':>8} {'mean k*':>8} {'med k*':>7} {'k*=0':>6} "
          f"{'reach|small':>12} {'reach|mid':>10} {'reach|large':>12}")
    for name, (g, l) in (("production", (0.65, 0.0)), ("SHIP 0.80/0.007", (0.80, 0.007)),
                         ("1.00/0.020", (1.00, 0.020))):
        k, nr = ks(g, l)
        r = k / np.maximum(nr, 1)
        t1, t2 = np.percentile(nr, [33, 67])
        lo, mid, hi = nr <= t1, (nr > t1) & (nr <= t2), nr > t2
        print(f"{name:22} {r.mean():8.5f} {k.mean():8.2f} {np.median(k):7.1f} "
              f"{(k==0).sum():6d} {r[lo].mean():12.5f} {r[mid].mean():10.5f} {r[hi].mean():12.5f}")

anatomy("../095-clean-harness/de_thinned.npz",
        "../104-context-heldout-credit-adapter/work/hct116/frozen.npz", "HCT116")
anatomy("de_hek293t_ref.npz",
        "../104-context-heldout-credit-adapter/work/hek293t/frozen.npz", "HEK293T")
