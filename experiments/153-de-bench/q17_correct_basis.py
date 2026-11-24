"""Q17 — the (gamma, lam) grid on the operator each surface's baseline was ACTUALLY built with.

`104/work/<surface>/manifest.json` records the operator per surface:
    hct116  -> 062/.../held_orion_hct116.npz   (sha 7cecbf7f)
    hek293t -> 062/.../held_orion_hek293t.npz  (sha 97e486e0)
Both are LEAKAGE-SAFE: the evaluated line is held out of the basis.  Sections 9-12 used the
095 operator for both surfaces.  For HCT116 that is the same fit (same donors, max eigenvalue
differs by 1.3e-5 of 8.036, and W = V f(ev) V^T is invariant to eigenvector sign convention),
so HCT116 stands.  For HEK293T it was the wrong operator and this supersedes it.

NOTE on which operator is "right": `data/donors/whitener_union{,_pruned}.npz` -- what production
ships -- lists donors ['k562_gwps','orion_hct116','orion_hek293t'], so it CONTAINS both
evaluation lines.  Using it to score HCT116 or HEK293T truth leaks.  Section 12's "production
operator" run was therefore leaky, and this file supersedes it too.  The held-safe operator is
the only valid one for measurement; the shipped operator is a different object and the gap
between them is a limit of this measurement, not something to paper over.
"""
import numpy as np, csv
from debench import Surface
from exact_pds import pds

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
wcols = np.load("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
                "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")
wgenes = [canon[i] for i in wcols]
SPAN = {"pds": 0.452, "reach": 0.885, "nmae": 0.60}
OP = "../062-cd4-post-whiten/cloud-output/whiteners/held_orion_%s.npz"

def build(de_npz, froz, real_npz, opname):
    W = np.load(OP % opname); ev, V = W["evals"], W["evecs"]
    z = np.load(froz, allow_pickle=True)
    g104 = {str(g): i for i, g in enumerate(z["genes"])}
    base = z["baseline"][:, [g104[g] for g in wgenes]].astype(np.float32)
    t104 = [str(t) for t in z["targets"]]
    R = np.load(real_npz, allow_pickle=True)
    rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
    pgi = np.array([i for i, g in enumerate(wgenes) if g in rg and g not in panel])
    pcol = np.array([rg[wgenes[i]] for i in pgi])
    prow_r = np.array([rt.get(t, -1) for t in t104]); keep_t = prow_r >= 0
    Dreal = R["delta"][np.ix_(prow_r[keep_t], pcol)].astype(np.float64)
    s = Surface(de_npz)
    gpos = {g: i for i, g in enumerate(s.genes)}; tpos = {t: i for i, t in enumerate(t104)}
    gsel = np.array([i for i, g in enumerate(wgenes) if g in gpos])
    gdst = np.array([gpos[wgenes[i]] for i in gsel])
    prow = np.array([tpos.get(p, -1) for p in s.perts]); pok = prow >= 0
    kp = np.zeros_like(s.real_sig); kp[np.ix_(pok, gdst)] = True
    s.real_sig &= kp; s.adj &= kp
    U_q = (base @ V) * (ev ** (0.65 / 2.0)); med0 = [None]
    def score(gamma, lam):
        P = ((U_q * ((ev + lam) ** (-gamma / 2.0))).astype(np.float32) @ V.T)
        m = float(np.median(np.linalg.norm(P, axis=1)))
        if med0[0] is None: med0[0] = m
        P = P * (med0[0] / m)
        F = np.zeros_like(s.pred_lfc); F[np.ix_(pok, gdst)] = P[np.ix_(prow[pok], gsel)]
        s.pred_lfc = F
        return (pds(P[np.ix_(keep_t, pgi)].astype(np.float64), Dreal),
                s.reach(key=[np.abs(F[p]) for p in range(len(s.perts))]), s.nmae())
    return score

A = build("../095-clean-harness/de_thinned.npz",
          "../104-context-heldout-credit-adapter/work/hct116/frozen.npz",
          "real_delta_hct116.npz", "hct116")
B = build("de_hek293t_ref.npz",
          "../104-context-heldout-credit-adapter/work/hek293t/frozen.npz",
          "real_delta_hek293t.npz", "hek293t")
bA, bB = A(0.65, 0.0), B(0.65, 0.0)
def comb(x, b):
    return ((x[0]-b[0])/SPAN["pds"] + (x[1]-b[1])/SPAN["reach"] - (x[2]-b[2])/SPAN["nmae"]) / 6.0
rows = []
for gamma in (0.70, 0.80, 0.90, 1.00, 1.10):
    for lam in (0.005, 0.007, 0.010, 0.015, 0.020, 0.030):
        a, b = A(gamma, lam), B(gamma, lam)
        ok = a[0] > bA[0] and b[0] > bB[0] and a[1] > bA[1] and b[1] > bB[1] and a[2] < bA[2] and b[2] < bB[2]
        rows.append((min(comb(a, bA), comb(b, bB)), gamma, lam, a, b, ok))
rows.sort(key=lambda r: -r[0])
print(f"{'gamma':>6} {'lam':>7} {'HCTΔpds':>9} {'HCTΔreach':>10} {'HEKΔpds':>9} {'HEKΔreach':>10} {'WORST':>9}  all6+")
for w, gamma, lam, a, b, ok in rows[:12]:
    print(f"{gamma:6.2f} {lam:7.3f} {a[0]-bA[0]:+9.5f} {a[1]-bA[1]:+10.5f} {b[0]-bB[0]:+9.5f} "
          f"{b[1]-bB[1]:+10.5f} {w:+9.5f}  {'YES' if ok else 'no'}")
print(f"\nbaselines: HCT pds {bA[0]:.5f} reach {bA[1]:.5f} nmae {bA[2]:.5f}")
print(f"           HEK pds {bB[0]:.5f} reach {bB[1]:.5f} nmae {bB[2]:.5f}")
