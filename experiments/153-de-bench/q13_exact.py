"""Q13 — the (gamma, lam) screen with EXACT scorer-space pds, on both surfaces."""
import numpy as np, csv, sys
from debench import Surface, overall_delta
from exact_pds import pds

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
wcols = np.load("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
                "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")
wgenes = [canon[i] for i in wcols]
W = np.load("../095-clean-harness/out/whitener_union.npz"); ev, V = W["evals"], W["evecs"]

def run(tag, de_npz, froz, real_npz):
    z = np.load(froz, allow_pickle=True)
    g104 = {str(g): i for i, g in enumerate(z["genes"])}
    base = z["baseline"][:, [g104[g] for g in wgenes]].astype(np.float32)
    t104 = [str(t) for t in z["targets"]]
    R = np.load(real_npz, allow_pickle=True)
    rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
    # exact-pds axis: whitener genes present in the real table, panel targets removed (BRIEF 1)
    pg = [i for i, g in enumerate(wgenes) if g in rg and g not in panel]
    pcol = np.array([rg[wgenes[i]] for i in pg]); pg = np.array(pg)
    prow_r = np.array([rt.get(t, -1) for t in t104]); keep_t = prow_r >= 0
    Dreal = R["delta"][np.ix_(prow_r[keep_t], pcol)].astype(np.float64)

    s = Surface(de_npz)
    gpos = {g: i for i, g in enumerate(s.genes)}; tpos = {t: i for i, t in enumerate(t104)}
    gsel = np.array([i for i, g in enumerate(wgenes) if g in gpos])
    gdst = np.array([gpos[wgenes[i]] for i in gsel])
    prow = np.array([tpos.get(p, -1) for p in s.perts]); pok = prow >= 0
    kp = np.zeros_like(s.real_sig); kp[np.ix_(pok, gdst)] = True
    s.real_sig &= kp; s.adj &= kp
    U_q = (base @ V) * (ev ** (0.65 / 2.0))
    med0 = [None]
    def score(gamma, lam):
        P = ((U_q * ((ev + lam) ** (-gamma / 2.0))).astype(np.float32) @ V.T)
        med = float(np.median(np.linalg.norm(P, axis=1)))
        if med0[0] is None: med0[0] = med
        P = P * (med0[0] / med)
        F = np.zeros_like(s.pred_lfc); F[np.ix_(pok, gdst)] = P[np.ix_(prow[pok], gsel)]
        s.pred_lfc = F
        r = s.reach(key=[np.abs(F[p]) for p in range(len(s.perts))])
        return pds(P[np.ix_(keep_t, pg)].astype(np.float64), Dreal), r, s.nmae()
    print(f"\n=== {tag} ===")
    print(f"{'gamma':>6} {'lam':>7} {'pds_EXACT':>10} {'Δpds':>9} {'reach':>9} {'Δreach':>9} {'nmae':>8} {'Δnmae':>9}")
    b = None
    for gamma, lam in [(0.65,0.0),(0.65,0.007),(0.65,0.020),(0.80,0.007),(0.80,0.020),
                       (0.80,0.040),(1.00,0.020),(1.00,0.040),(1.20,0.040),(1.00,0.080)]:
        p_, r_, n_ = score(gamma, lam)
        if b is None: b = (p_, r_, n_)
        mark = "  <- PRODUCTION" if (gamma,lam)==(0.65,0.0) else ""
        print(f"{gamma:6.2f} {lam:7.3f} {p_:10.5f} {p_-b[0]:+9.5f} {r_:9.5f} {r_-b[1]:+9.5f} {n_:8.5f} {n_-b[2]:+9.5f}{mark}")

run("HCT116", "../095-clean-harness/de_thinned.npz",
    "../104-context-heldout-credit-adapter/work/hct116/frozen.npz", "real_delta_hct116.npz")
run("HEK293T", "de_hek293t_ref.npz",
    "../104-context-heldout-credit-adapter/work/hek293t/frozen.npz", "real_delta_hek293t.npz")
