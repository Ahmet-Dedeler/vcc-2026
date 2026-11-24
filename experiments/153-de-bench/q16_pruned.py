"""Q14 — fine (gamma, lam) grid with exact pds, ranked by the WORST surface.

A cell that wins on one surface and loses on the other is how this repo has been burned
repeatedly.  Rank on min(HCT, HEK) of the combined overall delta, so only a cell that is good
on both can win.  Combined delta uses the member spans from the 780-entry anchor fit; nmae is
included at its span but flagged, since it is the one member my renorm convention can move.
"""
import numpy as np, csv
from debench import Surface
from exact_pds import pds

canon = [r["gene_name"] for r in csv.DictReader(open("../../data/challenge/gene_names.csv"))]
panel = {r["target_gene"] for r in csv.DictReader(open("../../data/challenge/pert_counts.csv"))}
wcols = np.load("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/"
                "aa534728-8641-48d7-a125-17e75d9c9399/scratchpad/wcols.npy")
wgenes = [canon[i] for i in wcols]
import sys as _s; W = np.load(_s.argv[1] if len(_s.argv)>1 else "../095-clean-harness/out/whitener_union.npz"); ev, V = W["evals"], W["evecs"]
SPAN = {"pds": 0.452, "reach": 0.885, "nmae": 0.60}

def build(de_npz, froz, real_npz):
    z = np.load(froz, allow_pickle=True)
    g104 = {str(g): i for i, g in enumerate(z["genes"])}
    base = z["baseline"][:, [g104[g] for g in wgenes]].astype(np.float32)
    t104 = [str(t) for t in z["targets"]]
    R = np.load(real_npz, allow_pickle=True)
    rg = {str(g): i for i, g in enumerate(R["genes"])}; rt = {str(t): i for i, t in enumerate(R["perts"])}
    pg = np.array([i for i, g in enumerate(wgenes) if g in rg and g not in panel])
    pcol = np.array([rg[wgenes[i]] for i in pg])
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
        return (pds(P[np.ix_(keep_t, pg)].astype(np.float64), Dreal),
                s.reach(key=[np.abs(F[p]) for p in range(len(s.perts))]), s.nmae())
    return score

A = build("../095-clean-harness/de_thinned.npz",
          "../104-context-heldout-credit-adapter/work/hct116/frozen.npz", "real_delta_hct116.npz")
B = build("de_hek293t_ref.npz",
          "../104-context-heldout-credit-adapter/work/hek293t/frozen.npz", "real_delta_hek293t.npz")
bA, bB = A(0.65, 0.0), B(0.65, 0.0)

def comb(x, b):
    return ((x[0]-b[0])/SPAN["pds"] + (x[1]-b[1])/SPAN["reach"] - (x[2]-b[2])/SPAN["nmae"]) / 6.0

rows = []
for gamma in (0.80, 0.90, 1.00, 1.10, 1.20, 1.40):
    for lam in (0.007, 0.010, 0.015, 0.020, 0.030):
        a, b = A(gamma, lam), B(gamma, lam)
        rows.append((min(comb(a, bA), comb(b, bB)), gamma, lam, a, b, comb(a, bA), comb(b, bB)))
rows.sort(key=lambda r: -r[0])
print(f"{'gamma':>6} {'lam':>7} {'HCT Δpds':>9} {'HCT Δreach':>11} {'HEK Δpds':>9} {'HEK Δreach':>11} "
      f"{'ΔovrHCT':>9} {'ΔovrHEK':>9} {'WORST':>9}")
for w, gamma, lam, a, b, ca, cb in rows[:14]:
    flag = "" if (a[0] > bA[0] and b[0] > bB[0] and a[1] > bA[1] and b[1] > bB[1]) else "   (a member regresses)"
    print(f"{gamma:6.2f} {lam:7.3f} {a[0]-bA[0]:+9.5f} {a[1]-bA[1]:+11.5f} {b[0]-bB[0]:+9.5f} "
          f"{b[1]-bB[1]:+11.5f} {ca:+9.5f} {cb:+9.5f} {w:+9.5f}{flag}")
