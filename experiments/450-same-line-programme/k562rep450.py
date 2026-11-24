"""450e: does the same-line effect replicate in a CANCER line? Truth = VIPerturb-seq K562 (Satija lab 2026, a different lab
and technology, median 48 cells per perturbation). Predictors: Replogle K562-GWPS (same line), Orion HCT116 and Orion
HEK293T (other lines, the production donors). Per target, on genes all four banks measure:
  cosine (raw, and centred = each bank's mean over the evaluated targets removed),
  sign agreement on the truth's top-100 |delta| genes.
Targets split by K562-essential membership (the challenge panel is essential-free, so the non-essential stratum is the
relevant one). Deterministic sample of up to 1,500 targets per stratum (seed 20260925).
"""
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data/donors"


def bank(name):
    d = D / name
    return (np.load(d / "delta.npy", mmap_mode="r"), [l.strip() for l in open(d / "perts.txt")],
            [l.strip() for l in open(d / "genes.txt")])


names = {"truth_viperturb": "viperturb_k562", "k562_gwps": "k562_gwps", "hct116": "orion_hct116", "hek293t": "orion_hek293t"}
B = {k: bank(v) for k, v in names.items()}
ess = set(l.strip() for l in open(D / "k562_essential/perts.txt"))
common_t = set(B["truth_viperturb"][1])
for k in ("k562_gwps", "hct116", "hek293t"):
    common_t &= set(B[k][1])
common_g = set(B["truth_viperturb"][2])
for k in ("k562_gwps", "hct116", "hek293t"):
    common_g &= set(B[k][2])
common_g = sorted(common_g)                     # each target's own gene is zeroed per row below (cis effect out)
rng = np.random.default_rng(20260925)
strata = {"non-essential": sorted(t for t in common_t if t not in ess), "essential": sorted(t for t in common_t if t in ess)}
print(f"common targets {len(common_t)} (non-essential {len(strata['non-essential'])}, essential {len(strata['essential'])}); "
      f"common genes {len(common_g)}")
for sname, ts in strata.items():
    ts = list(rng.choice(ts, size=min(1500, len(ts)), replace=False))
    M = {}
    for k, (delta, perts, genes) in B.items():
        pi = {p: i for i, p in enumerate(perts)}; gi = {g: i for i, g in enumerate(genes)}
        rows = np.array([pi[t] for t in ts]); cols = np.array([gi[g] for g in common_g])
        order = np.argsort(rows)
        X = np.empty((len(ts), len(cols)), dtype=np.float64)
        X[order] = np.asarray(delta[rows[order]], dtype=np.float64)[:, cols]
        own = {g: j for j, g in enumerate(common_g)}
        for r, t in enumerate(ts):
            if t in own:
                X[r, own[t]] = 0.0
        M[k] = X
    T = M["truth_viperturb"]
    strong = np.linalg.norm(T - T.mean(0), axis=1) >= np.quantile(np.linalg.norm(T - T.mean(0), axis=1), 0.75)
    Tc = T - T.mean(0)
    top = np.argsort(-np.abs(T), axis=1)[:, :100]

    def cos(A, Bm):
        return np.sum(A * Bm, 1) / np.maximum(np.linalg.norm(A, axis=1) * np.linalg.norm(Bm, axis=1), 1e-12)

    print(f"\n== {sname}: {len(ts)} targets")
    print(f"{'predictor':12s} {'cos raw':>8s} {'cos centred':>12s} {'sign@top100':>12s} {'centred cos > hct116':>21s}")
    ref = None
    for k in ("k562_gwps", "hct116", "hek293t"):
        P = M[k]; Pc = P - P.mean(0)
        cr, cc = cos(P, T), cos(Pc, Tc)
        sg = np.mean(np.sign(np.take_along_axis(P, top, 1)) == np.sign(np.take_along_axis(T, top, 1)), axis=1)
        if k == "hct116":
            ref = cc
        print(f"{k:12s} {cr.mean():8.4f} {cc.mean():12.4f} {sg.mean():12.4f}", end="")
        print("" if k != "k562_gwps" else f"   (vs hct116 per target: computed below)")
        if k == "k562_gwps":
            kc = cc
    for k in ("k562_gwps", "hct116"):
        P = M[k]; Pc = P - P.mean(0)
        print(f"  strongest-quartile truth rows ({int(strong.sum())}): {k} centred cos {cos(Pc, Tc)[strong].mean():.4f}")
    print(f"k562_gwps centred cosine beats hct116 on {int((kc > ref).sum())} / {len(ref)} targets "
          f"(mean diff {float((kc - ref).mean()):+.4f})")
