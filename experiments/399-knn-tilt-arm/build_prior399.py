#!/usr/bin/env python3
"""399 -- per-target neighbour-conditioned block sign prior for the challenge panel (397/398 construction).

For each panel target present in K562-GWPS: neighbours = top-K GWPS perturbations by Pearson correlation of GWPS response
rows (self excluded); P_down(g) = fraction of down votes over the neighbours' rows in the CLEAN4 bank family (GWPS,
K562-essential, Nadig Jurkat, Nadig HepG2; non-zero entries only), shrunk toward the same family's marginal P_down with
prior weight m; sign = down if P > 0.5 else up; genes with no vote -> the pinned CLEAN4 sign (which itself falls back to
down); targets absent from GWPS -> the CLEAN4 sign vector. Writes an artifact on the official 18,533-gene axis:
genes, psign (the default vector: CLEAN4 for the clean4 family, else the family marginal sign), psign_clean4, targets (panel order), psign_by_target (300 x 18533),
has_row (300,), P_by_target (float32, NaN where no vote), plus a JSON sidecar with input hashes and diagnostics."""
import numpy as np, json, hashlib, sys, time
from pathlib import Path
K = int(sys.argv[1]) if len(sys.argv) > 1 else 25; M = float(sys.argv[2]) if len(sys.argv) > 2 else 400.0
FAMILY = sys.argv[3] if len(sys.argv) > 3 else "orion"   # "clean4" = 4 banks (397/398 first pass); "orion" = + raw Orion HCT116/HEK293T (398c/398d)
REPO = Path("~/Code/virtual-cell-challenge-2026"); OUT = REPO / "experiments/399-knn-tilt-arm/work"; OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time()
off = [l.strip() for l in open(REPO / "data/challenge/gene_names.csv")]
if off[0] != "TSPAN6": off = off[1:]
panel = [l.strip() for l in open(REPO / "data/challenge/pert_counts.csv")][1:]
assert len(off) == 18533 and len(panel) == 300, (len(off), len(panel))
z = np.load(REPO / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=False)
assert list(map(str, z["genes"])) == off; psign4 = np.asarray(z["psign"], dtype=np.float64)
BANKS = {"k562_gwps": REPO / "data/donors/k562_gwps", "k562_essential": REPO / "data/donors/k562_essential",
         "nadig_jurkat": REPO / "data/donors/nadig_jurkat", "nadig_hepg2": REPO / "experiments/166-donor-acquisition/donors/nadig_hepg2"}
if FAMILY == "orion":
    BANKS.update({"orion_hct116": REPO / "data/donors/orion_hct116", "orion_hek293t": REPO / "data/donors/orion_hek293t"})
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
B = {}
for b, d in BANKS.items():
    D = np.asarray(np.load(d / "delta.npy", mmap_mode="r"), dtype=np.float32)
    P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
    gi = {g: i for i, g in enumerate(G)}; B[b] = (D, {p: i for i, p in enumerate(P)}, np.array([gi.get(g, -1) for g in off]), P)
    print(b, D.shape, "delta sha", sha(d / "delta.npy")[:16], flush=True)
GW, gwp, gwmap, GWP = B["k562_gwps"]
GWc = GW - GW.mean(1, keepdims=True); GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)
nG = len(off)
def votes(rows_gw):
    """down / total votes on the official axis over the given GWPS perturbation indices, across the four banks."""
    dn = np.zeros(nG); nn = np.zeros(nG); names = [GWP[j] for j in rows_gw]
    for b, (D, pi, gm, _) in B.items():
        ridx = [pi[p] for p in names if p in pi]
        if not ridx: continue
        cv = gm >= 0; sub = D[ridx][:, np.where(cv, gm, 0)]; nz = (sub != 0) & cv[None, :]
        dn += ((sub < 0) & nz).sum(0); nn += nz.sum(0)
    return dn, nn
dn_m, nn_m = votes(np.arange(len(GWP))); P_marg = np.where(nn_m > 0, dn_m / np.maximum(nn_m, 1), np.nan)
psign_marg = np.where(np.isfinite(P_marg), np.where(P_marg > 0.5, -1.0, 1.0), psign4)
print(f"marginal (4 banks): coverage {np.isfinite(P_marg).mean():.3f}; agreement with pinned CLEAN4 on covered genes "
      f"{(psign_marg == psign4)[np.isfinite(P_marg)].mean():.3f}", flush=True)
has = np.array([t in gwp for t in panel]); print(f"panel targets in GWPS: {has.sum()}/300", flush=True)
PS = np.tile(psign_marg, (300, 1)); PT = np.full((300, nG), np.nan, dtype=np.float32); top1 = np.full(300, np.nan); nflip_vs4 = np.zeros(300, int)
for i, t in enumerate(panel):
    if not has[i]: continue
    s = GWn @ GWn[gwp[t]]; s[gwp[t]] = -np.inf; nb = np.argsort(-s)[:K]; top1[i] = s[nb[0]]
    dn, nn = votes(nb); P = np.where(nn > 0, dn / np.maximum(nn, 1), np.nan)
    Ps = np.where(np.isfinite(P) & np.isfinite(P_marg), (nn * np.nan_to_num(P) + M * np.nan_to_num(P_marg)) / (nn + M),
                  np.where(np.isfinite(P), P, P_marg))
    PT[i] = Ps; PS[i] = np.where(np.isfinite(Ps), np.where(Ps > 0.5, -1.0, 1.0), psign_marg)
    nflip_vs4[i] = int((PS[i] != psign4).sum())
out = OUT / (f"knn{K}s{int(M)}_prior_challenge.npz" if FAMILY == "clean4" else f"knn{K}o{int(M)}_prior_challenge.npz")
np.savez(out, genes=np.array(off), psign=(psign4 if FAMILY == "clean4" else psign_marg), psign_clean4=psign4, targets=np.array(panel), psign_by_target=PS, has_row=has, P_by_target=PT,
         P_marginal=P_marg.astype(np.float32), K=np.array(K), shrink_m=np.array(M), banks=np.array(list(BANKS)))
meta = {"artifact_sha256": sha(out), "K": K, "shrink_m": M, "family": FAMILY, "default_vector": ("CLEAN4 psign" if FAMILY == "clean4" else "family marginal sign, CLEAN4 where no vote"), "banks": {b: sha(d / "delta.npy") for b, d in BANKS.items()},
        "clean4_source_sha256": sha(REPO / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz"),
        "targets_with_row": int(has.sum()), "top1_neighbour_corr_median": float(np.nanmedian(top1)),
        "genes_changed_vs_clean4_per_target_median": float(np.median(nflip_vs4[has])),
        "coverage_marginal": float(np.isfinite(P_marg).mean()), "rule": "see docstring", "elapsed_s": round(time.time() - t0, 1)}
json.dump(meta, open(str(out).replace(".npz", ".json"), "w"), indent=1); print(json.dumps(meta, indent=1), flush=True)
