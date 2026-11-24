#!/usr/bin/env python3
"""402 -- per-context per-target prior artifact: contexts B/C = 399's rows (six-bank family, K=25, shrink 400, family-marginal
default); context A = Jurkat-only votes (the 25 nearest GWPS perturbations' Nadig-Jurkat rows, shrink 50 toward the Jurkat
marginal; default = Jurkat marginal sign, CLEAN4 where no vote). Keys: genes, psign (399 default), psign_clean4, targets,
psign_by_target (399 rows), has_row, psign_A (A default), psign_by_target_A, has_row_A, P_by_target_A, source shas."""
import numpy as np, json, hashlib, time
from pathlib import Path
REPO = Path("~/Code/virtual-cell-challenge-2026"); OUT = REPO / "experiments/402-ctxA-jurkat-knn/work"; OUT.mkdir(parents=True, exist_ok=True)
t0 = time.time()
base = np.load(REPO / "experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.npz", allow_pickle=False)
off = [str(g) for g in base["genes"]]; panel = [str(t) for t in base["targets"]]; nG = len(off)
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
gd = REPO / "data/donors/k562_gwps"; GW = np.asarray(np.load(gd / "delta.npy", mmap_mode="r"), dtype=np.float32)
GWP = np.array([l.strip() for l in open(gd / "perts.txt")]); gwp = {p: i for i, p in enumerate(GWP)}
GWc = GW - GW.mean(1, keepdims=True); GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)
jd = REPO / "data/donors/nadig_jurkat"; JD = np.asarray(np.load(jd / "delta.npy", mmap_mode="r"), dtype=np.float32)
JP = np.array([l.strip() for l in open(jd / "perts.txt")]); JG = np.array([l.strip() for l in open(jd / "genes.txt")])
jp = {p: i for i, p in enumerate(JP)}; jg = {g: i for i, g in enumerate(JG)}; jmap = np.array([jg.get(g, -1) for g in off]); cvj = jmap >= 0
def votes_j(rows_gw):
    ridx = [jp[GWP[j]] for j in rows_gw if GWP[j] in jp]
    dn = np.zeros(nG); nn = np.zeros(nG)
    if ridx:
        sub = JD[ridx][:, np.where(cvj, jmap, 0)]; nz = (sub != 0) & cvj[None, :]; dn += ((sub < 0) & nz).sum(0); nn += nz.sum(0)
    return dn, nn
dn_m, nn_m = votes_j(np.arange(len(GWP))); PJm = np.where(nn_m > 0, dn_m / np.maximum(nn_m, 1), np.nan)
psign4 = np.asarray(base["psign_clean4"], dtype=np.float64)
psA = np.where(np.isfinite(PJm), np.where(PJm > 0.5, -1.0, 1.0), psign4)
K, M = 25, 50.0
has = np.array([t in gwp for t in panel]); PSA = np.tile(psA, (300, 1)); PTA = np.full((300, nG), np.nan, dtype=np.float32)
for i, t in enumerate(panel):
    if not has[i]: continue
    s = GWn @ GWn[gwp[t]]; s[gwp[t]] = -np.inf; nb = np.argsort(-s)[:K]; dn, nn = votes_j(nb); P = np.where(nn > 0, dn / np.maximum(nn, 1), np.nan)
    Ps = np.where(np.isfinite(P) & np.isfinite(PJm), (nn * np.nan_to_num(P) + M * np.nan_to_num(PJm)) / (nn + M), np.where(np.isfinite(P), P, PJm))
    PTA[i] = Ps; PSA[i] = np.where(np.isfinite(Ps), np.where(Ps > 0.5, -1.0, 1.0), psA)
out = OUT / "knn25o400_A_jurkat50_prior_challenge.npz"
np.savez(out, genes=base["genes"], psign=base["psign"], psign_clean4=base["psign_clean4"], targets=base["targets"], psign_by_target=base["psign_by_target"],
         has_row=base["has_row"], P_by_target=base["P_by_target"], P_marginal=base["P_marginal"], K=base["K"], shrink_m=base["shrink_m"], banks=base["banks"],
         psign_A=psA, psign_by_target_A=PSA, has_row_A=has, P_by_target_A=PTA, P_marginal_A=PJm.astype(np.float32), K_A=np.array(K), shrink_m_A=np.array(M), banks_A=np.array(["nadig_jurkat"]))
blk = np.load(REPO / "experiments/395-protected-signs-arm/work/jurkat_tilt_A_scratch.npz")["block_A"]
diff_vs_399 = (PSA[:, blk] != np.asarray(base["psign_by_target"])[:, blk]).mean(1)
meta = {"artifact_sha256": sha(out), "base_399_sha256": sha(REPO / "experiments/399-knn-tilt-arm/work/knn25o400_prior_challenge.npz"),
        "jurkat_delta_sha256": sha(jd / "delta.npy"), "gwps_delta_sha256": sha(gd / "delta.npy"), "K_A": K, "shrink_m_A": M,
        "targets_with_row_A": int(has.sum()), "coverage_A_marginal": float(np.isfinite(PJm).mean()),
        "A_block_signs_differing_from_399_per_target_median": float(np.median(diff_vs_399[has])), "elapsed_s": round(time.time() - t0, 1)}
json.dump(meta, open(str(out).replace(".npz", ".json"), "w"), indent=1); print(json.dumps(meta, indent=1))
