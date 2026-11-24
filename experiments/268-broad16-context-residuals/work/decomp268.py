"""268 stage 1, steps 2-3: variance decomposition per target (shared / line-specific-reproducible / noise) from
independent cell halves, guide halves, the OR-control empirical null, a synthetic known-answer, and the programme-level
(SVD) analysis of the reproducible residual.  Working lines only (reserved lines never loaded)."""
import json, sys, numpy as np, pandas as pd
from load268 import *
rng = np.random.default_rng(268)
MIN_CPM = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
MIN_CELLS = 6
H, genes = load_halves(WORKING); lines = WORKING; T = usable_targets(H, lines, MIN_CELLS)
print(f"targets with >= {MIN_CELLS} cells on every working line: {len(T)} of {len(COMPLETE)} complete", flush=True)
R, ok, avail = build_responses(H, genes, lines, T + CONTROLS, min_cpm=MIN_CPM)
G = int(ok.sum()); nT = len(T); nL = len(lines)
print(f"genes kept {G} (min mean control CPM {MIN_CPM}), targets {nT} complete + {len(CONTROLS)} OR controls, lines {nL}", flush=True)
assert avail[:nT].all() and R["avail"]["A"][:nT].all() and R["avail"]["B"][:nT].all(), "usable targets must have both halves on every working line"
# for OR-control 'targets' the control profile must exclude the KO itself: rebuild their rows with the other 9 controls
def ctrl_excl(k, l, c):
    tl = H[l]["targets"]; ci = [tl.index(x) for x in CONTROLS if x != c]; return logprof(np.nansum(H[l][k][ci], 0), ok)
AV = R.pop("avail")
for k in R:
    for i, c in enumerate(CONTROLS):
        for j, l in enumerate(lines):
            tl = H[l]["targets"]
            if c in tl and not np.isnan(H[l][k][tl.index(c)][ok]).any() and H[l][k][tl.index(c)][ok].sum() > 0:
                R[k][nT + i, j] = logprof(H[l][k][tl.index(c)], ok) - ctrl_excl(k, l, c)
            else: avail[nT + i, j] = False; AV[k][nT + i, j] = False

def decompose(RA, RB, RF=None):
    """RA, RB: (n, L, G) independent halves. Returns per-target energies (half level) and fractions."""
    mA, mB = RA.mean(1), RB.mean(1); rA, rB = RA - mA[:, None], RB - mB[:, None]
    L = RA.shape[1]
    E_tot = 0.5 * ((RA ** 2).sum((1, 2)) + (RB ** 2).sum((1, 2)))
    X_sh = L * (mA * mB).sum(1); X_ln = (rA * rB).sum((1, 2)); X_tot = X_sh + X_ln
    outd = {"E_tot_half": E_tot, "X_shared": X_sh, "X_line": X_ln, "noise_half": E_tot - X_tot,
            "f_shared": X_sh / E_tot, "f_line": X_ln / E_tot, "f_noise": 1 - X_tot / E_tot}
    # implied at full depth (noise halves): total_full = X_tot + noise_half/2
    tot_full = X_tot + 0.5 * (E_tot - X_tot)
    outd.update({"f_shared_full": X_sh / tot_full, "f_line_full": X_ln / tot_full, "f_noise_full": 0.5 * (E_tot - X_tot) / tot_full})
    if RF is not None:
        E_full = (RF ** 2).sum((1, 2)); outd["noise_ratio_full_over_half"] = (E_full - X_tot) / (E_tot - X_tot)
    # per (t,l) split-half cosine of the residual and of the full response
    nrm = lambda x: np.linalg.norm(x, axis=-1)
    outd["cos_res"] = (rA * rB).sum(-1) / np.maximum(nrm(rA) * nrm(rB), 1e-12)
    outd["cos_resp"] = (RA * RB).sum(-1) / np.maximum(nrm(RA) * nrm(RB), 1e-12)
    # null for the residual cosine: pair A of (t,l) with B of another target, same line
    perm = rng.permutation(RA.shape[0]); rBp = rB[perm]
    outd["cos_res_null"] = (rA * rBp).sum(-1) / np.maximum(nrm(rA) * nrm(rBp), 1e-12)
    # line offset d_l = mean_t r[t,l]: reproducibility and energy share
    dA, dB = rA.mean(0), rB.mean(0)
    outd["offset_cos"] = (dA * dB).sum(-1) / np.maximum(nrm(dA) * nrm(dB), 1e-12)
    # amplitude (parallel to the shared mean of the same half) vs direction (orthogonal) parts of the residual
    def split(r, m):
        u = m / np.maximum(nrm(m), 1e-12)[:, None]; par = (r * u[:, None]).sum(-1, keepdims=True) * u[:, None]; return par, r - par
    pA, qA = split(rA, mA); pB, qB = split(rB, mB)
    outd["X_line_par"] = (pA * pB).sum((1, 2)); outd["X_line_perp"] = (qA * qB).sum((1, 2))
    outd["f_line_par_full"] = outd["X_line_par"] / tot_full; outd["f_line_perp_full"] = outd["X_line_perp"] / tot_full
    outd["cos_res_perp"] = (qA * qB).sum(-1) / np.maximum(nrm(qA) * nrm(qB), 1e-12)
    permq = rng.permutation(RA.shape[0]); qBp = qB[permq]; outd["cos_res_perp_null"] = (qA * qBp).sum(-1) / np.maximum(nrm(qA) * nrm(qBp), 1e-12)
    outd["offset_energy"] = RA.shape[0] * (dA * dB).sum()  # cross-half energy of the offset summed over targets
    return outd, (mA, mB, rA, rB, dA, dB)

def summarize(dd, n, label):
    tot = dd["E_tot_half"][:n].sum(); xs = dd["X_shared"][:n].sum(); xl = dd["X_line"][:n].sum(); nz = dd["noise_half"][:n].sum()
    tot_full = xs + xl + nz / 2
    s = {"label": label, "n_targets": int(n), "pooled_half": {"shared": xs / tot, "line": xl / tot, "noise": nz / tot},
         "pooled_full_depth": {"shared": xs / tot_full, "line": xl / tot_full, "noise": (nz / 2) / tot_full},
         "median_half": {k: float(np.median(dd[f"f_{k}"][:n])) for k in ("shared", "line", "noise")},
         "median_full_depth": {k: float(np.median(dd[f"f_{k}_full"][:n])) for k in ("shared", "line", "noise")},
         "iqr_line_full": [float(np.quantile(dd["f_line_full"][:n], q)) for q in (0.25, 0.75)],
         "offset_share_of_line_energy": float(dd["offset_energy"] / max(xl, 1e-12)) if n > 1 else None,
         "line_par_share": float(dd["X_line_par"][:n].sum() / max(xl, 1e-12)), "line_perp_share": float(dd["X_line_perp"][:n].sum() / max(xl, 1e-12)),
         "pooled_full_depth_line_par": float(dd["X_line_par"][:n].sum() / tot_full), "pooled_full_depth_line_perp": float(dd["X_line_perp"][:n].sum() / tot_full),
         "median_full_depth_line_perp": float(np.median(dd["f_line_perp_full"][:n])),
         "cos_res_perp_median": float(np.median(dd["cos_res_perp"][:n])), "frac_tl_perp_above_null95": float((dd["cos_res_perp"][:n] > np.quantile(dd["cos_res_perp_null"][:n], 0.95)).mean()),
         "offset_cos_median": float(np.median(dd["offset_cos"])),
         "cos_res_median": float(np.median(dd["cos_res"][:n])), "cos_res_null_p95": float(np.quantile(dd["cos_res_null"][:n], 0.95)),
         "frac_tl_above_null95": float((dd["cos_res"][:n] > np.quantile(dd["cos_res_null"][:n], 0.95)).mean()),
         "cos_resp_median": float(np.median(dd["cos_resp"][:n]))}
    if "noise_ratio_full_over_half" in dd: s["noise_ratio_full_over_half_median"] = float(np.median(dd["noise_ratio_full_over_half"][:n]))
    return s

results = {"genes_kept": G, "min_cpm": MIN_CPM, "lines": lines, "n_targets": nT}
# --- known-answer 1: synthetic tensor with known energies
def synth(f_sh, f_ln, f_nz, n=87, L=12, Gs=2000, seed=1):
    r = np.random.default_rng(seed); mu = r.normal(size=(n, 1, Gs)) * np.sqrt(f_sh)
    rho = r.normal(size=(n, L, Gs)); rho -= rho.mean(1, keepdims=True); rho *= np.sqrt(f_ln / (rho ** 2).mean())
    nzA = r.normal(size=(n, L, Gs)) * np.sqrt(f_nz); nzB = r.normal(size=(n, L, Gs)) * np.sqrt(f_nz)
    return (mu + rho + nzA).astype(np.float32), (mu + rho + nzB).astype(np.float32)
for fr in ((0.6, 0.25, 0.15), (0.2, 0.1, 0.7), (0.0, 0.0, 1.0)):
    dA_, dB_ = synth(*fr); dd, _ = decompose(dA_, dB_); s = summarize(dd, 87, f"synthetic {fr}")
    print("KA synthetic", fr, "-> pooled_half", {k: round(v, 3) for k, v in s["pooled_half"].items()}, "frac(t,l) above null95", round(s["frac_tl_above_null95"], 3), flush=True)
    results.setdefault("ka_synthetic", []).append({"truth": fr, "estimate": s["pooled_half"], "frac_above_null95": s["frac_tl_above_null95"]})
# --- real data: cell halves
ddc, (mA, mB, rA, rB, dA, dB) = decompose(R["A"][:nT], R["B"][:nT], R["full"][:nT])
results["cell_halves"] = summarize(ddc, nT, f"cell halves, {nT} usable targets")
print("CELL HALVES", json.dumps(results["cell_halves"], indent=None, default=float), flush=True)
# --- OR-control empirical null (10 'targets', their own controls excluded)
cav = avail[nT:].all(1); Rc = {k: R[k][nT:][cav] for k in R}; nC = int(cav.sum())
ddn, _ = decompose(Rc["A"], Rc["B"], Rc["full"]); results["or_control_null"] = summarize(ddn, nC, f"OR-control null, {nC} controls")
print("OR-CONTROL NULL", json.dumps(results["or_control_null"], default=float), flush=True)
# --- guide halves
gmask = (AV["v1"][:nT] & AV["v2"][:nT]).all(1); nG_ = int(gmask.sum())
ddg, _ = decompose(R["v1"][:nT][gmask], R["v2"][:nT][gmask]); results["guide_halves"] = summarize(ddg, nG_, f"guide halves (v1 vs v2), {nG_} targets with both guides on all {nL} lines")
ddg_full = {k: np.full(nT, np.nan) for k in ("f_line_full", "f_shared_full")}
for k in ddg_full: ddg_full[k][gmask] = ddg[k]
print("GUIDE HALVES", json.dumps(results["guide_halves"], default=float), flush=True)
# --- per-target table
per = pd.DataFrame({"target": T, "shared_full": ddc["f_shared_full"], "line_full": ddc["f_line_full"], "noise_full": ddc["f_noise_full"],
                    "shared_half": ddc["f_shared"], "line_half": ddc["f_line"], "noise_half": ddc["f_noise"],
                    "line_full_guide": ddg_full["f_line_full"], "shared_full_guide": ddg_full["f_shared_full"],
                    "line_perp_full": ddc["f_line_perp_full"], "line_par_full": ddc["f_line_par_full"], "n_lines_perp_repro": (ddc["cos_res_perp"] > np.quantile(ddc["cos_res_perp_null"], 0.95)).sum(1),
                    "n_cells_min": [min(H[l]["ncell"][H[l]["targets"].index(t), 0] for l in lines) for t in T],
                    "n_lines_res_repro": (ddc["cos_res"] > np.quantile(ddc["cos_res_null"], 0.95)).sum(1),
                    "cos_resp_median": np.median(ddc["cos_resp"], 1), "E_tot_half": ddc["E_tot_half"], "X_shared": ddc["X_shared"], "X_line": ddc["X_line"]})
per.to_csv(W / "per_target_decomp.csv", index=False)
# --- batch structure of the residual: cross-half line x line similarity
batch = [META.loc[l, "batch"] for l in lines]
S = np.zeros((nL, nL))
for i in range(nL):
    for j in range(nL):
        S[i, j] = (rA[:, i] * rB[:, j]).sum() / np.sqrt(max((rA[:, i] * rB[:, i]).sum(), 1e-9) * max((rA[:, j] * rB[:, j]).sum(), 1e-9))
wb = [S[i, j] for i in range(nL) for j in range(nL) if i != j and batch[i] == batch[j]]; bb = [S[i, j] for i in range(nL) for j in range(nL) if i != j and batch[i] != batch[j]]
obs = np.mean(wb) - np.mean(bb); null = []
for _ in range(2000):
    pb_ = list(rng.permutation(batch)); w = [S[i, j] for i in range(nL) for j in range(nL) if i != j and pb_[i] == pb_[j]]; b = [S[i, j] for i in range(nL) for j in range(nL) if i != j and pb_[i] != pb_[j]]
    null.append(np.mean(w) - np.mean(b))
results["batch"] = {"within_batch_mean_sim": float(np.mean(wb)), "between_batch_mean_sim": float(np.mean(bb)), "n_within_pairs": len(wb) // 2, "perm_p": float((np.array(null) >= obs).mean()), "line_sim_matrix": S.round(3).tolist(), "batches": batch}
print("BATCH", {k: v for k, v in results["batch"].items() if k != "line_sim_matrix"}, flush=True)
# --- programme-level: SVD of the residual (fit on A, cross-half energy on B)
def svd_repro(A, B, label, kmax=60):
    U, s, Vt = np.linalg.svd(A, full_matrices=False); kmax = min(kmax, len(s))
    a = A @ Vt.T; b = B @ Vt.T; c = (a * b).sum(0); tot = (A * B).sum()
    cum = np.cumsum(c) / tot
    # null: permute the rows of B
    nulls = np.array([((a * b[rng.permutation(len(b))]).sum(0))[:kmax] for _ in range(200)])
    p99 = np.quantile(nulls, 0.99, axis=0); sig = c[:kmax] > p99
    corr = np.array([np.corrcoef(a[:, k], b[:, k])[0, 1] for k in range(kmax)])
    out_ = {"label": label, "n_rows": int(A.shape[0]), "sv_top10": s[:10].round(3).tolist(), "cross_energy_share_top": c[:kmax] / tot, "cum_share": cum[:kmax],
            "k_for_50_80_90": [int(np.searchsorted(cum, q) + 1) if cum.max() >= q else None for q in (0.5, 0.8, 0.9)],
            "n_sig_components_p99": int(sig.sum()), "first_nonsig": int(np.argmin(sig)) if not sig.all() else kmax, "corr_ab": corr, "Vt": Vt[:10], "scores_a": a[:, :10], "cross_energy": c[:kmax], "null_p99": p99}
    print(label, "cross-half energy share of top 1..10:", np.round(c[:10] / tot, 3).tolist(), "cum:", np.round(cum[[0, 2, 4, 9, 19, 39]], 3).tolist(), "k50/80/90", out_["k_for_50_80_90"], "sig comps", out_["n_sig_components_p99"], "first non-sig", out_["first_nonsig"], flush=True)
    return out_
resA = rA.reshape(nT * nL, G); resB = rB.reshape(nT * nL, G)
sv_full = svd_repro(resA, resB, "residual (with line offset)")
iA = (rA - dA[None]).reshape(nT * nL, G); iB = (rB - dB[None]).reshape(nT * nL, G)
sv_int = svd_repro(iA, iB, "interaction residual (line offset removed)")
sv_sh = svd_repro(mA, mB, "shared response (87 x G)")
sv_noise = svd_repro(*(x.reshape(-1, G) for x in (Rc["A"] - Rc["A"].mean(1, keepdims=True), Rc["B"] - Rc["B"].mean(1, keepdims=True))), "OR-control residual (noise null)")
results["svd"] = {k: {kk: (vv.tolist() if isinstance(vv, np.ndarray) and kk not in ("Vt", "scores_a") else vv) for kk, vv in v.items() if kk not in ("Vt", "scores_a")} for k, v in (("residual", sv_full), ("interaction", sv_int), ("shared", sv_sh), ("noise_null", sv_noise))}
# component interpretation: top genes, line/target eta^2, hallmark z
gs = pd.read_csv("~/Code/virtual-cell-challenge-2026/data/external/broad16line/data/all_genesets.csv"); gs = gs[gs.collection == "Hallmark"]
kept = genes[ok]; gidx = {g: i for i, g in enumerate(kept)}; sets = {t: np.array([gidx[g] for g in grp.gene if g in gidx]) for t, grp in gs.groupby("term")}
comp = []
for k in range(6):
    v = sv_full["Vt"][k]; sc_ = sv_full["scores_a"][:, k].reshape(nT, nL)
    eta_line = sc_.mean(0).var() / max(sc_.var(), 1e-12)
    ce_t = (sc_ * (resB @ v).reshape(nT, nL)).sum(1); eta_t = float(np.sort(np.abs(ce_t))[::-1][:5].sum() / max(np.abs(ce_t).sum(), 1e-12))
    top_pos = kept[np.argsort(-v)[:12]].tolist(); top_neg = kept[np.argsort(v)[:12]].tolist()
    zs = {}
    for t, ix in sets.items():
        if len(ix) < 15: continue
        m = v[ix].mean(); nullm = np.array([v[rng.choice(G, len(ix), replace=False)].mean() for _ in range(300)]); zs[t] = (m - nullm.mean()) / max(nullm.std(), 1e-12)
    top_sets = sorted(zs.items(), key=lambda x: -abs(x[1]))[:4]
    lines_top = [lines[i] for i in np.argsort(-np.abs(sc_.mean(0)))[:3]]; targets_top = [T[i] for i in np.argsort(-np.abs(sc_).mean(1))[:5]]
    comp.append({"k": k + 1, "cross_energy_share": float(sv_full["cross_energy"][k] / (resA * resB).sum()), "corr_ab": float(sv_full["corr_ab"][k]), "eta2_line": float(eta_line), "top5_target_share": float(eta_t),
                 "top_pos": top_pos, "top_neg": top_neg, "hallmark_z": [(t.replace("HALLMARK_", ""), round(z, 1)) for t, z in top_sets], "lines_top": lines_top, "targets_top": targets_top})
    print(f"C{k+1}: share {comp[-1]['cross_energy_share']:.3f} corr {comp[-1]['corr_ab']:.2f} eta2 line {eta_line:.2f} top5-target share {eta_t:.2f} lines {lines_top} targets {targets_top}\n   + {top_pos[:8]}\n   - {top_neg[:8]}\n   hallmark {comp[-1]['hallmark_z']}", flush=True)
results["components"] = comp
json.dump(results, open(W / f"decomp268_cpm{int(MIN_CPM)}.json", "w"), indent=1, default=float)
per.to_csv(W / f"per_target_decomp_cpm{int(MIN_CPM)}.csv", index=False)
np.savez_compressed(W / f"decomp268_arrays_cpm{int(MIN_CPM)}.npz", cos_res=ddc["cos_res"], cos_res_null=ddc["cos_res_null"], cos_resp=ddc["cos_resp"], f_line_full=ddc["f_line_full"], f_shared_full=ddc["f_shared_full"], Vt=sv_full["Vt"], scores_a=sv_full["scores_a"], kept_genes=kept)
print("done", flush=True)
