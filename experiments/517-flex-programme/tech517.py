#!/usr/bin/env python3
"""517 §3 — is part of a Flex "programme" a library-size (composition) effect, and does it predict DE signs?

Mechanism under test. The scorer's DE reads per-cell CPM (arithmetic means, Wilcoxon ranks). If a gene's per-cell CPM
scales with the cell's library as L^eps_g (eps_g = 0: proportional; > 0: larger share in larger cells), then a
perturbation that shifts its cells' libraries by delta_t = E[ln L | t] - E[ln L | ctrl] moves that gene's mean CPM by
~ eps_g * delta_t with no biology of its own for g. eps_g is measurable in any context's CONTROL cells, so if perturbations
mostly shift libraries one way, sign(eps_g * typical delta) is a zero-shot, recipient-side programme component.

H1 (Arc's own Flex data):
  1. delta_t for all 300 public perturbations (raw files, all cells) and for h1-test's 85 (its own 400 cells vs its controls);
     its distribution, and its relation to response strength (h1-test n_real; the bank row's norm);
  2. eps_g in the 38,176 public H1 controls: 10 library-decile bins, slope of ln(arithmetic-mean CPM) on the bin's mean
     ln L (the scorer's statistic), with split-half (random control halves) reliability;
  3. on h1-test's gate genes: sign accuracy of eps_g * delta_t (per-target oracle delta) and of eps_g * median delta
     (target-independent), overall and by |delta_t| tercile; its correlation with Arc's programme and CLEAN4; whether it
     adds to them (z-sum accuracy); the share of each target's truth lfc variance it explains.
A/B/C (the challenge's own controls, 46 NTC guides x 400 cells):
  4. eps_g per context; agreement across A/B/C/H1 (a chemistry-technical axis would be shared);
  5. NTC guides: between-guide variance of mean ln L vs the within-guide sampling expectation, and whether guide pseudobulk
     deviations line up with eps (a library-driven axis among cells that carry no perturbation at all).
Known answers: h1-test per-target mean CPM vs its controls reproduces the reference lfc (as halves517); the raw control
count equals 38,176; per-perturbation cell counts equal the official pert_counts CSVs.
Memory: stream517 (preallocated buffers).  usage: uv run python experiments/517-flex-programme/tech517.py
"""
from __future__ import annotations

import csv
import glob
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
from stream517 import CSRStream  # noqa: E402

FILES = [ROOT / "data/external/vcc2025/train/adata_Training.h5ad", ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
         ROOT / "data/external/vcc2025/test/adata_Test.h5ad"]
COUNTS = [ROOT / "data/external/vcc2025/train/pert_counts_Training.csv",
          ROOT / "data/external/vcc2025/validation/pert_counts_Validation.csv",
          ROOT / "data/external/vcc2025/test/pert_counts_Test.csv"]
H1T = ROOT / "data/harness/h1-test.h5ad"
NB = 10
t0 = time.time()
log = lambda *x: print("[517t]", *x, f"({time.time() - t0:.0f}s)", flush=True)   # noqa: E731


def cat(n):
    return n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)


def ids_of(f):
    ix = f["obs/_index"]
    return (ix[:] if not isinstance(ix, h5py.Group) else ix["values"][:]).astype(str)


def eps_from_controls(it, ctrl_rows, lib, cols, rng, n_bins=NB):
    """eps_g: OLS slope of ln(arithmetic-mean CPM in a library bin) on the bin's mean ln L, over n_bins library-quantile
    bins of the given control rows; also the two random-half estimates (reliability). Returns dict of arrays on `cols`."""
    lnL = np.log(lib[ctrl_rows])
    edges = np.quantile(lnL, np.linspace(0, 1, n_bins + 1))
    b = np.clip(np.searchsorted(edges, lnL, side="right") - 1, 0, n_bins - 1)
    half = rng.integers(0, 2, ctrl_rows.size)
    labels = np.full(it.n, -1, np.int64)
    labels[ctrl_rows] = b + n_bins * half                       # 2 * n_bins groups: (half, bin)
    S = it.group_sums(labels, 2 * n_bins, cols=cols, cpm=True, lib=lib)
    n = np.bincount(labels[labels >= 0], minlength=2 * n_bins).astype(float)
    xb = np.array([lnL[(b == k) & (half == h)].mean() for h in (0, 1) for k in range(n_bins)])
    A = S / n[:, None]                                          # arithmetic-mean CPM per (half, bin)

    def slope(Ab, x):
        y = np.log(np.maximum(Ab, 1e-3))
        xc = x - x.mean()
        return (xc[:, None] * (y - y.mean(0))).sum(0) / (xc * xc).sum()
    Aall = (S[:n_bins] + S[n_bins:]) / (n[:n_bins] + n[n_bins:])[:, None]
    xall = np.array([lnL[b == k].mean() for k in range(n_bins)])
    return dict(eps=slope(Aall, xall), eps_h0=slope(A[:n_bins], xb[:n_bins]), eps_h1=slope(A[n_bins:], xb[n_bins:]),
                mean_cpm=(S.sum(0) / n.sum()), bin_lnL=xall, n_ctrl=int(ctrl_rows.size))


def pearson(a, b, m=None):
    m = np.isfinite(a) & np.isfinite(b) if m is None else (m & np.isfinite(a) & np.isfinite(b))
    return float(np.corrcoef(a[m], b[m])[0, 1])


def main():
    out: dict = {}
    rng = np.random.default_rng(517)
    z = np.load(glob.glob(str(ROOT / "experiments/323-*/work/reference.npz"))[0], allow_pickle=False)
    genes, keep, targets = z["genes"].astype(str), np.asarray(z["keep"], bool), list(z["targets"].astype(str))
    q_ref, lfc_ref, bctrl = np.asarray(z["q"], float), np.asarray(z["lfc"], float), np.asarray(z["b_ctrl"], float)
    z.close()
    kg = genes[keep]
    cols_keep = np.full(genes.size, -1, np.int64); cols_keep[np.flatnonzero(keep)] = np.arange(int(keep.sum()))
    # ------------------------------------------------ 1. libraries of every public H1 cell; delta_t per perturbation
    per = {}
    ctrl_lnL = None
    eps = None
    seen: set[str] = set()
    for path in FILES:
        it = CSRStream(path)
        tv, ids = cat(it.f["obs/target_gene"]), ids_of(it.f)
        lib = it.libraries()
        is_c = tv == "non-targeting"
        new_c = is_c & ~np.isin(ids, list(seen)) if seen else is_c.copy()
        if eps is None:
            # every control is in the first (training) file: eps from its controls, once
            assert int(new_c.sum()) == 38176, int(new_c.sum())
            eps = eps_from_controls(it, np.flatnonzero(new_c), lib, cols_keep, rng)
            ctrl_lnL = np.log(lib[new_c])
            log(f"eps from {eps['n_ctrl']} H1 controls; split-half r(eps) {pearson(eps['eps_h0'], eps['eps_h1']):.3f}; "
                f"bin mean ln L {np.round(eps['bin_lnL'], 2).tolist()}")
        seen |= set(ids[is_c])
        for t in np.unique(tv[~is_c]):
            m = tv == t
            per[t] = dict(n=int(m.sum()), mean_lnL=float(np.log(lib[m]).mean()), median_L=float(np.median(lib[m])))
        it.close()
        log(f"{path.name}: {tv.size} cells")
    want = {}
    for c in COUNTS:
        for r_ in csv.DictReader(open(c)):
            want[r_["target_gene"]] = int(r_["n_cells"])
    assert all(per[t]["n"] == n for t, n in want.items()) and set(want) == set(per), "KA cell counts"
    mu_c = float(ctrl_lnL.mean())
    for t in per:
        per[t]["delta"] = per[t]["mean_lnL"] - mu_c
    d_all = np.array([per[t]["delta"] for t in per])
    out["h1_delta_all300"] = dict(mean=float(d_all.mean()), median=float(np.median(d_all)), sd=float(d_all.std()),
                                  frac_neg=float((d_all < 0).mean()), q10=float(np.quantile(d_all, .1)), q90=float(np.quantile(d_all, .9)),
                                  ctrl_median_L=float(np.exp(np.median(ctrl_lnL))), ctrl_sd_lnL=float(ctrl_lnL.std()),
                                  se_of_delta_at_1071_cells=float(ctrl_lnL.std() / np.sqrt(1071)))
    log(f"KA: per-perturbation cell counts == pert_counts CSVs. delta over 300 perturbations: mean {d_all.mean():+.4f}, "
        f"median {np.median(d_all):+.4f}, sd {d_all.std():.4f}, frac < 0 {(d_all < 0).mean():.3f} (control sd ln L "
        f"{ctrl_lnL.std():.3f}; se at 1,071 cells {ctrl_lnL.std() / np.sqrt(1071):.4f})")
    # strength of each perturbation's response: the bank row's norm over expressed genes
    d = ROOT / "data/donors/h1_2025_full"
    bdelta = np.load(d / "delta.npy", mmap_mode="r"); bp = [l.strip() for l in open(d / "perts.txt")]
    norms = {p: float(np.linalg.norm(np.asarray(bdelta[i], dtype=np.float64)[keep])) for i, p in enumerate(bp)}
    nm = np.array([norms[t] for t in per])
    from scipy.stats import spearmanr
    out["h1_delta_vs_strength_spearman"] = float(spearmanr(d_all, nm)[0])
    log(f"Spearman(delta_t, bank row norm) over 300: {out['h1_delta_vs_strength_spearman']:+.3f}")
    # ------------------------------------------------ 2. h1-test: its own delta_t; truth; accuracies
    it = CSRStream(H1T)
    tv = cat(it.f["obs/target_gene"])
    lib = it.libraries()
    is_c = tv == "non-targeting"
    mu = float(np.log(lib[is_c]).mean())
    d_t = np.array([np.log(lib[tv == t]).mean() - mu for t in targets])
    tid = {t: i for i, t in enumerate(targets)}
    lab = np.array([tid.get(t, len(targets) if t == "non-targeting" else -1) for t in tv])
    S = it.group_sums(lab, len(targets) + 1, cols=cols_keep, cpm=True, lib=lib)
    it.close()
    nn = np.bincount(lab[lab >= 0], minlength=len(targets) + 1)
    M = S / nn[:, None]
    L = np.log2((M[:-1] + 1e-9) / (M[-1] + 1e-9))
    fin = np.isfinite(lfc_ref) & np.isfinite(L)
    ka = float(np.abs(L[fin] - lfc_ref[fin]).max())
    log(f"KA h1-test: truth recomputed vs reference lfc max |d| {ka:.2e}; h1-test delta_t: mean {d_t.mean():+.4f}, "
        f"frac < 0 {(d_t < 0).mean():.3f}, sd {d_t.std():.4f}")
    assert ka < 1e-6
    nreal = (q_ref < 0.05).sum(1)
    out["h1test_delta"] = dict(mean=float(d_t.mean()), frac_neg=float((d_t < 0).mean()), sd=float(d_t.std()),
                               spearman_vs_nreal=float(spearmanr(d_t, nreal)[0]),
                               raw_vs_harness_delta_r=float(np.corrcoef(d_t, [per[t]["delta"] for t in targets])[0, 1]))
    log(f"h1-test: Spearman(delta_t, n_real) {out['h1test_delta']['spearman_vs_nreal']:+.3f}; r(harness delta, all-cell delta) "
        f"{out['h1test_delta']['raw_vs_harness_delta_r']:.3f}")
    e = eps["eps"]
    own = np.array([[g == t for g in kg] for t in targets])
    gate = (q_ref < 0.05) & np.isfinite(lfc_ref) & (lfc_ref != 0) & ~own
    e_ = np.expm1(bctrl[keep]); cpm = 1e6 * e_ / e_.sum()
    order = np.argsort(-cpm)
    block = np.zeros(kg.size, bool); block[[i for i in order if kg[i] not in set(targets)][:2200]] = True
    src = np.load(HERE / "work/sources517.npz", allow_pickle=False)
    assert list(src["genes_keep"]) == list(kg)
    P_arc, P_c4 = src["P_arc"], src["P_c4"]

    def acc_pred(pred, region, tsel=None):
        hit = n = 0
        for t in range(len(targets)) if tsel is None else tsel:
            m = gate[t] & region & np.isfinite(pred[t]) & (pred[t] != 0)
            hit += int((np.sign(pred[t][m]) == np.sign(lfc_ref[t][m])).sum()); n += int(m.sum())
        return round(hit / max(n, 1), 4), n
    pred_oracle = e[None, :] * d_t[:, None]
    pred_typ = np.broadcast_to(e * np.median(d_all), lfc_ref.shape)
    terc = np.argsort(np.argsort(np.abs(d_t))) * 3 // len(targets)
    acc = {}
    for name, pred in (("eps x delta_t (oracle)", pred_oracle), ("eps x median delta (target-independent)", pred_typ)):
        acc[name] = {"all": acc_pred(pred, np.ones(kg.size, bool)), "block": acc_pred(pred, block),
                     "nonblock": acc_pred(pred, ~block),
                     **{f"|delta| tercile {k}": acc_pred(pred, np.ones(kg.size, bool), np.flatnonzero(terc == k)) for k in range(3)}}
        log(f"{name}: " + "; ".join(f"{k} {v[0]:.4f} (n {v[1]})" for k, v in acc[name].items()))
    out["h1_accuracy"] = acc
    # relation to the programmes (programme sign: + = up)
    zA = 0.5 - P_arc; zC = 0.5 - P_c4
    s_typ = e * np.sign(np.median(d_all))
    exp_ = cpm >= 5
    out["h1_eps_vs_programmes"] = {"r(eps*sign(median delta), Arc 0.5-P)": pearson(s_typ, zA, exp_),
                                   "r(eps*sign(median delta), CLEAN4 0.5-P)": pearson(s_typ, zC, exp_),
                                   "r(eps, mean truth lfc over targets)": pearson(e, np.nanmean(np.where(gate, lfc_ref, np.nan), 0), exp_)}
    log("eps vs programmes: " + json.dumps({k: round(v, 3) for k, v in out["h1_eps_vs_programmes"].items()}))
    # does it add to a programme? z-sum accuracy on gate genes, target-independent and oracle forms
    def zs(v):
        return v / np.nanstd(v)
    add = {}
    for pname, zp in (("Arc", zA), ("CLEAN4", zC)):
        base_pred = np.broadcast_to(np.where(np.isfinite(zp), zp, np.nan), lfc_ref.shape)
        add[pname] = {"alone": acc_pred(base_pred, np.ones(kg.size, bool))[0]}
        for lam in (0.25, 0.5, 1.0):
            comb_t = np.broadcast_to(zs(zp) + lam * zs(s_typ), lfc_ref.shape)
            comb_o = zs(zp)[None, :] + lam * pred_oracle / np.nanstd(pred_oracle)
            add[pname][f"+{lam} eps (typ)"] = acc_pred(comb_t, np.ones(kg.size, bool))[0]
            add[pname][f"+{lam} eps x delta_t (oracle)"] = acc_pred(comb_o, np.ones(kg.size, bool))[0]
        log(f"{pname} + eps: " + json.dumps(add[pname]))
    out["h1_additivity_gate_acc"] = add
    # share of each target's gate-gene truth lfc variance explained by eps x delta_t (least squares through 0)
    r2 = []
    for t in range(len(targets)):
        m = gate[t]
        if m.sum() < 10:
            continue
        x, y = pred_oracle[t][m], lfc_ref[t][m]
        bta = (x * y).sum() / max((x * x).sum(), 1e-12)
        r2.append((targets[t], float(d_t[t]), float(1 - ((y - bta * x) ** 2).sum() / (y * y).sum()), float(bta)))
    r2a = np.array([x[2] for x in r2])
    out["h1_r2_eps_x_delta"] = dict(median=float(np.median(r2a)), q90=float(np.quantile(r2a, .9)), max=float(r2a.max()),
                                    top=[x for x in sorted(r2, key=lambda x: -x[2])[:8]])
    log(f"R2 of truth lfc on eps x delta_t per target (gate genes): median {np.median(r2a):.4f}, q90 {np.quantile(r2a, .9):.4f}, "
        f"max {r2a.max():.4f}; top {[(a, round(b, 3), round(c, 3)) for a, b, c, _ in out['h1_r2_eps_x_delta']['top'][:5]]}")
    np.savez_compressed(HERE / "work/tech517_h1.npz", eps=e, eps_h0=eps["eps_h0"], eps_h1=eps["eps_h1"], genes_keep=kg,
                        delta_h1test=d_t, targets=np.asarray(targets))
    # ------------------------------------------------ 3. A/B/C controls: eps per context; NTC guides
    official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    opos = {g: i for i, g in enumerate(official)}
    h1_on_off = np.array([opos.get(g, -1) for g in kg])
    ctx_eps = {}
    for c in "ABC":
        it = CSRStream(ROOT / f"data/challenge/context_{c}.h5ad")
        lib = it.libraries()
        guide = cat(it.f["obs/ntc_id"])
        eps_c = eps_from_controls(it, np.arange(it.n), lib, None, rng)
        lnL = np.log(lib)
        ug = np.unique(guide)
        gm = np.array([lnL[guide == g].mean() for g in ug]); gn = np.array([(guide == g).sum() for g in ug])
        within = np.mean([lnL[guide == g].var(ddof=1) for g in ug])
        F = float(gm.var(ddof=1) / (within / gn.mean()))
        # guide pseudobulk deviations (ln mean CPM per guide minus the context's), top PC, and alignment with eps
        glab = np.searchsorted(ug, guide)
        Sg = it.group_sums(glab, len(ug), cols=None, cpm=True, lib=lib)
        it.close()
        Ag = Sg / gn[:, None]
        expd = eps_c["mean_cpm"] >= 5
        Y = np.log(np.maximum(Ag[:, expd], 1e-3)); Y -= Y.mean(0)
        u, sv, vt = np.linalg.svd(Y, full_matrices=False)
        pc1 = vt[0]
        r_pc1_eps = pearson(pc1, eps_c["eps"][expd])
        r_pc1score_lib = float(np.corrcoef(u[:, 0] * sv[0], gm)[0, 1])
        # the eps direction's projection of the guide deviations vs the guides' library means
        proj = Y @ (eps_c["eps"][expd] / np.linalg.norm(eps_c["eps"][expd]))
        r_proj_lib = float(np.corrcoef(proj, gm)[0, 1])
        # null for the eps-projection correlation: permute guide labels 200 x (sampling noise alone)
        ctx_eps[c] = eps_c
        out[f"ctx_{c}"] = dict(n=int(it.n), median_L=float(np.exp(np.median(lnL))), sd_lnL=float(lnL.std()),
                               eps_split_half_r=pearson(eps_c["eps_h0"], eps_c["eps_h1"], expd),
                               guide_F_lnL=F, guide_mean_lnL_sd=float(gm.std()), expected_sd=float(np.sqrt(within / gn.mean())),
                               guide_pc1_var_share=float(sv[0] ** 2 / (sv ** 2).sum()), r_pc1_eps=r_pc1_eps,
                               r_pc1score_guide_lnL=r_pc1score_lib, r_epsproj_guide_lnL=r_proj_lib)
        log(f"{c}: median L {out[f'ctx_{c}']['median_L']:.0f}, eps split-half r {out[f'ctx_{c}']['eps_split_half_r']:.3f}; "
            f"NTC guides: F(ln L) {F:.2f} (sd of guide means {gm.std():.4f} vs sampling {np.sqrt(within / gn.mean()):.4f}); "
            f"guide PC1 share {out[f'ctx_{c}']['guide_pc1_var_share']:.3f}, r(PC1, eps) {r_pc1_eps:+.3f}, "
            f"r(PC1 score, guide ln L) {r_pc1score_lib:+.3f}, r(eps-projection, guide ln L) {r_proj_lib:+.3f}")
    # eps agreement across contexts (expressed >= 5 CPM in both), H1 mapped to the official axis
    E = {"H1": np.full(len(official), np.nan)}
    okh = h1_on_off >= 0
    E["H1"][h1_on_off[okh]] = e[okh]
    for c in "ABC":
        E[c] = np.where(ctx_eps[c]["mean_cpm"] >= 5, ctx_eps[c]["eps"], np.nan)
    Eh1 = np.full(len(official), np.nan); Eh1[h1_on_off[okh]] = np.where(eps["mean_cpm"][okh] >= 5, e[okh], np.nan)
    E["H1"] = Eh1
    agree = {f"{a}-{b}": pearson(E[a], E[b]) for i, a in enumerate(E) for b in list(E)[i + 1:]}
    out["eps_cross_context_r"] = agree
    log("eps agreement across contexts: " + json.dumps({k: round(v, 3) for k, v in agree.items()}))
    np.savez_compressed(HERE / "work/tech517_ctx.npz", **{f"eps_{k}": v for k, v in E.items()})
    (HERE / "work/tech517.json").write_text(json.dumps(out, indent=1, default=float))
    log("done")


if __name__ == "__main__":
    main()
