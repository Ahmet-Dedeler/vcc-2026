"""518 — H1 (h1-test, 85 targets) as the third surface: the mse identity EXACTLY, with per-gene jackknife noise correction.

Prediction: 435's harness stats `base_s0_bps_prodem_mde0.1cap0.1min0.3coarse` (410b-construction bulk realised by ballast,
coarse MDE, seed 0), the only scored H1 bulk on disk. Truth: 323's reference.npz (85 x 400 cells + 18,400 NTC).
KNOWN ANSWERS (abort): (i) my legs N_p, D_p == msekit's (cell_eval2 primitives) per target (1e-12 rel) and the raw ==
0.9253584049461503; (ii) the per-gene delete-1 jackknife I stream from the cells sums to the reference's jk per target and
for the control (1e-9 rel).
Decomposition (all in the metric's units, per-row weights w_p = 1/|G_p|):
    raw = 1 + [A - 2X - rho*claim + Kc] / sumD      (A, X over each row's scored genes; Kc = sum w_p C_ctrl)
    per target: s_p = |h_p|, x_p = <h_p, delta_p> (unbiased for <h_p, delta*_p>), m_p = x_p/s_p,
                noise var of m_p: v_p = sum_g hhat_pg^2 (C_pg + C_ctrl,g);  |delta*_p|^2 = |delta_p|^2 - C_p - C_ctrl.
    magnitude oracle at pinned spread energy: X* = sqrt(sum s^2 * sum (m+^2 - v)) (noise-corrected) vs X.
    gene classes: own target (excluded), the other 84 H1 targets, the 300 A/B/C panel genes, ungated (< 10 CPM in the H1
    control: the builder's bulk is zero there), expression deciles of the gated genes -> truth signal share (noise-
    corrected), our energy share, X share, class cosine.
    mean component: n|mean delta*|^2 / sum|delta*|^2 (noise-corrected, the control noise does not average out) and
    cos(mean h, mean delta*).

    uv run python experiments/518-bulk-mse-gap/h1_518.py        (~1.5 GB; streams the 676 MB h5ad once, caches the jk)
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import h5py  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))
import msekit  # noqa: E402
from cell_eval2.prep import bulk_lognorm_means  # noqa: E402

TS = 50_000.0
STATS = ROOT / "experiments/435-zs-on-coarse/work/h1/stats-base_s0_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
REF = ROOT / "experiments/323-current-family-validation/work/reference.npz"
H5 = ROOT / "data/harness/h1-test.h5ad"
RAW_KA = 0.9253584049461503
JK_CACHE = WORK / "h1_jk_pergene.npz"


def load():
    z = np.load(STATS)
    r = np.load(REF)
    targets = [str(t) for t in r["targets"]]
    assert [str(t) for t in z["targets"]] == targets and str(r["perts"][-1]) == "non-targeting"
    genes = np.asarray([str(g) for g in r["genes"]])
    return dict(targets=targets, genes=genes, pred_sums=z["sums"].astype(np.float64), pred_jk=z["jk"].astype(np.float64),
                ref_sums=r["sums"].astype(np.float64), ref_counts=r["counts"], ref_jk=r["jk"].astype(np.float64),
                perts=r["perts"], ref_sums_int=r["sums"], pred_sums_int=z["sums"])


def pergene_jk(D):
    """Delete-1 jackknife per (group, gene) from the cells, streamed; cached. Groups = 85 targets + control (last)."""
    if JK_CACHE.exists():
        c = np.load(JK_CACHE)
        return c["C"]
    t0 = time.time()
    with h5py.File(H5, "r") as h:
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in h["obs"]["target_gene"]["categories"][:]]
        codes = h["obs"]["target_gene"]["codes"][:]
        node = h["var"]["_index"]
        raw_g = node["values"][:] if isinstance(node, h5py.Group) else node[:]      # nullable-string layout
        genes = [g.decode() if isinstance(g, bytes) else str(g) for g in raw_g]
        assert np.array_equal(np.asarray(genes), D["genes"])
        grp_names = D["targets"] + ["non-targeting"]
        cat_to_grp = np.asarray([grp_names.index(c) for c in cats])
        grp = cat_to_grp[codes]
        indptr = h["X"]["indptr"][:]
        P = D["ref_sums"]                                  # (86, G) group sums
        S = P.sum(1)
        vfull = np.log1p(TS * P / S[:, None])
        G = P.shape[1]
        s1 = np.zeros_like(P)
        s2 = np.zeros_like(P)
        n = np.bincount(grp, minlength=len(grp_names)).astype(float)
        ncell = codes.size
        step = 2000
        for a in range(0, ncell, step):
            b = min(ncell, a + step)
            lo, hi = indptr[a], indptr[b]
            data = h["X"]["data"][lo:hi].astype(np.float64)
            idx = h["X"]["indices"][lo:hi]
            Y = np.zeros((b - a, G))
            rows = np.repeat(np.arange(b - a), np.diff(indptr[a:b + 1]))
            Y[rows, idx] = data
            gi = grp[a:b]
            lib = Y.sum(1)
            V = np.log1p(TS * (P[gi] - Y) / (S[gi] - lib)[:, None]) - vfull[gi]
            np.add.at(s1, gi, V)
            np.add.at(s2, gi, V * V)
            if (a // step) % 5 == 0:
                print(f"  jk stream {b}/{ncell} ({time.time() - t0:.0f}s)", flush=True)
        C = (n[:, None] - 1) / n[:, None] * (s2 - s1 * s1 / n[:, None])
    np.savez_compressed(JK_CACHE, C=C, n=n)
    return C


def main():
    t0 = time.time()
    D = load()
    targets, genes = D["targets"], D["genes"]
    n = len(targets)
    gpos = {g: i for i, g in enumerate(genes)}
    # ---------------- KA (i): the scorer's legs through msekit, and my own legs
    pred_sums = np.vstack([D["pred_sums_int"], D["ref_sums_int"][-1]])
    pred_jk = np.append(D["pred_jk"], D["ref_jk"][-1])
    mse, num, den = msekit.score(D["perts"], pred_sums, D["ref_counts"], pred_jk, D["perts"], D["ref_sums_int"],
                                 D["ref_counts"], D["ref_jk"], genes)
    assert abs(mse - RAW_KA) < 1e-12, f"KA raw {mse}"
    mu_p = bulk_lognorm_means(pred_sums.astype(np.float64), TS)
    mu_r = bulk_lognorm_means(D["ref_sums_int"].astype(np.float64), TS)
    b_ctrl = mu_r[-1]
    h = mu_p[:n] - b_ctrl
    dl = mu_r[:n] - b_ctrl
    own = np.asarray([gpos.get(t, -1) for t in targets])
    Gp = np.where(own >= 0, genes.size - 1, genes.size).astype(float)
    wp = 1.0 / Gp
    mask = np.ones((n, genes.size), bool)
    mask[np.arange(n)[own >= 0], own[own >= 0]] = False
    Cp, Cc, Chat = D["ref_jk"][:n], D["ref_jk"][-1], D["pred_jk"]
    # the scorer's budget B: sum_g sum_{p in R_g} w_p (mu_pg - wmean)^2
    W = mask * wp[:, None]
    wm = (W * mu_p[:n]).sum(0) / W.sum(0)
    B = float((W * (mu_p[:n] - wm[None, :]) ** 2).sum())
    term = np.minimum(Chat, 1.0 * Cp)
    claim = float((wp * term).sum())
    rho = min(1.0, B / claim)
    dist = ((h - dl) ** 2 * mask).sum(1)
    N = (dist - rho * term - Cp) / Gp
    Dd = ((dl ** 2) * mask).sum(1)
    Den = (Dd - Cp - Cc) / Gp
    dN = max(abs(N[i] - num[t]) / abs(num[t]) for i, t in enumerate(targets))
    dD = max(abs(Den[i] - den[t]) / abs(den[t]) for i, t in enumerate(targets))
    assert dN < 1e-9 and dD < 1e-9, f"KA legs {dN} {dD}"
    raw = N.sum() / Den.sum()
    print(f"KA (i) PASS: msekit raw {mse:.13f} == harness; my legs rel |d| N {dN:.1e} D {dD:.1e}; mine {raw:.13f}")
    # ---------------- identity terms
    A = float((wp * (h ** 2 * mask).sum(1)).sum())
    X = float((wp * (h * dl * mask).sum(1)).sum())
    Kc = float((wp * Cc).sum())
    sumD = float(Den.sum())
    raw_id = 1 + (A - 2 * X - rho * claim + Kc) / sumD
    assert abs(raw_id - raw) < 1e-10, f"identity {raw_id} vs {raw}"
    meanrow = A - B   # approx (the scorer centres per gene over R_g)
    ident = dict(A=A, B=B, claim=claim, rho=rho, Kc=Kc, X=X, sumD=sumD, raw=raw, A_over_sumD=A / sumD,
                 twoX_over_sumD=2 * X / sumD, credit_over_sumD=rho * claim / sumD, Kc_over_sumD=Kc / sumD,
                 meanrow_over_sumD=meanrow / sumD, cos_agg=X / np.sqrt(A * sumD), energy_opt_scale=X / A)
    print("H1 identity: " + " ".join(f"{k} {v:.4f}" for k, v in ident.items()))
    # the level: raw(s) with the prediction's effect scaled by s in log space (credit via B s^2 up to the claim)
    ss = np.linspace(0.2, 3.0, 281)
    rs = [1 + (s * s * A - 2 * s * X - min(claim, s * s * B) + Kc) / sumD for s in ss]
    i = int(np.argmin(rs))
    ident["level_scan"] = dict(best_s=float(ss[i]), best_raw=float(rs[i]), cliff_s=float(np.sqrt(claim / B)))
    print(f"H1 level: optimum s {ss[i]:.3f} raw {rs[i]:.4f} (cliff at s {np.sqrt(claim / B):.3f}); now s 1 raw {raw:.4f}")
    # ---------------- KA (ii): per-gene jackknife
    C = pergene_jk(D)
    dj = np.max(np.abs(C.sum(1) - D["ref_jk"]) / D["ref_jk"])
    assert dj < 1e-9, f"KA jk per gene {dj}"
    print(f"KA (ii) PASS: per-gene jackknife sums to the reference jk (max rel |d| {dj:.1e}) ({time.time() - t0:.0f}s)")
    Cpg, Ccg = C[:n], C[-1]
    # ---------------- per-target magnitude / direction
    hm = h * mask
    s = np.sqrt((hm ** 2).sum(1))
    x = (hm * dl).sum(1)
    m = x / s
    hh = hm / s[:, None]
    v = ((hh ** 2) * (Cpg + Ccg[None, :])).sum(1)
    tstar2 = Dd - Cp - Cc                                   # |delta*_p|^2 over G_p (noise-corrected)
    # at pinned (weighted) energy: X = sum w s m; oracle X* = sqrt(sum w s^2 * sum w m+^2) (w equal up to 1/G_p)
    mp2_raw = np.maximum(m, 0) ** 2
    mp2_cor = np.maximum(m ** 2 - v, 0) * (m > 0)
    Xn = float((wp * s * m).sum())
    Xstar_raw = float(np.sqrt((wp * s * s).sum() * (wp * mp2_raw).sum()))
    Xstar_cor = float(np.sqrt((wp * s * s).sum() * (wp * mp2_cor).sum()))
    cosp = m / np.sqrt(np.maximum(tstar2, 1e-12))
    per = dict(r_mag_raw=Xstar_raw / Xn - 1, r_mag_corrected=Xstar_cor / Xn - 1, frac_neg=float((m < 0).mean()),
               noise_share_of_m2=float(v.sum() / (m ** 2).sum()),
               energy_weighted_rms_cos=float(np.sqrt((np.maximum(tstar2, 0) * np.clip(cosp, -1, 1) ** 2).sum() / np.maximum(tstar2, 0).sum())))
    # energy deciles by the truth's noise-corrected energy
    order = np.argsort(-tstar2)
    top = order[: max(1, n // 10)]
    per["top_decile_truth_share"] = float(np.maximum(tstar2[top], 0).sum() / np.maximum(tstar2, 0).sum())
    per["top_decile_X_share"] = float(x[top].sum() / x.sum())
    per["top_decile_our_energy_share"] = float((s[top] ** 2).sum() / (s ** 2).sum())
    lam = np.sqrt((s * s).sum() / mp2_cor.sum())
    ratio = lam * np.sqrt(mp2_cor) / s
    per["ideal_ratio_q"] = np.quantile(ratio, [0.1, 0.25, 0.5, 0.75, 0.9]).tolist()
    per["ideal_ratio_top_decile_median"] = float(np.median(ratio[top]))
    per["spearman_ideal_vs_truth_energy"] = float(np.corrcoef(np.argsort(np.argsort(ratio)), np.argsort(np.argsort(tstar2)))[0, 1])
    per["spearman_ideal_vs_our_energy"] = float(np.corrcoef(np.argsort(np.argsort(ratio)), np.argsort(np.argsort(s)))[0, 1])
    per["spearman_our_energy_vs_truth_energy"] = float(np.corrcoef(np.argsort(np.argsort(s)), np.argsort(np.argsort(tstar2)))[0, 1])
    print("H1 per-target: " + json.dumps({k: (round(v, 4) if isinstance(v, float) else [round(u, 3) for u in v]) for k, v in per.items()}))
    # ---------------- gene classes
    cpm = 1e6 * np.expm1(b_ctrl) / np.expm1(b_ctrl).sum()
    gated = (cpm > 10.0) & (b_ctrl > 0)
    h1panel = np.zeros(genes.size, bool)
    h1panel[own[own >= 0]] = True
    abc_panel = np.zeros(genes.size, bool)
    for t in (l.split(",")[0] for l in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if l.strip()):
        if t in gpos:
            abc_panel[gpos[t]] = True
    sig = (dl ** 2 - Cpg - Ccg[None, :]) * mask              # noise-corrected truth signal per (p, g)
    sig_tot = sig.sum()
    ours = (h ** 2) * mask
    cross = h * dl * mask
    classes = {"other_H1_targets": h1panel, "ABC_panel_genes": abc_panel & ~h1panel, "ungated(<10cpm)": ~gated,
               "gated": gated}
    qs = np.quantile(b_ctrl[gated], np.linspace(0, 1, 11))
    dec = np.clip(np.searchsorted(qs, b_ctrl, side="right") - 1, 0, 9)
    for k in range(10):
        classes[f"gated_expr_decile_{k}"] = gated & (dec == k)
    gc = {}
    for name, gm in classes.items():
        cols = gm[None, :]
        ts = float((sig * cols).sum() / sig_tot)
        oe = float((ours * cols).sum() / ours.sum())
        xs = float((cross * cols).sum() / cross.sum())
        tn = float(np.sqrt(max((sig * cols).sum(), 1e-300)))
        on = float(np.sqrt((ours * cols).sum()))
        ccos = float((cross * cols).sum() / max(tn * on, 1e-300))
        gc[name] = dict(n_genes=int(gm.sum()), truth_signal_share=ts, our_energy_share=oe, X_share=xs, class_cos=ccos)
    print(f"{'gene class':28s}{'genes':>7}{'truth sig':>10}{'ours':>8}{'X':>8}{'cos':>8}")
    for k, g_ in gc.items():
        print(f"{k:28s}{g_['n_genes']:7d}{100 * g_['truth_signal_share']:9.1f}%{100 * g_['our_energy_share']:7.1f}%"
              f"{100 * g_['X_share']:7.1f}%{g_['class_cos']:8.3f}")
    # value of predicting a class we leave at zero, at pinned energy, if predicted at cosine c:  X ratio sqrt(1 + phi c^2/(X-part))
    # ---------------- mean component
    dbar = (dl * mask).sum(0) / mask.sum(0)
    hbar = (h * mask).sum(0) / mask.sum(0)
    noise_bar = (Cpg.sum(0) / n ** 2) + Ccg
    dbar_sig = float((dbar ** 2 - noise_bar).sum())
    mean_share = n * dbar_sig / float(tstar2.sum())
    mean = dict(truth_mean_share=mean_share, our_mean_share=float(n * (hbar ** 2).sum() / (s ** 2).sum()),
                cos_mean=float((hbar * dbar).sum() / np.sqrt((hbar ** 2).sum() * max(dbar_sig, 1e-300))),
                pot_raw=float(n * ((hbar - dbar) ** 2 - noise_bar).sum() * np.mean(wp) / sumD))
    print("H1 mean component: " + " ".join(f"{k} {v:.4f}" for k, v in mean.items()))
    out = dict(identity=ident, per_target=per, gene_classes=gc, mean_component=mean,
               arrays=dict(s=s.tolist(), m=m.tolist(), v=v.tolist(), tstar2=tstar2.tolist(), targets=targets))
    (WORK / "h1_518.json").write_text(json.dumps(out, indent=1, default=float))
    np.savez_compressed(WORK / "h1_518_rows.npz", h=h, dl=dl, mask=mask, Cpg=Cpg, Ccg=Ccg, b_ctrl=b_ctrl, wp=wp,
                        Chat=Chat, Cp=Cp, Cc=Cc)
    print(f"done ({time.time() - t0:.0f}s)")


if __name__ == "__main__":
    main()
