"""Does the dispersion the recipe adds land where nothing was callable anyway?

The recipe raises `C_hat_p` to the #247 cap by over-dispersing a TAIL of genes. That lowers
those genes' effect-to-dispersion ratio, which is exactly what sets the Wilcoxon p-values behind
`fid`, `jac` and `reach`'s ordering tier. This runs the REAL DE path -- production `emit_counts`,
then the exact scanpy-identical rank-sum against the REAL challenge control cells, then BH -- and
counts the call set under each policy. No proxy.

Two tail keys are compared, because they are NOT the same set:
    top-K by |predicted log-effect|   (the recipe as first written)
    control CPM >= T                  (Zhou Li's diagnosis: callability is set by expression)
"""
import os, sys, time
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix
from scipy.stats import false_discovery_control

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from de_kernel import ctrl_tiesums, wilcoxon_z          # noqa: E402
from cell_eval2.moments import jackknife_correction     # noqa: E402
import mech_alpha_chal as M                             # noqa: E402

CTXS = os.environ.get("CTXS", "A").split(",")
NP    = int(os.environ.get("NP", "15"))
MIN_CPM, ALPHA = 5.0, 0.05

def norm_log(x):
    s = x.sum(1); s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))

D = M.load(f"{ROOT}/experiments/096-mse-floor/work/v11-stats.npz")
cpn = dict(np.load(f"{ROOT}/experiments/096-mse-floor/work/cp-proxy-ntc.npz").items())

for ctx in CTXS:
    t0 = time.time()
    A = ad.read_h5ad(f"{ROOT}/data/challenge/context_{ctx}.h5ad")
    Xc = np.asarray(A.X.todense(), dtype=np.float64)
    C = norm_log(Xc)
    cpm_ctrl = np.expm1(C).mean(0)
    keep = cpm_ctrl > MIN_CPM
    cs = np.ascontiguousarray(np.sort(C[:, keep].T, axis=1)); tie = ctrl_tiesums(cs)
    del C, Xc
    G = keep.size; Gk = int(keep.sum())
    print(f"\n===== context {ctx}: {Gk}/{G} genes past the DE gate (control CPM > {MIN_CPM}) "
          f"[{time.time()-t0:.0f}s] =====", flush=True)

    d = D[ctx]
    idx = np.linspace(0, len(d["tgt"]) - 1, NP).astype(int)
    Cp = float(cpn[ctx].mean())

    def policy_vectors(e, name):
        sp = np.full(G, 0.15); ph = np.full(G, 16.0)
        if name == "a10":
            pass
        elif name.startswith("topK"):
            K = int(name[4:]); top = np.argpartition(-np.abs(e), K)[:K]
            sp[:] = 1.0; ph[:] = 4.0; sp[top] = 0.15; ph[top] = 16.0
        elif name.startswith("cpm"):
            T = float(name[3:]); hi = cpm_ctrl >= T
            sp[:] = 1.0; ph[:] = 4.0; sp[hi] = 0.15; ph[hi] = 16.0
        return sp, ph

    POL = ["a10", "topK8000", "cpm40", "cpm20", "cpm10"]
    res = {p: dict(npred=[], chat=[], sig=[]) for p in POL}
    for j in idx:
        e = d["dh"][j]
        prof = np.expm1(d["bctrl"] + e); np.clip(prof, 0, None, out=prof); prof /= prof.sum()
        libs = np.full(400, d["S"][j] / 400)
        for p in POL:
            sp, ph = policy_vectors(e, p)
            v, col, ip = M.emit_counts(prof, libs, 0.0, np.random.default_rng(7000 + j),
                                       mode="quasi", spread=sp, overdispersion=ph)
            Xm = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                            shape=(400, G))
            res[p]["chat"].append(
                jackknife_correction(Xm, np.zeros(400, np.intp), 1, M.TS, chunk=64)[0])
            P = norm_log(np.asarray(Xm.todense(), dtype=np.float64))
            z = wilcoxon_z(cs, tie, np.ascontiguousarray(P[:, keep].T))
            from scipy.stats import norm as _n
            pv = 2.0 * _n.sf(np.abs(z))
            padj = false_discovery_control(np.nan_to_num(pv, nan=1.0), method="bh")
            sig = padj < ALPHA
            res[p]["npred"].append(int(sig.sum())); res[p]["sig"].append(sig)

    base = np.stack(res["a10"]["sig"])
    print(f"  {'policy':>10}{'n_pred':>9}{'vs a10':>9}{'C_hat':>9}{'/C_p':>7}{'claim':>9}"
          f"{'kept calls':>12}")
    for p in POL:
        n = float(np.mean(res[p]["npred"])); n0 = float(np.mean(res["a10"]["npred"]))
        ch = float(np.mean(res[p]["chat"]))
        claim = 300 * float(np.mean(np.minimum(res[p]["chat"], Cp))) / 18532
        S = np.stack(res[p]["sig"])
        kept = float(np.mean([(S[i] & base[i]).sum() / max(base[i].sum(), 1)
                              for i in range(S.shape[0])]))
        print(f"  {p:>10}{n:>9.0f}{n/n0:>9.3f}{ch:>9.2f}{ch/Cp:>7.2f}{claim:>9.5f}{kept:>12.3f}",
              flush=True)

    # Zhou Li's diagnosis, on our own data: is callability set by control expression?
    S0 = base.mean(0) > 0.5
    full = np.zeros(G, bool); full[np.flatnonzero(keep)] = S0
    print(f"  callability vs control expression (a10 emission, genes past the gate):")
    print(f"    called   n={int(S0.sum()):6d}  median control CPM {np.median(cpm_ctrl[keep][S0]):8.1f}")
    print(f"    uncalled n={int((~S0).sum()):6d}  median control CPM {np.median(cpm_ctrl[keep][~S0]):8.1f}")
    for T in (10, 20, 40):
        hi = cpm_ctrl >= T
        print(f"    control CPM >= {T:3d}: {int(hi.sum()):6d} genes, holding "
              f"{100*np.mean(full[hi]):5.1f}% called; they are "
              f"{100*full[hi].sum()/max(full.sum(),1):5.1f}% of ALL called genes")
