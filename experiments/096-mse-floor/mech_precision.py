"""The tail policy's whole risk, reduced to one measurable scalar.

`fid` factorises as  F_p = k/max(n_pred, n_real) = (k/n_pred) * min(1, n_pred/n_real).
The COVERAGE leg carries `n_real` and is unreadable on our harnesses (RESEARCH-METHOD rule 1).
The PRECISION leg `k/n_pred` does not carry `n_real` at all -- it needs only our own call set and
the reference's signs on those genes, both of which the harnesses have. And on the challenge
`n_pred` 4772 >= `n_real` ~4757 already saturates coverage, and the tail policy raises `n_pred`
to 5417, so coverage STAYS saturated and precision is the only thing that can move `fid`.

So: emit the same profile under both policies, run the exact DE path against the real controls,
and measure the directional precision of the calls that are NEW under the tiered policy.
"""
import os, sys, time
import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix, issparse
from scipy.stats import false_discovery_control, norm as _n

ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src")
sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
sys.path.insert(0, f"{ROOT}/experiments/092-call-set-truncation")
from de_kernel import ctrl_tiesums, wilcoxon_z          # noqa: E402
from msekit import CTRL, TS, group_stats                # noqa: E402
from cell_eval2.prep import bulk_lognorm_means          # noqa: E402
from vcc2026.donor_transfer import emit_counts          # noqa: E402

MIN_CPM, ALPHA, NCELL = 5.0, 0.05, 400
NP = int(os.environ.get("NP", "25"))

def norm_log(x):
    s = x.sum(1); s[s == 0] = 1.0
    return np.log1p(x * (1e6 / s[:, None]))

def call_set(cs, tie, block):
    z = wilcoxon_z(cs, tie, np.ascontiguousarray(block.T))
    pv = 2.0 * _n.sf(np.abs(z))
    return false_discovery_control(np.nan_to_num(pv, nan=1.0), method="bh") < ALPHA

HARNESS = [("h1-test", f"{ROOT}/data/harness/h1-test.h5ad",
            f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
           ("orion-thinned", f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
            f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad")]

for name, real, pred in HARNESS:
    t0 = time.time()
    A = ad.read_h5ad(real)
    lab = np.asarray(A.obs["target_gene"].astype(str).values)
    genes = np.asarray([str(g) for g in A.var_names]); G = genes.size
    X = A.X.tocsr() if issparse(A.X) else csr_matrix(A.X)
    ci = np.flatnonzero(lab == CTRL)
    C = norm_log(np.asarray(X[ci].todense(), dtype=np.float64))
    cpm_ctrl = np.expm1(C).mean(0); keep = cpm_ctrl > MIN_CPM
    cs = np.ascontiguousarray(np.sort(C[:, keep].T, axis=1)); tie = ctrl_tiesums(cs)
    ctrl_lin = np.expm1(C[:, keep]).mean(0); del C
    Gk = int(keep.sum()); gk = genes[keep]; gpos = {g: i for i, g in enumerate(gk)}

    pp, ps, pc, pj, pg = group_stats(pred)
    pm = bulk_lognorm_means(ps.astype(np.float64), TS)
    rp, rs, rc, rj, rg = group_stats(real)
    bctrl = bulk_lognorm_means(rs[int(np.flatnonzero(rp == CTRL)[0])][None, :].astype(np.float64), TS)[0]
    pidx = {str(p): i for i, p in enumerate(pp)}
    targets = [t for t in sorted(set(lab) - {CTRL}) if t in pidx][:NP]
    print(f"\n===== {name}: {Gk}/{G} genes past the gate, {len(targets)} perturbations "
          f"[{time.time()-t0:.0f}s] =====", flush=True)

    acc = {k: [] for k in ("n0", "n1", "p0", "p1", "pnew", "pdrop", "nnew", "ndrop",
                           "hit0", "hitnew", "j0", "j1")}
    for t in targets:
        rows = np.flatnonzero(lab == t)
        Preal = norm_log(np.asarray(X[rows].todense(), dtype=np.float64))
        Rsig = call_set(cs, tie, Preal[:, keep])
        rlfc = np.log2(np.expm1(Preal[:, keep]).mean(0) + 1e-9) - np.log2(ctrl_lin + 1e-9)
        del Preal
        i = pidx[t]
        e = pm[i] - bctrl
        prof = np.expm1(bctrl + e); np.clip(prof, 0, None, out=prof); prof /= prof.sum()
        libs = np.full(NCELL, ps[i].sum() / NCELL)
        sets = {}
        for pol in ("a10", "tiered"):
            sp = np.full(G, 0.15); ph = np.full(G, 16.0)
            if pol == "tiered":
                top = np.argpartition(-np.abs(e), 8000)[:8000]
                sp[:] = 1.0; ph[:] = 4.0; sp[top] = 0.15; ph[top] = 16.0
            v, col, ip = emit_counts(prof, libs, 0.0, np.random.default_rng(8000 + i),
                                     mode="quasi", spread=sp, overdispersion=ph)
            M = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                           shape=(NCELL, G))
            Pn = norm_log(np.asarray(M.todense(), dtype=np.float64))
            sig = call_set(cs, tie, Pn[:, keep])
            plfc = np.log2(np.expm1(Pn[:, keep]).mean(0) + 1e-9) - np.log2(ctrl_lin + 1e-9)
            sets[pol] = (sig, plfc)
        ok = np.ones(Gk, bool)
        if t in gpos: ok[gpos[t]] = False                       # the scorer's own-gene exclusion
        ind = ok & (rlfc != 0) & ~np.isnan(rlfc)                # scored universe = reference's
        def prec(mask, plfc):
            m = mask & ind
            if m.sum() == 0: return np.nan, 0
            return float((np.sign(plfc[m]) == np.sign(rlfc[m])).mean()), int(m.sum())
        s0, l0 = sets["a10"]; s1, l1 = sets["tiered"]
        p0, n0 = prec(s0, l0); p1, n1 = prec(s1, l1)
        pn, nn = prec(s1 & ~s0, l1); pd, nd = prec(s0 & ~s1, l0)
        acc["n0"].append(n0); acc["n1"].append(n1)
        acc["p0"].append(p0); acc["p1"].append(p1)
        acc["pnew"].append(pn); acc["pdrop"].append(pd)
        acc["nnew"].append(nn); acc["ndrop"].append(nd)
        R = Rsig & ok
        acc["hit0"].append(float((s0 & R & ind).sum()) / max(n0, 1))
        acc["hitnew"].append(float(((s1 & ~s0) & R & ind).sum()) / max(nn, 1))
        for tag, s in (("j0", s0), ("j1", s1)):
            u = int((R | (s & ok)).sum())
            acc[tag].append(1.0 if u == 0 else int((R & s & ok).sum()) / u)
    m = {k: float(np.nanmean(v)) for k, v in acc.items()}
    print(f"  n_pred          a10 {m['n0']:7.0f}   tiered {m['n1']:7.0f}   ({m['n1']/m['n0']:.3f}x)")
    print(f"  precision k/n_pred    a10 {m['p0']:.4f}      tiered {m['p1']:.4f}     "
          f"delta {m['p1']-m['p0']:+.4f}")
    print(f"  >>> precision of the {m['nnew']:.0f} NEW calls   {m['pnew']:.4f}   "
          f"(chance = 0.5; a10's own = {m['p0']:.4f})")
    print(f"      precision of the {m['ndrop']:.0f} DROPPED calls {m['pdrop']:.4f}")
    print(f"  hit rate into the reference's DE set   a10 {m['hit0']:.4f}   new calls {m['hitnew']:.4f}")
    print(f"  jaccard         a10 {m['j0']:.4f}   tiered {m['j1']:.4f}   delta {m['j1']-m['j0']:+.4f}")
