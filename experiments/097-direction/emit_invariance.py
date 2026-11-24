"""Rule-3 check: is X really invariant to the emission policy?

MSE-MECHANISM.md asserts "`A`, `B`, `X`, `pds` and `nmae` are all invariant to (spread, phi) by
construction" because the quasi emitter repairs each gene's panel total exactly.  a12's entire
plan rests on that (it raises `claim` by 2x through emission and assumes X is untouched).  My
realisation of the orion base at spread 0.5/phi 4 returned X = 0.005072 against the cached
prediction's 0.008136 -- a 38 % gap -- so the assertion needs testing rather than citing.

Two things are confounded in that gap and this separates them:
  (1) the emission policy (spread, phi),
  (2) my reconstruction using UNIFORM library sizes S_i/n instead of the sampled real ones.
"""
from __future__ import annotations
import sys, time
import numpy as np
from scipy.sparse import csr_matrix
ROOT = "~/Code/virtual-cell-challenge-2026"
sys.path.insert(0, f"{ROOT}/src"); sys.path.insert(0, f"{ROOT}/experiments/096-mse-floor")
from cell_eval2.prep import bulk_lognorm_means
from decompose import terms
from msekit import CTRL, TS, group_stats
from vcc2026.donor_transfer import emit_counts
import anndata as ad

which = sys.argv[1] if len(sys.argv) > 1 else "orion"
REAL, PRED = {
    "h1": (f"{ROOT}/data/harness/h1-test.h5ad",
           f"{ROOT}/experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
    "orion": (f"{ROOT}/data/harness/orion-hct116-thinned.h5ad",
              f"{ROOT}/experiments/095-clean-harness/out/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad"),
}[which]
rp, rs, rc, rj, rg = group_stats(REAL)
pp, ps, pc, pj, pg = group_stats(PRED)
G = rg.size
ci = int(np.flatnonzero(rp == CTRL)[0])
bctrl = bulk_lognorm_means(rs[ci][None, :].astype(np.float64), TS)[0]
dh0 = bulk_lognorm_means(ps.astype(np.float64), TS) - bctrl[None, :]
S = ps.sum(1).astype(np.float64)
tb = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
print(f"# {which}: cached prediction X={tb['X']:.6f} A={tb['A']:.6f} "
      f"cells/group median={np.median(pc):.0f}")

# the real control library sizes, which is what 031-somatic's build() actually sampled
R = ad.read_h5ad(REAL, backed="r")
rl = np.asarray(R.obs["target_gene"].astype(str).values)
Xr = R.X[np.flatnonzero(rl == CTRL)]
libs_real = np.asarray(Xr.sum(1)).ravel().astype(np.float64)
R.file.close()
print(f"# real control libs: median {np.median(libs_real):.0f}  "
      f"mean {libs_real.mean():.0f}  sd {libs_real.std():.0f}")

NC = int(round(float(np.median(pc))))
for lib_mode in ("uniform", "sampled-real"):
    for spread, phi in ((0.25, 0.0), (0.5, 4.0)):
        t0 = time.time()
        sums = np.zeros_like(ps)
        rng0 = np.random.default_rng(5)
        for i in range(pp.size):
            prof = np.expm1(bctrl + dh0[i]); np.clip(prof, 0.0, None, out=prof)
            tot = prof.sum()
            if tot <= 0:
                continue
            prof = prof / tot
            libs = (np.full(NC, S[i] / NC) if lib_mode == "uniform"
                    else libs_real[rng0.choice(libs_real.size, NC,
                                               replace=libs_real.size < NC)] * (S[i] / (NC * libs_real.mean())))
            v, col, ip = emit_counts(prof, libs, 0.0, np.random.default_rng(7919 + i),
                                     mode="quasi", spread=spread, overdispersion=phi)
            M = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                           shape=(NC, G))
            sums[i] = np.asarray(M.sum(0)).ravel().astype(np.int64)
        t = terms(pp, sums, np.full(pp.size, float(NC)), pj, rp, rs, rc, rj, rg)
        # how exactly did the per-gene group totals survive?
        rel = np.abs(sums[:5].astype(float) - ps[:5].astype(float)).sum() / max(ps[:5].sum(), 1)
        print(f"  libs={lib_mode:<13} spread={spread:<5} phi={phi:<4} "
              f"X={t['X']:.6f} ({100*t['X']/tb['X']:.1f} % of cached)  A={t['A']:.6f} "
              f"({100*t['A']/tb['A']:.1f} %)  per-gene total drift={rel:.4f}  [{time.time()-t0:.0f}s]",
              flush=True)
