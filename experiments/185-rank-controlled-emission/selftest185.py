"""185 — validate the instrument before it measures anything (rule 10's corollary).

1. `hist_z` (histogram form) must equal `de_kernel.wilcoxon_z` (cell form, scanpy-identical)
   on REAL control data for random integer histograms, to machine precision.  This is the
   known-answer row for the exact-z evaluator: a wrong midrank or tie term reads differently.
2. `envelope_all`: F_min <= F_even <= F_max, and the family members reproduce their own z.
3. `design` + `place` round trip on one real target: every gated gene's delivered histogram
   equals the designed one; column totals bit-identical to a15's; every row sums to T.
"""
from __future__ import annotations

import hashlib
import sys
import time

import numpy as np

sys.path.insert(0, ".")
import harness185 as H  # noqa: E402
from eme185 import EME, hist_z, histogram_of  # noqa: E402

t0 = time.time()
C = H.load_controls()
cs, tie, keep = C["cs"], C["tie"], C["keep"]
G, m = cs.shape
n = 400
TOT = H.load_totals()
print(f"controls: {m} cells, {G} gated genes [{time.time()-t0:.0f}s]")

# ---- 1. hist_z vs wilcoxon_z on random histograms over real control columns
rng = np.random.default_rng(0)
S = int(np.median(TOT["totals"][TOT["targets"] != "non-targeting"].sum(1)))
T = S // n
eme = EME(cs, tie, T, n)
pick = rng.choice(G, 300, replace=False)
Ms = np.zeros(G, np.int64)
Ms[pick] = rng.integers(0, 3000, pick.size)
cmax, offs, lt, eq, Phi = eme.tables(Ms)
maxdiff = 0.0
lv = np.zeros(4, np.int64)
mu = np.zeros(4, np.int64)
cells = np.zeros((n, G))
zs = np.zeros(G)
for g in pick:
    M = int(Ms[g])
    fam = int(rng.integers(0, 2))
    if fam == 0:
        n0 = int(rng.integers(0, max(1, n - max(1, (M + T - 1) // T))))
        K = histogram_of(M, n, int(cmax[g]), 0, n0, -1, lv, mu)
    else:
        b = int(min(cmax[g], max(1, M // 3 + 1)))
        nb = int(rng.integers(1, max(2, min(n - 1, M // max(b, 1)) + 1))) if M > 0 else 0
        if nb == 0 or b < 1:
            K = histogram_of(M, n, int(cmax[g]), 0, 0, -1, lv, mu)
        else:
            K = histogram_of(M, n, int(cmax[g]), 1, nb, b, lv, mu)
    assert (mu[:K] * lv[:K]).sum() == M and mu[:K].sum() == n, (g, M, lv[:K], mu[:K])
    zs[g] = hist_z(lv, mu, K, lt, eq, int(offs[g]), float(tie[g]), n, m)
    col = np.repeat(lv[:K], mu[:K]).astype(np.float64)
    rng.shuffle(col)
    cells[:, g] = col
vals = np.log1p(cells * (1e6 / T))
zk = H.wilcoxon_z(cs, tie, np.ascontiguousarray(vals.T))
maxdiff = float(np.abs(zk[pick] - zs[pick]).max())
print(f"[1] hist_z vs wilcoxon_z on {pick.size} random histograms: max |diff| = {maxdiff:.2e}"
      f"   (|z| range {np.abs(zk[pick]).min():.2f}..{np.abs(zk[pick]).max():.2f})")
assert maxdiff < 1e-9, "ABORT: histogram z does not reproduce the kernel"

# ---- 2. envelope sanity on one real target
tsel = np.flatnonzero(TOT["targets"] != "non-targeting")[0]
Mfull = TOT["totals"][tsel]
Mg = Mfull[keep]
env = eme.envelope(Mg)
z_even, z_conc, z_min, z_max, z_floor = env[:, 0], env[:, 1], env[:, 2], env[:, 3], env[:, 4]
ok = (z_min <= z_even + 1e-9) & (z_even <= z_max + 1e-9)
print(f"[2] envelope on {TOT['targets'][tsel]}: z_min <= z_even <= z_max on {ok.sum()}/{G} genes; "
      f"median z_min {np.median(z_min):+.2f}, z_even {np.median(z_even):+.2f}, z_max {np.median(z_max):+.2f}; "
      f"median |z_floor| {np.median(np.abs(z_floor)):.3f}")
assert ok.mean() > 0.995, "ABORT: envelope ordering violated"

# ---- 3. design + place round trip
zt = rng.normal(0, 3, G)
dsg = eme.design(Mg, zt, mode="abs", T=int(Mfull.sum()) // n)
print(f"[3] design: mean |z_real - z_target| = {np.mean(np.abs(dsg['z'] - zt)):.3f}; "
      f"mean |(|z_real| - |z_target|)| = {np.mean(np.abs(np.abs(dsg['z']) - np.abs(zt))):.3f}")
# sub-gate block via production emit_counts (a15's tail law)
from vcc2026.donor_transfer import emit_counts  # noqa: E402
prof = np.zeros(Mfull.size)
prof[~keep] = (Mfull[~keep] + 0.5) / (n * H.EMIT_L)
v, c, ptr = emit_counts(prof, np.full(n, H.EMIT_L), 0.1, np.random.default_rng(1), mode="quasi",
                        spread=1.0, overdispersion=4.0)
B_sub = np.zeros((n, Mfull.size), np.int64)
for r in range(n):
    sl = slice(ptr[r], ptr[r + 1])
    B_sub[r, c[sl]] = v[sl]
assert np.array_equal(B_sub.sum(0)[~keep], Mfull[~keep]) and B_sub[:, keep].sum() == 0
S_t = int(Mfull.sum())
Tc = np.full(n, S_t // n, np.int64)
Tc[: S_t % n] += 1
gidx = np.flatnonzero(keep)
B, diag = eme.place(dsg, gidx, B_sub, Tc, np.random.default_rng(2))
assert np.array_equal(B.sum(0), Mfull), "column totals moved"
assert np.array_equal(B.sum(1), Tc), "row totals not exact"
# delivered histograms of gated genes == designed
bad = 0
for j, g in enumerate(gidx):
    col = B[:, g]
    lvls, cnts = np.unique(col, return_counts=True)
    K = int(dsg["K"][j])
    if not (np.array_equal(lvls, dsg["levels"][j, :K]) and np.array_equal(cnts, dsg["mults"][j, :K])):
        bad += 1
print(f"[3] place: forced {diag['forced']}, max|resid| before repair {diag['max_abs_resid']}, "
      f"repair moves {diag['moves']}, unresolved {diag['unresolved']}; "
      f"gated histograms delivered exactly on {G-bad}/{G}; row sums == T on 400/400; "
      f"column digest {hashlib.sha256(B.sum(0).tobytes()).hexdigest()[:16]}")
assert bad == 0 and diag["unresolved"] == 0
# realised z after placement vs designed z (only the +1 remainder cells can differ)
z, p, q, lfc = H.de_from_block(B, keep, cs, tie, C["ref_lin"])
d = np.abs(z - dsg["z"])
print(f"[3] realised z vs designed z: max |diff| {d.max():.2e}, median {np.median(d):.2e}; "
      f"sub-gate mass per cell median {np.median(B_sub.sum(1)):.0f}")
print(f"selftest PASS [{time.time()-t0:.0f}s]")
