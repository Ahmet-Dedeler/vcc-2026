"""170 — the h1-test scoring surface, shared by every measurement in this experiment.

RULE 1 / 168d: h1-test is the sizing surface. It is the only surface validated against a
board `reach` point, it carries a 3.6% target-independent blob against the Orion pair's
26.4% / 44.1%, and it is not an Orion line, so the PRODUCTION `union_pruned` operator can be
applied to it without leakage. `data/donors/h1_2025` covers 85/85 of its targets and is pure
oracle leakage; it is NOT in the donor set here (k562_gwps + orion_hct116 + orion_hek293t
only, which is also production's exact three-donor set).

Everything here calls a shipped function or a `cell_eval2` primitive:
  render / gate / mass fix   `make_submission_pm_cd4.profile_from_effect` (MASS_FIX=minnorm)
  emission                   `vcc2026.donor_transfer.emit_counts` + the builder's
                             `emission_spread_vector` (a15's tier-8000 policy)
  mse                        096's `msekit.score` + `decompose.terms` -- cell_eval2's own
                             `mse_unbiased_capped` / `distance_unbiased` / jackknife
  pds                        `cell_eval2.metrics.discrimination_score` on the vcc2026 profile
  DE                         092's `de_kernel` (scanpy-identical rank-sum) + BH, then 153's
                             `debench.Surface` for nmae / reach / fid / jac
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"

os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

sys.path.insert(0, str(REPO / "src"))
sys.path.insert(0, str(REPO / "experiments/096-mse-floor"))
sys.path.insert(0, str(REPO / "experiments/092-call-set-truncation"))
sys.path.insert(0, str(REPO / "experiments/153-de-bench"))

H1 = str(REPO / "data/harness/h1-test.h5ad")
TS = 50_000.0
NCELL = 400
CTRL = "non-targeting"

# a15's shipped emission policy (`submissions/a15-cd4f010.vcc.meta.json`).
MIN_CPM = 10.0
SPREAD, TAIL_SPREAD, KEEP_TOP, KEEP_SPREAD = 0.25, 1.0, 8000, 0.25
OVERDISPERSION, DISPERSION = 4.0, 0.1


def builder():
    path = REPO / "scripts/make_submission_pm_cd4.py"
    spec = importlib.util.spec_from_file_location("prod_builder_170k", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_170k"] = mod
    spec.loader.exec_module(mod)
    if mod.MASS_FIX != "minnorm":
        raise RuntimeError("VCC_MASS_FIX must be minnorm (a15 ships it)")
    return mod


def _labels(h, path):
    """h5ad string column, categorical or plain."""
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        c = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                        for x in o["categories"][:]])
        return c[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


class Surface170:
    """h1-test's reference side, cached once: group sums, jackknife, control profile, DE."""

    def __init__(self):
        from msekit import group_stats

        self.rp, self.rs, self.rc, self.rj, self.rg = group_stats(H1)
        self.ci = int(np.flatnonzero(self.rp == CTRL)[0])
        from cell_eval2.prep import bulk_lognorm_means

        self.rm = bulk_lognorm_means(self.rs.astype(np.float64), TS)
        self.b_ctrl = self.rm[self.ci]
        self.genes = self.rg
        # production's own gate, on the surface's own control (builder line ~1290)
        expr = np.expm1(self.b_ctrl)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        self.gate = (cpm > MIN_CPM) & (self.b_ctrl > 0)
        self.targets = np.asarray([str(p) for p in self.rp if str(p) != CTRL])
        self.trow = np.asarray([int(np.flatnonzero(self.rp == t)[0]) for t in self.targets])
        # 400 real control library sizes, drawn once (quasi uses only their median, so this
        # is stable, but it is drawn from the surface's own controls exactly as the builder
        # draws from each context's controls).
        self.libs = self._control_libraries()

    def _control_libraries(self):
        import h5py
        import scipy.sparse as sp

        with h5py.File(H1, "r") as h:
            n, g = (int(x) for x in h["X"].attrs["shape"])
            o = h["obs/target_gene"]
            cats = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                               for x in o["categories"][:]])
            lab = cats[o["codes"][:]]
            ip = h["X/indptr"][:]
            rows = np.flatnonzero(lab == CTRL)
            data = h["X/data"]
            tot = np.empty(rows.size, dtype=np.float64)
            for k, r in enumerate(rows):
                tot[k] = float(np.asarray(data[int(ip[r]):int(ip[r + 1])]).sum())
        rng = np.random.default_rng(0)
        return tot[rng.choice(tot.size, NCELL, replace=tot.size < NCELL)]


def to_h1_axis(effect_prod, prod_genes, h1_genes):
    """Map an 18,533-gene production effect onto h1-test's 18,080-gene axis.

    3 of h1-test's genes are absent from the production axis and get exactly zero effect.
    They stay inside production's CPM gate, so the mass fix may place a correction on them;
    that is 3 genes in 18,080 and is reported rather than special-cased.
    """
    pos = {g: i for i, g in enumerate(map(str, prod_genes))}
    idx = np.asarray([pos.get(str(g), -1) for g in h1_genes], dtype=np.int64)
    out = np.zeros((effect_prod.shape[0], idx.size), dtype=np.float64)
    have = idx >= 0
    out[:, have] = effect_prod[:, idx[have]]
    return out, have


def deliver(B, effect_h1, b_ctrl):
    """Render exactly as the builder does: gate -> mass fix -> profile; return (profile, dh).

    `dh` is what the SCORER reads back, `log1p(TS*p) - b_ctrl`, i.e. the DELIVERED geometry
    (156): every effect-level screen before experiment 156 scored the intended object.
    """
    P = np.empty_like(effect_h1)
    for i in range(effect_h1.shape[0]):
        P[i] = B.profile_from_effect(effect_h1[i], b_ctrl, min_cpm=MIN_CPM)
    dh = np.log1p(TS * P) - b_ctrl[None, :]
    return P, dh


class DEKit:
    """h1-test's control block, prepared once: the scorer's DE gate + the sorted rank base.

    The universe is the scorer's own reference-side gate (mean control CPM > 5,
    `distances.py:98` / `de_compute.py:414`), so it is IDENTICAL for every arm in the sweep
    and cannot be moved by gamma.
    """

    def __init__(self, gate_cpm=5.0, chunk=2000):
        import h5py
        from de_kernel import ctrl_tiesums

        # Read off h5py as a float32/int32 CSR (2.7 GB), never through `anndata.read_h5ad`
        # into dense: the naive route holds a dense 18,400 x 18,080 float64 control block,
        # a dense log-CPM copy of it and a fancy-indexed transpose at once, and takes the
        # process to 13.6 GiB. h1-test's rows are also INTERLEAVED (85/85 targets are
        # non-contiguous), so block-slicing the file directly would re-read almost all of it
        # per chunk.
        import scipy.sparse as sp

        with h5py.File(H1, "r") as h:
            n, g = (int(x) for x in h["X"].attrs["shape"])
            o = h["obs/target_gene"]
            cats = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                               for x in o["categories"][:]])
            self._lab = _labels(h, "obs/target_gene")
            gfull = _labels(h, "var/_index")
            self._X = sp.csr_matrix((h["X/data"][:], h["X/indices"][:], h["X/indptr"][:]),
                                    shape=(n, g))
        self._G = g
        ci = np.flatnonzero(self._lab == CTRL)

        acc = np.zeros(g, dtype=np.float64)             # pass 1: mean control CPM per gene
        for lo in range(0, ci.size, chunk):
            X = self._dense(ci[lo:min(lo + chunk, ci.size)])
            s = X.sum(1)
            s[s == 0] = 1.0
            X *= (1e6 / s[:, None])
            acc += X.sum(0)
        self.ref_lin_full = acc / ci.size
        self.keep = self.ref_lin_full > gate_cpm
        self.genes = gfull[self.keep]
        self.ref_lin = self.ref_lin_full[self.keep]

        cs = np.empty((int(self.keep.sum()), ci.size), dtype=np.float64)  # pass 2
        for lo in range(0, ci.size, chunk):
            hi = min(lo + chunk, ci.size)
            X = self._dense(ci[lo:hi])
            s = X.sum(1)
            s[s == 0] = 1.0
            cs[:, lo:hi] = np.log1p(X[:, self.keep] * (1e6 / s[:, None])).T
        cs.sort(axis=1)
        self.cs = cs
        self.tie = ctrl_tiesums(self.cs)

    def _dense(self, rows):
        return np.asarray(self._X[rows].todense(), dtype=np.float64)

    def real_table(self, targets):
        """The reference's own DE table. Gamma-independent -> computed once, cached."""
        cache = WORK / "de_real_h1.npz"
        if cache.exists():
            z = np.load(cache, allow_pickle=True)
            if (np.array_equal(z["genes"].astype(str), self.genes)
                    and np.array_equal(z["perts"].astype(str), targets)):
                return {k: z[k] for k in ("real_p", "real_padj", "real_lfc")}
        P, Q, L = self._alloc(targets.size)
        for i, t in enumerate(targets):
            P[i], Q[i], L[i] = self.row(self._dense(np.flatnonzero(self._lab == str(t))))
        np.savez_compressed(cache, genes=self.genes, perts=targets,
                            real_p=P, real_padj=Q, real_lfc=L)
        return {"real_p": P, "real_padj": Q, "real_lfc": L}

    def release_reference(self):
        """Drop the 2.7 GB reference CSR once the real DE table is built."""
        self._X = None

    def _alloc(self, n):
        G = int(self.keep.sum())
        return (np.full((n, G), np.nan), np.full((n, G), np.nan), np.full((n, G), np.nan))

    def row(self, X):
        from de_kernel import wilcoxon_z
        from scipy.stats import false_discovery_control, norm

        s = X.sum(1)
        s[s == 0] = 1.0
        Bm = np.log1p(X * (1e6 / s[:, None]))
        z = wilcoxon_z(self.cs, self.tie, np.ascontiguousarray(Bm[:, self.keep].T))
        p = 2 * norm.sf(np.abs(z))
        lfc = np.log2((np.expm1(Bm[:, self.keep]).mean(0) + 1e-9) / (self.ref_lin + 1e-9))
        return p, false_discovery_control(p, method="bh"), lfc


def emit_and_score(B, profiles, effect_h1, b_ctrl, libs, de, seed=0):
    """a15's emission policy, one target at a time; group sums, jackknife and DE in one pass.

    Streaming matters: holding all 85 emitted blocks costs ~2.7 GB of CSR, and this session
    is capped at 10 GB with two sibling agents running.
    """
    from cell_eval2.moments import jackknife_correction
    from vcc2026.donor_transfer import emit_counts

    n, G = profiles.shape
    sums = np.zeros((n, G), dtype=np.int64)
    jk = np.zeros(n, dtype=np.float64)
    P, Q, L = de._alloc(n)
    for i in range(n):
        spread = B.emission_spread_vector(
            effect_h1[i], b_ctrl, min_cpm=MIN_CPM, tail_spread=TAIL_SPREAD,
            keep_top=KEEP_TOP, keep_spread=KEEP_SPREAD, spread_map="block")
        v, col, ip = emit_counts(
            profiles[i], libs, DISPERSION,
            np.random.default_rng(seed * 7919 + i), mode="quasi",
            overdispersion=OVERDISPERSION, spread=spread)
        X = csr_matrix((v.astype(np.float64), col.astype(np.int64), ip.astype(np.int64)),
                       shape=(NCELL, G))
        sums[i] = np.asarray(X.sum(0)).ravel().astype(np.int64)
        jk[i] = jackknife_correction(X, np.zeros(NCELL, np.intp), 1, TS, chunk=64)[0]
        P[i], Q[i], L[i] = de.row(np.asarray(X.todense()))
        del X
    return sums, jk, (P, Q, L)


def pds_exact(pred_perts, pred_sums, S):
    """`pds_cosine` through cell_eval2's own `discrimination_score`, vcc2026 settings."""
    from cell_eval2.metrics.discrimination import discrimination_score
    from cell_eval2.prep import bulk_lognorm_means

    pm = bulk_lognorm_means(pred_sums.astype(np.float64), TS)
    out = discrimination_score(
        pred_bulk=(pred_perts, pm), real_bulk=(S.rp, S.rm),
        pert_col="target_gene", control=CTRL, distance="cosine",
        rank_denominator="n-1", tie_policy="midrank", exclude_target_gene=True,
        exclusion_scope="panel", control_source="real", genes=S.genes,
        target_gene_map=None)
    return float(np.mean(list(out.values())))
