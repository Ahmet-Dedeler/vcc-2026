"""Shared machinery for experiment 097: measure `W = X^2/A` for donor-combination arms.

The objective, restated from `experiments/096-mse-floor/MSE-MECHANISM.md` (b)+(c): once the
amplitude has been re-optimised onto the rho cliff, the only property of a model change that
survives is the scalar

        W = X^2 / A ,        gain in raw mse  =  2*sqrt(claim*W)/sumD

so a model change is worth exactly `sqrt(W_new/W_base)`.  `W` is invariant to a global rescale
of `dh` (X -> sX, A -> s^2 A), which is why the amplitude knob and the direction knob separate.

Everything here reconstructs `031-somatic/run.py`'s `main()`/`build()` EXACTLY for the shipped
config `g0.65 q1.0 a1.08 spread0.25 cpm5.0`, with two deliberate substitutions:

  * the whitener is LOADED from `095-clean-harness`'s cache (never refitted), and its `wcols`
    are re-derived here and checked bit-for-bit against the cached `wcols_sha256`;
  * emission is replaced by an EXACT reconstruction rather than a surrogate.  `emit_counts`
    in `mode="quasi"` sets `totals = floor(n*L*profile)` and then repairs every gene's total
    exactly, so no emission noise reaches the pseudobulk at all; the only stochastic input is
    `L = int(median(libs[pick]))`, one integer per perturbation.  `w_floor.py` recovers those
    integers from the shipped prediction's own group sums and reproduces the scorer's `dh`
    bit-for-bit (85/85 rows on h1-test; 131/132 on orion, the last off by one count in one
    gene from a floor tie).  Arms are then measured through the same path with `L` frozen.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("~/Code/virtual-cell-challenge-2026")
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments" / "096-mse-floor"))

from vcc2026.donor_transfer import TS, Donor  # noqa: E402
from vcc2026.multi_donor import MultiDonorModel  # noqa: E402
from vcc2026.whiten import Whitener  # noqa: E402

from decompose import terms  # noqa: E402
from msekit import group_stats  # noqa: E402

CONTROL, PERT_COL = "non-targeting", "target_gene"
# the shipped grid point, `experiments/095-clean-harness/grid.json`
GAMMA, Q, ALPHA, SPREAD, MIN_CPM = 0.65, 1.0, 1.08, 0.25, 5.0
N_EMIT = 400          # 031/run.py::build default
BASE_DONORS = ("k562_gwps", "orion_hek293t")

HARNESS = {
    "h1-test": dict(
        real=ROOT / "data/harness/h1-test.h5ad",
        pred=ROOT / "experiments/095-clean-harness/out-h1test/runs/"
                    "g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        whit=ROOT / "experiments/095-clean-harness/out-h1test/whitener_union.npz",
        banned="h1_2025",          # leak rule: h1-test IS h1_2025's panel
    ),
    "orion": dict(
        real=ROOT / "data/harness/orion-hct116-thinned.h5ad",
        pred=ROOT / "experiments/095-clean-harness/out/runs/"
                    "g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        whit=ROOT / "experiments/095-clean-harness/out/whitener_union.npz",
        banned="orion_hct116",     # leak rule: the harness is a thinning of this donor
    ),
}


# --------------------------------------------------------------------------- harness context
class Ctx:
    """Everything needed to turn a donor-weighting rule into an (X, A, W) triple."""

    def __init__(self, name: str, extra_donors: tuple[str, ...] = ()):
        h = HARNESS[name]
        self.name = name
        for d in extra_donors:
            if d == h["banned"]:
                raise ValueError(f"LEAK: {d} may not be used on harness {name}")
        self.donor_names = tuple(BASE_DONORS) + tuple(extra_donors)

        a = ad.read_h5ad(h["real"])
        X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
        labels = a.obs[PERT_COL].astype(str).to_numpy()
        self.genes = np.asarray(a.var_names)
        Xc = X[labels == CONTROL].tocsr()
        Pc = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
        self.b_ctrl = np.log1p(TS * Pc / Pc.sum())          # exactly build()'s b_ctrl
        self.libs = np.asarray(Xc.sum(1)).ravel()           # emitter's library-size pool
        self.perts = sorted(set(labels) - {CONTROL})
        del a, X, Xc

        donors = {n: Donor.load(ROOT / "data/donors" / n) for n in self.donor_names}
        self.model = MultiDonorModel(donors, self.genes, gene_norm=False, eb=False)
        del donors                                          # free the raw delta matrices

        # wcols exactly as 031's main() builds them for --whiten-axis union, --exclude-panel
        # auto -- but on the BASE two donors only, because that is what the cached whitener
        # was fitted with and we are forbidden from refitting.
        m = self.model
        shared: set[str] = set()
        for n in BASE_DONORS:
            shared |= set(m.target_genes[m.tcols[n]])
        chal = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")
                   ["target_gene"].astype(str))
        rt = set(self.perts)
        excluded = chal if rt <= chal else rt
        shared -= excluded
        ti = {g: i for i, g in enumerate(m.target_genes)}
        self.wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)

        meta = json.loads(Path(str(h["whit"]).replace(".npz", ".json")).read_text())
        got = hashlib.sha256(self.wcols.tobytes()).hexdigest()
        if got != meta["wcols_sha256"] or len(self.wcols) != meta["n_cols"]:
            raise SystemExit(
                f"wcols MISMATCH on {name}: rebuilt n={len(self.wcols)} sha={got[:16]} vs "
                f"cached n={meta['n_cols']} sha={meta['wcols_sha256'][:16]} -- STOP")
        self.wcols_sha = got
        self.W = Whitener.load(h["whit"])

        # the realized prediction, through the scorer's own primitives
        rp, rs, rc, rj, rg = group_stats(str(h["real"]))
        pp, ps, pc, pj, pg = group_stats(str(h["pred"]))
        assert np.array_equal(rg, pg) and np.array_equal(rg, self.genes)
        t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
        for k, c in t["cols"].items():
            t["dh"][k, c] = 0.0
            t["d"][k, c] = 0.0
        self.t = t
        self.targets = [str(x) for x in pp[t["keep"]]]
        assert self.targets == self.perts, "row order drift between build() and the scorer"
        self.w = t["w"]
        self.d = t["d"]
        self.dh_real = t["dh"]
        self.sumD = t["sumD"]
        self.claim = t["claim"]
        self.Kc = t["Kc"]
        # scorer's own control row, to confirm build()'s b_ctrl is the metric's origin
        self.bctrl_scorer = t["bctrl"]
        self.Lvec = None
        self.fit_L()

    # ------------------------------------------------------------------ the pipeline itself
    def effects(self, wfun) -> np.ndarray:
        """`E` for every perturbation, given a per-perturbation donor-weight rule.

        `wfun(pert) -> dict[str, float]` (or None for uniform).  Everything downstream --
        whitening, q, the alpha renormalisation -- is bit-identical to `031/run.py::build`.
        """
        rows = []
        for p in self.perts:
            w = wfun(p)
            rows.append(self.model.effect(p, 1.0, 0.0, weights=w))
        E = np.vstack(rows).astype(np.float64)
        if GAMMA > 0:
            before = float(np.median(np.linalg.norm(E[:, self.wcols], axis=1)))
            E[:, self.wcols] = self.W.apply(E[:, self.wcols], GAMMA, 0.0, renorm=False)
            E[:, self.wcols] *= before / (
                float(np.median(np.linalg.norm(E[:, self.wcols], axis=1))) or 1.0)
        if Q != 1.0:
            med = float(np.median(np.linalg.norm(E, axis=1))) or 1.0
            E = np.sign(E) * np.abs(E) ** Q
            E *= med / (float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
        E *= ALPHA * self.model.effect_scale / (
            float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
        return E

    def profiles(self, E: np.ndarray) -> np.ndarray:
        return np.vstack([self.model.profile_from_effect(E[i], self.b_ctrl, min_cpm=MIN_CPM,
                                                         gate_unexpressed=True)
                          for i in range(len(self.perts))])

    def dh_analytic(self, E: np.ndarray) -> np.ndarray:
        """`log1p(TS*profile) - b_ctrl`, the emission-free surrogate for the scored effect."""
        out = np.log1p(TS * self.profiles(E)) - self.b_ctrl[None, :]
        for k, c in self.t["cols"].items():
            out[k, c] = 0.0                      # the scorer never sees a row's own target
        return out

    def dh_exact(self, E: np.ndarray, Lvec: np.ndarray | None = None) -> np.ndarray:
        """The `dh` the scorer actually reads: `log1p(TS*floor(n*L*p)/Sigma) - b_ctrl`."""
        L = self.Lvec if Lvec is None else Lvec
        P = self.profiles(E)
        T = np.floor(N_EMIT * L[:, None] * P)
        out = np.log1p(TS * T / T.sum(axis=1, keepdims=True)) - self.b_ctrl[None, :]
        for k, c in self.t["cols"].items():
            out[k, c] = 0.0
        return out

    def XAW(self, dh: np.ndarray) -> tuple[float, float, float]:
        w = self.w
        X = float((np.einsum("ij,ij->i", dh, self.d) * w).sum())
        A = float((np.einsum("ij,ij->i", dh, dh) * w).sum())
        return X, A, (X * X / A if A > 0 else 0.0)

    def run(self, wfun, Lvec: np.ndarray | None = None) -> dict:
        E = self.effects(wfun)
        dh = self.dh_exact(E, Lvec)
        X, A, W = self.XAW(dh)
        return dict(X=X, A=A, W=W, E=E, dh=dh,
                    Enorm=float(np.median(np.linalg.norm(E, axis=1))),
                    dhnorm=float(np.median(np.linalg.norm(dh, axis=1))))

    def resample_L(self, seed: int) -> np.ndarray:
        """An alternative draw of the emitter's `L`, the only stochastic input to `dh`."""
        rng = np.random.default_rng(seed)
        return np.array([int(min(np.median(self.libs[rng.choice(len(self.libs), N_EMIT,
                                 replace=len(self.libs) < N_EMIT)]), 1_000_000))
                         for _ in self.perts], dtype=np.int64)

    # ------------------------------------------------------------ the emitter's `L` per row
    def fit_L(self) -> np.ndarray:
        """Recover `L = int(median(libs[pick]))` per perturbation from the shipped prediction.

        Cached: it depends only on the base config, and every arm holds it fixed (an arm
        changes the profile, not the control-cell library sizes the emitter samples).
        """
        cache = ROOT / f"experiments/097-direction/work/{self.name}-Lvec.npy"
        if cache.exists():
            self.Lvec = np.load(cache)
            return self.Lvec
        pp, ps, _, _, _ = group_stats(str(HARNESS[self.name]["pred"]))
        tot = {str(p): int(ps[i].sum()) for i, p in enumerate(pp)}
        P = self.profiles(self.effects(uniform))
        L = np.zeros(len(self.perts), dtype=np.int64)
        for i, p in enumerate(self.perts):
            lo, target = int(tot[p] / N_EMIT), tot[p]
            hit = -1
            for step in range(0, 3000):
                for cand in (lo + step, lo - step):
                    if cand > 0 and int(np.floor(N_EMIT * cand * P[i]).sum()) == target:
                        hit = cand
                        break
                if hit > 0:
                    break
            if hit < 0:      # a floor tie: take the nearest total instead of an exact match
                cands = np.arange(lo - 400, lo + 400)
                sums = np.array([int(np.floor(N_EMIT * cc * P[i]).sum()) for cc in cands])
                hit = int(cands[np.argmin(np.abs(sums - target))])
            L[i] = hit
        np.save(cache, L)
        self.Lvec = L
        return L


# ------------------------------------------------------------------------------- weight rules
def uniform(_p):
    return None


def fixed(weights: dict[str, float]):
    """A constant global weighting; donors absent for a perturbation drop out and renormalise."""
    def f(_p):
        return dict(weights)
    return f


def carriers(model, p) -> list[str]:
    return [n for n in model.names if p in model.pert_index[n]]


def donor_vec(model, n: str, p: str) -> np.ndarray:
    """The donor's renormalised gene-specific effect for `p`, on the target axis."""
    v = np.zeros(len(model.target_genes))
    v[model.tcols[n]] = model.s[n][model.pert_index[n][p]][model.dcols[n]]
    return model._renorm[n] * v
