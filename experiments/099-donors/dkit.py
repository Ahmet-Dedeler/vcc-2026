"""099 machinery: measure a DONOR SET's directional accuracy on a harness.

Adapted from `experiments/097-direction/wkit.py` (which is itself an exact reconstruction of
`031-somatic/run.py::build` at the shipped grid point) with one change: the donor set is a
parameter, donors may live outside `data/donors/`, and the frozen artifacts (whitener, wcols,
emitter `L`) are pinned to the BASE two-donor configuration so that the ONLY thing an arm
varies is the donor set entering `MultiDonorModel.effect()`.

Frozen across every arm, and asserted:
  * the whitening operator, loaded from `095-clean-harness`'s cache, never refitted;
  * `wcols`, re-derived from the BASE donors and sha256-checked against the cached operator;
  * the emitter's per-perturbation `L`, loaded from 097's cache;
  * the amplitude policy: `E *= ALPHA * effect_scale / median(||E||)`, so every arm emits the
    same median effect norm. `W = X^2/A` is scale-invariant anyway; `nmae` is not, which is why
    the amplitude POLICY rather than the amplitude VALUE is what is held fixed.

LEGALITY. `orion-hct116-thinned` is a thinning of the `orion_hct116` donor and `h1-test` is
`h1_2025`'s own panel, so a configuration containing either may not be scored on the
corresponding harness. `Ctx` refuses to construct such a pair.
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
sys.path.insert(0, str(ROOT / "experiments" / "098-nmae"))

from vcc2026.donor_transfer import TS, Donor            # noqa: E402
from vcc2026.multi_donor import MultiDonorModel         # noqa: E402
from vcc2026.whiten import Whitener                     # noqa: E402
from decompose import terms                             # noqa: E402
from msekit import group_stats                          # noqa: E402
from cell_eval2.metrics.delta import _across_pert_budget  # noqa: E402

CONTROL, PERT_COL = "non-targeting", "target_gene"
GAMMA, Q, ALPHA, SPREAD, MIN_CPM = 0.65, 1.0, 1.08, 0.25, 5.0
N_EMIT = 400
BASE_DONORS = ("k562_gwps", "orion_hek293t")
# The two donors a harness is ALLOWED to predict from, and the pair its whitening operator was
# fitted on. Both are the production donor set minus the one the harness's reference IS.
BASE_BY_HARNESS = {
    "h1-test": ("k562_gwps", "orion_hek293t"),
    "orion":   ("k562_gwps", "orion_hek293t"),
    "hek293t": ("k562_gwps", "orion_hct116"),
}

DONOR_PATH = {
    "k562_gwps":      ROOT / "data/donors/k562_gwps",
    "orion_hct116":   ROOT / "data/donors/orion_hct116",
    "orion_hek293t":  ROOT / "data/donors/orion_hek293t",
    "h1_2025":        ROOT / "data/donors/h1_2025",
    "k562_essential": ROOT / "data/donors/k562_essential",
    "nadig_jurkat":   ROOT / "data/donors/nadig_jurkat",
    "rpe1":           ROOT / "data/donors/rpe1",
    "kolf2":          ROOT / "data/external/kolf2",
    "perturbmatch":   ROOT / "experiments/099-donors/work/perturbmatch",
    "feng_ipsc":      ROOT / "experiments/099-donors/work/feng_ipsc",
    "cd4_rest":       ROOT / "data/external/cd4t/cloud-compact/cd4_compact/cd4_rest",
    "cd4_stim8hr":    ROOT / "data/external/cd4t/cloud-compact/cd4_compact/cd4_stim8hr",
    "cd4_stim48hr":   ROOT / "data/external/cd4t/cloud-compact/cd4_compact/cd4_stim48hr",
}

HARNESS = {
    "h1-test": dict(
        real=ROOT / "data/harness/h1-test.h5ad",
        pred=ROOT / "experiments/095-clean-harness/out-h1test/runs/"
                    "g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        whit=ROOT / "experiments/095-clean-harness/out-h1test/whitener_union.npz",
        lvec=ROOT / "experiments/097-direction/work/h1-test-Lvec.npy",
        refde=ROOT / "experiments/098-nmae/work/refde_h1.npz",
        banned=("h1_2025",),
    ),
    "orion": dict(
        real=ROOT / "data/harness/orion-hct116-thinned.h5ad",
        pred=ROOT / "experiments/095-clean-harness/out/runs/"
                    "g0.65_q1.0_a1.08_spread0.25/pred.h5ad",
        whit=ROOT / "experiments/095-clean-harness/out/whitener_union.npz",
        lvec=ROOT / "experiments/097-direction/work/orion-Lvec.npy",
        refde=ROOT / "experiments/098-nmae/work/refde_orion.npz",
        banned=("orion_hct116",),
    ),
    # NEW SURFACE, built by this experiment. `scripts/build_orion_harness.py --line HEK293T`
    # with the same greedy+thin policy and the same 132 x 150 shape as the HCT116 harness,
    # PASSing `check_harness_depth.py` at ratio 0.9999. There is no cached whitening operator
    # for it and none is fitted (that would be a new, unvalidated operator), so this surface
    # is used at `gamma = 0` only -- which is the setting at which the donor-set question is
    # asked without the frozen-operator confound anyway. `pred` is the reference itself: it is
    # read ONLY for `keep`/`cols`/`w`/`sumD`/`Kc`, all of which are reference-side; `claim`
    # and `dh_ship` are therefore meaningless here and nothing in 099 reads them on it.
    "hek293t": dict(
        real=ROOT / "data/harness/orion-hek293t-thinned.h5ad",
        pred=ROOT / "data/harness/orion-hek293t-thinned.h5ad",
        whit=ROOT / "experiments/099-donors/work/whit/whitener_hek293t.npz",
        lvec=None,
        refde=ROOT / "experiments/099-donors/work/refde_hek293t.npz",
        banned=("orion_hek293t",),
    ),
}
MIN_GATE = 10          # 098's minimum reference DE-gate size for a scorable row


class Ctx:
    """One harness, loaded once. Donor sets are swapped in via `set_donors`."""

    def __init__(self, name: str):
        h = HARNESS[name]
        self.name = name
        self.banned = set(h["banned"])

        a = ad.read_h5ad(h["real"])
        X = a.X.tocsr() if not sp.isspmatrix_csr(a.X) else a.X
        labels = a.obs[PERT_COL].astype(str).to_numpy()
        self.genes = np.asarray(a.var_names)
        Xc = X[labels == CONTROL].tocsr()
        Pc = np.asarray(Xc.sum(0)).ravel().astype(np.float64)
        self.b_ctrl = np.log1p(TS * Pc / Pc.sum())
        self.libs = np.asarray(Xc.sum(1)).ravel()
        self.perts = sorted(set(labels) - {CONTROL})
        del a, X, Xc

        # -- frozen whitener + wcols, derived from the BASE donors exactly as 097 does ------
        self.base_donors = BASE_BY_HARNESS[name]
        base = {n: Donor.load(DONOR_PATH[n]) for n in self.base_donors}
        bm = MultiDonorModel(base, self.genes, gene_norm=False, eb=False)
        shared: set[str] = set()
        for n in self.base_donors:
            shared |= set(bm.target_genes[bm.tcols[n]])
        chal = set(pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str))
        rt = set(self.perts)
        shared -= (chal if rt <= chal else rt)
        ti = {g: i for i, g in enumerate(bm.target_genes)}
        self.wcols = np.array(sorted(ti[g] for g in shared), dtype=np.int64)
        if h["whit"] is None:
            self.W = None
        else:
            meta = json.loads(Path(str(h["whit"]).replace(".npz", ".json")).read_text())
            if "wcols_sha256" in meta:
                got = hashlib.sha256(self.wcols.tobytes()).hexdigest()
                if got != meta["wcols_sha256"] or len(self.wcols) != meta["n_cols"]:
                    raise SystemExit(f"wcols MISMATCH on {name} -- STOP")
            self.W = Whitener.load(h["whit"])
            if self.W.evecs.shape[0] != len(self.wcols):
                raise SystemExit(f"operator/wcols shape mismatch on {name}: "
                                 f"{self.W.evecs.shape} vs {len(self.wcols)}")
        del base, bm

        # -- the scored reference, through the scorer's own primitives ---------------------
        rp, rs, rc, rj, rg = group_stats(str(h["real"]))
        pp, ps, pc, pj, pg = group_stats(str(h["pred"]))
        assert np.array_equal(rg, pg) and np.array_equal(rg, self.genes)
        t = terms(pp, ps, pc, pj, rp, rs, rc, rj, rg)
        for k, c in t["cols"].items():
            t["dh"][k, c] = 0.0
            t["d"][k, c] = 0.0
        self.t = t
        self.targets = [str(x) for x in pp[t["keep"]]]
        assert self.targets == self.perts, "row order drift"
        self.w, self.d, self.dh_ship = t["w"], t["d"], t["dh"]
        self.sumD, self.claim, self.Kc, self.cols = t["sumD"], t["claim"], t["Kc"], t["cols"]
        self.P, self.G = self.d.shape
        if h["lvec"] is None:
            # no shipped prediction to recover `L` from: use the emitter's own rule on this
            # harness's control library pool, once, and freeze it across every arm.
            rng = np.random.default_rng(0)
            self.Lvec = np.array(
                [int(min(np.median(self.libs[rng.choice(len(self.libs), N_EMIT,
                                   replace=len(self.libs) < N_EMIT)]), 1_000_000))
                 for _ in self.perts], dtype=np.int64)
        else:
            self.Lvec = np.load(h["lvec"])
        assert self.Lvec.size == self.P

        # -- 098's reference DE tables, for the exact `nmae` --------------------------------
        z = np.load(h["refde"], allow_pickle=True)
        assert np.array_equal(z["genes"].astype(str), self.genes), "refde gene order"
        self.gate = z["gate"]
        self.ref_lin = z["ref_lin"]
        rperts = z["perts"].astype(str)
        self.rlfc = z["lfc"].astype(np.float64)
        si, sp_ = z["sig_idx"], z["sig_ptr"]
        sig = [si[sp_[i]:sp_[i + 1]] for i in range(len(rperts))]
        rpos = {p: i for i, p in enumerate(rperts)}
        rows, rrows = [], []
        for i, p in enumerate(self.targets):
            j = rpos.get(p)
            if j is None or sig[j].size < MIN_GATE:
                continue
            rows.append(i); rrows.append(j)
        self.rows = np.asarray(rows); self.rrows = np.asarray(rrows); self.sig = sig

        self.model = None
        self.donor_names: tuple[str, ...] = ()

    # ------------------------------------------------------------------ donor set selection
    def set_donors(self, names: tuple[str, ...]):
        bad = self.banned & set(names)
        if bad:
            raise ValueError(f"ILLEGAL on harness {self.name}: {sorted(bad)}")
        donors = {n: Donor.load(DONOR_PATH[n]) for n in names}
        self.model = MultiDonorModel(donors, self.genes, gene_norm=False, eb=False)
        self.donor_names = tuple(names)
        del donors
        return self.model

    def coverage(self) -> dict[str, int]:
        m = self.model
        return {n: int(sum(p in m.pert_index[n] for p in self.perts)) for n in self.donor_names}

    # ------------------------------------------------------------------ the pipeline itself
    def effects(self, wfun=None, gamma: float | None = None) -> np.ndarray:
        gamma = GAMMA if gamma is None else gamma
        if gamma > 0 and self.W is None:
            raise ValueError(f"harness {self.name} has no cached whitener; gamma must be 0")
        rows = [self.model.effect(p, 1.0, 0.0, weights=(wfun(p) if wfun else None))
                for p in self.perts]
        E = np.vstack(rows).astype(np.float64)
        if gamma > 0:
            before = float(np.median(np.linalg.norm(E[:, self.wcols], axis=1)))
            E[:, self.wcols] = self.W.apply(E[:, self.wcols], gamma, 0.0, renorm=False)
            E[:, self.wcols] *= before / (
                float(np.median(np.linalg.norm(E[:, self.wcols], axis=1))) or 1.0)
        E *= ALPHA * self.model.effect_scale / (
            float(np.median(np.linalg.norm(E, axis=1))) or 1.0)
        return E

    def dh_exact(self, E: np.ndarray) -> np.ndarray:
        P = np.vstack([self.model.profile_from_effect(E[i], self.b_ctrl, min_cpm=MIN_CPM,
                                                      gate_unexpressed=True)
                       for i in range(len(self.perts))])
        T = np.floor(N_EMIT * self.Lvec[:, None] * P)
        out = np.log1p(TS * T / T.sum(axis=1, keepdims=True)) - self.b_ctrl[None, :]
        for k, c in self.cols.items():
            out[k, c] = 0.0
        return out

    # ----------------------------------------------------------------------------- scoring
    def XAW(self, dh):
        X = float((np.einsum("ij,ij->i", dh, self.d) * self.w).sum())
        A = float((np.einsum("ij,ij->i", dh, dh) * self.w).sum())
        return X, A, (X * X / A if A > 0 else 0.0)

    def kappa(self, dh):
        X = float((np.einsum("ij,ij->i", dh, self.d) * self.w).sum())
        B = float(_across_pert_budget(dh, self.cols, self.w))
        return X / np.sqrt(B * self.sumD)

    def pds(self, dh):
        tgt = np.asarray([self.cols[k] for k in sorted(self.cols)], int)
        keep = np.ones(self.G, bool); keep[tgt] = False
        Pm, Rm = dh[:, keep], self.d[:, keep]
        pn = np.linalg.norm(Pm, axis=1); rn = np.linalg.norm(Rm, axis=1)
        ok = (pn > 0) & (rn > 0)
        C = (Pm @ Rm.T) / np.maximum(pn[:, None] * rn[None, :], 1e-300)
        C[~ok, :] = -1.0; C[:, ~ok] = -1.0
        n = C.shape[0]
        own = np.diag(C).copy()
        k = (C > own[:, None]).sum(1) + 0.5 * ((C == own[:, None]).sum(1) - 1)
        return float(np.mean(1.0 - k / (n - 1)))

    def lfc_hat(self, eff_rows):
        prof = np.expm1(self.b_ctrl[None, :] + eff_rows)
        np.clip(prof, 0.0, None, out=prof)
        prof /= prof.sum(1, keepdims=True)
        return np.log2((1e6 * prof[:, self.gate] + 1e-9) / (self.ref_lin[None, :] + 1e-9))

    def nmae(self, dh, per_row=False):
        """098's exact analytic `de_wilcoxon_lfc_nmae` on the reference's own gate."""
        L = self.lfc_hat(dh[self.rows])
        out = np.empty(self.rows.size)
        for k, j in enumerate(self.rrows):
            gs = self.sig[j]; r = self.rlfc[j][gs]
            out[k] = np.abs(L[k][gs] - r).sum() / np.abs(r).sum()
        return out if per_row else float(out.mean())

    def signdot(self, dh):
        """`<sign(lfc_ref), lfc_hat>` -- 098's direction leg of `nmae`, normalised the same way.

        `NMAE = 1 - signdot` to within +0.0008 / +0.0036 in our regime (LOG 2026-08-29 §2).
        """
        L = self.lfc_hat(dh[self.rows])
        out = np.empty(self.rows.size)
        for k, j in enumerate(self.rrows):
            gs = self.sig[j]; r = self.rlfc[j][gs]
            out[k] = float(np.sign(r) @ L[k][gs]) / np.abs(r).sum()
        return float(out.mean())

    def per_pert_cos(self, dh):
        """True per-perturbation cosine between prediction and reference, panel targets out."""
        tgt = np.asarray([self.cols[k] for k in sorted(self.cols)], int)
        keep = np.ones(self.G, bool); keep[tgt] = False
        Pm, Rm = dh[:, keep], self.d[:, keep]
        pn = np.linalg.norm(Pm, axis=1); rn = np.linalg.norm(Rm, axis=1)
        return np.einsum("ij,ij->i", Pm, Rm) / np.maximum(pn * rn, 1e-300)


# ------------------------------------------------------------------------------- the keys
def donor_matrices(model, perts):
    """`(M[name] -> [P, G] renormalised effect, H[name] -> [P] carried)` for the key builders."""
    M, H = {}, {}
    G = len(model.target_genes)
    for n in model.names:
        Mx = np.zeros((len(perts), G)); h = np.zeros(len(perts), bool)
        for r, p in enumerate(perts):
            i = model.pert_index[n].get(p)
            if i is None:
                continue
            Mx[r, model.tcols[n]] = model.s[n][i][model.dcols[n]]
            h[r] = True
        M[n] = model._renorm[n] * Mx
        H[n] = h
    return M, H


def unanimity_key(M, H, dons, P, G):
    """`MultiDonorModel.donor_agreement` generalised: geometric mean gated on sign unanimity."""
    out = np.zeros((P, G))
    for i in range(P):
        have = [n for n in dons if H[n][i]]
        if len(have) < 2:
            continue
        st = np.vstack([M[n][i] for n in have])
        sg = np.sign(st)
        un = (np.abs(sg.sum(0)) == len(have)) & np.all(st != 0, axis=0)
        out[i] = np.where(un, np.exp(np.log(np.abs(st) + 1e-300).mean(0)), 0.0)
    return out


def pairwise_key(M, H, dons, P, G):
    out = np.zeros((P, G))
    for i in range(P):
        have = [n for n in dons if H[n][i]]
        if len(have) < 2:
            continue
        tot = np.zeros(G); cnt = 0
        for a in range(len(have)):
            for b in range(a + 1, len(have)):
                x, y = M[have[a]][i], M[have[b]][i]
                tot += np.sqrt(np.abs(x) * np.abs(y)) * (
                    (np.sign(x) == np.sign(y)) & (x != 0) & (y != 0))
                cnt += 1
        out[i] = tot / max(cnt, 1)
    return out


def donor_cos_key(M, H, dons, P):
    """Mean pairwise cosine between carrying donors' effect vectors -- the per-pert key."""
    out = np.zeros(P); nd = np.zeros(P, int)
    for i in range(P):
        have = [n for n in dons if H[n][i]]
        nd[i] = len(have)
        if len(have) < 2:
            continue
        vals = []
        for a in range(len(have)):
            for b in range(a + 1, len(have)):
                x, y = M[have[a]][i], M[have[b]][i]
                nx, ny = np.linalg.norm(x), np.linalg.norm(y)
                vals.append(float(x @ y) / max(nx * ny, 1e-12))
        out[i] = float(np.mean(vals))
    return out, nd
