"""518 — exact PLOLO arms on the dumped compositions (dump518-<fold>.npz), no operator (~1 GB).

An arm = a function of the dumped fold: per-target key multipliers on the unkeyed stack composition (key_scale restores the
non-panel Frobenius norm, as the builder's apply_pert_scale), the sharpening q, the bulk render gate (MIN_CPM), optional
per-gene multipliers, optional across-target centring. Rows are rendered and realised by arms540's route (rows518), scored by
comp420.score_rows (cos_agg, literal pds, board-emulated pds at lam 0.5).
KNOWN ANSWER (every run, abort unless |d| < 1e-12): the base arm V1_0.5_3 (and V1_0.5_6 when listed) reproduces
arms540-<fold>.json's realised cos_agg / pds / pds_lam0.5.

    from score518 import Fold;  F = Fold("k562_gwps");  F.score(w=..., q=1.75, min_cpm=10.0)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
sys.path.insert(0, str(HERE))
import json  # noqa: E402

import numpy as np  # noqa: E402

import lib420  # noqa: E402,F401  (repoints plolo.WORK to 420's read-only caches)
from lib420 import TS, frames, key_scale, profile_from_effect, realise_k, sharpen, weights  # noqa: E402
from comp420 import score_rows  # noqa: E402

WORK = HERE / "work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
TOL = 1e-12


def gate_at(b, min_cpm):
    expr = np.expm1(b)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    return (cpm > min_cpm) & (b > 0)


class Fold:
    def __init__(self, held: str, check: bool = True):
        self.held = held
        z = np.load(WORK / f"dump518-{held}.npz")
        self.targets = [str(t) for t in z["targets"]]
        self.fPP, self.fST = z["fPP"], z["fST"]
        self.kv1, self.carried = z["kv1"], z["carried"]
        self.scored = z["scored"]
        self.pp_pos = z["panel_positions"]
        self.truth_full = np.asarray(z["truth"], np.float64)
        self.truth = self.truth_full[:, self.scored]
        self.betas = {c: z[f"beta_{c}"] for c in "ABC"}
        self.nonpanel = np.ones(self.fST.shape[1], bool)
        self.nonpanel[self.pp_pos] = False
        self.fr = dict(frames())
        self.P_DE = {c: [profile_from_effect(self.fPP[i], b, g, "minnorm") for i in range(self.fPP.shape[0])]
                     for c, (b, g) in self.fr.items()}
        self.n = len(self.targets)
        self.w_base = weights(self.kv1, self.carried, 0.5, 3.0)[0]
        self.ref540 = json.loads((W540 / f"arms540-{held}.json").read_text())["arms"]
        if check:
            for tag, (lo, hi) in (("V1_0.5_3", (0.5, 3.0)), ("V1_0.5_6", (0.5, 6.0))):
                s, _ = self.score(w=weights(self.kv1, self.carried, lo, hi)[0])
                want = self.ref540[tag]["realised"]["16"]
                d = max(abs(s[k] - want[k]) for k in ("cos_agg", "pds", "pds_lam0.5"))
                assert d < TOL, f"KA FAILED {held} {tag}: {d}"
            self.ka = "PASS"

    def rows(self, w=None, q=1.75, min_cpm=None, gene_mult=None, centre=None, bulk=None, scale=None):
        """Realised delivered rows per frame. `bulk` overrides the keyed composition. centre: None | 'effect' (subtract the
        across-target mean of the keyed effect before sharpening, non-panel) | ('shrink', f) (subtract f x that mean)."""
        if bulk is None:
            bulk = key_scale(self.fST, self.w_base if w is None else w)[0]
        if centre is not None:
            f = 1.0 if centre == "effect" else float(centre[1])
            mu = bulk.mean(0, keepdims=True)
            bulk = bulk - f * mu
            bulk[:, self.pp_pos] = 0.0
        out, K = {}, {}
        for c, (b, g10) in self.fr.items():
            g = g10 if min_cpm is None else gate_at(b, min_cpm)
            F = sharpen(bulk, q, g, mode="row") if q is not None else bulk.copy()
            if scale is not None:
                F = F * float(scale[c])
            if gene_mult is not None:
                F = F * gene_mult[c][None, :]
            rea = np.empty_like(F)
            Kc = np.empty(self.n, int)
            for i in range(self.n):
                p_bk = profile_from_effect(F[i], b, g, "minnorm")
                p_re, _, Kc[i], _, _ = realise_k(self.P_DE[c][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            out[c], K[c] = rea, Kc
        return out, K

    def score(self, **kw):
        rea, K = self.rows(**kw)
        s, arr = score_rows(rea, self.scored, self.truth, self.betas, self.nonpanel)
        s["K_mean"] = float(np.mean(np.concatenate([K[c] for c in "ABC"])))
        return s, arr

    def energies(self, rea):
        """render-space A (mean row energy) per frame over non-panel genes (for energy-ratio sanity only)."""
        return {c: float((rea[c][:, self.nonpanel] ** 2).sum(1).mean()) for c in "ABC"}
