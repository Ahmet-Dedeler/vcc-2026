"""Combine several donor screens into one effect estimate.

Motivation, from two independent directions:

**Noise.** Nadig et al. 2025 measured that within K562, the raw log-fold-change correlation
between two *replicate* experiments is 0.16, while the noise-corrected true correlation is
0.90. Donor deltas are dominated by sampling noise, and averaging independent screens is the
cheapest way to remove it.

**Biology.** The same paper measured cross-line effect correlations of 0.40 (K562/RPE1) to
0.75 (HepG2/RPE1), with 44% of perturbations behaving cell-type-preferentially. A single
donor is therefore also a single point of *bias*, not just variance -- whatever is specific
to K562 is baked into every prediction. Averaging donors from different tissues cancels
line-specific components and keeps what is shared, which is exactly what transfers.

TransPert (3rd of 338 at VCC 2025) used similarity-weighted aggregation across ten reference
screens and posted the best DE score on the board, so this is the empirically validated shape.

Weighting schemes:
  uniform     equal weight to every donor carrying the perturbation
  basal       weight by similarity between the donor's control profile and the target's,
              on the shared gene axis -- the only context signal the Challenge gives us
  precision   weight by donor cell count for that perturbation (inverse-variance flavoured)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from .donor_transfer import Donor, TS, _norm_scale, eb_shrink


class MultiDonorModel:
    """Average gene-specific effects across donors, on the target's gene axis."""

    def __init__(self, donors: dict[str, Donor], target_genes: np.ndarray,
                 gene_norm: bool | dict[str, bool] = False, fallback: str = "mean",
                 eb: bool | dict[str, bool] = True,
                 cells_per_pert: dict[str, int] | None = None,
                 effect_scale: float = 6.817):
        """`gene_norm` may be per-donor. Whether it helps depends on how dense the donor's
        gene axis is: K562-gwps carries 8,248 genes all above 5 CPM and normalising gains
        ~1%, while Orion carries 38,606 including 3,372 unexpressed, where normalising
        costs ~32% even with a floor. Decide per donor, not globally."""
        self.names = list(donors)
        self.target_genes = np.asarray(target_genes)
        self.fallback = fallback
        G = len(self.target_genes)
        tindex = {g: i for i, g in enumerate(self.target_genes)}

        self.tcols: dict[str, np.ndarray] = {}
        self.dcols: dict[str, np.ndarray] = {}
        self.s: dict[str, np.ndarray] = {}
        self.pert_index: dict[str, dict[str, int]] = {}
        self.b_ctrl: dict[str, np.ndarray] = {}
        self.m_full = np.zeros(G, dtype=np.float32)

        gn = gene_norm if isinstance(gene_norm, dict) else {n: gene_norm for n in donors}
        ebd = eb if isinstance(eb, dict) else {n: eb for n in donors}
        cpp = cells_per_pert or {}
        self.gene_norm, self.eb = gn, ebd
        for name, d in donors.items():
            m, s = d.decompose()
            if ebd.get(name, False):
                s = eb_shrink(s, d.b_ctrl, cpp.get(name, 150))
            if gn.get(name, False):
                scale = _norm_scale(s)
                s_n = s / scale
                before = np.median(np.linalg.norm(s, axis=1))
                after = np.median(np.linalg.norm(s_n, axis=1)) or 1.0
                s = (s_n * (before / after)).astype(np.float32)

            cols, dcols = [], []
            for j, g in enumerate(d.genes):
                t = tindex.get(g)
                if t is not None:
                    cols.append(t); dcols.append(j)
            self.tcols[name] = np.asarray(cols, dtype=np.int64)
            self.dcols[name] = np.asarray(dcols, dtype=np.int64)
            self.s[name] = s
            self.pert_index[name] = {p: i for i, p in enumerate(d.perts)}

            b = np.zeros(G, dtype=np.float32)
            b[self.tcols[name]] = d.b_ctrl[self.dcols[name]]
            self.b_ctrl[name] = b

            mf = np.zeros(G, dtype=np.float32)
            mf[self.tcols[name]] = m[self.dcols[name]]
            self.m_full += mf / len(donors)

        # Two levels of normalisation, both so that `alpha` means one fixed thing.
        #
        # Per donor: equalise effect magnitude across donors, so a perturbation carried by a
        # large-effect donor does not silently get a bigger prediction than one carried by a
        # small-effect donor.
        self._scale = {n: float(np.median(np.linalg.norm(self.s[n][:, self.dcols[n]], axis=1)))
                       for n in self.names}
        ref = float(np.median(list(self._scale.values())))
        self._renorm = {n: ref / (self._scale[n] or 1.0) for n in self.names}
        self._m_norm = float(np.linalg.norm(self.m_full)) or 1.0

        # Global: pin the *combined* effect to a fixed reference magnitude. Without this,
        # `alpha` silently changes meaning whenever the pipeline changes -- EB shrinkage alone
        # drops the median effect norm from 5.9 to 2.0, so an alpha tuned before it would emit
        # effects 3x too small afterwards, moving nmae and fid for a reason that has nothing
        # to do with the model. `effect_scale` is that reference: alpha = 1 emits a typical
        # real effect magnitude, whatever the upstream processing.
        self.effect_scale = effect_scale
        probe = min(len(self.names), 1)
        n0 = self.names[0]
        raw_med = float(np.median(np.linalg.norm(self.s[n0][:, self.dcols[n0]], axis=1))) or 1.0
        self._global = effect_scale / (raw_med * (ref / (self._scale[n0] or 1.0)))
        self._median_s_norm = effect_scale

    def coverage(self, perts) -> dict[str, float]:
        out = {n: float(np.mean([p in self.pert_index[n] for p in perts])) for n in self.names}
        out["any"] = float(np.mean([any(p in self.pert_index[n] for n in self.names)
                                    for p in perts]))
        return out

    def basal_weights(self, b_ctrl_target: np.ndarray) -> dict[str, float]:
        """Cosine similarity between each donor's control profile and the target's.

        Computed on the genes the donor shares with the target, then softmaxed. This is the
        only signal the Challenge gives us about which line a context resembles.
        """
        sims = {}
        for n in self.names:
            c = self.tcols[n]
            a, b = self.b_ctrl[n][c], b_ctrl_target[c]
            na, nb = np.linalg.norm(a), np.linalg.norm(b)
            sims[n] = float(a @ b / (na * nb)) if na and nb else 0.0
        return sims

    def effect(self, pert: str, alpha: float, beta: float,
               weights: dict[str, float] | None = None,
               scheme: str = "uniform",
               b_ctrl_target: np.ndarray | None = None) -> np.ndarray:
        """Weighted mean of the donors carrying `pert`, on the target gene axis."""
        eff = beta * self.m_full
        have = [n for n in self.names if pert in self.pert_index[n]]

        if not have:
            if self.fallback == "mean" and alpha != 0.0:
                eff = eff + (alpha * self._median_s_norm / self._m_norm) * self.m_full
            return eff
        if alpha == 0.0:
            return eff

        if weights is not None:
            w = {n: weights.get(n, 0.0) for n in have}
        elif scheme == "basal" and b_ctrl_target is not None:
            sims = self.basal_weights(b_ctrl_target)
            # Sharpened softmax: cross-line similarities sit in a narrow band, so a plain
            # softmax over them is nearly uniform and buys nothing.
            ex = {n: np.exp(10.0 * sims[n]) for n in have}
            w = ex
        else:
            w = {n: 1.0 for n in have}
        tot = sum(w.values()) or 1.0

        acc = np.zeros(len(self.target_genes), dtype=np.float64)
        for n in have:
            v = np.zeros_like(acc)
            v[self.tcols[n]] = self.s[n][self.pert_index[n][pert]][self.dcols[n]]
            acc += (w[n] / tot) * self._renorm[n] * v
        return eff + (alpha * self._global) * acc.astype(np.float32)

    def donor_agreement(self, pert: str) -> np.ndarray:
        """Per-gene cross-donor confidence: geometric magnitude, gated on sign unanimity.

        Experiment 013 nominated `sqrt(|A|.|B|) x sign-agreement` as a confidence ordering
        and measured it at +24% reach over magnitude ordering on a LODO proxy. Experiment
        095 tested it against the real metric on a depth-repaired harness and it is the best
        buildable ordering key found: **+0.00678 overall** as a full ordering of the
        reference's DE pool, against +0.00448 for `|predicted lfc|`.

        The point is that it answers a different question from `|effect|`. Magnitude asks
        "how big do we think this is"; agreement asks "do independent cell lines tell the
        same story", which is what `reach`'s 0.9 purity bar actually rewards. A gene both
        donors move the same way is one we are far more likely to have the SIGN right on,
        and sign is the only thing purity counts.

        Generalises exp013's two-donor form to however many donors carry `pert`: the
        geometric mean of |delta| over the carrying donors, zeroed unless every one of them
        agrees on the sign. Donors are renormalised first, exactly as `effect()` does, so a
        large-effect donor cannot dominate the ranking.

        Returns a length-`len(target_genes)` score on the target axis, all-zero when fewer
        than two donors carry the perturbation (nothing to agree about).
        """
        have = [n for n in self.names if pert in self.pert_index[n]]
        if len(have) < 2:
            return np.zeros(len(self.target_genes), dtype=np.float64)
        mats = []
        for n in have:
            v = np.zeros(len(self.target_genes), dtype=np.float64)
            v[self.tcols[n]] = self.s[n][self.pert_index[n][pert]][self.dcols[n]]
            mats.append(self._renorm[n] * v)
        stack = np.vstack(mats)
        sign = np.sign(stack)
        unanimous = (np.abs(sign.sum(0)) == len(have)) & np.all(stack != 0, axis=0)
        geo = np.exp(np.log(np.abs(stack) + 1e-300).mean(0))
        return np.where(unanimous, geo, 0.0)

    def cpm_gate_mask(self, b_ctrl_target: np.ndarray, min_cpm: float) -> np.ndarray:
        if min_cpm <= 0:
            return np.ones(len(b_ctrl_target), dtype=bool)
        e = np.expm1(b_ctrl_target)
        return (1e6 * e / e.sum()) > min_cpm

    def profile_from_effect(self, eff: np.ndarray, b_ctrl_target: np.ndarray,
                            min_cpm: float = 20.0,
                            gate_unexpressed: bool = True) -> np.ndarray:
        """Gate and exponentiate an effect vector that was computed elsewhere.

        Split out from `profile` because spectral whitening (experiment 017) transforms the
        whole 300-perturbation effect matrix jointly, so the effect can no longer be derived
        one perturbation at a time.
        """
        eff = np.asarray(eff, dtype=np.float64)
        if min_cpm > 0:
            eff = np.where(self.cpm_gate_mask(b_ctrl_target, min_cpm), eff, 0.0)
        if gate_unexpressed:
            eff = np.where(b_ctrl_target > 0, eff, 0.0)
        p = np.clip(np.expm1(b_ctrl_target + eff), 0.0, None) / TS
        tot = p.sum()
        if tot <= 0:
            raise ValueError("degenerate profile")
        return p / tot

    def profile(self, pert: str, b_ctrl_target: np.ndarray, alpha: float, beta: float,
                scheme: str = "uniform", weights: dict[str, float] | None = None,
                min_cpm: float = 20.0, gate_unexpressed: bool = True) -> np.ndarray:
        eff = self.effect(pert, alpha, beta, weights=weights, scheme=scheme,
                          b_ctrl_target=b_ctrl_target)
        if min_cpm > 0:
            eff = np.where(self.cpm_gate_mask(b_ctrl_target, min_cpm), eff, 0.0)
        if gate_unexpressed:
            eff = np.where(b_ctrl_target > 0, eff, 0.0)
        p = np.clip(np.expm1(b_ctrl_target + eff), 0.0, None) / TS
        tot = p.sum()
        if tot <= 0:
            raise ValueError(f"degenerate profile for {pert!r}")
        return p / tot


def load_donors(spec: dict[str, str]) -> dict[str, Donor]:
    return {name: Donor.load(Path(path)) for name, path in spec.items()
            if Path(path).exists()}
