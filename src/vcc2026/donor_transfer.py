"""Donor-transfer model: predict a perturbation response in an unseen cell context.

The shape follows directly from experiment 003. A donor line's effect decomposes as

    delta_g = m + s_g

where `m` is the line's mean perturbation response and `s_g` the gene-specific deviation.
Measured K562 -> RPE1: `m` transfers at cos 0.49, `s` at cos 0.128 against a null of 0.003.
**`m` earns exactly zero score** -- it is the context-mean baseline, the 0 end of every
metric's scale -- so the two components have completely different jobs and are kept apart
here with separate coefficients.

    predicted effect  =  beta * m  +  alpha * s_g

`alpha` is the lever from experiments 001/002: `pds_cosine` is exactly invariant to it, while
`expr_mse` is quadratic and the DE metrics rise with it. `beta` controls how much of the
free-but-unscored mean response to add; it moves `mse` and the DE calls without touching
discrimination.

Gene axes are matched by symbol. The donor carries 8,248 genes against the target's 18,533,
so ~29% of the DE-scored genes get no donor signal and are predicted unchanged.

Effects are treated as log fold changes. `delta_j = log1p(TS*p_pert) - log1p(TS*p_ctrl)`, and
for a well-expressed gene log1p ~ log, so delta ~ log(p_pert/p_ctrl) -- a ratio, which is
invariant to the donor's gene-set normalisation and therefore portable to a different axis.
The approximation degrades for low-expressed genes; `--space foldchange` applies it as an
explicit ratio, `--space log1p` adds it in the space `expr_mse` actually measures.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

TS = 5e4


def eb_shrink(s: np.ndarray, b_ctrl: np.ndarray, n_cells: int,
              depth: float = 20000.0, noise_scale: float = 1.0) -> np.ndarray:
    """Empirical-Bayes shrinkage of a donor effect matrix.

    Nadig 2025: within K562, the raw log-fold-change correlation between two *replicate*
    experiments is 0.16, while the noise-corrected true correlation is 0.90. Donor deltas are
    overwhelmingly sampling noise.

    Per gene, the observed across-perturbation variance `V` decomposes into true-effect
    variance `tau2` and Poisson sampling variance `sigma2`, and the Bayes-optimal shrinkage
    is `tau2 / (tau2 + sigma2)`. For a log-pseudobulk delta the Poisson variance is ~1/N with
    `N ~ n_cells * depth * p_gene`, and `p_gene` comes from the donor's own control profile --
    so this is estimated **entirely from the donor**, with no target data. That is what
    separates it from the supervised per-gene weights that failed in experiment 004: there is
    nothing here to overfit.

    Measured on h1-test (full gene axis): 0.0706 -> 0.0750 alone, 0.0809 -> 0.0823 combined
    with multi-donor averaging. Supersedes `gene_norm`, which does overlapping work and is
    worse (0.0715 alone, and 0.0737 in the multi-donor setting where it actively hurts).
    """
    e = np.expm1(b_ctrl)
    p_gene = e / e.sum() if e.sum() > 0 else np.full_like(e, 1.0 / e.size)
    sigma2 = noise_scale / np.maximum(n_cells * depth * p_gene, 1e-3)
    V = s.var(axis=0)
    tau2 = np.maximum(V - sigma2, 0.0)
    return (s * (tau2 / (tau2 + sigma2 + 1e-12))).astype(np.float32)


def _norm_scale(s: np.ndarray, floor_pct: float = 25.0) -> np.ndarray:
    """Per-gene divisor for effect normalisation, floored.

    Dividing each gene by its mean |effect| equalises how much a gene can steer the
    predicted direction (experiment 006, +3.4%). But a donor whose axis includes
    unexpressed genes has genes with mean |effect| ~ 0, and dividing by that turns pure
    noise into the loudest signal in the vector. Orion carries 3,372 zero-expression genes
    on a 38,606-gene axis and unfloored normalisation dropped its transfer cosine from
    0.049 to 0.019. Flooring at a low percentile keeps the equalising effect and removes
    the blow-up.
    """
    sc = np.abs(s).mean(axis=0)
    live = sc[sc > 0]
    floor = np.percentile(live, floor_pct) if live.size else 1e-9
    return np.maximum(sc, floor)


@dataclass
class Donor:
    """Per-perturbation effect vectors from one donor screen."""

    delta: np.ndarray          # (n_perts, n_genes) log-pseudobulk effect vs control
    perts: np.ndarray          # (n_perts,) gene symbols
    genes: np.ndarray          # (n_genes,) gene symbols
    b_ctrl: np.ndarray         # (n_genes,) donor's own control profile

    @classmethod
    def load(cls, d: str | Path) -> "Donor":
        d = Path(d)
        return cls(
            delta=np.load(d / "delta.npy"),
            perts=np.asarray((d / "perts.txt").read_text().split("\n")),
            genes=np.asarray((d / "genes.txt").read_text().split("\n")),
            b_ctrl=np.load(d / "b_ctrl.npy"),
        )

    def decompose(self) -> tuple[np.ndarray, np.ndarray]:
        """Split into the shared mean response `m` and gene-specific deviations `s`."""
        m = self.delta.mean(axis=0)
        return m, self.delta - m


class TransferModel:
    """Maps donor effects onto a target context's gene axis and basal profile."""

    def __init__(self, donor: Donor, target_genes: np.ndarray, fallback: str = "mean",
                 gene_norm: bool = True):
        self.donor = donor
        self.target_genes = np.asarray(target_genes)
        m, s = donor.decompose()

        # Column map: donor gene -> target column, for the shared symbols only.
        tindex = {g: i for i, g in enumerate(self.target_genes)}
        cols, dcols = [], []
        for j, g in enumerate(donor.genes):
            t = tindex.get(g)
            if t is not None:
                cols.append(t)
                dcols.append(j)
        self.tcols = np.asarray(cols, dtype=np.int64)
        self.dcols = np.asarray(dcols, dtype=np.int64)

        self.m_full = np.zeros(len(self.target_genes), dtype=np.float32)
        self.m_full[self.tcols] = m[self.dcols]

        # Per-gene normalisation (experiment 006). Cosine is driven by whichever genes carry
        # the largest values, so a gene that is merely noisy in the donor steals direction
        # without contributing signal. Dividing by each gene's mean |effect| across the donor
        # panel equalises that. Measured +3.4% on h1-test and +5.4% on h1-trainval.
        # Rescaled to preserve the median effect norm so `alpha` keeps its meaning.
        if gene_norm:
            scale = _norm_scale(s)
            s_n = s / scale
            before = np.median(np.linalg.norm(s[:, self.dcols], axis=1))
            after = np.median(np.linalg.norm(s_n[:, self.dcols], axis=1)) or 1.0
            s = (s_n * (before / after)).astype(np.float32)
        self.gene_norm = gene_norm
        self.s = s
        self.pert_index = {p: i for i, p in enumerate(donor.perts)}
        self.fallback = fallback
        # Typical gene-specific effect magnitude, on the target axis, for scaling the
        # fallback so uncovered perturbations emit an effect of realistic size.
        norms = np.linalg.norm(s[:, self.dcols], axis=1)
        self._median_s_norm = float(np.median(norms))
        self._m_norm = float(np.linalg.norm(self.m_full)) or 1.0

    @property
    def coverage(self) -> float:
        return len(self.tcols) / len(self.target_genes)

    def covered(self, perts) -> np.ndarray:
        return np.asarray([p in self.pert_index for p in perts])

    def effect(self, pert: str, alpha: float, beta: float) -> np.ndarray:
        """Predicted log-fold-change vector on the target gene axis.

        A perturbation the donor does not carry cannot get a gene-specific effect. Leaving
        it at zero is the worst option available: the prediction then equals the control,
        the Wilcoxon test finds nothing, `n_pred` collapses and `fid` takes its floor --
        `F_p = k/max(n_pred, n_real) = 0` whenever the reference found anything. With
        `fallback="mean"` such a perturbation instead emits the donor's mean response,
        scaled to a typical effect magnitude: chance-level directions but full coverage,
        which is worth roughly the baseline rather than the floor.
        """
        eff = beta * self.m_full
        i = self.pert_index.get(pert)
        if i is not None:
            if alpha != 0.0:
                add = np.zeros_like(eff)
                add[self.tcols] = self.s[i][self.dcols]
                eff = eff + alpha * add
        elif self.fallback == "mean" and alpha != 0.0:
            eff = eff + (alpha * self._median_s_norm / self._m_norm) * self.m_full
        return eff

    def cpm_gate_mask(self, b_ctrl_target: np.ndarray, min_cpm: float) -> np.ndarray:
        """Genes above `min_cpm` in the target's controls.

        Experiment 006: gating the prediction at target CPM > 20 raised transfer cosine on
        both harnesses independently (h1-test +6.4%, h1-trainval +7.2%), with the optimum at
        20 in both. Below that a gene contributes noise to the predicted direction and is
        excluded from the DE metrics anyway (their filter is 5 CPM in control).
        """
        if min_cpm <= 0:
            return np.ones(len(b_ctrl_target), dtype=bool)
        e = np.expm1(b_ctrl_target)
        return (1e6 * e / e.sum()) > min_cpm

    def profile(
        self,
        pert: str,
        b_ctrl_target: np.ndarray,
        alpha: float,
        beta: float,
        space: str = "log1p",
        gate_unexpressed: bool = True,
        min_cpm: float = 20.0,
    ) -> np.ndarray:
        """Predicted relative expression profile (sums to 1) for one perturbation.

        `b_ctrl_target` is the target context's own control profile in log1p(TS*p) space,
        which is what identifies the context to the model -- the only thing we are given.
        """
        eff = self.effect(pert, alpha, beta)

        if min_cpm > 0:
            eff = np.where(self.cpm_gate_mask(b_ctrl_target, min_cpm), eff, 0.0)

        if gate_unexpressed:
            # A gene with no reads in the target's controls cannot be measured to change.
            # Predicting motion there is pure loss on every DE metric.
            eff = np.where(b_ctrl_target > 0, eff, 0.0)

        if space == "log1p":
            p = np.clip(np.expm1(b_ctrl_target + eff), 0.0, None) / TS
        elif space == "foldchange":
            p = np.clip(np.expm1(b_ctrl_target), 0.0, None) / TS * np.exp(eff)
        else:
            raise ValueError(f"unknown space {space!r}")

        tot = p.sum()
        if tot <= 0:
            raise ValueError(f"degenerate profile for {pert!r}")
        return p / tot


def control_profile(X, TS_: float = TS) -> np.ndarray:
    """log1p(TS * p) of a control cell matrix -- the context's basal signature."""
    P = np.asarray(X.sum(0)).ravel().astype(np.float64)
    return np.log1p(TS_ * P / P.sum()).astype(np.float64)


def _repair_totals(block, totals, rng):
    """Force each gene's column total to `totals` exactly, moving single counts between cells.

    This is the same repair the quasi emitter applies to its scattered part, and it exists for
    the same reason: any residual in a gene's total lands directly in the pseudobulk, which is
    what `pds`, `mse` and `nmae` read.
    """
    out = np.asarray(block, dtype=np.int64).copy()
    n = out.shape[0]
    resid = np.asarray(totals, dtype=np.int64) - out.sum(0)
    for j in np.flatnonzero(resid):
        d = int(resid[j])
        if d > 0:
            cells = rng.integers(0, n, d)
            np.add.at(out[:, j], cells, 1)
        else:
            for _ in range(-d):
                nz = np.flatnonzero(out[:, j] > 0)
                if nz.size == 0:
                    break
                out[rng.choice(nz), j] -= 1
    return out


def _stochastic_round(x, rng):
    """Round to integers preserving the expectation, so no systematic mass is destroyed.

    Flooring independently would silently lose up to one count per cell, which drives the
    emitted law off the control's and inverts the called direction.
    """
    fl = np.floor(x)
    return (fl + (rng.random(x.shape) < (x - fl))).astype(np.int64)


def _quasi_block(profile_p, n, G, L, spread, overdispersion, rng):
    """The shipped quasi law as a standalone block, byte-for-byte the same construction."""
    sp = np.clip(np.asarray(spread, dtype=np.float64), 0.0, 1.0)
    totals = np.floor(n * L * profile_p).astype(np.int64)
    t_rand = np.floor(totals * sp).astype(np.int64)
    t_even = totals - t_rand
    block = np.zeros((n, G), dtype=np.int64)
    nz = np.flatnonzero(t_even > 0)
    if nz.size:
        base, rem = np.divmod(t_even[nz], n)
        offsets = np.concatenate(([0], np.cumsum(rem)[:-1])) % n
        pos = (np.arange(n)[:, None] - offsets[None, :]) % n
        block[:, nz] = base[None, :] + (pos < rem[None, :])
    rz = np.flatnonzero(t_rand > 0)
    if rz.size:
        mu = np.broadcast_to(t_rand[rz] / n, (n, rz.size))
        od = np.asarray(overdispersion, dtype=np.float64)
        if od.ndim:
            if od.size != G:
                raise ValueError("per-gene overdispersion needs one value per gene")
            od = np.clip(od[rz], 0.0, None)[None, :]
        if np.any(od > 0.0):
            safe = np.where(od > 0.0, od, 1.0)
            lam = np.where(
                np.broadcast_to(od > 0.0, mu.shape),
                rng.gamma(shape=1.0 / safe, scale=safe * mu),
                mu,
            )
            draw = rng.poisson(lam)
        else:
            draw = rng.poisson(mu)
        resid = t_rand[rz] - draw.sum(0)
        for j in np.flatnonzero(resid):
            d = int(resid[j])
            if d > 0:
                cells = rng.integers(0, n, d)
                np.add.at(draw[:, j], cells, 1)
            else:
                for _ in range(-d):
                    nzc = np.flatnonzero(draw[:, j] > 0)
                    if nzc.size == 0:
                        break
                    draw[rng.choice(nzc), j] -= 1
        block[:, rz] += draw
    return block


def emit_counts(profile_p, lib_sizes, dispersion, rng, max_counts: int = 1_000_000,
                mode: str = "nb", tight_mask=None, spread=None, overdispersion=0.0,
                control_counts=None, control_profile=None, param_mask=None):
    """Sample count vectors from a relative-expression profile.

    Four of the six scored metrics are recomputed from these counts by a Wilcoxon test
    against the reference controls, so **the cell-to-cell spread we emit is a scoring
    decision, not a realism decision**. It sets how many genes clear the BH bar, which
    drives `fid`'s yield term `min(1, n_pred/n_real)`, `reach`'s depth, and `jac`.

    mode="nb"    gamma-Poisson. Realistic single-cell noise. The Wilcoxon then calls only
                 a few hundred genes, `n_pred << n_real`, and `fid` takes its coverage
                 penalty -- measured at -0.77 to -0.86 on the official scale (exp 005).
    mode="mixed" per-gene choice between the two, via `tight_mask`. This is a control
                 channel the other metrics cannot see: a gene's Wilcoxon p-value depends on
                 its mean effect **and** on its cell-to-cell spread, and we set them
                 independently. Genes we trust are emitted at the variance floor, so they
                 clear significance and rank first in `reach`'s ordering; genes we do not are
                 emitted with realistic overdispersion, so they fall down the p-value order
                 and out of our called set. **The mean profile is identical either way**, so
                 `pds`, `mse` and `nmae` should not move while `fid` and `reach` do.

                 The reason a minimum-variance emitter calls essentially every gene is not
                 the effect sizes -- it is that 400 near-identical values sit at the control
                 distribution's *mean*, which for skewed counts is above its median, so the
                 rank-sum against 18,400 real controls is significant even at zero effect.
                 Restoring realistic spread on a gene removes that artefact.

    mode="even"  minimum-variance allocation: every cell gets the same library size, and
                 each gene's panel total `round(n * L * p_g)` is spread as evenly as
                 possible over the n cells (floor everywhere, the remainder given one
                 extra each). Cell-to-cell variance is as low as integer counts allow, so
                 the test has maximum power and `n_pred` covers the DE universe.
                 Over-calling costs `jac`, which divides by the union -- the trade is
                 deliberate and must be checked against the scorer, not assumed.
    """
    import scipy.sparse as sp

    n = len(lib_sizes)
    G = profile_p.size

    if mode == "quasi":
        # Sub-Poisson, per gene, with the gene's panel total held fixed.
        #
        # mode="mixed" moves a gene all the way from the integer variance floor to full
        # Poisson, which does not reorder its p-value -- it removes it from the called set
        # entirely. `jac` and `fid` then pay for whatever `reach` gains (exp 027). What is
        # actually wanted is a *dial*: enough extra spread to push a gene down the p-value
        # ordering, not enough to lose its significance.
        #
        # Split each gene's total `T` into an evenly-allocated part and a randomly-scattered
        # part, `spread` being the fraction scattered. Variance per cell is then ~spread*T/n
        # against Poisson's T/n, so `spread` is literally the fraction of Poisson noise. The
        # gene's total is preserved, so the pseudobulk profile -- and therefore `pds`, `mse`
        # and `nmae` -- is untouched by construction.
        if spread is None:
            raise ValueError("mode='quasi' requires spread")
        sp = np.clip(np.asarray(spread, dtype=np.float64), 0.0, 1.0)
        L = int(min(np.median(lib_sizes), max_counts))
        totals = np.floor(n * L * profile_p).astype(np.int64)
        t_rand = np.floor(totals * sp).astype(np.int64)
        t_even = totals - t_rand

        block = np.zeros((n, G), dtype=np.int64)
        nz = np.flatnonzero(t_even > 0)
        if nz.size:
            base, rem = np.divmod(t_even[nz], n)
            offsets = np.concatenate(([0], np.cumsum(rem)[:-1])) % n
            pos = (np.arange(n)[:, None] - offsets[None, :]) % n
            block[:, nz] = base[None, :] + (pos < rem[None, :])
        rz = np.flatnonzero(t_rand > 0)
        if rz.size:
            # Scatter the random part, then repair each gene's total EXACTLY.
            #
            # A plain Poisson draw only preserves the gene total in expectation, and the
            # residual lands directly in the pseudobulk -- which is what `pds`, `mse` and
            # `nmae` read. Measured on h1-fast at spread 0.25, that alone moved `pds`
            # between 0.8404 and 0.8583 across two seeds, i.e. it put our best metric at
            # risk to buy `reach`. Repairing the totals makes the mean profile bit-identical
            # to the even emitter's, so only the DE metrics -- which read cell-to-cell
            # spread -- can move at all.
            mu = np.broadcast_to(t_rand[rz] / n, (n, rz.size))
            # `overdispersion` may be a scalar or one value per gene, mirroring `spread`.
            # A per-gene phi is the principled form -- each gene's own control cells have
            # their own dispersion, and a single global phi over- or under-disperses most of
            # them -- so the plumbing accepts it even though the shipped arms use a scalar.
            od = np.asarray(overdispersion, dtype=np.float64)
            if od.ndim:
                if od.size != G:
                    raise ValueError("per-gene overdispersion needs one value per gene")
                od = np.clip(od[rz], 0.0, None)[None, :]
            if np.any(od > 0.0):
                # Gamma-Poisson (negative binomial) scatter: Var = mu + phi*mu^2.
                #
                # `spread` caps at Poisson, and experiment 095 measured that this is far too
                # weak to reorder anything -- predicted p-values span hundreds of orders of
                # magnitude and Poisson noise moves them by a few, so a spread policy changes
                # which genes are CALLED long before it changes their ORDER. Over-dispersion
                # is the one emission lever that moves p-values enough to matter, and it is
                # also what the data says we are missing: measured against the real challenge
                # controls, our per-gene CV across 400 emitted cells is 0.38-0.40 against
                # 0.84-0.92 in the real control cells -- we are under-dispersed ~2.3x.
                #
                # Realised on the depth-repaired harness through this exact path, against a
                # 5e-4 same-policy noise floor: phi=0.30 takes reach 0.0838 -> 0.1888 and
                # phi=0.45 -> 0.1620, with `jac` up throughout and `nmae` invariant.
                #
                # The exact-total repair below is untouched, so `pds`, `mse` and `nmae` stay
                # bit-identical exactly as they do at any `spread`.
                safe = np.where(od > 0.0, od, 1.0)
                lam = np.where(
                    np.broadcast_to(od > 0.0, mu.shape),
                    rng.gamma(shape=1.0 / safe, scale=safe * mu),
                    mu,
                )
                draw = rng.poisson(lam)
            else:
                draw = rng.poisson(mu)
            resid = t_rand[rz] - draw.sum(0)
            for j in np.flatnonzero(resid):
                d = int(resid[j])
                if d > 0:
                    cells = rng.integers(0, n, d)
                    np.add.at(draw[:, j], cells, 1)
                else:
                    # remove only where there is something to remove, so counts stay >= 0
                    for _ in range(-d):
                        nz = np.flatnonzero(draw[:, j] > 0)
                        if nz.size == 0:
                            break
                        draw[rng.choice(nz), j] -= 1
            block[:, rz] += draw

        keep = block > 0
        indptr = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(keep.sum(1), out=indptr[1:])
        rows, cols = np.nonzero(keep)
        order = np.lexsort((cols, rows))
        return (block[rows[order], cols[order]].astype(np.float32),
                cols[order].astype(np.int32), indptr)

    if mode == "resample_mix":
        # Per-gene choice between the parametric quasi law and RESAMPLED real control cells.
        #
        # Why: the quasi emitter places all `n` cells at one fixed library size
        # (`L = median(lib_sizes)`), so emitted CPM lives on a lattice the real controls do
        # not. Against the scorer's Mann-Whitney that is not null under the null -- measured
        # on the challenge's own controls, ~19% of gated genes clear BH at ZERO predicted
        # effect and ~88% of those are called UP. Resampling real control cells and applying
        # the predicted fold change makes emitted and control exchangeable at zero effect, so
        # `Fbar = 0.5` identically rather than asymptotically, and it inherits the real
        # library-size distribution, skew and zero-inflation for free.
        #
        # It also raises the `#247` dispersion credit: `C_hat` moves 0.656 -> ~1.03 of the
        # real-cell value, on a cap measured to fire on 0 of 900 rows.
        #
        # `param_mask` keeps the top genes by |predicted effect| on the parametric law: those
        # are the calls we have signal on (top-decile sign accuracy 0.74 against 0.51 in the
        # tail), and holding them back keeps `n_pred` off the floor while the uninformative
        # tail is silenced.
        if control_counts is None or control_profile is None or param_mask is None:
            raise ValueError(
                "mode='resample_mix' requires control_counts, control_profile and param_mask"
            )
        if spread is None:
            raise ValueError("mode='resample_mix' requires spread (for the parametric part)")
        pm = np.asarray(param_mask, dtype=bool)
        if pm.size != G:
            raise ValueError("param_mask needs one value per gene")
        C = np.asarray(control_counts, dtype=np.float64)
        if C.shape != (n, G):
            raise ValueError(f"control_counts must be ({n}, {G}), got {C.shape}")
        cp = np.asarray(control_profile, dtype=np.float64)
        if cp.size != G:
            raise ValueError("control_profile needs one value per gene")

        L = int(min(np.median(lib_sizes), max_counts))
        totals = np.floor(n * L * profile_p).astype(np.int64)

        # -- parametric half: the shipped quasi law, unchanged --
        param_block = _quasi_block(profile_p, n, G, L, spread, overdispersion, rng)

        # -- resampled half: real cells carrying the predicted fold change --
        # fold = predicted relative expression / control relative expression, applied
        # multiplicatively, then each cell renormalised to ITS OWN library so the
        # composition shift does not smear across genes under the scorer's CPM step.
        fold = np.divide(profile_p, cp, out=np.ones_like(cp), where=cp > 0)
        scaled = C * fold[None, :]
        row_now = scaled.sum(1)
        row_want = C.sum(1)
        rescale = np.divide(row_want, row_now, out=np.ones_like(row_now), where=row_now > 0)
        scaled *= rescale[:, None]
        resamp_block = _stochastic_round(scaled, rng)

        # -- splice, THEN repair --
        # The two halves carry different per-cell library structure, so splicing genes
        # between them breaks the pseudobulk pinning. Repairing per-gene totals AFTER the
        # mix restores it: measured rel L2 against the predicted profile is 0.0155 without
        # this ordering and 0.0015 with it, and 1.55% profile drift is the exact channel
        # that cost a19 -0.0959 on `mse`.
        block = np.where(pm[None, :], param_block, resamp_block)
        block = _repair_totals(block, totals, rng)

        keep = block > 0
        indptr = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(keep.sum(1), out=indptr[1:])
        rows, cols = np.nonzero(keep)
        order = np.lexsort((cols, rows))
        return (block[rows[order], cols[order]].astype(np.float32),
                cols[order].astype(np.int32), indptr)

    if mode == "mixed":
        if tight_mask is None:
            raise ValueError("mode='mixed' requires tight_mask")
        tight = np.asarray(tight_mask, dtype=bool)
        L = int(min(np.median(lib_sizes), max_counts))
        block = np.zeros((n, G), dtype=np.int64)

        # tight genes: even allocation, as in mode="even"
        totals = np.floor(n * L * np.where(tight, profile_p, 0.0)).astype(np.int64)
        nz = np.flatnonzero(totals > 0)
        if nz.size:
            base, rem = np.divmod(totals[nz], n)
            offsets = np.concatenate(([0], np.cumsum(rem)[:-1])) % n
            pos = (np.arange(n)[:, None] - offsets[None, :]) % n
            block[:, nz] = base[None, :] + (pos < rem[None, :])

        # loose genes: gamma-Poisson at the same mean, so the pseudobulk is unchanged in
        # expectation and only the spread differs
        loose = np.flatnonzero(~tight & (profile_p > 0))
        if loose.size:
            mu = np.outer(np.minimum(lib_sizes, max_counts), profile_p[loose])
            if dispersion:
                mu = rng.gamma(shape=1.0 / dispersion, scale=mu * dispersion)
            block[:, loose] = rng.poisson(mu)

        keep = block > 0
        indptr = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(keep.sum(1), out=indptr[1:])
        rows, cols = np.nonzero(keep)
        order = np.lexsort((cols, rows))
        return (block[rows[order], cols[order]].astype(np.float32),
                cols[order].astype(np.int32), indptr)

    if mode == "even":
        L = int(min(np.median(lib_sizes), max_counts))
        totals = np.floor(n * L * profile_p).astype(np.int64)
        nz = np.flatnonzero(totals > 0)
        if nz.size == 0:
            return (np.zeros(0, np.float32), np.zeros(0, np.int32), np.zeros(n + 1, np.int64))
        base, rem = np.divmod(totals[nz], n)

        # Fully vectorised: the block is only (n_cells x n_nonzero_genes) -- ~400 x 7,700
        # here, i.e. ~12 MB -- so build it dense and convert, rather than appending per
        # (cell, gene). Each gene contributes `base` to every cell plus one extra to `rem`
        # of them; rotating which cells get the extra keeps library sizes equal.
        offsets = np.concatenate(([0], np.cumsum(rem)[:-1])) % n
        pos = (np.arange(n)[:, None] - offsets[None, :]) % n
        block = base[None, :].astype(np.int64) + (pos < rem[None, :])

        keep = block > 0
        counts_per_row = keep.sum(1)
        indptr = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(counts_per_row, out=indptr[1:])
        rows, cols = np.nonzero(keep)
        order = np.lexsort((nz[cols], rows))
        return (block[rows[order], cols[order]].astype(np.float32),
                nz[cols[order]].astype(np.int32),
                indptr)

    data, indices, indptr = [], [], [0]
    for lib in lib_sizes:
        mu = profile_p * float(min(lib, max_counts))
        if dispersion:
            mu = rng.gamma(shape=1.0 / dispersion, scale=mu * dispersion)
        c = rng.poisson(mu)
        idx = np.flatnonzero(c)
        vals = c[idx]
        if vals.sum() > max_counts:
            vals = np.floor(vals * (max_counts / vals.sum())).astype(np.int64)
            keep = vals > 0
            idx, vals = idx[keep], vals[keep]
        data.append(vals.astype(np.float32))
        indices.append(idx.astype(np.int32))
        indptr.append(indptr[-1] + idx.size)
    return (np.concatenate(data), np.concatenate(indices), np.asarray(indptr))
