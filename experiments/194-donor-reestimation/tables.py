"""194 — turn accumulated weighted sums into donor effect tables, format-compatible
with `data/donors/<donor>/`.

Every scheme goes through the incumbent's own final arithmetic

    b_ctrl  = log1p(TS * Pc / sum Pc)
    delta_T = log1p(TS * S_T / sum S_T) - b_ctrl

so the ONLY difference between a variant and the incumbent is which cells contributed to
`S_T` and with what weight.  Two post-hoc variants need no extra pass:

`debias`  The incumbent's estimator is a non-linear function of a finite count vector, and
          `log1p` is concave, so `E[log1p(TS p_hat)] < log1p(TS p)` by
          `TS^2 p / (2 (1+TS p)^2 N_eff)`.  The control pool has 10^2-10^3 times more
          counts than any single target, so this bias does NOT cancel in the difference:
          it leaves a systematic, target-independent negative offset on low-expressed
          genes in every row of the table -- exactly the shape of a shared blob.  The
          correction is the second-order delta-method term, applied to both sides.

`ctrlfix` The Orion builder's control regex `^neg` matches the real gene NEGR1, so NEGR1's
          perturbed cells are pooled into the control and NEGR1 gets no row at all.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
TS = 5e4

POST_HOC = ("debias", "debias2", "ctrlfix")


def _profile(S: np.ndarray) -> np.ndarray:
    tot = S.sum(-1, keepdims=True)
    return np.log1p(TS * S / np.maximum(tot, 1e-300))


def _debias(b: np.ndarray, neff: np.ndarray, k: float = 1.0) -> np.ndarray:
    """Add back the second-order log1p bias.  `neff` is the effective total counts."""
    u = np.expm1(b)
    return b + k * 0.5 * TS * u / ((1.0 + u) ** 2 * np.maximum(neff, 1e-12))


class SumBank:
    def __init__(self, donor: str, mode: str):
        self.dir = WORK / f"sums-{donor}-{mode}"
        self.meta = json.loads((self.dir / "meta.json").read_text())
        self.donor = donor
        self.mode = mode
        self.names = np.array(self.meta["names"])
        self.ctrl_code = int(self.meta["ctrl_code"])
        self.schemes = list(self.meta["schemes"])
        self.neff = np.load(self.dir / "_neff.npy")
        self.cells = np.load(self.dir / "_cells.npy")
        self.genes = (ROOT / "data/donors" / donor / "genes.txt").read_text().split("\n")

    def sums(self, scheme: str) -> np.ndarray:
        return np.load(self.dir / f"{scheme}.npy", mmap_mode="r")

    def _base_sums(self, variant: str) -> tuple[np.ndarray, np.ndarray, str]:
        """Group sums for a scheme, a post-hoc variant of `v0`, or a DOSE `scheme@lam`.

        A dose is exact, not an approximation: every scheme is a weighted count sum and
        `_norm_within_target` puts every target's weights on a mean of 1, so the sums
        under `w = (1-lam) * 1 + lam * w_scheme` are exactly
        `(1-lam) * S_v0 + lam * S_scheme`.  `lam > 1` sharpens the scheme beyond its
        nominal dose; `lam < 1` softens it.  This makes a full dose-response curve free.
        """
        if "@" in variant and not variant.startswith("kdrescale"):
            name, lam = variant.split("@")
            lam = float(lam)
            # w = (1-lam) + lam*w_scheme must stay non-negative, and normalised weights
            # reach 0, so lam > 1 (sharpening beyond the nominal dose) is not expressible
            # as a cell weighting at all.  Interpolation only.
            if not 0.0 <= lam <= 1.0:
                raise ValueError(f"dose {lam} outside [0,1]: not a valid cell weighting")
            k0, k1 = self.schemes.index("v0"), self.schemes.index(name)
            # put both schemes' weights on a mean of 1 per target first, otherwise the
            # blend is dominated by whichever scheme happens to have the larger raw
            # weights (`cellmean`'s are 1/depth, i.e. ~1e-4 of `v0`'s, so an unnormalised
            # 25% dose of it is numerically indistinguishable from v0)
            wsum = np.load(self.dir / "_wsum.npy")
            cells = np.load(self.dir / "_cells.npy").astype(np.float64)
            s0 = (cells / np.maximum(wsum[k0], 1e-300))[:, None]
            s1 = (cells / np.maximum(wsum[k1], 1e-300))[:, None]
            S = ((1.0 - lam) * np.asarray(self.sums("v0"), dtype=np.float64) * s0
                 + lam * np.asarray(self.sums(name), dtype=np.float64) * s1)
            neff = (1.0 - lam) * self.neff[k0] + lam * self.neff[k1]
            return S, neff, name
        base = variant if variant in self.schemes else "v0"
        k = self.schemes.index(base)
        return np.asarray(self.sums(base), dtype=np.float64), self.neff[k].copy(), base

    def table(self, variant: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """(delta, b_ctrl, names) for a scheme name, a dose, or a post-hoc variant."""
        S, neff, _ = self._base_sums(variant)

        if variant.endswith("ctrlfix"):
            Pc = np.load(self.dir / "_ctrl_clean.npy")
        else:
            Pc = S[self.ctrl_code]
        b_ctrl = np.log1p(TS * Pc / Pc.sum())

        rows = np.array([i for i in range(len(self.names))])
        tot = S[rows].sum(1)
        good = tot > 0
        prof = np.zeros((len(rows), S.shape[1]))
        prof[good] = _profile(S[rows][good])

        if variant.endswith("debias") or variant.endswith("debias2"):
            kk = 2.0 if variant.endswith("debias2") else 1.0
            prof[good] = _debias(prof[good], neff[rows][good, None], kk)
            b_ctrl = _debias(b_ctrl, np.array(neff[self.ctrl_code]), kk)

        delta = prof[good] - b_ctrl
        if variant.startswith("kdrescale"):
            delta = delta * self._kd_scale(variant, rows[good])[:, None]
        return (delta.astype(np.float32), b_ctrl.astype(np.float32),
                self.names[rows][good])

    def _kd_scale(self, variant: str, rows: np.ndarray) -> np.ndarray:
        """Per-target de-attenuation by the target's own measured repression.

        A CRISPRi guide that leaves a fraction `r` of the target transcript standing --
        whether because some cells escaped or because repression is partial in all of
        them -- leaves a pseudobulk effect diluted roughly in proportion to the
        repression actually achieved, `1 - r`.  Dividing each row by `1 - r` puts every
        donor's estimate of the same target on a common scale BEFORE `build_baseline`
        averages the donors, which matters because the three donors' median residual
        expression differs by 3.3x (K562 0.155, HCT116 0.300, HEK293T 0.510) and the
        equal-weight consensus is otherwise dominated, target by target, by whichever
        donor happened to repress hardest.

        `pds_cosine` is exactly scale-invariant per row, so this can only act through the
        cross-donor average -- which is precisely the claim being tested.
        `kdrescale_perm` shuffles `r` across targets and must not reproduce the effect.
        """
        kd = np.load(self.dir / "_kd_target.npy")
        power = float(variant.split("@")[1]) if "@" in variant else 1.0
        r = kd[0][rows]
        med = float(np.nanmedian(r))
        r = np.where(np.isfinite(r), r, med)
        if variant.split("@")[0].endswith("_perm"):
            rng = np.random.default_rng(20260903)
            r = rng.permutation(r)
        floor = 0.25          # never inflate a row by more than 4x
        return (1.0 / np.clip(1.0 - r, floor, 1.0)) ** power

    def write(self, variant: str, out: Path) -> dict:
        delta, b_ctrl, names = self.table(variant)
        out.mkdir(parents=True, exist_ok=True)
        np.save(out / "delta.npy", delta)
        np.save(out / "b_ctrl.npy", b_ctrl)
        (out / "perts.txt").write_text("\n".join(map(str, names)))
        (out / "genes.txt").write_text("\n".join(self.genes))
        norms = np.linalg.norm(delta, axis=1)
        meta = dict(source=f"194 variant {variant} of {self.donor}", variant=variant,
                    donor=self.donor, mode=self.mode, target_sum=TS,
                    n_perturbations=int(len(names)), n_genes=int(delta.shape[1]),
                    effect_norm=dict(min=float(norms.min()),
                                     median=float(np.median(norms)),
                                     max=float(norms.max())))
        (out / "meta.json").write_text(json.dumps(meta, indent=2))
        return meta
