"""272 -- the NOISE-CORRECTED whitening operator: the constructions, FROZEN before any fold is measured (this file's
sha256 is recorded in every output; analyse272.py asserts it did not change between folds).

Setting.  A PLOLO fold's operator W = V diag(lambda^(-gamma/2)) V^T is fitted (062 `fit_union_whitener`) on the two
donors' full-depth centred rows b_t (every perturbation of both banks, ~28k rows on p = 17,942 / 17,806 columns);
gamma 0.65, lam 0.  271 section 1.2: the donors' batch-parity HALF rows A_t, B_t whitened through it ANTI-correlate
(median cos(WA, WB) -0.13...-0.24) while the held line's halves read +0.03...+0.07 -- the operator has learnt the
sampling-noise realisation of its own training rows (the low-eigenvalue directions are the ones in which each fitted
row's halves cancel).

A. DIAGONAL WIENER DAMPING IN THE FITTED EIGENBASIS (`noise_spectrum`, `DampedOperator`).
   Noise covariance from the half differences of the donors' carried panel rows:
       d_t = HALF_FACTOR * (A_t - B_t) / 2,       N = sum_j w_j * mean_t d_{jt} d_{jt}^T
   with A_t, B_t the half rows in the production row's units (271 `work/half/rows-<line>.npz`, `real` == the
   production un-whitened donor row at float32) converted to the operator's BASIS units by ref / EFFECT_SCALE
   (062: basis row = (delta - mean) * ref / raw_median_j; production row = (delta - mean) * EFFECT_SCALE /
   raw_median_j; ref = median over the fold's donors of raw_median_j), and w_j = the donor's share of the basis rows
   (the covariance the operator was fitted on weights the donors by row count).
   HALF_FACTOR = 1.0: if the full row is the mean of two independent halves, real = (A + B)/2 (271 measured
   cos((A+B)/2, real) median 0.996, norm ratio 1.000-1.003), then Var(real - signal) = (Var e_A + Var e_B)/4 =
   Var((A - B)/2) exactly -- d_t is a noise realisation AT FULL DEPTH, no half-vs-full factor is needed.  271's
   1.37-1.39x is the HALF rows' norm over the full row's (each half carries twice the full row's noise variance),
   not a factor on d.  The identity is checked empirically per donor: mean ||d||^2 vs (mean ||A||^2 + mean ||B||^2)/2
   - mean ||real||^2 (both equal the full-row noise energy under the model; reported as `half_factor_check`).
   Batch-parity halves put any between-batch effect into d as well as the sampling noise (267 section 1); N is
   therefore the covariance of everything in a donor row that does NOT reproduce across its batches -- the
   conservative (larger) noise reading.
   Per direction i (eigenvector v_i, eigenvalue lambda_i):  n_i = v_i^T N v_i,  s_i = clip((lambda_i - n_i) /
   lambda_i, 0, 1) -- the Wiener signal share of the fitted variance in that direction.  The damped operator
   multiplies direction i's weight lambda_i^(-gamma/2) by s_i^kappa, kappa in KAPPAS (0.5, 1, 2; 0.25 and 4 as the
   dose curve); kappa = 0 is the incumbent, bit for bit (`DampedOperator.apply` is `plolo.Operator.apply`'s
   arithmetic with one extra multiply that is skipped at kappa = 0).  The preserve (median row norm on the operator
   columns restored by ONE scalar) is production's, so the arm is shape-only up to the render.
   Controls: `shuffled_s` (s_i permuted across directions, seed SHUF_SEED) and `mirrored_s` (the same s profile
   reversed along the eigenvalue rank: what the arm does to the bottom of the spectrum is done to the top);
   `hard_s` (s_i >= 0.5 -> 1 else 0: the truncation reading) is reported as a comparison, not an arm.
   A'. THE EIGENVALUE FLOOR AT THE NOISE LEVEL (added after the FIRST read of fold orion_hct116 with the v1 file, which
   found n_i > lambda_i on 77 % of the directions -- the in-sample eigenvalue below the independent noise estimate, the
   Marchenko-Pastur signature of the learnt noise -- so that s^kappa zeroes three quarters of the spectrum; every fold
   was then re-run with THIS file, one sha for all).  The soft correction of the same finding: whiten by the
   OUT-OF-SAMPLE variance instead of the in-sample one, lambda'_i = max(lambda_i, c * n_i), w_i = lambda'_i^(-gamma/2),
   c in FLOOR_C (1 = the measured noise level; 0.25...4 the dose curve); c = 0 is the incumbent.  It damps the
   noise-dominated directions by (lambda_i / (c n_i))^(gamma/2) instead of removing them.  Controls: n_i permuted
   across directions, and n_i mirrored along the eigenvalue rank.

B. CROSS-FITTED OPERATORS (`fit_operator`, `apply_factors`, `preserve`).
   The fold operator refitted with 062's recipe (every row of each bank, centred on the bank's column mean over the
   last-source columns of the operator axis, scaled by ref / raw_median_bank, C = B^T B / N accumulated in row
   blocks, eigh, floor max * 1e-12, eigenvectors stored float32) on (i) the FULL banks -- the KA-2 reproduction of
   the cached held operator -- and (ii) the half-A banks and (iii) the half-B banks (267's batch-parity halves,
   rebuilt from 267's sums with 271's `write_half_bank`, verbatim below).  Tests: cos(W_A A_t, W_A B_t) (A in-sample,
   B out) and its mirror, the held line's halves through each, against the fully in-sample cos(W A_t, W B_t)
   (KA-0).  Predictions: `xfit_avg` = the production pre through (W_A + W_B)/2 (the mandate's construction: the
   mean of the two preserved outputs, one preserve); `half:xfit` = mean over the carrying donors of
   (A_t W_B + B_t W_A)/2 -- each half whitened by the operator that never saw it -- against its matched baseline
   `half:inc` = mean_j (A_t + B_t)/2 through the cached operator, and `half:infit` = (A_t W_A + B_t W_B)/2 (each half
   through the operator fitted ON it, the maximally in-sample reading).

C. THE PRODUCTION OPERATOR, DIAGNOSTIC ONLY (`prod272.py`): the same s_i spectrum for
   data/donors/whitener_union_pruned.npz from the THREE donors' half rows (basis weights 7784 / 10807 / 10803 from
   its metadata; ref = median of the three raw medians), and the whitened energy of the kolf2 fold's cached `pre`
   (the three-donor production rows on 273 panel targets) in the noise-dominated directions.  No production file is
   touched.
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/272-noise-corrected-operator"
WORK = HERE / "work"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
E267 = ROOT / "experiments/267-replicate-validated-oracle"
HALF267 = E267 / "work/half"
EFFECT_SCALE = 6.817            # E062.EFFECT_SCALE
HALF_FACTOR = 1.0               # see the module docstring
SHUF_SEED = 20260906
DERANGE_SEED = 20260906
KAPPAS = (0.5, 1.0, 2.0)
KAPPAS_EXTRA = (0.25, 4.0)
FLOOR_C = (0.25, 0.5, 1.0, 2.0, 4.0)
TS = 5e4                        # halfrows271
MIN_HALF = 5                    # halfrows271


def sha256_self() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for blk in iter(lambda: fh.read(8 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def row_cos(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(1) / (na * nb)
    return np.where((na > 0) & (nb > 0), c, np.nan)


def derangement(n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    while True:
        p = rng.permutation(n)
        if n < 2 or not np.any(p == np.arange(n)):
            return p


# ------------------------------------------------------------------------------------------ the half rows (271)
def load_half_rows(line: str, targets) -> dict:
    """271's half rows for `line`, re-indexed to `targets` (production units on the official axis)."""
    z = np.load(E271 / f"work/half/rows-{line}.npz", allow_pickle=True)
    tz = [str(t) for t in z["targets"]]
    idx = np.asarray([tz.index(t) for t in targets], dtype=np.int64)
    meta = json.loads((E271 / f"work/half/rows-{line}.json").read_text())
    return dict(line=line, A=z["A"][idx].astype(np.float64), B=z["B"][idx].astype(np.float64),
                real=z["real"][idx].astype(np.float64), ok=(z["avail_A"][idx] & z["avail_B"][idx]),
                avail=z["avail"][idx].astype(bool), ncells_full=np.asarray(z["ncells_full"][idx], dtype=np.int64),
                ncells_A=np.asarray(z["ncells_A"][idx], dtype=np.int64), ncells_B=np.asarray(z["ncells_B"][idx], dtype=np.int64),
                raw_median_real=float(z["raw_median_real"]), raw_median_A=float(z["raw_median_A"]), raw_median_B=float(z["raw_median_B"]),
                meta=meta, sha256=sha256_file(E271 / f"work/half/rows-{line}.npz"))


def half_differences(H: dict) -> np.ndarray:
    """d_t = HALF_FACTOR * (A_t - B_t) / 2 on the rows both halves carry (production units, official axis)."""
    return HALF_FACTOR * 0.5 * (H["A"][H["ok"]] - H["B"][H["ok"]])


def half_factor_check(H: dict, wcols: np.ndarray) -> dict:
    """mean ||d||^2 against (mean ||A||^2 + mean ||B||^2)/2 - mean ||real||^2 on the operator columns: both are the
    full row's noise energy if real = (A + B)/2 with independent halves; the ratio is the empirical half factor^2."""
    ok = H["ok"]
    A, B, R = H["A"][ok][:, wcols], H["B"][ok][:, wcols], H["real"][ok][:, wcols]
    d = 0.5 * (A - B)
    e_d = float((d * d).sum(1).mean())
    e_half = float(0.5 * ((A * A).sum(1).mean() + (B * B).sum(1).mean()))
    e_real = float((R * R).sum(1).mean())
    hm = 0.5 * (A + B)
    return dict(n=int(ok.sum()), mean_sq_d=e_d, mean_sq_half=e_half, mean_sq_real=e_real, implied_noise_energy=e_half - e_real,
                ratio_d_over_implied=e_d / (e_half - e_real) if e_half > e_real else float("nan"),
                median_norm_ratio_half_over_real=float(np.median(np.linalg.norm(A, axis=1) / np.linalg.norm(R, axis=1))),
                median_cos_halfmean_real=float(np.median(row_cos(hm, R))),
                median_norm_ratio_halfmean_over_real=float(np.median(np.linalg.norm(hm, axis=1) / np.linalg.norm(R, axis=1))),
                signal_share_of_real_energy=float(1.0 - e_d / e_real) if e_real > 0 else float("nan"))


# ------------------------------------------------------------------------------------------ A: the noise spectrum
def basis_unit_factor(raw_medians: dict[str, float]) -> tuple[float, dict[str, float]]:
    """062 / build_union_whitener: ref = median over donors of raw_median; basis row = production row * ref /
    EFFECT_SCALE (donor-independent, since renorm_j * raw_median_j = ref)."""
    ref = float(np.median(list(raw_medians.values())))
    return ref / EFFECT_SCALE, {n: ref / v for n, v in raw_medians.items()}


def noise_spectrum(evals: np.ndarray, evecs: np.ndarray, D_list: list[np.ndarray], weights: list[float]) -> tuple[np.ndarray, np.ndarray, list[np.ndarray]]:
    """n_i = v_i^T N v_i with N = sum_j w_j mean_t d_jt d_jt^T (D_j in BASIS units on the operator columns);
    s_i = clip((lambda_i - n_i) / lambda_i, 0, 1).  Returns (n, s, per-donor n_j)."""
    p = evals.size
    n = np.zeros(p)
    per = []
    wsum = float(sum(weights))
    for D, wj in zip(D_list, weights):
        Z = np.asarray(D, dtype=np.float64) @ evecs
        nj = (Z * Z).mean(0)
        per.append(nj)
        n += (wj / wsum) * nj
    with np.errstate(divide="ignore", invalid="ignore"):
        s = np.where(evals > 0, (evals - n) / evals, 0.0)
    return n, np.clip(s, 0.0, 1.0), per


def whitened_energy(P_wcols: np.ndarray, evecs: np.ndarray, evals: np.ndarray, gamma: float) -> tuple[np.ndarray, np.ndarray]:
    """Per-direction energy of the rows before (Z^2) and after whitening (Z^2 lambda^-gamma), un-preserved."""
    Z = np.asarray(P_wcols, dtype=np.float64) @ evecs
    raw = (Z * Z).sum(0)
    return raw, raw * evals ** (-gamma)


def shuffled_s(s: np.ndarray, seed: int = SHUF_SEED) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return s[rng.permutation(s.size)]


def mirrored_s(s: np.ndarray, evals: np.ndarray) -> np.ndarray:
    """The s profile reversed along the eigenvalue rank: s_mirror[rank r] = s[rank p-1-r] (ranks by lambda descending)."""
    order = np.argsort(evals)[::-1]
    out = np.empty_like(s)
    out[order] = s[order[::-1]]
    return out


def hard_s(s: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    return (s >= threshold).astype(np.float64)


def spectrum_report(evals: np.ndarray, s: np.ndarray, n: np.ndarray, e_raw: np.ndarray, e_wht: np.ndarray, gamma: float) -> dict:
    """The s_i spectrum by eigenvalue-rank decile (rank 0 = the largest eigenvalue), the noise-dominated set, and the
    whitened energy it carries; the amplification lambda^(-gamma/2) relative to the top eigenvalue."""
    order = np.argsort(evals)[::-1]
    p = evals.size
    amp = (evals / evals.max()) ** (-gamma / 2.0)
    dec = np.minimum((np.arange(p) * 10) // p, 9)
    rows = []
    s_o, n_o, l_o, er_o, ew_o, amp_o = s[order], n[order], evals[order], e_raw[order], e_wht[order], amp[order]
    for k in range(10):
        m = dec == k
        rows.append(dict(decile=k, n=int(m.sum()), lambda_min=float(l_o[m].min()), lambda_max=float(l_o[m].max()),
                         mean_s=float(s_o[m].mean()), median_s=float(np.median(s_o[m])), share_s_lt_half=float((s_o[m] < 0.5).mean()),
                         share_s_zero=float((s_o[m] == 0).mean()), mean_noise_over_lambda=float((n_o[m] / l_o[m]).mean()),
                         median_noise_over_lambda=float(np.median(n_o[m] / l_o[m])),
                         raw_energy_share=float(er_o[m].sum() / e_raw.sum()), whitened_energy_share=float(ew_o[m].sum() / e_wht.sum()),
                         amplification_median=float(np.median(amp_o[m])), amplification_max=float(amp_o[m].max())))
    noise_dom = s < 0.5
    zero = s == 0
    # where along the rank the noise-dominated directions sit
    rank = np.empty(p, dtype=np.int64)
    rank[order] = np.arange(p)
    return dict(p=int(p), gamma=gamma, deciles=rows,
                n_noise_dominated=int(noise_dom.sum()), share_noise_dominated=float(noise_dom.mean()),
                n_s_zero=int(zero.sum()), share_s_zero=float(zero.mean()),
                whitened_energy_share_noise_dominated=float(e_wht[noise_dom].sum() / e_wht.sum()),
                whitened_energy_share_s_zero=float(e_wht[zero].sum() / e_wht.sum()),
                raw_energy_share_noise_dominated=float(e_raw[noise_dom].sum() / e_raw.sum()),
                rank_of_noise_dominated=dict(median=float(np.median(rank[noise_dom])) if noise_dom.any() else float("nan"),
                                             p10=float(np.percentile(rank[noise_dom], 10)) if noise_dom.any() else float("nan"),
                                             p90=float(np.percentile(rank[noise_dom], 90)) if noise_dom.any() else float("nan"),
                                             first=int(rank[noise_dom].min()) if noise_dom.any() else -1),
                mean_s=float(s.mean()), energy_weighted_mean_s=float((s * e_wht).sum() / e_wht.sum()),
                s_at_energy_quantiles={str(q): float(np.percentile(s, q)) for q in (1, 5, 10, 25, 50, 75, 90, 95, 99)},
                total_noise_over_total_lambda=float(n.sum() / evals.sum()),
                n_lambda_floored=int((evals <= evals.max() * 1e-12 * (1 + 1e-9)).sum()))


# ------------------------------------------------------------------------------------------ the operators
class DampedOperator:
    """W' = V diag(lambda^(-gamma/2) * s^kappa) V^T with production's on-axis preserve.  kappa = 0 (or s None) is
    `plolo.Operator.apply` bit for bit: the same expressions in the same order, the extra multiply skipped."""

    def __init__(self, evals: np.ndarray, evecs: np.ndarray, s: np.ndarray | None = None, kappa: float = 0.0, label: str = "incumbent",
                 *, n: np.ndarray | None = None, floor_c: float = 0.0):
        self.evals = np.asarray(evals, dtype=np.float64)
        self.evecs = evecs                      # shared with the loaded operator (no copy)
        self.s = None if s is None else np.asarray(s, dtype=np.float64)
        self.kappa = float(kappa)
        self.n = None if n is None else np.asarray(n, dtype=np.float64)
        self.floor_c = float(floor_c)
        self.label = label
        self.n_cols = int(self.evals.shape[0])

    @classmethod
    def from_operator(cls, op, s=None, kappa=0.0, label="incumbent", *, n=None, floor_c=0.0):
        return cls(op.evals, op.evecs, s, kappa, label, n=n, floor_c=floor_c)

    def weights(self, gamma: float, lam: float) -> np.ndarray:
        if self.floor_c != 0.0 and self.n is not None:
            return (np.maximum(self.evals, self.floor_c * self.n) + lam) ** (-gamma / 2.0)      # A': the noise floor
        w = (self.evals + lam) ** (-gamma / 2.0)
        if self.kappa != 0.0 and self.s is not None:
            w = w * self.s ** self.kappa                                                        # A: Wiener damping
        return w

    @property
    def is_incumbent(self) -> bool:
        return (self.floor_c == 0.0 or self.n is None) and (self.kappa == 0.0 or self.s is None)

    def transform(self, M_wcols: np.ndarray, gamma: float, lam: float) -> np.ndarray:
        """The un-preserved whitened rows on the operator columns."""
        return ((np.asarray(M_wcols, dtype=np.float64) @ self.evecs) * self.weights(gamma, lam)) @ self.evecs.T

    def apply(self, M, wcols, gamma, lam, *, check_axis: bool = True):
        out = np.asarray(M, dtype=np.float64).copy()
        w = np.asarray(wcols, dtype=np.int64)
        if w.shape[0] != self.n_cols:
            raise ValueError(f"{self.label}: applied on {w.shape[0]} columns, fitted on {self.n_cols}")
        norms_in = np.linalg.norm(out[:, w], axis=1)
        before = float(np.median(norms_in))
        if not before > 0.0:
            raise ValueError(f"{self.label}: on-axis preserve undefined (median row norm {before})")
        wt = self.weights(gamma, lam)
        t = ((out[:, w] @ self.evecs) * wt) @ self.evecs.T
        after = float(np.median(np.linalg.norm(t, axis=1)))
        if not after > 0.0:
            raise ValueError(f"{self.label}: whitened median row norm {after}")
        out[:, w] = t * (before / after)
        return out


def preserve(M: np.ndarray, T_wcols: np.ndarray, wcols: np.ndarray) -> np.ndarray:
    """Production's on-axis preserve for an arbitrary transformed block: M's rows with the operator columns replaced by
    T_wcols rescaled by ONE scalar so their median norm equals M's median norm on those columns."""
    out = np.asarray(M, dtype=np.float64).copy()
    before = float(np.median(np.linalg.norm(out[:, wcols], axis=1)))
    after = float(np.median(np.linalg.norm(T_wcols, axis=1)))
    assert before > 0 and after > 0
    out[:, wcols] = np.asarray(T_wcols, dtype=np.float64) * (before / after)
    return out


def apply_factors(parts: list[tuple[np.ndarray, DampedOperator]], wcols: np.ndarray, gamma: float, lam: float, ref_rows: np.ndarray) -> np.ndarray:
    """sum_k (rows_k through op_k, preserved to rows_k's own median norm) / len(parts), then preserved to ref_rows'
    median norm on the operator columns.  Used for the cross-fitted arms."""
    acc = None
    for rows, op in parts:
        o = op.apply(rows, wcols, gamma, lam)[:, wcols]
        acc = o if acc is None else acc + o
    acc /= len(parts)
    return preserve(ref_rows, acc, wcols)


# ------------------------------------------------------------------------------------------ B: the refit (062's recipe)
def fit_operator(banks: list, raw_medians: dict[str, float], wcols: np.ndarray, official, out_npz: Path, *, E062, chunk_rows: int = 2048,
                 max_rss_gib: float = 13.0, log=print) -> dict:
    """062 `fit_union_whitener`'s arithmetic, the Gram accumulated in row blocks (no basis file): every row of each
    bank, centred on the bank's column mean over its last-source columns of the axis, scaled by ref / raw_median_bank
    (ref = median over banks); C = B^T B / N; eigh; floor max * 1e-12; eigenvectors float32.  Writes out_npz and its
    json (self-describing: sha256, axis_sha256, wcols_sha256, basis rows per bank)."""
    t0 = time.time()
    p = int(len(wcols))
    ref = float(np.median([raw_medians[b.name] for b in banks]))
    C = np.zeros((p, p), dtype=np.float64)
    total = 0
    rows_by = {}
    for bank in banks:
        src_cols, op_cols = E062.last_source_columns(bank, official, np.asarray(wcols, dtype=np.int64))
        mean = E062.stream_column_mean(bank, src_cols, chunk_rows=chunk_rows, max_rss_gib=max_rss_gib)
        mult = ref / raw_medians[bank.name]
        delta = np.load(bank.delta_path, mmap_mode="r")
        n_rows = bank.shape[0]
        for lo in range(0, n_rows, chunk_rows):
            hi = min(lo + chunk_rows, n_rows)
            blk = np.zeros((hi - lo, p), dtype=np.float64)
            blk[:, op_cols] = (np.asarray(delta[lo:hi, src_cols], dtype=np.float64) - mean) * mult
            C += blk.T @ blk
            del blk
            E062.enforce_memory_limit(max_rss_gib, f"{bank.name} gram")
        del delta
        total += n_rows
        rows_by[bank.name] = int(n_rows)
        log(f"  gram: {bank.name} {n_rows} rows x {len(src_cols)} cols (mult {mult:.6f}), {time.time() - t0:.0f}s, peak {peak_gib():.2f} GiB")
    C /= total
    log(f"  covariance {C.shape} over {total} rows; eigh ...")
    evals, evecs = np.linalg.eigh(C)
    del C
    evals = np.maximum(evals, float(evals.max()) * 1e-12)
    evecs32 = np.asarray(evecs, dtype=np.float32)
    del evecs
    out_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez(out_npz, evals=evals, evecs=evecs32)
    axis_symbols = [str(official[i]) for i in wcols]
    meta = dict(axis="union", donors=[b.name for b in banks], n_cols=p,
                wcols_sha256=hashlib.sha256(np.ascontiguousarray(wcols, dtype=np.int64).tobytes()).hexdigest(),
                axis_sha256=hashlib.sha256("\n".join(axis_symbols).encode("utf-8")).hexdigest(),
                basis_rows_by_source=rows_by, basis_total_rows=int(total), reference_norm=ref,
                raw_medians={b.name: raw_medians[b.name] for b in banks},
                basis_normalization="production equal-source renorm: median source norm / donor norm", gamma=0.65,
                fit="272 op272.fit_operator (062 fit_union_whitener arithmetic, blocked Gram)", seconds=round(time.time() - t0, 1),
                peak_rss_gib=peak_gib(), max_eigenvalue=float(evals.max()), min_eigenvalue=float(evals.min()))
    meta["sha256"] = sha256_file(out_npz)
    out_npz.with_suffix(".json").write_text(json.dumps(meta, indent=1))
    log(f"  saved {out_npz.name} sha256 {meta['sha256'][:12]} max eig {meta['max_eigenvalue']:.6f}, {meta['seconds']}s, peak {meta['peak_rss_gib']:.2f} GiB")
    return meta


def load_fitted(npz: Path, label: str) -> DampedOperator:
    with np.load(npz) as z:
        ev = np.asarray(z["evals"], dtype=np.float64)
        V = np.asarray(z["evecs"], dtype=np.float64)          # float32 on disk -> float64, as plolo.Operator does
    return DampedOperator(ev, V, None, 0.0, label)


# ------------------------------------------------------------------------------------------ 271's half banks, verbatim
class Sums:
    """halfrows271.Sums (267 truth267.Sums), verbatim: 267's exact integer pseudobulk sums per (perturbation, half)."""

    def __init__(self, line: str):
        self.line = line
        d = HALF267 / line
        z = np.load(d / "sums.npz", allow_pickle=True)
        self.keep = np.asarray(z["keep"], dtype=np.int64)
        self.keep_names = [str(x) for x in z["keep_names"]]
        self.pert_names = [str(x) for x in z["pert_names"]]
        self.tot = np.asarray(z["tot"], dtype=np.int64)
        self.ncell = np.asarray(z["ncell"], dtype=np.int64)
        self.ctrl_sum = np.asarray(z["ctrl_sum"], dtype=np.float64)
        self.ctrl_n = np.asarray(z["ctrl_n"], dtype=np.int64)
        self.acc = {h: np.load(d / f"acc_{h}.npy", mmap_mode="r") for h in ("batchA", "batchB")}
        self.k562 = line == "k562_gwps"
        self.on_all_genes = self.k562

    def half(self, side: str):
        A, B = self.acc["batchA"], self.acc["batchB"]
        cs, cn, tot, nc = self.ctrl_sum, self.ctrl_n, self.tot, self.ncell
        if side == "A":
            return (lambda lo, hi: np.asarray(A[lo:hi], dtype=np.float64)), tot[:, 0], nc[:, 0], cs[0], cn[0]
        return (lambda lo, hi: np.asarray(B[lo:hi], dtype=np.float64)), tot[:, 1], nc[:, 1], cs[1], cn[1]


def b_ctrl_of(ctrl_sum_all: np.ndarray) -> np.ndarray:
    return np.log1p(TS * ctrl_sum_all / ctrl_sum_all.sum()).astype(np.float32)


def write_half_bank(sm: Sums, side: str, out: Path) -> dict:
    """halfrows271.write_half_bank (truth267.write_half_bank for the batch split, min 5 cells in the half), verbatim."""
    S, tot, nc, cs, cn = sm.half(side)
    b_all = b_ctrl_of(cs)
    keep = sm.keep
    b_keep = b_all[keep]
    bank_perts = set((ROOT / f"data/donors/{sm.line}/perts.txt").read_text().split("\n"))
    rows = [i for i, nm in enumerate(sm.pert_names) if nm in bank_perts and nc[i] >= MIN_HALF and tot[i] > 0]
    rows = sorted(rows, key=lambda i: sm.pert_names[i])
    names = [sm.pert_names[i] for i in rows]
    out.mkdir(parents=True, exist_ok=True)
    delta = np.lib.format.open_memmap(out / "delta.npy", mode="w+", dtype=np.float32, shape=(len(rows), len(keep)))
    rows_a = np.asarray(rows, dtype=np.int64)
    step = 2048
    for lo in range(0, len(rows), step):
        r = rows_a[lo:lo + step]
        blk = S(0, len(sm.pert_names))[r] if len(sm.pert_names) <= 12000 else S(int(r.min()), int(r.max()) + 1)[r - int(r.min())]
        if sm.on_all_genes:
            blk = blk[:, keep]
        prof = np.log1p(TS * blk / tot[r][:, None]).astype(np.float32)
        delta[lo:lo + len(r)] = prof - b_keep
    delta.flush()
    del delta
    np.save(out / "b_ctrl.npy", b_keep.astype(np.float32))
    (out / "perts.txt").write_text("\n".join(names))
    (out / "genes.txt").write_text("\n".join(sm.keep_names))
    meta = dict(line=sm.line, split="batch", side=side, n_perts=len(names), min_cells=MIN_HALF, n_control_cells=int(cn),
                n_genes=len(keep), cells_per_pert=dict(min=int(nc[rows_a].min()), median=float(np.median(nc[rows_a])), max=int(nc[rows_a].max())))
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    return meta
