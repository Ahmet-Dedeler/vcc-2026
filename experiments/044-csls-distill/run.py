"""Nested donor-held-out distillation of CSLS into submit-able query vectors.

The hidden reference gallery makes exact CSLS impossible at challenge inference.
This prototype learns low-capacity query-side approximations from public donor
screens, with both the outer donor and all VCC-panel perturbation identities held
out from fitting and hyperparameter selection.

The implementation is deliberately memory-bounded:

* donor ``delta.npy`` files are memory-mapped;
* statistics are accumulated in row chunks;
* only 600 training and 259 evaluation effects are cached;
* randomized range finding replaces a full gene covariance eigendecomposition;
* scoring operates on small cosine matrices, never on single-cell counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

# Avoid accidental BLAS thread storms on a laptop. These need to be set before
# NumPy loads the BLAS runtime.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "2")
os.environ.setdefault("OMP_NUM_THREADS", "2")

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
CACHE = HERE / "cache"

DONORS = {
    "k562_gwps": ROOT / "data/donors/k562_gwps",
    "orion_hct116": ROOT / "data/donors/orion_hct116",
    "orion_hek293t": ROOT / "data/donors/orion_hek293t",
    "rpe1": ROOT / "data/donors/rpe1",
    "nadig_jurkat": ROOT / "data/donors/nadig_jurkat",
}
OUTER_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
CACHE_VERSION = 2


def rss_gib() -> float:
    """Peak RSS on macOS, in GiB."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def read_lines(path: Path) -> list[str]:
    return [x for x in path.read_text().splitlines() if x]


def row_normalize(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    return x / (np.linalg.norm(x, axis=1, keepdims=True) + 1e-12)


def sharpen(x: np.ndarray, exponent: float) -> np.ndarray:
    if exponent == 1.0:
        return x
    return np.sign(x) * np.abs(x) ** exponent


def pds_from_scores(scores: np.ndarray) -> float:
    """PDS from a similarity matrix whose diagonal is the correct match."""
    own = np.diag(scores)[:, None]
    n = scores.shape[1]
    rank = (scores > own).sum(1) + 0.5 * ((scores == own).sum(1) - 1)
    return float(np.mean(1.0 - rank / (n - 1)))


def candidate_density(scores: np.ndarray, k: int) -> np.ndarray:
    """Average top-k query similarities for each candidate/gallery column."""
    k = min(k, scores.shape[0])
    return np.partition(scores, -k, axis=0)[-k:, :].mean(0)


def exact_csls_pds(q: np.ndarray, t: np.ndarray, k: int) -> float:
    """Exact oracle CSLS. This is diagnostic only and is never a fitted feature."""
    scores = q @ t.T
    density = candidate_density(scores, k)
    # 2*S-r_q-r_t has the same per-query ranking as S-0.5*r_t.
    return pds_from_scores(scores - 0.5 * density[None, :])


def stable_seed(names: Iterable[str], extra: int = 0) -> int:
    text = "|".join(sorted(names)) + f"|{extra}"
    return int.from_bytes(hashlib.sha256(text.encode()).digest()[:4], "little")


@dataclass(frozen=True)
class Candidate:
    method: str
    k: int
    rank: int
    ridge: float
    eta: float

    @property
    def key(self) -> str:
        return (
            f"{self.method}:k={self.k}:r={self.rank}:"
            f"ridge={self.ridge:g}:eta={self.eta:g}"
        )


@dataclass
class Bank:
    raw_train: np.ndarray
    src_train: np.ndarray
    raw_outer: np.ndarray | None
    src_outer: np.ndarray | None


@dataclass
class Fitted:
    basis: np.ndarray
    bias: np.ndarray | None = None
    bands: np.ndarray | None = None


def make_cache(
    train_count: int,
    seed: int,
    chunk_rows: int,
) -> tuple[list[str], list[str], list[str], np.ndarray]:
    """Create compact centered/raw and gene-normalized effect caches."""
    CACHE.mkdir(parents=True, exist_ok=True)
    gene_lists = {n: read_lines(p / "genes.txt") for n, p in DONORS.items()}
    pert_lists = {n: read_lines(p / "perts.txt") for n, p in DONORS.items()}
    gene_sets = {n: set(v) for n, v in gene_lists.items()}
    pert_sets = {n: set(v) for n, v in pert_lists.items()}

    official = list(
        pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str)
    )
    vcc_panel = list(
        pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
    )
    vcc_set = set(vcc_panel)

    # The outer scorer axis mirrors the established 3-line LODO and excludes all
    # VCC target genes exactly as pds_cosine does.
    master_genes = [
        g for g in official
        if g not in vcc_set and all(g in gene_sets[n] for n in OUTER_DONORS)
    ]

    common_train = sorted(set.intersection(*(pert_sets[n] for n in DONORS)))
    if set(common_train) & vcc_set:
        raise RuntimeError("training perturbations unexpectedly overlap the VCC panel")
    if len(common_train) < train_count:
        raise RuntimeError(
            f"only {len(common_train)} common training perturbations, need {train_count}"
        )
    rng = np.random.default_rng(seed)
    train_perts = sorted(rng.choice(common_train, train_count, replace=False).tolist())
    outer_perts = [
        p for p in vcc_panel if all(p in pert_sets[n] for n in OUTER_DONORS)
    ]

    # The learned correction lives only where all five training donors have data
    # and away from every target gene in either gallery. It is embedded as zero
    # on the remaining outer-score coordinates.
    train_target_set = set(train_perts)
    train_mask = np.asarray([
        g not in train_target_set and all(g in gene_sets[n] for n in DONORS)
        for g in master_genes
    ])

    signature = {
        "version": CACHE_VERSION,
        "train_count": train_count,
        "seed": seed,
        "donors": list(DONORS),
        "master_genes_sha256": hashlib.sha256(
            "\n".join(master_genes).encode()
        ).hexdigest(),
        "train_perts_sha256": hashlib.sha256(
            "\n".join(train_perts).encode()
        ).hexdigest(),
        "outer_perts_sha256": hashlib.sha256(
            "\n".join(outer_perts).encode()
        ).hexdigest(),
    }
    meta_path = CACHE / "meta.json"
    expected = [
        CACHE / f"{name}_{kind}.npy"
        for name in DONORS
        for kind in ("raw_train", "src_train")
    ] + [
        CACHE / f"{name}_{kind}.npy"
        for name in OUTER_DONORS
        for kind in ("raw_outer", "src_outer")
    ]
    if (
        meta_path.exists()
        and json.loads(meta_path.read_text()) == signature
        and all(p.exists() for p in expected)
    ):
        print(
            f"cache hit: {len(train_perts)} train perts, {len(outer_perts)} outer "
            f"perts, {len(master_genes)} genes, {int(train_mask.sum())} correction genes",
            flush=True,
        )
        return master_genes, train_perts, outer_perts, train_mask

    print(
        f"building cache: {len(train_perts)} disjoint train perts, "
        f"{len(outer_perts)} VCC perts, {len(master_genes)} outer genes, "
        f"{int(train_mask.sum())} correction genes",
        flush=True,
    )
    mpos = {g: i for i, g in enumerate(master_genes)}
    for name, path in DONORS.items():
        genes = gene_lists[name]
        gpos = {g: i for i, g in enumerate(genes)}
        available = [g for g in master_genes if g in gpos]
        dst = np.asarray([mpos[g] for g in available], dtype=np.int64)
        src = np.asarray([gpos[g] for g in available], dtype=np.int64)
        perts = pert_lists[name]
        ppos = {p: i for i, p in enumerate(perts)}
        delta = np.load(path / "delta.npy", mmap_mode="r")

        # Center exactly over the donor's full perturbation catalogue while
        # touching only the experiment's feature columns.
        mean = np.zeros(len(src), dtype=np.float64)
        for lo in range(0, delta.shape[0], chunk_rows):
            hi = min(delta.shape[0], lo + chunk_rows)
            mean += np.asarray(delta[lo:hi, src], dtype=np.float64).sum(0)
        mean /= delta.shape[0]

        mean_abs = np.zeros(len(src), dtype=np.float64)
        for lo in range(0, delta.shape[0], chunk_rows):
            hi = min(delta.shape[0], lo + chunk_rows)
            block = np.asarray(delta[lo:hi, src], dtype=np.float64)
            mean_abs += np.abs(block - mean).sum(0)
        mean_abs /= delta.shape[0]
        live = mean_abs[mean_abs > 0]
        floor = float(np.percentile(live, 25.0)) if live.size else 1e-9
        scale = np.maximum(mean_abs, floor)

        def extract(
            rows: list[int],
            suffix: str,
            source_delta: np.ndarray = delta,
        ) -> None:
            raw_path = CACHE / f"{name}_raw_{suffix}.npy"
            src_path = CACHE / f"{name}_src_{suffix}.npy"
            raw = np.lib.format.open_memmap(
                raw_path, mode="w+", dtype=np.float32,
                shape=(len(rows), len(master_genes)),
            )
            norm = np.lib.format.open_memmap(
                src_path, mode="w+", dtype=np.float32,
                shape=(len(rows), len(master_genes)),
            )
            raw[:] = 0
            norm[:] = 0
            for lo in range(0, len(rows), chunk_rows):
                hi = min(len(rows), lo + chunk_rows)
                rr = np.asarray(rows[lo:hi], dtype=np.int64)
                centered = np.asarray(
                    source_delta[np.ix_(rr, src)], dtype=np.float64
                ) - mean
                raw[lo:hi, dst] = centered.astype(np.float32)
                norm[lo:hi, dst] = (centered / scale).astype(np.float32)
            raw.flush()
            norm.flush()
            del raw, norm

        extract([ppos[p] for p in train_perts], "train")
        if name in OUTER_DONORS:
            extract([ppos[p] for p in outer_perts], "outer")
        print(
            f"  cached {name:16} features={len(src):5d} peak_rss={rss_gib():.2f} GiB",
            flush=True,
        )
        del delta

    meta_path.write_text(json.dumps(signature, indent=2))
    (CACHE / "master_genes.txt").write_text("\n".join(master_genes) + "\n")
    (CACHE / "train_perts.txt").write_text("\n".join(train_perts) + "\n")
    (CACHE / "outer_perts.txt").write_text("\n".join(outer_perts) + "\n")
    np.save(CACHE / "train_mask.npy", train_mask)
    return master_genes, train_perts, outer_perts, train_mask


def load_banks() -> dict[str, Bank]:
    banks = {}
    for name in DONORS:
        banks[name] = Bank(
            raw_train=np.load(CACHE / f"{name}_raw_train.npy", mmap_mode="r"),
            src_train=np.load(CACHE / f"{name}_src_train.npy", mmap_mode="r"),
            raw_outer=(
                np.load(CACHE / f"{name}_raw_outer.npy", mmap_mode="r")
                if name in OUTER_DONORS else None
            ),
            src_outer=(
                np.load(CACHE / f"{name}_src_outer.npy", mmap_mode="r")
                if name in OUTER_DONORS else None
            ),
        )
    return banks


def source_norms(
    banks: dict[str, Bank], names: Iterable[str], feat: np.ndarray
) -> dict[str, float]:
    out = {}
    for name in names:
        arr = np.asarray(banks[name].src_train[:, feat], dtype=np.float64)
        out[name] = float(np.median(np.linalg.norm(arr, axis=1))) or 1.0
    return out


def transfer_pair(
    banks: dict[str, Bank],
    target: str,
    sources: tuple[str, ...],
    panel: str,
    score_feat: np.ndarray,
    norm_feat: np.ndarray,
    sharpen_q: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return unit query/truth matrices on ``score_feat``."""
    norms = source_norms(banks, sources, norm_feat)
    ref = float(np.median(list(norms.values())))
    src_attr = "src_train" if panel == "train" else "src_outer"
    raw_attr = "raw_train" if panel == "train" else "raw_outer"
    parts = []
    for name in sources:
        arr = getattr(banks[name], src_attr)
        if arr is None:
            raise ValueError(f"{name} has no {panel} panel")
        parts.append((ref / norms[name]) * np.asarray(arr[:, score_feat], dtype=np.float64))
    q = sum(parts) / len(parts)
    truth_arr = getattr(banks[target], raw_attr)
    if truth_arr is None:
        raise ValueError(f"{target} has no {panel} truth")
    t = np.asarray(truth_arr[:, score_feat], dtype=np.float64)
    return row_normalize(sharpen(q, sharpen_q)), row_normalize(t)


def randomized_right_basis(
    x: np.ndarray,
    rank: int,
    seed: int,
    oversample: int = 8,
    power_iters: int = 1,
) -> np.ndarray:
    """Top right singular vectors without a full covariance/eigendecomposition."""
    x = np.asarray(x, dtype=np.float32)
    width = min(rank + oversample, min(x.shape))
    rng = np.random.default_rng(seed)
    omega = rng.standard_normal((x.shape[1], width), dtype=np.float32)
    y = x @ omega
    for _ in range(power_iters):
        q, _ = np.linalg.qr(y, mode="reduced")
        y = x @ (x.T @ q)
    q, _ = np.linalg.qr(y, mode="reduced")
    small = q.T @ x
    _, _, vt = np.linalg.svd(small, full_matrices=False)
    return vt[:rank].T.astype(np.float64)


def standardized_ridge(
    gram: np.ndarray,
    rhs: np.ndarray,
    observations: int,
    ridge: float,
) -> np.ndarray:
    """Ridge solve after unit-variance feature scaling."""
    gram = gram / max(observations, 1)
    rhs = rhs / max(observations, 1)
    scale = np.sqrt(np.maximum(np.diag(gram), 1e-12))
    gs = gram / scale[:, None] / scale[None, :]
    hs = rhs / scale
    coef_s = np.linalg.solve(gs + ridge * np.eye(len(gs)), hs)
    return coef_s / scale


def build_candidates(args: argparse.Namespace) -> list[Candidate]:
    return [
        Candidate(method, k, rank, ridge, eta)
        for method in args.methods
        for k in args.ks
        for rank in args.ranks
        for ridge in args.ridges
        for eta in args.etas
    ]


def fit_bundle(
    banks: dict[str, Bank],
    train_names: tuple[str, ...],
    train_feat: np.ndarray,
    sharpen_q: float,
    candidates: list[Candidate],
) -> dict[str, Fitted]:
    """Fit every candidate on pseudo-transfer folds within ``train_names``."""
    pairs = []
    truths = []
    for target in train_names:
        sources = tuple(n for n in train_names if n != target)
        q, t = transfer_pair(
            banks, target, sources, "train", train_feat, train_feat, sharpen_q
        )
        pairs.append((q, t))
        truths.append(t)

    max_rank = max(c.rank for c in candidates)
    basis = randomized_right_basis(
        np.vstack(truths),
        max_rank,
        stable_seed(train_names, extra=max_rank),
    )
    fitted: dict[str, Fitted] = {}

    for k in sorted({c.k for c in candidates}):
        labels = []
        projections = []
        for q, t in pairs:
            s = q @ t.T
            y = -0.5 * candidate_density(s, k)
            y -= y.mean()
            projections.append((q @ basis, t @ basis, y))
            labels.append(y)

        for rank in sorted({c.rank for c in candidates}):
            v = basis[:, :rank]

            # Translation q' = q + b, with b in the donor-trained spectral span.
            zb = np.vstack([zt[:, :rank] for _, zt, _ in projections])
            yb = np.concatenate(labels)
            gb = zb.T @ zb
            hb = zb.T @ yb

            # Diagonal spectral map q' = q + V diag(a) V^T q.
            ga = np.zeros((rank, rank), dtype=np.float64)
            ha = np.zeros(rank, dtype=np.float64)
            n_pair = 0
            # Joint affine + spectral normal equations.
            gj = np.zeros((2 * rank, 2 * rank), dtype=np.float64)
            hj = np.zeros(2 * rank, dtype=np.float64)
            for zq0, zt0, y in projections:
                zq = zq0[:, :rank]
                zt = zt0[:, :rank]
                nq = len(zq)
                gq = zq.T @ zq
                gt = zt.T @ zt
                sumq = zq.sum(0)
                zty = zt.T @ y
                ga += gq * gt
                ha += sumq * zty
                n_pair += nq * len(zt)

                gj[:rank, :rank] += nq * gt
                cross = gt * sumq[None, :]
                gj[:rank, rank:] += cross
                gj[rank:, :rank] += cross.T
                gj[rank:, rank:] += gq * gt
                hj[:rank] += nq * zty
                hj[rank:] += sumq * zty

            for ridge in sorted({c.ridge for c in candidates}):
                bias = standardized_ridge(gb, hb, len(zb), ridge)
                bands = standardized_ridge(ga, ha, n_pair, ridge)
                both = standardized_ridge(gj, hj, n_pair, ridge)
                for method in sorted({c.method for c in candidates}):
                    for eta in sorted({c.eta for c in candidates}):
                        candidate = Candidate(method, k, rank, ridge, eta)
                        if candidate not in candidates:
                            continue
                        if method == "bias":
                            fitted[candidate.key] = Fitted(v, bias=bias)
                        elif method == "bands":
                            fitted[candidate.key] = Fitted(v, bands=bands)
                        elif method == "combo":
                            fitted[candidate.key] = Fitted(
                                v, bias=both[:rank], bands=both[rank:]
                            )
                        else:
                            raise ValueError(f"unknown method {method!r}")
    print(
        f"  fit {','.join(train_names)}: pairs={len(pairs)} rank={max_rank} "
        f"models={len(fitted)} peak_rss={rss_gib():.2f} GiB",
        flush=True,
    )
    return fitted


def score_fitted(
    q_full: np.ndarray,
    t_full: np.ndarray,
    q_corr: np.ndarray,
    t_corr: np.ndarray,
    fit: Fitted,
    eta: float,
) -> float:
    """Score the exact cosine ranking induced by a low-rank query transform.

    Query normalization after the transform is unnecessary here: it divides an
    entire score row by one positive constant and therefore cannot change rank.
    Targets are already unit-normalized on the full score axis.
    """
    score = q_full @ t_full.T
    zq = q_corr @ fit.basis
    zt = t_corr @ fit.basis
    if fit.bias is not None:
        score += eta * (zt @ fit.bias)[None, :]
    if fit.bands is not None:
        score += eta * ((zq * fit.bands) @ zt.T)
    return pds_from_scores(score)


def evaluate_bundle(
    banks: dict[str, Bank],
    bundle: dict[str, Fitted],
    target: str,
    sources: tuple[str, ...],
    panel: str,
    score_feat: np.ndarray,
    corr_feat: np.ndarray,
    norm_feat: np.ndarray,
    sharpen_q: float,
    candidates: list[Candidate],
) -> tuple[float, dict[str, float], dict[int, float]]:
    q, t = transfer_pair(
        banks, target, sources, panel, score_feat, norm_feat, sharpen_q
    )
    if np.array_equal(score_feat, corr_feat):
        qc, tc = q, t
    else:
        pos = {int(g): i for i, g in enumerate(score_feat)}
        local = np.asarray([pos[int(g)] for g in corr_feat], dtype=np.int64)
        # These are slices of already full-axis-normalized vectors. That is the
        # correct coordinate block for an embedded query transform.
        qc, tc = q[:, local], t[:, local]

    scores = {}
    by_key = {c.key: c for c in candidates}
    for key, fit in bundle.items():
        scores[key] = score_fitted(q, t, qc, tc, fit, by_key[key].eta)
    oracle = {
        k: exact_csls_pds(q, t, k) for k in sorted({c.k for c in candidates})
    }
    return pds_from_scores(q @ t.T), scores, oracle


def choose_candidate(
    mean_scores: dict[str, float],
    baseline: float,
    candidates: list[Candidate],
    min_gain: float,
) -> str:
    by_key = {c.key: c for c in candidates}

    def complexity(key: str) -> tuple:
        c = by_key[key]
        # Prefer the smaller, more regularized model inside a tiny PDS tie.
        return (c.rank, -c.ridge, abs(c.eta - 1.0), c.k, c.method)

    best_score = max(mean_scores.values())
    near = [
        key for key, value in mean_scores.items()
        if value >= best_score - 5e-4
    ]
    selected = min(near, key=complexity)
    if mean_scores[selected] < baseline + min_gain:
        return "baseline"
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-count", type=int, default=600)
    parser.add_argument("--seed", type=int, default=44)
    parser.add_argument("--chunk-rows", type=int, default=256)
    parser.add_argument("--sharpen", type=float, default=1.0)
    parser.add_argument("--ks", nargs="+", type=int, default=[5, 10, 20])
    parser.add_argument("--ranks", nargs="+", type=int, default=[8, 16, 32])
    parser.add_argument("--ridges", nargs="+", type=float, default=[0.01, 0.1, 1.0])
    parser.add_argument("--etas", nargs="+", type=float, default=[0.5, 1.0])
    parser.add_argument(
        "--methods", nargs="+", choices=["bias", "bands", "combo"],
        default=["bias", "bands", "combo"],
    )
    parser.add_argument("--min-inner-gain", type=float, default=0.001)
    parser.add_argument("--output", type=Path, default=HERE / "results.json")
    args = parser.parse_args()

    master_genes, train_perts, outer_perts, train_mask = make_cache(
        args.train_count, args.seed, args.chunk_rows
    )
    banks = load_banks()
    all_feat = np.arange(len(master_genes), dtype=np.int64)
    train_feat = np.flatnonzero(train_mask).astype(np.int64)
    candidates = build_candidates(args)
    print(
        f"running nested LODO: {len(candidates)} candidates, "
        f"train={len(train_perts)} disjoint perts, outer={len(outer_perts)} VCC perts, "
        f"sharpen={args.sharpen:g}",
        flush=True,
    )

    bundle_cache: dict[tuple[str, ...], dict[str, Fitted]] = {}

    def bundle(names: Iterable[str]) -> dict[str, Fitted]:
        key = tuple(sorted(names))
        if key not in bundle_cache:
            bundle_cache[key] = fit_bundle(
                banks, key, train_feat, args.sharpen, candidates
            )
        return bundle_cache[key]

    nested = {}
    fixed_scores = {c.key: {} for c in candidates}
    for outer in OUTER_DONORS:
        outer_train = tuple(n for n in DONORS if n != outer)
        print(f"\nouter holdout {outer}; fit/tune donors={outer_train}", flush=True)
        inner_by_candidate = {c.key: [] for c in candidates}
        inner_baselines = []
        inner_detail = {}
        for valid in outer_train:
            fit_names = tuple(n for n in outer_train if n != valid)
            base, scores, oracle = evaluate_bundle(
                banks,
                bundle(fit_names),
                valid,
                fit_names,
                "train",
                train_feat,
                train_feat,
                train_feat,
                args.sharpen,
                candidates,
            )
            inner_baselines.append(base)
            for key, value in scores.items():
                inner_by_candidate[key].append(value)
            inner_detail[valid] = {
                "baseline": base,
                "oracle": oracle,
            }
            print(
                f"  inner valid={valid:16} baseline={base:.4f} "
                f"oracle_best={max(oracle.values()):.4f}",
                flush=True,
            )

        inner_mean = {
            key: float(np.mean(values)) for key, values in inner_by_candidate.items()
        }
        inner_base = float(np.mean(inner_baselines))
        selected = choose_candidate(
            inner_mean, inner_base, candidates, args.min_inner_gain
        )

        final_bundle = bundle(outer_train)
        outer_sources = tuple(n for n in OUTER_DONORS if n != outer)
        base, scores, oracle = evaluate_bundle(
            banks,
            final_bundle,
            outer,
            outer_sources,
            "outer",
            all_feat,
            train_feat,
            train_feat,
            args.sharpen,
            candidates,
        )
        for key, value in scores.items():
            fixed_scores[key][outer] = value
        selected_score = base if selected == "baseline" else scores[selected]
        nested[outer] = {
            "inner_baseline": inner_base,
            "inner_selected": selected,
            "inner_selected_score": (
                inner_base if selected == "baseline" else inner_mean[selected]
            ),
            "inner_best_score": max(inner_mean.values()),
            "outer_baseline": base,
            "outer_selected_score": selected_score,
            "outer_delta": selected_score - base,
            "outer_oracle_by_k": oracle,
            "outer_oracle_best": max(oracle.values()),
            "inner_detail": inner_detail,
        }
        print(
            f"  selected={selected}\n"
            f"  outer baseline={base:.4f} selected={selected_score:.4f} "
            f"delta={selected_score-base:+.4f} "
            f"oracle={max(oracle.values()):.4f}",
            flush=True,
        )

    nested_base = float(np.mean([v["outer_baseline"] for v in nested.values()]))
    nested_score = float(np.mean([v["outer_selected_score"] for v in nested.values()]))
    nested_oracle = float(np.mean([v["outer_oracle_best"] for v in nested.values()]))

    fixed_summary = {}
    for key, folds in fixed_scores.items():
        fixed_summary[key] = {
            **folds,
            "macro": float(np.mean(list(folds.values()))),
            "delta": float(np.mean(list(folds.values()))) - nested_base,
        }
    fixed_order = sorted(
        fixed_summary, key=lambda key: fixed_summary[key]["macro"], reverse=True
    )

    result = {
        "design": {
            "outer_donors": list(OUTER_DONORS),
            "all_training_donors": list(DONORS),
            "train_perturbations": len(train_perts),
            "outer_vcc_perturbations": len(outer_perts),
            "train_outer_overlap": len(set(train_perts) & set(outer_perts)),
            "outer_score_genes": len(master_genes),
            "correction_genes": int(train_mask.sum()),
            "sharpen": args.sharpen,
            "candidate_count": len(candidates),
            "full_eigendecompositions": 0,
        },
        "nested": nested,
        "nested_macro": {
            "baseline": nested_base,
            "selected": nested_score,
            "delta": nested_score - nested_base,
            "oracle_csls": nested_oracle,
            "oracle_headroom": nested_oracle - nested_base,
            "recovered_fraction": (
                (nested_score - nested_base) / (nested_oracle - nested_base)
                if nested_oracle > nested_base else None
            ),
        },
        "fixed_grid_lodo_optimistic": {
            "warning": (
                "Exploratory only: selecting from this table and reporting the same "
                "folds would be optimistic. Use nested_macro for the verdict."
            ),
            "top": [
                {"config": key, **fixed_summary[key]} for key in fixed_order[:20]
            ],
        },
        "peak_rss_gib": rss_gib(),
    }
    args.output.write_text(json.dumps(result, indent=2))

    print("\n== nested verdict ==", flush=True)
    print(
        f"baseline={nested_base:.5f} selected={nested_score:.5f} "
        f"delta={nested_score-nested_base:+.5f} "
        f"oracle_csls={nested_oracle:.5f} "
        f"headroom={nested_oracle-nested_base:+.5f}",
        flush=True,
    )
    print("\nexploratory fixed-grid top 5 (not unbiased):", flush=True)
    for key in fixed_order[:5]:
        row = fixed_summary[key]
        print(f"  {key:52} macro={row['macro']:.5f} delta={row['delta']:+.5f}")
    print(f"results: {args.output}", flush=True)
    print(f"peak RSS: {rss_gib():.2f} GiB", flush=True)


if __name__ == "__main__":
    main()
