"""Build a format-valid VCC 2026 submission .h5ad.

The submission contract is unforgiving and the failure modes are silent (a dense
matrix is 1.40x over the density cap; a log-normalised matrix is simply
unscoreable). This module owns that contract so model code never has to.

Memory note: a submission with ~6,000 stored entries per cell is ~2.2e9 nonzeros,
which is ~17 GB as float32 data + int32 indices. We therefore stream CSR rows
straight into the HDF5 file chunk by chunk instead of materialising the matrix.

Usage
-----
    from build_submission import ChallengeSpec, write_submission

    spec = ChallengeSpec.load("data/challenge")

    def predict(context: str, target_gene: str) -> np.ndarray:
        '''Return a length-18533 vector of expected *relative* expression
        (any non-negative scale; it is renormalised internally).'''
        ...

    write_submission(spec, predict, "submissions/prediction.h5ad")

Then:
    vcc prep submissions/prediction.h5ad -g data/challenge/gene_names.csv \
        --perts data/challenge/pert_counts.csv -o submissions/prediction.vcc --dry-run
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Sequence

import h5py
import numpy as np
import pandas as pd

# --- Contract constants, from the official spec. Do not soften these. -------
N_GENES = 18_533
CELLS_PER_PERT = 400
MAX_STORED_ENTRIES = 4_750_000_000
MAX_COUNTS_PER_CELL = 1_000_000
PERT_COL = "target_gene"
CONTEXT_COL = "context"


@dataclass
class ChallengeSpec:
    """The three reference files that define a round."""

    genes: np.ndarray  # (18533,) gene symbols, canonical order
    perts: np.ndarray  # (300,) target gene symbols
    contexts: list[str]  # ["A", "B", "C"] or ["D", "E", "F"]
    cells_per_pert: int
    manifest: dict

    @classmethod
    def load(cls, bundle_dir: str | Path) -> "ChallengeSpec":
        d = Path(bundle_dir)
        # Both CSVs ship WITH a header row despite `vcc prep --genes` documenting
        # "headerless". Read by column name and never by position.
        genes = pd.read_csv(d / "gene_names.csv")["gene_name"].to_numpy()
        perts = pd.read_csv(d / "pert_counts.csv")[PERT_COL].to_numpy()
        manifest = json.loads((d / "manifest.json").read_text())

        if len(genes) != N_GENES:
            raise ValueError(f"expected {N_GENES} genes, got {len(genes)}")
        if len(genes) != len(set(genes)):
            raise ValueError("duplicate gene symbols in gene_names.csv")

        return cls(
            genes=genes,
            perts=perts,
            contexts=list(manifest["contexts"]),
            cells_per_pert=int(manifest.get("cells_per_pert", CELLS_PER_PERT)),
            manifest=manifest,
        )

    @property
    def n_rows(self) -> int:
        return len(self.perts) * self.cells_per_pert * len(self.contexts)


def sample_counts(
    profile: np.ndarray,
    lib_sizes: np.ndarray,
    rng: np.random.Generator,
    dispersion: float | None = 0.5,
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Turn one expected expression profile into `len(lib_sizes)` count vectors.

    The scorer recomputes differential expression from your counts with a
    Wilcoxon rank-sum test, so the cell-to-cell spread you emit *is* what decides
    significance. Emitting 400 identical cells makes every gene trivially
    significant and destroys `de_wilcoxon_sig_jaccard`. Treat `dispersion` as a
    tuned hyperparameter, not a default.

    Parameters
    ----------
    profile
        Non-negative, length-N_GENES vector of expected relative expression.
    lib_sizes
        Target total counts per cell. Clipped to MAX_COUNTS_PER_CELL.
    dispersion
        Negative-binomial dispersion (variance = mu + dispersion * mu^2).
        `None` gives pure multinomial/Poisson sampling (under-dispersed relative
        to real scRNA-seq).

    Returns a list of (gene_indices, counts) pairs, both int-typed and sorted.
    """
    profile = np.asarray(profile, dtype=np.float64)
    if profile.shape != (N_GENES,):
        raise ValueError(f"profile must be ({N_GENES},), got {profile.shape}")
    if np.any(profile < 0) or not np.all(np.isfinite(profile)):
        raise ValueError("profile must be finite and non-negative")

    total = profile.sum()
    if total <= 0:
        raise ValueError("profile sums to zero")
    p = profile / total

    lib_sizes = np.minimum(np.asarray(lib_sizes), MAX_COUNTS_PER_CELL)
    out = []
    for lib in lib_sizes:
        mu = p * float(lib)
        if dispersion:
            # NB as a gamma-Poisson mixture.
            shape = 1.0 / dispersion
            mu = rng.gamma(shape=shape, scale=mu * dispersion)
        counts = rng.poisson(mu)
        idx = np.flatnonzero(counts)
        vals = counts[idx]
        # Enforce the per-cell cap after sampling; NB tails can overshoot.
        if vals.sum() > MAX_COUNTS_PER_CELL:
            scale = MAX_COUNTS_PER_CELL / vals.sum()
            vals = np.floor(vals * scale).astype(np.int64)
            keep = vals > 0
            idx, vals = idx[keep], vals[keep]
        out.append((idx.astype(np.int32), vals.astype(np.float32)))
    return out


def write_submission(
    spec: ChallengeSpec,
    predict: Callable[[str, str], np.ndarray],
    out_path: str | Path,
    *,
    lib_size_sampler: Callable[[str, np.random.Generator, int], np.ndarray] | None = None,
    dispersion: float | None = 0.5,
    seed: int = 0,
    verbose: bool = True,
) -> Path:
    """Stream a valid submission .h5ad to disk.

    `predict(context, target_gene) -> (18533,) non-negative array` is called once
    per (context, perturbation) pair — 900 times for a full round.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)

    if lib_size_sampler is None:
        # Match the reference's ~20,000 median UMI without touching the cap.
        def lib_size_sampler(context, rng, n):  # noqa: ARG001
            return rng.normal(20_000, 4_000, n).clip(3_000, 60_000).astype(np.int64)

    n_rows = spec.n_rows
    obs_pert: list[str] = []
    obs_ctx: list[str] = []

    with h5py.File(out_path, "w") as f:
        # --- .X as a growable CSR ---------------------------------------
        X = f.create_group("X")
        X.attrs["encoding-type"] = "csr_matrix"
        X.attrs["encoding-version"] = "0.1.0"
        X.attrs["shape"] = np.array([n_rows, N_GENES], dtype=np.int64)
        d_data = X.create_dataset(
            "data", shape=(0,), maxshape=(None,), dtype=np.float32, chunks=(1 << 20,)
        )
        d_idx = X.create_dataset(
            "indices", shape=(0,), maxshape=(None,), dtype=np.int32, chunks=(1 << 20,)
        )
        indptr = np.zeros(n_rows + 1, dtype=np.int64)

        nnz = 0
        row = 0
        for ctx in spec.contexts:
            for pert in spec.perts:
                profile = predict(ctx, pert)
                libs = lib_size_sampler(ctx, rng, spec.cells_per_pert)
                rows = sample_counts(profile, libs, rng, dispersion=dispersion)

                block_idx = np.concatenate([r[0] for r in rows])
                block_val = np.concatenate([r[1] for r in rows])
                lengths = np.fromiter((len(r[0]) for r in rows), dtype=np.int64)

                d_data.resize((nnz + len(block_val),))
                d_idx.resize((nnz + len(block_idx),))
                d_data[nnz:] = block_val
                d_idx[nnz:] = block_idx

                indptr[row + 1 : row + 1 + spec.cells_per_pert] = nnz + np.cumsum(lengths)
                nnz += len(block_val)
                row += spec.cells_per_pert

                obs_pert.extend([pert] * spec.cells_per_pert)
                obs_ctx.extend([ctx] * spec.cells_per_pert)

                if nnz > MAX_STORED_ENTRIES:
                    raise RuntimeError(
                        f"density cap exceeded at row {row}: {nnz:,} > "
                        f"{MAX_STORED_ENTRIES:,}. Emit fewer stored entries per cell."
                    )
            if verbose:
                print(f"  context {ctx}: {row:,} rows, {nnz:,} stored entries")

        X.create_dataset("indptr", data=indptr, dtype=np.int64)

        # --- .obs --------------------------------------------------------
        obs = f.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "_index"
        obs.attrs.create("column-order", [PERT_COL, CONTEXT_COL], dtype=_str_dtype())

        _write_str(obs, "_index", np.array([f"cell_{i}" for i in range(n_rows)], dtype=object))
        _write_categorical(obs, PERT_COL, np.asarray(obs_pert, dtype=object))
        _write_categorical(obs, CONTEXT_COL, np.asarray(obs_ctx, dtype=object))

        # --- .var --------------------------------------------------------
        var = f.create_group("var")
        var.attrs["encoding-type"] = "dataframe"
        var.attrs["encoding-version"] = "0.2.0"
        var.attrs["_index"] = "_index"
        var.attrs.create("column-order", np.empty(0, dtype=object), dtype=_str_dtype())
        _write_str(var, "_index", spec.genes.astype(object))

        for name in ("layers", "obsm", "varm", "obsp", "varp", "uns"):
            g = f.create_group(name)
            g.attrs["encoding-type"] = "dict"
            g.attrs["encoding-version"] = "0.1.0"

        f.attrs["encoding-type"] = "anndata"
        f.attrs["encoding-version"] = "0.1.0"

    if verbose:
        print(
            f"wrote {out_path} — {n_rows:,} cells x {N_GENES:,} genes, "
            f"{nnz:,} stored entries ({nnz / n_rows:.0f}/cell, "
            f"{100 * nnz / MAX_STORED_ENTRIES:.1f}% of cap)"
        )
    return out_path


def _str_dtype():
    return h5py.string_dtype(encoding="utf-8")


def _write_str(group: h5py.Group, name: str, values: np.ndarray) -> None:
    ds = group.create_dataset(name, data=values, dtype=_str_dtype())
    ds.attrs["encoding-type"] = "string-array"
    ds.attrs["encoding-version"] = "0.2.0"


def _write_categorical(group: h5py.Group, name: str, values: np.ndarray) -> None:
    cats, codes = np.unique(values, return_inverse=True)
    g = group.create_group(name)
    g.attrs["encoding-type"] = "categorical"
    g.attrs["encoding-version"] = "0.2.0"
    g.attrs["ordered"] = False
    _write_str(g, "categories", cats.astype(object))
    g.create_dataset("codes", data=codes.astype(np.int32), dtype=np.int32)


def check_submission(path: str | Path, spec: ChallengeSpec) -> None:
    """Re-read a written submission and assert every contract clause.

    `vcc prep --dry-run` is the authority; this is a fast local pre-flight that
    catches mistakes without loading the whole matrix.
    """
    import anndata as ad

    a = ad.read_h5ad(path, backed="r")
    assert a.shape == (spec.n_rows, N_GENES), f"shape {a.shape}"
    assert list(a.var_names) == list(spec.genes), "gene axis mismatch"
    assert set(a.obs[CONTEXT_COL]) == set(spec.contexts), "context labels mismatch"
    assert set(a.obs[PERT_COL]) == set(spec.perts), "perturbation set mismatch"
    assert "non-targeting" not in set(a.obs[PERT_COL]), "control rows present"
    counts = a.obs.groupby([CONTEXT_COL, PERT_COL], observed=True).size()
    assert (counts == spec.cells_per_pert).all(), "wrong cells per perturbation"
    print(f"OK: {path} satisfies the shape/label contract "
          f"({spec.n_rows:,} rows, {len(spec.perts)} perts x {len(spec.contexts)} contexts)")
