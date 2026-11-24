#!/usr/bin/env python3
"""Experiment 108: signed-regulon source-transfer pre-gate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
EXPERIMENT = "108-signed-regulon"

NETWORK = ROOT / "data/external/omnipath-dorothea/dorothea-abc.tsv"
EXPECTED_NETWORK_SHA256 = (
    "37f4f6959e9faf960837e8ad5a456514bcd081d12dd11161205dad1f3bc8f671"
)
SOURCE = ROOT / "data/donors/k562_gwps"
TRUTH = {
    "hct116": ROOT / "data/donors/orion_hct116",
    "hek293t": ROOT / "data/donors/orion_hek293t",
}
PANEL = tuple(
    pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
)
PANEL_SET = set(PANEL)
CONFIDENCE = {"A": 1.0, "B": 0.5, "C": 0.25}
GALLERY_SIZE = 300
MIN_REGULON_TARGETS = 5
RIDGE = 0.01
BETA = 0.10
BLEND = 0.10
N_NULL = 20
N_BOOT = 2_000
CHUNK_ROWS = 256


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(line.strip() for line in path.read_text().splitlines() if line.strip())
    if not values:
        raise ValueError(f"{path}: empty axis")
    if unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate axis")
    return values


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def stable_seed(*parts: object) -> int:
    return int.from_bytes(
        hashlib.sha256("\0".join(map(str, parts)).encode()).digest()[:8], "little"
    )


def atomic_npz(path: Path, **values: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with tmp.open("wb") as handle:
        np.savez_compressed(handle, **values)
    tmp.replace(path)


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def last_lookup(values: Sequence[str]) -> dict[str, int]:
    return {str(value): index for index, value in enumerate(values)}


def frozen_gallery() -> tuple[str, ...]:
    target_sets = [
        set(read_lines(path / "perts.txt"))
        for path in [SOURCE, *TRUTH.values()]
    ]
    pool = set.intersection(*target_sets) - PANEL_SET
    ordered = sorted(
        pool,
        key=lambda target: (
            hashlib.sha256(f"{EXPERIMENT}\0{target}".encode()).digest(),
            target,
        ),
    )
    if len(ordered) < GALLERY_SIZE:
        raise ValueError(f"only {len(ordered)} shared non-panel targets")
    return tuple(ordered[:GALLERY_SIZE])


def response_axis(gallery: Sequence[str]) -> tuple[str, ...]:
    gene_sets = [
        set(read_lines(path / "genes.txt", unique=False))
        for path in [SOURCE, *TRUTH.values()]
    ]
    allowed = set.intersection(*gene_sets) - PANEL_SET - set(gallery)
    official = tuple(
        pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str)
    )
    axis = tuple(gene for gene in official if gene in allowed)
    if len(axis) < 1_000:
        raise ValueError(f"only {len(axis)} response genes")
    return axis


def centered_source(
    path: Path, gallery: Sequence[str], axis: Sequence[str]
) -> tuple[np.ndarray, float]:
    genes = read_lines(path / "genes.txt", unique=False)
    targets = read_lines(path / "perts.txt")
    gene_lookup = last_lookup(genes)
    target_lookup = last_lookup(targets)
    cols = np.asarray([gene_lookup[gene] for gene in axis], dtype=np.int64)
    fit_rows = np.asarray(
        [index for index, target in enumerate(targets) if target not in PANEL_SET],
        dtype=np.int64,
    )
    gallery_rows = np.asarray(
        [target_lookup[target] for target in gallery], dtype=np.int64
    )
    delta = np.load(path / "delta.npy", mmap_mode="r")
    total = np.zeros(len(axis), dtype=np.float64)
    for start in range(0, len(fit_rows), CHUNK_ROWS):
        rows = fit_rows[start : start + CHUNK_ROWS]
        total += np.asarray(delta[np.ix_(rows, cols)], dtype=np.float64).sum(0)
    mean = total / len(fit_rows)
    norms = np.empty(len(fit_rows), dtype=np.float64)
    offset = 0
    for start in range(0, len(fit_rows), CHUNK_ROWS):
        rows = fit_rows[start : start + CHUNK_ROWS]
        block = np.asarray(delta[np.ix_(rows, cols)], dtype=np.float64)
        block -= mean
        norms[offset : offset + len(rows)] = np.linalg.norm(block, axis=1)
        offset += len(rows)
    scale = float(np.median(norms)) or 1.0
    selected = np.asarray(
        delta[np.ix_(gallery_rows, cols)], dtype=np.float64
    )
    selected = (selected - mean) / scale
    return selected.astype(np.float32), scale


def parse_bool(series: pd.Series) -> np.ndarray:
    return series.astype(str).str.lower().map({"true": True, "false": False}).fillna(
        False
    ).to_numpy(dtype=bool)


def build_regulon(
    network_path: Path, axis: Sequence[str]
) -> tuple[np.ndarray, tuple[str, ...], np.ndarray, dict[str, int]]:
    if sha256_file(network_path) != EXPECTED_NETWORK_SHA256:
        raise ValueError("DoRothEA export hash changed")
    table = pd.read_csv(network_path, sep="\t")
    stimulation = parse_bool(table["consensus_stimulation"])
    inhibition = parse_bool(table["consensus_inhibition"])
    signed = stimulation ^ inhibition
    table = table.loc[signed].copy()
    signs = np.where(stimulation[signed], 1.0, -1.0)
    levels = table["dorothea_level"].astype(str).str.split(";").str[0]
    weights = levels.map(CONFIDENCE).fillna(0.0).to_numpy(dtype=np.float64)
    axis_lookup = {gene: index for index, gene in enumerate(axis)}

    collapsed: dict[tuple[str, int], float] = {}
    for tf, target, sign, weight in zip(
        table["source_genesymbol"].astype(str),
        table["target_genesymbol"].astype(str),
        signs,
        weights,
        strict=True,
    ):
        column = axis_lookup.get(target)
        if column is None or weight <= 0:
            continue
        key = (tf, column)
        collapsed[key] = collapsed.get(key, 0.0) + float(sign * weight)

    by_tf: dict[str, list[tuple[int, float]]] = {}
    for (tf, column), weight in collapsed.items():
        if weight != 0:
            by_tf.setdefault(tf, []).append((column, weight))
    tf_names = tuple(
        sorted(
            tf
            for tf, edges in by_tf.items()
            if len({column for column, _ in edges}) >= MIN_REGULON_TARGETS
        )
    )
    matrix = np.zeros((len(tf_names), len(axis)), dtype=np.float64)
    for row, tf in enumerate(tf_names):
        for column, weight in by_tf[tf]:
            matrix[row, column] = weight
    row_norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    matrix /= np.maximum(row_norms, 1e-300)
    support = np.any(matrix != 0, axis=0)
    audit = {
        "input_rows": int(len(pd.read_csv(network_path, sep="\t"))),
        "signed_rows": int(signed.sum()),
        "collapsed_edges": int(sum(np.count_nonzero(row) for row in matrix)),
        "transcription_factors": len(tf_names),
        "supported_response_genes": int(support.sum()),
    }
    return matrix.astype(np.float32), tf_names, support, audit


def project_regulon(
    effects: np.ndarray, regulon: np.ndarray, ridge: float = RIDGE
) -> np.ndarray:
    x = np.asarray(effects, dtype=np.float64)
    r = np.asarray(regulon, dtype=np.float64)
    gram = r @ r.T
    gram.flat[:: len(gram) + 1] += ridge
    activities = np.linalg.solve(gram, r @ x.T).T
    return (activities @ r).astype(np.float32)


def sign_preserving_reweight(
    base: np.ndarray,
    reconstruction: np.ndarray,
    support: np.ndarray,
    beta: float = BETA,
) -> np.ndarray:
    base = np.asarray(base, dtype=np.float64)
    reconstruction = np.asarray(reconstruction, dtype=np.float64)
    support = np.asarray(support, dtype=bool)
    if base.shape != reconstruction.shape:
        raise ValueError("base and reconstruction shapes differ")
    out = base.copy()
    eps = 1e-12
    for row in range(len(out)):
        active = support & (base[row] != 0) & (reconstruction[row] != 0)
        before = float(np.linalg.norm(base[row]))
        if int(active.sum()) < 2 or before <= 0:
            continue
        ratio = np.minimum(
            np.abs(reconstruction[row, active])
            / (np.abs(base[row, active]) + eps),
            4.0,
        )
        score = np.sign(base[row, active] * reconstruction[row, active]) * ratio
        median = float(np.median(score))
        mad = float(np.median(np.abs(score - median)))
        z = np.clip((score - median) / max(1.4826 * mad, 1e-6), -2.0, 2.0)
        out[row, active] *= np.exp(beta * z)
        after = float(np.linalg.norm(out[row]))
        if after > 0:
            out[row] *= before / after
    if np.any(np.sign(out[base != 0]) != np.sign(base[base != 0])):
        raise AssertionError("sign-preserving regulon transform flipped a sign")
    return out.astype(np.float32)


def direction_blend(
    base: np.ndarray, reconstruction: np.ndarray, weight: float = BLEND
) -> np.ndarray:
    base = np.asarray(base, dtype=np.float64)
    reconstruction = np.asarray(reconstruction, dtype=np.float64)
    out = base.copy()
    for row in range(len(out)):
        base_norm = float(np.linalg.norm(base[row]))
        recon_norm = float(np.linalg.norm(reconstruction[row]))
        if base_norm <= 0 or recon_norm <= 0:
            continue
        matched = reconstruction[row] * (base_norm / recon_norm)
        out[row] = (1.0 - weight) * base[row] + weight * matched
        mixed_norm = float(np.linalg.norm(out[row]))
        if mixed_norm > 0:
            out[row] *= base_norm / mixed_norm
    return out.astype(np.float32)


def freeze() -> None:
    if WORK.exists():
        raise FileExistsError(f"refusing to overwrite {WORK}")
    gallery = frozen_gallery()
    axis = response_axis(gallery)
    base, scale = centered_source(SOURCE, gallery, axis)
    regulon, tf_names, support, regulon_audit = build_regulon(NETWORK, axis)
    reconstruction = project_regulon(base, regulon)
    reweight = sign_preserving_reweight(base, reconstruction, support)
    blend10 = direction_blend(base, reconstruction)
    rngs = [
        np.random.default_rng(stable_seed(EXPERIMENT, "null", draw))
        for draw in range(N_NULL)
    ]
    support_cols = np.flatnonzero(support)
    null_permutations = np.stack(
        [rng.permutation(support_cols) for rng in rngs]
    ).astype(np.int64)

    WORK.mkdir(parents=True)
    frozen_path = WORK / "frozen.npz"
    atomic_npz(
        frozen_path,
        genes=np.asarray(axis, dtype="U"),
        targets=np.asarray(gallery, dtype="U"),
        base=base,
        reconstruction=reconstruction,
        reweight=reweight,
        blend10=blend10,
        regulon_support=support,
        null_permutations=null_permutations,
    )
    manifest = {
        "experiment": EXPERIMENT,
        "sealed_before_held_truth": True,
        "held_truth_values_opened": False,
        "gallery_size": len(gallery),
        "gallery_sha256": sha256_lines(gallery),
        "response_genes": len(axis),
        "response_axis_sha256": sha256_lines(axis),
        "source_scale": scale,
        "network": {
            "path": str(NETWORK.relative_to(ROOT)),
            "sha256": sha256_file(NETWORK),
            "confidence_weights": CONFIDENCE,
            "minimum_targets_per_tf": MIN_REGULON_TARGETS,
            **regulon_audit,
        },
        "ridge": RIDGE,
        "reweight_beta": BETA,
        "direction_blend": BLEND,
        "transcription_factors_sha256": sha256_lines(tf_names),
        "source_hashes": {
            str((SOURCE / name).relative_to(ROOT)): sha256_file(SOURCE / name)
            for name in ("delta.npy", "genes.txt", "perts.txt", "meta.json")
        },
        "held_metadata_hashes": {
            str((path / name).relative_to(ROOT)): sha256_file(path / name)
            for path in TRUTH.values()
            for name in ("genes.txt", "perts.txt", "meta.json")
        },
        "frozen_sha256": sha256_file(frozen_path),
        "implementation_sha256": sha256_file(Path(__file__)),
        "readme_sha256": sha256_file(HERE / "README.md"),
    }
    atomic_json(WORK / "manifest.json", manifest)
    print(WORK / "manifest.json")


def load_truth(
    path: Path, targets: Sequence[str], genes: Sequence[str]
) -> np.ndarray:
    source_genes = read_lines(path / "genes.txt", unique=False)
    source_targets = read_lines(path / "perts.txt")
    gene_lookup = last_lookup(source_genes)
    target_lookup = last_lookup(source_targets)
    cols = np.asarray([gene_lookup[gene] for gene in genes], dtype=np.int64)
    rows = np.asarray([target_lookup[target] for target in targets], dtype=np.int64)
    delta = np.load(path / "delta.npy", mmap_mode="r")
    return np.asarray(delta[np.ix_(rows, cols)], dtype=np.float64)


def pds_per_target(predicted: np.ndarray, real: np.ndarray) -> np.ndarray:
    predicted = np.asarray(predicted, dtype=np.float64)
    real = np.asarray(real, dtype=np.float64)
    pn = np.linalg.norm(predicted, axis=1)
    rn = np.linalg.norm(real, axis=1)
    scores = (predicted @ real.T) / np.maximum(
        pn[:, None] * rn[None, :], 1e-300
    )
    scores[pn == 0] = -1
    scores[:, rn == 0] = -1
    own = np.diag(scores)[:, None]
    ranks = (scores > own).sum(axis=1)
    ties = (scores == own).sum(axis=1) - 1
    return 1.0 - (ranks + 0.5 * ties) / (len(predicted) - 1)


def bootstrap(values: np.ndarray, label: str) -> dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(stable_seed(EXPERIMENT, label))
    draws = np.empty(N_BOOT, dtype=np.float64)
    for index in range(N_BOOT):
        draws[index] = values[rng.integers(0, len(values), len(values))].mean()
    return {
        "estimate": float(values.mean()),
        "lo": float(np.percentile(draws, 2.5)),
        "hi": float(np.percentile(draws, 97.5)),
    }


def score() -> None:
    manifest = json.loads((WORK / "manifest.json").read_text())
    frozen_path = WORK / "frozen.npz"
    if manifest.get("sealed_before_held_truth") is not True:
        raise ValueError("candidate was not frozen before held truth")
    if sha256_file(frozen_path) != manifest.get("frozen_sha256"):
        raise ValueError("frozen candidate hash changed")
    output_path = WORK / "results.json"
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    with np.load(frozen_path, allow_pickle=False) as z:
        genes = z["genes"].astype(str)
        targets = z["targets"].astype(str)
        base = z["base"].astype(np.float64)
        reconstruction = z["reconstruction"].astype(np.float64)
        support = z["regulon_support"].astype(bool)
        candidates = {
            "base": base,
            "regulon_only": reconstruction,
            "reweight": z["reweight"].astype(np.float64),
            "blend10": z["blend10"].astype(np.float64),
        }
        null_permutations = z["null_permutations"].astype(np.int64)

    surface_results = {}
    for surface, path in TRUTH.items():
        truth = load_truth(path, targets, genes)
        pds = {
            name: pds_per_target(candidate, truth)
            for name, candidate in candidates.items()
        }
        delta = pds["reweight"] - pds["base"]
        null_means = []
        support_cols = np.flatnonzero(support)
        for permutation in null_permutations:
            null_reconstruction = reconstruction.copy()
            null_reconstruction[:, support_cols] = reconstruction[:, permutation]
            null_candidate = sign_preserving_reweight(
                base, null_reconstruction, support
            )
            null_means.append(
                float(
                    np.mean(
                        pds_per_target(null_candidate, truth) - pds["base"]
                    )
                )
            )
        surface_results[surface] = {
            "means": {
                name: float(values.mean()) for name, values in pds.items()
            },
            "reweight_minus_base": bootstrap(
                delta, f"{surface}:reweight-minus-base"
            ),
            "blend10_minus_base": bootstrap(
                pds["blend10"] - pds["base"],
                f"{surface}:blend10-minus-base",
            ),
            "regulon_only_minus_base": bootstrap(
                pds["regulon_only"] - pds["base"],
                f"{surface}:regulon-only-minus-base",
            ),
            "null": {
                "values": null_means,
                "max": float(np.max(null_means)),
                "empirical_p": float(
                    (1 + np.sum(np.asarray(null_means) >= delta.mean()))
                    / (1 + len(null_means))
                ),
            },
            "truth_hashes": {
                str((path / name).relative_to(ROOT)): sha256_file(path / name)
                for name in ("delta.npy", "genes.txt", "perts.txt", "meta.json")
            },
        }

    macro = float(
        np.mean(
            [
                value["reweight_minus_base"]["estimate"]
                for value in surface_results.values()
            ]
        )
    )
    gates = {
        "nonnegative_each_surface": all(
            value["reweight_minus_base"]["estimate"] >= 0
            for value in surface_results.values()
        ),
        "macro_at_least_0_003": macro >= 0.003,
        "bootstrap_positive_each_surface": all(
            value["reweight_minus_base"]["lo"] > 0
            for value in surface_results.values()
        ),
        "beats_null_each_surface": all(
            value["reweight_minus_base"]["estimate"] > value["null"]["max"]
            for value in surface_results.values()
        ),
    }
    result = {
        "experiment": EXPERIMENT,
        "truth_read_after_freeze": True,
        "surfaces": surface_results,
        "macro_reweight_minus_base": macro,
        "gates": gates,
        "verdict": "PROMOTE" if all(gates.values()) else "NO-GO",
        "submission_legal_status": (
            "mechanism evidence only; a supervised rules-safe model is required"
        ),
    }
    atomic_json(output_path, result)
    print(json.dumps(result, indent=2, sort_keys=True))


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("command", choices=("freeze", "score"))
    return p


def main() -> None:
    args = parser().parse_args()
    if args.command == "freeze":
        freeze()
    else:
        score()


if __name__ == "__main__":
    main()
