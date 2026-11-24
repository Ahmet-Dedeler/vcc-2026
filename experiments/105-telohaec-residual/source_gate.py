#!/usr/bin/env python3
"""Leakage-safe direct-source gate for the repaired TeloHAEC extraction.

Freeze builds TeloHAEC and K562 predictions using target names but never opens
HCT116 or HEK293T effect values. Score opens those two held contexts only after
the source predictions, gallery, masks, and null permutations are sealed.
"""

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
WORK = HERE / "work-source-gate"
EXPERIMENT = "105-telohaec-source-gate"

SOURCE_PATHS = {
    "telohaec": ROOT / "data/external/igvf-telohaec",
    "k562": ROOT / "data/donors/k562_gwps",
}
SOURCE_METADATA = {
    "telohaec": "extraction_report.json",
    "k562": "meta.json",
}
TRUTH_PATHS = {
    "hct116": ROOT / "data/donors/orion_hct116",
    "hek293t": ROOT / "data/donors/orion_hek293t",
}
PANEL = tuple(
    pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str)
)
PANEL_SET = set(PANEL)
GALLERY_SIZE = 300
N_BOOT = 2_000
N_NULL = 20
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
        for path in [*SOURCE_PATHS.values(), *TRUTH_PATHS.values()]
    ]
    pool = set.intersection(*target_sets) - PANEL_SET
    ordered = sorted(
        pool,
        key=lambda gene: (
            hashlib.sha256(f"{EXPERIMENT}\0{gene}".encode()).digest(),
            gene,
        ),
    )
    if len(ordered) < GALLERY_SIZE:
        raise ValueError(f"only {len(ordered)} eligible shared targets")
    return tuple(ordered[:GALLERY_SIZE])


def common_gene_axis(gallery: Sequence[str]) -> tuple[str, ...]:
    gene_sets = [
        set(read_lines(path / "genes.txt", unique=False))
        for path in [*SOURCE_PATHS.values(), *TRUTH_PATHS.values()]
    ]
    allowed = set.intersection(*gene_sets) - PANEL_SET - set(gallery)
    official = tuple(
        pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str)
    )
    axis = tuple(gene for gene in official if gene in allowed)
    if len(axis) < 1_000:
        raise ValueError(f"source gate has only {len(axis)} response genes")
    return axis


def source_prediction(
    path: Path, gallery: Sequence[str], axis: Sequence[str]
) -> tuple[np.ndarray, float]:
    genes = read_lines(path / "genes.txt", unique=False)
    targets = read_lines(path / "perts.txt")
    gene_lookup = last_lookup(genes)
    target_lookup = last_lookup(targets)
    cols = np.asarray([gene_lookup[gene] for gene in axis], dtype=np.int64)
    nonpanel_rows = np.asarray(
        [index for index, target in enumerate(targets) if target not in PANEL_SET],
        dtype=np.int64,
    )
    selected_rows = np.asarray(
        [target_lookup[target] for target in gallery], dtype=np.int64
    )
    delta = np.load(path / "delta.npy", mmap_mode="r")
    if delta.shape != (len(targets), len(genes)):
        raise ValueError(f"{path}: effect matrix shape changed")

    total = np.zeros(len(axis), dtype=np.float64)
    norms = np.empty(len(nonpanel_rows), dtype=np.float64)
    offset = 0
    for start in range(0, len(nonpanel_rows), CHUNK_ROWS):
        rows = nonpanel_rows[start : start + CHUNK_ROWS]
        block = np.asarray(delta[np.ix_(rows, cols)], dtype=np.float64)
        total += block.sum(axis=0)
        norms[offset : offset + len(rows)] = np.linalg.norm(block, axis=1)
        offset += len(rows)
    mean = total / len(nonpanel_rows)

    offset = 0
    for start in range(0, len(nonpanel_rows), CHUNK_ROWS):
        rows = nonpanel_rows[start : start + CHUNK_ROWS]
        block = np.asarray(delta[np.ix_(rows, cols)], dtype=np.float64)
        block -= mean
        norms[offset : offset + len(rows)] = np.linalg.norm(block, axis=1)
        offset += len(rows)
    scale = float(np.median(norms)) or 1.0
    prediction = np.asarray(
        delta[np.ix_(selected_rows, cols)], dtype=np.float64
    )
    prediction -= mean
    prediction /= scale
    return prediction.astype(np.float32), scale


def equal_source_blend(telohaec: np.ndarray, k562: np.ndarray) -> np.ndarray:
    if telohaec.shape != k562.shape:
        raise ValueError("source predictions do not align")
    return (0.5 * telohaec + 0.5 * k562).astype(np.float32)


def freeze() -> None:
    if WORK.exists():
        raise FileExistsError(f"refusing to overwrite {WORK}")
    gallery = frozen_gallery()
    axis = common_gene_axis(gallery)
    telohaec, telohaec_scale = source_prediction(
        SOURCE_PATHS["telohaec"], gallery, axis
    )
    k562, k562_scale = source_prediction(SOURCE_PATHS["k562"], gallery, axis)
    equal = equal_source_blend(telohaec, k562)
    rngs = [
        np.random.default_rng(stable_seed(EXPERIMENT, "null", draw))
        for draw in range(N_NULL)
    ]
    null_permutations = np.stack(
        [rng.permutation(len(gallery)) for rng in rngs]
    ).astype(np.int64)

    WORK.mkdir(parents=True)
    frozen_path = WORK / "frozen.npz"
    atomic_npz(
        frozen_path,
        genes=np.asarray(axis, dtype="U"),
        targets=np.asarray(gallery, dtype="U"),
        telohaec=telohaec,
        k562=k562,
        equal=equal,
        null_permutations=null_permutations,
    )
    manifest = {
        "experiment": EXPERIMENT,
        "stage": "source-freeze",
        "sealed_before_held_truth": True,
        "held_truth_values_opened": False,
        "gallery_size": len(gallery),
        "gallery_sha256": sha256_lines(gallery),
        "response_genes": len(axis),
        "response_axis_sha256": sha256_lines(axis),
        "mask": "VCC-300 union frozen 300-target gallery",
        "source_scales": {
            "telohaec": telohaec_scale,
            "k562": k562_scale,
        },
        "source_hashes": {
            str((path / name).relative_to(ROOT)): sha256_file(path / name)
            for source, path in SOURCE_PATHS.items()
            for name in (
                "delta.npy",
                "genes.txt",
                "perts.txt",
                SOURCE_METADATA[source],
            )
        },
        "held_metadata_hashes": {
            str((path / name).relative_to(ROOT)): sha256_file(path / name)
            for path in TRUTH_PATHS.values()
            for name in ("genes.txt", "perts.txt", "meta.json")
        },
        "frozen_sha256": sha256_file(frozen_path),
        "implementation_sha256": sha256_file(Path(__file__)),
    }
    atomic_json(WORK / "manifest.json", manifest)
    print(WORK / "manifest.json")


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


def cosine_per_target(predicted: np.ndarray, real: np.ndarray) -> np.ndarray:
    predicted = np.asarray(predicted, dtype=np.float64)
    real = np.asarray(real, dtype=np.float64)
    return np.einsum("ij,ij->i", predicted, real) / np.maximum(
        np.linalg.norm(predicted, axis=1) * np.linalg.norm(real, axis=1),
        1e-300,
    )


def load_truth(
    path: Path, gallery: Sequence[str], axis: Sequence[str]
) -> np.ndarray:
    genes = read_lines(path / "genes.txt", unique=False)
    targets = read_lines(path / "perts.txt")
    gene_lookup = last_lookup(genes)
    target_lookup = last_lookup(targets)
    cols = np.asarray([gene_lookup[gene] for gene in axis], dtype=np.int64)
    rows = np.asarray([target_lookup[target] for target in gallery], dtype=np.int64)
    delta = np.load(path / "delta.npy", mmap_mode="r")
    return np.asarray(delta[np.ix_(rows, cols)], dtype=np.float64)


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
        raise ValueError("source predictions were not sealed before held truth")
    if sha256_file(frozen_path) != manifest.get("frozen_sha256"):
        raise ValueError("frozen source bundle hash changed")
    output_path = WORK / "results.json"
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    with np.load(frozen_path, allow_pickle=False) as z:
        genes = z["genes"].astype(str)
        targets = z["targets"].astype(str)
        candidates = {
            "telohaec": z["telohaec"].astype(np.float64),
            "k562": z["k562"].astype(np.float64),
            "equal": z["equal"].astype(np.float64),
        }
        permutations = z["null_permutations"].astype(np.int64)

    surfaces = {}
    for surface, path in TRUTH_PATHS.items():
        truth = load_truth(path, targets, genes)
        pds = {
            name: pds_per_target(value, truth)
            for name, value in candidates.items()
        }
        cosine = {
            name: cosine_per_target(value, truth)
            for name, value in candidates.items()
        }
        equal_delta = pds["equal"] - pds["k562"]
        telo_delta = pds["telohaec"] - pds["k562"]
        null_values = []
        for permutation in permutations:
            shuffled = equal_source_blend(
                candidates["telohaec"][permutation], candidates["k562"]
            )
            null_values.append(
                float(
                    np.mean(
                        pds_per_target(shuffled, truth) - pds["k562"]
                    )
                )
            )
        surfaces[surface] = {
            "means": {
                f"{metric}_{name}": float(np.mean(values))
                for metric, by_name in (("pds", pds), ("cosine", cosine))
                for name, values in by_name.items()
            },
            "telohaec_minus_k562_pds": bootstrap(
                telo_delta, f"{surface}:telohaec-minus-k562"
            ),
            "equal_minus_k562_pds": bootstrap(
                equal_delta, f"{surface}:equal-minus-k562"
            ),
            "equal_null": {
                "values": null_values,
                "max": float(np.max(null_values)),
                "p95_conservative": float(np.max(null_values)),
                "empirical_p": float(
                    (1 + np.sum(np.asarray(null_values) >= equal_delta.mean()))
                    / (1 + len(null_values))
                ),
            },
            "truth_hashes": {
                str((path / name).relative_to(ROOT)): sha256_file(path / name)
                for name in ("delta.npy", "genes.txt", "perts.txt", "meta.json")
            },
        }

    macro_equal = float(
        np.mean(
            [
                value["equal_minus_k562_pds"]["estimate"]
                for value in surfaces.values()
            ]
        )
    )
    gates = {
        "equal_nonnegative_each_surface": all(
            value["equal_minus_k562_pds"]["estimate"] >= 0
            for value in surfaces.values()
        ),
        "equal_macro_at_least_0_003": macro_equal >= 0.003,
        "equal_bootstrap_positive_each_surface": all(
            value["equal_minus_k562_pds"]["lo"] > 0
            for value in surfaces.values()
        ),
        "equal_beats_null_each_surface": all(
            value["equal_minus_k562_pds"]["estimate"]
            > value["equal_null"]["p95_conservative"]
            for value in surfaces.values()
        ),
    }
    result = {
        "experiment": EXPERIMENT,
        "truth_read_after_freeze": True,
        "surfaces": surfaces,
        "macro_equal_minus_k562_pds": macro_equal,
        "gates": gates,
        "verdict": "PROMOTE" if all(gates.values()) else "NO-GO",
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
