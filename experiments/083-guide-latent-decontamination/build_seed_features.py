#!/usr/bin/env python3
"""Build provenance-checked sequence/seed features for K562 guide pairs."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import re
import urllib.request
from collections import defaultdict
from collections.abc import Iterable
from pathlib import Path
from typing import BinaryIO

import numpy as np

SOURCE_COMMIT = "198fe3d84263002e9cfa17e1a2178d4e8dd0a669"
SOURCE_FILENAME = "replogle_k562_2kb_filt_no_lfc_cutoff_off_target_candidates.csv.gz"
SOURCE_URL = (
    "https://raw.githubusercontent.com/AustinHartman/"
    f"perturb_seed_reproduce/{SOURCE_COMMIT}/data/{SOURCE_FILENAME}"
)
EXPECTED_BYTES = 14_460_423
EXPECTED_SHA256 = "820f9e5ac81bbdbc1700e4520d7f95b0da6034f7640da1f9c0577e535454db8d"
GUIDE_PART = re.compile(
    r"^(?P<target>.+)_(?P<strand>[+-])_(?P<position>[0-9]+)"
    r"(?:\.[0-9]+)?-(?P<library>P1P2|P1|P2)$"
)
BASES = "ACGT"
KMER3 = tuple(a + b + c for a in BASES for b in BASES for c in BASES)
KMER_POS = {value: index for index, value in enumerate(KMER3)}


def sha256_file(path: Path, block_bytes: int = 4 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


class HashingReader:
    def __init__(self, raw: BinaryIO):
        self.raw = raw
        self.count = 0
        self.digest = hashlib.sha256()

    def read(self, size: int = -1) -> bytes:
        value = self.raw.read(size)
        self.count += len(value)
        if self.count > EXPECTED_BYTES:
            raise ValueError("compressed source exceeds pinned byte count")
        self.digest.update(value)
        return value

    def readable(self) -> bool:
        return True


def open_source(local_path: Path | None) -> tuple[BinaryIO, str]:
    if local_path is not None:
        return local_path.open("rb"), str(local_path.resolve())
    request = urllib.request.Request(
        SOURCE_URL, headers={"User-Agent": "vcc-exp083/1.0"}
    )
    return urllib.request.urlopen(request, timeout=120), SOURCE_URL


def normalize_sequence(value: str) -> str | None:
    sequence = str(value).strip().upper().replace("U", "T")
    if not sequence or any(base not in BASES for base in sequence):
        return None
    return sequence


def kmer3_frequencies(sequence: str) -> np.ndarray:
    result = np.zeros(len(KMER3), dtype=np.float64)
    if len(sequence) < 3:
        return result
    windows = 0
    for i in range(len(sequence) - 2):
        kmer = sequence[i : i + 3]
        index = KMER_POS.get(kmer)
        if index is not None:
            result[index] += 1
            windows += 1
    if windows:
        result /= windows
    return result


def positional_one_hot(sequence: str, length: int = 20) -> np.ndarray:
    result = np.zeros(length * len(BASES), dtype=np.float64)
    sequence = sequence[-length:].rjust(length, "N")
    for position, base in enumerate(sequence):
        if base in BASES:
            result[position * 4 + BASES.index(base)] = 1.0
    return result


def structural_features(guide_id: str) -> tuple[np.ndarray, tuple[str, ...]]:
    parts = guide_id.split("|")
    parsed = [GUIDE_PART.match(part) for part in parts]
    feature_names = (
        "pair_members",
        "same_strand",
        "opposite_strand",
        "log1p_pair_gap",
        "member_plus_fraction",
        "member_P1_fraction",
        "member_P2_fraction",
        "member_P1P2_fraction",
        "parse_complete",
    )
    if not parts or any(match is None for match in parsed):
        return np.zeros(len(feature_names), dtype=np.float64), feature_names
    rows = [match.groupdict() for match in parsed if match is not None]
    strands = [row["strand"] for row in rows]
    positions = [int(row["position"]) for row in rows]
    libraries = [row["library"] for row in rows]
    gap = max(positions) - min(positions) if len(positions) > 1 else 0
    values = np.asarray(
        [
            len(rows),
            float(len(set(strands)) == 1),
            float(len(set(strands)) > 1),
            np.log1p(gap),
            sum(value == "+" for value in strands) / len(rows),
            sum(value == "P1" for value in libraries) / len(rows),
            sum(value == "P2" for value in libraries) / len(rows),
            sum(value == "P1P2" for value in libraries) / len(rows),
            1.0,
        ],
        dtype=np.float64,
    )
    return values, feature_names


def build_features(
    guide_ids: Iterable[str],
    *,
    local_source: Path | None,
) -> tuple[np.ndarray, np.ndarray, tuple[str, ...], dict[str, object]]:
    guide_ids = tuple(map(str, guide_ids))
    guide_pos = {guide: index for index, guide in enumerate(guide_ids)}
    if len(guide_pos) != len(guide_ids):
        raise ValueError("duplicate guide IDs")
    sequences: dict[int, set[str]] = defaultdict(set)
    seeds: dict[int, set[str]] = defaultdict(set)
    seed_lengths: dict[int, list[float]] = defaultdict(list)
    hamming: dict[int, list[float]] = defaultdict(list)
    directions: dict[int, list[str]] = defaultdict(list)
    candidate_count = np.zeros(len(guide_ids), dtype=np.int64)

    raw, source = open_source(local_source)
    hashing = HashingReader(raw)
    rows = 0
    matched_rows = 0
    try:
        with gzip.GzipFile(fileobj=hashing) as compressed:
            reader = csv.DictReader(
                io.TextIOWrapper(compressed, encoding="utf-8", newline="")
            )
            required = {
                "guide_id",
                "protospacer_sequence",
                "seed",
                "seed_match_len",
                "guide_hamming_dist",
                "direction",
            }
            if reader.fieldnames is None or not required <= set(reader.fieldnames):
                raise ValueError("source table schema mismatch")
            for row in reader:
                rows += 1
                index = guide_pos.get(row["guide_id"])
                if index is None:
                    continue
                matched_rows += 1
                candidate_count[index] += 1
                sequence = normalize_sequence(row["protospacer_sequence"])
                seed = normalize_sequence(row["seed"])
                if sequence:
                    sequences[index].add(sequence)
                if seed:
                    seeds[index].add(seed)
                try:
                    seed_lengths[index].append(float(row["seed_match_len"]))
                except (TypeError, ValueError):
                    pass
                try:
                    hamming[index].append(float(row["guide_hamming_dist"]))
                except (TypeError, ValueError):
                    pass
                directions[index].append(str(row["direction"]).strip().lower())
    finally:
        raw.close()
    digest = hashing.digest.hexdigest()
    if hashing.count != EXPECTED_BYTES or digest != EXPECTED_SHA256:
        raise ValueError(
            f"source provenance mismatch: bytes={hashing.count}, sha256={digest}"
        )

    _, structure_names = structural_features(guide_ids[0] if guide_ids else "")
    feature_names = (
        *structure_names,
        *(
            f"protospacer_pos_{position}_{base}"
            for position in range(20)
            for base in BASES
        ),
        *(f"protospacer_3mer_{kmer}" for kmer in KMER3),
        *(f"seed_3mer_{kmer}" for kmer in KMER3),
        "log1p_candidate_count",
        "unique_protospacers",
        "unique_seeds",
        "protospacer_gc_mean",
        "seed_gc_mean",
        "seed_length_mean",
        "seed_length_max",
        "hamming_mean",
        "hamming_min",
        "direction_forward_fraction",
        "direction_reverse_fraction",
        "direction_both_fraction",
    )
    matrix = np.zeros((len(guide_ids), len(feature_names)), dtype=np.float32)
    identifiable = np.zeros(len(guide_ids), dtype=bool)
    for index, guide_id in enumerate(guide_ids):
        structure, _ = structural_features(guide_id)
        protospacers = sorted(sequences.get(index, set()))
        seed_values = sorted(seeds.get(index, set()))
        offset = 0
        matrix[index, offset : offset + len(structure)] = structure
        offset += len(structure)
        if protospacers:
            matrix[index, offset : offset + 80] = np.mean(
                [positional_one_hot(value) for value in protospacers], axis=0
            )
        offset += 80
        if protospacers:
            matrix[index, offset : offset + 64] = np.mean(
                [kmer3_frequencies(value) for value in protospacers], axis=0
            )
        offset += 64
        if seed_values:
            matrix[index, offset : offset + 64] = np.mean(
                [kmer3_frequencies(value) for value in seed_values], axis=0
            )
        offset += 64
        lengths = seed_lengths.get(index, [])
        distances = hamming.get(index, [])
        direction_values = directions.get(index, [])
        scalars = (
            np.log1p(candidate_count[index]),
            len(protospacers),
            len(seed_values),
            np.mean(
                [
                    sum(base in "GC" for base in value) / len(value)
                    for value in protospacers
                ]
            )
            if protospacers
            else 0.0,
            np.mean(
                [
                    sum(base in "GC" for base in value) / len(value)
                    for value in seed_values
                ]
            )
            if seed_values
            else 0.0,
            np.mean(lengths) if lengths else 0.0,
            np.max(lengths) if lengths else 0.0,
            np.mean(distances) if distances else 0.0,
            np.min(distances) if distances else 0.0,
            (
                sum(value == "forward" for value in direction_values)
                / len(direction_values)
                if direction_values
                else 0.0
            ),
            (
                sum(value == "reverse" for value in direction_values)
                / len(direction_values)
                if direction_values
                else 0.0
            ),
            (
                sum(value == "both" for value in direction_values)
                / len(direction_values)
                if direction_values
                else 0.0
            ),
        )
        matrix[index, offset:] = scalars
        identifiable[index] = bool(protospacers and seed_values)
    provenance = {
        "source": source,
        "repository": "AustinHartman/perturb_seed_reproduce",
        "repository_commit": SOURCE_COMMIT,
        "compressed_bytes": hashing.count,
        "compressed_sha256": digest,
        "data_rows": rows,
        "rows_matching_local_guides": matched_rows,
        "local_guides": len(guide_ids),
        "identifiable_guides": int(identifiable.sum()),
        "feature_count": len(feature_names),
        "forbidden_expression_features": [
            "target_log2fc",
            "offtarget_log2fc",
            "offtarget_gene",
            "offtarget_ensembl",
        ],
        "uses_expression_or_named_offtarget_columns": False,
    }
    return matrix, identifiable, tuple(feature_names), provenance


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--guide-ids", type=Path, required=True)
    parser.add_argument("--source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    guide_ids = tuple(
        line.strip() for line in args.guide_ids.read_text().splitlines() if line.strip()
    )
    matrix, identifiable, names, provenance = build_features(
        guide_ids, local_source=args.source
    )
    args.output.mkdir(parents=True, exist_ok=True)
    np.save(args.output / "guide_features.npy", matrix, allow_pickle=False)
    np.save(args.output / "guide_identifiable.npy", identifiable, allow_pickle=False)
    (args.output / "feature_names.txt").write_text("\n".join(names) + "\n")
    # Experiment 109 consumes this cache jointly with the compact guide-effect
    # cache.  Preserve the ordered guide axis as a first-class receipt so a
    # same-shaped but reordered feature matrix cannot be silently joined.
    (args.output / "guide_ids.txt").write_text("\n".join(guide_ids) + "\n")
    provenance["outputs"] = {
        "guide_features.npy": sha256_file(args.output / "guide_features.npy"),
        "guide_identifiable.npy": sha256_file(args.output / "guide_identifiable.npy"),
        "feature_names.txt": sha256_file(args.output / "feature_names.txt"),
        "guide_ids.txt": sha256_file(args.output / "guide_ids.txt"),
    }
    (args.output / "seed-feature-provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n"
    )
    print(json.dumps(provenance, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
