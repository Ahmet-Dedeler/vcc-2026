#!/usr/bin/env python3
"""Build deterministic Orion control-only inputs for held-context STATE inference.

The selection procedure observes only rows labelled as controls. Perturbation
expression, perturbation depth, and perturbation identity never influence the
reservoir. This keeps the STATE query input independent of held perturbation
truth while matching the challenge's 18,400-control-cell budget.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DEFAULT_DATA = ROOT / "data/external/orion/data"
DEFAULT_META = ROOT / "data/external/orion/metadata/gene_metadata.parquet"
GENE_PATH = ROOT / "data/challenge/gene_names.csv"
CONTROL_LABEL = "non-targeting"
CONTROL_PATTERN = re.compile(
    r"^(?:non[-_\s]*targeting|ntc|control|safe[-_\s]*harbou?r|neg)" r"(?:[-_\s].*)?$",
    flags=re.IGNORECASE,
)
EXPECTED_LAYOUT = {
    "HCT116": {
        "labels": ("Non-Targeting",),
        "eligible_controls": 165_777,
        "row_groups": 109,
    },
    "HEK293T": {
        "labels": ("Non-Targeting",),
        "eligible_controls": 218_838,
        "row_groups": 223,
    },
}
DEFAULT_N_CONTROL = 18_400
DEFAULT_SEED = 106


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: object) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def is_control(name: object) -> bool:
    """Recognize only exact or separator-delimited control tokens."""
    return CONTROL_PATTERN.fullmatch(str(name).strip()) is not None


def batch_number(path: Path) -> int:
    match = re.search(r"Batch(\d+)", path.name)
    if match is None:
        raise ValueError(f"could not parse Orion batch number: {path}")
    return int(match.group(1))


def stable_seed(*parts: object) -> int:
    payload = "\0".join(map(str, parts)).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little")


def reservoir_slot(
    seen_controls: int,
    limit: int,
    rng: np.random.Generator,
) -> int | None:
    """Return the reservoir slot for the next control-only observation."""
    if seen_controls <= 0 or limit <= 0:
        raise ValueError("seen_controls and limit must be positive")
    if seen_controls <= limit:
        return seen_controls - 1
    candidate = int(rng.integers(0, seen_controls))
    return candidate if candidate < limit else None


def map_sparse_row(
    token_ids: Sequence[int],
    expression: Sequence[float],
    token_to_official: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    token_ids = np.asarray(token_ids, dtype=np.int64)
    expression = np.asarray(expression, dtype=np.float32)
    columns = token_to_official[token_ids]
    keep = columns >= 0
    return columns[keep].astype(np.int32), expression[keep]


@dataclass
class FileControlLayout:
    path: Path
    row_groups: int
    mixed_row_groups: int
    all_control_row_groups: int
    eligible_controls: int
    accepted_label_counts: dict[str, int]
    eligible_source_rows: np.ndarray


@dataclass
class ControlPreflight:
    line: str
    files: tuple[FileControlLayout, ...]
    row_groups: int
    mixed_row_groups: int
    all_control_row_groups: int
    eligible_controls: int
    accepted_label_counts: dict[str, int]


def inspect_control_layout(
    path: Path,
    *,
    allowed_labels: Sequence[str],
) -> FileControlLayout:
    """Inspect only label/filter columns and record the eligible source-row order."""
    parquet = pq.ParquetFile(path)
    allowed = set(map(str, allowed_labels))
    accepted_label_counts: dict[str, int] = {}
    eligible_source_rows: list[np.ndarray] = []
    mixed_row_groups = 0
    all_control_row_groups = 0
    row_offset = 0
    for row_group in range(parquet.num_row_groups):
        metadata = parquet.read_row_group(
            row_group,
            columns=["gene_target", "pass_guide_filter"],
        ).to_pandas()
        targets = metadata["gene_target"].astype(str).to_numpy()
        recognized = np.fromiter(
            (is_control(value) for value in targets),
            dtype=bool,
            count=len(targets),
        )
        unexpected = sorted(set(targets[recognized]) - allowed)
        if unexpected:
            raise ValueError(
                f"{path}: recognized unexpected control labels {unexpected}; "
                f"allowed labels are {sorted(allowed)}"
            )
        controls = np.fromiter(
            (value in allowed for value in targets),
            dtype=bool,
            count=len(targets),
        )
        eligible = controls & (metadata["pass_guide_filter"].to_numpy() == 1)
        for label, count in zip(*np.unique(targets[eligible], return_counts=True)):
            accepted_label_counts[str(label)] = accepted_label_counts.get(
                str(label), 0
            ) + int(count)
        if np.any(eligible):
            eligible_source_rows.append(row_offset + np.flatnonzero(eligible))
            if np.all(controls):
                all_control_row_groups += 1
            else:
                mixed_row_groups += 1
        row_offset += len(metadata)
    source_rows = (
        np.concatenate(eligible_source_rows).astype(np.int64, copy=False)
        if eligible_source_rows
        else np.empty(0, dtype=np.int64)
    )
    return FileControlLayout(
        path=path,
        row_groups=parquet.num_row_groups,
        mixed_row_groups=mixed_row_groups,
        all_control_row_groups=all_control_row_groups,
        eligible_controls=len(source_rows),
        accepted_label_counts=accepted_label_counts,
        eligible_source_rows=source_rows,
    )


def preflight_control_layout(
    line: str,
    files: Sequence[Path],
) -> ControlPreflight:
    """Validate the real mixed-row-group layout before any expression projection."""
    expected = EXPECTED_LAYOUT[line]
    allowed_labels = tuple(map(str, expected["labels"]))
    layouts = tuple(
        inspect_control_layout(path, allowed_labels=allowed_labels) for path in files
    )
    accepted_label_counts: dict[str, int] = {}
    for layout in layouts:
        for label, count in layout.accepted_label_counts.items():
            accepted_label_counts[label] = accepted_label_counts.get(label, 0) + count
    result = ControlPreflight(
        line=line,
        files=layouts,
        row_groups=sum(layout.row_groups for layout in layouts),
        mixed_row_groups=sum(layout.mixed_row_groups for layout in layouts),
        all_control_row_groups=sum(layout.all_control_row_groups for layout in layouts),
        eligible_controls=sum(layout.eligible_controls for layout in layouts),
        accepted_label_counts=accepted_label_counts,
    )
    if tuple(sorted(result.accepted_label_counts)) != tuple(sorted(allowed_labels)):
        raise ValueError(
            f"{line}: accepted control labels {sorted(result.accepted_label_counts)} "
            f"differ from expected {sorted(allowed_labels)}"
        )
    if result.eligible_controls != int(expected["eligible_controls"]):
        raise ValueError(
            f"{line}: found {result.eligible_controls} eligible controls; "
            f"expected {expected['eligible_controls']}"
        )
    if result.row_groups != int(expected["row_groups"]):
        raise ValueError(
            f"{line}: found {result.row_groups} row groups; "
            f"expected {expected['row_groups']}"
        )
    if (
        result.mixed_row_groups != result.row_groups
        or result.all_control_row_groups != 0
    ):
        raise ValueError(
            f"{line}: expected every row group to mix controls and perturbations; "
            f"found mixed={result.mixed_row_groups}, "
            f"all_control={result.all_control_row_groups}, "
            f"total={result.row_groups}"
        )
    return result


def read_predicate_filtered_controls(
    layout: FileControlLayout,
    *,
    accepted_labels: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Materialize only rows selected by exact control predicates in experiment code."""
    labels = tuple(map(str, accepted_labels))
    table = pq.read_table(
        layout.path,
        columns=[
            "gene_token_id",
            "gene_expression",
            "gene_target",
            "pass_guide_filter",
        ],
        filters=[
            ("pass_guide_filter", "=", 1),
            ("gene_target", "in", labels),
        ],
        use_threads=False,
    )
    returned_labels = table["gene_target"].to_pylist()
    returned_filter = table["pass_guide_filter"].to_numpy()
    if len(table) != layout.eligible_controls:
        raise ValueError(
            f"{layout.path}: predicate returned {len(table)} rows; "
            f"metadata preflight found {layout.eligible_controls}"
        )
    if any(label not in labels or not is_control(label) for label in returned_labels):
        raise ValueError(
            f"{layout.path}: filtered projection returned a non-control row"
        )
    if np.any(returned_filter != 1):
        raise ValueError(f"{layout.path}: filtered projection returned a failed guide")
    return (
        table["gene_token_id"].to_numpy(),
        table["gene_expression"].to_numpy(),
    )


def validate_existing(
    output_path: Path,
    *,
    line: str,
    n_control: int,
    seed: int,
) -> dict[str, object] | None:
    manifest_path = output_path.with_suffix(output_path.suffix + ".manifest.json")
    if not output_path.exists() and not manifest_path.exists():
        return None
    if not output_path.is_file() or not manifest_path.is_file():
        output_path.unlink(missing_ok=True)
        manifest_path.unlink(missing_ok=True)
        return None
    manifest = json.loads(manifest_path.read_text())
    if (
        manifest.get("experiment") != "106-state-baseline"
        or manifest.get("stage") != "held-control-only-input"
        or manifest.get("surface") != line.lower()
        or manifest.get("seed") != seed
        or manifest.get("control_rows") != n_control
        or manifest.get("contains_perturbation_rows") is not False
        or manifest.get("accepted_control_label_counts")
        != {
            label: count
            for label, count in zip(
                EXPECTED_LAYOUT[line]["labels"],
                (EXPECTED_LAYOUT[line]["eligible_controls"],),
            )
        }
        or manifest.get("expression_read_policy")
        != (
            "only predicate-filtered control rows are materialized by "
            "experiment code"
        )
        or manifest.get("output_sha256") != sha256_file(output_path)
    ):
        raise ValueError(f"existing control artifact failed validation: {output_path}")
    return manifest


def build_controls(
    *,
    line: str,
    data_dir: Path,
    metadata_path: Path,
    output_path: Path,
    n_control: int,
    seed: int,
) -> dict[str, object]:
    existing = validate_existing(
        output_path,
        line=line,
        n_control=n_control,
        seed=seed,
    )
    if existing is not None:
        print(output_path.with_suffix(output_path.suffix + ".manifest.json"))
        return existing
    if n_control <= 0:
        raise ValueError("n_control must be positive")

    official = tuple(pd.read_csv(GENE_PATH)["gene_name"].astype(str).tolist())
    official_index = {gene: index for index, gene in enumerate(official)}
    metadata = pd.read_parquet(metadata_path).sort_values("gene_token_id")
    n_tokens = int(metadata["gene_token_id"].max()) + 1
    token_names = np.full(n_tokens, "", dtype=object)
    token_names[metadata["gene_token_id"].to_numpy()] = metadata["gene_name"].astype(
        str
    )
    token_to_official = np.full(n_tokens, -1, dtype=np.int64)
    for token, gene in enumerate(token_names):
        column = official_index.get(str(gene))
        if column is not None:
            token_to_official[token] = column

    files = sorted(data_dir.glob(f"{line}_Batch*.parquet"), key=batch_number)
    if not files:
        raise FileNotFoundError(f"no Orion files for {line} under {data_dir}")
    preflight = preflight_control_layout(line, files)
    accepted_labels = tuple(EXPECTED_LAYOUT[line]["labels"])

    rows: list[tuple[np.ndarray, np.ndarray] | None] = [None] * n_control
    row_ids: list[str | None] = [None] * n_control
    source_batches: list[str | None] = [None] * n_control
    seen_controls = 0
    rng = np.random.default_rng(stable_seed("experiment-106-held-controls", line, seed))

    for file_index, layout in enumerate(preflight.files, 1):
        token_ids, expression = read_predicate_filtered_controls(
            layout,
            accepted_labels=accepted_labels,
        )
        for local_row, source_row in enumerate(layout.eligible_source_rows):
            seen_controls += 1
            slot = reservoir_slot(seen_controls, n_control, rng)
            if slot is None:
                continue
            rows[slot] = map_sparse_row(
                token_ids[local_row],
                expression[local_row],
                token_to_official,
            )
            row_ids[slot] = f"{layout.path.stem}:row{int(source_row)}"
            source_batches[slot] = layout.path.name
        if file_index % 20 == 0 or file_index == len(files):
            print(
                f"{line}: {file_index}/{len(files)} batches; "
                f"{seen_controls:,} eligible controls seen",
                flush=True,
            )

    if seen_controls != preflight.eligible_controls:
        raise AssertionError(
            f"{line}: reservoir saw {seen_controls} controls after metadata "
            f"preflight found {preflight.eligible_controls}"
        )
    if seen_controls < n_control or any(row is None for row in rows):
        raise ValueError(
            f"{line}: saw {seen_controls} eligible controls, need {n_control}"
        )
    selected = [row for row in rows if row is not None]
    indptr = np.zeros(n_control + 1, dtype=np.int64)
    np.cumsum([len(row[0]) for row in selected], out=indptr[1:])
    matrix = sp.csr_matrix(
        (
            np.concatenate([row[1] for row in selected]),
            np.concatenate([row[0] for row in selected]),
            indptr,
        ),
        shape=(n_control, len(official)),
        dtype=np.float32,
    )
    matrix.eliminate_zeros()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(
        tempfile.mkdtemp(prefix=f".{output_path.name}.", dir=output_path.parent)
    )
    staged_output = staging_dir / output_path.name
    try:
        ad.AnnData(
            X=matrix,
            obs=pd.DataFrame(
                {
                    "target_gene": [CONTROL_LABEL] * n_control,
                    "context": [line] * n_control,
                    "source_batch": source_batches,
                },
                index=row_ids,
            ),
            var=pd.DataFrame(index=official),
        ).write_h5ad(staged_output, compression="lzf")

        library_sizes = np.asarray(matrix.sum(axis=1)).ravel()
        manifest = {
            "experiment": "106-state-baseline",
            "stage": "held-control-only-input",
            "surface": line.lower(),
            "source": f"X-Atlas/Orion {line}",
            "selection": (
                "uniform reservoir over pass_guide_filter=1 and exact "
                "gene_target control-label predicates; perturbation expression "
                "and depth are not materialized by experiment code"
            ),
            "seed": seed,
            "control_rows": n_control,
            "eligible_controls_seen": seen_controls,
            "accepted_control_label_counts": preflight.accepted_label_counts,
            "genes": len(official),
            "gene_axis_sha256": hashlib.sha256(
                ("\n".join(official) + "\n").encode()
            ).hexdigest(),
            "source_files": [
                {
                    "name": path.name,
                    "bytes": path.stat().st_size,
                    "mtime_ns": path.stat().st_mtime_ns,
                }
                for path in files
            ],
            "metadata_sha256": sha256_file(metadata_path),
            "median_library_size": float(np.median(library_sizes)),
            "min_library_size": float(library_sizes.min()),
            "max_library_size": float(library_sizes.max()),
            "nnz": int(matrix.nnz),
            "contains_perturbation_rows": False,
            "expression_read_policy": (
                "only predicate-filtered control rows are materialized by "
                "experiment code"
            ),
            "metadata_preflight": {
                "row_groups": preflight.row_groups,
                "mixed_row_groups": preflight.mixed_row_groups,
                "all_control_row_groups": preflight.all_control_row_groups,
                "filtering_mandatory": preflight.mixed_row_groups > 0,
            },
        }
        manifest_path = output_path.with_suffix(output_path.suffix + ".manifest.json")
        manifest["output_sha256"] = sha256_file(staged_output)
        staged_manifest = staging_dir / manifest_path.name
        atomic_json(staged_manifest, manifest)
        os.replace(staged_output, output_path)
        os.replace(staged_manifest, manifest_path)
    finally:
        shutil.rmtree(staging_dir, ignore_errors=True)
    print(manifest_path)
    return manifest


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--line", choices=("HCT116", "HEK293T"), required=True)
    p.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    p.add_argument("--metadata", type=Path, default=DEFAULT_META)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--n-control", type=int, default=DEFAULT_N_CONTROL)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    return p


def main() -> None:
    args = parser().parse_args()
    build_controls(
        line=args.line,
        data_dir=args.data_dir.resolve(),
        metadata_path=args.metadata.resolve(),
        output_path=args.output.resolve(),
        n_control=args.n_control,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
