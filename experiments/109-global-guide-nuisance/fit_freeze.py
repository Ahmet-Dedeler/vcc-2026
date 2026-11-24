#!/usr/bin/env python3
"""Source-only fit/freeze scaffold for experiment 109.

This file intentionally does not extract expression or expose a real-data fit
command yet. It validates the preregistration and provides array-level helpers
for a later guarded runner built on experiment 083's compact caches.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import json
import os
import tempfile
import time
from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

os.environ["PYTHONDONTWRITEBYTECODE"] = "1"

import model
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
SOURCE_AXES_PATH = HERE / "source-axes.json"
EXPERIMENT = "109-global-guide-nuisance"
REQUIRED_SPLITS = ("fit", "tune", "test")
REQUIRED_FROZEN_ARTIFACTS = (
    "axes.json",
    "source-input-manifest.json",
    "source-cache-snapshot.json",
    "memory-guard.json",
    "source-test-model.npz",
    "source-test-holdouts.npz",
    "selected-model.npz",
    "tune-selection.json",
    "source-test.json",
    "source-bootstrap.npz",
    "source-nulls.npz",
    "source-null-selection.json",
    "singleton-corrections.npz",
    "held-label-nulls.npz",
)
REQUIRED_SOURCE_GATES = (
    "cosine_positive",
    "nmse_nonpositive",
    "bootstrap_lower_positive",
    "above_within_target_null_p95",
    "above_cross_target_null_p95",
    "direction_preserving_shrink_not_better",
)
SOURCE_FILES = {
    "model.py": HERE / "model.py",
    "fit_freeze.py": HERE / "fit_freeze.py",
    "source_only.py": HERE / "source_only.py",
    "evaluate.py": HERE / "evaluate.py",
    "inspect_source_axes.py": HERE / "inspect_source_axes.py",
    "experiment_083_model.py": model.EXP083_MODEL_PATH,
    "experiment_083_extract_guide_effects.py": (
        ROOT / "experiments/083-guide-latent-decontamination/extract_guide_effects.py"
    ),
    "experiment_083_build_seed_features.py": (
        ROOT / "experiments/083-guide-latent-decontamination/build_seed_features.py"
    ),
    "experiment_083_inspect_source.py": (
        ROOT / "experiments/083-guide-latent-decontamination/inspect_source.py"
    ),
}
SOURCE_INPUT_SCHEMA_VERSION = 4
SOURCE_TEST_HOLDOUT_SCHEMA_VERSION = 1
SINGLETON_CORRECTION_SCHEMA_VERSION = 2
HELD_LABEL_NULL_SCHEMA_VERSION = 1
SOURCE_LIFECYCLE_SCHEMA_VERSION = 1
SOURCE_CACHE_SNAPSHOT_SCHEMA_VERSION = 1
MEMORY_GUARD_RECEIPT_SCHEMA_VERSION = 1
FEATURE_ELIGIBILITY_PREDICATE = "normalized_valid_protospacer_and_seed"
HELD_LABEL_NULL_SEED = 1_092_027
GUIDE_CACHE_OUTPUTS = (
    "guide_delta.npy",
    "target_pooled_delta.npy",
    "guide_cells.npy",
    "target_cells.npy",
    "genes.txt",
    "guide_ids.txt",
    "guide_targets.txt",
    "target_names.txt",
)
FEATURE_CACHE_OUTPUTS = (
    "guide_features.npy",
    "guide_identifiable.npy",
    "feature_names.txt",
    "guide_ids.txt",
)
SOURCE_CACHE_FILES = {
    "guide": ("metadata.json", *GUIDE_CACHE_OUTPUTS),
    "feature": ("seed-feature-provenance.json", *FEATURE_CACHE_OUTPUTS),
}
LIFECYCLE_FILENAME = "source-lifecycle.json"
LOCK_FILENAME = ".source-fit.lock"
MEMORY_GUARD_PATH = ROOT / "scripts/run_with_memory_guard.py"
MEMORY_GUARD_SHA256 = "3be02ceaa44ea0c0ed4641827dd9305095255ec6cbacb82d667cefb4f7915198"


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_lines(values: Iterable[str]) -> str:
    return hashlib.sha256(("\n".join(map(str, values)) + "\n").encode()).hexdigest()


def sha256_array(values: np.ndarray) -> str:
    """Hash an array with its shape and dtype, independent of container format."""
    array = np.ascontiguousarray(np.asarray(values))
    digest = hashlib.sha256()
    digest.update(str(array.dtype).encode())
    digest.update(b"\0")
    digest.update(",".join(map(str, array.shape)).encode())
    digest.update(b"\0")
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


def assert_finite_allclose(
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    label: str,
    rtol: float,
    atol: float,
    block_values: int = 1 << 20,
) -> None:
    """Validate large arrays without NumPy testing's version-specific filtering."""
    observed = np.asarray(actual, dtype=np.float64)
    reference = np.asarray(expected, dtype=np.float64)
    if observed.shape != reference.shape:
        raise ValueError(
            f"{label} shape changed: {observed.shape} != {reference.shape}"
        )
    if block_values < 1:
        raise ValueError("allclose block size must be positive")
    observed_flat = observed.reshape(-1)
    reference_flat = reference.reshape(-1)
    for start in range(0, observed_flat.size, block_values):
        stop = min(start + block_values, observed_flat.size)
        observed_block = observed_flat[start:stop]
        reference_block = reference_flat[start:stop]
        if (
            not np.isfinite(observed_block).all()
            or not np.isfinite(reference_block).all()
        ):
            raise ValueError(f"{label} contains non-finite values")
        difference = np.abs(observed_block - reference_block)
        tolerance = atol + rtol * np.abs(reference_block)
        failed = difference > tolerance
        if np.any(failed):
            local = int(np.argmax(difference - tolerance))
            index = start + local
            raise ValueError(
                f"{label} differs at flat index {index}: "
                f"actual={observed_flat[index]!r} "
                f"expected={reference_flat[index]!r}"
            )


def sha256_json(payload: Mapping[str, object]) -> str:
    """Hash one JSON object in a canonical, path-independent representation."""
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def source_manifest_content(payload: Mapping[str, object]) -> dict[str, object]:
    """Return the cache-content contract without machine-local locations."""
    validated = validate_source_input_manifest(payload)
    return {
        key: value
        for key, value in validated.items()
        if key != "cache_locations"
    }


def source_cache_member_hashes(
    guide_cache: Path,
    feature_cache: Path,
    *,
    exact_members: bool = False,
) -> dict[str, dict[str, str]]:
    """Hash the exact cache members that may influence experiment 109."""
    result: dict[str, dict[str, str]] = {}
    for kind, directory in (("guide", guide_cache), ("feature", feature_cache)):
        if directory.is_symlink() or not directory.is_dir():
            raise ValueError(f"{kind} cache must be a real directory")
        expected = SOURCE_CACHE_FILES[kind]
        observed = tuple(sorted(path.name for path in directory.iterdir()))
        missing = sorted(set(expected) - set(observed))
        extras = sorted(set(observed) - set(expected))
        if missing or (exact_members and extras):
            raise ValueError(
                f"{kind} cache member set changed: "
                f"missing={missing} extras={extras}"
            )
        members: dict[str, str] = {}
        for filename in expected:
            path = directory / filename
            if path.is_symlink() or not path.is_file():
                raise ValueError(f"{kind} cache member is not a regular file: {filename}")
            members[filename] = sha256_file(path)
        result[kind] = members
    return result


def _atomic_write_json(path: Path, payload: Mapping[str, object]) -> None:
    """Durably replace one JSON marker without exposing a partial state."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def read_lifecycle(output: Path) -> dict[str, object] | None:
    path = output / LIFECYCLE_FILENAME
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    if (
        int(payload.get("schema_version", -1)) != SOURCE_LIFECYCLE_SCHEMA_VERSION
        or payload.get("experiment") != EXPERIMENT
    ):
        raise ValueError("source lifecycle marker is invalid")
    return payload


def write_lifecycle(
    output: Path,
    *,
    state: str,
    source_input_manifest_sha256: str | None = None,
) -> None:
    """Advance a one-way source/frozen/truth-open lifecycle marker."""
    allowed = {
        "source_open",
        "source_rejected",
        "source_frozen",
        "held_truth_opened",
    }
    if state not in allowed:
        raise ValueError("unknown source lifecycle state")
    previous = read_lifecycle(output)
    order = {
        "source_open": 0,
        "source_rejected": 1,
        "source_frozen": 1,
        "held_truth_opened": 2,
    }
    if previous is not None and order[state] < order[str(previous["state"])]:
        raise ValueError("source lifecycle cannot move backwards")
    if previous is not None and previous["state"] in {
        "source_rejected",
        "held_truth_opened",
    }:
        raise ValueError("terminal source lifecycle state forbids refit")
    payload: dict[str, object] = {
        "schema_version": SOURCE_LIFECYCLE_SCHEMA_VERSION,
        "experiment": EXPERIMENT,
        "state": state,
        "updated_unix_seconds": time.time(),
        "held_evaluation_enabled": False,
    }
    if source_input_manifest_sha256 is not None:
        if not is_sha256(source_input_manifest_sha256):
            raise ValueError("invalid source-input manifest hash for lifecycle")
        payload["source_input_manifest_sha256"] = source_input_manifest_sha256
    _atomic_write_json(output / LIFECYCLE_FILENAME, payload)


@contextlib.contextmanager
def exclusive_source_lock(output: Path):
    """Prevent simultaneous source fits or verification against moving files."""
    output.mkdir(parents=True, exist_ok=True)
    lock = output / LOCK_FILENAME
    try:
        descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError(f"source runner already active: {lock}") from exc
    try:
        os.write(descriptor, f"{os.getpid()}\n".encode())
        os.fsync(descriptor)
        yield
    finally:
        os.close(descriptor)
        lock.unlink(missing_ok=True)


def is_sha256(value: object) -> bool:
    text = str(value)
    return len(text) == 64 and all(char in "0123456789abcdef" for char in text)


def read_lines(path: Path, *, unique: bool = True) -> tuple[str, ...]:
    values = tuple(
        line.strip() for line in path.read_text().splitlines() if line.strip()
    )
    if unique and len(values) != len(set(values)):
        raise ValueError(f"{path}: duplicate labels")
    return values


def load_protocol(path: Path = PROTOCOL_PATH) -> dict[str, object]:
    protocol = json.loads(path.read_text())
    if protocol.get("experiment") != EXPERIMENT:
        raise ValueError("wrong experiment protocol")
    return protocol


def load_canonical_source_axes(
    protocol: Mapping[str, object] | None = None,
) -> dict[str, object]:
    frozen = load_protocol() if protocol is None else dict(protocol)
    contract = frozen["canonical_source_axes"]
    expected_path = ROOT / str(contract["path"])
    if expected_path.resolve() != SOURCE_AXES_PATH.resolve():
        raise ValueError("canonical source-axis path drift")
    if sha256_file(SOURCE_AXES_PATH) != str(contract["sha256"]):
        raise ValueError("canonical source-axis manifest changed")
    payload = json.loads(SOURCE_AXES_PATH.read_text())
    if payload.get("experiment") != EXPERIMENT:
        raise ValueError("canonical source axes name the wrong experiment")
    if payload.get("status") != contract["status"]:
        raise ValueError("canonical source-axis status drift")
    if contract.get("expression_values_opened") is not False:
        raise ValueError("canonical source axes are not metadata-only")
    return payload


def canonical_axis(
    name: str,
    source_axes: Mapping[str, object] | None = None,
) -> tuple[tuple[str, ...], model.AxisContract]:
    payload = load_canonical_source_axes() if source_axes is None else dict(source_axes)
    raw = payload["axes"][name]
    values = tuple(map(str, raw["values"]))
    contract = model.AxisContract(
        count=int(raw["count"]),
        sha256=str(raw["sha256"]),
    )
    model.validate_axis(
        values,
        contract,
        label=f"canonical {name}",
        require_sorted=name
        in {
            "panel_targets",
            "development_targets",
            "singleton_centering_targets",
            "perturbmatch_targets",
            "perturbmatch_genes",
        },
    )
    return values, contract


def validate_preregistration(
    protocol: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Fail if code constants drift from the frozen protocol."""
    frozen = dict(load_protocol() if protocol is None else protocol)
    if int(frozen.get("schema_version", -1)) != 5:
        raise ValueError("unsupported experiment-109 protocol schema")
    if (
        frozen.get("status")
        != "source_only_runner_ready_after_prefit_extraction_repair"
    ):
        raise ValueError("experiment-109 protocol status drift")
    split = frozen["fit_tune_test"]
    feature = frozen["feature_model"]
    grid = frozen["candidate_grid"]
    freeze = frozen["freeze_boundary"]
    source_axes = load_canonical_source_axes(frozen)
    if int(split["seed"]) != model.SPLIT_SEED:
        raise ValueError("split seed drift")
    if (
        float(feature["stable_singular_ratio_strictly_greater_than"])
        != model.STABLE_SINGULAR_RATIO
    ):
        raise ValueError("stable singular ratio drift")
    observed = model.candidate_grid()
    nonidentity = [candidate for candidate in observed if not candidate.is_identity]
    if sorted({candidate.response_rank for candidate in nonidentity}) != sorted(
        feature["response_ranks"]
    ):
        raise ValueError("response-rank grid drift")
    if sorted({candidate.ridge for candidate in nonidentity}) != sorted(
        feature["ridge_lambdas"]
    ):
        raise ValueError("ridge grid drift")
    if sorted({candidate.eta for candidate in observed}) != sorted(grid["eta"]):
        raise ValueError("eta grid drift")
    if tuple(grid["guide_cell_weight_power"]) != (0.0,):
        raise ValueError("guide-cell weight must be frozen to zero")
    if {candidate.weight_power for candidate in observed} != {0.0}:
        raise ValueError("guide-cell weight grid drift")
    if {candidate.huber_c for candidate in observed} != {float(grid["huber_c"])}:
        raise ValueError("Huber constant drift")
    if tuple(freeze["required_artifacts"]) != REQUIRED_FROZEN_ARTIFACTS:
        raise ValueError("required frozen artifact contract drift")
    if tuple(freeze["required_source_gates"]) != REQUIRED_SOURCE_GATES:
        raise ValueError("required source-gate contract drift")
    cross_null = frozen["nulls"]["cross_target_guide_row_permutation"]
    if (
        cross_null.get("split_local") is not True
        or cross_null.get("expression_or_effect_outcomes_used_for_strata") is not False
        or cross_null.get(
            "every_eligible_guide_row_and_target_must_move_in_every_replicate"
        )
        is not True
    ):
        raise ValueError("cross-target null isolation contract drift")
    guard = frozen["source_only_runner"]["memory_guard"]
    if (
        guard.get("script") != str(MEMORY_GUARD_PATH.relative_to(ROOT))
        or guard.get("sha256") != MEMORY_GUARD_SHA256
        or float(guard.get("max_rss_gib", -1)) <= 0
        or int(guard.get("min_free_percent", -1)) < 0
        or float(guard.get("poll_seconds", -1)) <= 0
    ):
        raise ValueError("memory-guard preregistration drift")
    reuse = frozen["reuse"]
    extractor_contract = reuse["experiment_083_extractor_pooling_contract"]
    if (
        extractor_contract.get("guide_effects")
        != "only guide categories with at least 20 cells"
        or extractor_contract.get("compact_target_reconstruction")
        != "all cells assigned to targets with at least 20 total cells"
        or float(extractor_contract.get("maximum_absolute_error", -1)) != 2e-5
    ):
        raise ValueError("experiment-083 extractor pooling contract drift")
    artifact_manifest = (
        ROOT / "experiments/083-guide-latent-decontamination/artifact_manifest.json"
    )
    if sha256_file(artifact_manifest) != str(
        reuse["experiment_083_artifact_manifest_sha256"]
    ):
        raise ValueError("experiment-083 artifact manifest changed")
    expected_exp083_sources = reuse["experiment_083_source_code_sha256"]
    observed_exp083_sources = {
        name: sha256_file(path)
        for name, path in SOURCE_FILES.items()
        if name.startswith("experiment_083_")
    }
    if observed_exp083_sources != expected_exp083_sources:
        raise ValueError("experiment-083 source dependency changed")
    if float(frozen["source_validation"]["sibling_center_weight_power"]) != 0.0:
        raise ValueError("source sibling-center weight must be frozen to zero")
    held_null = frozen["nulls"]["held_singleton_correction_label_permutation"]
    if int(held_null["seed"]) != HELD_LABEL_NULL_SEED:
        raise ValueError("held-label permutation seed drift")
    if "hct116_directional_credit_delta_minimum" in frozen["promotion_gates"]:
        raise ValueError("depth-confounded HCT116 directional gate is forbidden")
    if "feature_nullspace_prediction_change_abs_max" not in frozen["promotion_gates"]:
        raise ValueError("feature-nullspace invariance gate is missing")
    expectations = frozen["frozen_source_audit_expectations"]
    if int(source_axes["counts"]["development"]) != int(
        expectations["identifiable_nonpanel_multiguide_targets"]
    ):
        raise ValueError("canonical development target count drift")
    if int(source_axes["counts"]["singleton"]) != int(
        expectations["source_sequence_supported_singletons"]
    ):
        raise ValueError("canonical singleton target count drift")
    panel, _ = canonical_axis("panel_targets", source_axes)
    singletons, _ = canonical_axis("singleton_centering_targets", source_axes)
    if len(set(panel) & set(singletons)) != int(
        expectations["panel_sequence_supported_singletons"]
    ):
        raise ValueError("canonical panel singleton count drift")
    return frozen


def validate_source_input_manifest(
    payload_or_path: Mapping[str, object] | Path,
) -> dict[str, object]:
    """Validate the immutable receipt binding experiment-083 compact outputs."""
    payload = (
        json.loads(payload_or_path.read_text())
        if isinstance(payload_or_path, Path)
        else dict(payload_or_path)
    )
    required = {
        "schema_version",
        "experiment",
        "guide_axis",
        "gene_axis",
        "feature_axis",
        "cache_locations",
        "guide_cache",
        "feature_cache",
        "eligibility",
        "derived_arrays",
        "raw_sources",
        "source_code_sha256",
    }
    if set(payload) != required:
        raise ValueError("source-input manifest field set is invalid")
    if int(payload["schema_version"]) != SOURCE_INPUT_SCHEMA_VERSION:
        raise ValueError("unsupported source-input manifest schema")
    if payload["experiment"] != EXPERIMENT:
        raise ValueError("source-input manifest names the wrong experiment")
    locations = payload["cache_locations"]
    if set(locations) != {"guide_cache", "feature_cache"}:
        raise ValueError("source-input cache-location fields are invalid")
    for name, value in locations.items():
        location = Path(str(value))
        if not location.is_absolute():
            raise ValueError(f"{name} cache location must be absolute")

    source_axes = load_canonical_source_axes()
    raw_sources = payload["raw_sources"]
    expected_raw_sources = {
        "source_h5ad_sha256": source_axes["source"][
            "h5ad_sha256_from_pinned_source_audit"
        ],
        "source_h5ad_bytes": int(source_axes["source"]["h5ad_bytes"]),
        "seed_table_sha256": source_axes["source"]["seed_metadata"][
            "compressed_sha256"
        ],
        "panel_csv_sha256": source_axes["source"]["panel_csv_sha256"],
        "perturbmatch_perts_sha256": source_axes["source"]["perturbmatch_perts_sha256"],
        "perturbmatch_genes_sha256": source_axes["source"]["perturbmatch_genes_sha256"],
        "experiment_083_artifact_manifest_sha256": validate_preregistration()[
            "reuse"
        ]["experiment_083_artifact_manifest_sha256"],
    }
    if raw_sources != expected_raw_sources:
        raise ValueError("source-input raw provenance differs from pinned sources")
    genes, gene_contract = canonical_axis("gene_axis", source_axes)
    development, development_contract = canonical_axis(
        "development_targets",
        source_axes,
    )
    singletons, singleton_contract = canonical_axis(
        "singleton_centering_targets",
        source_axes,
    )
    panel, panel_contract = canonical_axis("panel_targets", source_axes)
    for name, expected in (
        ("gene_axis", gene_contract),
        ("development_target_axis", development_contract),
        ("singleton_target_axis", singleton_contract),
        ("panel_target_axis", panel_contract),
    ):
        observed = (
            payload["gene_axis"]
            if name == "gene_axis"
            else payload["eligibility"][name]
        )
        if (
            int(observed["count"]) != expected.count
            or observed["sha256"] != expected.sha256
        ):
            raise ValueError(
                f"source-input {name} contract differs from canonical axes"
            )

    guide_axis = payload["guide_axis"]
    feature_axis = payload["feature_axis"]
    if not is_sha256(guide_axis.get("sha256")) or int(guide_axis.get("count", 0)) < 1:
        raise ValueError("invalid guide-axis contract")
    if not is_sha256(feature_axis.get("sha256")):
        raise ValueError("invalid feature-axis hash")
    if int(feature_axis.get("count", 0)) != 229:
        raise ValueError("feature axis must contain exactly 229 features")

    guide_cache = payload["guide_cache"]
    feature_cache = payload["feature_cache"]
    if set(guide_cache.get("outputs", {})) != set(GUIDE_CACHE_OUTPUTS):
        raise ValueError("guide-cache output set is incomplete")
    if set(feature_cache.get("outputs", {})) != set(FEATURE_CACHE_OUTPUTS):
        raise ValueError("feature-cache output set is incomplete")
    hashes = (
        list(guide_cache["outputs"].values())
        + list(feature_cache["outputs"].values())
        + [
            guide_cache.get("metadata_sha256"),
            feature_cache.get("provenance_sha256"),
        ]
    )
    if not all(is_sha256(value) for value in hashes):
        raise ValueError("source-input cache hash is invalid")
    if guide_cache.get("guide_axis_sha256") != guide_axis["sha256"]:
        raise ValueError("guide cache is not bound to the shared guide axis")
    if feature_cache.get("guide_axis_sha256") != guide_axis["sha256"]:
        raise ValueError("feature cache is not bound to the shared guide axis")
    if guide_cache.get("gene_axis_sha256") != sha256_lines(genes):
        raise ValueError("guide cache is not bound to the canonical gene axis")
    if feature_cache.get("feature_axis_sha256") != feature_axis["sha256"]:
        raise ValueError("feature cache is not bound to the feature axis")

    eligibility = payload["eligibility"]
    if eligibility.get("predicate") != FEATURE_ELIGIBILITY_PREDICATE:
        raise ValueError(
            "feature eligibility predicate is not normalized-sequence exact"
        )
    expectations = validate_preregistration()["frozen_source_audit_expectations"]
    if int(eligibility.get("development_targets_count", -1)) != len(development):
        raise ValueError("feature-derived development target count changed")
    if int(eligibility.get("singleton_targets_count", -1)) != len(singletons):
        raise ValueError("feature-derived singleton target count changed")
    if int(eligibility.get("panel_singletons_count", -1)) != len(
        set(panel) & set(singletons)
    ):
        raise ValueError("feature-derived panel singleton count changed")
    if int(eligibility.get("identifiable_guides_count", 0)) < (
        2 * int(expectations["identifiable_nonpanel_multiguide_targets"])
        + int(expectations["source_sequence_supported_singletons"])
    ):
        raise ValueError("too few feature-identifiable guide rows")

    derived = payload["derived_arrays"]
    expected_derived = {
        "singleton_features",
        "panel_singleton_features",
        "development_features",
        "development_target_ids",
        "panel_singleton_source_index",
    }
    if set(derived) != expected_derived:
        raise ValueError("source-input derived-array set is incomplete")
    for name, contract in derived.items():
        if not is_sha256(contract.get("sha256")):
            raise ValueError(f"invalid derived-array hash: {name}")
        if not isinstance(contract.get("shape"), list) or not contract["shape"]:
            raise ValueError(f"invalid derived-array shape: {name}")
        if not isinstance(contract.get("dtype"), str) or not contract["dtype"]:
            raise ValueError(f"invalid derived-array dtype: {name}")
    if derived["singleton_features"]["shape"] != [len(singletons), 229]:
        raise ValueError("singleton feature matrix shape changed")
    if derived["panel_singleton_features"]["shape"] != [
        int(eligibility["panel_singletons_count"]),
        229,
    ]:
        raise ValueError("panel-singleton feature matrix shape changed")
    if derived["panel_singleton_source_index"]["shape"] != [len(panel)]:
        raise ValueError("panel singleton-index shape changed")
    if (
        derived["development_target_ids"]["shape"][0]
        != derived["development_features"]["shape"][0]
    ):
        raise ValueError("development feature/target rows differ")
    if derived["development_features"]["shape"][1] != 229:
        raise ValueError("development feature width changed")

    expected_sources = validate_preregistration()["reuse"][
        "experiment_083_source_code_sha256"
    ]
    if payload["source_code_sha256"] != expected_sources:
        raise ValueError("experiment-083 source dependency hashes changed")
    return payload


def build_source_input_manifest(
    guide_cache: Path,
    feature_cache: Path,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    """Build a hash-complete receipt and canonical derived feature arrays.

    The feature cache must contain a copied ``guide_ids.txt`` receipt. This is
    intentionally stricter than experiment 083 so same-shaped reordered
    feature matrices cannot be accepted.
    """
    guide_metadata_path = guide_cache / "metadata.json"
    feature_provenance_path = feature_cache / "seed-feature-provenance.json"
    guide_metadata = json.loads(guide_metadata_path.read_text())
    feature_provenance = json.loads(feature_provenance_path.read_text())
    guide_outputs = {
        name: sha256_file(guide_cache / name) for name in GUIDE_CACHE_OUTPUTS
    }
    feature_outputs = {
        name: sha256_file(feature_cache / name) for name in FEATURE_CACHE_OUTPUTS
    }
    if guide_metadata.get("outputs") != {
        name: guide_outputs[name] for name in GUIDE_CACHE_OUTPUTS
    }:
        raise ValueError("guide-cache metadata does not reproduce output hashes")
    expected_feature_outputs = {
        name: feature_outputs[name] for name in FEATURE_CACHE_OUTPUTS
    }
    if feature_provenance.get("outputs") != expected_feature_outputs:
        raise ValueError("feature provenance does not reproduce output hashes")
    if (
        feature_provenance.get("uses_expression_or_named_offtarget_columns")
        is not False
    ):
        raise ValueError(
            "feature cache used forbidden expression or named-gene columns"
        )
    if int(feature_provenance.get("feature_count", -1)) != 229:
        raise ValueError("feature cache width changed")

    guide_ids = read_lines(guide_cache / "guide_ids.txt")
    feature_guide_ids = read_lines(feature_cache / "guide_ids.txt")
    if feature_guide_ids != guide_ids:
        raise ValueError("feature and guide caches use different guide axes")
    guide_targets = read_lines(guide_cache / "guide_targets.txt", unique=False)
    genes = read_lines(guide_cache / "genes.txt")
    feature_names = read_lines(feature_cache / "feature_names.txt")
    guide_effects = np.load(
        guide_cache / "guide_delta.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    guide_cells = np.load(
        guide_cache / "guide_cells.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    features = np.load(
        feature_cache / "guide_features.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    identifiable = np.load(
        feature_cache / "guide_identifiable.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    if guide_effects.shape != (len(guide_ids), len(genes)):
        raise ValueError("guide-effect cache shape changed")
    if guide_cells.shape != (len(guide_ids),):
        raise ValueError("guide-cell cache shape changed")
    if features.shape != (len(guide_ids), len(feature_names)):
        raise ValueError("guide-feature cache shape changed")
    if identifiable.shape != (len(guide_ids),):
        raise ValueError("guide eligibility vector shape changed")
    if len(guide_targets) != len(guide_ids):
        raise ValueError("guide target axis length changed")

    source_axes = load_canonical_source_axes()
    canonical_genes, gene_contract = canonical_axis("gene_axis", source_axes)
    development, development_contract = canonical_axis(
        "development_targets",
        source_axes,
    )
    singletons, singleton_contract = canonical_axis(
        "singleton_centering_targets",
        source_axes,
    )
    panel, panel_contract = canonical_axis("panel_targets", source_axes)
    if genes != canonical_genes:
        raise ValueError("guide cache gene axis is not canonical")
    if (
        guide_metadata.get("source_h5ad_sha256")
        != source_axes["source"]["h5ad_sha256_from_pinned_source_audit"]
    ):
        raise ValueError("guide cache came from a different source H5AD")
    if int(guide_metadata.get("minimum_cells_per_guide", -1)) != int(
        source_axes["minimum_cells_per_guide"]
    ):
        raise ValueError("guide cache minimum-cell threshold changed")
    if guide_metadata.get("compact_reconstruction", {}).get("passes") is not True:
        raise ValueError("guide cache did not reproduce compact K562")
    if (
        feature_provenance.get("compressed_sha256")
        != source_axes["source"]["seed_metadata"]["compressed_sha256"]
    ):
        raise ValueError("feature cache came from different sequence metadata")
    grouped = model.EXP083.target_rows(guide_targets)
    identifiable = np.asarray(identifiable, dtype=bool)
    singleton_rows = []
    for target in singletons:
        rows = grouped.get(target)
        if rows is None or len(rows) != 1 or not identifiable[rows[0]]:
            raise ValueError(
                "canonical singleton axis disagrees with normalized feature eligibility"
            )
        singleton_rows.append(int(rows[0]))
    development_rows = []
    development_labels = []
    for target in development:
        rows = grouped.get(target)
        if rows is None:
            raise ValueError("canonical development target is absent from guide cache")
        usable = rows[identifiable[rows]]
        if len(usable) < 2:
            raise ValueError(
                "canonical development axis disagrees with feature eligibility"
            )
        development_rows.extend(map(int, usable))
        development_labels.extend([target] * len(usable))
    singleton_features = np.asarray(features[singleton_rows])
    development_features = np.asarray(features[development_rows])
    development_target_ids = np.asarray(development_labels)
    singleton_pos = {target: index for index, target in enumerate(singletons)}
    panel_index = np.asarray(
        [singleton_pos.get(target, -1) for target in panel],
        dtype=np.int64,
    )
    panel_singleton_features = np.asarray(
        singleton_features[panel_index[panel_index >= 0]]
    )
    derived = {
        "singleton_features": singleton_features,
        "panel_singleton_features": panel_singleton_features,
        "development_features": development_features,
        "development_target_ids": development_target_ids,
        "panel_singleton_source_index": panel_index,
    }
    payload = {
        "schema_version": SOURCE_INPUT_SCHEMA_VERSION,
        "experiment": EXPERIMENT,
        "guide_axis": {
            "count": len(guide_ids),
            "sha256": sha256_lines(guide_ids),
        },
        "gene_axis": {
            "count": gene_contract.count,
            "sha256": gene_contract.sha256,
        },
        "feature_axis": {
            "count": len(feature_names),
            "sha256": sha256_lines(feature_names),
        },
        "cache_locations": {
            "guide_cache": str(guide_cache.resolve()),
            "feature_cache": str(feature_cache.resolve()),
        },
        "guide_cache": {
            "metadata_sha256": sha256_file(guide_metadata_path),
            "guide_axis_sha256": sha256_lines(guide_ids),
            "gene_axis_sha256": sha256_lines(genes),
            "outputs": guide_outputs,
        },
        "feature_cache": {
            "provenance_sha256": sha256_file(feature_provenance_path),
            "guide_axis_sha256": sha256_lines(feature_guide_ids),
            "feature_axis_sha256": sha256_lines(feature_names),
            "outputs": feature_outputs,
        },
        "eligibility": {
            "predicate": FEATURE_ELIGIBILITY_PREDICATE,
            "identifiable_guides_count": int(np.sum(identifiable)),
            "development_targets_count": len(development),
            "singleton_targets_count": len(singletons),
            "panel_singletons_count": int(np.sum(panel_index >= 0)),
            "development_target_axis": {
                "count": development_contract.count,
                "sha256": development_contract.sha256,
            },
            "singleton_target_axis": {
                "count": singleton_contract.count,
                "sha256": singleton_contract.sha256,
            },
            "panel_target_axis": {
                "count": panel_contract.count,
                "sha256": panel_contract.sha256,
            },
        },
        "derived_arrays": {
            name: {
                "shape": list(values.shape),
                "dtype": str(values.dtype),
                "sha256": sha256_array(values),
            }
            for name, values in derived.items()
        },
        "raw_sources": {
            "source_h5ad_sha256": source_axes["source"][
                "h5ad_sha256_from_pinned_source_audit"
            ],
            "source_h5ad_bytes": int(source_axes["source"]["h5ad_bytes"]),
            "seed_table_sha256": source_axes["source"]["seed_metadata"][
                "compressed_sha256"
            ],
            "panel_csv_sha256": source_axes["source"]["panel_csv_sha256"],
            "perturbmatch_perts_sha256": source_axes["source"][
                "perturbmatch_perts_sha256"
            ],
            "perturbmatch_genes_sha256": source_axes["source"][
                "perturbmatch_genes_sha256"
            ],
            "experiment_083_artifact_manifest_sha256": validate_preregistration()[
                "reuse"
            ]["experiment_083_artifact_manifest_sha256"],
        },
        "source_code_sha256": dict(
            validate_preregistration()["reuse"][
                "experiment_083_source_code_sha256"
            ]
        ),
    }
    validate_source_input_manifest(payload)
    return payload, derived


def load_authenticated_source_cache(
    guide_cache: Path,
    feature_cache: Path,
    *,
    expected_manifest: Mapping[str, object] | Path | None = None,
) -> tuple[dict[str, object], dict[str, np.ndarray], dict[str, object]]:
    """Re-hash and load the only source inputs accepted by experiment 109.

    This is deliberately stricter than trusting a receipt written during an
    earlier run: every compact-cache output, metadata file, ordered axis and
    derived source array is rebuilt from the files that will control fitting.
    """
    guide_cache = guide_cache.resolve()
    feature_cache = feature_cache.resolve()
    payload, derived = build_source_input_manifest(guide_cache, feature_cache)
    if expected_manifest is not None:
        expected = validate_source_input_manifest(expected_manifest)
        if dict(expected) != payload:
            raise ValueError("authenticated compact caches differ from frozen receipt")

    guide_ids = read_lines(guide_cache / "guide_ids.txt")
    guide_targets = read_lines(guide_cache / "guide_targets.txt", unique=False)
    genes = read_lines(guide_cache / "genes.txt")
    guide_effects = np.load(
        guide_cache / "guide_delta.npy", mmap_mode="r", allow_pickle=False
    )
    guide_cells = np.load(
        guide_cache / "guide_cells.npy", mmap_mode="r", allow_pickle=False
    )
    features = np.load(
        feature_cache / "guide_features.npy", mmap_mode="r", allow_pickle=False
    )
    identifiable = np.asarray(
        np.load(feature_cache / "guide_identifiable.npy", mmap_mode="r", allow_pickle=False),
        dtype=bool,
    )
    if guide_effects.dtype.kind != "f" or features.dtype.kind != "f":
        raise ValueError("compact source matrices must use floating-point dtypes")
    if guide_cells.dtype.kind not in {"i", "u"}:
        raise ValueError("guide-cell counts must use an integer dtype")
    if (
        not np.isfinite(guide_effects).all()
        or not np.isfinite(features).all()
        or np.any(guide_cells < 1)
    ):
        raise ValueError("compact source cache contains invalid numerical values")
    cache: dict[str, object] = {
        "guide_ids": guide_ids,
        "guide_targets": guide_targets,
        "genes": genes,
        "guide_effects": guide_effects,
        "guide_cells": guide_cells,
        "features": features,
        "identifiable": identifiable,
    }
    return cache, derived, payload


def validate_source_cache_snapshot(
    frozen_dir: Path,
    *,
    expected_source_manifest: Mapping[str, object] | Path | None = None,
) -> tuple[dict[str, object], dict[str, np.ndarray], dict[str, object]]:
    """Authenticate and load the cache copy embedded in one frozen bundle."""
    receipt_path = frozen_dir / "source-cache-snapshot.json"
    receipt = json.loads(receipt_path.read_text())
    required = {
        "schema_version",
        "experiment",
        "guide_cache",
        "feature_cache",
        "members",
        "source_manifest_content_sha256",
        "original_member_hashes_before",
        "original_member_hashes_after",
    }
    if set(receipt) != required:
        raise ValueError("source-cache snapshot receipt field set is invalid")
    if int(receipt["schema_version"]) != SOURCE_CACHE_SNAPSHOT_SCHEMA_VERSION:
        raise ValueError("unsupported source-cache snapshot schema")
    if receipt["experiment"] != EXPERIMENT:
        raise ValueError("source-cache snapshot names the wrong experiment")
    expected_paths = {
        "guide_cache": "source-cache/guide",
        "feature_cache": "source-cache/feature",
    }
    for name, relative in expected_paths.items():
        if receipt[name] != relative:
            raise ValueError(f"source-cache snapshot path changed: {name}")
    guide_cache = (frozen_dir / expected_paths["guide_cache"]).resolve()
    feature_cache = (frozen_dir / expected_paths["feature_cache"]).resolve()
    root = frozen_dir.resolve()
    if root not in guide_cache.parents or root not in feature_cache.parents:
        raise ValueError("source-cache snapshot escapes the frozen directory")
    observed_members = source_cache_member_hashes(
        guide_cache,
        feature_cache,
        exact_members=True,
    )
    if receipt["members"] != observed_members:
        raise ValueError("source-cache snapshot member hashes changed")
    if receipt["original_member_hashes_before"] != observed_members:
        raise ValueError("source-cache snapshot differs from pre-compute source hashes")
    if receipt["original_member_hashes_after"] != observed_members:
        raise ValueError("source-cache snapshot differs from post-compute source hashes")
    expected_manifest = (
        validate_source_input_manifest(expected_source_manifest)
        if expected_source_manifest is not None
        else None
    )
    if expected_manifest is None:
        raise ValueError("source-cache snapshot requires its source-input manifest")
    content = source_manifest_content(expected_manifest)
    if sha256_json(content) != receipt["source_manifest_content_sha256"]:
        raise ValueError("source-cache snapshot manifest content hash changed")
    guide_ids = read_lines(guide_cache / "guide_ids.txt")
    feature_guide_ids = read_lines(feature_cache / "guide_ids.txt")
    guide_targets = read_lines(guide_cache / "guide_targets.txt", unique=False)
    feature_names = read_lines(feature_cache / "feature_names.txt")
    features = np.load(
        feature_cache / "guide_features.npy",
        mmap_mode="r",
        allow_pickle=False,
    )
    identifiable = np.asarray(
        np.load(
            feature_cache / "guide_identifiable.npy",
            mmap_mode="r",
            allow_pickle=False,
        ),
        dtype=bool,
    )
    if guide_ids != feature_guide_ids:
        raise ValueError("snapshot guide and feature axes differ")
    if (
        len(guide_targets) != len(guide_ids)
        or features.shape != (len(guide_ids), len(feature_names))
        or identifiable.shape != (len(guide_ids),)
    ):
        raise ValueError("snapshot null-design arrays do not align")
    if (
        sha256_lines(guide_ids) != expected_manifest["guide_axis"]["sha256"]
        or sha256_lines(feature_names) != expected_manifest["feature_axis"]["sha256"]
    ):
        raise ValueError("snapshot null-design axes differ from source receipt")
    if not np.isfinite(features).all():
        raise ValueError("snapshot sequence features contain non-finite values")
    cache: dict[str, object] = {
        "guide_ids": guide_ids,
        "guide_targets": guide_targets,
        "features": features,
        "identifiable": identifiable,
    }
    return cache, {}, expected_manifest


def validate_memory_guard_receipt(
    path: Path,
    *,
    expected_source_manifest_content_sha256: str,
    expected_source_cache_snapshot_sha256: str,
) -> dict[str, object]:
    """Verify that the complete heavy run was monitored by the pinned guard."""
    payload = json.loads(path.read_text())
    required = {
        "schema_version",
        "experiment",
        "status",
        "scope",
        "guard_script",
        "guard_script_sha256",
        "root_pid",
        "monitor_pid",
        "max_rss_gib",
        "min_free_percent",
        "poll_seconds",
        "peak_rss_gib",
        "samples",
        "last_process_count",
        "minimum_free_percent",
        "elapsed_seconds",
        "source_manifest_content_sha256",
        "source_cache_snapshot_sha256",
    }
    if set(payload) != required:
        raise ValueError("memory-guard receipt field set is invalid")
    protocol = validate_preregistration()
    expected = protocol["source_only_runner"]["memory_guard"]
    if int(payload["schema_version"]) != MEMORY_GUARD_RECEIPT_SCHEMA_VERSION:
        raise ValueError("unsupported memory-guard receipt schema")
    if payload["experiment"] != EXPERIMENT or payload["status"] != "complete":
        raise ValueError("memory guard did not complete successfully")
    if payload["scope"] != expected["scope"]:
        raise ValueError("memory-guard scope changed")
    if payload["guard_script"] != str(MEMORY_GUARD_PATH.relative_to(ROOT)):
        raise ValueError("memory-guard script path changed")
    if (
        payload["guard_script_sha256"] != MEMORY_GUARD_SHA256
        or sha256_file(MEMORY_GUARD_PATH) != MEMORY_GUARD_SHA256
    ):
        raise ValueError("memory-guard script hash changed")
    for name in ("max_rss_gib", "min_free_percent", "poll_seconds"):
        if float(payload[name]) != float(expected[name]):
            raise ValueError(f"memory-guard threshold changed: {name}")
    if int(payload["root_pid"]) <= 1 or int(payload["monitor_pid"]) <= 1:
        raise ValueError("memory-guard process identity is invalid")
    if int(payload["samples"]) < 1 or int(payload["last_process_count"]) < 1:
        raise ValueError("memory guard did not observe the source process")
    numeric = [
        float(payload["peak_rss_gib"]),
        float(payload["elapsed_seconds"]),
    ]
    if not np.isfinite(numeric).all() or min(numeric) < 0.0:
        raise ValueError("memory-guard measurements are invalid")
    if float(payload["peak_rss_gib"]) > float(payload["max_rss_gib"]):
        raise ValueError("memory-guard receipt records an RSS limit violation")
    minimum_free = payload["minimum_free_percent"]
    if minimum_free is not None and int(minimum_free) < int(payload["min_free_percent"]):
        raise ValueError("memory-guard receipt records a free-memory violation")
    if (
        payload["source_manifest_content_sha256"]
        != expected_source_manifest_content_sha256
    ):
        raise ValueError("memory guard is bound to a different source manifest")
    if payload["source_cache_snapshot_sha256"] != expected_source_cache_snapshot_sha256:
        raise ValueError("memory guard is bound to a different cache snapshot")
    return payload


def authenticated_guide_targets_from_manifest(
    payload_or_path: Mapping[str, object] | Path,
) -> tuple[str, ...]:
    """Read the guide-target axis only after authenticating its frozen receipt."""
    payload = validate_source_input_manifest(payload_or_path)
    guide_cache = Path(str(payload["cache_locations"]["guide_cache"]))
    target_path = guide_cache / "guide_targets.txt"
    if not target_path.is_file():
        raise FileNotFoundError(f"authenticated guide-target axis is missing: {target_path}")
    expected_hash = payload["guide_cache"]["outputs"]["guide_targets.txt"]
    if sha256_file(target_path) != expected_hash:
        raise ValueError("guide-target axis changed after source-input freeze")
    targets = read_lines(target_path, unique=False)
    if len(targets) != int(payload["guide_axis"]["count"]):
        raise ValueError("guide-target axis length differs from guide-axis contract")
    return targets


def development_split(
    identifiable_targets: Iterable[str],
    panel_targets: Iterable[str],
) -> dict[str, tuple[str, ...]]:
    """Create the immutable target split and assert complete isolation."""
    panel = set(map(str, panel_targets))
    development = tuple(sorted(set(map(str, identifiable_targets)) - panel))
    splits = model.deterministic_target_split(development, seed=model.SPLIT_SEED)
    model.assert_split_contract(splits, panel)
    combined = set().union(*(set(values) for values in splits.values()))
    if combined != set(development):
        raise AssertionError("split does not exactly cover development targets")
    return splits


def fit_candidate_from_arrays(
    guide_effects: np.ndarray,
    guide_features: np.ndarray,
    guide_targets: Sequence[str],
    cell_counts: np.ndarray,
    identifiable: np.ndarray,
    response_axis: model.AuthenticatedResponseAxis,
    fit_targets: Sequence[str],
    tune_targets: Sequence[str],
    singleton_target_features: np.ndarray,
    singleton_target_ids: Sequence[str],
    singleton_supported: np.ndarray,
    centering_axis_contract: model.AxisContract,
    candidate: model.Candidate,
    *,
    fixed_response_basis: np.ndarray | None = None,
) -> tuple[model.StableSubspaceModel, dict[str, float]]:
    """Fit on fit targets and score every guide holdout on tune targets."""
    fitted = model.fit_stable_subspace_model(
        guide_effects,
        guide_features,
        guide_targets,
        identifiable,
        fit_targets,
        response_axis,
        candidate,
        fixed_response_basis=fixed_response_basis,
    )
    nuisance_mean = model.singleton_nuisance_mean(
        fitted,
        singleton_target_features,
        singleton_target_ids,
        singleton_supported,
        centering_axis_contract=centering_axis_contract,
    )
    truth, baseline, predicted, labels, _ = model.held_guide_predictions(
        fitted,
        guide_effects,
        guide_features,
        guide_targets,
        cell_counts,
        identifiable,
        tune_targets,
        nuisance_mean,
    )
    return fitted, model.target_macro_metrics(truth, baseline, predicted, labels)


def recompute_tune_records(
    cache: Mapping[str, object],
    *,
    response_axis: model.AuthenticatedResponseAxis,
    splits: Mapping[str, Sequence[str]],
    singleton_features: np.ndarray,
    singleton_target_ids: Sequence[str],
    singleton_axis_contract: model.AxisContract,
) -> tuple[list[tuple[model.Candidate, dict[str, float]]], model.Candidate]:
    """Re-run the complete frozen tune grid from authenticated cache arrays."""
    required = {
        "guide_effects",
        "features",
        "guide_targets",
        "guide_cells",
        "identifiable",
    }
    if not required.issubset(cache):
        raise ValueError("source cache lacks arrays required for tune recomputation")
    support = np.ones(len(singleton_target_ids), dtype=bool)
    fitted_by_map: dict[tuple[int, float], model.StableSubspaceModel] = {}
    records: list[tuple[model.Candidate, dict[str, float]]] = []
    for candidate in model.candidate_grid():
        if candidate.is_identity:
            fitted, metrics = fit_candidate_from_arrays(
                cache["guide_effects"],
                cache["features"],
                cache["guide_targets"],
                cache["guide_cells"],
                cache["identifiable"],
                response_axis,
                splits["fit"],
                splits["tune"],
                singleton_features,
                singleton_target_ids,
                support,
                singleton_axis_contract,
                candidate,
            )
        else:
            key = (candidate.response_rank, candidate.ridge)
            if key not in fitted_by_map:
                template = model.Candidate(
                    candidate.response_rank,
                    candidate.ridge,
                    1.0,
                    candidate.weight_power,
                    candidate.huber_c,
                )
                fitted_by_map[key] = model.fit_stable_subspace_model(
                    cache["guide_effects"],
                    cache["features"],
                    cache["guide_targets"],
                    cache["identifiable"],
                    splits["fit"],
                    response_axis,
                    template,
                )
            fitted = model.clone_with_candidate(fitted_by_map[key], candidate)
            nuisance_mean = model.singleton_nuisance_mean(
                fitted,
                singleton_features,
                singleton_target_ids,
                support,
                centering_axis_contract=singleton_axis_contract,
            )
            truth, baseline, predicted, labels, _ = model.held_guide_predictions(
                fitted,
                cache["guide_effects"],
                cache["features"],
                cache["guide_targets"],
                cache["guide_cells"],
                cache["identifiable"],
                splits["tune"],
                nuisance_mean,
            )
            metrics = model.target_macro_metrics(truth, baseline, predicted, labels)
        records.append((candidate, metrics))
    selected, _ = select_candidate(records)
    return records, selected


def recompute_source_test_holdouts(
    cache: Mapping[str, object],
    *,
    response_axis: model.AuthenticatedResponseAxis,
    fit_tune_targets: Sequence[str],
    test_targets: Sequence[str],
    singleton_features: np.ndarray,
    singleton_target_ids: Sequence[str],
    singleton_axis_contract: model.AxisContract,
    candidate: model.Candidate,
) -> tuple[model.StableSubspaceModel, dict[str, object]]:
    """Refit once on fit+tune and derive every source-test value from cache."""
    support = np.ones(len(singleton_target_ids), dtype=bool)
    fitted = model.fit_stable_subspace_model(
        cache["guide_effects"],
        cache["features"],
        cache["guide_targets"],
        cache["identifiable"],
        fit_tune_targets,
        response_axis,
        candidate,
    )
    nuisance_mean = model.singleton_nuisance_mean(
        fitted,
        singleton_features,
        singleton_target_ids,
        support,
        centering_axis_contract=singleton_axis_contract,
    )
    truth, baseline, predicted, labels, held_rows = model.held_guide_predictions(
        fitted,
        cache["guide_effects"],
        cache["features"],
        cache["guide_targets"],
        cache["guide_cells"],
        cache["identifiable"],
        test_targets,
        nuisance_mean,
    )
    shrink_truth, _, shrink, shrink_labels, shrink_rows = (
        model.held_guide_direction_preserving_shrink(
            fitted,
            cache["guide_effects"],
            cache["features"],
            cache["guide_targets"],
            cache["guide_cells"],
            cache["identifiable"],
            test_targets,
            nuisance_mean,
        )
    )
    np.testing.assert_array_equal(shrink_truth, truth)
    if shrink_labels != labels:
        raise AssertionError("shrink holdouts changed source-test labels")
    np.testing.assert_array_equal(shrink_rows, held_rows)
    metrics = model.target_macro_metrics(truth, baseline, predicted, labels)
    shrink_metrics = model.target_macro_metrics(truth, baseline, shrink, labels)
    cosine_delta, nmse_delta = model.per_guide_metric_deltas(truth, predicted)
    return fitted, {
        "truth": truth,
        "predicted": predicted,
        "shrink_predicted": shrink,
        "labels": labels,
        "held_guide_rows": held_rows,
        "metrics": metrics,
        "shrink_metrics": shrink_metrics,
        "per_guide_cosine_delta": cosine_delta,
        "per_guide_nmse_delta": nmse_delta,
    }


def validate_preregistered_final_invariants(
    fitted: model.StableSubspaceModel,
    *,
    source_derived: Mapping[str, np.ndarray],
    source_input_manifest: Mapping[str, object] | Path,
) -> dict[str, float | int]:
    """Bind the final model to the exact frozen panel-singleton feature subset."""
    manifest = validate_source_input_manifest(source_input_manifest)
    required = {
        "singleton_features",
        "panel_singleton_features",
        "panel_singleton_source_index",
    }
    if not required.issubset(source_derived):
        raise ValueError("source-derived panel-singleton evidence is incomplete")
    singleton_features = np.asarray(source_derived["singleton_features"])
    panel_features = np.asarray(source_derived["panel_singleton_features"])
    panel_index = np.asarray(
        source_derived["panel_singleton_source_index"],
        dtype=np.int64,
    )
    for name, values in (
        ("singleton_features", singleton_features),
        ("panel_singleton_features", panel_features),
        ("panel_singleton_source_index", panel_index),
    ):
        contract = manifest["derived_arrays"][name]
        if list(values.shape) != contract["shape"]:
            raise ValueError(f"final invariant array shape differs: {name}")
        if str(values.dtype) != contract["dtype"]:
            raise ValueError(f"final invariant array dtype differs: {name}")
        if sha256_array(values) != contract["sha256"]:
            raise ValueError(f"final invariant array hash differs: {name}")
    valid = panel_index >= 0
    if np.any(panel_index[valid] >= len(singleton_features)):
        raise ValueError("panel-singleton source index is out of bounds")
    expected_panel_features = singleton_features[panel_index[valid]]
    np.testing.assert_array_equal(panel_features, expected_panel_features)

    expected = validate_preregistration()["frozen_source_audit_expectations"]
    varying = int(np.sum(fitted.feature_keep))
    stable_rank = fitted.stable_rank
    if varying != int(expected["varying_fit_features"]):
        raise ValueError(f"varying feature invariant failed: {varying}")
    if stable_rank != int(expected["stable_contrast_dimensions"]):
        raise ValueError(f"stable rank invariant failed: {stable_rank}")
    standardized = fitted.standardized_features(panel_features)
    projected = fitted.projected_standardized_features(panel_features)
    denominator = np.maximum(np.linalg.norm(standardized, axis=1), model.EPS)
    fraction = np.linalg.norm(projected, axis=1) / denominator
    median_fraction = float(np.median(fraction))
    registered = float(expected["median_panel_singleton_stable_subspace_norm_fraction"])
    if round(median_fraction, 3) != round(registered, 3):
        raise ValueError(
            "median singleton stable-subspace projection invariant failed: "
            f"{median_fraction:.12f}"
        )
    return {
        "varying_fit_features": varying,
        "stable_contrast_dimensions": stable_rank,
        "median_panel_singleton_stable_subspace_norm_fraction": median_fraction,
    }


def candidate_passes_tune(metrics: Mapping[str, float]) -> bool:
    return (
        int(metrics["targets"]) > 0
        and float(metrics["cosine_delta"]) > 0.0
        and float(metrics["nmse_delta"]) <= 0.0
    )


def selection_key(
    candidate: model.Candidate,
    metrics: Mapping[str, float],
) -> tuple[float, float, int, float, float]:
    """Frozen deterministic tie-break order from the protocol."""
    return (
        float(metrics["cosine_delta"]),
        -float(metrics["nmse_delta"]),
        -candidate.response_rank,
        -candidate.eta,
        -candidate.ridge,
    )


def select_candidate(
    records: Sequence[tuple[model.Candidate, Mapping[str, float]]],
) -> tuple[model.Candidate, dict[str, float]]:
    """Select once from one complete, duplicate-free frozen candidate grid."""
    expected = set(model.candidate_grid())
    observed = [candidate for candidate, _ in records]
    if len(observed) != len(set(observed)):
        raise ValueError("candidate selection records contain duplicates")
    if set(observed) != expected:
        missing = sorted(expected - set(observed))
        extra = sorted(set(observed) - expected)
        raise ValueError(
            f"candidate selection grid mismatch: missing={missing[:2]} "
            f"extra={extra[:2]}"
        )
    passing = [
        (candidate, dict(metrics))
        for candidate, metrics in records
        if candidate_passes_tune(metrics)
    ]
    if not passing:
        identity = next(candidate for candidate in expected if candidate.is_identity)
        metrics = next(
            dict(values) for candidate, values in records if candidate == identity
        )
        return identity, metrics
    return max(passing, key=lambda item: selection_key(item[0], item[1]))


def target_cluster_bootstrap(
    per_guide_delta: np.ndarray,
    labels: Sequence[str],
    *,
    repetitions: int,
    seed: int,
) -> np.ndarray:
    """Bootstrap target-macro deltas without overweighting guide-rich targets."""
    values = np.asarray(per_guide_delta, dtype=np.float64)
    target_labels = np.asarray(tuple(map(str, labels)), dtype=object)
    if values.ndim != 1 or len(values) != len(target_labels):
        raise ValueError("bootstrap value/label mismatch")
    if repetitions < 1:
        raise ValueError("bootstrap repetitions must be positive")
    targets = tuple(sorted(set(target_labels.tolist())))
    if len(targets) < 2:
        raise ValueError("target-cluster bootstrap needs at least two targets")
    target_means = np.asarray(
        [np.mean(values[target_labels == target]) for target in targets],
        dtype=np.float64,
    )
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(targets), size=(repetitions, len(targets)))
    return np.mean(target_means[draws], axis=1)


def within_target_permutation_indices(
    guide_targets: Sequence[str],
    *,
    seed: int,
    row_mask: np.ndarray | None = None,
    targets: Iterable[str] | None = None,
) -> np.ndarray:
    """Return a deterministic within-target guide-row permutation."""
    labels = np.asarray(tuple(map(str, guide_targets)), dtype=object)
    eligible = (
        np.ones(len(labels), dtype=bool)
        if row_mask is None
        else np.asarray(row_mask, dtype=bool)
    )
    if eligible.shape != (len(labels),):
        raise ValueError("within-target permutation row mask shape mismatch")
    selected_targets = (
        set(labels.tolist()) if targets is None else set(map(str, targets))
    )
    if not selected_targets.issubset(set(labels.tolist())):
        raise ValueError("within-target permutation names unknown targets")
    rng = np.random.default_rng(seed)
    permutation = np.arange(len(labels), dtype=np.int64)
    for target in sorted(selected_targets):
        rows = np.flatnonzero((labels == target) & eligible)
        if len(rows) < 2:
            raise ValueError(
                f"within-target permutation needs two eligible rows: {target}"
            )
        shift = int(rng.integers(1, len(rows)))
        permutation[rows] = np.roll(rows, shift)
        if np.any(permutation[rows] == rows):
            raise AssertionError("within-target permutation left a guide row fixed")
    return permutation


def permute_features_within_target(
    features: np.ndarray,
    guide_targets: Sequence[str],
    *,
    seed: int,
    row_mask: np.ndarray | None = None,
    targets: Iterable[str] | None = None,
) -> np.ndarray:
    """Permute guide-feature assignments independently inside every target."""
    values = np.asarray(features)
    if values.ndim != 2 or len(values) != len(guide_targets):
        raise ValueError("within-target permutation shape mismatch")
    permutation = within_target_permutation_indices(
        guide_targets,
        seed=seed,
        row_mask=row_mask,
        targets=targets,
    )
    return values[permutation]


def _cross_target_covariates(
    features: np.ndarray,
    guide_targets: Sequence[str],
    identifiable: np.ndarray,
    targets: Sequence[str],
) -> dict[str, tuple[int, float]]:
    """Build target covariates without reading expression/effect outcomes."""
    values = np.asarray(features, dtype=np.float64)
    labels = tuple(map(str, guide_targets))
    usable = np.asarray(identifiable, dtype=bool)
    if (
        values.ndim != 2
        or len(values) != len(labels)
        or usable.shape != (len(labels),)
    ):
        raise ValueError("cross-target covariate arrays do not align")
    if not np.isfinite(values).all():
        raise ValueError("cross-target covariates contain invalid values")
    grouped = model.EXP083.target_rows(labels)
    result: dict[str, tuple[int, float]] = {}
    for target in targets:
        rows = grouped.get(str(target))
        if rows is None:
            raise ValueError(f"cross-target null target is absent: {target}")
        rows = rows[usable[rows]]
        if len(rows) < 2:
            raise ValueError(f"cross-target null target lacks two guides: {target}")
        feature_norm = float(np.linalg.norm(values[rows], ord="fro"))
        result[str(target)] = (
            len(rows),
            feature_norm,
        )
    return result


def _bucket_has_row_derangement(
    targets: Sequence[str],
    grouped: Mapping[str, np.ndarray],
    identifiable: np.ndarray,
) -> bool:
    rows = np.concatenate(
        [grouped[target][identifiable[grouped[target]]] for target in sorted(targets)]
    )
    labels = np.asarray(
        [
            target
            for target in sorted(targets)
            for _ in grouped[target][identifiable[grouped[target]]]
        ],
        dtype=object,
    )
    for shift in range(1, len(rows)):
        if np.all(np.roll(labels, shift) != labels):
            return True
    return False


def build_split_local_cross_target_strata(
    features: np.ndarray,
    guide_targets: Sequence[str],
    identifiable: np.ndarray,
    splits: Mapping[str, Sequence[str]],
) -> tuple[dict[str, str], dict[str, object]]:
    """Create outcome-free strata, coarsening until every target can move."""
    labels = tuple(map(str, guide_targets))
    usable = np.asarray(identifiable, dtype=bool)
    grouped = model.EXP083.target_rows(labels)
    target_strata: dict[str, str] = {}
    split_audit: dict[str, object] = {}
    levels = (
        "guide_count_feature_norm_decile",
        "guide_count",
        "all_targets",
    )
    if set(splits) != set(REQUIRED_SPLITS):
        raise ValueError("cross-target null split set is incomplete")
    split_sets = [set(map(str, splits[name])) for name in REQUIRED_SPLITS]
    if any(left & right for i, left in enumerate(split_sets) for right in split_sets[i + 1 :]):
        raise ValueError("cross-target null splits overlap")
    for split_name in REQUIRED_SPLITS:
        targets = tuple(sorted(map(str, splits[split_name])))
        if not targets or len(targets) != len(set(targets)):
            raise ValueError(f"{split_name} cross-target null targets are invalid")
        covariates = _cross_target_covariates(
            features,
            labels,
            usable,
            targets,
        )
        norms = np.asarray([covariates[target][1] for target in targets])
        quantiles = np.quantile(norms, np.linspace(0.0, 1.0, 11))
        norm_bins = np.searchsorted(quantiles[1:-1], norms, side="right")
        chosen: dict[str, str] | None = None
        chosen_level = ""
        for level in levels:
            candidate: dict[str, str] = {}
            for target, norm_bin in zip(targets, norm_bins, strict=True):
                guide_count, _ = covariates[target]
                if level == "guide_count_feature_norm_decile":
                    key = (
                        f"guides={guide_count}|"
                        f"feature_norm_decile={int(norm_bin)}"
                    )
                elif level == "guide_count":
                    key = f"guides={guide_count}"
                else:
                    key = "all"
                candidate[target] = f"split={split_name}|{key}"
            buckets: dict[str, list[str]] = {}
            for target, key in candidate.items():
                buckets.setdefault(key, []).append(target)
            if all(
                len(bucket_targets) >= 2
                and _bucket_has_row_derangement(bucket_targets, grouped, usable)
                for bucket_targets in buckets.values()
            ):
                chosen = candidate
                chosen_level = level
                break
        if chosen is None:
            raise ValueError(
                f"{split_name} cross-target null cannot move every eligible target"
            )
        target_strata.update(chosen)
        split_audit[split_name] = {
            "targets": len(targets),
            "coarsening_level": chosen_level,
            "strata": len(set(chosen.values())),
            "all_targets_movable": True,
        }
    return target_strata, {
        "covariates": [
            "identifiable_guide_count",
            "sequence_feature_bundle_frobenius_norm_decile",
        ],
        "outcome_values_used": False,
        "expression_values_used": False,
        "split_local": True,
        "splits": split_audit,
    }


def cross_target_permutation_indices(
    guide_targets: Sequence[str],
    target_strata: Mapping[str, str],
    *,
    seed: int,
    row_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Permute guide rows across targets, split-locally and without fixed targets."""
    labels = np.asarray(tuple(map(str, guide_targets)), dtype=object)
    grouped = model.EXP083.target_rows(labels.tolist())
    requested = set(map(str, target_strata))
    if not requested or not requested.issubset(grouped):
        raise ValueError("target strata do not identify source feature rows")
    eligible = (
        np.ones(len(labels), dtype=bool)
        if row_mask is None
        else np.asarray(row_mask, dtype=bool)
    )
    if eligible.shape != (len(labels),):
        raise ValueError("cross-target permutation row mask shape mismatch")
    buckets: dict[str, list[str]] = {}
    for target, stratum in target_strata.items():
        rows = grouped[target][eligible[grouped[target]]]
        if len(rows) < 2:
            raise ValueError("cross-target permutation target has too few guides")
        buckets.setdefault(str(stratum), []).append(target)
    rng = np.random.default_rng(seed)
    permutation = np.arange(len(labels), dtype=np.int64)
    for stratum in sorted(buckets):
        targets = sorted(buckets[stratum])
        rows = np.concatenate([grouped[target][eligible[grouped[target]]] for target in targets])
        row_labels = labels[rows]
        valid_shifts = [
            shift
            for shift in range(1, len(rows))
            if np.all(np.roll(row_labels, shift) != row_labels)
        ]
        if not valid_shifts:
            raise ValueError(
                f"cross-target stratum cannot derange every guide row: {stratum}"
            )
        shift = valid_shifts[int(rng.integers(0, len(valid_shifts)))]
        donors = np.roll(rows, shift)
        permutation[rows] = donors
        if np.any(labels[permutation[rows]] == labels[rows]):
            raise AssertionError("cross-target permutation left a target fixed")
    return permutation


def permute_target_feature_rows(
    features: np.ndarray,
    guide_targets: Sequence[str],
    target_strata: Mapping[str, str],
    *,
    seed: int,
    row_mask: np.ndarray | None = None,
) -> np.ndarray:
    """Apply the split-local cross-target guide-row null."""
    values = np.asarray(features)
    if values.ndim != 2 or len(values) != len(guide_targets):
        raise ValueError("cross-target row permutation shape mismatch")
    permutation = cross_target_permutation_indices(
        guide_targets,
        target_strata,
        seed=seed,
        row_mask=row_mask,
    )
    return values[permutation]


def source_test_decision(
    metrics: Mapping[str, float],
    *,
    bootstrap_samples: np.ndarray,
    within_target_null: np.ndarray,
    cross_target_null: np.ndarray,
    direction_preserving_shrink_delta: float,
) -> dict[str, object]:
    """Apply every source-side promotion gate before a freeze can be built."""
    bootstrap = np.asarray(bootstrap_samples, dtype=np.float64)
    within = np.asarray(within_target_null, dtype=np.float64)
    cross = np.asarray(cross_target_null, dtype=np.float64)
    if min(len(bootstrap), len(within), len(cross)) < 1:
        raise ValueError("source bootstrap/null arrays must be nonempty")
    if not all(np.isfinite(values).all() for values in (bootstrap, within, cross)):
        raise ValueError("source bootstrap/null arrays contain non-finite values")
    observed = float(metrics["cosine_delta"])
    nmse = float(metrics["nmse_delta"])
    shrink = float(direction_preserving_shrink_delta)
    if not np.isfinite([observed, nmse, shrink]).all():
        raise ValueError("source metrics contain non-finite values")
    bootstrap_lower = float(np.quantile(bootstrap, 0.025))
    within_p95 = float(np.quantile(within, 0.95))
    cross_p95 = float(np.quantile(cross, 0.95))
    gates = {
        "cosine_positive": observed > 0.0,
        "nmse_nonpositive": nmse <= 0.0,
        "bootstrap_lower_positive": bootstrap_lower > 0.0,
        "above_within_target_null_p95": observed > within_p95,
        "above_cross_target_null_p95": observed > cross_p95,
        "direction_preserving_shrink_not_better": observed > shrink,
    }
    return {
        "passes": all(gates.values()),
        "gates": gates,
        "metrics": {
            "targets": int(metrics["targets"]),
            "guide_holdouts": int(metrics["guide_holdouts"]),
            "cosine_delta": observed,
            "nmse_delta": nmse,
        },
        "bootstrap_lower_95": bootstrap_lower,
        "within_target_null_p95": within_p95,
        "cross_target_null_p95": cross_p95,
        "direction_preserving_shrink_delta": shrink,
    }


def validate_source_test(source_test: Mapping[str, object]) -> None:
    """Reject a failed, partial, or self-contradictory source authorization."""
    if source_test.get("passes") is not True:
        raise ValueError("source test did not pass")
    gates = source_test.get("gates")
    if not isinstance(gates, Mapping):
        raise TypeError("source test has no gate map")
    if set(gates) != set(REQUIRED_SOURCE_GATES):
        raise ValueError("source test gate set is incomplete")
    if not all(gates[name] is True for name in REQUIRED_SOURCE_GATES):
        raise ValueError("one or more mandatory source gates failed")
    metrics = source_test.get("metrics")
    if not isinstance(metrics, Mapping):
        raise TypeError("source test metrics are missing")
    required_metrics = {"targets", "guide_holdouts", "cosine_delta", "nmse_delta"}
    if not required_metrics.issubset(metrics):
        raise ValueError("source test metrics are incomplete")
    expected_test_targets = int(
        validate_preregistration()["frozen_source_audit_expectations"]["split_counts"][
            "test"
        ]
    )
    if int(metrics["targets"]) != expected_test_targets:
        raise ValueError("source test does not cover the canonical test targets")
    if int(metrics["guide_holdouts"]) < expected_test_targets * 2:
        raise ValueError("source test has too few guide holdouts")
    numeric = [
        float(metrics["cosine_delta"]),
        float(metrics["nmse_delta"]),
        float(source_test["bootstrap_lower_95"]),
        float(source_test["within_target_null_p95"]),
        float(source_test["cross_target_null_p95"]),
        float(source_test["direction_preserving_shrink_delta"]),
    ]
    if not np.isfinite(numeric).all():
        raise ValueError("source test contains non-finite evidence")
    if float(metrics["cosine_delta"]) <= 0.0:
        raise ValueError("source cosine gate is inconsistent")
    if float(metrics["nmse_delta"]) > 0.0:
        raise ValueError("source NMSE gate is inconsistent")
    observed = float(metrics["cosine_delta"])
    if float(source_test["bootstrap_lower_95"]) <= 0.0:
        raise ValueError("source bootstrap gate is inconsistent")
    if observed <= float(source_test["within_target_null_p95"]):
        raise ValueError("within-target null gate is inconsistent")
    if observed <= float(source_test["cross_target_null_p95"]):
        raise ValueError("cross-target null gate is inconsistent")
    if observed <= float(source_test["direction_preserving_shrink_delta"]):
        raise ValueError("direction-preserving shrink gate is inconsistent")


def validate_source_null_selection(
    path: Path,
    *,
    cache: Mapping[str, object],
    splits: Mapping[str, Sequence[str]],
    within_target_values: np.ndarray,
    cross_target_values: np.ndarray,
    within_target_permutations: np.ndarray,
    cross_target_permutations: np.ndarray,
) -> dict[str, object]:
    """Recompute null assignments and independently reselect every winner."""
    payload = json.loads(path.read_text())
    if set(payload) != {
        "schema_version",
        "within_target",
        "cross_target",
        "cross_target_strata",
    }:
        raise ValueError("source-null selection field set is invalid")
    if int(payload["schema_version"]) != 2:
        raise ValueError("unsupported source-null selection schema")
    protocol = validate_preregistration()
    expected_repetitions = {
        "within_target": int(
            protocol["nulls"]["within_target_guide_feature_permutation"][
                "repetitions"
            ]
        ),
        "cross_target": int(
            protocol["nulls"]["cross_target_guide_row_permutation"][
                "repetitions"
            ]
        ),
    }
    observed_values = {
        "within_target": np.asarray(within_target_values, dtype=np.float64),
        "cross_target": np.asarray(cross_target_values, dtype=np.float64),
    }
    observed_permutations = {
        "within_target": np.asarray(within_target_permutations, dtype=np.int64),
        "cross_target": np.asarray(cross_target_permutations, dtype=np.int64),
    }
    required_cache = {"features", "guide_targets", "identifiable"}
    if not required_cache.issubset(cache):
        raise ValueError("source cache lacks arrays required for null verification")
    features = np.asarray(cache["features"])
    guide_targets = tuple(map(str, cache["guide_targets"]))
    identifiable = np.asarray(cache["identifiable"], dtype=bool)
    if features.ndim != 2 or len(features) != len(guide_targets):
        raise ValueError("source cache arrays do not align for null verification")
    target_strata, expected_strata_audit = build_split_local_cross_target_strata(
        features,
        guide_targets,
        identifiable,
        splits,
    )
    development_targets = tuple(
        sorted(target for name in REQUIRED_SPLITS for target in splits[name])
    )
    if len(development_targets) != len(set(development_targets)):
        raise ValueError("null verification development target axis overlaps")
    labels = np.asarray(guide_targets, dtype=object)
    grouped = model.EXP083.target_rows(guide_targets)
    for family in ("within_target", "cross_target"):
        records = payload[family]
        if not isinstance(records, list) or len(records) != expected_repetitions[family]:
            raise ValueError(f"{family} null selection repetition count changed")
        if observed_values[family].shape != (expected_repetitions[family],):
            raise ValueError(f"{family} null metric array shape changed")
        if observed_permutations[family].shape != (
            expected_repetitions[family],
            len(guide_targets),
        ):
            raise ValueError(f"{family} null assignment array shape changed")
        for replicate, record in enumerate(records):
            if int(record.get("replicate", -1)) != replicate:
                raise ValueError(f"{family} null replicate order changed")
            expected_seed = model.EXP083.stable_seed(
                model.SPLIT_SEED,
                f"exp109-{family}-null",
                replicate,
            )
            if int(record.get("seed", -1)) != expected_seed:
                raise ValueError(f"{family} null seed changed")
            expected_permutation = (
                within_target_permutation_indices(
                    guide_targets,
                    seed=expected_seed,
                    row_mask=identifiable,
                    targets=development_targets,
                )
                if family == "within_target"
                else cross_target_permutation_indices(
                    guide_targets,
                    target_strata,
                    seed=expected_seed,
                    row_mask=identifiable,
                )
            )
            np.testing.assert_array_equal(
                observed_permutations[family][replicate],
                expected_permutation,
            )
            if sha256_array(expected_permutation) != record.get(
                "permutation_indices_sha256"
            ):
                raise ValueError(f"{family} null assignment hash changed")
            if sha256_array(features[expected_permutation]) != record.get(
                "permuted_features_sha256"
            ):
                raise ValueError(f"{family} null feature hash changed")
            candidate_records = record.get("candidate_records")
            if not isinstance(candidate_records, list):
                raise TypeError(f"{family} null candidate records are missing")
            parsed_records = [
                (
                    candidate_from_payload(item["candidate"]),
                    dict(item["metrics"]),
                )
                for item in candidate_records
            ]
            candidate, selected_metrics = select_candidate(parsed_records)
            if candidate_from_payload(record["selected_candidate"]) != candidate:
                raise ValueError(f"{family} null selected candidate was not reproduced")
            claimed_selected_metrics = record.get("selected_tune_metrics")
            if not isinstance(claimed_selected_metrics, Mapping):
                raise TypeError(f"{family} null selected tune metrics are missing")
            if dict(claimed_selected_metrics) != selected_metrics:
                raise ValueError(f"{family} null selected tune metrics changed")
            metrics = record.get("source_test_metrics")
            if not isinstance(metrics, Mapping):
                raise TypeError(f"{family} null source-test metrics are missing")
            if not np.isclose(
                float(metrics["cosine_delta"]),
                float(observed_values[family][replicate]),
                rtol=0.0,
                atol=1e-15,
            ):
                raise ValueError(f"{family} null metric differs from frozen array")
            if family == "within_target":
                for target in splits["test"]:
                    rows = grouped[str(target)][identifiable[grouped[str(target)]]]
                    if np.any(expected_permutation[rows] == rows):
                        raise ValueError(
                            "within-target null left an eligible test guide fixed"
                        )
            else:
                for target in splits["test"]:
                    rows = grouped[str(target)][identifiable[grouped[str(target)]]]
                    if np.any(labels[expected_permutation[rows]] == str(target)):
                        raise ValueError(
                            "cross-target null left an eligible test target fixed"
                        )
    strata = payload["cross_target_strata"]
    if not isinstance(strata, Mapping) or dict(strata) != expected_strata_audit:
        raise ValueError("cross-target null strata audit is missing")
    return payload


def candidate_payload(candidate: model.Candidate) -> dict[str, float | int]:
    return {
        "response_rank": candidate.response_rank,
        "ridge": candidate.ridge,
        "eta": candidate.eta,
        "weight_power": candidate.weight_power,
        "huber_c": candidate.huber_c,
    }


def candidate_from_payload(payload: Mapping[str, object]) -> model.Candidate:
    return model.Candidate(
        response_rank=int(payload["response_rank"]),
        ridge=float(payload["ridge"]),
        eta=float(payload["eta"]),
        weight_power=float(payload["weight_power"]),
        huber_c=float(payload["huber_c"]),
    )


def validate_tune_selection(path: Path) -> model.Candidate:
    payload = json.loads(path.read_text())
    records = payload.get("candidate_records")
    if not isinstance(records, list):
        raise TypeError("tune selection has no candidate records")
    parsed = [
        (
            candidate_from_payload(record["candidate"]),
            dict(record["metrics"]),
        )
        for record in records
    ]
    selected, _ = select_candidate(parsed)
    claimed = candidate_from_payload(payload["selected_candidate"])
    if claimed != selected:
        raise ValueError("tune selection does not match the complete candidate grid")
    return selected


def tune_selection_payload(
    records: Sequence[tuple[model.Candidate, Mapping[str, float]]],
    selected: model.Candidate,
) -> dict[str, object]:
    """Canonical serialized evidence for a complete recomputed tune grid."""
    return {
        "candidate_records": [
            {
                "candidate": candidate_payload(candidate),
                "metrics": {
                    name: (
                        int(value)
                        if name in {"targets", "guide_holdouts"}
                        else float(value)
                    )
                    for name, value in metrics.items()
                },
            }
            for candidate, metrics in records
        ],
        "selected_candidate": candidate_payload(selected),
    }


def assert_models_identical(
    observed: model.StableSubspaceModel,
    expected: model.StableSubspaceModel,
    *,
    label: str,
) -> None:
    """Reject a serialized model not exactly regenerated from source arrays."""
    if observed.candidate != expected.candidate:
        raise ValueError(f"{label} candidate differs from regenerated model")
    scalar_fields = (
        "gene_axis_sha256",
        "panel_targets_sha256",
        "development_targets_sha256",
        "response_axis_sha256",
        "fit_targets_sha256",
    )
    if any(getattr(observed, field) != getattr(expected, field) for field in scalar_fields):
        raise ValueError(f"{label} source-axis binding differs from regenerated model")
    array_fields = (
        "feature_mean",
        "feature_scale",
        "feature_keep",
        "stable_feature_basis",
        "stable_singular_values",
        "response_mask",
        "response_basis",
        "coefficients",
    )
    for field in array_fields:
        np.testing.assert_array_equal(
            getattr(observed, field),
            getattr(expected, field),
            err_msg=f"{label} differs from regenerated source fit: {field}",
        )


def canonical_response_axis() -> model.AuthenticatedResponseAxis:
    source_axes = load_canonical_source_axes()
    genes, gene_contract = canonical_axis("gene_axis", source_axes)
    panel, panel_contract = canonical_axis("panel_targets", source_axes)
    development, development_contract = canonical_axis(
        "development_targets",
        source_axes,
    )
    _, response_contract = canonical_axis("response_axis", source_axes)
    return model.authenticate_response_axis(
        genes,
        panel,
        development,
        gene_axis_contract=gene_contract,
        panel_targets_contract=panel_contract,
        development_targets_contract=development_contract,
        response_axis_contract=response_contract,
    )


def validate_frozen_model_structure(
    fitted: model.StableSubspaceModel,
    *,
    response_axis: model.AuthenticatedResponseAxis,
    expected_fit_targets_sha256: str,
) -> None:
    """Validate numerical dimensions, axes, and fit-population identity."""
    if fitted.candidate not in set(model.candidate_grid()):
        raise ValueError("frozen model candidate is outside the preregistered grid")
    if fitted.candidate.weight_power != 0.0:
        raise ValueError("frozen model changed the fixed guide-cell weight")
    if (
        fitted.gene_axis_sha256 != response_axis.gene_axis_contract.sha256
        or fitted.panel_targets_sha256 != response_axis.panel_targets_contract.sha256
        or fitted.development_targets_sha256
        != response_axis.development_targets_contract.sha256
        or fitted.response_axis_sha256 != response_axis.response_axis_contract.sha256
    ):
        raise ValueError("frozen model is bound to different response axes")
    if fitted.fit_targets_sha256 != expected_fit_targets_sha256:
        raise ValueError("frozen model fit-target axis changed")
    if not np.array_equal(fitted.response_mask, response_axis.response_mask):
        raise ValueError("frozen model response mask differs from authenticated axes")
    if fitted.feature_keep.shape != (229,):
        raise ValueError("frozen feature mask width changed")
    stable_rank = fitted.stable_rank
    response_rank = fitted.candidate.response_rank
    kept_features = int(np.sum(fitted.feature_keep))
    if fitted.feature_mean.shape != (kept_features,):
        raise ValueError("frozen feature mean width changed")
    if fitted.feature_scale.shape != (kept_features,):
        raise ValueError("frozen feature scale width changed")
    if np.any(fitted.feature_scale <= 0):
        raise ValueError("frozen feature scales must be positive")
    if fitted.stable_feature_basis.shape != (kept_features, stable_rank):
        raise ValueError("frozen stable feature basis shape changed")
    if fitted.stable_singular_values.shape != (stable_rank,):
        raise ValueError("frozen stable singular-value count changed")
    if np.any(fitted.stable_singular_values <= 0):
        raise ValueError("frozen stable singular values must be positive")
    if fitted.response_basis.shape != (
        response_rank,
        int(response_axis.response_mask.sum()),
    ):
        raise ValueError("frozen model response basis width changed")
    if fitted.coefficients.shape != (stable_rank, response_rank):
        raise ValueError("frozen model response-rank dimensions differ")
    if stable_rank:
        np.testing.assert_allclose(
            fitted.stable_feature_basis.T @ fitted.stable_feature_basis,
            np.eye(stable_rank),
            atol=1e-8,
            rtol=1e-8,
        )
    if response_rank:
        np.testing.assert_allclose(
            fitted.response_basis @ fitted.response_basis.T,
            np.eye(response_rank),
            atol=1e-8,
            rtol=1e-8,
        )
    arrays = (
        fitted.feature_mean,
        fitted.feature_scale,
        fitted.stable_feature_basis,
        fitted.stable_singular_values,
        fitted.response_basis,
        fitted.coefficients,
    )
    if not all(np.isfinite(value).all() for value in arrays):
        raise ValueError("frozen model contains non-finite values")


def _npz_scalar_string(archive: Mapping[str, np.ndarray], name: str) -> str:
    value = np.asarray(archive[name])
    if value.shape != ():
        raise ValueError(f"{name} must be a scalar string")
    return str(value.item())


def validate_source_test_holdouts(
    path: Path,
    *,
    expected_test_targets: Sequence[str],
    expected_guide_targets: Sequence[str],
    expected_model_sha256: str,
    expected_source_input_sha256: str,
    expected_response_axis_sha256: str,
    expected_response_width: int,
    recomputed: Mapping[str, object] | None = None,
) -> dict[str, object]:
    """Validate frozen per-guide outputs and recompute all source metrics."""
    expected_targets = tuple(sorted(map(str, expected_test_targets)))
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != {
            "schema_version",
            "raw_residual",
            "recovered_nuisance",
            "shrink_recovered_nuisance",
            "target_labels",
            "held_guide_rows",
            "source_test_model_sha256",
            "source_input_manifest_sha256",
            "response_axis_sha256",
        }:
            raise ValueError("source-test holdout field set is invalid")
        if int(np.asarray(archive["schema_version"]).item()) != (
            SOURCE_TEST_HOLDOUT_SCHEMA_VERSION
        ):
            raise ValueError("unsupported source-test holdout schema")
        truth = np.asarray(archive["raw_residual"], dtype=np.float64)
        predicted = np.asarray(archive["recovered_nuisance"], dtype=np.float64)
        shrink_predicted = np.asarray(
            archive["shrink_recovered_nuisance"],
            dtype=np.float64,
        )
        labels = tuple(map(str, np.asarray(archive["target_labels"]).tolist()))
        held_rows = np.asarray(archive["held_guide_rows"], dtype=np.int64)
        if _npz_scalar_string(archive, "source_test_model_sha256") != (
            expected_model_sha256
        ):
            raise ValueError("source-test holdouts name a different test model")
        if _npz_scalar_string(archive, "source_input_manifest_sha256") != (
            expected_source_input_sha256
        ):
            raise ValueError("source-test holdouts name different source inputs")
        if _npz_scalar_string(archive, "response_axis_sha256") != (
            expected_response_axis_sha256
        ):
            raise ValueError("source-test holdouts use a different response axis")
    if (
        truth.ndim != 2
        or predicted.shape != truth.shape
        or shrink_predicted.shape != truth.shape
    ):
        raise ValueError("source-test holdout matrix shapes differ")
    if truth.shape[1] != expected_response_width:
        raise ValueError("source-test holdout response width changed")
    if len(labels) != len(truth) or held_rows.shape != (len(truth),):
        raise ValueError("source-test holdout row metadata differs")
    if not all(
        np.isfinite(values).all() for values in (truth, predicted, shrink_predicted)
    ):
        raise ValueError("source-test holdouts contain non-finite values")
    if np.any(held_rows < 0) or len(set(held_rows.tolist())) != len(held_rows):
        raise ValueError("held guide rows must be unique nonnegative indices")
    guide_targets = tuple(map(str, expected_guide_targets))
    if np.any(held_rows >= len(guide_targets)):
        raise ValueError("held guide row lies outside the authenticated guide axis")
    if tuple(guide_targets[index] for index in held_rows) != labels:
        raise ValueError("held guide rows do not identify their recorded targets")
    if tuple(sorted(set(labels))) != expected_targets:
        raise ValueError("source-test holdouts do not cover the canonical test targets")
    counts = {target: labels.count(target) for target in expected_targets}
    if any(count < 2 for count in counts.values()):
        raise ValueError("each source-test target requires at least two guide holdouts")
    baseline = np.zeros_like(truth)
    metrics = model.target_macro_metrics(truth, baseline, predicted, labels)
    shrink_metrics = model.target_macro_metrics(
        truth,
        baseline,
        shrink_predicted,
        labels,
    )
    cosine_delta, nmse_delta = model.per_guide_metric_deltas(truth, predicted)
    observed = {
        "truth": truth,
        "predicted": predicted,
        "shrink_predicted": shrink_predicted,
        "labels": labels,
        "held_guide_rows": held_rows,
        "metrics": metrics,
        "shrink_metrics": shrink_metrics,
        "per_guide_cosine_delta": cosine_delta,
        "per_guide_nmse_delta": nmse_delta,
    }
    if recomputed is not None:
        required = {
            "truth",
            "predicted",
            "shrink_predicted",
            "labels",
            "held_guide_rows",
        }
        if not required.issubset(recomputed):
            raise ValueError("recomputed source-test evidence is incomplete")
        for name in ("truth", "predicted", "shrink_predicted", "held_guide_rows"):
            np.testing.assert_array_equal(observed[name], recomputed[name])
        if tuple(recomputed["labels"]) != labels:
            raise ValueError("held source-test labels were not recomputed")
    return observed


def deterministic_stratified_permutations(
    strata: np.ndarray,
    movable: np.ndarray,
    *,
    repetitions: int,
    seed: int,
) -> np.ndarray:
    """Freeze target-label permutations while preserving deployment strata."""
    strata = np.asarray(strata, dtype=np.int64)
    movable = np.asarray(movable, dtype=bool)
    if strata.ndim != 1 or movable.shape != strata.shape:
        raise ValueError("held-label permutation strata shape mismatch")
    if repetitions < 1:
        raise ValueError("held-label permutation repetitions must be positive")
    result = np.tile(np.arange(len(strata), dtype=np.int64), (repetitions, 1))
    groups = [
        np.flatnonzero(movable & (strata == value))
        for value in sorted(set(strata[movable].tolist()))
    ]
    if not groups or any(len(group) < 2 for group in groups):
        raise ValueError("every movable permutation stratum needs at least two rows")
    rng = np.random.default_rng(seed)
    identity = np.arange(len(strata), dtype=np.int64)
    for repeat in range(repetitions):
        for group in groups:
            result[repeat, group] = rng.permutation(group)
        if np.array_equal(result[repeat], identity):
            group = max(groups, key=len)
            result[repeat, group] = np.roll(group, 1)
    return result


def validate_singleton_corrections(
    path: Path,
    *,
    fitted: model.StableSubspaceModel,
    expected_model_sha256: str,
    expected_source_test_model_sha256: str,
    expected_source_input_sha256: str,
) -> dict[str, object]:
    """Authenticate singleton features and reproduce deployed panel corrections."""
    source_axes = load_canonical_source_axes()
    singleton_targets, singleton_contract = canonical_axis(
        "singleton_centering_targets",
        source_axes,
    )
    panel_targets, _ = canonical_axis("panel_targets", source_axes)
    genes, _ = canonical_axis("gene_axis", source_axes)
    pm_targets, _ = canonical_axis("perturbmatch_targets", source_axes)
    pm_genes, _ = canonical_axis("perturbmatch_genes", source_axes)
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != {
            "schema_version",
            "singleton_target_ids",
            "singleton_features",
            "panel_singleton_features",
            "development_features",
            "development_target_ids",
            "panel_target_ids",
            "panel_singleton_source_index",
            "response_corrections",
            "source_slot_corrections",
            "source_test_model_sha256",
            "selected_model_sha256",
            "source_input_manifest_sha256",
            "response_axis_sha256",
            "centering_axis_sha256",
            "extrapolation_metric_names",
            "extrapolation_metric_values",
        }:
            raise ValueError("singleton-correction field set is invalid")
        if int(np.asarray(archive["schema_version"]).item()) != (
            SINGLETON_CORRECTION_SCHEMA_VERSION
        ):
            raise ValueError("unsupported singleton-correction schema")
        singleton_ids = tuple(
            map(str, np.asarray(archive["singleton_target_ids"]).tolist())
        )
        singleton_features_raw = np.asarray(archive["singleton_features"])
        panel_singleton_features_raw = np.asarray(
            archive["panel_singleton_features"]
        )
        development_features_raw = np.asarray(archive["development_features"])
        development_ids_raw = np.asarray(archive["development_target_ids"])
        singleton_features = np.asarray(singleton_features_raw, dtype=np.float64)
        development_features = np.asarray(
            development_features_raw,
            dtype=np.float64,
        )
        development_ids = tuple(map(str, development_ids_raw.tolist()))
        panel_ids = tuple(map(str, np.asarray(archive["panel_target_ids"]).tolist()))
        panel_index = np.asarray(
            archive["panel_singleton_source_index"],
            dtype=np.int64,
        )
        response_corrections = np.asarray(
            archive["response_corrections"],
            dtype=np.float64,
        )
        slot_corrections = np.asarray(
            archive["source_slot_corrections"],
            dtype=np.float64,
        )
        metric_names = tuple(
            map(str, np.asarray(archive["extrapolation_metric_names"]).tolist())
        )
        metric_values = np.asarray(
            archive["extrapolation_metric_values"],
            dtype=np.float64,
        )
        scalar_hashes = {
            name: _npz_scalar_string(archive, name)
            for name in (
                "source_test_model_sha256",
                "selected_model_sha256",
                "source_input_manifest_sha256",
                "response_axis_sha256",
                "centering_axis_sha256",
            )
        }
    if scalar_hashes["selected_model_sha256"] != expected_model_sha256:
        raise ValueError("singleton corrections name a different deployment model")
    if scalar_hashes["source_test_model_sha256"] != expected_source_test_model_sha256:
        raise ValueError("singleton corrections name a different source-test model")
    if scalar_hashes["source_input_manifest_sha256"] != expected_source_input_sha256:
        raise ValueError("singleton corrections name different source inputs")
    if scalar_hashes["response_axis_sha256"] != fitted.response_axis_sha256:
        raise ValueError("singleton corrections use a different response axis")
    if scalar_hashes["centering_axis_sha256"] != singleton_contract.sha256:
        raise ValueError("singleton corrections use a different centering axis")
    if not is_sha256(scalar_hashes["source_test_model_sha256"]):
        raise ValueError("invalid source-test model binding")
    if singleton_ids != singleton_targets:
        raise ValueError("singleton correction target axis is not canonical")
    if panel_ids != panel_targets:
        raise ValueError("panel correction target axis is not canonical")
    if singleton_features.shape != (len(singleton_targets), len(fitted.feature_keep)):
        raise ValueError("singleton feature matrix shape changed")
    if development_features.ndim != 2 or development_features.shape[1] != len(
        fitted.feature_keep
    ):
        raise ValueError("development feature matrix shape changed")
    if len(development_ids) != len(development_features):
        raise ValueError("development feature labels differ from matrix rows")
    if set(development_ids) != set(
        canonical_axis("development_targets", source_axes)[0]
    ):
        raise ValueError("development feature rows do not cover canonical targets")
    development_counts = {
        target: development_ids.count(target) for target in set(development_ids)
    }
    if any(count < 2 for count in development_counts.values()):
        raise ValueError("development features require at least two rows per target")
    if panel_index.shape != (len(panel_targets),):
        raise ValueError("panel singleton-index shape changed")
    valid = panel_index >= 0
    if int(np.sum(valid)) != len(set(panel_targets) & set(singleton_targets)):
        raise ValueError("panel singleton support count changed")
    if np.any(panel_index[valid] >= len(singleton_targets)):
        raise ValueError("panel singleton index is out of bounds")
    for panel_row in np.flatnonzero(valid):
        if singleton_targets[int(panel_index[panel_row])] != panel_targets[panel_row]:
            raise ValueError("panel singleton index maps to the wrong target")
    input_manifest_path = path.parent / "source-input-manifest.json"
    if sha256_file(input_manifest_path) != expected_source_input_sha256:
        raise ValueError("source-input manifest changed before correction validation")
    input_manifest = validate_source_input_manifest(input_manifest_path)
    derived = input_manifest["derived_arrays"]
    observed_derived = {
        "singleton_features": singleton_features_raw,
        "panel_singleton_features": panel_singleton_features_raw,
        "development_features": development_features_raw,
        "development_target_ids": development_ids_raw,
        "panel_singleton_source_index": panel_index,
    }
    for name, values in observed_derived.items():
        contract = derived[name]
        if list(values.shape) != contract["shape"]:
            raise ValueError(f"derived-array shape differs: {name}")
        if str(values.dtype) != contract["dtype"]:
            raise ValueError(f"derived-array dtype differs: {name}")
        if sha256_array(values) != contract["sha256"]:
            raise ValueError(f"derived-array hash differs: {name}")
    np.testing.assert_array_equal(
        panel_singleton_features_raw,
        singleton_features_raw[panel_index[valid]],
    )
    expected_response = model.centered_singleton_corrections(
        fitted.nuisance(singleton_features),
        singleton_targets,
        np.ones(len(singleton_targets), dtype=bool),
        eta=fitted.candidate.eta,
        centering_axis_contract=singleton_contract,
    )
    expected_panel = np.zeros(
        (len(panel_targets), fitted.response_width),
        dtype=np.float64,
    )
    expected_panel[valid] = expected_response[panel_index[valid]]
    if response_corrections.shape != expected_panel.shape:
        raise ValueError("panel response-correction shape changed")
    assert_finite_allclose(
        response_corrections,
        expected_panel,
        label="panel response corrections",
        rtol=0.0,
        atol=1e-12,
    )
    full = model.expand_response_corrections(fitted, expected_panel)
    expected_slot = model.apply_k562_slot_retention(
        full,
        np.asarray([target in set(pm_targets) for target in panel_targets]),
        np.asarray([gene in set(pm_genes) for gene in genes]),
    )
    if slot_corrections.shape != expected_slot.shape:
        raise ValueError("source-slot correction shape changed")
    assert_finite_allclose(
        slot_corrections,
        expected_slot,
        label="source-slot corrections",
        rtol=0.0,
        atol=1e-12,
    )
    np.testing.assert_array_equal(response_corrections[~valid], 0.0)
    np.testing.assert_array_equal(slot_corrections[:, ~fitted.response_mask], 0.0)
    source_mean_abs_max = float(np.max(np.abs(expected_response.sum(axis=0))))
    if source_mean_abs_max > 1e-12:
        raise ValueError("pre-PerturbMatch singleton correction is not zero mean")
    probe_features = singleton_features[: min(64, len(singleton_features))]
    projected = fitted.projected_standardized_features(probe_features)
    projected_features = probe_features.copy()
    projected_features[:, fitted.feature_keep] = (
        fitted.feature_mean + fitted.feature_scale * projected
    )
    nullspace_prediction_change = float(
        np.max(
            np.abs(
                fitted.nuisance(probe_features) - fitted.nuisance(projected_features)
            )
        )
    )
    nullspace_limit = float(
        validate_preregistration()["promotion_gates"][
            "feature_nullspace_prediction_change_abs_max"
        ]
    )
    if nullspace_prediction_change > nullspace_limit:
        raise ValueError("feature-nullspace perturbation changed nuisance prediction")

    diagnostics = model.stable_score_extrapolation_diagnostics(
        fitted,
        development_features,
        singleton_features,
    )
    if metric_names != tuple(diagnostics):
        raise ValueError("extrapolation diagnostic field order changed")
    np.testing.assert_allclose(
        metric_values,
        np.asarray(list(diagnostics.values()), dtype=np.float64),
        rtol=0.0,
        atol=1e-12,
    )
    final_invariants = validate_preregistered_final_invariants(
        fitted,
        source_derived=observed_derived,
        source_input_manifest=input_manifest,
    )
    return {
        "supported_panel_targets": int(np.sum(valid)),
        "source_mean_abs_max_pre_perturbmatch": source_mean_abs_max,
        "source_mean_abs_max_post_perturbmatch": float(
            np.max(np.abs(expected_slot.sum(axis=0)))
        ),
        "feature_nullspace_prediction_change_abs_max": nullspace_prediction_change,
        "extrapolation": diagnostics,
        "final_model_invariants": final_invariants,
        "panel_support": valid,
    }


def validate_held_label_nulls(
    path: Path,
    *,
    panel_support: np.ndarray,
    singleton_corrections_sha256: str,
) -> np.ndarray:
    """Require exact deterministic, support-preserving label permutations."""
    source_axes = load_canonical_source_axes()
    panel_targets, _ = canonical_axis("panel_targets", source_axes)
    pm_targets, _ = canonical_axis("perturbmatch_targets", source_axes)
    support = np.asarray(panel_support, dtype=bool)
    if support.shape != (len(panel_targets),):
        raise ValueError("held-label support shape mismatch")
    strata = np.zeros(len(panel_targets), dtype=np.int64)
    strata[support] = 1
    strata[
        support
        & np.asarray([target not in set(pm_targets) for target in panel_targets])
    ] = 2
    expected = deterministic_stratified_permutations(
        strata,
        support,
        repetitions=int(
            validate_preregistration()["nulls"][
                "held_singleton_correction_label_permutation"
            ]["repetitions"]
        ),
        seed=HELD_LABEL_NULL_SEED,
    )
    with np.load(path, allow_pickle=False) as archive:
        if set(archive.files) != {
            "schema_version",
            "seed",
            "permutations",
            "strata",
            "movable",
            "panel_targets_sha256",
            "singleton_corrections_sha256",
        }:
            raise ValueError("held-label-null field set is invalid")
        if int(np.asarray(archive["schema_version"]).item()) != (
            HELD_LABEL_NULL_SCHEMA_VERSION
        ):
            raise ValueError("unsupported held-label-null schema")
        if int(np.asarray(archive["seed"]).item()) != HELD_LABEL_NULL_SEED:
            raise ValueError("held-label-null seed changed")
        observed = np.asarray(archive["permutations"], dtype=np.int64)
        observed_strata = np.asarray(archive["strata"], dtype=np.int64)
        observed_movable = np.asarray(archive["movable"], dtype=bool)
        if _npz_scalar_string(archive, "panel_targets_sha256") != sha256_lines(
            panel_targets
        ):
            raise ValueError("held-label nulls use a different panel axis")
        if (
            _npz_scalar_string(
                archive,
                "singleton_corrections_sha256",
            )
            != singleton_corrections_sha256
        ):
            raise ValueError("held-label nulls name different singleton corrections")
    np.testing.assert_array_equal(observed_strata, strata)
    np.testing.assert_array_equal(observed_movable, support)
    np.testing.assert_array_equal(observed, expected)
    identity = np.arange(len(panel_targets), dtype=np.int64)
    for permutation in observed:
        if not np.array_equal(np.sort(permutation), identity):
            raise ValueError("held-label null row is not a bijection")
        if np.array_equal(permutation, identity):
            raise ValueError("held-label null contains the identity")
        if np.any(permutation[~support] != identity[~support]):
            raise ValueError("held-label null moves unsupported targets")
        if np.any(strata[permutation] != strata):
            raise ValueError("held-label null crosses deployment strata")
    return observed


def validate_source_evidence(
    output: Path,
    source_test: Mapping[str, object],
    splits: Mapping[str, Sequence[str]],
) -> model.StableSubspaceModel:
    """Recompute source gates and bind every source artifact to one model."""
    protocol = validate_preregistration()
    source_input_path = output / "source-input-manifest.json"
    source_input = validate_source_input_manifest(source_input_path)
    source_input_sha256 = sha256_file(source_input_path)
    cache, _, _ = validate_source_cache_snapshot(
        output,
        expected_source_manifest=source_input,
    )
    expected_guide_targets = tuple(map(str, cache["guide_targets"]))
    snapshot_sha256 = sha256_file(output / "source-cache-snapshot.json")
    validate_memory_guard_receipt(
        output / "memory-guard.json",
        expected_source_manifest_content_sha256=sha256_json(
            source_manifest_content(source_input)
        ),
        expected_source_cache_snapshot_sha256=snapshot_sha256,
    )
    selected = validate_tune_selection(output / "tune-selection.json")
    if selected.is_identity:
        raise ValueError("identity candidate cannot authorize held evaluation")
    source_test_model_path = output / "source-test-model.npz"
    source_test_model = load_frozen_model(
        source_test_model_path,
        expected_role="source_test_model",
    )
    fitted = load_frozen_model(
        output / "selected-model.npz",
        expected_role="deployment_model",
    )
    if fitted.candidate != selected or source_test_model.candidate != selected:
        raise ValueError("selected model differs from tune-selection winner")
    if candidate_from_payload(source_test["selected_candidate"]) != selected:
        raise ValueError("source test names a different selected candidate")
    fit_tune_targets = tuple(sorted((*splits["fit"], *splits["tune"])))
    response_axis = canonical_response_axis()
    validate_frozen_model_structure(
        source_test_model,
        response_axis=response_axis,
        expected_fit_targets_sha256=sha256_lines(fit_tune_targets),
    )
    validate_frozen_model_structure(
        fitted,
        response_axis=response_axis,
        expected_fit_targets_sha256=response_axis.development_targets_contract.sha256,
    )
    if set(fit_tune_targets) & set(splits["test"]):
        raise ValueError("source-test model training targets overlap test targets")

    holdouts = validate_source_test_holdouts(
        output / "source-test-holdouts.npz",
        expected_test_targets=splits["test"],
        expected_guide_targets=expected_guide_targets,
        expected_model_sha256=sha256_file(source_test_model_path),
        expected_source_input_sha256=source_input_sha256,
        expected_response_axis_sha256=source_test_model.response_axis_sha256,
        expected_response_width=source_test_model.response_width,
    )
    metrics = holdouts["metrics"]
    claimed_metrics = source_test.get("metrics")
    if not isinstance(claimed_metrics, Mapping):
        raise TypeError("source test metrics are missing")
    for name in ("targets", "guide_holdouts"):
        if int(claimed_metrics[name]) != int(metrics[name]):
            raise ValueError(f"source-test metric was not reproduced: {name}")
    for name in ("cosine_delta", "nmse_delta"):
        if not np.isclose(
            float(claimed_metrics[name]),
            float(metrics[name]),
            rtol=0.0,
            atol=1e-15,
        ):
            raise ValueError(f"source-test metric was not reproduced: {name}")

    with np.load(output / "source-bootstrap.npz", allow_pickle=False) as archive:
        if set(archive.files) != {
            "schema_version",
            "seed",
            "cosine_delta",
            "source_test_holdouts_sha256",
        }:
            raise ValueError("source bootstrap field set is invalid")
        if int(np.asarray(archive["schema_version"]).item()) != 1:
            raise ValueError("unsupported source bootstrap schema")
        if int(np.asarray(archive["seed"]).item()) != model.SPLIT_SEED:
            raise ValueError("source bootstrap seed changed")
        bootstrap = np.asarray(archive["cosine_delta"], dtype=np.float64)
        if _npz_scalar_string(
            archive,
            "source_test_holdouts_sha256",
        ) != sha256_file(output / "source-test-holdouts.npz"):
            raise ValueError("source bootstrap names different holdout outputs")
    with np.load(output / "source-nulls.npz", allow_pickle=False) as archive:
        if set(archive.files) != {
            "schema_version",
            "within_target_cosine_delta",
            "cross_target_cosine_delta",
            "within_target_permutation_indices",
            "cross_target_permutation_indices",
            "direction_preserving_shrink_delta",
        }:
            raise ValueError("source null field set is invalid")
        if int(np.asarray(archive["schema_version"]).item()) != 2:
            raise ValueError("unsupported source null schema")
        within = np.asarray(archive["within_target_cosine_delta"], dtype=np.float64)
        cross = np.asarray(archive["cross_target_cosine_delta"], dtype=np.float64)
        within_permutations = np.asarray(
            archive["within_target_permutation_indices"],
            dtype=np.int64,
        )
        cross_permutations = np.asarray(
            archive["cross_target_permutation_indices"],
            dtype=np.int64,
        )
        shrink_values = np.asarray(
            archive["direction_preserving_shrink_delta"], dtype=np.float64
        )
    repetitions = int(
        protocol["source_validation"]["test_target_cluster_bootstrap_repetitions"]
    )
    expected_bootstrap = target_cluster_bootstrap(
        holdouts["per_guide_cosine_delta"],
        holdouts["labels"],
        repetitions=repetitions,
        seed=model.SPLIT_SEED,
    )
    if len(bootstrap) != repetitions:
        raise ValueError("source bootstrap repetition count changed")
    np.testing.assert_array_equal(bootstrap, expected_bootstrap)
    null_repetitions = int(
        protocol["nulls"]["within_target_guide_feature_permutation"]["repetitions"]
    )
    if len(within) != null_repetitions or len(cross) != int(
        protocol["nulls"]["cross_target_guide_row_permutation"]["repetitions"]
    ):
        raise ValueError("source null repetition count changed")
    validate_source_null_selection(
        output / "source-null-selection.json",
        cache=cache,
        splits=splits,
        within_target_values=within,
        cross_target_values=cross,
        within_target_permutations=within_permutations,
        cross_target_permutations=cross_permutations,
    )
    if shrink_values.shape != (1,):
        raise ValueError("direction-preserving shrink artifact must be scalar")
    if not np.isclose(
        float(shrink_values[0]),
        float(holdouts["shrink_metrics"]["cosine_delta"]),
        rtol=0.0,
        atol=1e-15,
    ):
        raise ValueError("direction-preserving shrink metric was not reproduced")
    recomputed = source_test_decision(
        metrics,
        bootstrap_samples=bootstrap,
        within_target_null=within,
        cross_target_null=cross,
        direction_preserving_shrink_delta=float(shrink_values[0]),
    )
    for key in (
        "passes",
        "gates",
        "metrics",
        "bootstrap_lower_95",
        "within_target_null_p95",
        "cross_target_null_p95",
        "direction_preserving_shrink_delta",
    ):
        if recomputed[key] != source_test.get(key):
            raise ValueError(f"source-test field was not reproduced: {key}")
    validate_source_test(source_test)

    bindings = source_test.get("bindings")
    if not isinstance(bindings, Mapping):
        raise TypeError("source test has no artifact bindings")
    expected_bindings = {
        "source_input_manifest_sha256": source_input_sha256,
        "source_cache_snapshot_sha256": snapshot_sha256,
        "memory_guard_sha256": sha256_file(output / "memory-guard.json"),
        "source_test_model_sha256": sha256_file(source_test_model_path),
        "source_test_holdouts_sha256": sha256_file(output / "source-test-holdouts.npz"),
        "selected_model_sha256": sha256_file(output / "selected-model.npz"),
        "tune_selection_sha256": sha256_file(output / "tune-selection.json"),
        "source_bootstrap_sha256": sha256_file(output / "source-bootstrap.npz"),
        "source_nulls_sha256": sha256_file(output / "source-nulls.npz"),
        "source_null_selection_sha256": sha256_file(
            output / "source-null-selection.json"
        ),
        "axes_sha256": sha256_file(output / "axes.json"),
        "singleton_corrections_sha256": sha256_file(
            output / "singleton-corrections.npz"
        ),
        "held_label_nulls_sha256": sha256_file(output / "held-label-nulls.npz"),
        "selected_candidate": candidate_payload(selected),
        "split_hashes": {name: sha256_lines(splits[name]) for name in REQUIRED_SPLITS},
    }
    if dict(bindings) != expected_bindings:
        raise ValueError("source-test artifact bindings are inconsistent")
    correction_audit = validate_singleton_corrections(
        output / "singleton-corrections.npz",
        fitted=fitted,
        expected_model_sha256=sha256_file(output / "selected-model.npz"),
        expected_source_test_model_sha256=sha256_file(source_test_model_path),
        expected_source_input_sha256=source_input_sha256,
    )
    if source_test.get("final_model_invariants") != correction_audit[
        "final_model_invariants"
    ]:
        raise ValueError("final-model invariant evidence was not reproduced")
    validate_held_label_nulls(
        output / "held-label-nulls.npz",
        panel_support=correction_audit["panel_support"],
        singleton_corrections_sha256=sha256_file(output / "singleton-corrections.npz"),
    )
    return fitted


def axis_payload(
    response_axis: model.AuthenticatedResponseAxis,
    centering_target_ids: Sequence[str],
    centering_axis_contract: model.AxisContract,
) -> dict[str, object]:
    centering = model.validate_axis(
        centering_target_ids,
        centering_axis_contract,
        label="singleton centering axis",
        require_sorted=True,
    )
    response_genes = tuple(
        gene
        for gene, keep in zip(
            response_axis.gene_axis,
            response_axis.response_mask,
            strict=True,
        )
        if keep
    )
    source_axes = load_canonical_source_axes()
    perturbmatch_targets, _ = canonical_axis("perturbmatch_targets", source_axes)
    perturbmatch_genes, _ = canonical_axis("perturbmatch_genes", source_axes)
    return {
        "centering_targets": list(centering),
        "gene_axis": list(response_axis.gene_axis),
        "panel_targets": list(response_axis.panel_targets),
        "development_targets": list(response_axis.development_targets),
        "response_genes": list(response_genes),
        "perturbmatch_targets": list(perturbmatch_targets),
        "perturbmatch_genes": list(perturbmatch_genes),
    }


def write_axis_payload(
    path: Path,
    response_axis: model.AuthenticatedResponseAxis,
    centering_target_ids: Sequence[str],
    centering_axis_contract: model.AxisContract,
) -> None:
    payload = axis_payload(
        response_axis,
        centering_target_ids,
        centering_axis_contract,
    )
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def save_frozen_model(
    path: Path,
    fitted: model.StableSubspaceModel,
    *,
    role: str = "deployment_model",
) -> None:
    """Serialize only numerical model state; no pickle is used."""
    if role not in {"source_test_model", "deployment_model"}:
        raise ValueError("invalid frozen model role")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        np.savez_compressed(
            handle,
            schema_version=np.asarray(1),
            model_role=np.asarray(role),
            candidate=np.asarray(
                [
                    fitted.candidate.response_rank,
                    fitted.candidate.ridge,
                    fitted.candidate.eta,
                    fitted.candidate.weight_power,
                    fitted.candidate.huber_c,
                ],
                dtype=np.float64,
            ),
            feature_mean=fitted.feature_mean,
            feature_scale=fitted.feature_scale,
            feature_keep=fitted.feature_keep,
            stable_feature_basis=fitted.stable_feature_basis,
            stable_singular_values=fitted.stable_singular_values,
            response_mask=fitted.response_mask,
            response_axis_hashes=np.asarray(
                [
                    fitted.gene_axis_sha256,
                    fitted.panel_targets_sha256,
                    fitted.development_targets_sha256,
                    fitted.response_axis_sha256,
                    fitted.fit_targets_sha256,
                ],
                dtype="<U64",
            ),
            response_basis=fitted.response_basis,
            coefficients=fitted.coefficients,
        )


def load_frozen_model(
    path: Path,
    *,
    expected_role: str | None = None,
) -> model.StableSubspaceModel:
    with np.load(path, allow_pickle=False) as archive:
        expected_fields = {
            "schema_version",
            "model_role",
            "candidate",
            "feature_mean",
            "feature_scale",
            "feature_keep",
            "stable_feature_basis",
            "stable_singular_values",
            "response_mask",
            "response_axis_hashes",
            "response_basis",
            "coefficients",
        }
        if set(archive.files) != expected_fields:
            raise ValueError("frozen model field set is invalid")
        if int(np.asarray(archive["schema_version"]).item()) != 1:
            raise ValueError("unsupported frozen model schema")
        role = _npz_scalar_string(archive, "model_role")
        if role not in {"source_test_model", "deployment_model"}:
            raise ValueError("invalid frozen model role")
        if expected_role is not None and role != expected_role:
            raise ValueError(f"expected {expected_role}, found {role}")
        candidate_values = np.asarray(archive["candidate"], dtype=np.float64)
        if candidate_values.shape != (5,):
            raise ValueError("invalid frozen candidate")
        candidate = model.Candidate(
            response_rank=int(candidate_values[0]),
            ridge=float(candidate_values[1]),
            eta=float(candidate_values[2]),
            weight_power=float(candidate_values[3]),
            huber_c=float(candidate_values[4]),
        )
        response_axis_hashes = np.asarray(archive["response_axis_hashes"])
        if response_axis_hashes.shape != (5,):
            raise ValueError("invalid frozen response-axis hashes")
        if not all(is_sha256(value) for value in response_axis_hashes.tolist()):
            raise ValueError("invalid frozen response-axis hash")
        return model.StableSubspaceModel(
            candidate=candidate,
            feature_mean=np.asarray(archive["feature_mean"], dtype=np.float64),
            feature_scale=np.asarray(archive["feature_scale"], dtype=np.float64),
            feature_keep=np.asarray(archive["feature_keep"], dtype=bool),
            stable_feature_basis=np.asarray(
                archive["stable_feature_basis"], dtype=np.float64
            ),
            stable_singular_values=np.asarray(
                archive["stable_singular_values"], dtype=np.float64
            ),
            response_mask=np.asarray(archive["response_mask"], dtype=bool),
            gene_axis_sha256=str(response_axis_hashes[0]),
            panel_targets_sha256=str(response_axis_hashes[1]),
            development_targets_sha256=str(response_axis_hashes[2]),
            response_axis_sha256=str(response_axis_hashes[3]),
            fit_targets_sha256=str(response_axis_hashes[4]),
            response_basis=np.asarray(archive["response_basis"], dtype=np.float64),
            coefficients=np.asarray(archive["coefficients"], dtype=np.float64),
        )


def build_freeze_manifest(
    output: Path,
    *,
    splits: Mapping[str, Sequence[str]],
    frozen_files: Sequence[str],
    source_test: Mapping[str, object],
    response_axis: model.AuthenticatedResponseAxis,
    centering_target_ids: Sequence[str],
    centering_axis_contract: model.AxisContract,
) -> dict[str, object]:
    """Build the future pre-held authorization manifest."""
    protocol = validate_preregistration()
    if set(splits) != set(REQUIRED_SPLITS):
        raise ValueError("freeze splits are incomplete")
    expected_split_counts = protocol["frozen_source_audit_expectations"]["split_counts"]
    for name in REQUIRED_SPLITS:
        values = tuple(map(str, splits[name]))
        if values != tuple(sorted(values)) or len(set(values)) != len(values):
            raise ValueError(f"{name} split must be sorted and unique")
        if len(values) != int(expected_split_counts[name]):
            raise ValueError(f"{name} split count differs from preregistration")
    model.assert_split_contract(splits, response_axis.panel_targets)
    expected_splits = development_split(
        response_axis.development_targets,
        response_axis.panel_targets,
    )
    if {
        name: tuple(map(str, splits[name])) for name in REQUIRED_SPLITS
    } != expected_splits:
        raise ValueError("freeze splits differ from the canonical deterministic split")

    source_axes = load_canonical_source_axes(protocol)
    canonical_genes, canonical_gene_contract = canonical_axis("gene_axis", source_axes)
    canonical_panel, canonical_panel_contract = canonical_axis(
        "panel_targets", source_axes
    )
    canonical_development, canonical_development_contract = canonical_axis(
        "development_targets", source_axes
    )
    canonical_response, canonical_response_contract = canonical_axis(
        "response_axis", source_axes
    )
    canonical_centering, canonical_centering_contract = canonical_axis(
        "singleton_centering_targets", source_axes
    )
    _, canonical_pm_targets_contract = canonical_axis(
        "perturbmatch_targets", source_axes
    )
    _, canonical_pm_genes_contract = canonical_axis("perturbmatch_genes", source_axes)
    if response_axis.gene_axis != canonical_genes:
        raise ValueError("response gene axis is not canonical")
    if response_axis.panel_targets != canonical_panel:
        raise ValueError("panel target axis is not canonical")
    if response_axis.development_targets != canonical_development:
        raise ValueError("development target axis is not canonical")
    if response_axis.gene_axis_contract != canonical_gene_contract:
        raise ValueError("gene-axis contract is not canonical")
    if response_axis.panel_targets_contract != canonical_panel_contract:
        raise ValueError("panel-axis contract is not canonical")
    if response_axis.development_targets_contract != canonical_development_contract:
        raise ValueError("development-axis contract is not canonical")
    if response_axis.response_axis_contract != canonical_response_contract:
        raise ValueError("response-axis contract is not canonical")
    observed_response = tuple(
        gene
        for gene, keep in zip(
            response_axis.gene_axis,
            response_axis.response_mask,
            strict=True,
        )
        if keep
    )
    if observed_response != canonical_response:
        raise ValueError("derived response axis is not canonical")
    if tuple(map(str, centering_target_ids)) != canonical_centering:
        raise ValueError("singleton centering target axis is not canonical")
    if centering_axis_contract != canonical_centering_contract:
        raise ValueError("singleton centering contract is not canonical")
    files = tuple(map(str, frozen_files))
    if len(files) != len(set(files)) or set(files) != set(REQUIRED_FROZEN_ARTIFACTS):
        raise ValueError("frozen artifact set is incomplete or contains extras")
    expected_axes = axis_payload(
        response_axis,
        centering_target_ids,
        centering_axis_contract,
    )
    axes_path = output / "axes.json"
    if json.loads(axes_path.read_text()) != expected_axes:
        raise ValueError("axes.json differs from authenticated axes")
    source_test_path = output / "source-test.json"
    if json.loads(source_test_path.read_text()) != dict(source_test):
        raise ValueError("source-test.json differs from the passing source test")
    validate_source_evidence(output, source_test, splits)
    source_hashes = {name: sha256_file(path) for name, path in SOURCE_FILES.items()}
    return {
        "experiment": EXPERIMENT,
        "status": "frozen_before_held_truth",
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "candidate_grid_frozen_before_expression": True,
        "held_context_truth_opened": False,
        "post_held_refit_possible": False,
        "source_code_sha256": source_hashes,
        "splits": {
            name: {
                "count": len(values),
                "targets": list(values),
                "targets_sha256": sha256_lines(values),
            }
            for name, values in splits.items()
        },
        "source_test": dict(source_test),
        "axes": {
            "centering_targets": {
                "count": centering_axis_contract.count,
                "sha256": centering_axis_contract.sha256,
            },
            "gene_axis": {
                "count": response_axis.gene_axis_contract.count,
                "sha256": response_axis.gene_axis_contract.sha256,
            },
            "panel_targets": {
                "count": response_axis.panel_targets_contract.count,
                "sha256": response_axis.panel_targets_contract.sha256,
            },
            "development_targets": {
                "count": response_axis.development_targets_contract.count,
                "sha256": response_axis.development_targets_contract.sha256,
            },
            "response_axis": {
                "count": response_axis.response_axis_contract.count,
                "sha256": response_axis.response_axis_contract.sha256,
            },
            "perturbmatch_targets": {
                "count": canonical_pm_targets_contract.count,
                "sha256": canonical_pm_targets_contract.sha256,
            },
            "perturbmatch_genes": {
                "count": canonical_pm_genes_contract.count,
                "sha256": canonical_pm_genes_contract.sha256,
            },
        },
        "frozen_artifacts": {
            filename: sha256_file(output / filename) for filename in files
        },
        "promotion_gates": protocol["promotion_gates"],
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check-contract",
        action="store_true",
        help="validate only source code against the frozen preregistration",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    if not args.check_contract:
        raise SystemExit(
            "source-only scaffold: no real-data fit command exists; "
            "use --check-contract"
        )
    protocol = validate_preregistration()
    print(
        json.dumps(
            {
                "experiment": EXPERIMENT,
                "status": protocol["status"],
                "candidate_count": len(model.candidate_grid()),
                "protocol_sha256": sha256_file(PROTOCOL_PATH),
                "real_data_opened": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
