#!/usr/bin/env python3
"""Build Experiment 151's source-only held-HCT selector candidate.

The selector reads only the frozen Experiment 149/150 effect archives. During
candidate emission, the held HCT AnnData is opened in backed mode and expression
is read only from the preflighted contiguous non-targeting rows.
"""

from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import os
import resource
import shutil
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PREREGISTRATION_PATH = HERE / "PREREGISTRATION.v1.json"
PREREGISTRATION_DIGEST_PATH = HERE / "PREREGISTRATION.v1.sha256"
BUILD_AUTHORIZATION_PATH = HERE / "SOURCE_BUILD_AUTHORIZATION.v1.json"
OUTPUT_DIR = HERE / "work/source-build-v1"
CONTROL = "non-targeting"
PERT_COL = "target_gene"
EXPECTED_CONTROL_START = 19_800
EXPECTED_CONTROL_STOP = 38_200
CONTROL_CHUNK_ROWS = 256
GIB = 1024**3
RUNTIME_SOURCE_PATHS = (
    "experiments/149-dual-moment-hek293t/build_predictions.py",
    "experiments/147-dual-moment-emitter/build_candidates.py",
    "experiments/147-dual-moment-emitter/dual_moment.py",
    "scripts/make_submission_pm_cd4.py",
)


@dataclass(frozen=True)
class SelectorPlan:
    ranked_targets: tuple[str, ...]
    top_targets: tuple[str, ...]
    bottom_targets: tuple[str, ...]
    primary_targets_hct_order: tuple[str, ...]
    hct_targets_order: tuple[str, ...]
    hct_only_targets: tuple[str, ...]
    scores: dict[str, float]
    eligible_counts: dict[str, int]
    lambda_by_target: dict[str, float]
    hashes: dict[str, str]

    def to_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["lambda_counts"] = {
            str(value): sum(
                assigned == value for assigned in self.lambda_by_target.values()
            )
            for value in sorted(set(self.lambda_by_target.values()))
        }
        return payload


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_strings(values: tuple[str, ...] | list[str] | np.ndarray) -> str:
    payload = "\n".join(map(str, values)) + "\n"
    return hashlib.sha256(payload.encode()).hexdigest()


def stable_seed(base: int, target: str, blend_lambda: float) -> int:
    payload = f"{base}:{target}:{blend_lambda:.2f}".encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "little")


def peak_rss_gib() -> float:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value * 1024 / GIB


def load_preregistration() -> dict[str, Any]:
    digest_line = PREREGISTRATION_DIGEST_PATH.read_text().strip().split()
    if len(digest_line) != 2 or digest_line[1] != PREREGISTRATION_PATH.name:
        raise ValueError("malformed immutable preregistration digest file")
    observed = sha256(PREREGISTRATION_PATH)
    if observed != digest_line[0]:
        raise ValueError(
            f"immutable preregistration hash changed: {observed} != {digest_line[0]}"
        )
    protocol = json.loads(PREREGISTRATION_PATH.read_text())
    if protocol.get("experiment") != "151-relaxation-consensus-selector":
        raise ValueError("wrong preregistration experiment")
    if protocol.get("status") != "preregistered_before_hct_outcome_access":
        raise ValueError("Experiment 151 is not in its frozen preregistered state")
    return protocol


def verify_build_authorization() -> dict[str, Any]:
    if not BUILD_AUTHORIZATION_PATH.is_file():
        raise PermissionError("Experiment 151 source build is not authorized")
    authorization = json.loads(BUILD_AUTHORIZATION_PATH.read_text())
    if authorization.get("experiment") != "151-relaxation-consensus-selector":
        raise ValueError("wrong source-build authorization experiment")
    if authorization.get("status") != "authorized_once_source_only":
        raise PermissionError("Experiment 151 source build is not authorized once")
    builder = authorization.get("builder")
    if not isinstance(builder, dict):
        raise ValueError("source-build authorization has no builder binding")
    if builder.get("sha256") != sha256(Path(__file__).resolve()):
        raise ValueError("builder source does not match source-build authorization")
    bindings = authorization.get("bindings_sha256")
    if not isinstance(bindings, dict) or not bindings:
        raise ValueError("source-build authorization has no bindings")
    observed = {
        str((ROOT / relative).relative_to(ROOT)): sha256(ROOT / relative)
        for relative in bindings
    }
    if observed != bindings:
        raise ValueError("source-build authorization bindings do not match")
    return authorization


def verify_source_artifacts(protocol: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Verify source-only files without scanning held perturbation expression."""
    report: dict[str, dict[str, Any]] = {}
    for name, record in protocol["source_artifacts"].items():
        path = ROOT / record["path"]
        if not path.is_file():
            raise FileNotFoundError(path)
        if name == "held_hct_harness":
            report[name] = {
                "path": record["path"],
                "expected_sha256": record["sha256"],
                "verification": (
                    "existence_only; full hashing intentionally deferred because it "
                    "would scan sealed perturbation-expression bytes"
                ),
            }
            continue
        observed = sha256(path)
        if observed != record["sha256"]:
            raise ValueError(
                f"{name} hash changed: {observed} != {record['sha256']}"
            )
        report[name] = {
            "path": record["path"],
            "sha256": observed,
            "verification": "sha256",
        }
    return report


def _load_effect_archive(
    path: Path,
    required_keys: tuple[str, ...],
) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        missing = sorted(set(required_keys) - set(archive.files))
        if missing:
            raise ValueError(f"{path}: missing required arrays {missing}")
        return {key: np.asarray(archive[key]) for key in archive.files}


def compute_selector(
    protocol: dict[str, Any],
    *,
    verify_sources: bool = True,
) -> tuple[SelectorPlan, dict[str, np.ndarray], dict[str, np.ndarray]]:
    if verify_sources:
        verify_source_artifacts(protocol)
    required = tuple(protocol["required_effect_keys"])
    hek_path = ROOT / protocol["source_artifacts"]["hek_effects"]["path"]
    hct_path = ROOT / protocol["source_artifacts"]["hct_effects"]["path"]
    hek = _load_effect_archive(hek_path, required)
    hct = _load_effect_archive(hct_path, required)

    hek_genes = np.asarray(hek["genes"]).astype(str)
    hct_genes = np.asarray(hct["genes"]).astype(str)
    if not np.array_equal(hek_genes, hct_genes):
        raise ValueError("HEK and HCT frozen effect gene axes differ")

    hek_targets = np.asarray(hek["targets"]).astype(str)
    hct_targets = np.asarray(hct["targets"]).astype(str)
    hek_index = {target: index for index, target in enumerate(hek_targets)}
    hct_index = {target: index for index, target in enumerate(hct_targets)}
    common = set(hek_index) & set(hct_index)
    selector = protocol["selector"]
    if len(common) != int(selector["common_target_count"]):
        raise ValueError(
            f"common target count changed: {len(common)} != "
            f"{selector['common_target_count']}"
        )

    hct_only = tuple(sorted(set(hct_index) - set(hek_index)))
    expected_hct_only = tuple(
        sorted(selector["hct_only_targets_excluded_from_primary_analysis"])
    )
    if hct_only != expected_hct_only:
        raise ValueError(f"HCT-only targets changed: {hct_only} != {expected_hct_only}")

    hek_relaxation = (
        np.asarray(hek["direction"], dtype=np.float64)
        - np.asarray(hek["incumbent"], dtype=np.float64)
    )
    hct_relaxation = (
        np.asarray(hct["direction"], dtype=np.float64)
        - np.asarray(hct["incumbent"], dtype=np.float64)
    )
    scored: list[tuple[str, float]] = []
    eligible_counts: dict[str, int] = {}
    expected_eligible = int(selector["eligible_gene_count_per_target"])
    for target in common:
        hek_values = hek_relaxation[hek_index[target]]
        hct_values = hct_relaxation[hct_index[target]]
        weights = np.minimum(np.abs(hek_values), np.abs(hct_values))
        eligible = (
            np.isfinite(hek_values)
            & np.isfinite(hct_values)
            & np.isfinite(weights)
            & (weights > 0)
        )
        count = int(eligible.sum())
        if count != expected_eligible:
            raise ValueError(
                f"{target}: eligible gene count {count} != {expected_eligible}"
            )
        denominator = float(weights[eligible].sum())
        if denominator <= 0:
            raise ValueError(f"{target}: selector has zero eligible weight")
        agreement = hek_values[eligible] * hct_values[eligible] > 0
        score = float(
            np.sum(weights[eligible] * agreement, dtype=np.float64) / denominator
        )
        scored.append((target, score))
        eligible_counts[target] = count

    scored.sort(key=lambda item: (-item[1], item[0]))
    ranked = tuple(target for target, _ in scored)
    top_count = int(selector["top_count"])
    bottom_count = int(selector["bottom_count"])
    top = ranked[:top_count]
    bottom = ranked[top_count:]
    if len(bottom) != bottom_count:
        raise ValueError(f"bottom target count {len(bottom)} != {bottom_count}")
    primary_hct_order = tuple(target for target in hct_targets if target in common)
    hashes = {
        "top_names_sha256_newline_terminated": sha256_strings(top),
        "bottom_names_sha256_newline_terminated": sha256_strings(bottom),
        "ranked_names_sha256_newline_terminated": sha256_strings(ranked),
        "primary_targets_in_hct_order_sha256_newline_terminated": sha256_strings(
            primary_hct_order
        ),
    }
    for name, observed in hashes.items():
        expected = selector[name]
        if observed != expected:
            raise ValueError(f"selector hash {name} changed: {observed} != {expected}")

    scores = {target: score for target, score in scored}
    if top[-1] != selector["last_top_target"]:
        raise ValueError("last top target changed")
    if bottom[0] != selector["first_bottom_target"]:
        raise ValueError("first bottom target changed")
    if bool(scores[top[-1]] == scores[bottom[0]]) != bool(selector["boundary_tie"]):
        raise ValueError("selector boundary tie state changed")
    for target, expected_score in (
        (top[-1], float(selector["last_top_score"])),
        (bottom[0], float(selector["first_bottom_score"])),
    ):
        if not np.isclose(scores[target], expected_score, rtol=0.0, atol=1e-9):
            raise ValueError(
                f"{target}: boundary score changed: {scores[target]} != {expected_score}"
            )

    top_lambda = float(selector["top_lambda"])
    bottom_lambda = float(selector["bottom_lambda"])
    hct_only_lambda = float(selector["hct_only_lambda"])
    lambda_by_target = {target: top_lambda for target in top}
    lambda_by_target.update({target: bottom_lambda for target in bottom})
    lambda_by_target.update({target: hct_only_lambda for target in hct_only})
    if set(lambda_by_target) != set(hct_targets):
        raise ValueError("lambda assignment does not cover the full HCT target set")

    return (
        SelectorPlan(
            ranked_targets=ranked,
            top_targets=top,
            bottom_targets=bottom,
            primary_targets_hct_order=primary_hct_order,
            hct_targets_order=tuple(map(str, hct_targets)),
            hct_only_targets=hct_only,
            scores=scores,
            eligible_counts=eligible_counts,
            lambda_by_target=lambda_by_target,
            hashes=hashes,
        ),
        hek,
        hct,
    )


def _load_build_runtime() -> ModuleType:
    return load_module(
        "exp151_build_base",
        ROOT / "experiments/149-dual-moment-hek293t/build_predictions.py",
    )


def runtime_source_hashes() -> dict[str, str]:
    return {relative: sha256(ROOT / relative) for relative in RUNTIME_SOURCE_PATHS}


def build_candidate(
    protocol: dict[str, Any],
    plan: SelectorPlan,
    hct: dict[str, np.ndarray],
    staging: Path,
) -> dict[str, Any]:
    import anndata as ad
    import h5py

    base = _load_build_runtime()
    dm = base.DM
    build = base.BUILD
    real_path = ROOT / protocol["source_artifacts"]["held_hct_harness"]["path"]
    output_path = staging / "candidate.h5ad"
    metadata_path = staging / f".metadata-{os.getpid()}.h5ad"

    genes = np.asarray(hct["genes"]).astype(str)
    targets = np.asarray(hct["targets"]).astype(str)
    incumbent_effects = np.asarray(hct["incumbent"], dtype=np.float64)
    direction_effects = np.asarray(hct["direction"], dtype=np.float64)
    bctrl = np.asarray(hct["bctrl"], dtype=np.float64)
    fixed = protocol["fixed_emission"]
    if incumbent_effects.shape != direction_effects.shape:
        raise ValueError("HCT incumbent and direction effect shapes differ")
    if incumbent_effects.shape != (len(targets), len(genes)):
        raise ValueError("HCT effect arrays do not match their frozen axes")
    if int(fixed["cells_per_target"]) != base.CELLS_PER_TARGET:
        raise ValueError("Experiment 151 cell count differs from Experiment 150")

    real = ad.read_h5ad(real_path, backed="r")
    temporary_output = output_path.with_name(
        f".{output_path.name}.tmp-{os.getpid()}"
    )
    try:
        labels = real.obs[PERT_COL].astype(str).to_numpy()
        control_indices = np.flatnonzero(labels == CONTROL)
        expected_indices = np.arange(
            EXPECTED_CONTROL_START,
            EXPECTED_CONTROL_STOP,
            dtype=np.int64,
        )
        if not np.array_equal(control_indices, expected_indices):
            raise ValueError("HCT control-row preflight contract changed")
        if not np.array_equal(np.asarray(real.var_names).astype(str), genes):
            raise ValueError("HCT harness and frozen effect gene axes differ")
        real_targets = set(labels) - {CONTROL}
        if real_targets != set(targets):
            raise ValueError("HCT metadata target set differs from frozen effects")

        libraries = base.control_libraries(
            real,
            EXPECTED_CONTROL_START,
            EXPECTED_CONTROL_STOP,
            chunk_rows=CONTROL_CHUNK_ROWS,
        )
        control_sums = np.zeros(len(genes), dtype=np.float64)
        for start in range(
            EXPECTED_CONTROL_START,
            EXPECTED_CONTROL_STOP,
            CONTROL_CHUNK_ROWS,
        ):
            stop = min(start + CONTROL_CHUNK_ROWS, EXPECTED_CONTROL_STOP)
            control_sums += np.asarray(
                build.read_rows(real, slice(start, stop)).sum(axis=0),
                dtype=np.float64,
            ).ravel()
        observed_bctrl = np.log1p(50_000.0 * control_sums / control_sums.sum())
        if not np.array_equal(
            observed_bctrl.astype(np.float32),
            bctrl.astype(np.float32),
        ):
            raise ValueError("HCT controls differ from the frozen control profile")

        build.write_metadata_template(
            metadata_path,
            target_order=list(map(str, targets)),
            cells_per_target=base.CELLS_PER_TARGET,
            control_cells=len(control_indices),
            var=real.var,
        )
        inc_pick_rng, inc_emission_rng = base.PROD.make_rngs(0)
        dir_pick_rng, dir_emission_rng = base.PROD.make_rngs(0)
        per_target: dict[str, dict[str, Any]] = {}
        with (
            h5py.File(metadata_path, "r") as metadata_h5,
            h5py.File(temporary_output, "w") as candidate_h5,
        ):
            build.clone_non_x_metadata(metadata_h5, candidate_h5)
            n_rows = (
                len(targets) * base.CELLS_PER_TARGET + len(control_indices)
            )
            writer = build.CSRStreamWriter(
                candidate_h5,
                n_rows=n_rows,
                n_cols=len(genes),
                data_dtype=np.dtype("float32"),
            )
            for index, target_value in enumerate(targets):
                target = str(target_value)
                incumbent = base.emitted_block(
                    incumbent_effects[index],
                    bctrl,
                    libraries,
                    inc_pick_rng,
                    inc_emission_rng,
                )
                direction = base.emitted_block(
                    direction_effects[index],
                    bctrl,
                    libraries,
                    dir_pick_rng,
                    dir_emission_rng,
                )
                blend_lambda = float(plan.lambda_by_target[target])
                candidate, diagnostics = dm.emit_dual_moment(
                    incumbent,
                    direction,
                    blend_lambda=blend_lambda,
                    depth_contrast=float(fixed["depth_contrast"]),
                    quasi_spread=float(fixed["within_group_quasi_spread"]),
                    seed=stable_seed(
                        int(fixed["dual_moment_seed"]),
                        target,
                        blend_lambda,
                    ),
                )
                incumbent_columns = np.asarray(
                    incumbent.sum(axis=0), dtype=np.int64
                ).ravel()
                candidate_columns = np.asarray(
                    candidate.sum(axis=0), dtype=np.int64
                ).ravel()
                if not np.array_equal(incumbent_columns, candidate_columns):
                    raise RuntimeError(f"{target}: candidate column totals changed")
                if not np.array_equal(
                    dm.scorer_pseudobulk(incumbent),
                    dm.scorer_pseudobulk(candidate),
                ):
                    raise RuntimeError(f"{target}: candidate pseudobulk changed")
                row_totals = np.asarray(
                    candidate.sum(axis=1), dtype=np.int64
                ).ravel()
                if (
                    np.any(row_totals <= 0)
                    or row_totals.max(initial=0) > dm.MAX_COUNTS_PER_CELL
                ):
                    raise RuntimeError(f"{target}: invalid candidate row total")
                writer.append(candidate.astype(np.float32))
                per_target[target] = {
                    "lambda": blend_lambda,
                    "selector_group": (
                        "top"
                        if target in plan.top_targets
                        else "bottom"
                        if target in plan.bottom_targets
                        else "hct_only"
                    ),
                    "selector_score": plan.scores.get(target),
                    "columns_exact": True,
                    "pseudobulk_exact": True,
                    "incumbent_columns_sha256": hashlib.sha256(
                        np.ascontiguousarray(incumbent_columns, dtype="<i8").tobytes()
                    ).hexdigest(),
                    "candidate_nnz": int(candidate.nnz),
                    "incumbent_nnz": int(incumbent.nnz),
                    "max_row_total": int(row_totals.max(initial=0)),
                    "dual_moment": diagnostics.to_dict(),
                }
                if (index + 1) % 10 == 0:
                    print(f"[emit] {index + 1}/{len(targets)} targets", flush=True)
                del incumbent, direction, candidate
                gc.collect()

            for start in range(
                EXPECTED_CONTROL_START,
                EXPECTED_CONTROL_STOP,
                CONTROL_CHUNK_ROWS,
            ):
                stop = min(start + CONTROL_CHUNK_ROWS, EXPECTED_CONTROL_STOP)
                control_block = build.read_rows(real, slice(start, stop))
                writer.append(control_block.astype(np.float32))
                del control_block
            writer.finish()
            candidate_h5.flush()
            output_nnz = writer.nnz

        os.replace(temporary_output, output_path)
        metadata_path.unlink(missing_ok=True)
        built = ad.read_h5ad(output_path, backed="r")
        try:
            built_labels = built.obs[PERT_COL].astype(str).to_numpy()
            if built.n_vars != len(genes):
                raise RuntimeError("completed candidate gene count differs")
            if not np.array_equal(np.asarray(built.var_names).astype(str), genes):
                raise RuntimeError("completed candidate gene axis differs")
            if np.count_nonzero(built_labels == CONTROL) != len(control_indices):
                raise RuntimeError("completed candidate control count differs")
            for target in targets:
                if (
                    np.count_nonzero(built_labels == str(target))
                    != base.CELLS_PER_TARGET
                ):
                    raise RuntimeError(f"{target}: completed candidate row count differs")
        finally:
            built.file.close()

        return {
            "candidate": str((OUTPUT_DIR / "candidate.h5ad").relative_to(ROOT)),
            "candidate_sha256": sha256(output_path),
            "candidate_nnz": int(output_nnz),
            "target_order_sha256_newline_terminated": sha256_strings(targets),
            "targets": len(targets),
            "controls": len(control_indices),
            "real_control_rows_read": len(control_indices),
            "real_perturbation_rows_read": 0,
            "all_column_totals_exact": all(
                record["columns_exact"] for record in per_target.values()
            ),
            "all_pseudobulks_exact": all(
                record["pseudobulk_exact"] for record in per_target.values()
            ),
            "per_target": per_target,
        }
    finally:
        real.file.close()
        metadata_path.unlink(missing_ok=True)
        temporary_output.unlink(missing_ok=True)


def main() -> int:
    started = time.monotonic()
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"immutable build output already exists: {OUTPUT_DIR}")
    protocol = load_preregistration()
    authorization = verify_build_authorization()
    source_verification = verify_source_artifacts(protocol)
    plan, _, hct = compute_selector(protocol, verify_sources=False)
    staging = OUTPUT_DIR.with_name(f".{OUTPUT_DIR.name}.tmp-{os.getpid()}")
    if staging.exists():
        raise FileExistsError(f"staging directory already exists: {staging}")
    staging.mkdir(parents=True)
    try:
        (staging / "selector.json").write_text(
            json.dumps(plan.to_dict(), indent=2, sort_keys=True) + "\n"
        )
        prediction = build_candidate(protocol, plan, hct, staging)
        report = {
            "experiment": "151-relaxation-consensus-selector",
            "stage": "source_only_build",
            "preregistration_sha256": sha256(PREREGISTRATION_PATH),
            "preregistration_digest_file_sha256": sha256(
                PREREGISTRATION_DIGEST_PATH
            ),
            "source_build_authorization_sha256": sha256(
                BUILD_AUTHORIZATION_PATH
            ),
            "source_build_authorization_status": authorization["status"],
            "builder_sha256": sha256(Path(__file__).resolve()),
            "runtime_source_hashes": runtime_source_hashes(),
            "source_verification": source_verification,
            "selector": plan.to_dict(),
            "prediction": prediction,
            "peak_rss_gib": peak_rss_gib(),
            "elapsed_seconds": time.monotonic() - started,
            "held_hct_expression_access": {
                "control_rows_read": prediction["real_control_rows_read"],
                "perturbation_rows_read": 0,
            },
            "cell_eval2_invoked": False,
            "challenge_build_created": False,
            "challenge_submission_created": False,
        }
        (staging / "build-report.json").write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        os.replace(staging, OUTPUT_DIR)
        print(
            json.dumps(
                {
                    "candidate_sha256": prediction["candidate_sha256"],
                    "lambda_counts": plan.to_dict()["lambda_counts"],
                    "all_column_totals_exact": prediction[
                        "all_column_totals_exact"
                    ],
                    "all_pseudobulks_exact": prediction[
                        "all_pseudobulks_exact"
                    ],
                    "real_perturbation_rows_read": 0,
                },
                indent=2,
            )
        )
        return 0
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
