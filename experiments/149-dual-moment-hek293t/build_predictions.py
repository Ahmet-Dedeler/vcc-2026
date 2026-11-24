#!/usr/bin/env python3
"""Build held-HEK incumbent and dual-moment predictions without opening truth."""

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
from pathlib import Path
from types import ModuleType
from typing import Any

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
OUTPUT_DIR = HERE / "work/source-build-v1"
REAL_PATH = ROOT / "data/harness/orion-hek293t-thinned.h5ad"
REFERENCE_104 = ROOT / (
    "experiments/104-context-heldout-credit-adapter/"
    "work-v2/hek293t/frozen.npz"
)
REFERENCE_106 = ROOT / (
    "experiments/106-state-baseline/work/runs/20260830T064800Z/"
    "incumbent-a14/hek293t/frozen.npz"
)
UNPRUNED_OPERATOR = ROOT / (
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/"
    "held_orion_hek293t.npz"
)
PRUNED_OPERATOR = ROOT / (
    "experiments/106-state-baseline/operators/"
    "held_orion_hek293t_pruned.npz"
)
HARNESS_METADATA = ROOT / "data/harness/orion-hek293t-thinned.h5ad.meta.json"
CONTROL = "non-targeting"
PERT_COL = "target_gene"
CELLS_PER_TARGET = 400
MAX_RSS_GIB = 24.0
GIB = 1024**3


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


E104 = load_module(
    "exp149_e104",
    ROOT / "experiments/104-context-heldout-credit-adapter/run.py",
)
E106 = load_module(
    "exp149_e106",
    ROOT / "experiments/106-state-baseline/freeze_a14_endpoints.py",
)
PROD = load_module(
    "exp149_production",
    ROOT / "scripts/make_submission_pm_cd4.py",
)
DM = load_module(
    "exp149_dual_moment",
    ROOT / "experiments/147-dual-moment-emitter/dual_moment.py",
)
BUILD = load_module(
    "exp149_stream_writer",
    ROOT / "experiments/147-dual-moment-emitter/build_candidates.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def array_sha256(value: np.ndarray) -> str:
    array = np.ascontiguousarray(value, dtype="<f8")
    return hashlib.sha256(array.tobytes()).hexdigest()


def peak_rss_gib() -> float:
    value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value * 1024 / GIB


def memory_checkpoint(label: str) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {label}: {peak:.3f} GiB", flush=True)
    if peak >= MAX_RSS_GIB:
        raise MemoryError(f"{label}: peak RSS {peak:.3f} GiB exceeded limit")


def verify_protocol() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL_PATH.read_text())
    if protocol.get("status") != "preregistered_source_only":
        raise ValueError("experiment is not source-only preregistered")
    expected = protocol["source_bindings_sha256"]
    observed = {
        str((ROOT / relative).relative_to(ROOT)): sha256(ROOT / relative)
        for relative in expected
    }
    if observed != expected:
        raise ValueError("source bindings differ from preregistration")
    builder = protocol["builder"]
    if builder["sha256"] != sha256(Path(__file__).resolve()):
        raise ValueError("builder source differs from preregistration")
    return protocol


def verify_biological_inputs() -> dict[str, str]:
    if sha256(E106.A14_REFERENCE_META) != E106.A14_REFERENCE_META_SHA256:
        raise ValueError("a15 reference manifest changed")
    reference = json.loads(E106.A14_REFERENCE_META.read_text())
    expected = reference["hashes"]["inputs"]
    observed: dict[str, str] = {}
    for key, path in E106.a14_biological_inputs("hek293t").items():
        actual = sha256(path)
        if expected.get(key) != actual:
            raise ValueError(f"a15 biological input changed: {key}")
        observed[str(path.relative_to(ROOT))] = actual
    return observed


def load_reference(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        return {key: np.asarray(archive[key]) for key in archive.files}


def load_operator(
    path: Path,
    banks: list[object],
    official: tuple[str, ...],
) -> tuple[object, np.ndarray]:
    return E104.E064.load_frozen_operator(
        path,
        banks,
        official,
        E104.PANEL_SET,
        expected_sha256=sha256(path),
    )


def apply_operator_gamma(
    matrix: np.ndarray,
    operator: object,
    columns: np.ndarray,
    gamma: float,
) -> tuple[np.ndarray, dict[str, float]]:
    output = np.asarray(matrix, dtype=np.float64).copy()
    before = float(np.median(np.linalg.norm(output[:, columns], axis=1))) or 1.0
    weights = np.asarray(operator.evals, dtype=np.float64) ** (-gamma / 2.0)
    transformed = (
        (output[:, columns] @ operator.evecs) * weights
    ) @ operator.evecs.T
    after = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    scale = before / after
    output[:, columns] = transformed * scale
    return output, {
        "gamma": gamma,
        "before_operator_median_norm": before,
        "after_operator_median_norm": after,
        "applied_scale": scale,
    }


def prepare_source_components(
    targets: np.ndarray,
) -> tuple[np.ndarray, list[object], object, dict[str, Any]]:
    official = E104.OFFICIAL
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[gene] for gene in E104.PANEL if gene in official_pos],
        dtype=np.int64,
    )
    source_names = E104.SURFACES["hek293t"]["sources"]
    banks = {
        name: E104.E062.Bank.load(name, E104.REAL_DONOR_PATHS[name])
        for name in source_names
    }
    prepared = {
        name: E104.E062.prepare_main_source(
            bank,
            targets,
            official,
            panel_positions,
            chunk_rows=256,
            max_rss_gib=MAX_RSS_GIB,
        )
        for name, bank in banks.items()
    }
    pm_bank = E104.E064.load_perturbmatch(E104.PM_PATH)
    pm_prepared = E104.E062.prepare_main_source(
        pm_bank,
        targets,
        official,
        panel_positions,
        chunk_rows=256,
        max_rss_gib=MAX_RSS_GIB,
    )
    shared = E104.E064.shared_k562_pm_columns(
        banks["k562_gwps"],
        pm_bank,
        official,
        E104.PANEL_SET,
    )
    pm_slot, pm_supported = E104._partial_pm_slot(
        prepared["k562_gwps"],
        pm_prepared,
        shared,
    )
    source_prepared = [
        pm_slot if name == "k562_gwps" else prepared[name]
        for name in source_names
    ]
    pre, votes = E104.E062.build_baseline(source_prepared)
    pre *= E104.ALPHA
    cd4_banks = [
        E104.E062.Bank.load(
            name,
            E104.CD4_ROOT / name,
            require_cd4_scale=True,
        )
        for name in E104.CD4_NAMES
    ]
    cd4 = E104.E062.prepare_cd4_mean(
        cd4_banks,
        targets,
        official,
        E104.PANEL_SET,
        "gene_norm",
        chunk_rows=256,
        max_rss_gib=MAX_RSS_GIB,
    )
    audit = {
        "sources": list(source_names),
        "pm_supported_targets": int(pm_supported.sum()),
        "cd4_supported_targets": int(np.sum(cd4.target_condition_count > 0)),
        "vote_histogram": {
            str(value): int(np.sum(votes == value))
            for value in np.unique(votes)
        },
    }
    return pre, [banks[name] for name in source_names], cd4, audit


def build_effects(
    *,
    pre: np.ndarray,
    banks: list[object],
    cd4: object,
    genes: np.ndarray,
    operator_path: Path,
    gamma: float,
    cd4_fraction: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    operator, columns = load_operator(operator_path, banks, E104.OFFICIAL)
    whitened, whitening = apply_operator_gamma(pre, operator, columns, gamma)
    final, cd4_audit = E104.E062.add_post_whiten_residual(
        whitened,
        cd4,
        cd4_fraction,
        columns,
    )
    official = E104.OFFICIAL
    official_pos = {gene: index for index, gene in enumerate(official)}
    panel_positions = np.asarray(
        [official_pos[gene] for gene in E104.PANEL if gene in official_pos],
        dtype=np.int64,
    )
    final[:, panel_positions] = 0.0
    held_columns = np.asarray(
        [official_pos.get(str(gene), -1) for gene in genes],
        dtype=np.int64,
    )
    if np.any(held_columns < 0):
        raise ValueError("held HEK gene axis is not contained in official genes")
    held = final[:, held_columns]
    del operator
    gc.collect()
    return held, {
        "operator": str(operator_path.relative_to(ROOT)),
        "operator_sha256": sha256(operator_path),
        "gamma": gamma,
        "cd4_fraction": cd4_fraction,
        "whitening": whitening,
        "cd4": cd4_audit,
        "effect_sha256": array_sha256(held),
    }


def control_libraries(
    handle: ad.AnnData,
    start_row: int,
    stop_row: int,
    chunk_rows: int = 256,
) -> np.ndarray:
    libraries = np.zeros(stop_row - start_row, dtype=np.int64)
    for start in range(start_row, stop_row, chunk_rows):
        stop = min(start + chunk_rows, stop_row)
        block = BUILD.read_rows(handle, slice(start, stop))
        libraries[start - start_row : stop - start_row] = np.asarray(
            block.sum(axis=1),
            dtype=np.int64,
        ).ravel()
    if np.any(libraries <= 0):
        raise ValueError("isolated controls contain an empty cell")
    return libraries


def emitted_block(
    effect: np.ndarray,
    bctrl: np.ndarray,
    libraries: np.ndarray,
    pick_rng: np.random.Generator,
    emission_rng: np.random.Generator,
) -> sp.csr_matrix:
    profile = PROD.profile_from_effect(
        effect,
        bctrl,
        min_cpm=10.0,
    )
    picks = pick_rng.choice(
        len(libraries),
        CELLS_PER_TARGET,
        replace=len(libraries) < CELLS_PER_TARGET,
    )
    spread = PROD.emission_spread_vector(
        effect,
        bctrl,
        min_cpm=10.0,
        tail_spread=1.0,
        keep_top=8_000,
        keep_spread=0.25,
        spread_map="block",
        spread_lo=0.0,
        spread_hi=0.60,
    )
    values, columns, pointers = PROD.emit_counts(
        profile,
        libraries[picks],
        0.1,
        emission_rng,
        mode="quasi",
        spread=spread,
        overdispersion=4.0,
    )
    return DM.as_integer_csr(
        sp.csr_matrix(
            (values, columns, pointers),
            shape=(CELLS_PER_TARGET, len(effect)),
        )
    )


def write_predictions(
    output_dir: Path,
    genes: np.ndarray,
    targets: np.ndarray,
    bctrl: np.ndarray,
    incumbent_effects: np.ndarray,
    direction_effects: np.ndarray,
) -> dict[str, Any]:
    incumbent_path = output_dir / "incumbent.h5ad"
    candidate_path = output_dir / "candidate.h5ad"
    if incumbent_path.exists() or candidate_path.exists():
        raise FileExistsError("immutable prediction output already exists")
    controls = ad.read_h5ad(REAL_PATH, backed="r")
    temporary_paths: list[Path] = []
    metadata_path = output_dir / f".metadata-{os.getpid()}.h5ad"
    try:
        labels = controls.obs[PERT_COL].astype(str).to_numpy()
        control_indices = np.flatnonzero(labels == CONTROL)
        if (
            control_indices.size != 18_400
            or not np.array_equal(
                control_indices,
                np.arange(19_800, 38_200, dtype=np.int64),
            )
        ):
            raise ValueError("HEK control-row preflight contract changed")
        control_start = int(control_indices[0])
        control_stop = int(control_indices[-1] + 1)
        if not np.array_equal(np.asarray(controls.var_names).astype(str), genes):
            raise ValueError("HEK harness and frozen gene axes differ")
        libraries = control_libraries(
            controls,
            control_start,
            control_stop,
        )
        sums = np.zeros(len(genes), dtype=np.float64)
        for start in range(control_start, control_stop, 256):
            stop = min(start + 256, control_stop)
            sums += np.asarray(
                BUILD.read_rows(controls, slice(start, stop)).sum(axis=0),
                dtype=np.float64,
            ).ravel()
        observed_bctrl = np.log1p(50_000.0 * sums / sums.sum())
        if not np.array_equal(
            observed_bctrl.astype(np.float32),
            bctrl.astype(np.float32),
        ):
            raise ValueError("harness controls differ from frozen control profile")

        BUILD.write_metadata_template(
            metadata_path,
            target_order=list(map(str, targets)),
            cells_per_target=CELLS_PER_TARGET,
            control_cells=control_indices.size,
            var=controls.var,
        )
        temporary_paths.append(metadata_path)
        incumbent_tmp = incumbent_path.with_name(
            f".{incumbent_path.name}.tmp-{os.getpid()}"
        )
        candidate_tmp = candidate_path.with_name(
            f".{candidate_path.name}.tmp-{os.getpid()}"
        )
        temporary_paths.extend((incumbent_tmp, candidate_tmp))
        for path in temporary_paths:
            if path != metadata_path:
                path.unlink(missing_ok=True)

        inc_pick_rng, inc_emission_rng = PROD.make_rngs(0)
        dir_pick_rng, dir_emission_rng = PROD.make_rngs(0)
        per_target: dict[str, dict[str, Any]] = {}
        with (
            h5py.File(metadata_path, "r") as metadata_h5,
            h5py.File(incumbent_tmp, "w") as incumbent_h5,
            h5py.File(candidate_tmp, "w") as candidate_h5,
        ):
            BUILD.clone_non_x_metadata(metadata_h5, incumbent_h5)
            BUILD.clone_non_x_metadata(metadata_h5, candidate_h5)
            n_rows = len(targets) * CELLS_PER_TARGET + control_indices.size
            incumbent_writer = BUILD.CSRStreamWriter(
                incumbent_h5,
                n_rows=n_rows,
                n_cols=len(genes),
                data_dtype=np.dtype("float32"),
            )
            candidate_writer = BUILD.CSRStreamWriter(
                candidate_h5,
                n_rows=n_rows,
                n_cols=len(genes),
                data_dtype=np.dtype("float32"),
            )
            for index, target_value in enumerate(targets):
                target = str(target_value)
                incumbent = emitted_block(
                    incumbent_effects[index],
                    bctrl,
                    libraries,
                    inc_pick_rng,
                    inc_emission_rng,
                )
                direction = emitted_block(
                    direction_effects[index],
                    bctrl,
                    libraries,
                    dir_pick_rng,
                    dir_emission_rng,
                )
                candidate, diagnostics = DM.emit_dual_moment(
                    incumbent,
                    direction,
                    blend_lambda=1.0,
                    depth_contrast=0.5,
                    quasi_spread=0.25,
                    seed=int.from_bytes(
                        hashlib.sha256(
                            f"14700:{target}:1.00".encode()
                        ).digest()[:4],
                        "little",
                    ),
                )
                incumbent_columns = np.asarray(
                    incumbent.sum(axis=0),
                    dtype=np.int64,
                ).ravel()
                candidate_columns = np.asarray(
                    candidate.sum(axis=0),
                    dtype=np.int64,
                ).ravel()
                if not np.array_equal(incumbent_columns, candidate_columns):
                    raise RuntimeError(f"{target}: candidate column totals changed")
                if not np.array_equal(
                    DM.scorer_pseudobulk(incumbent),
                    DM.scorer_pseudobulk(candidate),
                ):
                    raise RuntimeError(f"{target}: candidate pseudobulk changed")
                incumbent_writer.append(incumbent.astype(np.float32))
                candidate_writer.append(candidate.astype(np.float32))
                per_target[target] = {
                    "rows": CELLS_PER_TARGET,
                    "pseudobulk_exact": True,
                    "incumbent_nnz": int(incumbent.nnz),
                    "candidate_nnz": int(candidate.nnz),
                    "dual_moment": diagnostics.to_dict(),
                }
                if (index + 1) % 10 == 0:
                    print(
                        f"[emit] {index + 1}/{len(targets)} targets",
                        flush=True,
                    )
                del incumbent, direction, candidate
                gc.collect()

            for start in range(control_start, control_stop, 256):
                stop = min(start + 256, control_stop)
                block = BUILD.read_rows(controls, slice(start, stop))
                incumbent_writer.append(block.astype(np.float32))
                candidate_writer.append(block.astype(np.float32))
            incumbent_writer.finish()
            candidate_writer.finish()
            incumbent_h5.flush()
            candidate_h5.flush()
            incumbent_nnz = incumbent_writer.nnz
            candidate_nnz = candidate_writer.nnz

        os.replace(incumbent_tmp, incumbent_path)
        os.replace(candidate_tmp, candidate_path)
        metadata_path.unlink(missing_ok=True)
        temporary_paths.clear()
        return {
            "incumbent_sha256": sha256(incumbent_path),
            "candidate_sha256": sha256(candidate_path),
            "incumbent_nnz": int(incumbent_nnz),
            "candidate_nnz": int(candidate_nnz),
            "targets": len(per_target),
            "controls": int(control_indices.size),
            "real_control_rows_read": int(control_indices.size),
            "real_perturbation_rows_read": 0,
            "all_pseudobulk_exact": all(
                row["pseudobulk_exact"] for row in per_target.values()
            ),
            "per_target": per_target,
        }
    finally:
        controls.file.close()
        for path in temporary_paths:
            path.unlink(missing_ok=True)


def main() -> int:
    started = time.monotonic()
    if OUTPUT_DIR.exists():
        raise FileExistsError(f"immutable build directory exists: {OUTPUT_DIR}")
    verify_protocol()
    biological_inputs = verify_biological_inputs()
    staging = OUTPUT_DIR.with_name(f".{OUTPUT_DIR.name}.tmp-{os.getpid()}")
    if staging.exists():
        raise FileExistsError(f"staging directory exists: {staging}")
    staging.mkdir(parents=True)

    try:
        reference_104 = load_reference(REFERENCE_104)
        reference_106 = load_reference(REFERENCE_106)
        genes = np.asarray(reference_104["genes"]).astype(str)
        targets = np.asarray(reference_104["targets"]).astype(str)
        bctrl = np.asarray(reference_104["bctrl"], dtype=np.float64)
        if not np.array_equal(
            genes,
            np.asarray(reference_106["genes"]).astype(str),
        ):
            raise ValueError("frozen HEK reference gene axes differ")
        if not np.array_equal(
            targets,
            np.asarray(reference_106["targets"]).astype(str),
        ):
            raise ValueError("frozen HEK reference target axes differ")
        metadata = json.loads(HARNESS_METADATA.read_text())
        if set(map(str, metadata["perts"])) != set(targets):
            raise ValueError("safe harness metadata target set differs")

        pre, banks, cd4, source_audit = prepare_source_components(targets)
        unpruned_065, unpruned_audit = build_effects(
            pre=pre,
            banks=banks,
            cd4=cd4,
            genes=genes,
            operator_path=UNPRUNED_OPERATOR,
            gamma=0.65,
            cd4_fraction=E104.CD4_FRACTION,
        )
        unpruned_match = np.array_equal(
            unpruned_065.astype(np.float32),
            np.asarray(reference_104["baseline"], dtype=np.float32),
        )
        if not unpruned_match:
            raise RuntimeError(
                "failed Experiment 104 exact reconstruction bridge"
            )
        del unpruned_065
        gc.collect()

        pruned_065_cd4f005, pruned_005_audit = build_effects(
            pre=pre,
            banks=banks,
            cd4=cd4,
            genes=genes,
            operator_path=PRUNED_OPERATOR,
            gamma=0.65,
            cd4_fraction=0.05,
        )
        pruned_005_match = np.array_equal(
            pruned_065_cd4f005.astype(np.float32),
            np.asarray(reference_106["baseline"], dtype=np.float32),
        )
        if not pruned_005_match:
            raise RuntimeError(
                "failed Experiment 106 exact reconstruction bridge"
            )
        del pruned_065_cd4f005
        gc.collect()

        incumbent_effects, incumbent_audit = build_effects(
            pre=pre,
            banks=banks,
            cd4=cd4,
            genes=genes,
            operator_path=PRUNED_OPERATOR,
            gamma=0.65,
            cd4_fraction=E104.CD4_FRACTION,
        )
        direction_effects, direction_audit = build_effects(
            pre=pre,
            banks=banks,
            cd4=cd4,
            genes=genes,
            operator_path=PRUNED_OPERATOR,
            gamma=0.45,
            cd4_fraction=E104.CD4_FRACTION,
        )
        memory_checkpoint("effects ready")
        effects_path = staging / "effects.npz"
        with effects_path.open("wb") as handle:
            np.savez_compressed(
                handle,
                genes=genes.astype("U"),
                targets=targets.astype("U"),
                bctrl=bctrl.astype(np.float32),
                incumbent=incumbent_effects.astype(np.float32),
                direction=direction_effects.astype(np.float32),
            )

        predictions = write_predictions(
            staging,
            genes,
            targets,
            bctrl,
            incumbent_effects,
            direction_effects,
        )
        predictions["incumbent"] = str(
            (OUTPUT_DIR / "incumbent.h5ad").relative_to(ROOT)
        )
        predictions["candidate"] = str(
            (OUTPUT_DIR / "candidate.h5ad").relative_to(ROOT)
        )
        memory_checkpoint("predictions written")
        report = {
            "experiment": "149-dual-moment-hek293t",
            "stage": "source_only_build",
            "provenance_utc": "2026-08-31T18:45:00+00:00",
            "protocol_sha256": sha256(PROTOCOL_PATH),
            "builder_sha256": sha256(Path(__file__).resolve()),
            "real_harness_opened_for_control_rows_only": True,
            "real_harness_control_rows_read": 18_400,
            "real_harness_perturbation_rows_read": 0,
            "real_harness_path_reserved_for_scoring": (
                "data/harness/orion-hek293t-thinned.h5ad"
            ),
            "source_biological_inputs": biological_inputs,
            "source_audit": source_audit,
            "bridge_checks": {
                "experiment_104_unpruned_gamma065_cd4_solved_exact": (
                    unpruned_match
                ),
                "experiment_106_pruned_gamma065_cd4_f005_exact": (
                    pruned_005_match
                ),
            },
            "bridge_audits": {
                "experiment_104": unpruned_audit,
                "experiment_106": pruned_005_audit,
            },
            "incumbent_effect": incumbent_audit,
            "direction_effect": direction_audit,
            "effects": str((OUTPUT_DIR / "effects.npz").relative_to(ROOT)),
            "effects_sha256": sha256(effects_path),
            "predictions": predictions,
            "peak_rss_gib": peak_rss_gib(),
            "elapsed_seconds": time.monotonic() - started,
            "challenge_build_created": False,
            "challenge_submission_created": False,
        }
        report_path = staging / "build-report.json"
        report_path.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        os.replace(staging, OUTPUT_DIR)
        print(
            json.dumps(
                {
                    "bridge_checks": report["bridge_checks"],
                    "incumbent_sha256": predictions["incumbent_sha256"],
                    "candidate_sha256": predictions["candidate_sha256"],
                    "all_pseudobulk_exact": (
                        predictions["all_pseudobulk_exact"]
                    ),
                    "peak_rss_gib": report["peak_rss_gib"],
                    "elapsed_seconds": report["elapsed_seconds"],
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
