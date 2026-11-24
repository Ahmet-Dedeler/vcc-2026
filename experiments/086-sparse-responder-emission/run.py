#!/usr/bin/env python3
"""Guarded HCT116 count-level test of exact-margin sparse responders.

This runner authenticates a terminal, passing experiment 077 result before it
loads donor data, builds count files, or invokes cell-eval2. It reconstructs
the selected exp077 low-alpha profile with the exact source archive from that
run, preserves every original quasi-emitter column total, and changes only the
within-perturbation cell layout.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import resource
import subprocess
import sys
import time
from pathlib import Path
from typing import Mapping

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CONTROL = "non-targeting"
PERT_COL = "target_gene"
N_EMIT = 400
MAX_COUNTS = 1_000_000
FROZEN_FRACTIONS = (0.35, 0.50, 0.65)
FROZEN_MAXIMUM_SCORERS = 2
FROZEN_SEED = 8_602_026
TRUSTED_EXP077_PROVENANCE = HERE / "exp077-trusted-provenance.json"
EXP077_LAUNCH_METADATA_SNAPSHOT = HERE / "exp077-launch-metadata.snapshot.json"
EXPECTED_EXP077_SOURCE_SHA256 = (
    "e1ec0cdb2b5647224431ce334931d9210d1d038494a92660cbff63bcecee3292"
)
EXPECTED_EXP077_EXP031_SHA256 = (
    "ce4c41c715e76a0dc9490a72cbd9924460f747f2e6cbc8b3d542652aa8144ba3"
)
EXPECTED_EXP077_DONOR_TRANSFER_SHA256 = (
    "cecb27035a6c301f96e7fe364032216b1211bb2a92fddd1f9439d9cf8730cb3e"
)
EXPECTED_TEMPLATE_EMITTER_SHA256 = (
    "6333d07ee6f18d62f23422b6af84fb267eb2d5614774c6819109e1f438a1d29d"
)
EXPECTED_GRID_SHA256 = (
    "75a4af2d08ad413ed17dcebfb925b8bbd2af44eceb2990307e08311e3c7a32ac"
)
EXPECTED_REAL_SHA256 = (
    "d20a1e9da59debca8db7d83b955d12c3e1384b530ae6154191d95039b1655863"
)
EXPECTED_WHITENER_SHA256 = (
    "a46264ebb6cb4b877032a6e4c68e76069eb9112c6f1bb7f15e1133fae04a8a77"
)
EXPECTED_WHITENER_METADATA_SHA256 = (
    "cf5c1de31d03875b15f6fd11d3f628aadfccf04e3346493c831c3a6f441b4793"
)
EXPECTED_PANEL_SHA256 = (
    "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"
)
EXPECTED_DONOR_SHA256 = {
    "k562_gwps": {
        "b_ctrl.npy": "d06cf51d51cde1ae999e86bbd3be0bdc15471392bae153a1aa1e21dc5c81b957",
        "delta.npy": "0c9422280fed7afb81ca474e8b9024af788c9ffbde32e11f96cb52422e296aa9",
        "genes.txt": "f499ffc3077e4f06c86cd4504bfeefd34741bc2fd515be93f7e436c5c08fa2ba",
        "perts.txt": "b784821e4a662aa8ab8d55da357fc712f109dda5a464c3e7e047b15b61408bec",
    },
    "orion_hek293t": {
        "b_ctrl.npy": "e04234cffbe102731cc562e8fb06bddb8ef0644b2af894b8a0ee6ce079ca4984",
        "delta.npy": "d95202d984d9322dc3830f23b6ef55001ba8562dea1d72c43fb0c9d5064c7d1b",
        "genes.txt": "db8ba814bd30dcc7fe48be659f4dee3f6955b626d7d7e385f4630232733f893e",
        "perts.txt": "3846edf57cc7a0bd91a499ce9be4ed68c9c8ccc1f36ba062c3bb1301695ef547",
    },
}
ANCHORS = {
    "pds_cosine": (0.500, 0.9525),
    "expr_mse_unbiased_capped_norm": (0.989, 0.0365),
    "de_wilcoxon_lfc_nmae": (1.0013, 0.400),
    "de_wilcoxon_direction_fidelity_yield_raw": (0.5135, 0.8135),
    "de_wilcoxon_direction_reach_raw": (0.072, 0.968),
    "de_wilcoxon_sig_jaccard": (0.029, 0.399),
}


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def sha256_arrays(records: list[tuple[str, np.ndarray]]) -> str:
    digest = hashlib.sha256()
    for label, values in records:
        encoded = label.encode()
        array = np.ascontiguousarray(values)
        digest.update(len(encoded).to_bytes(4, "little"))
        digest.update(encoded)
        digest.update(array.view(np.uint8))
    return digest.hexdigest()


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / 1024**3 if sys.platform == "darwin" else value / 1024**2


def load_local_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def official(raw: Mapping[str, float]) -> dict[str, float]:
    output: dict[str, float] = {}
    for metric, (baseline, replicate) in ANCHORS.items():
        value = (float(raw[metric]) - baseline) / (replicate - baseline)
        if metric == "expr_mse_unbiased_capped_norm":
            value = max(0.0, min(1.0, value))
        output[metric] = float(value)
    output["avg6"] = float(np.mean(list(output.values())))
    return output


def tag_for_config(config: Mapping[str, object]) -> str:
    return "_".join(f"{key}{value}" for key, value in config.items())


def validate_frozen_protocol(args: argparse.Namespace) -> None:
    if tuple(args.fractions) != FROZEN_FRACTIONS:
        raise ValueError(
            f"fractions {tuple(args.fractions)} != frozen {FROZEN_FRACTIONS}"
        )
    if (
        isinstance(args.maximum_scorers, bool)
        or args.maximum_scorers != FROZEN_MAXIMUM_SCORERS
    ):
        raise ValueError(
            f"maximum scorers {args.maximum_scorers} != frozen {FROZEN_MAXIMUM_SCORERS}"
        )
    if isinstance(args.seed, bool) or args.seed != FROZEN_SEED:
        raise ValueError(f"seed {args.seed} != frozen {FROZEN_SEED}")


def read_selected_config(exp077_repo: Path, selected_tag: str) -> dict[str, object]:
    grid = exp077_repo / "experiments/077-low-alpha-union/grid.json"
    if sha256_file(grid) != EXPECTED_GRID_SHA256:
        raise ValueError("exp077 grid hash changed")
    configs = json.loads(grid.read_text())
    matches = [config for config in configs if tag_for_config(config) == selected_tag]
    if len(matches) != 1:
        raise ValueError(
            "selected exp077 tag does not identify exactly one frozen config"
        )
    config = matches[0]
    expected = {
        "g": 0.65,
        "q": 1.0,
        "spread": 0.25,
        "cpm": 5.0,
        "seed": 0,
    }
    for key, value in expected.items():
        if float(config[key]) != value:
            raise ValueError(f"exp077 selected config changed at {key}")
    if float(config["a"]) not in (0.35, 0.40):
        raise ValueError("exp077 selected alpha is not 0.35 or 0.40")
    return config


def verify_frozen_inputs(
    *,
    real_path: Path,
    donor_specs: list[str],
    whitener_path: Path,
    whitener_metadata_path: Path,
    panel_csv: Path,
) -> None:
    fixed_files = {
        real_path: EXPECTED_REAL_SHA256,
        whitener_path: EXPECTED_WHITENER_SHA256,
        whitener_metadata_path: EXPECTED_WHITENER_METADATA_SHA256,
        panel_csv: EXPECTED_PANEL_SHA256,
    }
    for path, expected in fixed_files.items():
        observed = sha256_file(path)
        if observed != expected:
            raise ValueError(f"frozen input {path} SHA-256 {observed} != {expected}")

    parsed: dict[str, Path] = {}
    for specification in donor_specs:
        if "=" not in specification:
            raise ValueError(f"malformed donor specification {specification!r}")
        name, raw_path = specification.split("=", 1)
        if name in parsed:
            raise ValueError(f"duplicate donor specification {name!r}")
        parsed[name] = Path(raw_path)
    if tuple(parsed) != tuple(EXPECTED_DONOR_SHA256):
        raise ValueError(
            f"donor order/set {tuple(parsed)} != {tuple(EXPECTED_DONOR_SHA256)}"
        )
    for donor, files in EXPECTED_DONOR_SHA256.items():
        for filename, expected in files.items():
            path = parsed[donor] / filename
            observed = sha256_file(path)
            if observed != expected:
                raise ValueError(
                    f"frozen donor input {donor}/{filename} SHA-256 "
                    f"{observed} != {expected}"
                )


def control_mean_cpm(matrix: sp.csr_matrix, *, chunk_rows: int = 256) -> np.ndarray:
    output = np.zeros(matrix.shape[1], dtype=np.float64)
    for start in range(0, matrix.shape[0], chunk_rows):
        block = matrix[start : start + chunk_rows].tocsr().astype(np.float64)
        libraries = np.asarray(block.sum(axis=1), dtype=np.float64).ravel()
        if np.any(libraries <= 0):
            raise ValueError("control contains a zero-library cell")
        row_of_nnz = np.repeat(np.arange(block.shape[0]), np.diff(block.indptr))
        block.data *= 1_000_000.0 / libraries[row_of_nnz]
        output += np.asarray(block.sum(axis=0), dtype=np.float64).ravel()
    return output / matrix.shape[0]


def prepare_exact_exp077_geometry(
    *,
    exp077_repo: Path,
    real_path: Path,
    donor_specs: list[str],
    whitener_path: Path,
    whitener_metadata_path: Path,
    panel_csv: Path,
    config: Mapping[str, object],
) -> dict[str, object]:
    exp031_path = exp077_repo / "experiments/031-somatic/run.py"
    donor_transfer_path = exp077_repo / "src/vcc2026/donor_transfer.py"
    if sha256_file(exp031_path) != EXPECTED_EXP077_EXP031_SHA256:
        raise ValueError("exp077 experiment-031 source hash changed")
    if sha256_file(donor_transfer_path) != EXPECTED_EXP077_DONOR_TRANSFER_SHA256:
        raise ValueError("exp077 donor-transfer source hash changed")
    sys.path.insert(0, str(exp077_repo / "src"))
    from vcc2026.donor_transfer import Donor, TS, emit_counts
    from vcc2026.multi_donor import MultiDonorModel
    from vcc2026.whiten import Whitener

    real = ad.read_h5ad(real_path)
    matrix = real.X.tocsr() if not sp.isspmatrix_csr(real.X) else real.X
    labels = real.obs[PERT_COL].astype(str).to_numpy()
    genes = np.asarray(real.var_names, dtype=str)
    control_mask = labels == CONTROL
    controls = matrix[control_mask].tocsr()
    control_totals = np.asarray(controls.sum(axis=0), dtype=np.float64).ravel()
    b_ctrl = np.log1p(TS * control_totals / control_totals.sum())
    control_profile = np.expm1(b_ctrl) / TS
    control_profile /= control_profile.sum()
    libraries = np.asarray(controls.sum(axis=1), dtype=np.int64).ravel()
    targets = tuple(sorted(set(labels) - {CONTROL}))

    donors = {}
    for specification in donor_specs:
        name, path = specification.split("=", 1)
        donors[name] = Donor.load(Path(path))
    if tuple(donors) != ("k562_gwps", "orion_hek293t"):
        raise ValueError("exp077 donor order/set changed")
    model = MultiDonorModel(donors, genes, gene_norm=False, eb=False)

    shared: set[str] = set()
    for name in model.names:
        shared |= set(model.target_genes[model.tcols[name]])
    challenge_targets = set(pd.read_csv(panel_csv)[PERT_COL].astype(str))
    shared -= challenge_targets
    target_index = {gene: index for index, gene in enumerate(model.target_genes)}
    wcols = np.asarray(sorted(target_index[gene] for gene in shared), dtype=np.int64)
    metadata = json.loads(whitener_metadata_path.read_text())
    expected_metadata = {
        "axis": "union",
        "donors": list(model.names),
        "n_cols": int(len(wcols)),
        "wcols_sha256": hashlib.sha256(wcols.tobytes()).hexdigest(),
        "excluded_targets": sorted(challenge_targets),
    }
    if metadata != expected_metadata:
        raise ValueError("exp077 whitener metadata does not match reconstructed axis")
    whitener = Whitener.load(whitener_path)

    effects = np.vstack(
        [
            model.effect(target, 1.0, 0.0, weights={name: 1.0 for name in model.names})
            for target in targets
        ]
    ).astype(np.float64)
    gamma = float(config["g"])
    before = float(np.median(np.linalg.norm(effects[:, wcols], axis=1)))
    effects[:, wcols] = whitener.apply(effects[:, wcols], gamma, 0.0, renorm=False)
    after = float(np.median(np.linalg.norm(effects[:, wcols], axis=1))) or 1.0
    effects[:, wcols] *= before / after
    q = float(config["q"])
    if q != 1.0:
        median = float(np.median(np.linalg.norm(effects, axis=1))) or 1.0
        effects = np.sign(effects) * np.abs(effects) ** q
        effects *= median / (float(np.median(np.linalg.norm(effects, axis=1))) or 1.0)
    effects *= (
        float(config["a"])
        * model.effect_scale
        / (float(np.median(np.linalg.norm(effects, axis=1))) or 1.0)
    )
    return {
        "real": real,
        "matrix": matrix,
        "controls": controls,
        "genes": genes,
        "targets": targets,
        "libraries": libraries,
        "control_profile": control_profile,
        "control_mean_cpm": control_mean_cpm(controls),
        "model": model,
        "b_ctrl": b_ctrl,
        "effects": effects,
        "emit_counts": emit_counts,
    }


def build_candidate(
    *,
    state: Mapping[str, object],
    config: Mapping[str, object],
    fraction: float,
    template_module,
    emitter_module,
    output_path: Path,
    seed: int,
) -> dict[str, object]:
    controls = state["controls"]
    genes = state["genes"]
    targets = state["targets"]
    libraries = state["libraries"]
    control_profile = state["control_profile"]
    control_mean = state["control_mean_cpm"]
    model = state["model"]
    b_ctrl = state["b_ctrl"]
    effects = state["effects"]
    emit_counts = state["emit_counts"]
    assert isinstance(controls, sp.csr_matrix)
    assert isinstance(genes, np.ndarray)
    assert isinstance(libraries, np.ndarray)
    assert isinstance(control_profile, np.ndarray)
    assert isinstance(control_mean, np.ndarray)
    assert isinstance(effects, np.ndarray)

    rng = np.random.default_rng(int(config["seed"]))
    blocks: list[sp.csr_matrix] = []
    observations: list[str] = []
    diagnostics: list[dict[str, object]] = []
    margin_records: list[tuple[str, np.ndarray]] = []
    jackknife_records: list[dict[str, float | str]] = []
    for index, target in enumerate(targets):
        profile = model.profile_from_effect(
            effects[index],
            b_ctrl,
            min_cpm=float(config["cpm"]),
            gate_unexpressed=True,
        )
        picked = rng.choice(len(libraries), N_EMIT, replace=len(libraries) < N_EMIT)
        picked_libraries = libraries[picked]

        # Reproduce the exact original quasi call to preserve its RNG trajectory,
        # then authenticate its integer margins against the closed-form totals.
        data, indices, indptr = emit_counts(
            profile,
            picked_libraries,
            0.0,
            rng,
            mode="quasi",
            spread=np.full(len(genes), float(config["spread"])),
        )
        quasi = sp.csr_matrix((data, indices, indptr), shape=(N_EMIT, len(genes)))
        quasi_totals = np.asarray(quasi.sum(axis=0), dtype=np.int64).ravel()
        library = int(min(np.median(picked_libraries), MAX_COUNTS))
        expected_totals = np.floor(N_EMIT * library * profile).astype(np.int64)
        if not np.array_equal(quasi_totals, expected_totals):
            raise RuntimeError(f"{target}: reconstructed exp077 margins differ")
        grand_total = int(quasi_totals.sum())
        row_margins = template_module.integer_column_totals(
            picked_libraries.astype(np.float64),
            grand_total,
            seed=emitter_module.stable_seed(seed, target, "row-margins"),
        )
        baseline_totals = template_module.integer_column_totals(
            control_profile,
            grand_total,
            seed=emitter_module.stable_seed(seed, target, "control-columns"),
        )
        baseline, baseline_diag = template_module.emit_template_fixed_margins(
            controls[picked],
            column_totals=baseline_totals,
            row_margins=row_margins,
            template_library_sizes=picked_libraries,
            seed=emitter_module.stable_seed(seed, target, "baseline"),
        )
        scored = control_mean > 5.0
        target_position = np.flatnonzero(genes == target)
        if target_position.size:
            scored[int(target_position[0])] = False
        emitted, emission_diag = emitter_module.emit_sparse_responders(
            baseline,
            column_totals=quasi_totals,
            responder_fraction=fraction,
            key=target,
            seed=seed,
            control_reference_totals=baseline_totals,
            scored_gene_mask=scored,
            max_counts_per_cell=MAX_COUNTS,
            minimum_counts_per_cell=1,
        )
        if not np.array_equal(
            np.asarray(emitted.sum(axis=1), dtype=np.int64).ravel(),
            row_margins,
        ):
            raise RuntimeError(f"{target}: row margins changed")
        if not np.array_equal(
            np.asarray(emitted.sum(axis=0), dtype=np.int64).ravel(),
            quasi_totals,
        ):
            raise RuntimeError(f"{target}: exp077 column margins changed")
        quasi_jackknife = emitter_module.jackknife_correction(quasi, chunk_genes=32)
        sparse_jackknife = emitter_module.jackknife_correction(emitted, chunk_genes=32)
        jackknife_tolerance = max(1e-12, abs(quasi_jackknife) * 1e-12)
        jackknife_records.append(
            {
                "target": target,
                "quasi": quasi_jackknife,
                "sparse": sparse_jackknife,
                "delta": sparse_jackknife - quasi_jackknife,
                "non_decreasing": bool(
                    sparse_jackknife + jackknife_tolerance >= quasi_jackknife
                ),
            }
        )
        margin_records.append((target, quasi_totals))
        diagnostics.append(
            {
                "target": target,
                "baseline_template": baseline_diag.to_dict(),
                "sparse": emission_diag.to_dict(),
            }
        )
        blocks.append(emitted)
        observations.extend([target] * N_EMIT)

    blocks.append(controls)
    observations.extend([CONTROL] * controls.shape[0])
    combined = sp.vstack(blocks, format="csr", dtype=np.float32)
    if combined.data.min(initial=0) < 0 or not np.isfinite(combined.data).all():
        raise RuntimeError("candidate matrix contains invalid counts")
    if combined.data.max(initial=0) > np.iinfo(np.int32).max:
        raise RuntimeError("candidate matrix exceeds int32 count capacity")
    row_sums = np.asarray(combined.sum(axis=1), dtype=np.int64).ravel()
    if row_sums.max(initial=0) > MAX_COUNTS:
        raise RuntimeError("candidate matrix exceeds VCC row-count cap")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(
        X=combined,
        obs=pd.DataFrame(
            {PERT_COL: observations},
            index=[f"c{index}" for index in range(len(observations))],
        ),
        var=pd.DataFrame(index=genes),
    ).write_h5ad(output_path, compression="lzf")

    sparse_diags = [record["sparse"] for record in diagnostics]
    max_scored_clip = max(
        float(record["scored_nonresponder_clipped_mass_fraction"])
        for record in sparse_diags
    )
    min_concentration = min(
        float(record["absolute_change_mass_in_responders_fraction"])
        for record in sparse_diags
    )
    all_exact = all(
        int(record["max_abs_column_error"]) == 0
        and bool(record["baseline_row_margins_preserved"])
        for record in sparse_diags
    )
    jackknife_non_decreasing = all(
        bool(record["non_decreasing"]) for record in jackknife_records
    )
    jackknife_deltas = np.asarray(
        [float(record["delta"]) for record in jackknife_records],
        dtype=np.float64,
    )
    construction_gate = {
        "all_target_column_margins_exact": all_exact,
        "maximum_scored_nonresponder_clipped_mass_at_most_0_01": (
            max_scored_clip <= 0.01
        ),
        "minimum_response_mass_concentration_at_least_0_95": (
            min_concentration >= 0.95
        ),
        "all_rows_between_1_and_1000000": bool(
            row_sums.min() >= 1 and row_sums.max() <= MAX_COUNTS
        ),
        "all_target_jackknife_corrections_non_decreasing": (jackknife_non_decreasing),
    }
    construction_gate["passed"] = all(construction_gate.values())
    return {
        "fraction": float(fraction),
        "path": str(output_path),
        "sha256": sha256_file(output_path),
        "shape": list(combined.shape),
        "nnz": int(combined.nnz),
        "margin_sha256": sha256_arrays(margin_records),
        "max_scored_nonresponder_clipped_mass_fraction": max_scored_clip,
        "mean_scored_nonresponder_clipped_mass_fraction": float(
            np.mean(
                [
                    float(record["scored_nonresponder_clipped_mass_fraction"])
                    for record in sparse_diags
                ]
            )
        ),
        "minimum_response_mass_concentration": min_concentration,
        "jackknife_all_targets": {
            "count": len(jackknife_records),
            "non_decreasing_count": int(
                sum(bool(record["non_decreasing"]) for record in jackknife_records)
            ),
            "minimum_delta": float(jackknife_deltas.min()),
            "median_delta": float(np.median(jackknife_deltas)),
            "maximum_delta": float(jackknife_deltas.max()),
            "records": jackknife_records,
        },
        "construction_gate": construction_gate,
        "per_target": diagnostics,
    }


def read_aggregate(path: Path) -> dict[str, float]:
    frame = pd.read_csv(path).set_index("statistic")
    return {metric: float(frame.loc["mean", metric]) for metric in ANCHORS}


def score_candidate(
    *,
    prediction: Path,
    real: Path,
    output: Path,
) -> dict[str, object]:
    output.mkdir(parents=True, exist_ok=False)
    process = subprocess.run(
        [
            "uv",
            "run",
            "cell-eval2",
            "run",
            "-ap",
            str(prediction),
            "-ar",
            str(real),
            "-o",
            str(output),
            "--preset",
            "vcc2026",
            "--pert-col",
            PERT_COL,
            "--control",
            CONTROL,
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    (output / "stdout.log").write_text(process.stdout)
    (output / "stderr.log").write_text(process.stderr)
    aggregate = output / "agg_results.csv"
    if process.returncode != 0 or not aggregate.is_file():
        raise RuntimeError(
            f"cell-eval2 failed ({process.returncode}): "
            f"{(process.stderr or process.stdout)[-1000:]}"
        )
    raw = read_aggregate(aggregate)
    return {"raw": raw, "scaled": official(raw)}


def run(args: argparse.Namespace) -> dict[str, object]:
    validate_frozen_protocol(args)
    if args.output.is_symlink() or args.output.exists():
        raise FileExistsError(f"refusing existing output path {args.output}")
    emitter = load_local_module("exp086_emitter_run", HERE / "emitter.py")
    gate_module = load_local_module("exp086_gate_run", HERE / "gate.py")
    gate = gate_module.validate_terminal_exp077(
        args.exp077_output,
        run_id=args.exp077_run_id,
        trusted_provenance_path=TRUSTED_EXP077_PROVENANCE,
        launch_metadata_snapshot_path=EXP077_LAUNCH_METADATA_SNAPSHOT,
        launch_command_evidence=args.exp077_launch_command_evidence,
        terminal_command_evidence=args.exp077_terminal_command_evidence,
        terminal_invocation_evidence=args.exp077_terminal_invocation_evidence,
    )
    if sha256_file(args.exp077_source_archive) != EXPECTED_EXP077_SOURCE_SHA256:
        raise ValueError("exp077 source archive hash changed")
    if gate.source_archive_sha256 != EXPECTED_EXP077_SOURCE_SHA256:
        raise ValueError("SSM-attested exp077 source archive hash changed")
    if (
        sha256_file(ROOT / "experiments/043-template-emitter/template_emitter.py")
        != EXPECTED_TEMPLATE_EMITTER_SHA256
    ):
        raise ValueError("experiment 043 template emitter hash changed")
    template = load_local_module(
        "exp043_template_for_086",
        ROOT / "experiments/043-template-emitter/template_emitter.py",
    )
    verify_frozen_inputs(
        real_path=args.real,
        donor_specs=args.donors,
        whitener_path=args.whitener,
        whitener_metadata_path=args.whitener_metadata,
        panel_csv=args.panel_csv,
    )
    config = read_selected_config(args.exp077_repo, gate.selected_tag)
    if float(config["a"]) != gate.selected_alpha:
        raise ValueError("exp077 gate alpha and source config disagree")

    state = prepare_exact_exp077_geometry(
        exp077_repo=args.exp077_repo,
        real_path=args.real,
        donor_specs=args.donors,
        whitener_path=args.whitener,
        whitener_metadata_path=args.whitener_metadata,
        panel_csv=args.panel_csv,
        config=config,
    )
    args.output.mkdir(parents=True)
    builds = []
    for fraction in args.fractions:
        tag = f"responders-{fraction:g}"
        build = build_candidate(
            state=state,
            config=config,
            fraction=fraction,
            template_module=template,
            emitter_module=emitter,
            output_path=args.output / tag / "pred.h5ad",
            seed=args.seed,
        )
        builds.append(build)

    margin_hashes = {str(build["margin_sha256"]) for build in builds}
    if len(margin_hashes) != 1:
        raise RuntimeError("candidate fractions reconstructed different exp077 margins")
    eligible = [build for build in builds if build["construction_gate"]["passed"]]
    to_score = sorted(eligible, key=lambda value: float(value["fraction"]))[
        : args.maximum_scorers
    ]
    selected_summary = json.loads((args.exp077_output / "summary.json").read_text())
    exp077_candidate = next(
        candidate
        for candidate in selected_summary["candidates"]
        if candidate["tag"] == gate.selected_tag
    )
    baseline_raw = {
        metric: float(exp077_candidate["raw"][metric]) for metric in ANCHORS
    }
    scored = []
    for build in to_score:
        score = score_candidate(
            prediction=Path(str(build["path"])),
            real=args.real,
            output=Path(str(build["path"])).parent / "score",
        )
        raw = score["raw"]
        pds_delta = raw["pds_cosine"] - baseline_raw["pds_cosine"]
        score["delta_vs_exp077_raw"] = {
            metric: raw[metric] - baseline_raw[metric] for metric in ANCHORS
        }
        score["pseudobulk_pds_invariant_passed"] = abs(pds_delta) <= 1e-12
        score["corrected_mse_non_regression_passed"] = bool(
            raw["expr_mse_unbiased_capped_norm"]
            <= baseline_raw["expr_mse_unbiased_capped_norm"] + 1e-12
        )
        if not score["pseudobulk_pds_invariant_passed"]:
            raise RuntimeError(
                "PDS changed despite exact exp077 column margins; refusing result"
            )
        if not score["corrected_mse_non_regression_passed"]:
            raise RuntimeError(
                "corrected MSE regressed despite the jackknife non-regression gate"
            )
        scored.append(
            {
                "fraction": build["fraction"],
                "prediction_sha256": build["sha256"],
                **score,
            }
        )
    result = {
        "status": "completed" if scored else "construction_gate_failed",
        "experiment": "086-sparse-responder-emission",
        "exp077_gate": gate.to_dict(),
        "exp077_selected_config": config,
        "exp077_uniform_quasi_raw": baseline_raw,
        "fractions_built": list(args.fractions),
        "selection_rule": (
            "score at most two smallest responder fractions passing prediction-only "
            "construction gates; no held truth enters this choice"
        ),
        "builds": builds,
        "scored": scored,
        "submission_authorized": False,
        "peak_rss_gib": peak_rss_gib(),
        "elapsed_seconds": time.time() - args.started,
    }
    (args.output / "results.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    state["real"].file.close() if getattr(state["real"], "isbacked", False) else None
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--exp077-output", type=Path, required=True)
    parser.add_argument("--exp077-run-id", required=True)
    parser.add_argument("--exp077-source-archive", type=Path, required=True)
    parser.add_argument("--exp077-repo", type=Path, required=True)
    parser.add_argument("--exp077-launch-command-evidence", type=Path, required=True)
    parser.add_argument("--exp077-terminal-command-evidence", type=Path, required=True)
    parser.add_argument(
        "--exp077-terminal-invocation-evidence", type=Path, required=True
    )
    parser.add_argument("--real", type=Path, required=True)
    parser.add_argument("--donors", nargs="+", required=True)
    parser.add_argument("--whitener", type=Path, required=True)
    parser.add_argument("--whitener-metadata", type=Path, required=True)
    parser.add_argument("--panel-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--fractions", nargs="+", type=float, default=list(FROZEN_FRACTIONS)
    )
    parser.add_argument("--maximum-scorers", type=int, default=FROZEN_MAXIMUM_SCORERS)
    parser.add_argument("--seed", type=int, default=FROZEN_SEED)
    parser.add_argument("--guarded", action="store_true")
    args = parser.parse_args()
    if not args.guarded:
        raise SystemExit("refusing real-data execution without --guarded")
    try:
        validate_frozen_protocol(args)
    except ValueError as error:
        raise SystemExit(str(error)) from error
    args.started = time.time()
    result = run(args)
    print(
        json.dumps(
            {
                "status": result["status"],
                "exp077_selected_alpha": result["exp077_gate"]["selected_alpha"],
                "eligible_fractions": [
                    build["fraction"]
                    for build in result["builds"]
                    if build["construction_gate"]["passed"]
                ],
                "scored_fractions": [record["fraction"] for record in result["scored"]],
                "submission_authorized": False,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
