#!/usr/bin/env python3
"""Build fixed held-HCT dual-moment predictions from Experiment 149 code."""

from __future__ import annotations

import gc
import hashlib
import importlib.util
import json
import os
import shutil
import sys
import time
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PROTOCOL_PATH = HERE / "preregistered_protocol.json"
OUTPUT_DIR = HERE / "work/source-build-v1"
REAL_PATH = ROOT / "data/harness/orion-hct116-thinned.h5ad"
HARNESS_METADATA = REAL_PATH.with_suffix(".h5ad.meta.json")
REFERENCE_104 = ROOT / (
    "experiments/104-context-heldout-credit-adapter/"
    "work-v2/hct116/frozen.npz"
)
REFERENCE_106 = ROOT / (
    "experiments/106-state-baseline/work/runs/20260830T064800Z/"
    "incumbent-a14/hct116/frozen.npz"
)
UNPRUNED_OPERATOR = ROOT / (
    "experiments/062-cd4-post-whiten/cloud-output/whiteners/"
    "held_orion_hct116.npz"
)
PRUNED_OPERATOR = ROOT / (
    "experiments/106-state-baseline/operators/"
    "held_orion_hct116_pruned.npz"
)


def load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


BASE = load_module(
    "exp150_base",
    ROOT / "experiments/149-dual-moment-hek293t/build_predictions.py",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_protocol() -> dict[str, Any]:
    protocol = json.loads(PROTOCOL_PATH.read_text())
    if protocol.get("status") != "preregistered_source_only":
        raise ValueError("experiment is not source-only preregistered")
    if protocol["builder"]["sha256"] != sha256(Path(__file__).resolve()):
        raise ValueError("builder source differs from preregistration")
    expected = protocol["source_bindings_sha256"]
    observed = {
        str((ROOT / relative).relative_to(ROOT)): sha256(ROOT / relative)
        for relative in expected
    }
    if observed != expected:
        raise ValueError("source bindings differ from preregistration")
    return protocol


def verify_biological_inputs() -> dict[str, str]:
    if (
        sha256(BASE.E106.A14_REFERENCE_META)
        != BASE.E106.A14_REFERENCE_META_SHA256
    ):
        raise ValueError("a15 reference manifest changed")
    reference = json.loads(BASE.E106.A14_REFERENCE_META.read_text())
    expected = reference["hashes"]["inputs"]
    observed: dict[str, str] = {}
    for key, path in BASE.E106.a14_biological_inputs("hct116").items():
        actual = sha256(path)
        if expected.get(key) != actual:
            raise ValueError(f"a15 biological input changed: {key}")
        observed[str(path.relative_to(ROOT))] = actual
    return observed


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

    original_alias = dict(BASE.E104.SURFACES["hek293t"])
    BASE.E104.SURFACES["hek293t"] = dict(BASE.E104.SURFACES["hct116"])
    BASE.REAL_PATH = REAL_PATH
    try:
        reference_104 = BASE.load_reference(REFERENCE_104)
        reference_106 = BASE.load_reference(REFERENCE_106)
        genes = np.asarray(reference_104["genes"]).astype(str)
        targets = np.asarray(reference_104["targets"]).astype(str)
        bctrl = np.asarray(reference_104["bctrl"], dtype=np.float64)
        if not np.array_equal(
            genes,
            np.asarray(reference_106["genes"]).astype(str),
        ):
            raise ValueError("frozen HCT reference gene axes differ")
        if not np.array_equal(
            targets,
            np.asarray(reference_106["targets"]).astype(str),
        ):
            raise ValueError("frozen HCT reference target axes differ")
        metadata = json.loads(HARNESS_METADATA.read_text())
        if set(map(str, metadata["perts"])) != set(targets):
            raise ValueError("safe harness metadata target set differs")

        pre, banks, cd4, source_audit = BASE.prepare_source_components(targets)
        unpruned_065, unpruned_audit = BASE.build_effects(
            pre=pre,
            banks=banks,
            cd4=cd4,
            genes=genes,
            operator_path=UNPRUNED_OPERATOR,
            gamma=0.65,
            cd4_fraction=BASE.E104.CD4_FRACTION,
        )
        unpruned_match = np.array_equal(
            unpruned_065.astype(np.float32),
            np.asarray(reference_104["baseline"], dtype=np.float32),
        )
        if not unpruned_match:
            raise RuntimeError(
                "failed Experiment 104 HCT reconstruction bridge"
            )
        del unpruned_065
        gc.collect()

        pruned_065_cd4f005, pruned_005_audit = BASE.build_effects(
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
                "failed Experiment 106 HCT reconstruction bridge"
            )
        del pruned_065_cd4f005
        gc.collect()

        incumbent_effects, incumbent_audit = BASE.build_effects(
            pre=pre,
            banks=banks,
            cd4=cd4,
            genes=genes,
            operator_path=PRUNED_OPERATOR,
            gamma=0.65,
            cd4_fraction=BASE.E104.CD4_FRACTION,
        )
        direction_effects, direction_audit = BASE.build_effects(
            pre=pre,
            banks=banks,
            cd4=cd4,
            genes=genes,
            operator_path=PRUNED_OPERATOR,
            gamma=0.45,
            cd4_fraction=BASE.E104.CD4_FRACTION,
        )
        BASE.memory_checkpoint("effects ready")
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
        predictions = BASE.write_predictions(
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
        BASE.memory_checkpoint("predictions written")
        report = {
            "experiment": "150-dual-moment-score-level-hct116",
            "stage": "source_only_build",
            "provenance_utc": "2026-08-31T18:45:00+00:00",
            "protocol_sha256": sha256(PROTOCOL_PATH),
            "builder_sha256": sha256(Path(__file__).resolve()),
            "real_harness_opened_for_control_rows_only": True,
            "real_harness_control_rows_read": 18_400,
            "real_harness_perturbation_rows_read": 0,
            "real_harness_path_reserved_for_scoring": str(
                REAL_PATH.relative_to(ROOT)
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
            "peak_rss_gib": BASE.peak_rss_gib(),
            "elapsed_seconds": time.monotonic() - started,
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
    finally:
        BASE.E104.SURFACES["hek293t"] = original_alias


if __name__ == "__main__":
    raise SystemExit(main())
