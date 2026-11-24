#!/usr/bin/env python3
"""Export the staged a15 pre-emission effect matrix in its pinned build runtime."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np


EXPECTED_PREWHITEN_SHA256 = (
    "0601336fdf1ec5c533dcad69f6cb6c09d6368182bc8c69d10d516a1ecb19bf65"
)
EXPECTED_FINAL_SHA256 = (
    "2980d70655f2a9d4a81aa0d8cdb06f2389fca20dddf586d64c132cb05885e5f2"
)


def load_module(path: Path):
    os.environ["VCC_ALPHA"] = "0.982"
    os.environ["VCC_CD4_FRACTION"] = "0.05"
    os.environ["VCC_CD4_TARGET_F"] = "0.10"
    spec = importlib.util.spec_from_file_location("exp103_a15_export", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.assert_frozen_contract()
    return module


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production-script", type=Path, required=True)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--k562", type=Path, required=True)
    parser.add_argument("--hct116", type=Path, required=True)
    parser.add_argument("--hek293t", type=Path, required=True)
    parser.add_argument("--perturbmatch", type=Path, required=True)
    parser.add_argument("--cd4-root", type=Path, required=True)
    parser.add_argument("--whitener", type=Path, required=True)
    parser.add_argument("--whitener-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    module = load_module(args.production_script.resolve())
    spec = module.ChallengeSpec.load(args.bundle.resolve())
    donors = {
        "k562_gwps": module.Donor.load(args.k562.resolve()),
        "orion_hct116": module.Donor.load(args.hct116.resolve()),
        "orion_hek293t": module.Donor.load(args.hek293t.resolve()),
    }
    model = module.MultiDonorModel(
        donors,
        spec.genes,
        gene_norm=False,
        eb=False,
    )
    weights = module.donor_weights(model)
    pm = module.EXP061.prepare_perturbmatch(
        args.perturbmatch.resolve(),
        tuple(map(str, spec.genes)),
        chunk_rows=128,
        max_rss_gib=120.0,
    )
    operator = module.load_production_operator(
        model,
        tuple(map(str, spec.perts)),
        args.whitener.resolve(),
        args.whitener_metadata.resolve(),
        expected_sha256=module.WHITENER_REGISTRY["union_pruned"][0],
    )
    positions = {str(gene): index for index, gene in enumerate(spec.genes)}
    panel_positions = np.asarray(
        [positions[str(target)] for target in spec.perts], dtype=np.int64
    )
    _, cd4 = module.prepare_frozen_cd4(
        args.cd4_root.resolve(),
        tuple(map(str, spec.perts)),
        tuple(map(str, spec.genes)),
        tuple(map(str, spec.perts)),
        chunk_rows=128,
        max_rss_gib=120.0,
    )
    effects, audit = module.build_production_effects(
        model,
        tuple(map(str, spec.perts)),
        operator.whitener,
        operator.columns,
        panel_positions,
        pm=pm,
        cd4=cd4,
        alpha=0.982,
        weights=weights,
    )
    if audit["candidate_pre_sha256"] != EXPECTED_PREWHITEN_SHA256:
        raise ValueError("a15 pre-whitening hash mismatch")
    if audit["final_sha256"] != EXPECTED_FINAL_SHA256:
        raise ValueError(
            f"a15 final hash {audit['final_sha256']} != {EXPECTED_FINAL_SHA256}"
        )

    output = args.output.resolve()
    metadata = Path(str(output) + ".json")
    if output.exists() or metadata.exists():
        raise FileExistsError(f"refusing to overwrite {output} or {metadata}")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("wb") as handle:
        np.save(handle, np.asarray(effects, dtype="<f8"), allow_pickle=False)
    record = {
        "shape": list(effects.shape),
        "dtype": "<f8",
        "prewhiten_sha256": audit["candidate_pre_sha256"],
        "final_sha256": audit["final_sha256"],
        "file_sha256": module.file_sha256(output),
        "runtime": {
            "python": sys.version,
            "numpy": np.__version__,
        },
    }
    metadata.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps(record, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
