"""247 item C, the byte-inertness rows — stage the production EFFECT MATRIX, without emitting.

222's `stage222.py`, unchanged in construction (build210's flags verbatim, the production
operator/banks/frozen table), so its two staged digests are the known answers:

    flag OFF, VCC_H1_UNION_G unset    effects_array_sha256 87a4b6d1...  (222 stage-off.json)
    flag OFF, VCC_H1_UNION_G=0.28     effects_array_sha256 abeb0749...  (222 stage-on028.json)
                                      == the shipped h1union's `final_sha256`
    flag ON at g_cd4, G=0.28          effects.cd4.fraction == 0.32207201067578145 to <= 1e-9

Env is read at the builder's IMPORT time, so one process = one configuration; the driver
`stage247.sh` runs the three, each alone under `/usr/bin/time -l`, gated on memory.

    VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 [VCC_H1_UNION_G=0.28] [VCC_CD4_SIZING=carrying VCC_CD4_TARGET_G=...] \
        uv run python experiments/247-oct22-guards/stage247.py --tag <name>
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work" / "stage"
WORK.mkdir(parents=True, exist_ok=True)

BUILD_FLAGS = [
    "--whitener", "data/donors/whitener_union_pruned.npz",
    "--whitener-metadata", "data/donors/whitener_union_pruned.json",
    "--whitener-id", "union_pruned",
    "--min-cpm", "10",
    "--tail-spread", "1.0",
    "--keep-top", "8000",
    "--keep-spread", "0.25",
    "--overdispersion", "4.0",
    "--max-rss-gib", "24",
    "--emitter", "rc_callset",
]


def peak_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def load_builder(path: Path):
    spec = importlib.util.spec_from_file_location("staged_builder_247", path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / path.name)
    sys.modules["staged_builder_247"] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--builder", default=str(ROOT / "scripts/make_submission_pm_cd4.py"))
    ap.add_argument("--tag", required=True)
    ap.add_argument("--save", action="store_true", help="also save the 44 MB effects array as .npy")
    args_local = ap.parse_args()

    t0 = time.time()
    builder_path = Path(args_local.builder).resolve()
    mod = load_builder(builder_path)
    os.chdir(ROOT)

    args = mod.build_parser().parse_args(BUILD_FLAGS)
    bundle = args.bundle.expanduser().resolve()
    donor_paths = {
        mod.K562: args.k562.expanduser().resolve(),
        "orion_hct116": args.hct116.expanduser().resolve(),
        "orion_hek293t": args.hek293t.expanduser().resolve(),
    }
    perturbmatch = args.perturbmatch.expanduser().resolve()
    cd4_root = args.cd4_root.expanduser().resolve()
    operator_path = args.whitener.expanduser().resolve()
    operator_metadata_path = args.whitener_metadata.expanduser().resolve()

    frozen_bundle = mod.load_frozen_bundle(args.frozen_inputs.expanduser().resolve())
    bundle_contexts = mod.read_bundle_contexts(bundle)
    mod.collect_input_hashes(
        bundle, bundle_contexts, donor_paths, perturbmatch, cd4_root,
        operator_path, operator_metadata_path,
        whitener_id=args.whitener_id, frozen_bundle=frozen_bundle,
    )
    spec = mod.ChallengeSpec.load(bundle)
    targets = tuple(map(str, spec.perts))
    genes = tuple(map(str, spec.genes))
    print(f"[{args_local.tag}] contexts {list(bundle_contexts)}  {len(targets)} targets  "
          f"{len(genes)} genes  sizing={mod.CD4_SIZING!r} g={mod.CD4_TARGET_G!r} "
          f"f={mod.CD4_TARGET_F!r} h1={getattr(mod, 'H1_UNION_G', None)!r}", flush=True)

    donors = {name: mod.Donor.load(path) for name, path in donor_paths.items()}
    model = mod.MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False)
    weights = mod.donor_weights(model)
    pm = mod.EXP061.prepare_perturbmatch(
        perturbmatch, genes, chunk_rows=args.chunk_rows, max_rss_gib=args.max_rss_gib
    )
    operator = mod.load_production_operator(
        model, targets, operator_path, operator_metadata_path,
        expected_sha256=(args.expected_whitener_sha256
                         or mod.WHITENER_REGISTRY[args.whitener_id][0]),
        whitener_id=args.whitener_id,
    )
    gene_index = {gene: index for index, gene in enumerate(genes)}
    panel_positions = np.asarray([gene_index[t] for t in targets], dtype=np.int64)
    cd4_banks, cd4 = mod.prepare_frozen_cd4(
        cd4_root, targets, genes, targets,
        chunk_rows=args.chunk_rows, max_rss_gib=args.max_rss_gib,
    )
    print(f"[{args_local.tag}] inputs staged ({time.time() - t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)

    extra = {}
    if getattr(mod, "H1_UNION_G", None) is not None:
        h1_union, details = mod.prepare_h1_union(targets, genes, targets, max_rss_gib=args.max_rss_gib)
        audit = mod.h1_union_support_audit(h1_union, details)
        extra = {"h1_union": h1_union, "h1_union_audit": audit}

    effects, audit = mod.build_production_effects(
        model, targets, operator.whitener, operator.columns, panel_positions,
        pm=pm, cd4=cd4, weights=weights, **extra,
    )
    payload = {
        "tag": args_local.tag,
        "builder": str(builder_path),
        "builder_sha256": mod.file_sha256(builder_path),
        "env": {k: v for k, v in sorted(os.environ.items()) if k.startswith("VCC_")},
        "alpha": mod.ALPHA,
        "cd4_target_f": mod.CD4_TARGET_F,
        "cd4_sizing": mod.CD4_SIZING,
        "cd4_target_g": mod.CD4_TARGET_G,
        "h1_union_g": getattr(mod, "H1_UNION_G", None),
        "effects_array_sha256": mod.array_sha256(effects),
        "operator_sha256": operator.file_sha256,
        "operator_axis_sha256": operator.axis_sha256,
        "cd4_fraction": audit["cd4"]["fraction"],
        "cd4_residual_frame": audit["cd4"]["residual_frame"],
        "cd4_post_add_preserve": audit["cd4"]["post_add_preserve"],
        "final_sha256": audit["final_sha256"],
        "post_whiten_sha256": audit["post_whiten_sha256"],
        "peak_rss_gib": peak_gib(),
        "seconds": round(time.time() - t0, 1),
    }
    if args_local.save:
        np.save(WORK / f"effects-{args_local.tag}.npy", np.ascontiguousarray(effects, dtype=np.float64))
    (WORK / f"stage-{args_local.tag}.json").write_text(json.dumps(payload, indent=2))
    print(f"[{args_local.tag}] array_sha256 {payload['effects_array_sha256']}")
    print(f"[{args_local.tag}] cd4.fraction {payload['cd4_fraction']!r}")
    print(f"[{args_local.tag}] done in {payload['seconds']}s, peak {payload['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
