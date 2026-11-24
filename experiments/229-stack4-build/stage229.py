"""229 — stage the production EFFECT MATRIX from a chosen builder, without emitting.

222's `stage222.py`, copied byte-for-byte except this paragraph and the two lines that name
it.  `WORK` is `Path(__file__).parent/"work"`, so the copy writes into 229 and never into 222.
It reads `mod.H1_UNION_G` and `mod.prepare_h1_union` off the imported module, so it drives the
registry-driven builder unchanged; `VCC_UNION_SOURCES` is read at the builder's import time,
like every other flag, so one process is one configuration.

The effect matrix is 300 x 18,533 float64 = 44 MB; the 12.6 GB h5ad that follows it is 45
minutes of emission.  Every question this experiment has about the arm's GEOMETRY -- is the
default path bit-identical to the pre-change builder, is the residual non-vacuous, does the
union reproduce 218 -- is a question about that 44 MB object, so it is answered here first and
the two long builds are only spent once the geometry is proved.

KNOWN-ANSWER ROW 2 (mandate).  `--builder` takes a PATH, so the same driver runs the builder
resolved out of **git** (`work/make_submission_pm_cd4.PRE.py`, `git show HEAD:...`) and the
builder in the working tree, in two separate processes, and the two effect matrices are
compared byte-for-byte.  A test that only exercises the new code is not a control.

Env is read at the builder's IMPORT time, so one process = one configuration.

    uv run python stage222.py --builder <path> --tag <name> [--h1-union-g 0.28]
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
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

# build210.sh's flags, verbatim -- the incumbent `cd4f034`'s configuration.
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
    """Import a builder from an arbitrary path, with its `ROOT` still the repository root.

    The builder computes `ROOT = Path(__file__).resolve().parents[1]`, which is correct only
    for a file sitting in `scripts/`.  The pre-change control is a byte copy of
    `git show HEAD:scripts/make_submission_pm_cd4.py` living in this experiment's `work/`, so
    `__file__` is pinned to the path the bytes came FROM before the module body runs.  Only
    `ROOT` and `collect_implementation_hashes` read it; the executed bytes are git's.
    """
    spec = importlib.util.spec_from_file_location("staged_builder_229", path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / path.name)
    sys.modules["staged_builder_229"] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--builder", required=True)
    ap.add_argument("--tag", required=True)
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
    print(f"[{args_local.tag}] contexts {list(bundle_contexts)}  "
          f"{len(targets)} targets  {len(genes)} genes", flush=True)

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
    )
    gene_index = {gene: index for index, gene in enumerate(genes)}
    panel_positions = np.asarray([gene_index[t] for t in targets], dtype=np.int64)
    cd4_banks, cd4 = mod.prepare_frozen_cd4(
        cd4_root, targets, genes, targets,
        chunk_rows=args.chunk_rows, max_rss_gib=args.max_rss_gib,
    )
    print(f"[{args_local.tag}] inputs staged ({time.time() - t0:.0f}s, "
          f"peak {peak_gib():.2f} GiB)", flush=True)

    extra = {}
    if getattr(mod, "H1_UNION_G", None) is not None:
        h1_union, details = mod.prepare_h1_union(
            targets, genes, targets, max_rss_gib=args.max_rss_gib
        )
        audit = mod.h1_union_support_audit(h1_union, details)
        extra = {"h1_union": h1_union, "h1_union_audit": audit}
        print(f"[{args_local.tag}] h1 union g={mod.H1_UNION_G:g}  "
              f"{audit['carrying_targets']}/{len(targets)} targets carried, "
              f"{audit['targets_carried_by_both_sources']} by both", flush=True)
        # The union residual on its own, saved so `repro229` and the falsifiers can read
        # exactly the object the build used rather than rebuilding it.
        np.savez(
            WORK / f"union-{args_local.tag}.npz",
            mean_effect=np.asarray(h1_union.mean_effect, dtype=np.float64),
            gene_votes=np.asarray(h1_union.gene_votes, dtype=np.uint8),
            target_condition_count=np.asarray(
                h1_union.target_condition_count, dtype=np.uint8),
        )

    effects, audit = mod.build_production_effects(
        model, targets, operator.whitener, operator.columns, panel_positions,
        pm=pm, cd4=cd4, weights=weights, **extra,
    )
    out_npy = WORK / f"effects-{args_local.tag}.npy"
    np.save(out_npy, np.ascontiguousarray(effects, dtype=np.float64))
    payload = {
        "tag": args_local.tag,
        "builder": str(builder_path),
        "builder_sha256": mod.file_sha256(builder_path),
        "env": {k: v for k, v in sorted(os.environ.items()) if k.startswith("VCC_")},
        "alpha": mod.ALPHA,
        "cd4_target_f": mod.CD4_TARGET_F,
        "h1_union_g": getattr(mod, "H1_UNION_G", None),
        "effects_npy": str(out_npy),
        "effects_array_sha256": mod.array_sha256(effects),
        "effects_file_sha256": mod.file_sha256(out_npy),
        "operator_sha256": operator.file_sha256,
        "audit": audit,
        "peak_rss_gib": peak_gib(),
        "seconds": round(time.time() - t0, 1),
    }
    (WORK / f"stage-{args_local.tag}.json").write_text(json.dumps(payload, indent=2))
    print(f"[{args_local.tag}] final_sha256 {audit['final_sha256']}")
    print(f"[{args_local.tag}] array_sha256 {payload['effects_array_sha256']}")
    print(f"[{args_local.tag}] done in {payload['seconds']}s, peak "
          f"{payload['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
