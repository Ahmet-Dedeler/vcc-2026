"""248 -- stage the production EFFECT MATRIX from a chosen builder, without emitting.

222's `stage222.py` (= 229's `stage229.py` = 249's `stage249.py`), with two additions: it
drives the SECOND union when the imported builder carries `UNION2_SOURCES`, and it saves the
union residual(s) and the CD4 residual next to the effect matrix so `repro248` / the profile
falsifiers read exactly the objects the build used.  `WORK` is this experiment's, so nothing
is written into 222/229/249.

Env is read at the builder's IMPORT time, so one process = one configuration.

    uv run python stage248.py --builder <path> --tag <name>
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

# build210.sh's flags, verbatim -- the incumbent's configuration (222 / 229 / 249 unchanged).
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
    spec = importlib.util.spec_from_file_location("staged_builder_248", path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / path.name)
    sys.modules["staged_builder_248"] = module
    spec.loader.exec_module(module)
    return module


def save_residual(path: Path, res, targets) -> None:
    np.savez(path,
             mean_effect=np.asarray(res.mean_effect, dtype=np.float64),
             gene_votes=np.asarray(res.gene_votes, dtype=np.uint8),
             target_condition_count=np.asarray(res.target_condition_count, dtype=np.uint8),
             targets=np.asarray(targets))


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
    save_residual(WORK / f"cd4-{args_local.tag}.npz", cd4, targets)
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
        save_residual(WORK / f"union-{args_local.tag}.npz", h1_union, targets)
    # 248: the second union, driven exactly as `run()` drives it.
    union2_sources = getattr(mod, "UNION2_SOURCES", None)
    if union2_sources is not None:
        union2, details2 = mod.prepare_h1_union(
            targets, genes, targets, max_rss_gib=args.max_rss_gib, sources=union2_sources
        )
        audit2 = mod.h1_union_support_audit(
            union2, details2, sources=union2_sources, target_g=mod.UNION2_G
        )
        extra.update({"union2": union2, "union2_audit": audit2})
        print(f"[{args_local.tag}] union2 g={mod.UNION2_G:g}  sources "
              f"{[n for n, _ in union2_sources]}  {audit2['carrying_targets']}/{len(targets)} "
              f"targets carried, by source count {audit2['targets_by_source_count']}", flush=True)
        save_residual(WORK / f"union2-{args_local.tag}.npz", union2, targets)

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
        "union2_g": getattr(mod, "UNION2_G", None),
        "union2_sources": None if union2_sources is None else [n for n, _ in union2_sources],
        "effects_npy": str(out_npy),
        "effects_array_sha256": mod.array_sha256(effects),
        "effects_file_sha256": mod.file_sha256(out_npy),
        "operator_sha256": operator.file_sha256,
        "audit": audit,
        "peak_rss_gib": peak_gib(),
        "seconds": round(time.time() - t0, 1),
    }
    (WORK / f"stage-{args_local.tag}.json").write_text(json.dumps(payload, indent=2, default=str))
    print(f"[{args_local.tag}] final_sha256 {audit['final_sha256']}")
    print(f"[{args_local.tag}] array_sha256 {payload['effects_array_sha256']}")
    print(f"[{args_local.tag}] cd4.fraction {audit['cd4']['fraction']!r}")
    if "h1_union" in audit:
        print(f"[{args_local.tag}] h1_union.solved_fraction "
              f"{audit['h1_union']['residual_frame']['solved_fraction']!r}")
    if "union2" in audit:
        print(f"[{args_local.tag}] union2.solved_fraction "
              f"{audit['union2']['residual_frame']['solved_fraction']!r}  carrying "
              f"{audit['union2']['residual_frame']['carrying_rows']}")
    print(f"[{args_local.tag}] done in {payload['seconds']}s, peak "
          f"{payload['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
