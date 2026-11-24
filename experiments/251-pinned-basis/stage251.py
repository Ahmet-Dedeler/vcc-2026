"""251 -- stage the production EFFECT MATRIX from a chosen builder, without emitting.

249's `stage249.py` (= 222's `stage222.py`: build210's flags verbatim, the production operator /
banks / frozen table) plus the pinned-basis verification `run()` performs.  Rows (`stage251.sh`,
each alone, memory-gated):

    pre     the PRE-251 byte copy with VCC_CD4_DEFLATE_BASIS=V5_abc  -> abeb0749... (the flag is unknown to it)
    off     the working-tree builder, flags unset                     -> abeb0749... (222 stage-on028.json = the
                                                                        shipped h1union `final_sha256`)
    k5      the working-tree builder, VCC_CD4_DEFLATE_K=5             -> a28d5be2... (249's staged arm: the re-fit
                                                                        path is untouched)
    pinned  the working-tree builder, VCC_CD4_DEFLATE_BASIS=V5_abc    -> compared with a28d5be2... by compare251.py

Saved per tag: `effects-<tag>.npy`, `cd4-<tag>.npz` (the CD4 residual BEFORE deflation),
`union-<tag>.npz`, `stage-<tag>.json` (env, digests, the whole builder audit incl.
`effects.cd4.deflation`).  Env is read at the builder's IMPORT time: one process = one configuration.

    VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28 [VCC_CD4_DEFLATE_BASIS=V5_abc] \
        uv run python experiments/251-pinned-basis/stage251.py --builder <path> --tag <name>
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
    spec = importlib.util.spec_from_file_location("staged_builder_251", path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / "make_submission_pm_cd4.py")
    sys.modules["staged_builder_251"] = module
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
    basis_name = getattr(mod, "CD4_DEFLATE_BASIS", None)
    if basis_name is not None:
        # what `run()` does before the model load: hash the pins, verify the gene axis
        _v, info = mod.load_cd4_deflation_basis(basis_name, official_genes=genes)
        print(f"[{args_local.tag}] cd4-basis {basis_name}: {info['basis_path']} sha256 {info['basis_sha256']} "
              f"V {info['basis_v_sha256'][:16]} U_k {info['basis_u_k_sha256'][:16]} k={info['k']} "
              f"axis verified {info['gene_axis_verified']}", flush=True)
        del _v

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
    np.savez(
        WORK / f"cd4-{args_local.tag}.npz",
        mean_effect=np.asarray(cd4.mean_effect, dtype=np.float64),
        gene_votes=np.asarray(cd4.gene_votes, dtype=np.uint8),
        target_condition_count=np.asarray(cd4.target_condition_count, dtype=np.uint8),
        targets=np.asarray(targets),
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
        "cd4_deflate_k": getattr(mod, "CD4_DEFLATE_K", None),
        "cd4_deflate_basis": basis_name,
        "effects_npy": str(out_npy),
        "effects_array_sha256": mod.array_sha256(effects),
        "effects_file_sha256": mod.file_sha256(out_npy),
        "operator_sha256": operator.file_sha256,
        "audit": audit,
        "peak_rss_gib": peak_gib(),
        "seconds": round(time.time() - t0, 1),
    }
    (WORK / f"stage-{args_local.tag}.json").write_text(json.dumps(payload, indent=2))
    c = audit["cd4"]
    print(f"[{args_local.tag}] final_sha256 {audit['final_sha256']}")
    print(f"[{args_local.tag}] array_sha256 {payload['effects_array_sha256']}")
    print(f"[{args_local.tag}] cd4.fraction {c['fraction']!r}  residual_median_row_norm_at_unit_fraction "
          f"{c['residual_frame']['residual_median_row_norm_at_unit_fraction']!r}  f_per_unit "
          f"{c['residual_frame']['f_per_unit_fraction']!r}  preserve "
          f"{c['post_add_preserve']['post_add_global_preserve_scale']!r}")
    if "deflation" in c:
        print(f"[{args_local.tag}] cd4.deflation {json.dumps(c['deflation'])}")
    else:
        print(f"[{args_local.tag}] cd4.deflation ABSENT")
    print(f"[{args_local.tag}] done in {payload['seconds']}s, peak "
          f"{payload['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
