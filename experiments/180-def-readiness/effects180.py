"""180 gate 1 — the rule-10 bit-identity proof for the context adapter, on the REAL A/B/C run.

The item-1 change is entirely upstream of the effect matrix: it moves the bundle half of the
frozen-input table into `scripts/frozen_inputs.json`, builds the bundle file list from the
manifest's contexts instead of a hardcoded A/B/C, and reads the PerturbMatch / CD4 coverage
pair out of the same table. If ANY of that reached the model, the effects audit would move.

So this drives the production `run()` end to end at a23's exact configuration with only
`write_submission_from_effects` stubbed out -- the 3-hour, 19 GB emission is the one thing
being skipped, and it is covered separately, byte for byte, by
`scripts/test_make_submission_pm_cd4.py::test_default_emitter_reproduces_the_pre_change_builder_byte_for_byte`
which drives the real writer against the last pre-change revision out of git.

Everything else runs for real: `load_frozen_bundle`, the three-way context cross-check,
`collect_input_hashes` over all 34 frozen inputs, the donor model, PerturbMatch, the pruned
operator, the whitening, the CD4 residual, and the panel-support assertion.

Compared field by field against `submissions/a23-massfix.h5ad.meta.json`:
`effects.*` (including `final_sha256`), `operator.*`, `hashes.inputs`, `hashes.implementation`
and `config`. The expected diff is EXACTLY: the builder's own file hash (it was patched) and
nothing else.

It also writes the final effect matrix to `work/effects-a23.npy` for `knee180.py`, so the knee
is solved on a matrix whose sha256 has been proven equal to the incumbent's.

    uv run python effects180.py
"""
from __future__ import annotations

import argparse
import contextlib
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)

# a23's flags, set BEFORE the builder module is imported (read at module scope).
os.environ.setdefault("VCC_ALPHA", "0.982")
os.environ.setdefault("VCC_CD4_TARGET_F", "0.10")
os.environ.setdefault("VCC_MASS_FIX", "minnorm")

A23_META = REPO / "submissions/a23-massfix.h5ad.meta.json"
BUILDER = REPO / "scripts/make_submission_pm_cd4.py"


def load_builder():
    spec = importlib.util.spec_from_file_location("prod_builder_180", BUILDER)
    module = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_180"] = module
    spec.loader.exec_module(module)
    return module


def flatten(node, prefix=""):
    out = {}
    if isinstance(node, dict):
        for key, value in node.items():
            out.update(flatten(value, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(node, list):
        out[prefix] = json.dumps(node)
    else:
        out[prefix] = node
    return out


def main() -> int:
    # The builder's donor / PerturbMatch / CD4 defaults are repo-relative paths.
    os.chdir(REPO)
    module = load_builder()
    if module.ALPHA != 0.982 or module.MASS_FIX != "minnorm" or module.CD4_TARGET_F != 0.10:
        raise RuntimeError("a23's env was not applied")
    if module.ALPHA_MODE != "fixed":
        raise RuntimeError("this gate must run at a FIXED alpha, not at the knee")

    captured: dict[str, object] = {}
    real_writer = module.write_submission_from_effects

    def stub(spec, effects, bundle, output, **kwargs):
        captured["effects"] = np.asarray(effects, dtype=np.float64)
        captured["context_effects"] = kwargs.get("context_effects")
        captured["writer_kwargs"] = {
            key: (None if key == "context_effects" else value)
            for key, value in kwargs.items()
            if not isinstance(value, np.ndarray)
        }
        return {
            "path": str(output), "bytes": 0, "sha256": "stubbed", "n_rows": spec.n_rows,
            "nnz": 0, "nnz_per_cell": 0.0, "pct_of_density_cap": 0.0, "contexts": {},
        }

    module.write_submission_from_effects = stub
    args = module.build_parser().parse_args([
        "--bundle", str(REPO / "data/challenge"),
        "-o", str(WORK / "effects-only-never-written.h5ad"),
        "--min-cpm", "10.0",
        "--tail-spread", "1.0",
        "--keep-top", "8000",
        "--keep-spread", "0.25",
        "--overdispersion", "4.0",
        "--whitener-id", "union_pruned",
        "--whitener", str(REPO / "data/donors/whitener_union_pruned.npz"),
        "--whitener-metadata", str(REPO / "data/donors/whitener_union_pruned.json"),
        "--max-rss-gib", "26.0",
    ])
    report = module.run(args)
    module.write_submission_from_effects = real_writer

    np.save(WORK / "effects-a23.npy", captured["effects"])
    (WORK / "effects180-meta.json").write_text(json.dumps(report, indent=2, sort_keys=True))

    incumbent = json.loads(A23_META.read_text())
    fields = 0
    diffs = []
    for block in ("effects", "operator", "config", "hashes"):
        got = flatten(report.get(block, {}), block)
        want = flatten(incumbent.get(block, {}), block)
        for key in sorted(set(got) | set(want)):
            fields += 1
            if got.get(key, "<absent>") != want.get(key, "<absent>"):
                diffs.append((key, got.get(key, "<absent>"), want.get(key, "<absent>")))

    print(f"\n[gate1] fields compared: {fields}")
    effects_fields = sum(1 for k in flatten(report["effects"], "effects"))
    operator_fields = sum(1 for k in flatten(report["operator"], "operator"))
    print(f"[gate1]   effects.*  {effects_fields}   operator.*  {operator_fields}")
    for key, got, want in diffs:
        print(f"[gate1] DIFF {key}\n           built     {got}\n           incumbent {want}")
    if not diffs:
        print("[gate1] differences: 0")
    print(f"\n[gate1] effects.final_sha256 {report['effects']['final_sha256']}")
    print(f"[gate1] a23 final_sha256     {incumbent['effects']['final_sha256']}")
    print(f"[gate1] peak RSS {report['peak_rss_gib']:.3f} GiB")

    expected_only = {"hashes.implementation.scripts/make_submission_pm_cd4.py"}
    unexpected = [key for key, _, _ in diffs if key not in expected_only]
    verdict = "PASS" if not unexpected else "FAIL"
    print(f"\n[gate1] {verdict}: {len(diffs)} difference(s), "
          f"{len(unexpected)} outside the patched-file hash")
    (WORK / "gate1.json").write_text(json.dumps({
        "fields_compared": fields,
        "effects_fields": effects_fields,
        "operator_fields": operator_fields,
        "differences": [{"key": k, "built": g, "incumbent": w} for k, g, w in diffs],
        "unexpected": unexpected,
        "verdict": verdict,
        "final_sha256": report["effects"]["final_sha256"],
        "a23_final_sha256": incumbent["effects"]["final_sha256"],
        "frozen_bundle": report.get("frozen_bundle"),
        "peak_rss_gib": report["peak_rss_gib"],
    }, indent=1))
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
