#!/usr/bin/env python3
"""Build the frozen v9 + PerturbMatch + post-whitening CD4 submission.

This is a production adapter, not another model-selection harness.  It keeps the
v9 challenge-generation contract fixed:

* K562-GWPS, Orion HCT116, and Orion HEK293T have equal source weights;
* alpha=1.0796, beta=0, gamma=0.65, union-axis preserve whitening;
* CPM>5 gating and exact-total quasi emission at spread=0.25;
* independent RNG streams for control-library selection and count emission.

Two frozen changes are applied to the common 300 x 18,533 effect matrix:

1. Replace the K562 source slot 50:50 with PerturbMatch only when both
   analyses contain the target.  A target present in only one of the two is
   preserved exactly.  Genes absent from PerturbMatch retain the incumbent
   K562-slot value.
2. After the existing three-donor production whitener, add the complete
   Rest/Stim8hr/Stim48hr CD4 mean residual using experiment 062's gene
   normalization at fraction 0.05.  Targets missing any condition receive no
   CD4 residual.  Every VCC panel coordinate is excluded from that residual.

The script intentionally imports the audited experiment-061 and experiment-062
implementations for PerturbMatch preparation and CD4 preparation rather than
duplicating those transformations.

This file does not run on import.  A full build is intentionally expensive and
should be launched under ``scripts/run_with_memory_guard.py``.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import resource
import sys
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Iterable, Mapping, Sequence

os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")

import anndata as ad
import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
EXP061_PATH = ROOT / "experiments/061-perturbmatch-hct-count/build_candidate.py"
EXP062_PATH = ROOT / "experiments/062-cd4-post-whiten/run.py"

FROZEN_IMPLEMENTATION_SHA256 = {
    "scripts/build_submission.py": (
        "3e0a7643f604c59b5ff8494af2fa36dc08c37303dae19bc7557ac95a8dde4e44"
    ),
    "scripts/make_submission.py": (
        "e85879effd33334cd6b1ddd189b02e1373ea6f918add30298e9cf9bebb6a0ab1"
    ),
    # Updated 2026-08-27 for the ADDITION of emit_counts(..., overdispersion=). The
    # parameter defaults to 0.0, accepts a scalar or one value per gene, and the Poisson
    # branch is preserved verbatim, so every
    # existing call path is bit-identical -- which the v11-reproduction test asserts.
    # Previous hash cecb27035a6c301f96e7fe364032216b1211bb2a92fddd1f9439d9cf8730cb3e.
    #
    # Updated 2026-08-31 for the ADDITION of emit_counts(mode="resample_mix") plus the
    # helpers _quasi_block / _repair_totals / _stochastic_round. The new mode is reachable
    # only when control_counts, control_profile and param_mask are all supplied, which the
    # builder does only when VCC_RESAMPLE_X is set; unset, the mode is unreachable. The
    # quasi, mixed, even and nb branches are untouched. VERIFIED bit-identical against the
    # pre-change module at phi=0 and phi=4 over 400 cells x 500 genes with a fixed seed --
    # 161,937 and 152,647 nonzeros matching exactly, values, indices and indptr.
    # Previous hash c91cb05e1204138900058be2c82a227f132e373493a3b3e15914cd18f807f3b9.
    "src/vcc2026/donor_transfer.py": (
        "9ffca7ecac591954fdc5f4f4fe18646a4e67be4f1edf7d57a73b79f4865c7dfc"
    ),
    # Updated 2026-08-27 for the ADDITION of `MultiDonorModel.donor_agreement`, which is
    # a new read-only method. Nothing on the v11 build path changed: `effect()`,
    # `basal_weights()`, `cpm_gate_mask()` and the constructor are untouched, and
    # `test_frozen_constants_match_reused_helpers` still passes. Previous hash
    # 3bb75320ed54e438b0a1ea10acff413e9738fc3294056db06ac337ed4d51cd8f.
    "src/vcc2026/multi_donor.py": (
        "1daf35d97b50fb234a37255cf4136d34580f98f9d54faf175dca48461a33129a"
    ),
    "src/vcc2026/whiten.py": (
        "5dbddefde17192cda4bc7f8f7eceb202faf672de1833c42251366078e7e18bc1"
    ),
    "experiments/061-perturbmatch-hct-count/build_candidate.py": (
        "77e6fbf7edf1aeea6d39be4fe801591d0b0c2867f95042b697f0fc4042a2fae4"
    ),
    "experiments/062-cd4-post-whiten/run.py": (
        "16417a9b985d9aad8993bf2f24ce550ea1e6871d9401bb4672d11ec6d0046176"
    ),
}

# The BUNDLE half of the frozen-input table -- `manifest.json`, `gene_names.csv`,
# `pert_counts.csv` and one `context_<id>.h5ad` per context -- lives in
# `scripts/frozen_inputs.json`, NOT here. It is the only part of the table that moves when a
# new panel lands (the final round ships contexts D/E/F and a different 300-target panel),
# and everything below is ours and never moves. Keeping the panel half in a file makes a
# re-pin a reviewable diff (`--repin-frozen-inputs`) rather than a source edit, and keeping
# it under `scripts/` rather than `data/challenge/` keeps it in version control -- a freeze
# whose recorded value is gitignored is not a freeze.
FROZEN_BUNDLE_PATH = ROOT / "scripts/frozen_inputs.json"

FROZEN_INPUT_SHA256 = {
    "donors/k562_gwps/delta.npy": (
        "0c9422280fed7afb81ca474e8b9024af788c9ffbde32e11f96cb52422e296aa9"
    ),
    "donors/k562_gwps/b_ctrl.npy": (
        "d06cf51d51cde1ae999e86bbd3be0bdc15471392bae153a1aa1e21dc5c81b957"
    ),
    "donors/k562_gwps/genes.txt": (
        "f499ffc3077e4f06c86cd4504bfeefd34741bc2fd515be93f7e436c5c08fa2ba"
    ),
    "donors/k562_gwps/perts.txt": (
        "b784821e4a662aa8ab8d55da357fc712f109dda5a464c3e7e047b15b61408bec"
    ),
    "donors/k562_gwps/meta.json": (
        "0bd3478aa6b63164bd4816cd8bbb439b805b6b93b2a267ddf8114b196c5e19aa"
    ),
    "donors/orion_hct116/delta.npy": (
        "72c7a01dc2869bd3e94c37e89b158c3aaf4048a0dff63da6542afa2de4ec2b33"
    ),
    "donors/orion_hct116/b_ctrl.npy": (
        "b3044d439cfac932024e0a48152196387714e4b460a339fb9a03e32fcb2fe069"
    ),
    "donors/orion_hct116/genes.txt": (
        "db8ba814bd30dcc7fe48be659f4dee3f6955b626d7d7e385f4630232733f893e"
    ),
    "donors/orion_hct116/perts.txt": (
        "fd2fbbbebdd43c360cd274a246c1eeeb9a2139711ee12a475a8e116021c5caa8"
    ),
    "donors/orion_hct116/meta.json": (
        "7a1ec5a9323c6ebd4b586aa300b28bcd9ddf271e67f04bdaf353cb15b8b2d6c8"
    ),
    "donors/orion_hek293t/delta.npy": (
        "d95202d984d9322dc3830f23b6ef55001ba8562dea1d72c43fb0c9d5064c7d1b"
    ),
    "donors/orion_hek293t/b_ctrl.npy": (
        "e04234cffbe102731cc562e8fb06bddb8ef0644b2af894b8a0ee6ce079ca4984"
    ),
    "donors/orion_hek293t/genes.txt": (
        "db8ba814bd30dcc7fe48be659f4dee3f6955b626d7d7e385f4630232733f893e"
    ),
    "donors/orion_hek293t/perts.txt": (
        "3846edf57cc7a0bd91a499ce9be4ed68c9c8ccc1f36ba062c3bb1301695ef547"
    ),
    "donors/orion_hek293t/meta.json": (
        "a7310e69e2e09ad42d88f9192a938b14bf3db44be92575123b86508f0818710f"
    ),
    "perturbmatch/delta.npy": (
        "efa20944f6f6485a7cc7f484c7ab70bd2ba491a6e07b231f306e104b2e1e4cf2"
    ),
    "perturbmatch/genes.txt": (
        "e889994e2d0462d064018408ca456fc8be428e660523fb1bd9a7bcd7b3ab5959"
    ),
    "perturbmatch/perts.txt": (
        "1131758c3ab0d0c781a17580894511986d9baeadc1ad1304251f5a133086998e"
    ),
    "perturbmatch/meta.json": (
        "0b8224c2e04b1d2e4a331ecb26c5d1733eb32c143f949424b1cd518299066a41"
    ),
    "cd4/cd4_rest/delta.npy": (
        "b2bafd6e24c2aa34b03ec17836db34fcd9512f78920aed207413a9d727f9e765"
    ),
    "cd4/cd4_rest/b_ctrl.npy": (
        "5c47ef6b15a00a867b7b24034273bf3aba118d7cdf87ce58ec58427751313d81"
    ),
    "cd4/cd4_rest/genes.txt": (
        "cbada1fd0eaeeb0385fa9c905c08ea90f1d472f5f10e832c47836ed9d7fe45fd"
    ),
    "cd4/cd4_rest/perts.txt": (
        "ddcd8ff07e2f38d93fc035ce2aadebe4bcaeb097ef95e1f26e8e4860b641374b"
    ),
    "cd4/cd4_rest/meta.json": (
        "cece9564be470b133ac6399012af42bcfa1934776c1820838c6aa8ab6a15b9cc"
    ),
    "cd4/cd4_stim8hr/delta.npy": (
        "eea7996d6694450fb1167b9ebbb6b709a3477e6d0347b49a4475e13ff7800c99"
    ),
    "cd4/cd4_stim8hr/b_ctrl.npy": (
        "75e35371ec33677462acdbc56938c2cc266d583dc2e75ed50b6b825c4a20161e"
    ),
    "cd4/cd4_stim8hr/genes.txt": (
        "cbada1fd0eaeeb0385fa9c905c08ea90f1d472f5f10e832c47836ed9d7fe45fd"
    ),
    "cd4/cd4_stim8hr/perts.txt": (
        "639feb3fba1d6cf3b5bb5e2ef9dfdc0ff3044adb4c7477b4588c82ff096d828e"
    ),
    "cd4/cd4_stim8hr/meta.json": (
        "9ca1e9b8ccbf785d88b2a31bc47818479327c205145faa78d9b8eaa215bff4cc"
    ),
    "cd4/cd4_stim48hr/delta.npy": (
        "568119ac25c8d16116d8953d3ee61762fe3c45d9ee8c22a6e7b2669c39837a9a"
    ),
    "cd4/cd4_stim48hr/b_ctrl.npy": (
        "a35f6c639937332cc483c82b9717dec9ff8a3fc70e6030e1b9f6df5fcc6658fd"
    ),
    "cd4/cd4_stim48hr/genes.txt": (
        "cbada1fd0eaeeb0385fa9c905c08ea90f1d472f5f10e832c47836ed9d7fe45fd"
    ),
    "cd4/cd4_stim48hr/perts.txt": (
        "4bf6468927fc11f0709bb0e63d14852f5a49f4f577455f0fa92e88719777d618"
    ),
    "cd4/cd4_stim48hr/meta.json": (
        "62390f52782a592b2a144349bebebffc2675309a614bff0eef99c0b315668ef8"
    ),
    "whitener/whitener_union.npz": (
        "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3"
    ),
    "whitener/whitener_union.json": (
        "a77f9425ac524cdcf3a53a1744520e1e63cb8b925ef77fae2864c4d8b72337b6"
    ),
}

# Registered whitening operators, by `--whitener-id`. The frozen-input table above pins
# ONE operator, which is right for a single production artefact and wrong the moment a
# second operator has to be compared against it -- experiment 091's pruned basis could
# not be built at all without either weakening the guard or registering it here.
#
# Registering keeps the guard fail-closed in the way that matters: an operator must be
# named AND match both of its recorded hashes. An unlisted id, or a listed id whose file
# has drifted by one byte, still aborts before any model input is read.
WHITENER_REGISTRY: dict[str, tuple[str, str]] = {
    # v11 production: the full union basis, 33,787 rows.
    "union": (
        "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3",
        "a77f9425ac524cdcf3a53a1744520e1e63cb8b925ef77fae2864c4d8b72337b6",
    ),
    # Experiment 091: the same basis with the 4,393 common-essential rows dropped.
    # +0.002409 avg6 (t = 19.85) against a SIZE-MATCHED random prune, so the effect is
    # composition and not basis size. Its metadata records `dropped_in_panel: []`, and the
    # common-essential list has zero overlap with the 300-gene challenge panel, so the
    # prune cannot remove signal the submission needs.
    "union_pruned": (
        "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
        "838f263c6a4e1bed069d5faeb3bc3bb2b9234c5a5a933219a94241037a801f4d",
    ),
}

# --- THE OCT-22 GUARD: the operator's AXIS, keyed by the operator's own sha256 (exp 235).
#
# `DEF-RUNBOOK` section 0.1: the whitening operator's axis is
# `sorted_by_official_index(union(donor coverage) - panel)`, so it is A FUNCTION OF THE PANEL,
# and pointing the A/B/C operator at a new 300-target panel fails SILENTLY -- 98.8% of the
# coordinates pair with the wrong gene while every check passes. 212 priced that at **-0.0429
# overall on `pds` alone**, on the only entry that counts, with the board dark.
#
# The reason every check passed is that there was only one, `n_cols`, and 235 measured it at
# **17,942 on all four panels on disk** (A/B/C and 199's three synthetic D/E/F bundles). It
# cannot detect the failure it exists to detect. `columns_sha256` in the meta is closer but is
# only RECORDED, never asserted, and it hashes the integer column indices -- 235 measured
# 199's `def-bundle-badaxis` (two symbols transposed in `gene_names.csv`) as having the same
# index digest as `def-bundle` and a different SYMBOL digest.
#
# So the assertion is on the ordered column SYMBOLS. Operators fitted from 2026-09-04 carry
# their own `axis_sha256` in their metadata json (written by `build_union_whitener.py`, and
# that json's sha256 is itself pinned in `WHITENER_REGISTRY`), which is where the Oct-22 refit
# will get it. The two operators fitted BEFORE that predate the field, so their axis digests
# are pinned here -- by OPERATOR SHA256, never by id or filename (rule 11), so that every
# caller of `load_production_operator` is covered whether or not it passes an id.
#
# Both entries share one axis because the prune drops basis ROWS, not columns.
LEGACY_OPERATOR_AXIS_SHA256: dict[str, str] = {
    # data/donors/whitener_union.npz        ("union", 33,787 basis rows)
    "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3": (
        "edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f"
    ),
    # data/donors/whitener_union_pruned.npz ("union_pruned", 29,394 rows; production)
    "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123": (
        "edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f"
    ),
}


HEX64 = "0123456789abcdef"


def load_frozen_bundle(path: Path = FROZEN_BUNDLE_PATH) -> dict[str, object]:
    """Read and validate the bundle half of the frozen-input table.

    Fail-closed on the shape as well as the values: a table that has drifted out of sync
    with its own `contexts` list (an extra `context_D.h5ad` entry, a missing one, a
    truncated digest) must abort here rather than turn into a confusing path-set mismatch
    later.
    """
    table = json.loads(path.read_text())
    contexts = table.get("contexts")
    if (
        not isinstance(contexts, list)
        or not contexts
        or not all(isinstance(value, str) and value for value in contexts)
        or len(set(contexts)) != len(contexts)
    ):
        raise ValueError(f"{path}: 'contexts' must be a non-empty list of unique strings")
    files = table.get("files")
    if not isinstance(files, dict):
        raise ValueError(f"{path}: 'files' must be an object")
    expected_keys = {f"bundle/{name}" for name in bundle_file_names(contexts)}
    if set(files) != expected_keys:
        raise ValueError(
            f"{path}: 'files' does not match 'contexts': "
            f"missing={sorted(expected_keys - set(files))}, "
            f"unexpected={sorted(set(files) - expected_keys)}"
        )
    for name, digest in files.items():
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or any(ch not in HEX64 for ch in digest)
        ):
            raise ValueError(f"{path}: {name} is not a 64-character sha256")
    for key in ("n_constructs", "cells_per_pert"):
        if not isinstance(table.get(key), int) or table[key] <= 0:
            raise ValueError(f"{path}: '{key}' must be a positive integer")
    support = table.get("panel_support", None)
    if support is not None:
        if set(support) != {"perturbmatch_replacement_targets", "cd4_complete_targets"}:
            raise ValueError(
                f"{path}: 'panel_support' must hold exactly "
                "perturbmatch_replacement_targets and cd4_complete_targets, or be null"
            )
        if not all(isinstance(value, int) for value in support.values()):
            raise ValueError(f"{path}: 'panel_support' values must be integers")
    return table


def bundle_file_names(contexts: Sequence[str]) -> list[str]:
    """The bundle files a round ships, for any context set (not just A/B/C)."""
    return [
        "manifest.json",
        "gene_names.csv",
        "pert_counts.csv",
        *[f"context_{context}.h5ad" for context in contexts],
    ]


def read_bundle_contexts(bundle: Path) -> tuple[str, ...]:
    """The context ids `manifest.json` declares, read BEFORE anything is hashed.

    The manifest has to be read to know which `context_*.h5ad` files exist at all, so it is
    read here untrusted and its own hash is verified in the same pass as the rest of the
    bundle. A tampered manifest that names a different context set is caught by the
    cross-check against the frozen table, which happens before a single byte is hashed.
    """
    manifest = json.loads((bundle / "manifest.json").read_text())
    contexts = manifest.get("contexts")
    if not isinstance(contexts, list) or not contexts:
        raise ValueError(f"{bundle}/manifest.json declares no contexts")
    return tuple(map(str, contexts))


def whitener_expected_hashes(
    whitener_id: str, frozen_bundle: Mapping[str, object] | None = None
) -> dict[str, str]:
    """The frozen-input table with the whitener pair swapped for the named operator."""
    if whitener_id not in WHITENER_REGISTRY:
        raise ValueError(
            f"unknown --whitener-id {whitener_id!r}; "
            f"registered: {sorted(WHITENER_REGISTRY)}"
        )
    npz_sha, json_sha = WHITENER_REGISTRY[whitener_id]
    if frozen_bundle is None:
        frozen_bundle = load_frozen_bundle()
    table = {**dict(frozen_bundle["files"]), **FROZEN_INPUT_SHA256}
    table["whitener/whitener_union.npz"] = npz_sha
    table["whitener/whitener_union.json"] = json_sha
    return table


def file_sha256(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def verify_expected_hashes(
    paths: Mapping[str, Path],
    expected: Mapping[str, str],
    *,
    label: str,
) -> dict[str, str]:
    """Hash every named path and reject omissions, additions, and drift."""
    if set(paths) != set(expected):
        raise ValueError(
            f"{label} path set mismatch: "
            f"missing={sorted(set(expected) - set(paths))}, "
            f"unexpected={sorted(set(paths) - set(expected))}"
        )
    observed: dict[str, str] = {}
    for name in sorted(paths):
        path = paths[name]
        if not path.is_file():
            raise FileNotFoundError(f"{label} input is missing: {name} -> {path}")
        digest = file_sha256(path)
        observed[name] = digest
        if digest != expected[name]:
            raise ValueError(
                f"{label} hash drift for {name}: {digest} != {expected[name]}"
            )
    return observed


def verify_frozen_implementation_hashes() -> dict[str, str]:
    paths = {name: ROOT / name for name in FROZEN_IMPLEMENTATION_SHA256}
    return verify_expected_hashes(
        paths,
        FROZEN_IMPLEMENTATION_SHA256,
        label="frozen implementation",
    )


# Local helpers are executable model inputs. Verify them before importing any of
# them, so a drifted helper cannot run and merely report its hash afterward.
verify_frozen_implementation_hashes()

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from build_submission import (  # noqa: E402
    CONTEXT_COL,
    MAX_STORED_ENTRIES,
    N_GENES,
    PERT_COL,
    ChallengeSpec,
    _str_dtype,
    _write_categorical,
    _write_str,
)
from vcc2026 import rc_emit  # noqa: E402
from vcc2026.donor_transfer import Donor, TS, emit_counts  # noqa: E402
from vcc2026.multi_donor import MultiDonorModel  # noqa: E402
from vcc2026.whiten import Whitener  # noqa: E402


def _load_module(name: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


EXP061 = _load_module("production_pm_cd4_exp061", EXP061_PATH)
EXP062 = _load_module("production_pm_cd4_exp062", EXP062_PATH)

K562 = "k562_gwps"
PRODUCTION_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")

# `VCC_ALPHA` takes a number (the shipped behaviour, default 1.0796) or the literal string
# `knee`, which solves ONE alpha PER CONTEXT at that context's `mse` optimum.
#
# Experiment 170 established that raw `mse` is minimised where the across-perturbation budget
# `B_c` meets the claim `claim_c = sum_p w_p min(Chat_p, C_p)`: below that point the metric
# refunds the whole differential part of the amplitude through the #348 budget, above it the
# refund saturates at `claim` and the extra amplitude is pure cost. The optimum is therefore
# at `rho_c = min(1, B_c/claim_c) = 1`, i.e. at the amplitude `s* = sqrt(claim_c / B_c)`.
#
# a15/a23's single global 0.982 was fitted against A/B/C BOARD points, and 170 measured that
# it leaves the three contexts at rho = 0.992 / 0.940 / 0.882 -- A on the knee, B and C below
# it. A global alpha can only ever put the panel MEAN somewhere sensible; the per-context
# optimum is a per-context number.
#
# Everything the knee needs is computable with no board and no answers:
#   `B_c`      from our own rendered prediction (096 (a))
#   `Chat_p`   from our own emitted cells
#   `C_p`      from the CONTEXT'S OWN NTC control cells -- 096's delete-1 jackknife of one
#              ntc construct's 400 cells, which ships inside `context_<id>.h5ad`
# so it transfers to a D/E/F bundle on the day it lands, which a board-fitted constant does
# not.
#
# `VCC_ALPHA_REF` is the amplitude the effect matrix is FIRST built at, and the anchor the
# knee is expressed against: `alpha_c = alpha_ref * sqrt(claim_c / B_c(alpha_ref))`. If every
# context were already on its knee this returns `alpha_ref` exactly, which is the sense in
# which the rule reproduces the shipped scale. It is close to anchor-free in practice --
# `B` is very nearly quadratic in the amplitude and `claim` is saturated at the `C_p` cap --
# and `knee180.py` measures that residual rather than assuming it.
ALPHA_ENV = os.environ.get("VCC_ALPHA", "1.0796").strip()
ALPHA_MODE = "knee" if ALPHA_ENV.lower() == "knee" else "fixed"
ALPHA_REF_ENV = os.environ.get("VCC_ALPHA_REF", "").strip()
if ALPHA_MODE == "knee":
    ALPHA = float(ALPHA_REF_ENV) if ALPHA_REF_ENV else 1.0796
else:
    if ALPHA_REF_ENV:
        raise ValueError(
            "VCC_ALPHA_REF is only meaningful with VCC_ALPHA=knee; got "
            f"VCC_ALPHA={ALPHA_ENV!r}"
        )
    ALPHA = float(ALPHA_ENV)
if not np.isfinite(ALPHA) or ALPHA < 0.0:
    raise ValueError(f"VCC_ALPHA must be a non-negative number or 'knee'; got {ALPHA_ENV!r}")
BETA = 0.0
# DE-GATED AND CLOSED, 2026-08-28. Experiment 097 swept this the way the shrinkage was swept
# -- 031-somatic's build path reproduced exactly, the cached operator never refitted, a13's
# shrinkage on top, every arm re-optimised onto its own cliff, six members, re-seeded twins --
# and NO gamma is net positive on both surfaces. `kappa` is monotone DECREASING in gamma and
# `pds` monotone INCREASING, on both, across the whole range: the operator is exactly the trade
# exp017 bought it for and both halves are real. Going 0.65 -> 0, each +0.01 of kappa costs
# 0.0181 raw pds on h1-test and 0.0517 on orion, so pds alone more than eats the directional
# gain before any DE member is counted (net -0.0033 / -0.0080). Do not re-open this on an mse
# measurement; that is the surface that recommended it.
# Overridable the same way VCC_ALPHA is. The whitener is `W = (C + lam I)^(-gamma/2)`
# (`src/vcc2026/whiten.py:79`), and `lam` has sat at its 0.0 default in every build this
# project has ever shipped -- it was threaded through `apply_v9_whitening` and recorded in
# the meta as `whiten_lambda`, but never swept. With lam = 0 the union whitener's smallest
# eigenvalues (~1e-11) are amplified enormously, which is why AGENTS.md 4.2 records that
# whitening makes |predicted effect| a WORSE per-gene confidence signal. A ridge floors
# those near-null directions, which decouples the two jobs gamma is doing at once:
# suppressing dominant shared programmes (which buys pds) and amplifying the tail (which
# costs per-gene sign). Defaults are the shipped values, so omitting both variables
# reproduces every prior build bit-for-bit.
WHITEN_GAMMA = float(os.environ.get("VCC_WHITEN_GAMMA", "0.65"))
WHITEN_LAM = float(os.environ.get("VCC_WHITEN_LAM", "0.0"))
# The CD4 helper's frozen whitening gamma. NOT overridable: sweeping the main whitener
# must leave the CD4 residual path untouched.
CD4_WHITEN_GAMMA = 0.65
# Fraction of genes (by |predicted effect|) kept on the parametric quasi law when the
# resampling emitter is active. Unset -> the mode is OFF and every prior build reproduces
# bit-for-bit. See `emit_counts(mode="resample_mix")`.
# Invert the mask: resample the TOP genes by |predicted effect| and keep the tail
# parametric. a21 did the opposite and lost -0.0642 on `nmae`, because the cell-channel
# payoff scales as (2*acc - 1) and it was applied to the tail at acc ~0.51 (multiplier
# 0.02). Top-decile acc is 0.7359 (multiplier 0.47), ~20x better. Coverage is settled
# saturated by a21, so silencing head calls no longer costs `fid`.
RESAMPLE_INVERT = os.environ.get("VCC_RESAMPLE_INVERT", "").strip() not in ("", "0", "false")
RESAMPLE_X = os.environ.get("VCC_RESAMPLE_X", "").strip()
RESAMPLE_X = None if not RESAMPLE_X else float(RESAMPLE_X)
if RESAMPLE_X is not None and not (0.0 <= RESAMPLE_X <= 1.0):
    raise ValueError(f"VCC_RESAMPLE_X must lie in [0, 1]; got {RESAMPLE_X}")
if not (0.0 <= WHITEN_GAMMA <= 3.0):
    raise ValueError(f"VCC_WHITEN_GAMMA must lie in [0, 3]; got {WHITEN_GAMMA}")
if not (0.0 <= WHITEN_LAM <= 1.0):
    raise ValueError(f"VCC_WHITEN_LAM must lie in [0, 1]; got {WHITEN_LAM}")
# `VCC_DONOR_WEIGHTS` — the per-donor mixing weights, in `PRODUCTION_DONORS` order
# (experiment 175, priced composed in 178).
#
# `donor_weights()` has hard-coded `{name: 1.0}` since v9 because experiment 019 measured
# equal weighting as the BEST configuration for `pds` under whitening -- it flipped 010's
# result, and the mechanism is that averaging's smoothing cost is precisely what whitening
# removes. Experiment 175 then measured the OTHER side of the same trade on the h1-test DE
# pool: per-donor sign accuracy is K562 0.5670 / HCT116 0.5375 / HEK293T 0.5308, the Orion
# pair carries two-thirds of the weight while adding only +0.0027 over K562 alone, and a
# K562-heavy (2,1,1) mix buys +0.004 pooled sign accuracy. Same profile, two objectives,
# opposite optima -- so this is a rider that costs `pds` and buys every other member, and it
# has to be priced composed rather than alone.
#
# Format: as many comma-separated non-negative floats as `PRODUCTION_DONORS`, e.g.
# `VCC_DONOR_WEIGHTS=2,1,1`. `MultiDonorModel.effect` renormalises over the donors that
# actually carry a perturbation, so the values are relative, not absolute. Unset reproduces
# every prior build bit-for-bit (the default is the literal `{name: 1.0}` dict, not a
# reconstruction of it).
DONOR_WEIGHTS_ENV = os.environ.get("VCC_DONOR_WEIGHTS", "").strip()
DONOR_WEIGHTS: dict[str, float] | None = None
if DONOR_WEIGHTS_ENV:
    _parts = [chunk.strip() for chunk in DONOR_WEIGHTS_ENV.split(",")]
    if len(_parts) != len(PRODUCTION_DONORS):
        raise ValueError(
            f"VCC_DONOR_WEIGHTS needs {len(PRODUCTION_DONORS)} comma-separated values in "
            f"the order {PRODUCTION_DONORS}; got {DONOR_WEIGHTS_ENV!r}"
        )
    _values = [float(chunk) for chunk in _parts]
    if any(value < 0.0 for value in _values) or not any(value > 0.0 for value in _values):
        raise ValueError(
            f"VCC_DONOR_WEIGHTS must be non-negative with at least one positive value; "
            f"got {DONOR_WEIGHTS_ENV!r}"
        )
    DONOR_WEIGHTS = dict(zip(PRODUCTION_DONORS, _values))
# `VCC_SIGN_REPAIR` — undo the sign flips whitening introduces (experiment 173).
#
# 172 measured that whitening costs 0.011 pooled sign accuracy on the reference's DE pool
# while buying `pds`. 173 decomposed it: on the 15.7% of pool entries where the whitened
# effect disagrees in sign with the un-whitened donor consensus, the whitened sign is right
# **45.6%** of the time -- worse than chance -- and the consensus is right 54.4%. Those
# entries concentrate where the consensus is small (35% of its bottom quartile, 1.4% of its
# top decile) and carry only 1.6% of the whitened effect's energy.
#
# `cons` sets the sign of those entries back to the consensus's and keeps the magnitude, so
# it is a surgical partial un-whitening: it lands on the same `pds`-vs-accuracy exchange
# curve as `WHITEN_GAMMA` (-0.054 raw `pds` per +100% transfer cosine against gamma's
# -0.068), but it is a MECHANISM rather than a fitted exponent, so it should carry to an
# unseen line. 173's random-flip null (same count of entries per row) collapses `pds` by
# -0.034 and worsens `nmae` by +0.011, which is the evidence the targeting is real.
#
# The repair is applied to the whitened matrix BEFORE the CD4 residual is added, against the
# pre-whitening effect the operator was handed -- i.e. against the same consensus whose signs
# the operator moved, and nothing else. Unset is a no-op and reproduces every prior build
# bit-for-bit.
SIGN_REPAIR = os.environ.get("VCC_SIGN_REPAIR", "").strip().lower() or None
if SIGN_REPAIR is not None and SIGN_REPAIR != "cons":
    raise ValueError(f"VCC_SIGN_REPAIR must be 'cons'; got {SIGN_REPAIR!r}")
PM_BLEND = 0.50
CD4_MODE = "gene_norm"
# Overridable the same way VCC_ALPHA is, and for the same reason: it is a build-time scalar
# whose optimum is measured per arm rather than fixed. The default is the shipped 0.05, so
# omitting the variable reproduces every prior build bit-for-bit.
#
# Experiment 099 measured the post-whitening CD4 residual at x1.077 transferable on directional
# accuracy with Delta pds +0.0065 / +0.0109 -- the only added-signal route in that study whose
# readable members were all non-negative on both surfaces. NOTE its own control: adding back the
# base mixture's OWN unwhitened effect reproduces most of the direction curve (x1.122 / x1.123
# against CD4's x1.134 / x1.077), so most of this is the whitening-relaxation lever experiment
# 097 closed. What is CD4-specific is the pds sign, and that is +0.0018 transferable. The study
# pre-registered this arm AND stated in writing that it expects it to under-deliver.
CD4_FRACTION = float(os.environ.get("VCC_CD4_FRACTION", "0.05"))
# `CD4_FRACTION` is NOT experiment 099's `f`, and conflating them is the specific error 099
# §5.4 wrote a warning about: 062 rescales each condition to EFFECT_SCALE before averaging and
# never rescales the mean back up, so the fraction multiplies a differently-normalised object
# and production's effective `f` lands somewhere in 0.033-0.049 -- a range, because the
# challenge-panel mixture norm was never measured. Two codebases with a knob called `fraction`
# is what cost this project a 5x error on `gamma`.
#
# `f` is defined in 099 (`postwhiten.py`) as the added residual's median row norm over the
# whitened base's median row norm, both over the FULL gene axis, medians taken over every row
# including the ones carrying no residual. That is measurable here exactly, and because the
# residual is linear in the fraction, one measurement solves the fraction in closed form.
# Setting `VCC_CD4_TARGET_F` does that solve at build time and ignores `VCC_CD4_FRACTION`;
# leaving it unset changes nothing and every prior build reproduces bit-for-bit.
CD4_TARGET_F = os.environ.get("VCC_CD4_TARGET_F", "").strip()
CD4_TARGET_F = float(CD4_TARGET_F) if CD4_TARGET_F else None
if CD4_TARGET_F is not None and not 0.0 < CD4_TARGET_F < 1.0:
    raise ValueError(f"VCC_CD4_TARGET_F must be in (0, 1); got {CD4_TARGET_F}")

# --- `VCC_H1_UNION_G` — the COVERAGE-AWARE sizing rule, and the two H1 sources it unlocks.
#
# THE BUG (experiment 218 section 7.1, and the same failure class 200 and 207 found).
# `solve_cd4_fraction` above sizes a residual by the ratio of its MEDIAN ROW NORM OVER EVERY
# PANEL ROW to the whitened base's.  A source that covers <= 50% of the 300-target panel has a
# median row norm of EXACTLY ZERO, so `f_per_unit` is 0 and the solved fraction is 0.  The
# residual is not small, it is identically absent -- with no error, no warning, and every
# downstream check passing, because production's only guard fires on a FULLY zero residual and
# a partial-coverage source is not fully zero.  Two sources this repo has held for weeks are
# therefore absent from every build ever made: Pan-2026 H1 (64/300) and `h1_2025` (20/300).
#
# THE FIX.  Size a partial-coverage residual on the rows it ACTUALLY CARRIES --
# `median(||added_i|| : ||added_i|| > 0)` -- so `g` means "on the rows this source covers, the
# added norm is `g` x the base's typical row norm".  For a source above 50% coverage the two
# definitions coincide; CD4 sits at 93-98% and ITS OWN SOLVE IS NEVER TOUCHED by this flag
# (216/217's territory).  `solve_cd4_fraction` is not modified.
#
# THE ARM.  Unset (the default) this is inert and every prior build reproduces bit-for-bit.
# Set to a dose `g`, the builder adds the UNION of the two H1 sources -- 77/300 of the panel,
# the two sharing only 7 targets -- as a second un-whitened post-whitening residual, solved on
# its carrying rows, merged with the CD4 residual into ONE add followed by production's single
# global preserve rescale.  That is 218 section 7.2's construction exactly: `compose218.py`'s
# `merged()` + `add_post_whiten_residual(..., fraction=1.0)`.  218 measured macro cosine
# +7.680% and `pds` +0.004887 at `g = 0.28` across the three cancer-line PLOLO folds, positive
# on cosine AND `pds` on 12 of 12 (fold x dose) cells.
H1_UNION_G = os.environ.get("VCC_H1_UNION_G", "").strip()
H1_UNION_G = float(H1_UNION_G) if H1_UNION_G else None
if H1_UNION_G is not None and not 0.0 < H1_UNION_G < 5.0:
    raise ValueError(f"VCC_H1_UNION_G must be in (0, 5); got {H1_UNION_G}")
# --- THE SOURCE REGISTRY (experiment 229).  Rule 11: a path AND a sha, never a description.
#
# Everything `VCC_H1_UNION_G` does is already source-agnostic -- `prepare_h1_union` takes a
# list of `(name, path)` banks in the standard compact-donor layout, unions them, sizes the
# result on its CARRYING rows and merges it with the CD4 residual into one add.  The only
# H1-specific thing was this list, which used to be a hardcoded two-entry tuple.  It is now a
# name -> (path, per-file sha256) table plus `VCC_UNION_SOURCES`, an ordered comma-separated
# list of which entries are active.
#
#   VCC_UNION_SOURCES unset  ->  the H1 pair, so `VCC_H1_UNION_G=0.28` alone still reproduces
#                                the shipped `h1union` package BIT-FOR-BIT (known-answer row 2).
#   VCC_UNION_SOURCES=a,b,c  ->  those banks, IN THAT ORDER.
#
# ORDER IS PART OF THE ARM.  `union_residual` accumulates a float64 sum per coordinate, so a
# permutation of the same members can differ in the last bit; the env var's order is the
# order the union is built in and it is recorded in the meta.
#
# The shas are VERIFIED at build time (`resolve_union_sources`), not merely recorded: three of
# these banks live under `experiments/*/work/`, which other experiments write to, and a build
# that silently picked up a re-derived table would be a different arm wearing this one's name.
UNION_SOURCE_REGISTRY: dict[str, dict[str, object]] = {
    # --- 222's shipped pair (do not reorder: `VCC_UNION_SOURCES` unset selects these two) ---
    "gse295214_h1": {
        "path": ROOT / "experiments/115-gse295214-h1/work/compact-donor",
        "sha256": {
            "delta.npy": "d0fa6555f9a62ce503d2ad55c8bfda954322657b52dcd48a24e29c330e718c7a",
            "b_ctrl.npy": "c88816a57c57c3a088183ae24bca406d8ded374871e3d330aa444fcb13811b4f",
            "genes.txt": "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723",
            "perts.txt": "7060076272d411d4ee00b898cc26eedcb1eaaac5f5e94e2c60f3b45978320747",
        },
    },
    "h1_2025": {
        "path": ROOT / "data/donors/h1_2025",
        "sha256": {
            "delta.npy": "269628f6eae749857e619bad4556ce9ac8b8810882e7af8a9925e3178eec0cc8",
            "b_ctrl.npy": "bd527b8373d56dab8bfa4d6a1a54560987977e37e22630709fcbbd8f7ce42d7a",
            "genes.txt": "da2b383939d0939e911d273539b136c7a10f0482b2ee509d294e1632d5629822",
            "perts.txt": "977b808e723c73eea5dc4a60cec6e71512387dc520512bf133d3b18fce7287dd",
        },
    },
    # --- the sources experiment 225 swept and priced (experiment 223 built the first) -------
    "dld1_lowmoi": {
        "path": ROOT / "experiments/223-silent-zero-audit/work/tables/dld1_lowmoi",
        "sha256": {
            "delta.npy": "fd7bd0637cf9949971040f6c184b28f288391191f2e956c7cb16bb92d266aaee",
            "b_ctrl.npy": "6e59a43ab777d8c42a5ad2cf2db69af31b2909e65dfa1c7844339a980ca6d0d9",
            "genes.txt": "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723",
            "perts.txt": "db4847da3fd2312d66bc7e4808dcea597ed56af3fb2dee4bb91d8f25a801867d",
        },
    },
    "jiang_5line_noK562": {
        "path": ROOT / "experiments/225-source-sweep/work/tables/jiang_5line_noK562",
        "sha256": {
            "delta.npy": "11c4e9ec149419ec6e49d3ea8334900da51b9ceb370042fcb9fdd286b52eed97",
            "b_ctrl.npy": "6cd2fc18ea809d4a1611217b6976dc80cdf750330225789001398f6955f95df6",
            "genes.txt": "df49bf9411353f762a93b01c38535b68858a5d11b2931705dfddbebed29c5643",
            "perts.txt": "254d8388adad410ad40ab20923aefb738751d703165a2580f7c125980e63363e",
        },
    },
    "igvf_cardiomyocyte": {
        "path": ROOT / "experiments/225-source-sweep/work/tables/igvf_cardiomyocyte",
        "sha256": {
            "delta.npy": "d30b4930671c5f4b18c11ac88bb40b3c27a4d9c0a03875f48ba4ed2258736935",
            "b_ctrl.npy": "9eaf6088d2c264f0d7d02db719082456daeaf06182cc2fd4b6654e840e531da3",
            "genes.txt": "696c1f76d7ef6d69fa591e8ad0e940055bc77eacbce6a2c0beef42728048d382",
            "perts.txt": "3f272e2b2c0d635f2bdaf1c143f923909fc23cf588f23728d2f2f42510845cff",
        },
    },
    "myllia_a549": {
        "path": ROOT / "experiments/120-myllia-a549/work/compact-donor",
        "sha256": {
            "delta.npy": "e81436089579914637bfcb54df454e5bf4b3759fba44b8d8e48e5204ed1c9ed6",
            "b_ctrl.npy": "6a2c252df8d7b2ae0a1b555f74d15f1678af4ad21e5661b2e2e24973e4f39421",
            "genes.txt": "c2b8dec58c02810d64018164155b655c4e7c7be451a3ae579965db193af33723",
            "perts.txt": "56906d4a9a91925ac0a29194caf88434880d026e041ed0626d448acf9b325f19",
        },
    },
    # --- THE OCTOBER FOUR.  Registered, never selected. ------------------------------------
    # These four carry **0/300** of today's panel, so PLOLO cannot price them at all -- it
    # scores today's panel and they touch none of it.  225 §8.3 projected every source onto
    # `DEF-RUNBOOK` §5's two synthetic draws and found each of them carries ~30/300 of the
    # PESSIMISTIC October panel, which is where the incumbent H1 pair falls to 36.  They are
    # pre-registered here, with paths and shas, so that on 2026-10-22 selecting one is an env
    # var rather than a bank build under time pressure.
    #
    # Registry membership is INERT: nothing reads an entry that `VCC_UNION_SOURCES` does not
    # name, and `verify_union_source_hashes` only hashes the selected banks.  229 §5.2 proves
    # that by rebuilding the shipped package with these four present and re-deriving its hash.
    "nadig_hepg2": {
        "path": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2",
        "sha256": {
            "delta.npy": "ca4565e84045b91dbfca35280c01ae73f9a67112b3de2c8ba799cfda3fbce567",
            "b_ctrl.npy": "d8a986250735a35d1d5bc3011a00eade9c17963d8582718bb10c2f58d3720cb7",
            "genes.txt": "e808cb642deca07480ca9a19a139730c75fc4abfa644b96a32261026951e12a9",
            "perts.txt": "cea71d4d5d860f12cc2c3d77342159c1169dfa3ad9e836680eb66e3f12050866",
        },
    },
    "nadig_jurkat": {
        "path": ROOT / "data/donors/nadig_jurkat",
        "sha256": {
            "delta.npy": "1fa16024c2d7115fc4b7c06f0ab7d9afe8ba6ebd0ae6dcef62c4955636e29129",
            "b_ctrl.npy": "d9e59318b785ba1d32189975298b1af8cfb61ebe772f5c97e4b47f2074f53ac9",
            "genes.txt": "3071f4366db6896b2dac2cb0d59fbfff18f89840f021d556e42ccee083c39830",
            "perts.txt": "f5a5c1bbbbdf97bcedb789277b13b0901447c28a5f091d798c37a93ef86197d7",
        },
    },
    "rpe1": {
        "path": ROOT / "data/donors/rpe1",
        "sha256": {
            "delta.npy": "253cc16e51f4635f23e8e597d91bd38a84162a1102486d4a5764c48c379c2b4a",
            "b_ctrl.npy": "2c736c3f755df0c717428f387e65907dfa4f530816dc7855bd177da8d24a65f3",
            "genes.txt": "ac7c7f7b141a74004a5363b2e1aba52f5769c02ea769d9a5b72fdb5eacdbd2f0",
            "perts.txt": "4377d2dc2fd57c523d6632416ca5b9f9e7c192192e32b0eebd91a84fabab8c67",
        },
    },
    "k562_essential": {
        "path": ROOT / "data/donors/k562_essential",
        "sha256": {
            "delta.npy": "4daf73cb932fdb089857c8629a70e83d8d3090570aa1f7d678c831b620ea281b",
            "b_ctrl.npy": "72e57c224d232adc3e423ac2c0d05f053f3be3314821d1eb8fce6fd973f69b8b",
            "genes.txt": "7f1617ef778eb0a51a2b7998d42bc78ce1c58f6fdf5e08f2b8aa196d90c819d7",
            "perts.txt": "b1e5279534eb240f0f3487f0814d379e05b12799b21f6cb3818f5c7f55557e2d",
        },
    },
}
# Unset selects 222's pair, so the shipped arm is reproducible from `VCC_H1_UNION_G` alone.
UNION_SOURCES_DEFAULT: tuple[str, ...] = ("gse295214_h1", "h1_2025")
UNION_SOURCES_SELECTED = os.environ.get("VCC_UNION_SOURCES", "").strip()
UNION_SOURCE_NAMES: tuple[str, ...] = (
    tuple(part.strip() for part in UNION_SOURCES_SELECTED.split(",") if part.strip())
    if UNION_SOURCES_SELECTED
    else UNION_SOURCES_DEFAULT
)
if len(set(UNION_SOURCE_NAMES)) != len(UNION_SOURCE_NAMES):
    raise ValueError(f"VCC_UNION_SOURCES repeats a source: {UNION_SOURCE_NAMES}")
for _name in UNION_SOURCE_NAMES:
    if _name not in UNION_SOURCE_REGISTRY:
        raise ValueError(
            f"VCC_UNION_SOURCES names {_name!r}, which is not in UNION_SOURCE_REGISTRY "
            f"({sorted(UNION_SOURCE_REGISTRY)}).  Add it there, with its shas, first."
        )
H1_UNION_SOURCES: tuple[tuple[str, Path], ...] = tuple(
    (name, UNION_SOURCE_REGISTRY[name]["path"])
    for name in UNION_SOURCE_NAMES
)
H1_UNION_FILES = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt")
MIN_CPM = 5.0
# `VCC_MASS_FIX` — make the rendered profile need no renormalisation (experiment 156).
#
# `profile_from_effect` renders a log-space effect as `expm1(b_ctrl + gated)` and then divides
# by its own sum. The scorer reads `log1p(TS * p) - b_ctrl`, so that division is the identity
# if and only if `sum_g expm1(b_ctrl + gated) == TS`. The control profile satisfies that
# exactly by construction (it IS a `log1p(TS * p)` object), and so does every reference
# perturbation. Our effect does not: measured on the frozen 104 artifacts the mass ratio is
# 1.0018, and the renderer's implicit correction lands as a shift concentrated on the
# highest-expressed genes that is nearly the SAME direction for every target -- 4.0-4.5% of the
# row norm at row-row cosine 0.303, against an intended row-row cosine of 0.0040. Whitening's
# whole job is to make those rows mutually dissimilar, and the renderer was putting the shared
# blob back.
#
# expm1 is convex and `S(0) = TS`, so `S` has no second root along the ray `{g * effect}`:
# scaling cannot satisfy the constraint and the correction must move off the ray. But it is
# only ONE scalar, so the direction is ours to choose. `dS/dh_g = exp(b_g + h_g)`, so the
# minimum-L2 correction is `eps = -(S - TS) * w / ||w||^2` with `w_g = exp(b_g + h_g)`,
# restricted to the genes the CPM gate already lets move. Iterating that is Newton's method
# with a minimum-norm step.
#
# The mode is deliberately explicit rather than a boolean: experiment 156 measured that HOW the
# constraint is met dominates the result. `const` (a uniform log shift on the gate, satisfying
# the same scalar) LOSES 0.0056/0.0031 pds and 4.7%/18.8% transfer cosine, while `minnorm`
# GAINS 0.00035/0.00278 pds and 14.6%/67.2% cosine on the two held Orion lines. Unset leaves
# the renderer untouched and every prior build reproduces bit-for-bit.
MASS_FIX = os.environ.get("VCC_MASS_FIX", "").strip().lower() or None
if MASS_FIX is not None and MASS_FIX not in {"minnorm", "const"}:
    raise ValueError(f"VCC_MASS_FIX must be 'minnorm' or 'const'; got {MASS_FIX!r}")
DISPERSION = 0.1
SPREAD = 0.25
# The emission LAWS this builder can drive. `quasi` is the shipped one and the default, so
# an unflagged build writes no `emitter` key and reproduces every prior build bit-for-bit.
EMITTERS = ("quasi", "resample_mix", "rc_hybrid", "rc_hist_a", "rc_callset")
# The two rank-controlled laws of experiment 185. They share every fixed choice below and
# differ in EXACTLY one thing, which is `rc_emit.rc_hybrid_block`'s `call_set_first`:
# `rc_hybrid` (185 section 3.4) re-keys the incumbent's |z| multiset inside the incumbent's
# own BH call set and inside its complement, so the call SET is preserved; `rc_hist_a`
# (185 section 3a) re-keys the whole gene axis, so the set is re-drawn by `|H|` at the same
# size. Experiment 201 exists to price that difference on the board, so it is carried by
# the `emitter` key alone and adds NO second config key -- one variable, one difference.
RC_EMITTERS = ("rc_hybrid", "rc_hist_a", "rc_callset")
RC_CALL_SET_FIRST = {"rc_hybrid": True, "rc_hist_a": False, "rc_callset": True}
# `rc_callset` (experiment 202, riding 192 section 5.3) is `rc_hybrid` with ONE variable: the
# block that goes in front is a PRESCRIBED set instead of the incumbent's own --
#   top 400 by |H|  UNION  top 2,200 by control CPM   (n_pred 2,600, vs `rc_hybrid`'s 4,480)
# -- with the ORDER, the totals, the weld, the sign policy, the tail law and the placement
# all `rc_hybrid`'s, unchanged. `control CPM` is the per-gene mean of per-cell CPM over that
# context's own NTC control cells, which is `rc_emit.RcControlBlock.ref_lin`: the same object
# the scorer gates on, measured and not modelled. Both rankings run over the scorer's own
# 5-CPM DE universe, because a gene outside it can never be called.
# Unlike `rc_hist_a`, this arm CANNOT be carried by the `emitter` key alone -- its three
# numbers are the arm -- so they are written into the config and pre-declared to the rule-10
# diff. See `rc_emit.CallSet` for why `z_cap` exists and why it is the only place the
# incumbent's p-value multiset is allowed to move.
RC_CALLSET_HEAD = 400
RC_CALLSET_CPM = 2200
RC_CALLSET_Z_CAP = 1.90
# The rank-controlled emitters' three fixed choices, recorded in the meta so the rule-10
# diff names them:
# the scorer's own DE gate (mean per-cell CPM > 5 on the context's controls, NOT the
# builder's `--min-cpm`), the `abs` sign policy of 185 section 3a, and the ordering key.
RC_SCORER_MIN_CPM = 5.0
RC_SIGN_POLICY = "abs"
RC_KEY = "gated_effect"
# Ties in `|H|` -- every gene the CPM gate zeroed shares the key value 0 -- are broken by a
# deterministic jitter, seeded per block so a rebuild reproduces the same assignment.
RC_JITTER_SEED_BASE = 1000
SEED = 0
EMISSION_SEED_OFFSET = 1_000_003
GIB = 1024**3

# `--shrink-key agreement`. Both numbers are experiment 096's sweep, not a guess:
# `experiments/096-mse-floor/mech_miller.py` tried five alphas on both harnesses (hard
# top-K at 500/2000/8000 and a soft exponential in rank at gamma 0.5/1.0) and the soft
# exponential won on both. Its rank denominator is 2000 over an ~18k gene axis, which is the
# challenge panel's 18,533 to within 3%.
#
# THE DEFAULT IS 0.25, NOT THE mse OPTIMUM, AND IT IS THE ONE VALUE IN THIS FAMILY THAT HAS
# SURVIVED A DE GATE. a13 was built with the literal flag `--shrink-gamma 0.25` on this
# function, and experiment 097 gated THAT ARTIFACT at +0.0215 (h1-test) / +0.0040 (orion) net
# overall against a12, positive on both surfaces across all six members. (An earlier gate
# reported +0.0095 / +0.0025 for a row labelled "a13"; that row was 097's own power-law form at
# gg=0.25 -- a different object with the same name -- and its error was conservative by luck,
# not design.) This is
# also the ONLY knob in the shrinkage family that `pds` constrains, and on orion `pds` falls
# MONOTONICALLY from gamma = 0, so no positive value is free. Raising it is a follow-up
# experiment, not a tuning knob.
# STANDING RULE FOR EVERY NUMBER IN THIS FAMILY, earned four times over in one session:
# WHEN A NUMBER CROSSES A BOUNDARY, TRANSFER THE DEFINITION, NOT THE VALUE; WHERE THE
# DEFINITION CANNOT TRAVEL, TRANSFER THE RATIO WITHIN ONE FRAME AND SAY SO. The four crossings
# and the check that catches each, at increasing distance:
#   two columns in one table   analytic vs realised multipliers  -> grep the superseded digits
#   two codebases              097's gg vs this --shrink-gamma   -> ask for the formula written
#   two pipelines, one number  reconstructed vs shipped kappa    -> check against the artifact
#                                                                   and report the ratio
#   two panels, one policy     harness vs challenge rank of a    -> state which panel a rank
#                              0.0 cosine                           was computed on
#   code differs, artifact     `damp` vs 097's gate: two          -> run the REAL population
#   does not                   policies, zero differing rows         through both branches
#                              at FLOOR 0.30                         before believing a diff
# Every one of those numbers was CORRECT in its own frame and wrong the moment it crossed. No
# amount of extra precision would have helped; precision was not what was missing. Note the
# third: a per-row cosine of 0.99922 reads as a faithful reconstruction and is not one for a
# quantity living in 2% of the vector's norm -- A VALIDATION STATISTIC CHOSEN AT THE WRONG
# SCALE WILL PASS.
AGREEMENT_SHRINK_RANK_SCALE = 2000.0
AGREEMENT_SHRINK_GAMMA = 0.25
# ...and toward WHAT. The prior was that `expr_mse_unbiased_capped_norm`'s scored 0-anchor
# is the context's mean perturbation response assigned identically to every perturbation, NOT
# the control -- the board has both ends, team "Kun Lab"'s explicit zero-effect submission at
# raw 1.0424 against a scored baseline of raw ~0.9879 -- so a gene shrunk toward the control
# would walk toward the worse anchor. MEASURED, that prior is wrong, and the reason is that
# our own grand mean is not the metric's anchor: `mech_shrink_target.py` puts
# cos(our dbar, the reference's mean effect) at +0.016 (h1-test) and -0.066 (orion), i.e. no
# capture and not even a consistent sign, so shrinking toward it buys none of the baseline's
# advantage while inflating (A - B)/A from ~0.012 to ~0.055 -- and A - B is exactly the term
# the master formula drops. Scored through the identity with B recomputed, 'zero' beats
# 'mean' on both harnesses at every gamma and with both agreement keys, 12/12 -- reproduced
# by a second, independent implementation (`mech_shrink_exact.py`) that disagrees with the
# first on LEVELS and still gives the same ordering. The default was specified as 'mean' and
# was flipped to 'zero' against that spec once the measurement came in; 'mean' stays reachable
# only so the comparison can be rerun.
AGREEMENT_SHRINK_TOWARD = "zero"

# `--shrink-pert-key donor-cos`: the PER-PERTURBATION half, from experiment 097. Its key is
# the mean pairwise cosine between the carrying donors' effect vectors for that target, and
# the weight is `z_p = rank01(key)^gamma`. REALISED through experiment 097's full
# re-emit/re-jackknife/re-score, this half is worth 1.246 (h1-test) / 1.176 (orion) alone and
# 1.473 / 1.580 stacked on the per-gene key, against a permutation null at p <= 0.001. (The
# analytic route reads 1.261 / 1.220 and 1.478 / 1.610; quote the realised figures and keep the
# analytic ones as a cross-check only. Never mix the two within one comparison -- they agree to
# 0.3-1.9%, which is exactly the condition under which nobody notices which column a number
# came from.) Neither half should ship alone; see `--shrink-pert-key`.
#
# gamma 1.0. This was 2.0 -- 097's mse-optimal setting, 1.478 (h1-test) / 1.610 (orion) on a
# 33-cell grid through the exact identity -- and 2.0 IS better on mse. It is worse on the ARM.
# 097's DE gate (production emit_counts, real rank-sum, BH 0.05, 40 perturbations, every arm
# with a re-seeded twin) prices all six members:
#
#     arm            n_pred   reach    fid      jac      nmae     net overall vs a12
#     h1  a13 gp=1     3744   0.1442  0.5041  0.0745   0.9678     +0.0095
#     h1  a13 gp=2     3648   0.1101  0.5022  0.0750   0.9768     +0.0068
#     ori a13 gp=1     3516   0.0743  0.5742  0.0143   0.9369     +0.0025
#     ori a13 gp=2     3456   0.0639  0.5783  0.0145   0.9585     -0.0037
#
# gp=2 is NEGATIVE on orion and `nmae` alone accounts for -0.0165 of it. Rule 1 binds on the
# worse surface, so gp=1. THE MECHANISM IS THE ONE NONE OF THIS EXPERIMENT WAS PRICING: `nmae`
# is an absolute fold-change error scored over the REFERENCE's DE genes, which are not the genes
# the donors agree on, so any shrinkage toward zero makes it monotonically worse (h1 0.9630 ->
# 0.9678 -> 0.9768; orion 0.9189 -> 0.9369 -> 0.9585). a13 at gp=1 pays that and the other
# members more than cover it; at gp=2 the nmae cost overtakes them.
PERT_SHRINK_GAMMA = 1.0
# Applied to the RANK, INSIDE the exponent: `z = max(rank01, floor)**gamma`. Unfloored,
# `rank01` puts the weakest target at 0.5/n, so at n=300 the bottom perturbation keeps 0.17 %%
# of its effect and emits essentially the control profile -- not a crash (the emitter raises
# only on an all-zero row and a control profile is dense) but a perturbation with no DE calls
# is an `n_pred` loss the board reads directly, and it leaves pds a near-zero-norm row. At
# gamma=2 unfloored that was 66 of 300 targets under z=0.05. Flooring the RANK at 0.3 puts a
# hard bound of 0.3**gamma = 0.09 under every target, and 097 measured it as nearly free on
# mse (1.467 -> 1.478 at gamma=2) -- the guard costs nothing and removes the failure mode.
PERT_SHRINK_FLOOR = 0.3
# What to do with a target fewer than two donors carry, i.e. one whose donor cosine is
# UNDEFINED rather than zero. 'keep' (z = 1) was specified on the grounds that it matches the
# per-gene policy; it does not, and that is the whole point. alpha = 1 on a degenerate per-GENE
# row is a genuine no-op -- that row is emitted exactly as an unflagged build emits it. z = 1 on
# an evidence-free PERTURBATION is not a no-op, because z is RELATIVE: every other row is being
# pushed down, so the rows we know least about become the highest-weighted rows on the panel.
# Experiment 097 then measured the cost on h1-test, the only surface where it varies at all
# (orion is 132/132 two-donor): those 4 targets are 13.1% of total predicted amplitude A and
# +0.5% of X, mean per-perturbation cosine -0.0012 against +0.0369 for the other 81, and 'keep'
# takes the arm from 1.239 to 0.958 at gamma=2. See --shrink-pert-missing for why that number
# is n=4 on one surface and should not be read as a magnitude.
PERT_SHRINK_MISSING = "damp"

# The exact cached three-donor operator used to build v9.
FROZEN_PRODUCTION_WHITENER_SHA256 = (
    "5608e807230dea04f93f5e804ed26e8f7285d9487d16a2791283793ffb80fed3"
)
EXPECTED_PM_REPLACEMENT_TARGETS = 248
EXPECTED_CD4_COMPLETE_TARGETS = 276


@dataclass(frozen=True)
class LoadedOperator:
    whitener: Whitener
    columns: np.ndarray
    metadata: dict[str, object]
    file_sha256: str
    metadata_sha256: str
    # The axis this build recomputed from ITS OWN panel, and where the expectation it was
    # checked against came from (experiment 235). Defaulted so that a caller constructing a
    # `LoadedOperator` by hand still type-checks; the production path always fills them.
    axis_sha256: str = ""
    axis_sha256_source: str = ""


def peak_rss_gib() -> float:
    value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return value / GIB if sys.platform == "darwin" else value / 1024**2


def memory_checkpoint(stage: str, max_rss_gib: float) -> None:
    peak = peak_rss_gib()
    print(f"[rss] {stage}: {peak:.3f} GiB", flush=True)
    if max_rss_gib > 0 and peak >= max_rss_gib:
        raise MemoryError(
            f"{stage}: peak RSS {peak:.3f} GiB reached limit {max_rss_gib:.3f} GiB"
        )


def array_sha256(array: np.ndarray) -> str:
    canonical = np.ascontiguousarray(array, dtype="<f8")
    return hashlib.sha256(canonical.tobytes()).hexdigest()


def make_rngs(seed: int) -> tuple[np.random.Generator, np.random.Generator]:
    """Return v9's independent library-pick and count-emission streams."""
    return (
        np.random.default_rng(seed),
        np.random.default_rng(seed + EMISSION_SEED_OFFSET),
    )


def assert_frozen_contract() -> None:
    observed = {
        "pm_blend": float(EXP061.BLEND),
        "pm_min_cpm": float(EXP061.MIN_CPM),
        "pm_spread": float(EXP061.SPREAD),
        "cd4_effect_scale": float(EXP062.EFFECT_SCALE),
        "cd4_gamma": float(EXP062.WHITEN_GAMMA),
    }
    expected = {
        "pm_blend": PM_BLEND,
        "pm_min_cpm": MIN_CPM,
        "pm_spread": SPREAD,
        "cd4_effect_scale": 6.817,
        # EXP062 whitens the CD4 residual with its OWN frozen gamma
        # (`experiments/062-cd4-post-whiten/run.py:40` = 0.65), which is a different
        # object from the main production whitener. This guard exists to detect drift in
        # the imported helper, and it previously compared against WHITEN_GAMMA only
        # because the two happened to be equal in every shipped build. Sweeping the main
        # whitener's gamma must not silently move the CD4 residual too -- that would be a
        # second, unmeasured variable -- so the expected value is pinned to the helper's
        # own frozen constant instead.
        "cd4_gamma": CD4_WHITEN_GAMMA,
    }
    if observed != expected:
        raise ValueError(
            "frozen helper contract drifted:\n"
            + json.dumps({"expected": expected, "observed": observed}, indent=2)
        )


def donor_weights(model: MultiDonorModel) -> dict[str, float]:
    if tuple(model.names) != PRODUCTION_DONORS:
        raise ValueError(
            f"production donor order must be {PRODUCTION_DONORS}, got "
            f"{tuple(model.names)}"
        )
    if DONOR_WEIGHTS is None:
        return {name: 1.0 for name in model.names}
    # `DONOR_WEIGHTS` is keyed by `PRODUCTION_DONORS`, and the order was asserted above, so
    # the lookup cannot silently mis-assign a weight to the wrong donor.
    return {name: float(DONOR_WEIGHTS[name]) for name in model.names}


def main_source_support_audit(
    model: MultiDonorModel, targets: Sequence[str]
) -> dict[str, object]:
    per_donor = {
        name: int(sum(str(target) in model.pert_index[name] for target in targets))
        for name in model.names
    }
    vote_counts = np.asarray(
        [
            sum(str(target) in model.pert_index[name] for name in model.names)
            for target in targets
        ],
        dtype=np.int64,
    )
    return {
        "per_donor_target_counts": per_donor,
        "source_vote_count_histogram": {
            str(value): int(np.sum(vote_counts == value))
            for value in sorted(set(vote_counts.tolist()))
        },
        "targets_with_any_source": int(np.sum(vote_counts > 0)),
        "targets_without_source": int(np.sum(vote_counts == 0)),
    }


def build_v9_pre_whiten_effects(
    model: MultiDonorModel,
    targets: Sequence[str],
    *,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
) -> np.ndarray:
    """Construct the exact pre-whitening effect matrix used by v9."""
    active_weights = dict(weights) if weights is not None else donor_weights(model)
    return np.vstack(
        [
            model.effect(
                str(target),
                alpha,
                beta,
                weights=active_weights,
            )
            for target in targets
        ]
    ).astype(np.float64)


def apply_perturbmatch_k562_replacement(
    baseline: np.ndarray,
    model: MultiDonorModel,
    pm: Mapping[str, object],
    targets: Sequence[str],
    *,
    alpha: float = ALPHA,
    weights: Mapping[str, float] | None = None,
    blend: float = PM_BLEND,
) -> tuple[np.ndarray, dict[str, object]]:
    """Replace one K562 slot while preserving target-level source availability.

    The replacement uses experiment 061's broad gene policy: all official genes
    observed by PerturbMatch may represent the K562 slot.  A PerturbMatch-missing
    gene retains its incumbent value.  At the target level, replacement occurs
    only when both K562-GWPS and PerturbMatch contain the perturbation.
    """
    baseline = np.asarray(baseline, dtype=np.float64)
    if baseline.shape != (len(targets), len(model.target_genes)):
        raise ValueError("baseline shape does not match target/model axes")
    if K562 not in model.names:
        raise ValueError("production model has no K562 slot")

    active_weights = (
        dict(weights) if weights is not None else {name: 1.0 for name in model.names}
    )
    candidate = baseline.copy()
    k_scale = float(model._scale[K562]) or 1.0
    pm_scale = float(pm["scale"])
    pm_target_columns = np.asarray(pm["target_columns"], dtype=np.int64)
    pm_source_columns = np.asarray(pm["source_columns"], dtype=np.int64)
    pm_mean = np.asarray(pm["mean"], dtype=np.float64)
    pm_delta = pm["delta"]
    pm_pert_index = pm["pert_index"]

    if len(pm_target_columns) != len(pm_source_columns):
        raise ValueError("PerturbMatch source/target column maps differ in length")
    if np.any(pm_target_columns < 0) or np.any(
        pm_target_columns >= len(model.target_genes)
    ):
        raise ValueError("PerturbMatch target column lies outside official gene axis")

    support = {
        "both_supported_modified": [],
        "k562_supported_pm_missing_preserved": [],
        "pm_supported_k562_missing_preserved": [],
        "neither_supported_preserved": [],
    }
    source_vote_histogram: dict[str, int] = {}
    for row, target_value in enumerate(targets):
        target = str(target_value)
        has_k562 = target in model.pert_index[K562]
        has_pm = target in pm_pert_index
        if not (has_k562 and has_pm):
            if has_k562:
                key = "k562_supported_pm_missing_preserved"
            elif has_pm:
                key = "pm_supported_k562_missing_preserved"
            else:
                key = "neither_supported_preserved"
            support[key].append(target)
            if not np.array_equal(candidate[row], baseline[row]):
                raise AssertionError(f"{target}: unsupported target was modified")
            continue

        have = [
            name
            for name in model.names
            if target in model.pert_index[name] and active_weights.get(name, 0.0) > 0
        ]
        total_weight = float(sum(active_weights[name] for name in have))
        if total_weight <= 0:
            raise ValueError(f"{target}: no positive incumbent source weight")
        k562_share = float(active_weights[K562]) / total_weight
        source_vote_histogram[str(len(have))] = (
            source_vote_histogram.get(str(len(have)), 0) + 1
        )

        incumbent = np.zeros(len(model.target_genes), dtype=np.float64)
        k_row = model.pert_index[K562][target]
        incumbent[model.tcols[K562]] = model.s[K562][k_row][model.dcols[K562]] / k_scale
        pm_row = pm_pert_index[target]
        pm_standardized = (
            np.asarray(pm_delta[pm_row, pm_source_columns], dtype=np.float64) - pm_mean
        ) / pm_scale

        before = candidate[row].copy()
        candidate[row, pm_target_columns] += (
            alpha
            * k562_share
            * float(model.effect_scale)
            * blend
            * (pm_standardized - incumbent[pm_target_columns])
        )
        missing_gene_mask = np.ones(len(model.target_genes), dtype=bool)
        missing_gene_mask[pm_target_columns] = False
        if not np.array_equal(
            candidate[row, missing_gene_mask], before[missing_gene_mask]
        ):
            raise AssertionError(
                f"{target}: PerturbMatch-missing gene changed before whitening"
            )
        support["both_supported_modified"].append(target)

    count_sum = sum(len(values) for values in support.values())
    if count_sum != len(targets):
        raise AssertionError("PerturbMatch target-support partition is incomplete")
    audit: dict[str, object] = {
        "blend": blend,
        "gene_policy": (
            "experiment-061 broad official-gene support; PM-missing genes retain "
            "the incumbent K562 slot"
        ),
        "targets_total": len(targets),
        **{f"{name}_count": len(values) for name, values in support.items()},
        **support,
        "modified_target_source_vote_histogram": source_vote_histogram,
        "genes_observed_by_perturbmatch": int(len(pm_target_columns)),
        "incumbent_k562_scale": k_scale,
        "perturbmatch_scale": pm_scale,
    }
    return candidate, audit


def operator_axis_sha256(symbols: Iterable[str]) -> str:
    """sha256 over the operator's ORDERED column/gene axis (experiment 235).

    Byte-identical to `scripts/build_union_whitener.py: axis_sha256`, which writes this digest
    into every operator metadata json it produces. The two are separate copies on purpose --
    the builder verifies the sha256 of every helper it imports, so importing the fitting script
    would couple two independently-pinned files -- and 235's known-answer row 3 measures the
    two implementations against each other on a real refit.

    Canonical form: newline-joined UTF-8 of the symbols in axis order, no trailing newline.
    """
    return hashlib.sha256("\n".join(map(str, symbols)).encode("utf-8")).hexdigest()


def production_operator_columns(
    model: MultiDonorModel, panel: Iterable[str]
) -> np.ndarray:
    union: set[str] = set()
    for name in model.names:
        union |= set(map(str, model.target_genes[model.tcols[name]]))
    union -= set(map(str, panel))
    target_index = {str(gene): index for index, gene in enumerate(model.target_genes)}
    return np.asarray(
        sorted(target_index[gene] for gene in union),
        dtype=np.int64,
    )


def load_production_operator(
    model: MultiDonorModel,
    panel: Sequence[str],
    operator_path: Path,
    metadata_path: Path,
    *,
    expected_sha256: str | None = FROZEN_PRODUCTION_WHITENER_SHA256,
    whitener_id: str | None = None,
) -> LoadedOperator:
    columns = production_operator_columns(model, panel)
    panel_set = set(map(str, panel))
    if panel_set & set(map(str, model.target_genes[columns])):
        raise AssertionError("a VCC panel coordinate entered the production operator")

    operator_hash = file_sha256(operator_path)
    metadata_hash = file_sha256(metadata_path)
    if expected_sha256 and operator_hash != expected_sha256:
        raise ValueError(
            f"{operator_path}: SHA-256 {operator_hash} != frozen {expected_sha256}"
        )
    metadata = json.loads(metadata_path.read_text())
    expected_fields = {
        "n_cols": len(columns),
        "donors": list(model.names),
    }
    observed_fields = {key: metadata.get(key) for key in expected_fields}
    if observed_fields != expected_fields:
        raise ValueError(
            "production whitener metadata mismatch:\n"
            + json.dumps(
                {"expected": expected_fields, "observed": observed_fields},
                indent=2,
            )
        )
    basis = metadata.get("basis")
    if not isinstance(basis, list) or len(basis) != 2 or int(basis[1]) != len(columns):
        raise ValueError("production whitener basis metadata is invalid")

    # --- THE PANEL GUARD (experiment 235). Everything above this line passes on the WRONG
    # panel: `n_cols` is 17,942 on every panel measured, the basis shape is derived from it,
    # and the file sha256 only says the operator is the operator. The axis digest is the one
    # quantity that is a function of the panel, and it is checked HERE, before the 1.3 GB
    # eigenbasis is read and 21 minutes are spent on a build that would be wrong in 98.8% of
    # its coordinates with nothing to read it back (`DEF-RUNBOOK` section 0.1; -0.0429 on
    # `pds` alone, 212).
    axis_digest = operator_axis_sha256(map(str, model.target_genes[columns]))
    recorded = metadata.get("axis_sha256")
    axis_source = f"metadata:{metadata_path.name}"
    if not isinstance(recorded, str) or not recorded:
        recorded = LEGACY_OPERATOR_AXIS_SHA256.get(operator_hash)
        axis_source = "LEGACY_OPERATOR_AXIS_SHA256[operator sha256]"
    if not isinstance(recorded, str) or not recorded:
        raise ValueError(
            "REFUSING TO BUILD: this operator carries no axis digest, so the panel it was "
            "fitted for cannot be checked, and a panel mismatch here is SILENT.\n"
            f"  operator      {operator_path} ({operator_hash})\n"
            f"  metadata      {metadata_path}\n"
            f"  this panel's axis_sha256 {axis_digest} over {len(columns)} columns\n"
            "  Fix: re-fit it with scripts/build_union_whitener.py, which writes "
            "'axis_sha256' into the metadata json, or -- only if you have INDEPENDENTLY "
            "established which panel it was fitted for -- add its operator sha256 to "
            "LEGACY_OPERATOR_AXIS_SHA256."
        )
    if recorded != axis_digest:
        raise ValueError(
            "REFUSING TO BUILD: the whitening operator was fitted for a DIFFERENT PANEL.\n"
            f"  operator            {operator_path}\n"
            f"  operator sha256     {operator_hash}\n"
            f"  --whitener-id       {whitener_id!r}\n"
            f"  axis_sha256 recorded ({axis_source}) {recorded}\n"
            f"  axis_sha256 of THIS panel                {axis_digest}\n"
            f"  columns             {len(columns)} (this check is NOT a column count; "
            "n_cols is 17,942 on every panel this project has measured)\n"
            "  The operator's axis is union(donor coverage) - panel, so it is a function of "
            "the panel: using it against another one pairs ~98.8% of its coordinates with the "
            "wrong gene, silently. Re-fit it for this panel "
            "(scripts/build_union_whitener.py --bundle <this bundle> --basis-prune essential), "
            "re-pin it (scripts/def_pin_panel_artifacts.py --apply), and rebuild. "
            "DEF-RUNBOOK section 0.1 and section 2 step 5."
        )

    whitener = Whitener.load(operator_path)
    if whitener.evecs.shape != (len(columns), len(columns)):
        raise ValueError(
            f"operator shape {whitener.evecs.shape} != {(len(columns), len(columns))}"
        )
    print(
        f"[operator] axis_sha256 {axis_digest} over {len(columns)} columns "
        f"MATCHES {axis_source}",
        flush=True,
    )
    return LoadedOperator(
        whitener=whitener,
        columns=columns,
        metadata=metadata,
        file_sha256=operator_hash,
        metadata_sha256=metadata_hash,
        axis_sha256=axis_digest,
        axis_sha256_source=axis_source,
    )


def apply_v9_whitening(
    effects: np.ndarray,
    whitener,
    columns: np.ndarray,
    *,
    gamma: float = WHITEN_GAMMA,
    lam: float = WHITEN_LAM,
) -> tuple[np.ndarray, dict[str, float]]:
    """Apply v9's candidate-specific union-axis preserve normalization."""
    output = np.asarray(effects, dtype=np.float64).copy()
    columns = np.asarray(columns, dtype=np.int64)
    before_all = float(np.median(np.linalg.norm(output, axis=1)))
    before_operator = (
        float(np.median(np.linalg.norm(output[:, columns], axis=1))) or 1.0
    )
    transformed = whitener.apply(output[:, columns], gamma, lam, renorm=False)
    after_operator = float(np.median(np.linalg.norm(transformed, axis=1))) or 1.0
    scale = before_operator / after_operator
    output[:, columns] = transformed * scale
    after_all = float(np.median(np.linalg.norm(output, axis=1)))
    return output, {
        "gamma": gamma,
        "lambda": lam,
        "before_operator_median_norm": before_operator,
        "after_operator_median_norm": after_operator,
        "operator_preserve_scale": scale,
        "before_all_median_norm": before_all,
        "after_all_median_norm": after_all,
    }


def sign_repair_to_consensus(
    whitened: np.ndarray,
    pre_whiten: np.ndarray,
    columns: np.ndarray,
) -> tuple[np.ndarray, dict[str, object]]:
    """Put back the signs the whitening operator flipped (experiment 173).

    `flip = sign(whitened) == -sign(pre_whiten)`, both non-zero. On those entries the
    whitened sign is measurably ANTI-informative (45.6% correct on the h1-test reference DE
    pool against the consensus's 54.4%), so the magnitude is kept and only the sign is taken
    from the consensus. See the `SIGN_REPAIR` block above for the measurement.

    Restricted to `columns` on purpose. `apply_v9_whitening` writes only `output[:, columns]`
    and passes every other coordinate through untouched, so no off-column sign can have been
    flipped by the operator; confining the repair there makes that explicit and means the
    switch can never move a coordinate the operator did not.
    """
    out = np.asarray(whitened, dtype=np.float64).copy()
    columns = np.asarray(columns, dtype=np.int64)
    after = out[:, columns]
    before = np.asarray(pre_whiten, dtype=np.float64)[:, columns]
    flip = (np.sign(after) * np.sign(before)) < 0.0
    out[:, columns] = np.where(flip, -after, after)
    total = float(np.square(out).sum()) or 1.0
    off_axis = np.ones(out.shape[1], dtype=bool)
    off_axis[columns] = False
    return out, {
        "mode": "cons",
        "entries_on_axis": int(flip.size),
        "flipped_entries": int(flip.sum()),
        "flipped_fraction_on_axis": float(flip.mean()),
        "flipped_energy_share": float(np.square(after[flip]).sum() / total),
        "off_axis_columns_untouched": int(off_axis.sum()),
    }


def freeze_complete_mean3(prepared):
    """Keep CD4 residual rows only when all three frozen conditions exist."""
    complete = prepared.target_condition_count == len(CD4_NAMES)
    mean_effect = np.asarray(prepared.mean_effect, dtype=np.float64).copy()
    gene_votes = np.asarray(prepared.gene_votes, dtype=np.uint8).copy()
    mean_effect[~complete] = 0.0
    gene_votes[~complete] = 0
    return EXP062.CD4Prepared(
        mode=prepared.mode,
        mean_effect=mean_effect,
        gene_votes=gene_votes,
        target_condition_count=np.asarray(
            prepared.target_condition_count, dtype=np.uint8
        ).copy(),
        condition_details=prepared.condition_details,
    )


def prepare_frozen_cd4(
    root: Path,
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
):
    banks = tuple(
        EXP062.Bank.load(name, root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    prepared = EXP062.prepare_cd4_mean(
        banks,
        targets,
        official_genes,
        set(map(str, panel)),
        CD4_MODE,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    return banks, freeze_complete_mean3(prepared)


def cd4_support_audit(cd4, fraction: float = CD4_FRACTION) -> dict[str, object]:
    counts = np.asarray(cd4.target_condition_count, dtype=np.uint8)
    complete = counts == len(CD4_NAMES)
    return {
        "conditions": list(CD4_NAMES),
        "mode": str(cd4.mode),
        "fraction": float(fraction),
        "target_condition_count_histogram": {
            str(value): int(np.sum(counts == value))
            for value in sorted(set(counts.tolist()))
        },
        "fully_supported_targets": int(np.sum(complete)),
        "incomplete_targets_preserved_before_global_scale": int(np.sum(~complete)),
        "supported_entries": int(np.count_nonzero(cd4.gene_votes)),
    }


def solve_cd4_fraction(whitened: np.ndarray, cd4) -> tuple[float, dict[str, object]]:
    """Report experiment 099's `f` for this build, and solve the fraction when one is asked for.

    099 (`postwhiten.py`) normalises the candidate residual so its median row norm equals the
    whitened base's, then adds `f x` that.  So `f` IS the ratio of the two median row norms,
    over the full gene axis, with both medians taken over every row -- including rows that
    carry no residual.  The builder adds `CD4_FRACTION * mean_effect` masked by `gene_votes`
    instead, which is a differently-normalised object, so the two knobs do not share a scale.

    `norm(fraction * residual)` is linear in `fraction`, so measuring `f` once at the current
    fraction pins the whole line and the target solves in closed form.  Reported either way.
    """
    added = np.where(cd4.gene_votes > 0, cd4.mean_effect, 0.0)
    residual_norm = float(np.median(np.linalg.norm(added, axis=1)))
    base_norm = float(np.median(np.linalg.norm(np.asarray(whitened), axis=1)))
    # `f` per unit fraction: the ratio the fraction is multiplied into.
    f_per_unit = residual_norm / base_norm if base_norm else 0.0
    frame = {
        "definition": (
            "median row norm of the added residual over median row norm of the whitened base, "
            "full gene axis, medians over all rows (experiment 099 postwhiten.py)"
        ),
        "residual_median_row_norm_at_unit_fraction": residual_norm,
        "whitened_base_median_row_norm": base_norm,
        "f_per_unit_fraction": f_per_unit,
        "cd4_fraction_constant": CD4_FRACTION,
        "effective_f_at_cd4_fraction": CD4_FRACTION * f_per_unit,
        "target_f_requested": CD4_TARGET_F,
    }
    if CD4_TARGET_F is None:
        frame["fraction_source"] = "VCC_CD4_FRACTION"
        return CD4_FRACTION, frame
    if f_per_unit <= 0.0:
        raise ValueError(
            "VCC_CD4_TARGET_F was set but the CD4 residual has zero norm; nothing to solve"
        )
    solved = CD4_TARGET_F / f_per_unit
    frame["fraction_source"] = "VCC_CD4_TARGET_F"
    frame["solved_fraction"] = solved
    return solved, frame


def solve_carrying_row_fraction(
    whitened: np.ndarray, residual, target_g: float
) -> tuple[float, dict[str, object]]:
    """THE COVERAGE-AWARE SIZING RULE.  218 section 7.1's `solve_fraction_carrying`, verbatim.

    `solve_cd4_fraction` takes the median row norm over EVERY panel row, so a source covering
    <= 50% of the panel has median exactly 0, `f_per_unit` is 0, and the solved fraction is 0 --
    the residual is identically absent and nothing catches it.  For a partial-coverage source
    the meaningful denominator is the rows it actually carries, so `g` reads "on the rows this
    source covers, the added norm is `g` x the base's typical row norm".

    Above 50% coverage the two rules coincide.  **This function is never applied to CD4** --
    `apply_frozen_cd4_residual` still sizes CD4 through `solve_cd4_fraction` unchanged.
    """
    added = np.where(
        np.asarray(residual.gene_votes) > 0, np.asarray(residual.mean_effect), 0.0
    )
    norms = np.linalg.norm(added, axis=1)
    live = norms[norms > 0]
    carrying_median = float(np.median(live)) if live.size else 0.0
    base_norm = float(np.median(np.linalg.norm(np.asarray(whitened), axis=1)))
    per_unit = carrying_median / base_norm if base_norm else 0.0
    median_all_rows = float(np.median(norms))
    frame = {
        "definition": (
            "median row norm of the added residual over the rows it CARRIES, divided by the "
            "median row norm of the whitened base over all rows (experiment 218 section 7.1)"
        ),
        "residual_median_row_norm_over_carrying_rows": carrying_median,
        "residual_median_row_norm_over_all_rows": median_all_rows,
        "whitened_base_median_row_norm": base_norm,
        "g_per_unit_fraction": per_unit,
        "target_g_requested": float(target_g),
        "carrying_rows": int(live.size),
        "coverage_frac": float(live.size / norms.size) if norms.size else 0.0,
        # The number this whole change exists to make visible: what production's median-over-
        # all-rows rule would have solved for this source.  It is 0 for anything below 50%.
        "shipped_median_rule_would_solve": (
            float(target_g / (median_all_rows / base_norm))
            if median_all_rows > 0.0 and base_norm
            else 0.0
        ),
    }
    if per_unit <= 0.0:
        raise ValueError(
            "the coverage-aware solve found a residual with no carrying rows at all; "
            "nothing to size"
        )
    frame["solved_fraction"] = float(target_g / per_unit)
    return float(target_g / per_unit), frame


def single_source_residual(
    bank,
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel_targets: set[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
):
    """218's `multi_source_residual` for one bank -- `prepare_cd4_mean` generalised.

    Construction-identical to production's CD4 path (centre on the bank's own row mean, the
    `gene_norm` per-gene scale, rescale the median centred row norm to `EFFECT_SCALE`) except
    that the official-axis mapping goes through `map_bank_to_axis`, which keeps the LAST source
    column for a duplicated symbol instead of raising.  218's KA-4 asserted the two agree
    bit-for-bit on the three CD4 banks (max abs diff 0.0e+00, votes identical, 4/4 folds).
    """
    amap = EXP062.map_bank_to_axis(bank, official_genes)
    keep = np.asarray(
        [
            str(official_genes[position]) not in panel_targets
            for position in amap.unique_target_positions.tolist()
        ],
        dtype=bool,
    )
    columns = amap.unique_source_columns[keep]
    positions = amap.unique_target_positions[keep]
    mean = EXP062.stream_column_mean(
        bank, columns, chunk_rows=chunk_rows, max_rss_gib=max_rss_gib
    )
    raw_median, gene_scale, norm_median = EXP062.centered_statistics(
        bank,
        columns,
        mean,
        need_gene_scale=True,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    multiplier = EXP062.EFFECT_SCALE / norm_median
    value_sum = np.zeros((len(targets), len(official_genes)), dtype=np.float64)
    votes = np.zeros(value_sum.shape, dtype=np.uint8)
    counts = np.zeros(len(targets), dtype=np.uint8)
    delta = np.load(bank.delta_path, mmap_mode="r")
    covered = 0
    for index, target in enumerate(targets):
        row = bank.pert_pos.get(str(target))
        if row is None:
            continue
        centered = np.asarray(delta[row, columns], dtype=np.float64) - mean
        centered = centered / gene_scale
        value_sum[index, positions] += centered * multiplier
        votes[index, positions] += 1
        counts[index] += 1
        covered += 1
    available = votes > 0
    mean_effect = np.zeros_like(value_sum)
    mean_effect[available] = value_sum[available] / votes[available]
    detail = {
        "condition": bank.name,
        "path": str(Path(bank.path).relative_to(ROOT)),
        "background_perturbations": int(bank.shape[0]),
        "mapped_non_panel_genes": int(columns.size),
        "duplicate_symbols": len(amap.duplicate_symbols),
        "covered_targets": int(covered),
        # Experiment 235. The panel is disjoint from A/B/C's on Oct 22 and 225 measured two
        # registered sources going from 273 and 182 targets to EXACTLY ZERO on it, so a
        # source's own coverage is a first-class output, not a footnote.
        "coverage_frac": float(covered) / len(targets) if len(targets) else 0.0,
        "raw_median_centered_norm": float(raw_median),
        "gene_normalized_median_before_rescale": float(norm_median),
    }
    return (
        EXP062.CD4Prepared(
            mode=f"gene_norm:{bank.name}",
            mean_effect=mean_effect,
            gene_votes=votes,
            target_condition_count=counts,
            condition_details=[detail],
        ),
        detail,
    )


def union_residual(sources):
    """218's `union_rows`: per-gene mean over whichever sources carry a coordinate."""
    total = np.zeros_like(np.asarray(sources[0].mean_effect), dtype=np.float64)
    votes = np.zeros(total.shape, dtype=np.uint16)
    counts = np.zeros(
        np.asarray(sources[0].target_condition_count).shape, dtype=np.uint16
    )
    for source in sources:
        carried = np.asarray(source.gene_votes) > 0
        total[carried] += np.asarray(source.mean_effect)[carried]
        votes[carried] += 1
        counts += (np.asarray(source.target_condition_count) > 0).astype(np.uint16)
    mean_effect = np.zeros_like(total)
    live = votes > 0
    mean_effect[live] = total[live] / votes[live]
    return EXP062.CD4Prepared(
        mode="union:" + "+".join(str(source.mode) for source in sources),
        mean_effect=mean_effect,
        gene_votes=np.minimum(votes, 255).astype(np.uint8),
        target_condition_count=np.minimum(counts, 255).astype(np.uint8),
        condition_details=[],
    )


def verify_union_source_hashes() -> dict[str, dict[str, str]]:
    """Rule 11, ENFORCED: every selected bank's four files must hash to the registry's value.

    Three of the registered banks live under `experiments/*/work/`, a tree other experiments
    write to.  Recording a sha catches a changed file only if something compares it; this
    compares it, before a single byte of a 35-minute build is spent, and refuses to run
    otherwise.  Returns the measured table so the caller does not hash twice.
    """
    measured: dict[str, dict[str, str]] = {}
    for name, path in H1_UNION_SOURCES:
        expected = UNION_SOURCE_REGISTRY[name]["sha256"]
        got = {filename: file_sha256(path / filename) for filename in H1_UNION_FILES}
        bad = {k: (expected[k], v) for k, v in got.items() if expected[k] != v}
        if bad:
            raise ValueError(
                f"union source {name!r} at {path} does not match UNION_SOURCE_REGISTRY: "
                + "; ".join(f"{k}: expected {e[:16]}… got {g[:16]}…" for k, (e, g) in bad.items())
            )
        measured[name] = got
    return measured


def prepare_h1_union(
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel: Sequence[str],
    *,
    max_rss_gib: float,
):
    """The union of the selected registry banks on the official axis, panel columns excluded.

    `chunk_rows` is pinned at 218's 256 rather than taken from `--chunk-rows`: the streamed
    column mean accumulates chunk by chunk, so the chunk size is a (16th-decimal) numerical
    choice and this residual has to reproduce the object 218 measured.
    """
    verify_union_source_hashes()
    panel_targets = set(map(str, panel))
    prepared = []
    details = []
    for name, path in H1_UNION_SOURCES:
        bank = EXP062.Bank.load(name, path)
        source, detail = single_source_residual(
            bank,
            targets,
            official_genes,
            panel_targets,
            chunk_rows=256,
            max_rss_gib=max_rss_gib,
        )
        # THE PER-SOURCE SILENT ZERO (experiment 235). `solve_carrying_row_fraction` raises
        # only when the WHOLE union carries nothing; a single selected source that carries
        # nothing contributes exactly zero to the union and every check downstream still
        # passes, so the meta would describe an arm that is not the arm that ran. 225
        # measured `kolf2` 273 -> 0 and `feng_ipsc` 182 -> 0 against October's panel, which is
        # disjoint from A/B/C's by construction -- so this is the expected October failure,
        # not a hypothetical one. It is loud now.
        if int(detail["covered_targets"]) == 0:
            raise ValueError(
                f"union source {name!r} carries 0 of {len(targets)} panel targets, so it "
                "would contribute EXACTLY NOTHING to this build while still being named in "
                f"VCC_UNION_SOURCES and recorded in the meta.\n  path {path}\n"
                "  Drop it from VCC_UNION_SOURCES, or select a source that covers this "
                "panel. On a new panel, check coverage first: 225 section 8.3 measured two "
                "registered sources falling from 273 and 182 targets to zero."
            )
        prepared.append(source)
        details.append(detail)
        del bank
    union = union_residual(prepared)
    del prepared
    gc.collect()
    return union, details


def h1_union_support_audit(union, details) -> dict[str, object]:
    counts = np.asarray(union.target_condition_count, dtype=np.uint8)
    return {
        # Rule 11: a path AND a sha, never a description.  These two directories are not in
        # the frozen-input table (that table is the panel bundle plus the three production
        # donors plus PerturbMatch plus CD4, and adding to it is a separate, reviewable
        # re-pin), so the build records their content hashes here instead of asserting them.
        "sources": [
            {
                "name": name,
                "path": str(path.relative_to(ROOT)),
                "sha256": {
                    filename: file_sha256(path / filename)
                    for filename in H1_UNION_FILES
                },
            }
            for name, path in H1_UNION_SOURCES
        ],
        "source_details": details,
        "mode": str(union.mode),
        "sizing_rule": "median over CARRYING rows (coverage-aware)",
        "target_g": float(H1_UNION_G),
        "carrying_targets": int(np.sum(counts > 0)),
        # Correct at any arity, unlike `targets_carried_by_both_sources` below, which counts
        # `== 2` and is only meaningful for the shipped pair. Kept because it is in the
        # shipped `h1union` meta. (experiment 235)
        "targets_by_source_count": {
            str(k): int(np.sum(counts == k)) for k in range(0, len(H1_UNION_SOURCES) + 1)
        },
        "source_coverage": {
            str(detail["condition"]): int(detail["covered_targets"]) for detail in details
        },
        "targets_carried_by_both_sources": int(np.sum(counts == 2)),
        "supported_entries": int(np.count_nonzero(union.gene_votes)),
    }


def apply_frozen_cd4_residual(
    whitened: np.ndarray,
    cd4,
    operator_columns: np.ndarray,
    panel_positions: np.ndarray,
    *,
    extra=None,
    extra_audit: dict[str, object] | None = None,
) -> tuple[np.ndarray, dict[str, object], dict[str, object] | None]:
    """Apply frozen post-whitening CD4 and enforce zero panel support.

    `extra` is the optional second residual of `VCC_H1_UNION_G`.  When it is given, CD4 and the
    second source are solved against the SAME denominator (each by its own sizing rule), merged
    into one residual carrying `frac_cd4 * cd4 + frac_extra * extra`, and added in ONE step
    followed by production's single global preserve rescale -- 218 section 7.2's `merged()`.
    CD4's own solved fraction is bit-identical either way; the only thing the second residual
    can move in CD4's own audit is the shared preserve scalar, which is what "one global
    rescale" means.
    """
    panel_positions = np.asarray(panel_positions, dtype=np.int64)
    if np.any(cd4.gene_votes[:, panel_positions] != 0):
        raise AssertionError("a VCC panel coordinate entered CD4 gene support")
    if np.any(cd4.mean_effect[:, panel_positions] != 0):
        raise AssertionError("a VCC panel coordinate entered the CD4 residual")

    fraction, frame = solve_cd4_fraction(whitened, cd4)
    if extra is None:
        added = cd4
        applied = fraction
        extra_frame = None
        support = np.asarray(cd4.gene_votes) > 0
    else:
        if np.any(np.asarray(extra.gene_votes)[:, panel_positions] != 0):
            raise AssertionError("a VCC panel coordinate entered the second residual support")
        if np.any(np.asarray(extra.mean_effect)[:, panel_positions] != 0):
            raise AssertionError("a VCC panel coordinate entered the second residual")
        extra_fraction, extra_frame = solve_carrying_row_fraction(
            whitened, extra, float(H1_UNION_G)
        )
        cd4_votes = np.asarray(cd4.gene_votes) > 0
        extra_votes = np.asarray(extra.gene_votes) > 0
        merged_effect = np.where(cd4_votes, np.asarray(cd4.mean_effect), 0.0) * fraction
        merged_effect += (
            np.where(extra_votes, np.asarray(extra.mean_effect), 0.0) * extra_fraction
        )
        support = cd4_votes | extra_votes
        added = EXP062.CD4Prepared(
            mode="merged",
            mean_effect=merged_effect,
            gene_votes=support.astype(np.uint8),
            target_condition_count=np.maximum(
                np.asarray(cd4.target_condition_count),
                np.asarray(extra.target_condition_count),
            ).astype(np.uint8),
            condition_details=[],
        )
        applied = 1.0
        extra_frame["applied_fraction"] = float(extra_fraction)
        extra_frame["entries_outside_cd4_support"] = int(
            np.count_nonzero(extra_votes & ~cd4_votes)
        )
        extra_frame["targets_outside_cd4_support"] = int(
            np.sum(
                (np.asarray(extra.target_condition_count) > 0)
                & (np.asarray(cd4.target_condition_count) == 0)
            )
        )

    output, preserve = EXP062.add_post_whiten_residual(
        whitened,
        added,
        applied,
        np.asarray(operator_columns, dtype=np.int64),
    )
    scale = float(preserve["post_add_global_preserve_scale"])
    np.testing.assert_allclose(
        output[:, panel_positions],
        np.asarray(whitened)[:, panel_positions] * scale,
        rtol=0.0,
        atol=0.0,
        err_msg="panel coordinates changed by more than the preserve scalar",
    )
    unavailable = ~support
    np.testing.assert_allclose(
        output[unavailable],
        np.asarray(whitened)[unavailable] * scale,
        rtol=0.0,
        atol=0.0,
        err_msg="unsupported residual entries changed by more than the preserve scalar",
    )
    if extra is None:
        reported_preserve: dict[str, object] = preserve
    else:
        # `add_post_whiten_residual` echoes the fraction it was handed, which is 1.0 for a
        # merged add.  Report CD4's OWN solved fraction here instead, so this key keeps the
        # meaning it has in every prior build's meta and a rule-10 diff cannot read the merge
        # as a change in CD4's dose.  The two norm scalars below it DO move -- one global
        # preserve over one merged add is exactly what 218 measured.
        reported_preserve = {**preserve, "fraction": fraction}
    cd4_audit = {
        **cd4_support_audit(cd4, fraction),
        "residual_frame": frame,
        "post_add_preserve": reported_preserve,
        "panel_residual_nonzero_entries": int(
            np.count_nonzero(cd4.mean_effect[:, panel_positions])
        ),
    }
    if extra is None:
        return output, cd4_audit, None
    union_audit = {
        **(extra_audit or {}),
        "residual_frame": extra_frame,
        "merged_add": True,
        "merged_preserve": preserve,
        "panel_residual_nonzero_entries": int(
            np.count_nonzero(np.asarray(extra.mean_effect)[:, panel_positions])
        ),
    }
    return output, cd4_audit, union_audit


def build_production_effects(
    model: MultiDonorModel,
    targets: Sequence[str],
    whitener,
    operator_columns: np.ndarray,
    panel_positions: np.ndarray,
    *,
    pm: Mapping[str, object] | None,
    cd4=None,
    h1_union=None,
    h1_union_audit: dict[str, object] | None = None,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Build one context-independent production effect matrix."""
    active_weights = dict(weights) if weights is not None else donor_weights(model)
    baseline_pre = build_v9_pre_whiten_effects(
        model,
        targets,
        alpha=alpha,
        beta=beta,
        weights=active_weights,
    )
    if pm is None:
        candidate_pre = baseline_pre.copy()
        pm_audit: dict[str, object] = {"enabled": False}
    else:
        candidate_pre, pm_audit = apply_perturbmatch_k562_replacement(
            baseline_pre,
            model,
            pm,
            targets,
            alpha=alpha,
            weights=active_weights,
        )
        pm_audit = {"enabled": True, **pm_audit}

    whitened, whitening_audit = apply_v9_whitening(
        candidate_pre,
        whitener,
        operator_columns,
    )
    if SIGN_REPAIR is None:
        sign_repair_audit: dict[str, object] = {"enabled": False}
    else:
        whitened, sign_repair_audit = sign_repair_to_consensus(
            whitened, candidate_pre, operator_columns
        )
        sign_repair_audit = {"enabled": True, **sign_repair_audit}
    union_audit: dict[str, object] | None = None
    if cd4 is None:
        final = whitened
        cd4_audit: dict[str, object] = {"enabled": False}
    else:
        final, cd4_audit, union_audit = apply_frozen_cd4_residual(
            whitened,
            cd4,
            operator_columns,
            panel_positions,
            extra=h1_union,
            extra_audit=h1_union_audit,
        )
        cd4_audit = {"enabled": True, **cd4_audit}

    if not np.all(np.isfinite(final)):
        raise ValueError("final effect matrix contains non-finite values")
    return final, {
        "baseline_pre_sha256": array_sha256(baseline_pre),
        "candidate_pre_sha256": array_sha256(candidate_pre),
        "post_whiten_sha256": array_sha256(whitened),
        "final_sha256": array_sha256(final),
        "perturbmatch": pm_audit,
        "whitening": whitening_audit,
        # Emitted only when the switch is on, exactly as `mass_fix` is in the config block:
        # an always-present `{"enabled": false}` would make every default build's meta differ
        # from every shipped one by a key, which is the diff rule 10 relies on.
        **({} if SIGN_REPAIR is None else {"sign_repair": sign_repair_audit}),
        "cd4": cd4_audit,
        # Emitted only when `VCC_H1_UNION_G` is set, exactly as `mass_fix` and `sign_repair`
        # are: an always-present block would make every default build's audit differ from
        # every shipped one by a key, which is the diff rule 10 relies on.
        **({} if union_audit is None else {"h1_union": union_audit}),
    }


def _required_files(directory: Path, names: Sequence[str]) -> list[Path]:
    paths = [directory / name for name in names]
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("missing required inputs: " + ", ".join(missing))
    return paths


def frozen_input_paths(
    bundle: Path,
    contexts: Sequence[str],
    donor_paths: Mapping[str, Path],
    perturbmatch: Path,
    cd4_root: Path,
    operator_path: Path,
    operator_metadata_path: Path,
) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for name in bundle_file_names(contexts):
        paths[f"bundle/{name}"] = bundle / name
    for donor_name, directory in donor_paths.items():
        for name in ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json"):
            paths[f"donors/{donor_name}/{name}"] = directory / name
    for name in ("delta.npy", "genes.txt", "perts.txt", "meta.json"):
        paths[f"perturbmatch/{name}"] = perturbmatch / name
    for name in CD4_NAMES:
        for filename in (
            "delta.npy",
            "b_ctrl.npy",
            "genes.txt",
            "perts.txt",
            "meta.json",
        ):
            paths[f"cd4/{name}/{filename}"] = cd4_root / name / filename
    paths["whitener/whitener_union.npz"] = operator_path
    paths["whitener/whitener_union.json"] = operator_metadata_path
    return paths


def collect_input_hashes(
    bundle: Path,
    contexts: Sequence[str],
    donor_paths: Mapping[str, Path],
    perturbmatch: Path,
    cd4_root: Path,
    operator_path: Path,
    operator_metadata_path: Path,
    whitener_id: str = "union",
    frozen_bundle: Mapping[str, object] | None = None,
) -> dict[str, str]:
    """Verify and return every frozen data hash before any model input is read."""
    return verify_expected_hashes(
        frozen_input_paths(
            bundle,
            contexts,
            donor_paths,
            perturbmatch,
            cd4_root,
            operator_path,
            operator_metadata_path,
        ),
        whitener_expected_hashes(whitener_id, frozen_bundle),
        label="frozen model input",
    )


def repin_frozen_inputs(
    bundle: Path,
    path: Path = FROZEN_BUNDLE_PATH,
    *,
    now: str | None = None,
) -> dict[str, object]:
    """Write the frozen bundle table for the bundle currently on disk. Never builds.

    The table exists to freeze inputs, so re-pinning has to be a deliberate, separately
    invoked, reviewable act -- never something a build does for itself when a hash fails to
    match. This prints the old -> new diff and returns the new table.

    `panel_support` (the PerturbMatch / CD4 target-coverage pair the build asserts) is a
    property of the PANEL, not of the files, so it is preserved only when `pert_counts.csv`
    and `gene_names.csv` are byte-identical to the pinned ones and is set to `null`
    otherwise. `null` is not a bypass: the build then prints the observed values and
    records `panel_support_pinned: false` in its meta, and `--repin-panel-support` freezes
    them from one effects-only run.
    """
    contexts = read_bundle_contexts(bundle)
    manifest = json.loads((bundle / "manifest.json").read_text())
    previous = json.loads(path.read_text()) if path.is_file() else {}
    previous_files = dict(previous.get("files", {}))
    files = {
        f"bundle/{name}": file_sha256(bundle / name)
        for name in bundle_file_names(contexts)
    }
    panel_unchanged = all(
        previous_files.get(key) == files[key]
        for key in ("bundle/pert_counts.csv", "bundle/gene_names.csv")
    )
    import pandas as pd  # local: only the re-pin path needs it

    table = {
        "_comment": previous.get("_comment", []),
        "panel_id": str(manifest.get("panel_id", "")),
        "contexts": list(contexts),
        "n_constructs": int(len(pd.read_csv(bundle / "pert_counts.csv")[PERT_COL])),
        "cells_per_pert": int(manifest.get("cells_per_pert", 400)),
        "files": {key: files[key] for key in sorted(files)},
        "panel_support": previous.get("panel_support") if panel_unchanged else None,
        "repinned_at": now or __import__("datetime").datetime.now().isoformat(timespec="seconds"),
        "supersedes": file_sha256(path) if path.is_file() else None,
    }
    changed = [
        key
        for key in sorted(set(previous_files) | set(files))
        if previous_files.get(key) != files.get(key)
    ]
    print(f"[repin] {path}", flush=True)
    print(f"[repin] contexts {list(previous.get('contexts', []))} -> {list(contexts)}",
          flush=True)
    for key in changed:
        print(f"[repin]   {key}: {previous_files.get(key, '<absent>')} -> "
              f"{files.get(key, '<absent>')}", flush=True)
    if not changed:
        print("[repin]   no file hash changed", flush=True)
    print(f"[repin] panel_support: {previous.get('panel_support')} -> "
          f"{table['panel_support']}"
          f"{'' if panel_unchanged else '  (panel files moved; re-pin it from a build)'}",
          flush=True)
    atomic_write_text(path, json.dumps(table, indent=2) + "\n")
    load_frozen_bundle(path)  # the written table must pass its own validator
    return table


def collect_implementation_hashes() -> dict[str, str]:
    paths = (
        Path(__file__).resolve(),
        ROOT / "scripts/build_submission.py",
        ROOT / "scripts/make_submission.py",
        ROOT / "src/vcc2026/donor_transfer.py",
        ROOT / "src/vcc2026/multi_donor.py",
        ROOT / "src/vcc2026/whiten.py",
        EXP061_PATH,
        EXP062_PATH,
    )
    return {str(path): file_sha256(path) for path in paths}


def mass_consistent_effect(
    gated: np.ndarray,
    b_ctrl: np.ndarray,
    gate: np.ndarray,
    mode: str = "minnorm",
    iters: int = 64,
    rtol: float = 1e-12,
) -> np.ndarray:
    """Add the cheapest correction making `sum(expm1(b_ctrl + h)) == TS` (experiment 156).

    With the constraint met, `profile_from_effect`'s final division by `total` is exactly the
    identity and the effect the scorer reads back is exactly the effect we designed. See the
    `MASS_FIX` block above for why the constraint exists, why scaling cannot satisfy it, and
    why the direction of the correction is the whole result.

    `minnorm` takes the minimum-L2 step, `eps ∝ -exp(b + h)` on the gate. `const` takes a
    uniform log shift on the gate and exists only as the measured control -- it satisfies the
    identical scalar constraint and is decisively worse, which is the evidence that this is a
    geometry effect and not a normalisation bookkeeping detail.

    The correction is confined to `gate`, so a gene the CPM gate zeroed stays at exactly zero
    and no effect is put on a gene no metric can see.
    """
    if mode not in {"minnorm", "const"}:
        raise ValueError(f"unknown mass-fix mode {mode!r}")
    h = np.asarray(gated, dtype=np.float64).copy()
    b = np.asarray(b_ctrl, dtype=np.float64)
    on = np.asarray(gate, dtype=bool)
    for _ in range(iters):
        # The renderer CLIPS the profile at zero, so the mass being constrained is the clipped
        # sum, and a clipped gene has zero derivative -- it cannot carry any of the correction
        # and its read-back is -b_ctrl no matter what we do. Matching the clip here is what
        # makes the read-back exact on every gene that is not clipped.
        live = on & (b + h > 0.0)
        residual = float(np.clip(np.expm1(b + h), 0.0, None).sum() - TS)
        if abs(residual) <= rtol * TS:
            break
        w = np.where(live, np.exp(b + h), 0.0)
        denom = float(w @ w) if mode == "minnorm" else float(w.sum())
        if not np.isfinite(denom) or denom <= 0.0:
            break
        step = -residual / denom
        h = h + (step * w if mode == "minnorm" else np.where(live, step, 0.0))
    return h


def profile_from_effect(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    min_cpm: float = MIN_CPM,
    effect_scale: float = 1.0,
    gene_scale: np.ndarray | None = None,
    shrink_center: np.ndarray | None = None,
) -> np.ndarray:
    """Apply the v9 CPM gate and convert a log-pseudobulk effect to abundance.

    `min_cpm` was 5.0 through v10. Experiment 088 measured 10.0 as a sharp
    optimum on h1-test (+0.00455 avg6 over 6 seeds, sd 0.00058), improving
    reach, pds and nmae together. Commit 269f268 had previously settled this at
    5.0, but that was on the g=0.85 inter-axis operator, before the union axis
    and PerturbMatch+CD4 existed.
    """
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gate = (cpm > min_cpm) & (b_ctrl > 0)
    gated = np.where(gate, effect, 0.0)
    # `effect_scale` shrinks the predicted effect toward the control profile in LOG space,
    # which is where the metric's own distance is measured. 1.0 reproduces every prior build
    # bit-for-bit; 0.0 emits the context's control profile unchanged (the board's measured
    # `expr_mse_unbiased_capped_norm` floor of ~1.042, per team "Kun Lab"). See
    # experiments/096-mse-floor/.
    if effect_scale != 1.0:
        gated = gated * float(effect_scale)
    # `gene_scale` is the same shrinkage with one factor per GENE instead of one for the
    # panel; `--shrink-key agreement` supplies `agreement_shrink_vector(...)`. None (the
    # default) is bit-identical to every prior build. It is applied AFTER `effect_scale`,
    # and both compose multiplicatively with the VCC_ALPHA env var, which scales the same
    # effect earlier (before whitening) -- three knobs on one quantity, so do not
    # double-apply an amplitude cut.
    #
    # `shrink_center` is what a fully shrunk gene converges to. None means zero, i.e. the
    # context's control profile. `--shrink-toward mean` supplies the grand mean predicted
    # effect `mbar`, which is the metric's own 0-anchor, and the interpolation becomes the
    # James-Stein form `mbar + a*(e - mbar)`. It is written as `a*e + (1-a)*mbar` so a gene
    # at `a == 1` is bit-identical to an unflagged build rather than merely close.
    #
    # NOTE the centre is re-gated with the SAME mask. A gene the CPM gate zeroed must stay
    # at exactly zero; centring it without re-gating would resurrect it at `(1-a)*mbar_g`
    # and put effect on a gene no metric can see. Because the mask is a per-gene multiply
    # and the interpolation is affine per gene, gating before or after the shrink is the
    # same operation, and the same holds for `effect_scale` -- which is why it is applied to
    # the centre as well, keeping the two knobs commutative.
    if gene_scale is None:
        if shrink_center is not None:
            raise ValueError("shrink_center requires gene_scale")
    else:
        factors = np.asarray(gene_scale, dtype=np.float64)
        if factors.shape != gated.shape:
            raise ValueError(
                f"gene_scale must be one value per gene, got {factors.shape} for "
                f"{gated.shape}"
            )
        if shrink_center is None:
            gated = gated * factors
        else:
            center = np.asarray(shrink_center, dtype=np.float64)
            if center.shape != gated.shape:
                raise ValueError(
                    f"shrink_center must be one value per gene, got {center.shape} for "
                    f"{gated.shape}"
                )
            center = np.where(gate, center, 0.0)
            if effect_scale != 1.0:
                center = center * float(effect_scale)
            gated = factors * gated + (1.0 - factors) * center
    if MASS_FIX is not None:
        gated = mass_consistent_effect(gated, b_ctrl, gate, mode=MASS_FIX)
    profile = np.clip(np.expm1(b_ctrl + gated), 0.0, None) / TS
    total = float(profile.sum())
    if total <= 0:
        raise ValueError("degenerate profile")
    return profile / total



def donor_agreement_pairwise(
    model: "MultiDonorModel",
    pert: str,
) -> np.ndarray:
    """`donor_agreement` generalised over donor PAIRS instead of all-donor unanimity.

        agree_g = mean over unordered pairs (i, j) of carrying donors of
                  sqrt(|s_i[g]| * |s_j[g]|) * 1[sign(s_i[g]) == sign(s_j[g])] * 1[both nonzero]

    WHY THIS EXISTS, and why it is a separate function rather than an edit. `src/vcc2026` is
    frozen by `FROZEN_IMPLEMENTATION_SHA256`, and `MultiDonorModel.donor_agreement` is
    board-validated in a shipped arm (`--keep-key agreement`, +0.00678 overall as a full
    ordering). It is not touched. This is an alternative key for `--shrink-key` only.

    THE PROBLEM IT SOLVES. `donor_agreement` demands sign unanimity across EVERY carrying donor
    and zeroes any gene missing from any of their axes, so its graded fraction collapses as
    donors are added. Measured on the real challenge panel (experiment 097's `chal_alpha.py`,
    and reproducible from `mech_shrink_coverage.py`): the 259 three-donor targets carry a
    nonzero agreement on a median 10.6 % of the 18,533 panel genes, while the 35 two-donor
    targets carry one on 54.5 %. So a target we have LESS evidence about keeps roughly five
    times as many genes near full amplitude -- backwards on its face, and working directly
    against the per-perturbation weight, which is trying to damp exactly those rows.

    WHAT IT FIXES AND WHAT IT DOES NOT. Both harnesses use two donors, and WITH EXACTLY TWO
    DONORS THESE TWO FUNCTIONS ARE IDENTICAL -- one pair, and unanimity over that pair is the
    pair's own sign test (asserted with `assert_array_equal` in the test suite). So experiment
    097's 1.101x / 1.478x apply to BOTH and cannot distinguish them; the two separate only at
    production's third donor.

    MEASURED at three donors, and THE FROZEN `agreement` WINS -- prefer it, and treat this
    function as the documented fallback rather than the recommendation.

    Median graded (nonzero) fraction of the 18,533 panel genes, on the real challenge panel
    (`mech_shrink_coverage.py`), reproduced on h1-test to three decimals by experiment 097:

        3-donor targets (259)   unanimity 0.106    pairwise 0.756
        2-donor targets  (35)   unanimity 0.545    pairwise 0.545   (identical, as above)
        harness, 2 donors       0.208 (h1-test) / 0.206 (orion)

    A gene needs only ONE pair to agree, so this key is a UNION over pairs where unanimity is an
    intersection: adding a donor RAISES its coverage. Both keys are therefore extrapolations
    away from the 0.208 that was measured, in opposite directions.

    Experiment 097's `pairwise3.py` then measured which extrapolation is better, on h1-test with
    `orion_hct116` added as a legal third donor -- which makes that triple EXACTLY production's
    donor set. Multiplier at eta=1: 2-donor 1.101, 3-donor unanimity **1.278**, 3-donor pairwise
    1.213; stacked with the per-perturbation weight, 1.478 / **1.612** / 1.504. So both 3-donor
    forms beat the 2-donor estimator -- the production key is STRONGER than what was measured,
    not weaker -- and the AGGRESSIVE extrapolation is the good one. Corroborated independently
    and on BOTH surfaces by `graded_fraction.py`: holding the key fixed and dialling the graded
    fraction from 0.21 down to 0.05 raises the multiplier on h1-test (1.101 -> 1.138) and orion
    (1.525 -> 1.543). Grading fewer genes is mildly better wherever it can be probed.

    TWO CAVEATS THAT TRAVEL WITH THAT. First, n = 1 SURFACE: orion has no legal third donor (its
    only candidate covers 12 of 132 targets), so rule 1's two-surface test cannot be run on this
    comparison and never will be with what we hold. What makes it actionable anyway is that it
    is not a harness-specific quantity -- it is the same three donor screens production ships,
    scored against a held-out reference. Second, the structural objection this function was
    written for is NOT refuted, only outvoted: under unanimity a 2-donor target still keeps
    0.545 of its genes graded against a 3-donor target's 0.106, so the donor-count ordering is
    still backwards. It simply does not cost score on the one surface that can see it. If a
    second surface ever contradicts, this is the fallback and it is already built and tested.

    It is also the same generalisation `donor_cosine_key` already makes (mean of the pairwise
    cosines), so the two keys are built the same way rather than one pairwise and one unanimous.

    Returns all-zero when fewer than two donors carry `pert`, exactly as `donor_agreement` does.
    """
    have = [n for n in model.names if str(pert) in model.pert_index[n]]
    n_genes = len(model.target_genes)
    if len(have) < 2:
        return np.zeros(n_genes, dtype=np.float64)
    vectors = []
    for name in have:
        vector = np.zeros(n_genes, dtype=np.float64)
        vector[model.tcols[name]] = model.s[name][
            model.pert_index[name][str(pert)]
        ][model.dcols[name]]
        vectors.append(model._renorm[name] * vector)
    total = np.zeros(n_genes, dtype=np.float64)
    pairs = 0
    for i in range(len(vectors)):
        for j in range(i + 1, len(vectors)):
            stack = np.vstack([vectors[i], vectors[j]])
            sign = np.sign(stack)
            unanimous = (np.abs(sign.sum(0)) == 2) & np.all(stack != 0, axis=0)
            # the same geometric mean `donor_agreement` uses, epsilon included, so the
            # two-donor case reproduces it bit for bit
            geo = np.exp(np.log(np.abs(stack) + 1e-300).mean(0))
            total += np.where(unanimous, geo, 0.0)
            pairs += 1
    return total / pairs


def _descending_midrank(scores: np.ndarray) -> np.ndarray:
    """0-based descending ranks, EQUAL SCORES SHARING THE MEAN of the positions they span.

    The obvious `argsort(argsort(-x))` gives tied entries consecutive distinct ranks decided by
    their position in the array, which is not a property of the key at all. That is not a
    cosmetic issue here: `donor_agreement` is exactly 0 on every gene the donors do not agree
    about -- a median of 16,545 of 18,533 genes per target -- so under `argsort` the great
    majority of a row's alphas are set by array index. Midranks make equal keys give equal
    alphas, which is the invariant the caller actually needs.
    """
    n = int(scores.size)
    order = np.argsort(-scores, kind="stable")
    ordered = scores[order]
    edges = np.flatnonzero(
        np.concatenate(([True], ordered[1:] != ordered[:-1], [True]))
    )
    starts, ends = edges[:-1], edges[1:]
    # mean of the consecutive integers start .. end-1
    group_rank = (starts + ends - 1) / 2.0
    rank = np.empty(n, dtype=np.float64)
    rank[order] = np.repeat(group_rank, ends - starts)
    return rank


def agreement_shrink_vector(
    agreement: np.ndarray,
    gamma: float,
    *,
    rank_scale: float = AGREEMENT_SHRINK_RANK_SCALE,
    protect: int | None = None,
) -> np.ndarray:
    """Per-gene amplitude factor in (0, 1], monotone increasing in cross-donor agreement.

        alpha_g = exp(-gamma * rank_g / rank_scale),  rank_g = 0 at the highest score

    `agreement` is one target's `MultiDonorModel.donor_agreement` row: the geometric mean
    of |delta| over the donors that carry the target, zeroed unless every one of them
    agrees on the sign. `rank_g` is that row's DESCENDING rank over the whole gene axis,
    ties broken by gene index with a stable full sort so the vector is a deterministic
    function of the scores rather than of a sort pivot -- which matters here, because the
    key is exactly 0 on every gene the donors do not agree about and those are usually the
    majority.

    WHY THIS KEY AND NOT A FITTED COEFFICIENT. Experiment 096 measured both, and they are
    not close. A per-gene regression SLOPE fitted on one cell line and applied to the other
    gains `1.02x`: the cross-cell-line optimal ridge is lambda ~ 100, which shrinks the
    per-gene coefficients ~99% of the way back to a single global scalar, against a 4-6.5x
    oracle. Keying the same shrinkage on `donor_agreement` gains `1.17x` (orion) to `1.38x`
    (h1-test) on `sqrt(W)`, positive on BOTH surfaces and at every setting swept, and it
    beats the self-referential `|effect|`-rank key everywhere too. The rule that generalises:
    a per-gene slope is a CELL-LINE property and does not transfer; a per-gene confidence
    score is a MODEL property and does. See MSE-MECHANISM.md section (c) for the slope
    measurement, intervention 5 for this one, and `mech_miller.py` for the sweep.

    THIS FUNCTION RETURNS ALPHA ONLY. What alpha interpolates *toward* is the caller's
    choice, and it is not a free one. Miller et al. (bioRxiv 2025.10.20.683304), and
    `mech_miller.py` after them, interpolate toward the PANEL MEAN,
    `dh' = alpha*dh + (1-alpha)*dh_bar`. The argument for it is that
    `expr_mse_unbiased_capped_norm`'s scored 0-anchor is the context's mean perturbation
    response assigned identically to every perturbation, not the control -- the board has
    both ends, team "Kun Lab"'s explicit zero-effect submission at raw 1.0424 against a
    scored baseline of raw ~0.9879.

    MEASURED, that argument does not survive, for two reasons, both in
    `experiments/096-mse-floor/mech_shrink_target.py`. First, the anchor the metric rewards
    is the REFERENCE's mean effect, and our own panel mean is not an estimate of it:
    cos is +0.016 on h1-test and -0.066 on orion, so the capture fraction is 0.0% / 0.4% and
    the two surfaces disagree in sign. Second, shrinking toward the panel mean preserves the
    common mode while collapsing the deviations, which inflates `(A - B)/A` from ~0.012 to
    ~0.055 at gamma=1 -- and `A - B` is precisely the term the master formula
    `MSE_min = 1 + (Kc - 2*sqrt(claim*W))/sumD` drops, so the published `sqrt(W)` multiplier
    flatters it. Scored through the scorer's own identity with `B` recomputed rather than
    assumed, `zero` beats `mean` on both harnesses at every gamma with both keys (12/12), and
    `mean` is WORSE THAN NOT SHRINKING at all on h1-test from gamma=1 upward.

    THE SHRINKAGE IS A FUNCTION OF THE INFORMATION IN THE KEY, and that is a correctness
    requirement rather than a nicety. Two guards enforce it, both added after a shipped-path
    bug that all of the measurement above was blind to:

      * A row with FEWER THAN TWO DISTINCT key values expresses no ordering, so it returns all
        ones. The all-zero row (fewer than two donors carry the target) is the common case --
        `mse-mechanism` counted **111 of 900 rows on the real challenge panel**, 12% of what we
        would ship -- but a constant NON-zero row has the identical defect and nobody had
        counted those. Before this, such a row returned `exp(-gamma * g / 2000)` indexed by
        POSITION IN THE GENE AXIS: gene 0 kept full amplitude and gene 18,000 was shrunk to
        nothing, for no reason but where it sat in the array. Worse, those rows were then
        scored as though the shrinking had helped, which is what inflated the family's
        estimated gain; removing the artefact took it from +0.0102 to +0.0055 over a12.
      * Within a row, ties share a MIDRANK (`_descending_midrank`), so the ~16,545 genes whose
        agreement is exactly 0 all receive the SAME alpha instead of an arbitrary ramp across
        them. This is the same defect at the other scale and it was not separately reported;
        it is fixed structurally rather than by special-casing the zero row, because a
        special case would leave the constant-row and the tied-block variants standing.

    `protect` forces one gene's alpha to exactly 1.0 and is used for each row's own target
    coordinate. `cell_eval2` excludes that coordinate from all six scored members, so this
    costs exactly nothing -- but without it we would ship a knockdown shrunk up to 9x at
    gamma=0.5 and 99x at gamma=1.0 on the targets whose own gene ranks poorly, a defect that
    passes every number we check.

    Substituting a self-referential key on an uninformative row would contradict the finding
    this whole intervention rests on -- the two scores computed from our own prediction
    correlate `-0.001/+0.253` and `-0.337/-0.102` with the true per-perturbation cosine across
    the two harnesses, and the second is actively harmful.

    AND THE HARNESSES UNDERSTATE THE PENALTY. The `(A - B)` term enters the identity as
    `claim * eps / sumD` with `eps = (A - B)/B`, so its price is set by `claim/B`. That ratio is
    0.183 on h1-test and 0.369 on orion (both recomputed here from this script's own terms) but
    **0.879 on the challenge panel at the #247 cap** (experiment 097's `exact.py`). A family
    that inflates `eps` is therefore charged 2.4x more on the board than orion reports and 4.8x
    more than h1-test reports -- and a12's design raises `claim` TO that cap, which raises the
    price further. Whatever gamma looked survivable on a harness is more expensive after a12.

    WHAT IS AND IS NOT ESTABLISHED. The ORDERING above is checked three times, by three
    implementations that could have disagreed and did disagree on levels: this agent's
    `mech_shrink_target.py`, experiment 096's `mech_shrink_exact.py`, and experiment 097's
    `exact.py` differ by up to 0.027 in `d(mse)` on the same harness and even on whether the
    mean target is worse than not shrinking, yet all three put `zero` ahead of `mean` on both
    surfaces at every gamma (12/12, then 4/4). The
    MAGNITUDES are not established, and both routes share one omission: `claim` is held at its
    UNSHRUNK value. `claim = sum_p w_p min(C_hat_p, K * C_hat_p^ref)` is a functional of the
    EMITTED cells, so shrinking moves it, and the optimum `s* = sqrt(claim/B')` sits on it. The
    two targets shrink by the same alpha and differ only in the centre, so the ordering should
    survive; the numbers may not. And this omission is shared by ALL THREE routes, so the
    three-way concordance does not test it -- rule 10, a concordance is only a check on what
    the routes could have disagreed about, and none of them re-emitted cells.

    DISCHARGED, by `experiments/097-direction/capmargin.py` -- real emission through the
    production `emit_counts`, real `jackknife_correction` on the emitted cells, each arm at its
    OWN solved cliff. Note the panel total `claim` could NOT have settled it: `claim` sits on the
    #247 cap on h1-test (85/85 rows pinned), and a saturated `claim` is set by the REFERENCE's
    jackknife, so a flat total is equally consistent with "nothing moved" and with "our margin is
    being eaten and has not run out". The per-row diagnostic is what answers it:

        harness   arm     rows at cap    margin C_hat/(K*C_p) min / p05 / med
        h1-test   base       85/85           1.223 / 1.541 / 1.65
        h1-test   STACK      85/85           1.222 / 1.535 / 1.65
        orion     base        0/132          0.173 / 0.388 / 0.47
        orion     STACK       0/132          0.173 / 0.389 / 0.47

    No row changes side and the whole distribution is unchanged to three decimals: `C_hat_p` is
    invariant to the shrinkage ROW BY ROW, not merely in the panel total. The two harnesses
    happen to bracket both regimes, and the informative one is orion, where 0/132 rows are capped
    so `claim` is free to respond to our emission -- it moved +0.04 %. The mechanism is that
    `C_hat_p` is set by the emission and by the profile's second moment, and every profile is a
    small perturbation of the same control composition, so `C_hat_p` is near-constant in `p` and
    the margin varies almost entirely through the reference's `C_p`. The most-shrunk rows are
    therefore among the SAFEST (h1-test low-z quartile minimum margin 1.560 against a panel
    minimum of 1.222), not the closest to the edge.

    And the challenge panel is in the SAFE regime by construction: a12 deliberately raises
    `claim` TO the #247 cap, i.e. into the saturated h1-test case where rows are pinned and
    `claim` is reference-set. The concern is smallest exactly where we ship.

    NOT REALISED THROUGH THE DE PATH. Everything above is `mse`-only accounting. The per-gene
    damping changes which genes clear the significance gate, so `fid`, `reach`, `jac` and
    `nmae` all move and none of them has been measured under it. Run it through
    `experiments/095-clean-harness/realize.py` before spending a slot.

    NOTE WHAT THIS DOES TO GLOBAL AMPLITUDE, AND RE-SOLVE IT. At `gamma=1.0` over 18,533 genes
    only `2000*ln2 ~ 1,386` genes keep more than half their amplitude, so the panel-wide effect
    norm falls sharply. The framework that measured the gain re-optimises the global scale AFTER
    shrinking (`s* = sqrt(claim/B')`) and reports only the direction gain, so an arm that ships
    this at the previous VCC_ALPHA is NOT the arm that was measured.

    MEASURED compensation for the full stack (per-gene `agreement` at gamma=1 toward zero, times
    per-perturbation `donor-cos` at gamma=2 floor 0.3), solved from the emitted cells rather than
    computed from a mean: `s` goes 0.5576 -> 1.6956 on h1-test and 0.8330 -> 2.3748 on orion,
    i.e. **~2.9-3.0x**, so a12's VCC_ALPHA 1.012 becomes ~3.0. Take the RATIO and not the level:
    the levels are harness-specific, but both arms were solved at the same emission so the ratio
    is the transferable part. Solving `s` from a realisation is also the one form that cannot be
    got wrong by taking the mean of the wrong power, which happened three times in this
    experiment's history.
    """
    scores = np.asarray(agreement, dtype=np.float64)
    if scores.ndim != 1:
        raise ValueError(
            f"agreement must be one value per gene, got shape {scores.shape}"
        )
    decay = float(gamma)
    if not np.isfinite(decay) or decay < 0.0:
        raise ValueError(f"--shrink-gamma must be finite and >= 0, got {gamma!r}")
    scale = float(rank_scale)
    if not np.isfinite(scale) or scale <= 0.0:
        raise ValueError(f"rank_scale must be finite and > 0, got {rank_scale!r}")
    n_genes = int(scores.size)
    scores = np.nan_to_num(scores)
    if protect is not None and not 0 <= int(protect) < n_genes:
        raise ValueError(f"protect must index a gene, got {protect!r} for {n_genes}")
    # A key with fewer than two distinct values carries no ordering to express.
    if n_genes == 0 or np.unique(scores).size < 2:
        return np.ones(n_genes, dtype=np.float64)
    alpha = np.exp(-decay * _descending_midrank(scores) / scale)
    if protect is not None:
        alpha[int(protect)] = 1.0
    return alpha


def donor_cosine_key(
    model: "MultiDonorModel",
    perts: Sequence[str],
) -> tuple[np.ndarray, np.ndarray]:
    """Mean pairwise cosine between the carrying donors' effect vectors, per target.

    Returns `(score, carried)`, both length `len(perts)`; `carried` is False where fewer than
    two donors carry the target and `score` is 0.0 there.

    Experiment 097 found this transfers, REALISED, at 1.246x (h1-test) / 1.176x (orion) as a
    per-perturbation weight, against 1.035x for experiment 096's best per-perturbation key,
    and it survives a permutation null over the same grid and argmax at p <= 0.001. It asks a
    different question from `donor_agreement`: agreement is a per-GENE magnitude gated on sign,
    this is a whole-vector DIRECTION agreement for the perturbation as a unit. Both are
    computed from donors agreeing with each other rather than from our own prediction, which is
    the property experiment 096 measured as the one that transfers.

    Built from exactly the vectors `donor_agreement` scores -- `_renorm[n] * s[n][p]` mapped to
    the target axis with zeros off it, i.e. the gene-specific DEVIATIONS after each donor's own
    mean response is removed, which is also what experiment 097's `wkit.donor_vec` and
    `perpert.donor_rows` build. The renormalisation cancels out of a cosine and is applied only
    so the two are literally the same vector.

    Zeros off a donor's gene axis are deliberate and are 097's construction: a gene one donor
    does not carry contributes nothing to the inner product rather than restricting the axis
    per pair, so all pairs are scored on one axis and the cosines are comparable across targets.
    """
    names = list(model.names)
    n_genes = len(model.target_genes)
    score = np.zeros(len(perts), dtype=np.float64)
    carried = np.zeros(len(perts), dtype=bool)
    for row, pert in enumerate(perts):
        have = [n for n in names if str(pert) in model.pert_index[n]]
        if len(have) < 2:
            continue
        vectors = []
        for name in have:
            vector = np.zeros(n_genes, dtype=np.float64)
            vector[model.tcols[name]] = model.s[name][
                model.pert_index[name][str(pert)]
            ][model.dcols[name]]
            vectors.append(model._renorm[name] * vector)
        cosines = []
        for i in range(len(vectors)):
            for j in range(i + 1, len(vectors)):
                denominator = float(
                    np.linalg.norm(vectors[i]) * np.linalg.norm(vectors[j])
                )
                cosines.append(
                    float(vectors[i] @ vectors[j]) / denominator
                    if denominator > 1e-12 else 0.0
                )
        if cosines:
            score[row] = float(np.mean(cosines))
            carried[row] = True
    return score, carried


def pert_shrink_vector(
    score: np.ndarray,
    carried: np.ndarray,
    gamma: float,
    floor: float = PERT_SHRINK_FLOOR,
    missing: str = PERT_SHRINK_MISSING,
) -> np.ndarray:
    """Per-perturbation weight `z_p = max(floor, rank01(score)^gamma)`, in (0, 1].

        rank01(x) = (argsort(argsort(x)) + 0.5) / n      ascending, so the best target -> ~1

    `rank01` is experiment 097's, verbatim. It is taken over the CARRIED targets only, so those
    span the full (0, 1) rather than being compressed by uncarried rows sitting at the bottom;
    097 ranked over all rows because on its harnesses nearly all were carried, and on the real
    300-target panel the difference is 6 rows, 2% of the rank range. That is a judgement call,
    recorded rather than hidden.

    A target fewer than two donors carry gets `z = 1`: no cross-donor evidence, no shrink, the
    same policy `agreement_shrink_vector` applies per gene. Note this makes those targets the
    HIGHEST-weighted rows, which is the conservative direction (they are emitted exactly as an
    unflagged build emits them) but is worth knowing. 097's own scripts instead drop them to
    the bottom of the rank scale; that is a self-referential penalty on targets we have no
    information about, and the same reasoning that rejects an `|effect|` fallback rejects it.

    THE FLOOR IS APPLIED TO THE RANK, INSIDE THE EXPONENT: `z = max(rank01, floor)**gamma`, so
    `--shrink-pert-floor 0.3` at gamma=2 leaves a floored target at 0.09 of full amplitude, not
    0.3. That is experiment 097's construction and its grid is reported in those units; applying
    the clamp after the exponent instead would silently disagree with every number they quote.

    `missing` decides what happens to a target fewer than two donors carry, whose cosine is
    UNDEFINED rather than zero:

      * `damp`   -- score it 0.0 and rank it with everything else. 097's recommendation. Note
                    what this actually does on the real panel: challenge cosines run
                    [-0.285, +0.439] (min/max over the 294 CARRIED targets; p05/p95 are
                    -0.043/+0.115) with median +0.0284, so 0.0 is NOT the bottom of the scale,
                    it is a mid-rank. On h1-test, where most cosines are positive, it does sit
                    low, which is why 097 describes it as damping.
      * `neutral`-- rank over the carried targets only and give the rest the MEDIAN carried
                    weight. The policy that follows from the principle alone without leaning on
                    097's n=4: absence of evidence is neither confidence nor its opposite.
                    Nobody has measured this one.
      * `keep`   -- z = 1. Measured at 1.023 (gamma=1) and 0.958 (gamma=2) on h1-test against
                    1.207 / 1.239 for `damp`, i.e. worse than not shrinking. Kept reachable
                    because that measurement is 4 perturbations on ONE surface.

    `gamma = 0` is an exact no-op under every policy.
    """
    scores = np.nan_to_num(np.asarray(score, dtype=np.float64))
    keep = np.asarray(carried, dtype=bool)
    if scores.ndim != 1:
        raise ValueError(f"score must be one value per target, got {scores.shape}")
    if keep.shape != scores.shape:
        raise ValueError(f"carried shape {keep.shape} != score shape {scores.shape}")
    decay = float(gamma)
    if not np.isfinite(decay) or decay < 0.0:
        raise ValueError(f"--shrink-pert-gamma must be finite and >= 0, got {gamma!r}")
    bottom = float(floor)
    if not np.isfinite(bottom) or not 0.0 <= bottom <= 1.0:
        raise ValueError(f"--shrink-pert-floor must lie in [0, 1], got {floor!r}")
    if missing not in ("damp", "neutral", "keep"):
        raise ValueError(
            f"--shrink-pert-missing must be damp, neutral or keep, got {missing!r}"
        )

    def rank01(values: np.ndarray) -> np.ndarray:
        """097's rank01, ascending, ties sharing a midrank."""
        descending = _descending_midrank(values)
        return (((values.size - 1) - descending) + 0.5) / values.size

    weight = np.ones(scores.size, dtype=np.float64)
    if missing == "damp":
        # score 0.0 for an uncarried target and rank the whole panel together
        ranked = np.where(keep, scores, 0.0)
        # Same structural guard as `agreement_shrink_vector`: fewer than two distinct values
        # is no ordering at all, so nothing is shrunk.
        if scores.size and np.unique(ranked).size >= 2:
            weight = np.maximum(rank01(ranked), bottom) ** decay
        return weight
    index = np.flatnonzero(keep)
    if index.size and np.unique(scores[index]).size >= 2:
        carried_weight = np.maximum(rank01(scores[index]), bottom) ** decay
        weight[index] = carried_weight
        if missing == "neutral":
            weight[~keep] = float(np.median(carried_weight))
    return weight


def emission_spread_vector(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    *,
    min_cpm: float = MIN_CPM,
    tail_spread: float | None = None,
    keep_top: int = 0,
    keep_top_frac: float | None = None,
    keep_spread: float | None = None,
    rank_score: np.ndarray | None = None,
    spread_map: str = "block",
    spread_lo: float = 0.0,
    spread_hi: float = 0.60,
) -> np.ndarray:
    """Per-gene `spread` for one target's emission.

    `emit_counts(mode="quasi")` splits each gene's panel total into an evenly-allocated
    part and a Poisson-scattered part, `spread` being the scattered fraction, and then
    repairs the total EXACTLY. Per-cell variance is ~`spread * T / n` against Poisson's
    `T / n`, while the gene's panel total -- and therefore the pseudobulk that `pds` and
    `expr_mse` are computed from -- is untouched. Raising a gene's `spread` raises its
    Wilcoxon p-value without touching its effect, which is the only lever experiment 092
    found that shrinks our over-large significant call set without paying `nmae`: nulling
    the tail's effect instead degrades `nmae` from 0.9595 to 1.0154 and wipes out half the
    gain.

    Experiment 093 measured this channel end to end on the HCT116 harness, whose geometry is
    production's exactly (400 cells/target, 18,400 controls, 18,533 genes). It works, but it
    SATURATES: `spread` is clipped to 1.0 inside `emit_counts` -- Poisson is the noisiest a
    total-preserving scatter can be -- and even there median `n_pred` only falls 5340 ->
    2317 against a reference median of 304, so all 132 targets still over-call and the
    counterfactual's 850 is out of reach. Two thirds of experiment 092's +0.0299 arrives
    anyway, because `fid` is `k / n_pred` once we over-call and the calls that leave first
    are the weakest ones.

    Recommended: `tail_spread = 1.0, keep_top = 850` -- +0.0202 +- 0.0002 avg6 over four
    seeds, against +0.0196 for `keep_top = 0`. The `keep_top` curve is a broad plateau
    (400 -> +0.0193, 600 -> +0.0205, 850 -> +0.0203, 1200 -> +0.0199, 2000 -> +0.0192).
    |predicted effect| is a weak proxy for our own p-value order in isolation, but it beats
    ranking by that p-value here: an oracle arm ranked on the production run's own p_adj
    scored +0.0193, BELOW the no-keep arm, because it spends its protection on genes the
    spread was never going to de-significate. See
    experiments/093-spread-truncation/RESULTS.md.

    Experiment 094 reopened this. Two board facts retire the v14/v15 reading above.

    First, `fid` rose +0.0179 (v14 -> v15) and then +0.0176 (v15 -> v11) as the call set
    grew, i.e. LINEARLY and with no sign of saturating. `fid`'s coverage term
    `min(1, n_pred/n_real)` saturates the moment `n_pred >= n_real`, so a still-linear
    response is direct evidence that **v11 itself still undershoots `n_real`** on the
    challenge contexts. Every arm in 093 moved `n_pred` the wrong way.

    Second, and separately, the scorer's `reach` walks the reference's significant genes in
    OUR confidence order -- significant calls first, then ascending predicted adjusted p.
    That order is a deterministic function of `|F_bar_g - 1/2| / sqrt(tc_g)`, the mean of
    the control CDF over our emitted count law, which correlates with `|predicted effect|`
    at Spearman +0.07. So our ranking is very nearly independent of our own model's signal,
    and is instead dominated by each gene's control-distribution skewness: a minimum-variance
    emitter puts 400 cells near the predicted mean, which for right-skewed counts sits above
    the median, so the rank-sum reports "significantly up" at zero true effect.

    `keep_spread` decouples the two ends of the channel, which were previously locked to
    `tail_spread` and the constant `SPREAD`. Setting `keep_spread < tail_spread <= SPREAD`
    LOWERS spread on the genes carrying real predicted effect rather than raising it on the
    tail: `n_pred` goes UP (paying the coverage term that v14 and v15 both underpaid) while
    the informative genes still gain rank against the artefact calls. Both effects are free
    on `pds`, `expr_mse` and `nmae`, which the exact-total repair leaves bit-identical --
    board-verified in v14, where `pds` moved by exactly 0.000000.

    `rank_score` replaces `|effect|` as the ranking key. Experiment 095 measured
    `MultiDonorModel.donor_agreement` -- geometric magnitude gated on cross-donor sign
    unanimity -- at **+0.00678** overall as a full ordering against `|lfc|`'s **+0.00448**,
    the best buildable key found. It does NOT help under `spread_map="block"`, where it
    scores +0.00341 against `|effect|`'s +0.00583, because a block promotion preserves the
    emitted p-order *inside* the block and agreement selects a set whose internal order is
    worse. Pair it with `spread_map="graded"`.

    `spread_map="graded"` assigns spread by rank across the whole gated axis, linearly from
    `spread_lo` at the most confident gene to `spread_hi` at the least, which lets the
    realised p-value order follow the key rather than merely putting one block in front.
    That is the mechanism the agreement key needs. **Its effect on `n_pred` is unmeasured** --
    the block map was probed at +3.9% on the real challenge controls, and a graded map has no
    such measurement, so probe it with `experiments/094-reach-ordering/probe_de.py` before
    spending a submission on it.

    With `tail_spread` unset this returns the shipped constant vector, so an unflagged
    build is bit-identical to v11.
    """
    n_genes = int(np.asarray(effect).size)
    if tail_spread is None:
        return np.full(n_genes, SPREAD)
    tail = float(tail_spread)
    if not np.isfinite(tail) or not 0.0 <= tail <= 1.0:
        raise ValueError(f"--tail-spread must lie in [0, 1], got {tail_spread!r}")
    if keep_top and keep_top_frac is not None:
        raise ValueError("pass at most one of --keep-top and --keep-top-frac")
    kept = SPREAD if keep_spread is None else float(keep_spread)
    if not np.isfinite(kept) or not 0.0 <= kept <= 1.0:
        raise ValueError(f"--keep-spread must lie in [0, 1], got {keep_spread!r}")
    if keep_spread is not None and not (keep_top or keep_top_frac is not None):
        raise ValueError("--keep-spread requires --keep-top or --keep-top-frac")

    gated, score = emission_rank_axis(
        effect, b_ctrl, min_cpm=min_cpm, rank_score=rank_score
    )
    keep = resolve_keep_count(gated, keep_top=keep_top, keep_top_frac=keep_top_frac)

    if spread_map not in ("block", "graded"):
        raise ValueError(f"--spread-map must be 'block' or 'graded', got {spread_map!r}")

    if spread_map == "graded":
        lo, hi = float(spread_lo), float(spread_hi)
        for nm, v in (("--spread-lo", lo), ("--spread-hi", hi)):
            if not np.isfinite(v) or not 0.0 <= v <= 1.0:
                raise ValueError(f"{nm} must lie in [0, 1], got {v!r}")
        if lo > hi:
            raise ValueError("--spread-lo must not exceed --spread-hi")
        vector = np.full(n_genes, tail)
        gi = np.flatnonzero(gated)
        if gi.size:
            order = np.lexsort((gi, -score[gi]))
            rank = np.empty(gi.size, dtype=np.float64)
            rank[order] = np.arange(gi.size, dtype=np.float64)
            frac = rank / max(gi.size - 1, 1)
            vector[gi] = lo + (hi - lo) * frac
        return vector

    vector = np.full(n_genes, tail)
    vector[keep_block_indices(gated, score, keep)] = kept
    return vector


def emission_rank_axis(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    *,
    min_cpm: float = MIN_CPM,
    rank_score: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """The gated DE axis and the per-gene score `--keep-top` competes on.

    The DE universe: the scorer filters at 5 CPM on the reference control, and our own
    gate zeroes the effect outside `min_cpm`. A gene outside it can never be called, so
    its spread is irrelevant to every metric; the fraction form of `--keep-top` is taken
    over the gated genes because that is the population `keep_top` is competing with.
    """
    n_genes = int(np.asarray(effect).size)
    expression = np.expm1(b_ctrl)
    cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
    gated = (cpm > min_cpm) & (b_ctrl > 0)
    score = np.where(
        gated,
        np.abs(np.nan_to_num(np.asarray(effect, dtype=np.float64)))
        if rank_score is None
        else np.nan_to_num(np.asarray(rank_score, dtype=np.float64)),
        -np.inf,
    )
    if rank_score is not None and np.asarray(rank_score).size != n_genes:
        raise ValueError("rank_score must be one value per gene")
    return gated, score


def resolve_keep_count(
    gated: np.ndarray,
    *,
    keep_top: int = 0,
    keep_top_frac: float | None = None,
) -> int:
    """`--keep-top` / `--keep-top-frac` resolved against the gated axis."""
    if keep_top_frac is not None:
        frac = float(keep_top_frac)
        if not np.isfinite(frac) or not 0.0 <= frac <= 1.0:
            raise ValueError(f"--keep-top-frac must lie in [0, 1], got {keep_top_frac!r}")
        return int(round(frac * int(gated.sum())))
    keep = int(keep_top)
    if keep < 0:
        raise ValueError(f"--keep-top must be >= 0, got {keep_top!r}")
    return keep


def keep_block_indices(gated: np.ndarray, score: np.ndarray, keep: int) -> np.ndarray:
    """Positions of the top-`keep` gated genes, in the shipped stable order.

    Factored out of `emission_spread_vector` so the resampling emitter's `param_mask` can
    be the SAME set the block spread map promotes, rather than a second ranking that
    happens to agree. With `--resample-keep-top` equal to `--keep-top` the parametric half
    of `emit_counts(mode="resample_mix")` is then exactly a15's keep block and only the
    tail's emission law changes (experiment 177 section 4).
    """
    n_genes = int(np.asarray(score).size)
    # Capped at the gated count: a gene the gate zeroed can never be called, so promoting
    # one would spend a keep slot on a gene no metric can see.
    keep = min(int(keep), int(np.asarray(gated).sum()))
    if keep <= 0:
        return np.empty(0, dtype=np.int64)
    # Ranked with a stable full sort so the chosen set is a deterministic function of
    # the scores and not of an argpartition pivot.
    return np.lexsort((np.arange(n_genes), -np.asarray(score)))[:keep]


def emission_param_mask(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    *,
    min_cpm: float = MIN_CPM,
    resample_keep_top: int,
    rank_score: np.ndarray | None = None,
) -> np.ndarray:
    """`param_mask` for `emit_counts(mode="resample_mix")`: True = parametric quasi law.

    Experiment 177 measured the whole family on h1-test at a verified-identical delivered
    pseudobulk. The shipped `spread 1.0 / phi 4` tail calls 38-45% of the DE universe
    significant at ZERO predicted effect and ~100% of those DOWN; the resampled law calls
    0.0-0.2% (1 / 18 / 19 genes on challenge A/B/C). Replacing a15's tail block with
    resampled real control cells realises 69% of INBOX-4's full |H| reordering
    counterfactual, +0.0352 raw `reach` on h1-test, at `n_pred` 6066 -> 3304.
    """
    gated, score = emission_rank_axis(
        effect, b_ctrl, min_cpm=min_cpm, rank_score=rank_score
    )
    mask = np.zeros(int(np.asarray(effect).size), dtype=bool)
    mask[keep_block_indices(gated, score, int(resample_keep_top))] = True
    return mask


def rc_diagnostics_summary(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    """Median / worst of `rc_hybrid`'s per-block diagnostics, for the meta.

    These are the numbers 185 section 2d reports for the two coupling constraints, and they
    are the cheapest evidence that the weld and the repair did what they claim: `forced`,
    `max_abs_resid`, `moves` and `unresolved` all zero means the greedy placement balanced
    every row on its own and the sub-gate slack was never spent, and
    `z_realised_vs_designed_max` is how far the delivered Wilcoxon z sits from the designed
    one (the `T + 1` remainder cells are the only source of that difference).
    """
    keys = sorted({key for row in rows for key in row})
    summary: dict[str, object] = {"blocks": len(rows)}
    for key in keys:
        values = np.asarray([float(row[key]) for row in rows if key in row])
        summary[key] = {
            "median": float(np.median(values)),
            "max": float(values.max()),
            "min": float(values.min()),
        }
    return summary


def assert_raw_control_counts(
    control_counts: np.ndarray,
    lib_sizes: np.ndarray,
    *,
    where: str,
) -> None:
    """Fail closed on the confound experiment 177 section 0.1 diagnosed.

    `emit_counts(mode="resample_mix")` repairs each gene's panel total EXACTLY to
    `floor(n * L * profile)`, `L = int(median(lib_sizes))`. If the donated control cells
    are rescaled to some nominal depth while `lib_sizes` still carries the real one, the
    repair has to peel the difference back out one count per cell per round -- a
    maximally variance-DESTROYING operation that called 10,658 of 10,779 genes
    significant in 177's reproduction and inverts the sign of the whole arm.

    So the call-site contract is: RAW counts, the cells' OWN library sizes, no rescaling.
    This asserts it on the block actually handed to the emitter. It reads differently when
    the thing it tests is absent -- a block rescaled to a nominal 100k depth fails the row
    -sum equality by a factor of ~5, and a CPM-normalised block fails integrality.
    """
    block = np.asarray(control_counts, dtype=np.float64)
    libs = np.asarray(lib_sizes, dtype=np.float64)
    if block.ndim != 2 or block.shape[0] != libs.size:
        raise ValueError(f"{where}: control block {block.shape} does not match {libs.size} libraries")
    if not np.all(np.isfinite(block)):
        raise ValueError(f"{where}: control counts are not finite")
    if np.any(block < 0.0):
        raise ValueError(f"{where}: control counts contain negatives")
    if not np.array_equal(block, np.floor(block)):
        raise ValueError(
            f"{where}: control counts are not raw integer counts -- a normalised or "
            "rescaled control block breaks the exact-total repair (177 section 0.1)"
        )
    rows = block.sum(axis=1)
    if not np.array_equal(rows, libs):
        worst = float(np.max(np.abs(rows - libs)))
        raise ValueError(
            f"{where}: resampled control cells were rescaled -- row sums differ from "
            f"their own library sizes by up to {worst:g}. Pass RAW counts with the "
            "cells' own lib_sizes (177 section 0.1)"
        )


@contextmanager
def atomic_h5_output(output: Path):
    """Yield a new HDF5 handle and publish it only after a clean close."""
    temporary = output.with_name(f".{output.name}.partial-{os.getpid()}")
    if temporary.exists():
        raise FileExistsError(f"refusing to overwrite temporary file {temporary}")
    try:
        with h5py.File(temporary, "x") as handle:
            yield handle
        os.replace(temporary, output)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def atomic_write_text(path: Path, text: str) -> None:
    temporary = path.with_name(f".{path.name}.partial-{os.getpid()}")
    if temporary.exists():
        raise FileExistsError(f"refusing to overwrite temporary file {temporary}")
    try:
        temporary.write_text(text)
        os.replace(temporary, path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


# --- VCC_ALPHA=knee -------------------------------------------------------------------
# 170's `chal.py` measured `A`, `B`, `claim` and `rho` on the challenge panel with exactly
# these primitives and exactly these constants; the numbers below have to reproduce its
# table at the shipped alpha, which is this feature's known-answer row.
KNEE_SUBSAMPLE = 24            # targets per context used for `Chat`, 170's chal.py value
KNEE_JACKKNIFE_CHUNK = 64      # chunk-invariant to ~1 ULP; a cost knob, not an answer knob
KNEE_CHAT_SEED = 3000          # 170's chal.py seed base, so the two agree cell for cell


@dataclass(frozen=True)
class ControlFrame:
    """One context's control side: everything the knee needs that is not our prediction."""

    context: str
    b_ctrl: np.ndarray
    median_library: float
    cp: float
    cp_values: np.ndarray
    cp_source: str


def ntc_construct_cp(
    control_matrix, ntc_ids: np.ndarray, *, chunk: int = KNEE_JACKKNIFE_CHUNK
) -> np.ndarray:
    """096's `C_p` proxy: the delete-1 jackknife of each ntc construct's own cells.

    A reference perturbation group is 400 cells of one construct at ~20k median depth, and
    an ntc construct is 400 cells of one construct at ~20k median depth, so the two earn the
    same order of sampling correction. It is a mild UNDER-estimate (real perturbation groups
    are slightly more dispersed than controls), which makes `min(Chat, C_p)` an
    under-estimate too -- conservative in the direction that matters.

    096 quoted 14.544 / 14.772 / 16.372 on A / B / C from 16 constructs each. This uses ALL
    of them (46 on the current bundle) because there is no principled "first 16" rule to
    carry to a bundle nobody has seen; `cpproxy180.py` measures the two against each other.
    """
    from cell_eval2.moments import jackknife_correction

    ids = np.asarray(ntc_ids).astype(str)
    values = []
    for name in dict.fromkeys(ids.tolist()):
        rows = np.flatnonzero(ids == name)
        block = sp.csr_matrix(control_matrix[rows])
        values.append(
            float(
                jackknife_correction(
                    block, np.zeros(rows.size, np.intp), 1, TS, chunk=chunk
                )[0]
            )
        )
    return np.asarray(values, dtype=np.float64)


def random_draw_cp(
    control_matrix,
    *,
    n_cells: int,
    n_draw: int = 12,
    seed: int = 0,
    chunk: int = KNEE_JACKKNIFE_CHUNK,
) -> np.ndarray:
    """096's `rand` proxy -- the fallback when a bundle ships no `ntc_id` column."""
    from cell_eval2.moments import jackknife_correction

    rng = np.random.default_rng(seed)
    rows_total = int(control_matrix.shape[0])
    take = min(n_cells, rows_total)
    values = []
    for _ in range(n_draw):
        rows = np.sort(rng.choice(rows_total, take, replace=False))
        block = sp.csr_matrix(control_matrix[rows])
        values.append(
            float(
                jackknife_correction(block, np.zeros(take, np.intp), 1, TS, chunk=chunk)[0]
            )
        )
    return np.asarray(values, dtype=np.float64)


def load_control_frame(
    bundle: Path, context: str, genes: Sequence[str], *, cells_per_pert: int
) -> ControlFrame:
    """`b_ctrl`, the median control library and `C_p`, from one context's own control file."""
    controls = ad.read_h5ad(bundle / f"context_{context}.h5ad")
    matrix = controls.X.tocsr() if not sp.isspmatrix_csr(controls.X) else controls.X
    if list(map(str, controls.var_names)) != list(map(str, genes)):
        raise ValueError(f"context {context} gene order differs from official")
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(np.float64)
    b_ctrl = np.log1p(TS * totals / totals.sum())
    median_library = float(np.median(np.asarray(matrix.sum(axis=1)).ravel()))
    if "ntc_id" in controls.obs.columns:
        values = ntc_construct_cp(matrix, controls.obs["ntc_id"].to_numpy())
        source = f"ntc_construct_jackknife[{values.size}]"
    else:
        values = random_draw_cp(matrix, n_cells=cells_per_pert)
        source = f"random_{cells_per_pert}_cell_draw[{values.size}]"
    del controls, matrix, totals
    gc.collect()
    return ControlFrame(
        context=str(context),
        b_ctrl=b_ctrl,
        median_library=median_library,
        cp=float(values.mean()),
        cp_values=values,
        cp_source=source,
    )


def knee_context_terms(
    effects: np.ndarray,
    frame: ControlFrame,
    cols: Mapping[int, int],
    *,
    min_cpm: float,
    effect_scale: float,
    shrink_scales: np.ndarray | None,
    shrink_center: np.ndarray | None,
    rank_scores: np.ndarray | None,
    spread_kwargs: Mapping[str, object],
    overdispersion: float,
    cells_per_pert: int,
    subsample: int = KNEE_SUBSAMPLE,
) -> dict[str, float]:
    """`A`, `B`, `Chat`, `claim` and `rho` for one context, in the scorer's own space.

    `pm = log1p(TS * profile)` is the pseudobulk the scorer reads: the quasi emitter's
    exact-total repair makes the delivered per-gene panel totals proportional to `profile`
    (179 gate 4 verified that on two real 19 GB artifacts, 900/900 blocks), so this is the
    delivered object and not the intended one -- which is the distinction that invalidated
    every effect-level screen before a23.
    """
    from cell_eval2.metrics.delta import PRED_TRACE_CAP_K, _across_pert_budget, _row_weights
    from cell_eval2.moments import jackknife_correction

    n_targets, n_genes = effects.shape
    b_ctrl = frame.b_ctrl
    profiles = np.empty_like(effects)
    for index in range(n_targets):
        profiles[index] = profile_from_effect(
            effects[index],
            b_ctrl,
            min_cpm=min_cpm,
            effect_scale=effect_scale,
            gene_scale=None if shrink_scales is None else shrink_scales[index],
            shrink_center=shrink_center,
        )
    predicted = np.log1p(TS * profiles)
    delta = predicted - b_ctrl[None, :]
    weights = _row_weights(n_targets, n_genes, dict(cols))
    per_row = np.einsum("ij,ij->i", delta, delta)
    for row, column in cols.items():
        per_row[row] -= delta[row, column] ** 2
    amplitude = float((per_row * weights).sum())
    budget = float(_across_pert_budget(predicted, dict(cols), weights))
    del predicted, delta

    libraries = np.full(cells_per_pert, frame.median_library)
    chosen = np.linspace(0, n_targets - 1, subsample).astype(int)
    chats = []
    for index in chosen:
        spread = emission_spread_vector(
            effects[index],
            b_ctrl,
            min_cpm=min_cpm,
            rank_score=None if rank_scores is None else rank_scores[index],
            **spread_kwargs,
        )
        values, columns, pointers = emit_counts(
            profiles[index],
            libraries,
            DISPERSION,
            np.random.default_rng(KNEE_CHAT_SEED + int(index)),
            mode="quasi",
            overdispersion=overdispersion,
            spread=spread,
        )
        block = sp.csr_matrix(
            (
                values.astype(np.float64),
                columns.astype(np.int64),
                pointers.astype(np.int64),
            ),
            shape=(cells_per_pert, n_genes),
        )
        chats.append(
            float(
                jackknife_correction(
                    block,
                    np.zeros(cells_per_pert, np.intp),
                    1,
                    TS,
                    chunk=KNEE_JACKKNIFE_CHUNK,
                )[0]
            )
        )
    del profiles
    chats = np.asarray(chats, dtype=np.float64)
    cap_per_row = float(PRED_TRACE_CAP_K) * frame.cp
    weight_sum = float(weights.sum())
    claim = float(np.mean(np.minimum(chats, cap_per_row)) * weight_sum)
    cap = float(cap_per_row * weight_sum)
    return {
        "A": amplitude,
        "B": budget,
        "claim": claim,
        "claim_cap": cap,
        "claim_fraction_of_cap": (claim / cap) if cap else 0.0,
        "C_hat_mean": float(chats.mean()),
        "C_hat_min": float(chats.min()),
        # the per-row values, so `claim` can be recomputed under a different `C_p` estimate
        # without re-emitting: `claim = mean(min(C_hat, K * C_p)) * weight_sum`
        "C_hat_values": [float(value) for value in chats],
        "weight_sum": weight_sum,
        "C_p": frame.cp,
        "C_p_source": frame.cp_source,
        "C_p_sd": float(frame.cp_values.std()),
        "cap_saturated_rows": int(np.sum(chats >= cap_per_row)),
        "cap_rows_measured": int(chats.size),
        # `rho` is the scorer's own clamped ratio and it reads 1.0 both ON the knee and
        # ABOVE it, so "rho == 1 after the correction" cannot tell the two apart. The
        # UNCLAMPED ratio is the quantity that has to land at 1.000: below it we are
        # forfeiting refundable amplitude, above it we are paying for amplitude the metric
        # will not refund. Every check in this feature reads `rho_unclamped`.
        "rho": float(min(1.0, budget / claim)) if claim else 1.0,
        "rho_unclamped": float(budget / claim) if claim else float("inf"),
        "median_control_library": frame.median_library,
    }


def solve_knee_alphas(
    spec: ChallengeSpec,
    effects: np.ndarray,
    frames: Mapping[str, ControlFrame],
    *,
    alpha_ref: float,
    **terms_kwargs: object,
) -> dict[str, dict[str, float]]:
    """One `alpha` per context, at that context's `rho = 1` knee."""
    cols = panel_gene_columns(spec)
    return solve_knee_alphas_with_columns(
        spec, effects, frames, cols, alpha_ref=alpha_ref, **terms_kwargs
    )


def panel_gene_columns(spec: ChallengeSpec) -> dict[int, int]:
    """`{row: gene column}` for the targets that resolve on the official gene axis.

    The metric never scores a perturbation's own knocked-down gene (#172), so these are the
    columns both `A` and the #348 budget drop. Identical to 170's `chal.py`.
    """
    positions = {str(gene): index for index, gene in enumerate(spec.genes)}
    return {
        row: positions[str(target)]
        for row, target in enumerate(spec.perts)
        if str(target) in positions
    }


def solve_knee_alphas_with_columns(
    spec: ChallengeSpec,
    effects: np.ndarray,
    frames: Mapping[str, ControlFrame],
    cols: Mapping[int, int],
    *,
    alpha_ref: float,
    **terms_kwargs: object,
) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    for context in map(str, spec.contexts):
        measured = knee_context_terms(effects, frames[context], cols, **terms_kwargs)
        if measured["B"] <= 0.0 or measured["claim"] <= 0.0:
            raise ValueError(
                f"context {context}: degenerate knee terms B={measured['B']} "
                f"claim={measured['claim']}"
            )
        scale = float(np.sqrt(measured["claim"] / measured["B"]))
        measured["alpha_ref"] = float(alpha_ref)
        measured["knee_scale"] = scale
        measured["alpha"] = float(alpha_ref) * scale
        out[context] = measured
        print(
            f"[knee] {context}: A {measured['A']:.6f}  B {measured['B']:.6f}  "
            f"claim {measured['claim']:.6f} ({100 * measured['claim_fraction_of_cap']:.2f}% "
            f"of cap)  C_p {measured['C_p']:.3f} [{measured['C_p_source']}]  "
            f"rho {measured['rho']:.4f}  ->  alpha {measured['alpha']:.6f}",
            flush=True,
        )
    return out


def write_submission_from_effects(
    spec: ChallengeSpec,
    effects: np.ndarray,
    bundle: Path,
    output: Path,
    *,
    max_rss_gib: float,
    min_cpm: float = MIN_CPM,
    effect_scale: float = 1.0,
    shrink_scales: np.ndarray | None = None,
    shrink_center: np.ndarray | None = None,
    tail_spread: float | None = None,
    keep_top: int = 0,
    keep_top_frac: float | None = None,
    keep_spread: float | None = None,
    rank_scores: np.ndarray | None = None,
    spread_map: str = "block",
    spread_lo: float = 0.0,
    spread_hi: float = 0.60,
    overdispersion: float = 0.0,
    emitter: str = "quasi",
    resample_keep_top: int | None = None,
    context_effects: Mapping[str, np.ndarray] | None = None,
) -> dict[str, object]:
    """Stream every context's raw counts with the exact v9 RNG and emission semantics.

    `context_effects` is `VCC_ALPHA=knee`'s only entry point into the writer: one effect
    matrix per context instead of one shared matrix. `None` (the default) uses `effects` for
    every context and is bit-identical to every prior build.
    """
    if emitter not in EMITTERS:
        raise ValueError(f"--emitter must be one of {EMITTERS}, got {emitter!r}")
    if emitter != "resample_mix" and resample_keep_top is not None:
        raise ValueError("--resample-keep-top requires --emitter resample_mix")
    if emitter in RC_EMITTERS:
        # The key the order is delivered on is the incumbent's own GATED effect, which is
        # the object `profile_from_effect` renders. Every knob that would make the rendered
        # effect something other than `where(gate, effect, 0)` is refused rather than
        # silently keyed around -- an ordering emitter pointed at the wrong vector would
        # look exactly as reasonable as one pointed at the right one (rule 10's corollary).
        if effect_scale != 1.0 or shrink_scales is not None or shrink_center is not None:
            raise ValueError(
                f"--emitter {emitter} needs the unshrunk effect: --effect-scale 1.0 and no "
                "--shrink-key / --shrink-pert-key"
            )
        if MASS_FIX is not None:
            raise ValueError(f"--emitter {emitter} does not support VCC_MASS_FIX")
        if RESAMPLE_X is not None:
            raise ValueError(f"--emitter {emitter} and VCC_RESAMPLE_X are two emission laws")
    if emitter == "resample_mix":
        if RESAMPLE_X is not None:
            raise ValueError(
                "VCC_RESAMPLE_X and --emitter resample_mix are two spellings of the same "
                "channel; set exactly one"
            )
        resample_keep_top = keep_top if resample_keep_top is None else int(resample_keep_top)
        if resample_keep_top < 0:
            raise ValueError(f"--resample-keep-top must be >= 0, got {resample_keep_top!r}")
    effects = np.asarray(effects, dtype=np.float64)
    if effects.shape != (len(spec.perts), len(spec.genes)):
        raise ValueError(
            f"effect shape {effects.shape} != {(len(spec.perts), len(spec.genes))}"
        )
    if len(spec.genes) != N_GENES:
        raise ValueError(f"official gene dimension is {len(spec.genes)}, not {N_GENES}")
    if context_effects is not None:
        missing = [c for c in map(str, spec.contexts) if c not in context_effects]
        extra = [c for c in context_effects if c not in set(map(str, spec.contexts))]
        if missing or extra:
            raise ValueError(
                f"context_effects must cover exactly {list(map(str, spec.contexts))}; "
                f"missing={missing}, unexpected={extra}"
            )
        for name, matrix in context_effects.items():
            if np.asarray(matrix).shape != effects.shape:
                raise ValueError(
                    f"context_effects[{name!r}] shape {np.asarray(matrix).shape} != "
                    f"{effects.shape}"
                )
    if shrink_scales is not None:
        shrink_scales = np.asarray(shrink_scales, dtype=np.float64)
        if shrink_scales.shape != effects.shape:
            raise ValueError(
                f"shrink_scales shape {shrink_scales.shape} != {effects.shape}"
            )
    if shrink_center is not None:
        shrink_center = np.asarray(shrink_center, dtype=np.float64)
        if shrink_center.shape != (effects.shape[1],):
            raise ValueError(
                f"shrink_center shape {shrink_center.shape} != {(effects.shape[1],)}"
            )
        if shrink_scales is None:
            raise ValueError("shrink_center requires shrink_scales")

    output = output.expanduser().resolve()
    metadata_path = Path(str(output) + ".meta.json")
    if output.exists() or metadata_path.exists():
        raise FileExistsError(f"refusing to overwrite {output} or {metadata_path}")
    if output.suffix != ".h5ad":
        raise ValueError(f"challenge output must end in .h5ad, got {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    pick_rng, emission_rng = make_rngs(SEED)
    n_rows = spec.n_rows
    obs_pert: list[str] = []
    obs_context: list[str] = []
    indptr = np.zeros(n_rows + 1, dtype=np.int64)
    nnz = 0
    row = 0
    context_audit: dict[str, object] = {}
    rc_audit: list[dict[str, object]] = []

    with atomic_h5_output(output) as handle:
        matrix_group = handle.create_group("X")
        matrix_group.attrs["encoding-type"] = "csr_matrix"
        matrix_group.attrs["encoding-version"] = "0.1.0"
        matrix_group.attrs["shape"] = np.asarray([n_rows, N_GENES], dtype=np.int64)
        data_dataset = matrix_group.create_dataset(
            "data",
            shape=(0,),
            maxshape=(None,),
            dtype=np.float32,
            chunks=(1 << 20,),
        )
        index_dataset = matrix_group.create_dataset(
            "indices",
            shape=(0,),
            maxshape=(None,),
            dtype=np.int32,
            chunks=(1 << 20,),
        )

        for context in spec.contexts:
            # `VCC_ALPHA=knee` hands the writer one effect matrix per context; every other
            # configuration hands it one shared matrix and this is `effects` itself.
            context_effect = (
                effects
                if context_effects is None
                else np.asarray(context_effects[str(context)], dtype=np.float64)
            )
            control_path = bundle / f"context_{context}.h5ad"
            controls = ad.read_h5ad(control_path)
            control_matrix = (
                controls.X.tocsr() if not sp.isspmatrix_csr(controls.X) else controls.X
            )
            if list(map(str, controls.var_names)) != list(map(str, spec.genes)):
                raise ValueError(f"context {context} gene order differs from official")
            totals = np.asarray(control_matrix.sum(axis=0)).ravel().astype(np.float64)
            b_ctrl = np.log1p(TS * totals / totals.sum())
            libraries = np.asarray(control_matrix.sum(axis=1)).ravel()
            control_rows = int(control_matrix.shape[0])
            median_library = float(np.median(libraries))
            # The resampling emitter draws real control cells, so it needs the count matrix
            # itself rather than only its column totals and row sums.
            #
            # `totals / totals.sum()` is the POOLED control profile, and it is exactly the
            # space the rendered profile lives in: `b_ctrl = log1p(TS * totals/sum)`, so
            # `expm1(b_ctrl) / sum(expm1(b_ctrl))` IS this vector. 177 section 7 (b): the
            # mean-per-cell-CPM profile differs from it by an sd of 0.020 in log2 and makes
            # `fold != 1` at zero predicted effect, which breaks the zero-effect null the
            # whole arm rests on.
            resampling = RESAMPLE_X is not None or emitter == "resample_mix"
            resample_ctrl = control_matrix if resampling else None
            resample_profile = (totals / totals.sum()) if resampling else None
            # `rc_hybrid` needs the control cells too, but only to build the DE-space block
            # the Wilcoxon reads -- the sorted log1p(CPM) rows over the SCORER's own 5-CPM
            # gate -- so the count matrix is released as soon as that is built.
            rc_control = None
            if emitter in RC_EMITTERS:
                rc_control = rc_emit.control_block_from_counts(
                    control_matrix, min_cpm=RC_SCORER_MIN_CPM
                )
                print(
                    f"[rc] {context}: {rc_control.n_controls} controls, "
                    f"{int(rc_control.keep.sum())}/{N_GENES} genes past the "
                    f"{RC_SCORER_MIN_CPM}-CPM scorer gate",
                    flush=True,
                )
                memory_checkpoint(f"rc control block {context}", max_rss_gib)
            del controls, totals
            if resample_ctrl is None:
                del control_matrix
                gc.collect()

            context_start_nnz = nnz
            for target_index, target_value in enumerate(spec.perts):
                target = str(target_value)
                profile = profile_from_effect(
                    context_effect[target_index],
                    b_ctrl,
                    min_cpm=min_cpm,
                    effect_scale=effect_scale,
                    gene_scale=(
                        None if shrink_scales is None else shrink_scales[target_index]
                    ),
                    shrink_center=shrink_center,
                )
                picks = pick_rng.choice(
                    len(libraries),
                    spec.cells_per_pert,
                    replace=len(libraries) < spec.cells_per_pert,
                )
                emit_kwargs = {}
                emit_mode = "quasi"
                if emitter == "resample_mix":
                    # a15's keep block stays on the parametric quasi law and a15's tail is
                    # replaced by resampled real control cells carrying the predicted fold
                    # change. `param_mask` is the SAME index set `emission_spread_vector`
                    # promotes when `resample_keep_top == keep_top`, so at K = 8000 the
                    # parametric half is exactly the incumbent's keep block and the only
                    # thing that moved is the tail's emission law.
                    pmask = emission_param_mask(
                        context_effect[target_index],
                        b_ctrl,
                        min_cpm=min_cpm,
                        resample_keep_top=int(resample_keep_top),
                        rank_score=(
                            None if rank_scores is None else rank_scores[target_index]
                        ),
                    )
                    # RAW counts with the cells' OWN library sizes -- never rescaled to a
                    # nominal depth. `assert_raw_control_counts` is the guard; 177 section
                    # 0.1 is why it exists.
                    control_block = np.asarray(
                        resample_ctrl[picks].todense(), dtype=np.float64
                    )
                    assert_raw_control_counts(
                        control_block, libraries[picks], where=f"{context}/{target}"
                    )
                    emit_mode = "resample_mix"
                    emit_kwargs = {
                        "control_counts": control_block,
                        "control_profile": resample_profile,
                        "param_mask": pmask,
                    }
                elif RESAMPLE_X is not None:
                    # Top RESAMPLE_X of genes by |predicted effect| stay on the parametric
                    # law -- those are the calls we have signal on (top-decile sign accuracy
                    # 0.74 vs 0.51 in the tail) -- and the uninformative tail is resampled,
                    # which removes the fixed-library lattice artefact and raises `C_hat`.
                    eff_abs = np.abs(np.asarray(context_effect[target_index], dtype=np.float64))
                    k = int(round(float(RESAMPLE_X) * eff_abs.size))
                    pmask = np.zeros(eff_abs.size, dtype=bool)
                    if k > 0:
                        pmask[np.argsort(-eff_abs, kind="stable")[:k]] = True
                    if RESAMPLE_INVERT:
                        # top-x RESAMPLED, remainder parametric
                        pmask = ~pmask
                    emit_mode = "resample_mix"
                    emit_kwargs = {
                        "control_counts": np.asarray(
                            resample_ctrl[picks].todense(), dtype=np.float64
                        ),
                        "control_profile": resample_profile,
                        "param_mask": pmask,
                    }
                spread_vector = emission_spread_vector(
                    context_effect[target_index],
                    b_ctrl,
                    min_cpm=min_cpm,
                    tail_spread=tail_spread,
                    keep_top=keep_top,
                    keep_top_frac=keep_top_frac,
                    keep_spread=keep_spread,
                    rank_score=(
                        None if rank_scores is None else rank_scores[target_index]
                    ),
                    spread_map=spread_map,
                    spread_lo=spread_lo,
                    spread_hi=spread_hi,
                )
                if emitter in RC_EMITTERS:
                    # 185 section 7. Pass 1 emits the incumbent's own block under the law
                    # above and reads its per-gene Wilcoxon z against THIS context's real
                    # control cells; pass 2 re-assigns that |z| multiset by descending |H|
                    # -- inside the incumbent's own BH call set and inside its complement
                    # (`rc_hybrid`), or across the whole gene axis (`rc_hist_a`) -- and
                    # realises it exactly. The delivered per-gene panel totals are pass 1's,
                    # count for count, and every emitted cell carries the same library
                    # total -- so `pds`, `mse`'s profile half and `nmae` are pinned under
                    # both, and only the p-value ORDER (and, for `rc_hist_a`, which genes
                    # end up called at the incumbent's own call COUNT) moves.
                    gate = emission_rank_axis(
                        context_effect[target_index], b_ctrl, min_cpm=min_cpm
                    )[0]
                    key = np.where(
                        gate,
                        np.asarray(context_effect[target_index], dtype=np.float64),
                        0.0,
                    )
                    block, rc_diag = rc_emit.rc_hybrid_block(
                        profile,
                        libraries[picks],
                        DISPERSION,
                        emission_rng,
                        spread=spread_vector,
                        overdispersion=overdispersion,
                        control=rc_control,
                        key=key,
                        jitter_seed=RC_JITTER_SEED_BASE + row // spec.cells_per_pert,
                        call_set_first=RC_CALL_SET_FIRST[emitter],
                        call_set=(
                            rc_emit.CallSet(RC_CALLSET_HEAD, RC_CALLSET_CPM,
                                            RC_CALLSET_Z_CAP)
                            if emitter == "rc_callset" else None
                        ),
                        emit_counts=emit_counts,
                    )
                    values, columns, pointers = rc_emit.to_csr(block)
                    rc_audit.append(rc_diag)
                    del block
                else:
                    values, columns, pointers = emit_counts(
                        profile,
                        libraries[picks],
                        DISPERSION,
                        emission_rng,
                        mode=emit_mode,
                        **emit_kwargs,
                        overdispersion=overdispersion,
                        spread=spread_vector,
                    )
                row_totals = np.add.reduceat(
                    values.astype(np.float64, copy=False),
                    pointers[:-1],
                )
                empty_rows = pointers[1:] == pointers[:-1]
                row_totals[empty_rows] = 0.0
                if np.any(row_totals <= 0) or np.any(row_totals > 1_000_000):
                    raise ValueError(f"{context}/{target}: invalid emitted row total")
                data_dataset.resize((nnz + values.size,))
                index_dataset.resize((nnz + columns.size,))
                data_dataset[nnz:] = values
                index_dataset[nnz:] = columns
                indptr[row + 1 : row + 1 + spec.cells_per_pert] = nnz + pointers[1:]
                nnz += values.size
                row += spec.cells_per_pert
                obs_pert.extend([target] * spec.cells_per_pert)
                obs_context.extend([str(context)] * spec.cells_per_pert)
                if nnz > MAX_STORED_ENTRIES:
                    raise RuntimeError(f"density cap exceeded at row {row}: {nnz:,}")

            del rc_control
            rc_control = None
            gc.collect()
            context_audit[str(context)] = {
                "control_rows": control_rows,
                "median_control_library": median_library,
                "emitted_rows": len(spec.perts) * spec.cells_per_pert,
                "stored_entries": nnz - context_start_nnz,
            }
            memory_checkpoint(f"context {context} complete", max_rss_gib)

        if row != n_rows:
            raise AssertionError(f"emitted {row} rows, expected {n_rows}")
        if int(indptr[-1]) != nnz:
            raise AssertionError(f"CSR indptr ends at {indptr[-1]}, expected {nnz}")
        matrix_group.create_dataset("indptr", data=indptr, dtype=np.int64)

        obs = handle.create_group("obs")
        obs.attrs["encoding-type"] = "dataframe"
        obs.attrs["encoding-version"] = "0.2.0"
        obs.attrs["_index"] = "_index"
        obs.attrs.create(
            "column-order",
            [PERT_COL, CONTEXT_COL],
            dtype=_str_dtype(),
        )
        _write_str(
            obs,
            "_index",
            np.asarray([f"cell_{index}" for index in range(n_rows)], dtype=object),
        )
        _write_categorical(obs, PERT_COL, np.asarray(obs_pert, dtype=object))
        _write_categorical(
            obs,
            CONTEXT_COL,
            np.asarray(obs_context, dtype=object),
        )

        var = handle.create_group("var")
        var.attrs["encoding-type"] = "dataframe"
        var.attrs["encoding-version"] = "0.2.0"
        var.attrs["_index"] = "_index"
        var.attrs.create(
            "column-order",
            np.empty(0, dtype=object),
            dtype=_str_dtype(),
        )
        _write_str(var, "_index", spec.genes.astype(object))

        for name in ("layers", "obsm", "varm", "obsp", "varp", "uns"):
            group = handle.create_group(name)
            group.attrs["encoding-type"] = "dict"
            group.attrs["encoding-version"] = "0.1.0"
        handle.attrs["encoding-type"] = "anndata"
        handle.attrs["encoding-version"] = "0.1.0"

    return {
        "path": str(output),
        "bytes": output.stat().st_size,
        "sha256": file_sha256(output),
        "n_rows": n_rows,
        "nnz": int(nnz),
        "nnz_per_cell": float(nnz / n_rows),
        "pct_of_density_cap": float(100 * nnz / MAX_STORED_ENTRIES),
        "contexts": context_audit,
        **({} if not rc_audit else {emitter: rc_diagnostics_summary(rc_audit)}),
    }


def runtime_metadata() -> dict[str, object]:
    packages = {}
    for name in ("anndata", "h5py", "numpy", "scipy"):
        packages[name] = importlib.metadata.version(name)
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "packages": packages,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    implementation_hashes = {
        **verify_frozen_implementation_hashes(),
        "scripts/make_submission_pm_cd4.py": file_sha256(Path(__file__).resolve()),
    }
    assert_frozen_contract()
    # v10 shipped 5.0; experiment 088 measured 10.0 as a sharp optimum. The
    # frozen contract above still pins the PerturbMatch helper's own MIN_CPM
    # against the module constant, so this override cannot silently desync it.
    min_cpm = float(args.min_cpm)
    effect_scale = float(args.effect_scale)
    # Per-gene amplitude shrinkage keyed on cross-donor agreement (experiment 096,
    # intervention 5). Off by default, so an unflagged build stays bit-identical to v11.
    shrink_key = str(args.shrink_key)
    shrink_pert_key = str(args.shrink_pert_key)
    shrink_gamma = None if args.shrink_gamma is None else float(args.shrink_gamma)
    shrink_toward = None if args.shrink_toward is None else str(args.shrink_toward)
    shrink_pert_gamma = (
        None if args.shrink_pert_gamma is None else float(args.shrink_pert_gamma)
    )
    shrink_pert_floor = (
        None if args.shrink_pert_floor is None else float(args.shrink_pert_floor)
    )
    shrink_pert_missing = (
        None if args.shrink_pert_missing is None else str(args.shrink_pert_missing)
    )
    # Each modifier is rejected against ITS OWN key, and --shrink-toward against either, so a
    # flag that cannot act is an abort rather than a silent no-op. Same contract as
    # --keep-spread without --keep-top.
    if shrink_key == "none" and shrink_gamma is not None:
        raise ValueError("--shrink-gamma does nothing without --shrink-key agreement")
    if shrink_pert_key == "none":
        for name, value in (("--shrink-pert-gamma", shrink_pert_gamma),
                            ("--shrink-pert-floor", shrink_pert_floor),
                            ("--shrink-pert-missing", shrink_pert_missing)):
            if value is not None:
                raise ValueError(
                    f"{name} does nothing without --shrink-pert-key donor-cos"
                )
    if shrink_key == "none" and shrink_pert_key == "none":
        if shrink_toward is not None:
            raise ValueError(
                "--shrink-toward does nothing without --shrink-key agreement or "
                "--shrink-pert-key donor-cos"
            )
    else:
        if shrink_toward is None:
            shrink_toward = AGREEMENT_SHRINK_TOWARD
        if shrink_toward not in ("mean", "zero"):
            raise ValueError(f"--shrink-toward must be mean or zero, got {shrink_toward!r}")
    if shrink_key != "none" and shrink_gamma is None:
        shrink_gamma = AGREEMENT_SHRINK_GAMMA
    if shrink_pert_key != "none":
        if shrink_pert_gamma is None:
            shrink_pert_gamma = PERT_SHRINK_GAMMA
        if shrink_pert_floor is None:
            shrink_pert_floor = PERT_SHRINK_FLOOR
        if shrink_pert_missing is None:
            shrink_pert_missing = PERT_SHRINK_MISSING
    # Emission-spread policy (experiment 093). All three default to the shipped scalar
    # SPREAD, so an unflagged build stays bit-identical to v11.
    tail_spread = None if args.tail_spread is None else float(args.tail_spread)
    keep_top = int(args.keep_top)
    keep_top_frac = None if args.keep_top_frac is None else float(args.keep_top_frac)
    keep_spread = None if args.keep_spread is None else float(args.keep_spread)
    spread_map = str(args.spread_map)
    if args.keep_key == "agreement" and tail_spread is None:
        raise ValueError("--keep-key agreement does nothing without --tail-spread")
    if tail_spread is None and (keep_top or keep_top_frac is not None):
        raise ValueError("--keep-top/--keep-top-frac do nothing without --tail-spread")
    if keep_spread is not None and tail_spread is None:
        raise ValueError("--keep-spread does nothing without --tail-spread")
    # Emission LAW (experiment 177). `quasi` is the shipped law and the default, so an
    # unflagged build writes no new meta key and is bit-identical to every prior build.
    emitter = str(args.emitter)
    resample_keep_top = (
        None if args.resample_keep_top is None else int(args.resample_keep_top)
    )
    if emitter != "resample_mix" and resample_keep_top is not None:
        raise ValueError("--resample-keep-top requires --emitter resample_mix")
    if emitter == "resample_mix":
        if RESAMPLE_X is not None:
            raise ValueError(
                "VCC_RESAMPLE_X and --emitter resample_mix are two spellings of the same "
                "channel; set exactly one"
            )
        if resample_keep_top is None:
            resample_keep_top = keep_top
        if resample_keep_top < 0:
            raise ValueError(f"--resample-keep-top must be >= 0, got {resample_keep_top!r}")
    bundle = args.bundle.expanduser().resolve()
    donor_paths = {
        K562: args.k562.expanduser().resolve(),
        "orion_hct116": args.hct116.expanduser().resolve(),
        "orion_hek293t": args.hek293t.expanduser().resolve(),
    }
    perturbmatch = args.perturbmatch.expanduser().resolve()
    cd4_root = args.cd4_root.expanduser().resolve()
    operator_path = args.whitener.expanduser().resolve()
    operator_metadata_path = args.whitener_metadata.expanduser().resolve()
    if args.output is None:
        if not args.repin_panel_support:
            raise ValueError("-o/--output is required for a build")
        output = metadata_path = None
    else:
        output = args.output.expanduser().resolve()
        metadata_path = Path(str(output) + ".meta.json")
        if output.exists() or metadata_path.exists():
            raise FileExistsError(f"refusing to overwrite {output} or {metadata_path}")

    # --- context adapter -------------------------------------------------------------
    # The bundle declares its own contexts and the frozen table pins them. Three
    # statements of the same fact have to agree before a byte is hashed: the manifest's,
    # the frozen table's, and -- when given -- the operator's `--expect-contexts`. The
    # default for the last one is the manifest's list, so an unflagged build on a correctly
    # re-pinned bundle just works, while a bundle/table mismatch still aborts loudly.
    frozen_path = args.frozen_inputs.expanduser().resolve()
    frozen_bundle = load_frozen_bundle(frozen_path)
    bundle_contexts = read_bundle_contexts(bundle)
    frozen_contexts = tuple(map(str, frozen_bundle["contexts"]))
    expect_contexts = (
        bundle_contexts
        if args.expect_contexts is None
        else tuple(
            part.strip() for part in str(args.expect_contexts).split(",") if part.strip()
        )
    )
    if bundle_contexts != frozen_contexts:
        raise ValueError(
            f"bundle manifest declares contexts {list(bundle_contexts)} but "
            f"{frozen_path} pins {list(frozen_contexts)}; re-pin the frozen inputs "
            f"with --repin-frozen-inputs before building on a new panel"
        )
    if bundle_contexts != expect_contexts:
        raise ValueError(
            f"--expect-contexts {list(expect_contexts)} != the bundle's "
            f"{list(bundle_contexts)}"
        )
    print(f"[contexts] {list(bundle_contexts)} (manifest, frozen table and "
          f"--expect-contexts all agree)", flush=True)

    input_hashes = collect_input_hashes(
        bundle,
        bundle_contexts,
        donor_paths,
        perturbmatch,
        cd4_root,
        operator_path,
        operator_metadata_path,
        whitener_id=args.whitener_id,
        frozen_bundle=frozen_bundle,
    )
    memory_checkpoint("input hashing", args.max_rss_gib)

    spec = ChallengeSpec.load(bundle)
    if tuple(map(str, spec.contexts)) != bundle_contexts:
        raise ValueError(
            f"ChallengeSpec read contexts {list(spec.contexts)} from a manifest that "
            f"declares {list(bundle_contexts)}"
        )
    expected_perts = int(frozen_bundle["n_constructs"])
    if len(spec.perts) != expected_perts:
        raise ValueError(
            f"expected {expected_perts} perturbations (frozen table), got {len(spec.perts)}"
        )
    expected_cells = int(frozen_bundle["cells_per_pert"])
    if int(spec.cells_per_pert) != expected_cells:
        raise ValueError(
            f"expected {expected_cells} cells per perturbation (frozen table), got "
            f"{spec.cells_per_pert}"
        )

    donors = {name: Donor.load(path) for name, path in donor_paths.items()}
    model = MultiDonorModel(
        donors,
        spec.genes,
        gene_norm=False,
        eb=False,
    )
    weights = donor_weights(model)
    memory_checkpoint("three-donor model", args.max_rss_gib)

    pm = EXP061.prepare_perturbmatch(
        perturbmatch,
        tuple(map(str, spec.genes)),
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    operator = load_production_operator(
        model,
        tuple(map(str, spec.perts)),
        operator_path,
        operator_metadata_path,
        expected_sha256=(
            args.expected_whitener_sha256
            or WHITENER_REGISTRY[args.whitener_id][0]
        ),
        whitener_id=args.whitener_id,
    )

    panel_positions = np.asarray(
        [
            {str(gene): index for index, gene in enumerate(spec.genes)}[str(target)]
            for target in spec.perts
        ],
        dtype=np.int64,
    )
    cd4_banks, cd4 = prepare_frozen_cd4(
        cd4_root,
        tuple(map(str, spec.perts)),
        tuple(map(str, spec.genes)),
        tuple(map(str, spec.perts)),
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if H1_UNION_G is None:
        h1_union = None
        h1_union_audit = None
    else:
        h1_union, h1_union_details = prepare_h1_union(
            tuple(map(str, spec.perts)),
            tuple(map(str, spec.genes)),
            tuple(map(str, spec.perts)),
            max_rss_gib=args.max_rss_gib,
        )
        h1_union_audit = h1_union_support_audit(h1_union, h1_union_details)
        print(
            f"[h1-union] g={H1_UNION_G:g}  "
            f"{h1_union_audit['carrying_targets']}/{len(spec.perts)} targets carried "
            f"({h1_union_audit['targets_carried_by_both_sources']} by both sources), "
            f"{h1_union_audit['supported_entries']} supported entries",
            flush=True,
        )
        memory_checkpoint("h1 union residual", args.max_rss_gib)
    effects, effect_audit = build_production_effects(
        model,
        tuple(map(str, spec.perts)),
        operator.whitener,
        operator.columns,
        panel_positions,
        pm=pm,
        cd4=cd4,
        h1_union=h1_union,
        h1_union_audit=h1_union_audit,
        weights=weights,
    )
    rank_scores = None
    shrink_agreement = None
    shrink_scales = None
    shrink_center = None
    agreement = None
    blank = None
    if args.keep_key == "agreement" or shrink_key != "none":
        # (n_targets, n_genes) on the target axis. 300 x 18,533 float64 is 44 MB, so this
        # is cheap next to the 15 GiB the effect matrix already occupies. Built ONCE and
        # shared by --keep-key and --shrink-key so the two can never key on differently
        # computed agreement.
        # --keep-key always uses the frozen, board-validated unanimity method. --shrink-key
        # may ask for the pairwise one; when both are on and they differ, build both.
        agreement = np.vstack(
            [model.donor_agreement(str(t)) for t in spec.perts]
        ).astype(np.float64)
        blank = ~agreement.any(axis=1)
        shrink_agreement = agreement
        if shrink_key == "agreement-pairwise":
            shrink_agreement = np.vstack(
                [donor_agreement_pairwise(model, str(t)) for t in spec.perts]
            ).astype(np.float64)
        graded = (shrink_agreement > 0).sum(axis=1)
        print(f"[shrink-key] key={shrink_key}: graded (nonzero) gene fraction per target, "
              f"median {float(np.median(graded)) / len(spec.genes):.3f}, "
              f"p5 {float(np.quantile(graded, 0.05)) / len(spec.genes):.3f}, "
              f"p95 {float(np.quantile(graded, 0.95)) / len(spec.genes):.3f}", flush=True)
    if args.keep_key == "agreement":
        rank_scores = np.array(agreement, dtype=np.float64, copy=True)
        # A target carried by only ONE donor has nothing to agree about, so its whole row is
        # zero -- and an all-tied key makes the graded map rank by gene index, which is
        # arbitrary. Experiment 095 measured a random ordering at -0.0048 overall against
        # production, so leaving those rows tied actively costs. Fall back to |effect| for
        # them, which is the key the block map has always used.
        if blank.any():
            rank_scores[blank] = np.abs(
                np.nan_to_num(np.asarray(effects, dtype=np.float64)[blank])
            )
        print(f"[keep-key] cross-donor agreement for {int((~blank).sum())}/"
              f"{len(spec.perts)} targets; {int(blank.sum())} fell back to |effect| "
              f"(carried by a single donor)", flush=True)
    if shrink_key != "none":
        # The shrinkage takes the RAW agreement -- deliberately NOT the |effect| fallback
        # above. A target no two donors carry has no cross-donor evidence, so
        # `agreement_shrink_vector` returns all ones for it and it is emitted exactly as an
        # unflagged build emits it. Substituting a self-referential key there would shrink
        # on a signal experiment 096 measured as useless-to-harmful (rho with the true
        # per-perturbation cosine: -0.001/+0.253 for our own magnitude, -0.337/-0.102 for
        # prediction uniqueness, across the two harnesses).
        # `panel_positions[i]` is the gene column of target i, forced to alpha = 1: the
        # metric never sees that coordinate, so protecting it is free, and shipping a
        # silently-shrunk knockdown is a defect no number we check would catch.
        shrink_scales = np.vstack([
            agreement_shrink_vector(row, shrink_gamma, protect=int(position))
            for row, position in zip(shrink_agreement, panel_positions)
        ])
        degenerate = int(sum(1 for row in shrink_agreement if np.unique(row).size < 2))
    if shrink_pert_key == "donor-cos":
        # A whole-vector direction agreement for the perturbation as a unit, orthogonal to the
        # per-gene magnitude key above; experiment 097 measured them stacking, realised
        # (1.246 / 1.176 alone, 1.473 / 1.580 on top of the per-gene key). Applied as `z[:, None] *
        # alpha` into the same matrix, so a per-row scalar cannot reorder genes within a row --
        # which is also why it leaves the DE membership of a row untouched.
        cosine_score, cosine_carried = donor_cosine_key(
            model, [str(t) for t in spec.perts]
        )
        pert_weight = pert_shrink_vector(
            cosine_score, cosine_carried, shrink_pert_gamma, shrink_pert_floor,
            missing=shrink_pert_missing,
        )
        if shrink_scales is None:
            shrink_scales = np.ones(
                (len(spec.perts), len(spec.genes)), dtype=np.float64
            )
        shrink_scales = shrink_scales * pert_weight[:, None]
        zq = np.quantile(pert_weight, [0.0, 0.25, 0.5, 0.75, 1.0])
        sq = np.quantile(cosine_score[cosine_carried], [0.0, 0.5, 1.0]) if (
            cosine_carried.any()) else np.zeros(3)
        print(f"[shrink-pert] donor-cos gamma={shrink_pert_gamma} floor={shrink_pert_floor} "
              f"missing={shrink_pert_missing} on {int(cosine_carried.sum())}/"
              f"{len(spec.perts)} targets; {int((~cosine_carried).sum())} have fewer than two "
              f"donors and got z={np.round(pert_weight[~cosine_carried], 4).tolist()[:4]}. "
              f"key min/med/max = {sq[0]:+.4f}/{sq[1]:+.4f}/{sq[2]:+.4f}; "
              f"z min/q25/med/q75/max = {zq[0]:.4f}/{zq[1]:.4f}/{zq[2]:.4f}/{zq[3]:.4f}/"
              f"{zq[4]:.4f}; mean z^2 = {float((pert_weight ** 2).mean()):.4f}", flush=True)
        thin = int((pert_weight < 0.05).sum())
        if thin:
            print(f"[shrink-pert] WARNING {thin} targets carry z < 0.05 and will emit close "
                  f"to the bare control profile. That is an n_pred loss the board reads, and "
                  f"holding the cliff raises the global amplitude by ~1/sqrt(mean z^2) = "
                  f"{1.0 / max(float(np.sqrt((pert_weight ** 2).mean())), 1e-12):.2f}x. "
                  f"experiments/097-direction recommends --shrink-pert-floor 0.2.", flush=True)

    if shrink_scales is not None:
        if shrink_toward == "mean":
            # The James-Stein centre: the grand mean predicted log effect per gene, taken
            # over all 300 perturbations of our OWN effect matrix. Nothing from the
            # reference enters it, so there is no oracle leak by rule 4.
            #
            # The case for it was that this is the metric's own 0-anchor -- the context's
            # mean perturbation response assigned identically to every perturbation, which
            # the board prices at raw ~0.9879 against 1.0424 for the control. MEASURED, it
            # is not: the anchor is the REFERENCE's mean effect and ours does not estimate
            # it (cos +0.016 / -0.066 across the two harnesses), while shrinking toward it
            # inflates the (A - B) term the master formula drops. `mech_shrink_target.py`
            # has `zero` beating `mean` 12/12. See --shrink-toward's help.
            #
            # Each gene's own panel row also contributes 1/300 of its own on-target
            # knockdown to `mbar` at that gene; the metric drops each row's own target
            # coordinate and this does not. That is a 0.3% contamination of one coordinate
            # per gene, left in for simplicity.
            shrink_center = np.asarray(effects, dtype=np.float64).mean(axis=0)
        # Rule 3: print the DISTRIBUTION of the intervening quantity, not a summary, so a
        # null on this arm can be told apart from the shrinkage never having happened.
        shrunk = shrink_scales
        if shrunk.size:
            q = np.quantile(shrunk, [0.0, 0.25, 0.5, 0.75, 1.0])
            half = float((shrunk > 0.5).sum(axis=1).mean())
        else:
            q, half = np.ones(5), 0.0
        print(f"[shrink-key] agreement shrinkage gamma={shrink_gamma} toward "
              f"{shrink_toward} on {len(spec.perts) - degenerate}/{len(spec.perts)} targets; "
              f"{degenerate} left at alpha=1 (key has fewer than two distinct values, of "
              f"which {int(blank.sum())} are all-zero). alpha over "
              f"shrunk targets, min/q25/med/q75/max = {q[0]:.3e}/{q[1]:.3e}/{q[2]:.3e}/"
              f"{q[3]:.3e}/{q[4]:.3e}; {half:.0f} genes per target keep alpha > 0.5",
              flush=True)
        if shrink_center is not None:
            print("[shrink-key] WARNING --shrink-toward mean is the NON-default option and "
                  "was measured WORSE than 'zero' on both harnesses at every gamma "
                  "(experiments/096-mse-floor/mech_shrink_target.py and mech_shrink_exact.py, "
                  "12/12), and worse than not shrinking at all on h1-test. Drop the flag "
                  "unless this arm IS the comparison.", flush=True)
            cq = np.quantile(shrink_center, [0.0, 0.25, 0.5, 0.75, 1.0])
            print(f"[shrink-key] centre mbar over {len(spec.genes)} genes, "
                  f"min/q25/med/q75/max = {cq[0]:+.5f}/{cq[1]:+.5f}/{cq[2]:+.5f}/"
                  f"{cq[3]:+.5f}/{cq[4]:+.5f}; ||mbar||={float(np.linalg.norm(shrink_center)):.5f}, "
                  f"mean row ||effect||="
                  f"{float(np.linalg.norm(effects, axis=1).mean()):.5f}", flush=True)

    # PerturbMatch / CD4 coverage is a property of (panel x donor tables), so a new
    # 300-target panel moves it and the A/B/C values would abort every D/E/F build. The
    # expected pair therefore lives in the frozen bundle table next to the file hashes.
    # `null` means "this panel has never been pinned": the observed values are printed
    # loudly and recorded in the meta, and `--repin-panel-support` freezes them.
    pm_modified = int(effect_audit["perturbmatch"]["both_supported_modified_count"])
    cd4_complete = int(effect_audit["cd4"]["fully_supported_targets"])
    observed_support = {
        "perturbmatch_replacement_targets": pm_modified,
        "cd4_complete_targets": cd4_complete,
    }
    pinned_support = frozen_bundle.get("panel_support")
    if pinned_support is None:
        print("[panel-support] UNPINNED for this panel. Observed "
              f"perturbmatch_replacement_targets={pm_modified}, "
              f"cd4_complete_targets={cd4_complete}. Re-run with --repin-panel-support to "
              "freeze them before the round is built for real.", flush=True)
    if args.repin_panel_support:
        table = json.loads(frozen_path.read_text())
        table["panel_support"] = dict(sorted(observed_support.items()))
        atomic_write_text(frozen_path, json.dumps(table, indent=2) + "\n")
        load_frozen_bundle(frozen_path)
        print(f"[repin] panel_support pinned to {observed_support} in "
              f"{frozen_path}; no submission was built", flush=True)
        return {"status": "panel_support_repinned", "panel_support": observed_support,
                "submission": False}
    if pinned_support is not None:
        if pm_modified != int(pinned_support["perturbmatch_replacement_targets"]):
            raise ValueError(
                f"PerturbMatch replacement support changed: {pm_modified} != "
                f"{pinned_support['perturbmatch_replacement_targets']}"
            )
        if cd4_complete != int(pinned_support["cd4_complete_targets"]):
            raise ValueError(
                f"complete CD4 mean3 support changed: {cd4_complete} != "
                f"{pinned_support['cd4_complete_targets']}"
            )
    memory_checkpoint("final production effects", args.max_rss_gib)

    coverage = model.coverage(spec.perts)
    source_support = main_source_support_audit(model, spec.perts)
    operator_audit = {
        "path": str(operator_path),
        "metadata_path": str(operator_metadata_path),
        "sha256": operator.file_sha256,
        "metadata_sha256": operator.metadata_sha256,
        "metadata": operator.metadata,
        "columns": len(operator.columns),
        "columns_sha256": hashlib.sha256(
            np.ascontiguousarray(operator.columns, dtype="<i8").tobytes()
        ).hexdigest(),
        # The panel guard's own record (experiment 235). `columns_sha256` above hashes the
        # integer indices and cannot see a re-ordered gene axis; this hashes the ordered
        # SYMBOLS, and unlike `columns_sha256` it was ASSERTED against the operator's own
        # recorded axis before the build ran. Deliberately here and not in `hashes.inputs`:
        # every key there is resolved back to a file on disk by
        # `experiments/199-def-runbook/rule10_199.py`, a DEF-RUNBOOK section 4 known-answer
        # row, and this digest is not a file.
        "axis_sha256": operator.axis_sha256,
        "axis_sha256_checked_against": operator.axis_sha256_source,
        "panel_sha256": operator_axis_sha256(map(str, spec.perts)),
        "all_panel_coordinates_excluded": not bool(
            np.intersect1d(operator.columns, panel_positions).size
        ),
    }

    # --- VCC_ALPHA=knee -------------------------------------------------------------
    # Solve one alpha per context from the effect matrix just built, then REBUILD the
    # effect matrix at that alpha through the same production path rather than scaling the
    # finished one. The pipeline happens to be homogeneous of degree 1 in alpha under
    # a15/a23's flags (`model.effect` is linear in it, whitening renormalises by a norm
    # RATIO, and `VCC_CD4_TARGET_F` solves a fraction that is itself scale-relative), so
    # scaling would give the same answer today -- but that is a property of the current
    # configuration, not of the knob, and `VCC_CD4_FRACTION` alone already breaks it.
    # `effects` and its audit are left exactly as built at `alpha_ref`, so the rule-10
    # effects-audit diff against the incumbent stays clean.
    context_effects = None
    knee_audit: dict[str, object] | None = None
    if ALPHA_MODE == "knee":
        if emitter != "quasi" or RESAMPLE_X is not None:
            raise ValueError(
                "VCC_ALPHA=knee currently measures `C_hat` through the quasi emitter only; "
                "a resampling emission law needs its own `C_hat` probe before the knee can "
                "be solved for it"
            )
        spread_kwargs = {
            "tail_spread": tail_spread,
            "keep_top": keep_top,
            "keep_top_frac": keep_top_frac,
            "keep_spread": keep_spread,
            "spread_map": spread_map,
            "spread_lo": float(args.spread_lo),
            "spread_hi": float(args.spread_hi),
        }
        terms_kwargs = dict(
            min_cpm=min_cpm,
            effect_scale=effect_scale,
            shrink_scales=shrink_scales,
            shrink_center=shrink_center,
            rank_scores=rank_scores,
            spread_kwargs=spread_kwargs,
            overdispersion=float(args.overdispersion),
            cells_per_pert=int(spec.cells_per_pert),
        )
        frames = {
            str(context): load_control_frame(
                bundle, str(context), spec.genes, cells_per_pert=spec.cells_per_pert
            )
            for context in spec.contexts
        }
        memory_checkpoint("knee control frames", args.max_rss_gib)
        solved = solve_knee_alphas(
            spec, effects, frames, alpha_ref=ALPHA, **terms_kwargs
        )
        context_effects = {}
        for context in map(str, spec.contexts):
            per_context, per_audit = build_production_effects(
                model,
                tuple(map(str, spec.perts)),
                operator.whitener,
                operator.columns,
                panel_positions,
                pm=pm,
                cd4=cd4,
                h1_union=h1_union,
                h1_union_audit=h1_union_audit,
                alpha=float(solved[context]["alpha"]),
                weights=weights,
            )
            context_effects[context] = per_context
            solved[context]["effect_sha256"] = per_audit["final_sha256"]
            solved[context]["cd4_fraction"] = per_audit["cd4"]["fraction"]
            # Non-vacuity for the whole feature: if `alpha` never reached the effect
            # matrix, this ratio would read 1.0 for every context and the verification
            # below would read the ORIGINAL rho back.
            ratio = float(
                np.linalg.norm(per_context) / max(float(np.linalg.norm(effects)), 1e-300)
            )
            solved[context]["frobenius_ratio_to_reference"] = ratio
            solved[context]["homogeneity_residual"] = float(
                np.max(np.abs(per_context - solved[context]["knee_scale"] * effects))
                / max(float(np.max(np.abs(effects))), 1e-300)
            )
            memory_checkpoint(f"knee effects {context}", args.max_rss_gib)
        # Verification pass: re-measure `rho` on the REBUILT matrices. On the knee it must
        # read 1.000. This is not an identity -- `profile_from_effect` is nonlinear in the
        # effect (log1p/expm1, the CPM gate, and a23's minimum-norm mass fix), so `B` is
        # only approximately quadratic in alpha, and a mis-wired per-context alpha would
        # leave `rho` at its pre-correction value.
        panel_cols = panel_gene_columns(spec)
        for context in map(str, spec.contexts):
            after = knee_context_terms(
                context_effects[context], frames[context], panel_cols, **terms_kwargs
            )
            solved[context]["after"] = {
                key: after[key]
                for key in ("A", "B", "claim", "rho", "claim_fraction_of_cap", "C_hat_mean")
            }
            print(
                f"[knee] {context}: after rebuild at alpha "
                f"{solved[context]['alpha']:.6f} -> rho {after['rho']:.4f} "
                f"(was {solved[context]['rho']:.4f}), B {after['B']:.6f} vs claim "
                f"{after['claim']:.6f}",
                flush=True,
            )
        knee_audit = {
            "mode": "knee",
            "alpha_ref": float(ALPHA),
            "subsample_targets": KNEE_SUBSAMPLE,
            "jackknife_chunk": KNEE_JACKKNIFE_CHUNK,
            "chat_seed_base": KNEE_CHAT_SEED,
            "per_context": solved,
        }
        del frames
        gc.collect()

    del donors, pm, cd4_banks, cd4, h1_union, operator, model
    gc.collect()
    memory_checkpoint("before challenge emission", args.max_rss_gib)

    output_audit = write_submission_from_effects(
        spec,
        effects,
        bundle,
        output,
        context_effects=context_effects,
        max_rss_gib=args.max_rss_gib,
        min_cpm=min_cpm,
        effect_scale=effect_scale,
        shrink_scales=shrink_scales,
        shrink_center=shrink_center,
        tail_spread=tail_spread,
        keep_top=keep_top,
        keep_top_frac=keep_top_frac,
        keep_spread=keep_spread,
        rank_scores=rank_scores,
        spread_map=spread_map,
        spread_lo=float(args.spread_lo),
        spread_hi=float(args.spread_hi),
        overdispersion=float(args.overdispersion),
        emitter=emitter,
        resample_keep_top=resample_keep_top,
    )
    report = {
        "model": "v9 + frozen PerturbMatch K562-slot replacement + post-whiten CD4",
        "status": "challenge_h5ad_built",
        "submission": False,
        "config": {
            "contexts": list(map(str, spec.contexts)),
            "targets": len(spec.perts),
            "genes": len(spec.genes),
            "cells_per_target_context": spec.cells_per_pert,
            "donors": list(PRODUCTION_DONORS),
            "donor_weights": weights,
            "alpha": ALPHA,
            # Written only when the amplitude is NOT a single global constant, so a default
            # build's config block stays byte-identical to every prior build's -- the same
            # convention `mass_fix`, `sign_repair` and `emitter` follow.
            **(
                {}
                if ALPHA_MODE == "fixed"
                else {
                    "alpha_mode": ALPHA_MODE,
                    "alpha_ref": ALPHA,
                    "alpha_per_context": {
                        context: knee_audit["per_context"][context]["alpha"]
                        for context in map(str, spec.contexts)
                    },
                }
            ),
            **({} if args.expect_contexts is None else
               {"expect_contexts": str(args.expect_contexts)}),
            "beta": BETA,
            "whiten_gamma": WHITEN_GAMMA,
            "whiten_lambda": WHITEN_LAM,
            "resample_x": RESAMPLE_X,
            "resample_invert": RESAMPLE_INVERT if RESAMPLE_X is not None else None,
            "whiten_axis": "union",
            "whiten_renorm": "candidate-specific on-axis preserve, as v9",
            "perturbmatch_blend": PM_BLEND,
            "perturbmatch_target_policy": (
                "replace K562 slot only when K562 and PerturbMatch both "
                "support the target"
            ),
            "cd4_conditions": list(CD4_NAMES),
            "cd4_mode": CD4_MODE,
            "cd4_fraction": CD4_FRACTION,
            # The fraction actually applied, the measured experiment-099 `f`, and the frame it
            # was measured in are under `effects.cd4` -- they differ from the constant above
            # whenever `VCC_CD4_TARGET_F` is set.
            "cd4_target_f": CD4_TARGET_F,
            "cd4_target_policy": "complete mean3 only",
            # Written only when the coverage-aware second residual is ON, so a default
            # build's config block stays byte-identical to every prior build's -- the same
            # convention `mass_fix`, `sign_repair` and `emitter` follow.
            **(
                {}
                if H1_UNION_G is None
                else {
                    "h1_union_g": H1_UNION_G,
                    "h1_union_sources": [name for name, _ in H1_UNION_SOURCES],
                    "h1_union_sizing_rule": "median over carrying rows",
                }
            ),
            "min_cpm": min_cpm,
            "gate_unexpressed": True,
            **({} if MASS_FIX is None else {"mass_fix": MASS_FIX}),
            **({} if SIGN_REPAIR is None else {"sign_repair": SIGN_REPAIR}),
            # Written only when the emission law is NOT the shipped one, so a default
            # build's config block stays byte-identical to every prior build's rather than
            # differing by two null keys (the convention `mass_fix` and `sign_repair` set).
            **(
                {}
                if emitter == "quasi"
                else {"emitter": emitter}
            ),
            **(
                {}
                if emitter != "resample_mix"
                else {"resample_keep_top": resample_keep_top}
            ),
            **(
                {}
                if emitter not in RC_EMITTERS
                else {
                    "rc_key": RC_KEY,
                    "rc_sign_policy": RC_SIGN_POLICY,
                    "rc_scorer_min_cpm": RC_SCORER_MIN_CPM,
                    "rc_jitter_seed_base": RC_JITTER_SEED_BASE,
                }
            ),
            # The prescribed call set's three numbers, written ONLY under `rc_callset`, so
            # `rc_hybrid`'s and `rc_hist_a`'s config blocks are untouched and the rule-10
            # diff between this arm and `rc_hybrid` is exactly `emitter` plus these three.
            **(
                {}
                if emitter != "rc_callset"
                else {
                    "rc_callset_head": RC_CALLSET_HEAD,
                    "rc_callset_cpm": RC_CALLSET_CPM,
                    "rc_callset_z_cap": RC_CALLSET_Z_CAP,
                }
            ),
            "emission": "quasi exact-total",
            "dispersion_argument": DISPERSION,
            "spread": SPREAD,
            "whitener_id": args.whitener_id,
            "tail_spread": tail_spread,
            "keep_top": keep_top,
            "keep_top_frac": keep_top_frac,
            "keep_spread": keep_spread,
            "keep_key": args.keep_key,
            "spread_map": spread_map,
            "effect_scale": float(args.effect_scale),
            "shrink_key": shrink_key,
            "shrink_gamma": shrink_gamma,
            "shrink_toward": shrink_toward,
            "shrink_rank_scale": AGREEMENT_SHRINK_RANK_SCALE,
            "shrink_pert_key": shrink_pert_key,
            "shrink_pert_gamma": shrink_pert_gamma,
            "shrink_pert_floor": shrink_pert_floor,
            "shrink_pert_missing": shrink_pert_missing,
            "spread_lo": float(args.spread_lo),
            "spread_hi": float(args.spread_hi),
            "overdispersion": float(args.overdispersion),
            "seed": SEED,
            "library_pick_seed": SEED,
            "emission_seed": SEED + EMISSION_SEED_OFFSET,
        },
        "donor_coverage": coverage,
        "main_source_support": source_support,
        "effects": effect_audit,
        # New top-level blocks, deliberately OUTSIDE `config` and `hashes` so that a default
        # A/B/C build's config block and input-hash table stay byte-identical to a23's and
        # the rule-10 diff stays readable.
        "frozen_bundle": {
            "path": str(frozen_path),
            "sha256": file_sha256(frozen_path),
            "panel_id": frozen_bundle.get("panel_id"),
            "contexts": list(frozen_contexts),
            "n_constructs": int(frozen_bundle["n_constructs"]),
            "cells_per_pert": int(frozen_bundle["cells_per_pert"]),
            "panel_support_pinned": pinned_support is not None,
            "panel_support": pinned_support,
            "panel_support_observed": observed_support,
        },
        **({} if knee_audit is None else {"knee": knee_audit}),
        "operator": operator_audit,
        "hashes": {
            "inputs": input_hashes,
            "implementation": implementation_hashes,
        },
        "output": output_audit,
        "runtime": runtime_metadata(),
        "peak_rss_gib": peak_rss_gib(),
    }
    atomic_write_text(
        metadata_path,
        json.dumps(report, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(report, indent=2, sort_keys=True), flush=True)
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, default=Path("data/challenge"))
    parser.add_argument("--k562", type=Path, default=Path("data/donors/k562_gwps"))
    parser.add_argument(
        "--hct116",
        type=Path,
        default=Path("data/donors/orion_hct116"),
    )
    parser.add_argument(
        "--hek293t",
        type=Path,
        default=Path("data/donors/orion_hek293t"),
    )
    parser.add_argument(
        "--perturbmatch",
        type=Path,
        default=Path("data/external/perturbmatch"),
    )
    parser.add_argument(
        "--cd4-root",
        type=Path,
        default=Path("data/external/cd4t/cloud-compact/cd4_compact"),
    )
    parser.add_argument(
        "--whitener",
        type=Path,
        default=Path("data/donors/whitener_union.npz"),
    )
    parser.add_argument(
        "--whitener-metadata",
        type=Path,
        default=Path("data/donors/whitener_union.json"),
    )
    parser.add_argument(
        "--expected-whitener-sha256",
        default=None,
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=None,
        help=(
            "Output .h5ad. Required for every build; omit it only for "
            "--repin-frozen-inputs and --repin-panel-support, which never write one."
        ),
    )
    parser.add_argument(
        "--expect-contexts",
        default=None,
        help=(
            "Comma-separated context ids the operator BELIEVES this bundle carries, e.g. "
            "'A,B,C' or 'D,E,F'. Defaults to whatever `manifest.json` declares, which is "
            "also cross-checked against scripts/frozen_inputs.json -- so an unflagged build "
            "on a correctly re-pinned bundle just works, while a bundle that is not the one "
            "the frozen table pins aborts before a byte is hashed. Pass it when you want the "
            "build to fail rather than proceed if the bundle on disk is not the round you "
            "think it is."
        ),
    )
    parser.add_argument(
        "--frozen-inputs",
        type=Path,
        default=FROZEN_BUNDLE_PATH,
        help=(
            "The frozen bundle table this build verifies against. Defaults to the "
            "production one. Point it at a throwaway copy to REHEARSE a D/E/F re-pin "
            "without moving the table a concurrent A/B/C build is using -- the code path "
            "is identical, only the file differs, and the path and its sha256 are recorded "
            "in the meta under `frozen_bundle`."
        ),
    )
    parser.add_argument(
        "--repin-frozen-inputs",
        action="store_true",
        help=(
            "Hash the bundle on disk and rewrite scripts/frozen_inputs.json, then EXIT "
            "without building. This is the documented way to move the frozen-input table "
            "onto a new panel (contexts D/E/F, a new 300-target `pert_counts.csv`); a build "
            "never re-pins itself. Prints the old -> new diff. `panel_support` is preserved "
            "only when pert_counts.csv and gene_names.csv are unchanged, and set to null "
            "otherwise -- pin it with --repin-panel-support once."
        ),
    )
    parser.add_argument(
        "--repin-panel-support",
        action="store_true",
        help=(
            "Run the effects build, write the observed PerturbMatch / CD4 target-coverage "
            "pair into scripts/frozen_inputs.json, and EXIT before emission. Needed once "
            "per new panel, because that pair is a property of (panel x donor tables) and "
            "the A/B/C values (248 / 276) would abort every D/E/F build."
        ),
    )
    parser.add_argument("--chunk-rows", type=int, default=128)
    parser.add_argument("--max-rss-gib", type=float, default=32.0)
    parser.add_argument(
        "--min-cpm",
        type=float,
        default=MIN_CPM,
        help="CPM gate threshold. v10 shipped 5.0; experiment 088 measured 10.0.",
    )
    parser.add_argument(
        "--tail-spread",
        type=float,
        default=None,
        help=(
            "Per-gene emission spread for genes OUTSIDE the kept top-K, in [0, 1]. "
            "Unset (the default) emits every gene at the shipped scalar SPREAD=0.25, "
            "so the build is bit-identical to v11. Raising it raises those genes' "
            "Wilcoxon p-values while preserving each gene's panel total exactly, which "
            "shrinks our significant call set without moving pds/mse/nmae "
            "(experiment 093)."
        ),
    )
    parser.add_argument(
        "--keep-top",
        type=int,
        default=0,
        help=(
            "Number of genes per target held at SPREAD=0.25 when --tail-spread is set, "
            "ranked by |predicted effect| over the CPM-gated axis. 0 (the default) "
            "applies --tail-spread to every gene. Experiment 093 recommends "
            "--tail-spread 1.0 --keep-top 850 (+0.0202 avg6 on the HCT116 harness)."
        ),
    )
    parser.add_argument(
        "--overdispersion",
        type=float,
        default=0.0,
        help=(
            "Gamma-Poisson dispersion phi for the emitted scatter, Var = mu + phi*mu^2. "
            "0.0 (the default) is the shipped Poisson scatter and is bit-identical to v11. "
            "Our emitter is under-dispersed ~2.3x against the real challenge controls "
            "(per-gene CV 0.38-0.40 against 0.84-0.92), and experiment 095 measured this "
            "as the ONLY emission lever that moves the DE ordering: realised on a repaired "
            "harness, phi=0.30 takes reach 0.0838 -> 0.1888 with jac up and nmae invariant, "
            "against a 5e-4 noise floor. It shrinks n_pred sharply, so probe it with "
            "experiments/094-reach-ordering/probe_de.py before spending a slot."
        ),
    )
    parser.add_argument(
        "--keep-key",
        default="effect",
        choices=["effect", "agreement"],
        help=(
            "What ranks genes for --keep-top / --spread-map graded. 'effect' is "
            "|predicted effect|. 'agreement' is MultiDonorModel.donor_agreement, the "
            "geometric magnitude gated on cross-donor sign unanimity, which experiment 095 "
            "measured at +0.00678 overall as a full ordering against |lfc|'s +0.00448 -- "
            "the best buildable ordering key found. It needs --spread-map graded; under "
            "the block map it LOSES to 'effect' (+0.00341 vs +0.00583)."
        ),
    )
    parser.add_argument(
        "--spread-map",
        default="block",
        choices=["block", "graded"],
        help=(
            "'block' is the two-level v14/v15 shape: --tail-spread everywhere, "
            "--keep-spread on the top-K. 'graded' assigns spread by rank across the whole "
            "gated axis, from --spread-lo at the most confident gene to --spread-hi at the "
            "least, so the realised p-value order can follow the key rather than merely "
            "putting one block in front. Its effect on n_pred is UNMEASURED -- probe it "
            "with experiments/094-reach-ordering/probe_de.py before spending a slot."
        ),
    )
    parser.add_argument(
        "--effect-scale", type=float, default=1.0,
        help=(
            "Multiply the gated log-space effect by this factor before it is exponentiated. "
            "1.0 (default) reproduces every prior build bit-for-bit. Values below 1 shrink "
            "the prediction toward the context's own control profile, which is where the "
            "expression-error member measures distance; 0.0 emits the control unchanged. "
            "COMPOSES MULTIPLICATIVELY with the VCC_ALPHA env var (default 1.0796), which "
            "scales the effect EARLIER, before whitening. Setting both multiplies them. "
            "For a board-calibrated amplitude change prefer VCC_ALPHA, whose scale the "
            "two-point board solve in experiments/096-mse-floor/ is expressed in; use "
            "--effect-scale only for a post-gate change VCC_ALPHA cannot express -- for "
            "instance re-raising the global amplitude after the shrinkage flags, whose "
            "compensation is SOLVED from the emitted cells at ~2.9-3.0x for the full stack "
            "(s 0.5576 -> 1.6956 on h1-test, 0.8330 -> 2.3748 on orion) and which VCC_ALPHA "
            "cannot undo linearly because the "
            "PerturbMatch and CD4 terms are added at fixed magnitude AFTER it. "
            "--shrink-key is the THIRD knob on this same quantity and multiplies with both "
            "of these; do not double-apply an amplitude cut. "
            "See experiments/096-mse-floor/."
        ),
    )
    parser.add_argument(
        "--shrink-key",
        default="none",
        choices=["none", "agreement", "agreement-pairwise"],
        help=(
            "Per-gene shrinkage of the gated log-space effect, keyed on a confidence "
            "score. 'none' (the default) applies no per-gene factor and is bit-identical "
            "to every prior build. "
            "SHIP THE PAIR OR NEITHER. The two halves are exposed as independent flags, "
            "but REALISED rule-1 minimums (min over the two harnesses, realised "
            "throughout -- never mixed with the analytic column) are 1.098 for "
            "--shrink-key alone, 1.176 for --shrink-pert-key alone and 1.473 for the "
            "stack: the stack is worth more than the sum of what either half can be "
            "trusted for, because each surface's weak half is the other's strong one "
            "(h1-test per-pert 1.246 > per-gene 1.098; orion per-gene 1.556 > per-pert "
            "1.176). A second, independent reason: realised eps RISES for the per- "
            "perturbation half alone (+0.00222, where the analytic route had it "
            "falling), so the claim*eps penalty works against it -- only the STACK "
            "keeps eps below the base on both surfaces. Do not read the stack as "
            "decomposable; a later arm that drops one half may be justified by the "
            "wrong harness. "
            "AND THE BINDING MEMBER IS nmae, WHICH NONE OF THIS EXPERIMENT WAS "
            "PRICING: nmae is an absolute fold-change error over the REFERENCE's DE "
            "genes, which are not the genes the donors agree on, so every shrinkage "
            "toward zero makes it monotonically worse. AND FOR THE SHIPPED ARM ITS SIGN "
            "DOES NOT TRANSFER: gated on a13 as built, nmae is +0.0026 on h1-test and "
            "-0.0145 on orion, so its contribution is not reliably signed under rule 1, "
            "and a13's worst-member floor (-0.0092) is WORSE than a milder arm's "
            "(-0.0012) even though its net is better on both surfaces. a13 is the better "
            "bet and NOT the safer one; if it underperforms on the board, nmae is the "
            "member to look at first, and its n_pred also drops 7.7%% on h1-test, the "
            "largest of any arm still net positive. The older monotonic trend reads "
            "h1-test "
            "0.9630 -> 0.9678 -> 0.9768 and orion 0.9189 -> 0.9369 -> 0.9585 across "
            "a12, a13 and a13 at gamma_pert=2. The recommended arm pays that and the "
            "other members cover it; a more aggressive one does not, so do NOT tune "
            "these flags on mse alone. From the same gate: n_pred stays within 2.4%% "
            "and fid rises on both surfaces against a 0.0001-0.0003 seed floor, while "
            "reach is UNREADABLE there -- re-seeded same-policy twins differ by 0.0207 "
            "(h1) and 0.0245 (orion), an order of magnitude above every between-arm "
            "difference. "
            "'agreement' keys on MultiDonorModel.donor_agreement -- "
            "the geometric mean of |delta| over the donors carrying the target, zeroed "
            "unless all of them agree on the sign, the same quantity --keep-key agreement "
            "ranks by -- and multiplies each gene's effect by exp(-gamma * rank / 2000), "
            "rank 0 being the gene the donors agree most strongly about. Experiment 096 "
            "measured this family at 1.17x (orion) to 1.38x (h1-test) on sqrt(W), positive "
            "on BOTH surfaces, where the same shrinkage with a per-gene slope FITTED on a "
            "donor gains 1.02x, i.e. nothing: a confidence score is a model property and "
            "transfers, a slope is a cell-line property and does not. "
            "COMPOSES MULTIPLICATIVELY with --effect-scale and with the VCC_ALPHA env var "
            "(default 1.0796). All three scale the same predicted effect -- VCC_ALPHA "
            "before whitening, the other two after the CPM gate -- so setting more than "
            "one MULTIPLIES them; do not double-apply an amplitude cut. This one cuts the "
            "panel-wide effect norm hard (at gamma=1.0 only ~1,386 of 18,533 genes keep "
            "more than half their amplitude) and the measurement re-optimises the global "
            "scale AFTER shrinking, so re-solve VCC_ALPHA rather than reusing the unshrunk "
            "arm's. Targets carried by a single donor have no cross-donor evidence and are "
            "left unshrunk. See experiments/096-mse-floor/ intervention 5, mech_miller.py "
            "and mech_shrink_target.py."
        ),
    )
    parser.add_argument(
        "--shrink-toward",
        default=None,
        choices=["mean", "zero"],
        help=(
            "UNVALIDATED DEFAULT -- read this before trusting it. Four knobs in this "
            "family were tuned against a SINGLE member, mse, and their defaults then "
            "read as though someone had checked them. AND WITHIN THIS FAMILY THE mse "
            "OBJECTIVE IS ANTI-CORRELATED WITH THE TOTAL AT THE AGGRESSIVE END, so "
            "optimising it without a gate actively selects the loser rather than merely "
            "failing to see it: --shrink-gamma 1.0 carries the HIGHEST sqrt(W) anywhere "
            "in this experiment (1.515 / 2.102) and a raw nmae of 1.0667 against that "
            "member's baseline anchor of ~1.0013 -- worse than predicting nothing on it. "
            "Three arms have now been recommended by mse and refused by a DE gate: "
            "--shrink-gamma's own value, --shrink-pert-gamma, and WHITEN_GAMMA. When "
            "experiment 097 finally ran a "
            "DE gate it covered exactly one of the four (--shrink-pert-gamma, whose 2.0 "
            "default was mse-optimal and came out NEGATIVE on orion) and it found that "
            "nmae, priced by nobody -- not this builder, not MSE-MECHANISM.md, not "
            "097 -- is the binding member on the whole family. This flag is one of the "
            "three the gate did NOT cover. Its default was chosen against mse alone (12/12 on the exact identity, both harnesses) and has "
            "never faced a DE gate; it is not known to be neutral on nmae, reach, fid, "
            "jac or n_pred. It is left where it is because 097's gate is the only DE "
            "evidence we have and it does not reach here -- changing it now would swap "
            "one unvalidated choice for another. "
            "What a fully shrunk gene converges to under --shrink-key agreement. 'zero' "
            "(the default when the key is on) shrinks toward the context's control profile, "
            "i.e. e' = a*e. 'mean' is the James-Stein form e' = mbar + a*(e - mbar), mbar "
            "being the grand mean predicted log effect per gene over all 300 perturbations of "
            "our own effect matrix -- no reference, no oracle leak, but measured worse; it is "
            "kept reachable only so the comparison can be rerun. "
            "MEASURED, PREFER 'zero': experiments/096-mse-floor/mech_shrink_target.py scores "
            "both through the scorer's own identity with the #348 budget B recomputed rather "
            "than assumed, and 'zero' wins on BOTH harnesses at every gamma with both "
            "agreement keys, 12/12 -- at gamma=1 with the shipped key, d(mse) is -0.00856 vs "
            "+0.00336 on h1-test and -0.05625 vs -0.01181 on orion (RE-MEASURED after the "
            "midrank/degenerate-key fix; the pre-fix run read -0.0089/+0.0029 and "
            "-0.0552/-0.0115, and every number moved by less than 0.001 without changing an "
            "ordering). 'mean' is worse than not "
            "shrinking at all on h1-test. The argument FOR 'mean' -- that the metric's 0-"
            "anchor is the context's mean perturbation response (board raw ~0.9879) rather "
            "than the control (team \"Kun Lab\"'s zero-effect submission, raw 1.0424) -- "
            "fails because our own panel mean does not estimate the REFERENCE's mean effect: "
            "cos is +0.016 (h1-test) and -0.066 (orion), no capture and not even a shared "
            "sign. It also inflates (A - B)/A from ~0.012 to ~0.055, and A - B is the term "
            "the master formula drops, which is why mech_miller.py's sqrt(W) headline "
            "PREFER 'agreement'. 'agreement-pairwise' averages the same key over donor PAIRS "
            "instead of demanding unanimity across all three; with exactly two donors the two "
            "are IDENTICAL, so no two-donor harness number separates them and only production's "
            "third donor does. Measured there -- experiment 097's pairwise3.py, on h1-test with "
            "orion_hct116 added as a legal third donor, which makes that triple exactly "
            "production's donor set -- unanimity WINS: multiplier 1.278 against pairwise's "
            "1.213 at eta=1, and 1.612 against 1.504 stacked with --shrink-pert-key. Both beat "
            "the 2-donor estimator's 1.101, so the production key is stronger than what the "
            "harnesses measured, not weaker. Corroborated on BOTH surfaces by graded_fraction."
            "py: dialling the graded gene fraction down from 0.21 to 0.05 raises the multiplier "
            "on each. CAVEAT, n=1 SURFACE: orion has no legal third donor, so rule 1's "
            "two-surface test cannot be run on this comparison and never will be with what we "
            "hold. 'agreement-pairwise' remains the documented fallback because its donor-count "
            "ordering is not backwards (unanimity grades 0.106 of genes on a 3-donor target "
            "against 0.545 on a 2-donor one, i.e. less evidence is damped less) -- that "
            "objection is outvoted on score, not refuted. NOTE --keep-key is deliberately NOT "
            "given this choice: it is board-validated on the frozen donor_agreement and is left "
            "alone. "
            "The --shrink-toward default was originally specified as 'mean' and was flipped to "
            "'zero' against that spec once the measurement came in -- reproduced independently "
            "by mech_shrink_exact.py, which disagrees on levels and agrees on the ordering. "
            "What is NOT established is the magnitude: both routes hold `claim` at its "
            "unshrunk value, and claim is a functional of the emitted cells, so re-solve it on "
            "the built .vcc with mech_challenge.py before quoting a board number. "
            "Rejected without --shrink-key agreement."
        ),
    )
    parser.add_argument(
        "--shrink-gamma",
        type=float,
        default=None,
        help=(
            "Decay exponent for --shrink-key agreement: alpha = exp(-gamma * rank / 2000). "
            "Unset with --shrink-key agreement uses 0.25 -- the value a13 was BUILT and GATED "
            "with, gated as the ARTIFACT at +0.0215 (h1-test) / +0.0040 (orion) net overall "
            "against a12 across all six members. It is NOT the mse optimum (096's sweep liked 1.0) and the gap is "
            "large: at gamma=1.0 the ~87%% of each row whose cross-donor agreement is exactly "
            "zero is damped to alpha 0.0054, against 0.271 at gamma=0.25 -- mean alpha 0.080 "
            "against 0.347. THIS IS THE ONLY KNOB IN THE SHRINKAGE FAMILY THAT pds CONSTRAINS: "
            "on orion pds falls MONOTONICALLY from gamma=0, so no positive value is free, and "
            "rule 1 binds on that surface even though h1-test peaks at 1.0. NOTE the family is "
            "exponential in an absolute RANK COUNT with a 2000-gene decay length, so it is NOT "
            "scale-free in the gene axis -- halving G roughly doubles mean alpha (measured: "
            "0.166 at G=9,000 against 0.080 at G=18,533 for gamma=1). A figure quoted from any "
            "other parameterisation, e.g. a power law on rank/G, does not transfer by any "
            "scalar. "
            "SO QUOTE MEAN ALPHA, NOT GAMMA. It is the one number that survives a "
            "reparameterisation, and 097 measured that it is also what actually drives the "
            "score: their power law at mean alpha 0.834 and this exponential at mean alpha "
            "0.800 -- different families entirely -- land within 0.0007 of each other on BOTH "
            "surfaces. The amount of shrinkage matters; the functional form barely does. "
            "MEASURED on the real 300-target panel, mean alpha by gamma: 0.796 at 0.05, "
            "0.347 at 0.25 (the default), 0.160 at 0.5, 0.077 at 1.0, 0.047 at 2.0. "
            "Quote a row of that table, not a gamma, whenever a number leaves this file. "
            "gamma=0 is an exact no-op. Rejected without --shrink-key agreement."
        ),
    )
    parser.add_argument(
        "--shrink-pert-key",
        default="none",
        choices=["none", "donor-cos"],
        help=(
            "PER-PERTURBATION shrinkage, multiplied into the same matrix as --shrink-key: "
            "z_p = rank01(mean pairwise cosine between the carrying donors' effect vectors "
            "for that target)^gamma, applied as z[:, None] * alpha. 'none' (the default) is "
            "bit-identical to every prior build. "
            "SHIP THE PAIR OR NEITHER. The two halves are exposed as independent flags, "
            "but REALISED rule-1 minimums (min over the two harnesses, realised "
            "throughout -- never mixed with the analytic column) are 1.098 for "
            "--shrink-key alone, 1.176 for --shrink-pert-key alone and 1.473 for the "
            "stack: the stack is worth more than the sum of what either half can be "
            "trusted for, because each surface's weak half is the other's strong one "
            "(h1-test per-pert 1.246 > per-gene 1.098; orion per-gene 1.556 > per-pert "
            "1.176). A second, independent reason: realised eps RISES for the per- "
            "perturbation half alone (+0.00222, where the analytic route had it "
            "falling), so the claim*eps penalty works against it -- only the STACK "
            "keeps eps below the base on both surfaces. Do not read the stack as "
            "decomposable; a later arm that drops one half may be justified by the "
            "wrong harness. "
            "AND THE BINDING MEMBER IS nmae, WHICH NONE OF THIS EXPERIMENT WAS "
            "PRICING: nmae is an absolute fold-change error over the REFERENCE's DE "
            "genes, which are not the genes the donors agree on, so every shrinkage "
            "toward zero makes it monotonically worse. AND FOR THE SHIPPED ARM ITS SIGN "
            "DOES NOT TRANSFER: gated on a13 as built, nmae is +0.0026 on h1-test and "
            "-0.0145 on orion, so its contribution is not reliably signed under rule 1, "
            "and a13's worst-member floor (-0.0092) is WORSE than a milder arm's "
            "(-0.0012) even though its net is better on both surfaces. a13 is the better "
            "bet and NOT the safer one; if it underperforms on the board, nmae is the "
            "member to look at first, and its n_pred also drops 7.7%% on h1-test, the "
            "largest of any arm still net positive. The older monotonic trend reads "
            "h1-test "
            "0.9630 -> 0.9678 -> 0.9768 and orion 0.9189 -> 0.9369 -> 0.9585 across "
            "a12, a13 and a13 at gamma_pert=2. The recommended arm pays that and the "
            "other members cover it; a more aggressive one does not, so do NOT tune "
            "these flags on mse alone. From the same gate: n_pred stays within 2.4%% "
            "and fid rises on both surfaces against a 0.0001-0.0003 seed floor, while "
            "reach is UNREADABLE there -- re-seeded same-policy twins differ by 0.0207 "
            "(h1) and 0.0245 (orion), an order of magnitude above every between-arm "
            "difference. "
            "Experiment 097 measured the "
            "family, REALISED, at 1.246 (h1-test) / 1.176 (orion) alone and 1.473 / 1.580 "
            "stacked on --shrink-key "
            "agreement, against a permutation null over the same grid and argmax at "
            "p <= 0.001, and against 1.035 for experiment 096's best per-perturbation key. It "
            "asks a different question from --shrink-key: that one is a per-GENE magnitude "
            "gated on sign, this is a whole-vector DIRECTION agreement for the perturbation "
            "as a unit, so a per-row scalar cannot reorder genes within a row and the DE "
            "membership of a row is untouched. Buildable on 294/300 challenge targets "
            "(k562_gwps carries 269, orion_hct116 290, orion_hek293t 294); the other 6 have "
            "fewer than two donors and are left at z = 1. COMPOSES MULTIPLICATIVELY with "
            "--shrink-key, --effect-scale and VCC_ALPHA. See --shrink-pert-floor for the "
            "amplitude hazard, and experiments/097-direction/."
        ),
    )
    parser.add_argument(
        "--shrink-pert-gamma",
        type=float,
        default=None,
        help=(
            "Exponent for --shrink-pert-key: z_p = rank01(key)^gamma. Unset with the key on "
            "uses 1.0. NOT the mse optimum -- 2.0 is, at 1.478 (h1-test) / 1.610 (orion) -- but "
            "097's DE gate prices all six members and gp=2 comes out NEGATIVE on orion "
            "(-0.0037 overall against gp=1's +0.0025), with nmae alone accounting for -0.0165. "
            "Rule 1 binds on the worse surface. Raising this past 1.0 buys mse and spends more "
            "than it buys on nmae; see the SHIP THE PAIR warning for the mechanism. gamma=3 is "
            "worse again and also leaves 111 of 300 targets under z=0.05 even floored. "
            "gamma=0 is an exact no-op. "
            "THIS KNOB IS EXACTLY FREE FOR pds, PROVIDED --shrink-pert-floor > 0: a "
            "per-perturbation weight is a per-row scalar and pds_cosine is scale-invariant per "
            "row, so 097 measured pds bit-identical across gamma_pert = 0, 1, 2 at every "
            "gamma_gene on both surfaces. Only the per-gene gamma is constrained by pds. That "
            "guarantee depends on the floor: at floor 0 the weakest rows go near-zero-norm and "
            "stop being a scalar rescale of anything. Rejected "
            "without --shrink-pert-key donor-cos."
        ),
    )
    parser.add_argument(
        "--shrink-pert-floor",
        type=float,
        default=None,
        help=(
            "UNVALIDATED DEFAULT -- read this before trusting it. Four knobs in this "
            "family were tuned against a SINGLE member, mse, and their defaults then "
            "read as though someone had checked them. AND WITHIN THIS FAMILY THE mse "
            "OBJECTIVE IS ANTI-CORRELATED WITH THE TOTAL AT THE AGGRESSIVE END, so "
            "optimising it without a gate actively selects the loser rather than merely "
            "failing to see it: --shrink-gamma 1.0 carries the HIGHEST sqrt(W) anywhere "
            "in this experiment (1.515 / 2.102) and a raw nmae of 1.0667 against that "
            "member's baseline anchor of ~1.0013 -- worse than predicting nothing on it. "
            "Three arms have now been recommended by mse and refused by a DE gate: "
            "--shrink-gamma's own value, --shrink-pert-gamma, and WHITEN_GAMMA. When "
            "experiment 097 finally ran a "
            "DE gate it covered exactly one of the four (--shrink-pert-gamma, whose 2.0 "
            "default was mse-optimal and came out NEGATIVE on orion) and it found that "
            "nmae, priced by nobody -- not this builder, not MSE-MECHANISM.md, not "
            "097 -- is the binding member on the whole family. This flag is one of the "
            "three the gate did NOT cover. Its default was chosen against mse (097 measured it near-free, 1.467 -> 1.478) plus an n_pred/pds argument, never a scored arm and has "
            "never faced a DE gate; it is not known to be neutral on nmae, reach, fid, "
            "jac or n_pred. It is left where it is because 097's gate is the only DE "
            "evidence we have and it does not reach here -- changing it now would swap "
            "one unvalidated choice for another. "
            "Lower clamp applied to the RANK, INSIDE the exponent: z = max(rank01, floor)"
            "**gamma, so 0.3 at gamma=2 bounds every target at 0.09 of full amplitude, not "
            "0.3. That is experiment 097's construction and its grid is reported in those "
            "units; clamping after the exponent instead would silently disagree with every "
            "number they quote. Default 0.3. Unfloored, rank01 puts the weakest target at "
            "0.5/n -- at 300 targets it keeps 0.17%% of its effect and emits essentially the "
            "bare control profile, which is not a crash but is a perturbation with no DE calls "
            "(an n_pred loss the board reads) and a near-zero-norm row for pds; at gamma=2 "
            "that was 66 of 300 targets under z=0.05. 097 measured the floor as nearly free on "
            "mse (1.467 -> 1.478 at gamma=2), so the guard costs nothing. "
            "DO NOT LOWER IT BELOW 0.18 WITHOUT A NEW MEASUREMENT. That is the rank a 0.0 "
            "donor cosine earns on the real panel, and above it the --shrink-pert-missing "
            "policy is unobservable -- every branch floors to the same weight, which is why "
            "this builder and 097's gate produce bit-identical artifacts at 0.30 despite "
            "implementing different policies. Below 0.18 the policy choice becomes live and is "
            "unmeasured on every surface we hold, so a floor cut is two changes, not one. "
            "IT IS ALSO LOAD-BEARING FOR pds, WHICH IS THE STRONGER REASON. A near-zero-norm "
            "predicted effect ties its whole distance row, and metrics.discrimination's own "
            "docstring records that such a row scores exactly 0.5 under the midrank tie policy "
            "against a per-target baseline near 0.833. Unfloored at gamma=2 the panel minimum "
            "z is 2.8e-06, i.e. the bottom-ranked CARRIED targets are the ones driven to a "
            "degenerate row -- NOT the fewer-than-two-donor targets, which under 'damp' land at "
            "a MID rank on the challenge panel (z = 0.111 at gamma=2) because challenge cosines "
            "run [-0.285, +0.439] -- min/max over the 294 CARRIED targets, p05/p95 being "
            "-0.043/+0.115 -- with median +0.0284. floor 0.3 bounds every row at 0.09 and "
            "097 verified pds is then bit-identical across gamma_pert. Lowering this to 0 to be "
            "'more aggressive' buys nothing on mse and silently spends pds. Rejected without "
            "--shrink-pert-key."
        ),
    )
    parser.add_argument(
        "--shrink-pert-missing",
        default=None,
        choices=["damp", "neutral", "keep"],
        help=(
            "UNVALIDATED DEFAULT -- read this before trusting it. Four knobs in this "
            "family were tuned against a SINGLE member, mse, and their defaults then "
            "read as though someone had checked them. AND WITHIN THIS FAMILY THE mse "
            "OBJECTIVE IS ANTI-CORRELATED WITH THE TOTAL AT THE AGGRESSIVE END, so "
            "optimising it without a gate actively selects the loser rather than merely "
            "failing to see it: --shrink-gamma 1.0 carries the HIGHEST sqrt(W) anywhere "
            "in this experiment (1.515 / 2.102) and a raw nmae of 1.0667 against that "
            "member's baseline anchor of ~1.0013 -- worse than predicting nothing on it. "
            "Three arms have now been recommended by mse and refused by a DE gate: "
            "--shrink-gamma's own value, --shrink-pert-gamma, and WHITEN_GAMMA. When "
            "experiment 097 finally ran a "
            "DE gate it covered exactly one of the four (--shrink-pert-gamma, whose 2.0 "
            "default was mse-optimal and came out NEGATIVE on orion) and it found that "
            "nmae, priced by nobody -- not this builder, not MSE-MECHANISM.md, not "
            "097 -- is the binding member on the whole family. This flag is one of the "
            "three the gate did NOT cover. Its default was chosen against a structural argument plus n=4 perturbations on one surface and has "
            "never faced a DE gate; it is not known to be neutral on nmae, reach, fid, "
            "jac or n_pred. It is left where it is because 097's gate is the only DE "
            "evidence we have and it does not reach here -- changing it now would swap "
            "one unvalidated choice for another. "
            "What to do with a target fewer than two donors carry, whose donor cosine is "
            "UNDEFINED rather than zero. 'damp' (the default when the key is on) scores it 0.0 "
            "and ranks it with everything else -- experiment 097's recommendation. 'keep' "
            "gives it z = 1; it was the original spec on the grounds of parity with the "
            "per-gene policy, and that parity is FALSE: alpha = 1 on a degenerate per-gene row "
            "is a true no-op, but z = 1 is not, because z is RELATIVE -- every other row is "
            "pushed down, so the rows we know least about become the highest-weighted on the "
            "panel. 097 measured 'keep' at 1.023 (gamma=1) and 0.958 (gamma=2) against 'damp' "
            "at 1.207 / 1.239, i.e. worse than not shrinking; those 4 h1-test targets are "
            "13.1%% of predicted amplitude A and +0.5%% of X. TREAT THAT MAGNITUDE WITH CARE: "
            "it is 4 perturbations on ONE surface, orion is 132/132 two-donor and cannot see "
            "it at all, and the quoted rho +0.269 between n_donors and the true cosine is a "
            "function of the same 4 rows, so it is not independent corroboration. The "
            "PRINCIPLE is what justifies the default, not the number. 'neutral' gives them the "
            "median carried weight -- what the principle alone implies, absence of evidence "
            "being neither confidence nor its opposite -- and nobody has measured it. Note "
            "'damp' is not 'bottom' on the real panel: challenge cosines run [-0.285, +0.439] "
            "(min/max over the 294 CARRIED targets; quote p05/p95 -0.043/+0.115 if you want a "
            "robust spread, and never percentile a 0-filled 300) "
            "with median +0.0284, so a score of 0.0 is a mid-rank there. "
            "IT DIVERGES FROM 097'S GATE IN CODE AND NOT IN THE ARTIFACT, WHICH IS WORTH "
            "KNOWING BEFORE ANYONE LOWERS THE FLOOR. This function scores an uncarried target "
            "0.0 and ranks it with the panel; 097's gate instead discards that rank and forces "
            "the uncarried WEIGHT to 0.0 before flooring. Two different policies -- but "
            "MEASURED on the real 300-target panel they are bit-identical on all 300 rows at "
            "a13's settings, because only 17.3%% of donor cosines are negative, so a 0.0 cosine "
            "ranks at 0.180, and max(0.180, FLOOR) == max(0.0, FLOOR) for any FLOOR >= 0.18. "
            "a13 ships at 0.30. The divergence is LATENT: it goes live only below FLOOR 0.18, "
            "and is unmeasured on every surface we hold. See --shrink-pert-floor. NOTE this "
            "was first reported here as a live 294/6 split; that was computed on a SIMULATED "
            "cosine distribution with far more negative mass than the real donors have, and "
            "running the actual population through both branches returned zero differing rows. "
            "Rejected without --shrink-pert-key."
        ),
    )
    parser.add_argument("--spread-lo", type=float, default=0.0,
                        help="graded map: spread at the most confident gene.")
    parser.add_argument("--spread-hi", type=float, default=0.60,
                        help="graded map: spread at the least confident gene.")
    parser.add_argument(
        "--whitener-id",
        default="union",
        choices=sorted(WHITENER_REGISTRY),
        help=(
            "Which registered whitening operator this build uses. 'union' is v11's "
            "production basis; 'union_pruned' is experiment 091's, with the 4,393 "
            "common-essential rows dropped. The id selects which sha256 pair the "
            "frozen-input guard demands, so an unregistered or drifted operator still "
            "aborts before any model input is read."
        ),
    )
    parser.add_argument(
        "--keep-spread",
        type=float,
        default=None,
        help=(
            "Emission spread applied to the kept top-K genes, in [0, 1]. Unset uses the "
            "shipped scalar SPREAD=0.25, which is what v14 and v15 both used. Setting it "
            "BELOW --tail-spread lowers those genes' Wilcoxon p-values instead of only "
            "raising the tail's: n_pred rises rather than falls, which pays fid's "
            "coverage term min(1, n_pred/n_real) that v14 and v15 both underpaid, while "
            "still separating the informative genes from the skewness artefact in the "
            "confidence order that reach walks. Requires --keep-top or --keep-top-frac "
            "(experiment 094)."
        ),
    )
    parser.add_argument(
        "--keep-top-frac",
        type=float,
        default=None,
        help=(
            "Scale-adaptive form of --keep-top: a fraction in [0, 1] of the genes that "
            "pass the CPM gate. Mutually exclusive with --keep-top."
        ),
    )
    parser.add_argument(
        "--emitter",
        default="quasi",
        choices=list(EMITTERS),
        help=(
            "The emission LAW. 'quasi' (the default) is the shipped exact-total quasi "
            "emitter and writes no meta key, so an unflagged build is bit-identical to "
            "every prior build. "
            "'rc_hybrid' is experiment 185 section 7's rank-controlled hybrid: a TWO-PASS "
            "emitter that emits the incumbent's own block under the quasi law, reads its "
            "per-gene Wilcoxon z against the context's real control cells, and then "
            "re-assigns THAT |z| multiset inside the incumbent's own BH call set by "
            "descending |H| and inside its complement by descending |H| (the `abs` sign "
            "policy), realising it with an exact-marginal emitter that holds every gene's "
            "panel total and welds every cell's library to T = S // 400. So the delivered "
            "pseudobulk is the incumbent's count for count (`pds` and `mse`'s profile half "
            "are identities), the call COUNT is the incumbent's, `nmae` is pinned to five "
            "decimals, and the p-value ORDER is the only thing that moves: 185 measured "
            "Spearman(-log p, |H|) inside the call set going -0.391 -> +0.984 and h1-test "
            "`reach` 0.18283 -> 0.21681 at `jac` 0.06824 -> 0.06827. It needs the unshrunk "
            "effect (--effect-scale 1.0, no --shrink-key) and refuses VCC_MASS_FIX. "
            "'rc_hist_a' is the SAME emitter with the re-key applied to the whole gene "
            "axis instead of inside the call set and its complement (185 section 3a), so "
            "the call set is re-drawn by |H| at the incumbent's own call count: h1-test "
            "`reach` 0.23278 at `jac` 0.05932, i.e. it buys `reach` and spends `jac`. "
            "That trade is the one variable between the two and it is what experiment 201 "
            "prices; every other setting, including the four rc_* meta keys, is shared. "
            "'rc_callset' is 'rc_hybrid' with ONE variable, the call SET (experiment 202, "
            "riding experiment 192 section 5.3): the block that goes in front is the "
            "PRESCRIBED top 400 by |H| UNION top 2,200 by control CPM (n_pred 2,600 "
            "against `rc_hybrid`'s 4,480) instead of the incumbent's own BH set, with the "
            "inside it and inside its complement unchanged. Control CPM is the per-gene "
            "mean of per-cell CPM over that context's own NTC control cells -- measured, "
            "and the same object the scorer gates on. 192 measured, on h1-test at the "
            "retrodicted n_real level, `fid` +0.00097, `jac` +0.01288, `reach` +0.01123 "
            "raw marginal to `rc_hybrid`: the detectability half buys membership (`jac`), "
            "the protected |H| head keeps the order the reach walk reads. Because BH's "
            "count is a function of the p-value MULTISET, a smaller set needs the multiset "
            "to move; it moves in exactly one place, a single common factor putting the "
            "complement's largest target |z| on z_cap 1.90 (p = 0.0574 > alpha), which is "
            "strictly monotone and so leaves |H|'s order in the complement exactly where "
            "`rc_hybrid` put it. The three numbers are the arm, so they are written into "
            "the config as rc_callset_head / rc_callset_cpm / rc_callset_z_cap. "
            "'resample_mix' keeps the top --resample-keep-top genes on "
            "the quasi law and emits the REST as resampled real control cells carrying "
            "the predicted fold change. "
            "WHY (experiment 177): the shipped tail at --tail-spread 1.0 --overdispersion "
            "4.0 calls 43.3%% / 44.8%% / 38.4%% of the DE universe significant at ZERO "
            "predicted effect on challenge A/B/C, ~100%% of them DOWN, because 400 cells "
            "at one fixed library size sit on a CPM lattice the real controls do not; the "
            "resampled law calls 0.0-0.2%% (1 / 18 / 19 genes) because real control cells "
            "carrying a fold change are EXCHANGEABLE with the reference controls, so "
            "Fbar = 1/2 identically rather than asymptotically. Realised on h1-test at "
            "K = 8000: `reach` +0.0352 raw (69%% of INBOX-4's full |H| ordering "
            "counterfactual), rho(p_adj, |H|) inside our own call set -0.395 -> +0.485, "
            "rho(p_adj, expression) +0.484 -> +0.249, precision +0.0052, `nmae` -0.0002. "
            "The delivered per-gene panel totals are UNCHANGED -- `emit_counts` repairs "
            "them exactly after the splice -- so `pds` and `mse`'s profile half are "
            "pinned by identity, verified as a sha256 on 85/85 targets in 177 and on "
            "900/900 (context, perturbation) pairs in 179. "
            "WHAT IT COSTS: `n_pred` falls 6066 -> 3304 on h1-test (board ~4424 -> "
            "~2410), so `fid`'s yield term min(1, n_pred/n_real) and `jac`'s union both "
            "come off saturation and the arm's sign depends on the challenge's n_real, "
            "which AGENTS.md 11.1 rule 9 says no harness we hold can read. 177 section 6 "
            "brackets it: +0.006 overall at n_real ~240, -0.173 at 093's ~2670. "
            "TWO CALL-SITE CONTRACTS, or the arm INVERTS: the donated cells must be RAW "
            "counts with their OWN library sizes (never rescaled to a nominal depth -- "
            "the exact-total repair then peels ~80%% of every gene's column and called "
            "10,658 of 10,779 genes significant in 177's reproduction), and the control "
            "profile must be the POOLED one, expm1(b_ctrl)/sum, not the mean-per-cell "
            "CPM profile. Both are enforced here; see assert_raw_control_counts."
        ),
    )
    parser.add_argument(
        "--resample-keep-top",
        type=int,
        default=None,
        help=(
            "Number of genes per target held on the parametric quasi law under --emitter "
            "resample_mix, ranked by the SAME score --keep-top uses over the same CPM-"
            "gated axis. Unset defaults to --keep-top, which makes the parametric half "
            "exactly the incumbent's keep block so the ONLY thing that changes is the "
            "tail's emission law. Experiment 177 section 4 measured the K sweep on "
            "h1-test: 2000 -> `reach` +0.0523 at n_pred 1579, 4000 -> +0.0559 at 2307, "
            "6000 -> +0.0350 at 3020, 8000 -> +0.0352 at 3304. Lower K buys more `reach` "
            "and cuts `n_pred` harder, i.e. it bets more of the arm on a small n_real; "
            "K = 8000 is the most n_real-robust point, staying net positive to "
            "n_real ~1100. Requires --emitter resample_mix."
        ),
    )
    return parser


def main(argv: Sequence[str] | None = None) -> dict[str, object]:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.repin_frozen_inputs and args.repin_panel_support:
        parser.error(
            "--repin-frozen-inputs and --repin-panel-support are separate steps: re-pin the "
            "file hashes first, then pin panel_support from one effects build"
        )
    if args.repin_frozen_inputs:
        repin_frozen_inputs(
            args.bundle.expanduser().resolve(), args.frozen_inputs.expanduser().resolve()
        )
        return {"status": "frozen_inputs_repinned", "submission": False}
    if args.output is None and not args.repin_panel_support:
        parser.error("-o/--output is required for a build")
    return run(args)


if __name__ == "__main__":
    main()
