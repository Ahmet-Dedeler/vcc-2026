"""Apply the 410b BALLAST path to the builder (run ONLY after run410.sh logs "builder after restore").

Four exact edits, each asserted to match once:
  1. module-level knob + helpers (work/ballast_block.py.txt) before DE_SIGN_PRIOR_BY_TARGET
  2. write_submission_from_effects: `bulk_effects` kwarg + validation
  3. the DE-sign-prior RC branch: carriers on n - K cells, then ballast rows
  4. run(): with VCC_PERT_SCALE_BALLAST the writer gets the UNSCALED matrix as `effects` and the scaled one as `bulk_effects`;
     the config records pert_scale_ballast
usage: uv run python patch_ballast410.py [--write]
"""
import difflib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
import os
P = Path(os.environ.get("PATCH_TARGET", ROOT / "scripts/make_submission_pm_cd4.py"))
BLOCK = (Path(__file__).resolve().parent / "work/ballast_block.py.txt").read_text()
src = P.read_text()
out = src


def sub(old, new, label):
    global out
    n = out.count(old)
    if n != 1:
        raise SystemExit(f"[{label}] expected exactly one match, found {n}")
    out = out.replace(old, new)


# 1. helpers
sub("DE_SIGN_PRIOR_BY_TARGET: dict[str, np.ndarray] = {}   # filled by de_sign_prior_full for the per-target knobs",
    BLOCK + "DE_SIGN_PRIOR_BY_TARGET: dict[str, np.ndarray] = {}   # filled by de_sign_prior_full for the per-target knobs",
    "helpers")

# 2. kwarg + validation
sub("""    context_effects: Mapping[str, np.ndarray] | None = None,
) -> dict[str, object]:
    \"\"\"Stream every context's raw counts with the exact v9 RNG and emission semantics.""",
    """    context_effects: Mapping[str, np.ndarray] | None = None,
    bulk_effects: np.ndarray | None = None,
) -> dict[str, object]:
    \"\"\"Stream every context's raw counts with the exact v9 RNG and emission semantics.""",
    "kwarg")
sub("""    effects = np.asarray(effects, dtype=np.float64)
    if effects.shape != (len(spec.perts), len(spec.genes)):""",
    """    effects = np.asarray(effects, dtype=np.float64)
    if bulk_effects is not None:
        # 410b: the BALLAST path (VCC_PERT_SCALE_BALLAST). `effects` is the DE-channel design (carriers), `bulk_effects` the
        # per-target scaled design the block's per-gene totals are allocated to.
        bulk_effects = np.asarray(bulk_effects, dtype=np.float64)
        if bulk_effects.shape != effects.shape:
            raise ValueError(f"bulk_effects shape {bulk_effects.shape} != {effects.shape}")
        if emitter not in RC_EMITTERS or DE_SIGN_PRIOR is None or context_effects is not None:
            raise ValueError("bulk_effects (ballast) needs an RC emitter, VCC_DE_SIGN_PRIOR and one shared effect matrix")
    if effects.shape != (len(spec.perts), len(spec.genes)):""",
    "validate")

# 3. the RC + DE-sign-prior branch
sub("""                    if DE_SIGN_PRIOR is not None:
                        from vcc2026 import rc_emit_k  # noqa: PLC0415 -- only on the arm path
                        block, rc_diag = rc_emit_k.rc_hybrid_block_k(
                            profile,
                            libraries[picks],""",
    """                    if DE_SIGN_PRIOR is not None:
                        from vcc2026 import rc_emit_k  # noqa: PLC0415 -- only on the arm path
                        carrier_libs = libraries[picks]
                        if bulk_effects is not None:
                            profile_bulk = profile_from_effect(
                                bulk_effects[target_index],
                                b_ctrl,
                                min_cpm=min_cpm,
                                effect_scale=effect_scale,
                                gene_scale=(
                                    None if shrink_scales is None else shrink_scales[target_index]
                                ),
                                shrink_center=shrink_center,
                            )
                            ballast_k, ballast_ratio = ballast_count(
                                profile, profile_bulk, libraries[picks], int(spec.cells_per_pert))
                            carrier_libs = libraries[picks][: int(spec.cells_per_pert) - ballast_k]
                        block, rc_diag = rc_emit_k.rc_hybrid_block_k(
                            profile,
                            carrier_libs,""",
    "branch-head")
sub("""                        if rc_emit_k.LAST_PRIOR_REC:
                            rc_diag = {**rc_diag, **rc_emit_k.LAST_PRIOR_REC.pop()}""",
    """                        if rc_emit_k.LAST_PRIOR_REC:
                            rc_diag = {**rc_diag, **rc_emit_k.LAST_PRIOR_REC.pop()}
                        if bulk_effects is not None:
                            block, ballast_diag = ballast_rows(block, profile, profile_bulk, ballast_k,
                                                               seed=RC_JITTER_SEED_BASE + row // spec.cells_per_pert)
                            rc_diag = {**rc_diag, **ballast_diag, "ballast_ratio": float(ballast_ratio)}
                            if block.shape[0] != int(spec.cells_per_pert):
                                raise AssertionError(f"{context}/{target}: ballast block has {block.shape[0]} rows")""",
    "branch-tail")

# 4. run(): the scaled matrix goes to the bulk, the unscaled one to the carriers
sub("""        effects, effect_audit["pert_scale"] = apply_pert_scale(effects, spec.perts, spec.genes)""",
    """        effects_unscaled = effects
        effects, effect_audit["pert_scale"] = apply_pert_scale(effects, spec.perts, spec.genes)
        if PERT_SCALE_BALLAST is not None:
            effect_audit["pert_scale"]["ballast"] = {"lb": BALLAST_LB, "kmax": BALLAST_KMAX}""",
    "run-scale")
sub("""    output_audit = write_submission_from_effects(
        spec,
        effects,
        bundle,
        output,
        context_effects=context_effects,""",
    """    output_audit = write_submission_from_effects(
        spec,
        effects_unscaled if PERT_SCALE_BALLAST is not None else effects,
        bundle,
        output,
        context_effects=context_effects,
        bulk_effects=effects if PERT_SCALE_BALLAST is not None else None,""",
    "run-writer")
sub("""            **({} if PERT_SCALE is None else {"pert_scale": PERT_SCALE, "pert_scale_sha256": PERT_SCALE_ARTIFACTS[PERT_SCALE][1]}),""",
    """            **({} if PERT_SCALE is None else {"pert_scale": PERT_SCALE, "pert_scale_sha256": PERT_SCALE_ARTIFACTS[PERT_SCALE][1]}),
            **({} if PERT_SCALE_BALLAST is None else {"pert_scale_ballast": {"lb": BALLAST_LB, "kmax": BALLAST_KMAX}}),""",
    "config")

# 5. register the uncapped artifact (410c candidate): same key, clip(0.3, inf)
sub("""        "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e",
    ),
}""",
    """        "d7f76a22f62ad6991d7220ce8a65218db38b2d71dcfe3c3781ec6a8f70409e0e",
    ),
    "wshare_lo0.3_inf": (
        ROOT / "experiments/410-pertscale-revisit/work/pert_scale_410_lo0.3_hi1e+09.npz",
        "94c048e99aecd40004e70f559fe26bdd3b8b9aef9cb8c8cd0efbf8a72d2fe69e",
    ),
}""",
    "artifact-inf")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
