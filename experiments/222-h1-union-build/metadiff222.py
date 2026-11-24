"""222 falsifier 1 — the rule-10 audit diff, `h1union` against the incumbent `cd4f034`.

`experiments/210-cd4-build/metadiff210.py` with the two pre-declared sets rewritten for this
arm and **declared before the diff was read**.  210's arm changed the CD4 dose, so its
`effects.cd4.*` had to move in exactly the dose's consequences.  THIS arm must not move CD4's
dose at all: `cd4_target_f` stays 0.34, `effects.cd4.fraction` stays bit-identical, and the
only thing in CD4's own block that may move is the **shared preserve scalar** -- because 218's
construction merges the two residuals into ONE add followed by production's ONE global
rescale, so the scalar is by definition a property of both.  Everything else in `effects.cd4`
is pinned by name below.

Every difference is printed, not only the expected ones (mandate).

    uv run python metadiff222.py <arm.meta.json> <cd4f034.meta.json>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# ---------------------------------------------------------------------------------------
# PRE-DECLARED BEFORE THE DIFF WAS READ (rule 10).
#
# CONFIG: the arm sets one environment variable, `VCC_H1_UNION_G`, and the builder writes
# that key plus the two it names alongside it (the source list and the sizing rule) only when
# the flag is on -- the `mass_fix` / `sign_repair` / `emitter` convention.  So exactly three
# config keys may APPEAR and none may change.
EXPECTED_CONFIG_ADDED = {"h1_union_g", "h1_union_sources", "h1_union_sizing_rule"}
EXPECTED_CONFIG_CHANGED: set[str] = set()
#
# EFFECTS/OPERATOR: `VCC_H1_UNION_G` is consumed in exactly one place -- a second residual
# solved by `solve_carrying_row_fraction` and merged into `add_post_whiten_residual`'s single
# add.  So the only audit fields that may MOVE are the sha of the final effect matrix (which
# must move, or the residual never applied) and the two scalars of the shared preserve.  The
# whole `effects.h1_union` block is new.
EXPECTED_AUDIT_CHANGED = {
    "effects.final_sha256",
    "effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm",
    "effects.cd4.post_add_preserve.post_add_global_preserve_scale",
}
EXPECTED_AUDIT_ADDED_PREFIX = "effects.h1_union."
#
# The scalars that must NOT move, named individually so a silent drift cannot hide in the
# aggregate count (AGENTS.md rule 10: the cheap form of a15's alpha falsifier).
PINNED = [
    "config.alpha",
    "config.cd4_target_f",
    "config.cd4_fraction",
    "config.cd4_mode",
    "config.cd4_conditions",
    "config.cd4_target_policy",
    "config.whitener_id",
    "config.whiten_gamma",
    "config.whiten_lambda",
    "config.emitter",
    "config.min_cpm",
    "config.keep_top",
    "config.keep_spread",
    "config.tail_spread",
    "config.overdispersion",
    "config.spread_map",
    "config.effect_scale",
    "config.emission_seed",
    "config.rc_key",
    "config.rc_callset_head",
    "config.rc_callset_cpm",
    "config.rc_callset_z_cap",
    "operator.sha256",
    "operator.metadata_sha256",
    "operator.columns_sha256",
    "operator.columns",
    "operator.all_panel_coordinates_excluded",
    "effects.baseline_pre_sha256",
    "effects.candidate_pre_sha256",
    "effects.post_whiten_sha256",
    "effects.whitening.after_all_median_norm",
    "effects.whitening.before_all_median_norm",
    "effects.whitening.after_operator_median_norm",
    "effects.whitening.before_operator_median_norm",
    "effects.whitening.operator_preserve_scale",
    "effects.whitening.gamma",
    "effects.whitening.lambda",
    # CD4's own dose, form and SOLVED FRACTION -- the single-variable constraint.
    "effects.cd4.fraction",
    "effects.cd4.mode",
    "effects.cd4.fully_supported_targets",
    "effects.cd4.supported_entries",
    "effects.cd4.panel_residual_nonzero_entries",
    "effects.cd4.residual_frame.solved_fraction",
    "effects.cd4.residual_frame.target_f_requested",
    "effects.cd4.residual_frame.f_per_unit_fraction",
    "effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction",
    "effects.cd4.residual_frame.whitened_base_median_row_norm",
    "effects.cd4.residual_frame.cd4_fraction_constant",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction",
    "effects.cd4.post_add_preserve.fraction",
    "effects.cd4.post_add_preserve.baseline_output_median_norm",
    "frozen_bundle.sha256",
    "frozen_bundle.panel_id",
    "output.n_rows",
]
# Reported, never asserted: these are properties of THIS build rather than of the arm.
VOLATILE_PREFIX = ("output.", "peak_rss_gib", "runtime.")
# The one implementation hash a code change is allowed to move.
EXPECTED_IMPLEMENTATION_CHANGED = {
    "hashes.implementation.scripts/make_submission_pm_cd4.py"
}


def flatten(obj, prefix=""):
    out = {}
    if isinstance(obj, dict):
        for key, value in obj.items():
            out.update(flatten(value, f"{prefix}.{key}" if prefix else str(key)))
    elif isinstance(obj, list) and obj and isinstance(obj[0], (dict, list)):
        for index, value in enumerate(obj):
            out.update(flatten(value, f"{prefix}[{index}]"))
    else:
        out[prefix] = obj
    return out


def main() -> None:
    arm = flatten(json.loads(Path(sys.argv[1]).read_text()))
    inc = flatten(json.loads(Path(sys.argv[2]).read_text()))
    keys = sorted(set(arm) | set(inc))

    added, removed, changed = [], [], []
    for key in keys:
        if key not in inc:
            added.append(key)
        elif key not in arm:
            removed.append(key)
        elif arm[key] != inc[key]:
            changed.append(key)

    def bucket(key: str) -> str:
        if key.startswith(VOLATILE_PREFIX):
            return "E-output/runtime"
        if key.startswith("config."):
            return "A-config"
        if key.startswith(("effects.", "operator.", "donor_coverage.",
                           "main_source_support.")):
            return "B-effects/operator"
        if key.startswith(("hashes.inputs", "frozen_bundle.")):
            return "C-frozen inputs"
        if key.startswith("hashes.implementation"):
            return "D-implementation"
        return "F-other"

    print("=" * 100)
    print("EVERY DIFFERENCE, not only the expected ones")
    print("=" * 100)
    for label, items in (("ADDED", added), ("REMOVED", removed), ("CHANGED", changed)):
        for key in items:
            tag = bucket(key)
            if label == "CHANGED":
                print(f"[{tag:20s}] {label:7s} {key}")
                print(f"{'':32s}  cd4f034: {inc[key]!r}"[:200])
                print(f"{'':32s}  h1union: {arm[key]!r}"[:200])
            else:
                value = arm.get(key, inc.get(key))
                print(f"[{tag:20s}] {label:7s} {key} = {value!r}"[:200])
    print()

    failures = []
    cfg_added = {k[len("config."):] for k in added if k.startswith("config.")}
    cfg_added = {k.split("[")[0] for k in cfg_added}
    cfg_changed = {k[len("config."):] for k in changed if k.startswith("config.")}
    cfg_removed = {k[len("config."):] for k in removed if k.startswith("config.")}
    print(f"A. CONFIG added   {sorted(cfg_added)}\n   expected       "
          f"{sorted(EXPECTED_CONFIG_ADDED)}")
    print(f"   CONFIG changed {sorted(cfg_changed)}   removed {sorted(cfg_removed)}")
    if cfg_added != EXPECTED_CONFIG_ADDED:
        failures.append(f"config additions {sorted(cfg_added)} != pre-declared")
    if cfg_changed != EXPECTED_CONFIG_CHANGED or cfg_removed:
        failures.append(f"config changed {sorted(cfg_changed)} removed {sorted(cfg_removed)}")

    aud_changed = {k for k in changed if bucket(k) == "B-effects/operator"}
    aud_added = {k for k in added if bucket(k) == "B-effects/operator"}
    aud_removed = {k for k in removed if bucket(k) == "B-effects/operator"}
    unexpected_added = {k for k in aud_added
                        if not k.startswith(EXPECTED_AUDIT_ADDED_PREFIX)}
    print(f"\nB. AUDIT changed {sorted(aud_changed)}\n   expected      "
          f"{sorted(EXPECTED_AUDIT_CHANGED)}")
    print(f"   AUDIT added: {len(aud_added)} keys, "
          f"{len(aud_added) - len(unexpected_added)} of them under "
          f"`{EXPECTED_AUDIT_ADDED_PREFIX}`; unexpected {sorted(unexpected_added)}")
    if aud_changed != EXPECTED_AUDIT_CHANGED:
        failures.append(f"audit changes {sorted(aud_changed)} != pre-declared")
    if unexpected_added or aud_removed:
        failures.append(f"unexpected audit added {sorted(unexpected_added)} "
                        f"removed {sorted(aud_removed)}")

    frozen = [k for k in added + removed + changed if bucket(k) == "C-frozen inputs"]
    print(f"\nC. FROZEN INPUTS differences: {len(frozen)} (must be 0)  {frozen}")
    if frozen:
        failures.append(f"frozen input differences {frozen}")

    impl = {k for k in added + removed + changed if bucket(k) == "D-implementation"}
    print(f"\nD. IMPLEMENTATION differences: {sorted(impl)}")
    print(f"   expected (the one file this arm edits): "
          f"{sorted(EXPECTED_IMPLEMENTATION_CHANGED)}")
    if impl != EXPECTED_IMPLEMENTATION_CHANGED:
        failures.append(f"implementation differences {sorted(impl)} != pre-declared")

    other = [k for k in added + removed + changed if bucket(k) == "F-other"]
    print(f"\nF. UNCLASSIFIED differences: {other} (must be empty)")
    if other:
        failures.append(f"unclassified differences {other}")

    print("\nPINNED scalars -- each must be byte-identical:")
    for key in PINNED:
        a, b = arm.get(key, "<missing>"), inc.get(key, "<missing>")
        ok = a == b and a != "<missing>"
        print(f"  {'OK  ' if ok else 'FAIL'} {key:64s} {b!r}"[:150])
        if not ok:
            failures.append(f"pinned {key}: {b!r} -> {a!r}")

    print(f"\nE. build output / runtime (reported, not asserted): "
          f"{sorted(k for k in changed if bucket(k) == 'E-output/runtime')[:8]} ...")

    print()
    if failures:
        print(f"FALSIFIER 1 **FIRES** -- {len(failures)} problems:")
        for f in failures:
            print(f"  - {f}")
        raise SystemExit(1)
    print("FALSIFIER 1 PASS: the config diff is EXACTLY the pre-declared coverage-aware "
          "sizing change plus its deterministic consequences.")


if __name__ == "__main__":
    main()
