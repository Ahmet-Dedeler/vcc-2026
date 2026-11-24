"""249 falsifier 1 -- the rule-10 audit diff, `cd4svd5` against the incumbent `h1union`.

222's `metadiff222.py` with the two pre-declared sets rewritten for THIS arm and **declared
before the diff was read** (PREREGISTRATION.md section 4, falsifier 4).  The arm sets ONE
environment variable, `VCC_CD4_DEFLATE_K=5`; the builder writes `config.cd4_deflate_k` and the
`effects.cd4.deflation` block only when it is on.  The deflated residual has a lower median row
norm, so the SAME sizing rule (`f` = 0.34, untouched) re-solves a higher fraction: the CD4
frame's `residual_median_row_norm_at_unit_fraction`, `f_per_unit_fraction`, `solved_fraction`
(and `effects.cd4.fraction`, `post_add_preserve.fraction`, `effective_f_at_cd4_fraction`, which
are the same number or its multiple) move, and so does the ONE merged preserve scalar (which
218's construction makes a property of both residuals; it is written twice, under
`effects.cd4.post_add_preserve` and `effects.h1_union.merged_preserve`).  Nothing else may.

Every difference is printed, not only the expected ones (mandate).

    uv run python metadiff249.py <arm.meta.json> <h1union.meta.json>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"

# ---------------------------------------------------------------------------------------
# PRE-DECLARED BEFORE THE DIFF WAS READ (rule 10).
EXPECTED_CONFIG_ADDED = {"cd4_deflate_k"}
EXPECTED_CONFIG_CHANGED: set[str] = set()
EXPECTED_AUDIT_CHANGED = {
    "effects.final_sha256",
    # the CD4 frame: the same rule re-solved on the deflated object
    "effects.cd4.fraction",
    "effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction",
    "effects.cd4.residual_frame.f_per_unit_fraction",
    "effects.cd4.residual_frame.solved_fraction",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction",
    "effects.cd4.post_add_preserve.fraction",
    # the ONE merged preserve scalar and the norm it is solved from, written in both blocks
    "effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm",
    "effects.cd4.post_add_preserve.post_add_global_preserve_scale",
    "effects.h1_union.merged_preserve.candidate_pre_preserve_median_norm",
    "effects.h1_union.merged_preserve.post_add_global_preserve_scale",
}
EXPECTED_AUDIT_ADDED_PREFIX = "effects.cd4.deflation."
EXPECTED_DEFLATION_KEYS = {"k", "n_carried", "energy_removed_frac",
                           "top_k_energy_frac_of_centred", "basis_u_k_sha256"}
# The direction each moved CD4 scalar must take (a deflated residual has LESS norm, so the
# solve for the same target `f` returns a LARGER fraction).
EXPECTED_DIRECTION = {
    "effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction": "down",
    "effects.cd4.residual_frame.f_per_unit_fraction": "down",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction": "down",
    "effects.cd4.residual_frame.solved_fraction": "up",
    "effects.cd4.fraction": "up",
    "effects.cd4.post_add_preserve.fraction": "up",
}
PINNED = [
    "config.alpha",
    "config.cd4_target_f",
    "config.cd4_fraction",
    "config.cd4_mode",
    "config.cd4_conditions",
    "config.cd4_target_policy",
    "config.h1_union_g",
    "config.h1_union_sources",
    "config.h1_union_sizing_rule",
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
    # CD4's support, form and requested dose: the single-variable constraint.
    "effects.cd4.mode",
    "effects.cd4.fully_supported_targets",
    "effects.cd4.incomplete_targets_preserved_before_global_scale",
    "effects.cd4.supported_entries",
    "effects.cd4.panel_residual_nonzero_entries",
    "effects.cd4.target_condition_count_histogram.0",
    "effects.cd4.target_condition_count_histogram.1",
    "effects.cd4.target_condition_count_histogram.2",
    "effects.cd4.target_condition_count_histogram.3",
    "effects.cd4.residual_frame.target_f_requested",
    "effects.cd4.residual_frame.whitened_base_median_row_norm",
    "effects.cd4.residual_frame.cd4_fraction_constant",
    "effects.cd4.residual_frame.fraction_source",
    "effects.cd4.post_add_preserve.baseline_output_median_norm",
    "frozen_bundle.sha256",
    "frozen_bundle.panel_id",
    "output.n_rows",
]
VOLATILE_PREFIX = ("output.", "peak_rss_gib", "runtime.")
EXPECTED_IMPLEMENTATION_CHANGED = {"hashes.implementation.scripts/make_submission_pm_cd4.py"}


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
    # `--same-builder`: the reference was built by the SAME builder bytes (the unflagged
    # known-answer rebuild), so the implementation bucket must show ZERO differences and no
    # post-222 metadata field can appear as an addition.  `--tag <name>` names the outputs.
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    same_builder = "--same-builder" in sys.argv
    out_tag = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else "metadiff249"
    expected_impl = set() if same_builder else EXPECTED_IMPLEMENTATION_CHANGED
    arm = flatten(json.loads(Path(args[0]).read_text()))
    inc = flatten(json.loads(Path(args[1]).read_text()))
    keys = sorted(set(arm) | set(inc))
    # EVERY `effects.h1_union.*` field of the incumbent is pinned except the merged preserve pair
    # (the mandate: "every effects.h1_union.* field unchanged").
    pinned = list(PINNED) + sorted(k for k in inc if k.startswith("effects.h1_union.")
                                   and k not in EXPECTED_AUDIT_CHANGED)

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
        if key.startswith(("effects.", "operator.", "donor_coverage.", "main_source_support.")):
            return "B-effects/operator"
        if key.startswith(("hashes.inputs", "frozen_bundle.")):
            return "C-frozen inputs"
        if key.startswith("hashes.implementation"):
            return "D-implementation"
        return "F-other"

    lines = []

    def P(s=""):
        print(s)
        lines.append(s)

    P("=" * 100)
    P("EVERY DIFFERENCE, not only the expected ones")
    P("=" * 100)
    for label, items in (("ADDED", added), ("REMOVED", removed), ("CHANGED", changed)):
        for key in items:
            tag = bucket(key)
            if label == "CHANGED":
                P(f"[{tag:20s}] {label:7s} {key}")
                P(f"{'':32s}  h1union: {inc[key]!r}"[:200])
                P(f"{'':32s}  cd4svd5: {arm[key]!r}"[:200])
            else:
                value = arm.get(key, inc.get(key))
                P(f"[{tag:20s}] {label:7s} {key} = {value!r}"[:200])
    P()

    failures = []
    cfg_added = {k.split("[")[0][len("config."):] for k in added if k.startswith("config.")}
    cfg_changed = {k[len("config."):] for k in changed if k.startswith("config.")}
    cfg_removed = {k[len("config."):] for k in removed if k.startswith("config.")}
    P(f"A. CONFIG added   {sorted(cfg_added)}\n   expected       {sorted(EXPECTED_CONFIG_ADDED)}")
    P(f"   CONFIG changed {sorted(cfg_changed)}   removed {sorted(cfg_removed)}")
    if cfg_added != EXPECTED_CONFIG_ADDED:
        failures.append(f"config additions {sorted(cfg_added)} != pre-declared")
    if cfg_changed != EXPECTED_CONFIG_CHANGED or cfg_removed:
        failures.append(f"config changed {sorted(cfg_changed)} removed {sorted(cfg_removed)}")
    if arm.get("config.cd4_deflate_k") != 5:
        failures.append(f"config.cd4_deflate_k is {arm.get('config.cd4_deflate_k')!r}, not 5")

    aud_changed = {k for k in changed if bucket(k) == "B-effects/operator"}
    aud_added = {k for k in added if bucket(k) == "B-effects/operator"}
    aud_removed = {k for k in removed if bucket(k) == "B-effects/operator"}
    unexpected_added = {k for k in aud_added if not k.startswith(EXPECTED_AUDIT_ADDED_PREFIX)}
    defl_keys = {k[len(EXPECTED_AUDIT_ADDED_PREFIX):] for k in aud_added
                 if k.startswith(EXPECTED_AUDIT_ADDED_PREFIX)}
    P(f"\nB. AUDIT changed {sorted(aud_changed)}\n   expected      {sorted(EXPECTED_AUDIT_CHANGED)}")
    P(f"   AUDIT added: {len(aud_added)} keys, {len(aud_added) - len(unexpected_added)} of them under "
      f"`{EXPECTED_AUDIT_ADDED_PREFIX}` = {sorted(defl_keys)}; unexpected {sorted(unexpected_added)}")
    if aud_changed != EXPECTED_AUDIT_CHANGED:
        failures.append(f"audit changes {sorted(aud_changed)} != pre-declared "
                        f"(extra {sorted(aud_changed - EXPECTED_AUDIT_CHANGED)}, "
                        f"missing {sorted(EXPECTED_AUDIT_CHANGED - aud_changed)})")
    if unexpected_added or aud_removed:
        failures.append(f"unexpected audit added {sorted(unexpected_added)} removed {sorted(aud_removed)}")
    if defl_keys != EXPECTED_DEFLATION_KEYS:
        failures.append(f"deflation block keys {sorted(defl_keys)} != {sorted(EXPECTED_DEFLATION_KEYS)}")

    P("\n   the moved CD4 scalars, with direction:")
    for key, want in EXPECTED_DIRECTION.items():
        a, b = arm.get(key), inc.get(key)
        ok = (a > b) if want == "up" else (a < b)
        P(f"   {'OK  ' if ok else 'FAIL'} {key:70s} {b!r} -> {a!r}  (x{a / b:.6f}, expected {want})")
        if not ok:
            failures.append(f"{key} moved the wrong way: {b!r} -> {a!r}")
    same_scalar = (arm["effects.cd4.post_add_preserve.post_add_global_preserve_scale"]
                   == arm["effects.h1_union.merged_preserve.post_add_global_preserve_scale"])
    P(f"   {'OK  ' if same_scalar else 'FAIL'} the merged preserve scalar is ONE number written twice: "
      f"{inc['effects.cd4.post_add_preserve.post_add_global_preserve_scale']!r} -> "
      f"{arm['effects.cd4.post_add_preserve.post_add_global_preserve_scale']!r}")
    if not same_scalar:
        failures.append("the two preserve scalars differ")
    P(f"   the deflation block: " + ", ".join(
        f"{k}={arm[EXPECTED_AUDIT_ADDED_PREFIX + k]!r}" for k in sorted(defl_keys)))

    frozen = [k for k in added + removed + changed if bucket(k) == "C-frozen inputs"]
    P(f"\nC. FROZEN INPUTS differences: {len(frozen)} (must be 0)  {frozen}")
    if frozen:
        failures.append(f"frozen input differences {frozen}")
    impl = {k for k in added + removed + changed if bucket(k) == "D-implementation"}
    P(f"\nD. IMPLEMENTATION differences: {sorted(impl)}")
    P(f"   expected ({'same builder bytes: none' if same_builder else 'the one file this arm edits'}): {sorted(expected_impl)}")
    if impl != expected_impl:
        failures.append(f"implementation differences {sorted(impl)} != pre-declared {sorted(expected_impl)}")
    other = [k for k in added + removed + changed if bucket(k) == "F-other"]
    P(f"\nF. UNCLASSIFIED differences: {other} (must be empty)")
    if other:
        failures.append(f"unclassified differences {other}")

    P(f"\nPINNED scalars ({len(pinned)}) -- each must be byte-identical:")
    n_ok = 0
    for key in pinned:
        a, b = arm.get(key, "<missing>"), inc.get(key, "<missing>")
        ok = a == b and a != "<missing>"
        n_ok += ok
        if not ok or not key.startswith("effects.h1_union."):
            P(f"  {'OK  ' if ok else 'FAIL'} {key:64s} {b!r}"[:150])
        if not ok:
            failures.append(f"pinned {key}: {b!r} -> {a!r}")
    n_h1 = sum(k.startswith("effects.h1_union.") for k in pinned)
    P(f"  ... and {n_h1} `effects.h1_union.*` fields (every one of the incumbent's except the "
      f"merged preserve pair): {sum(arm.get(k) == inc.get(k) for k in pinned if k.startswith('effects.h1_union.'))}/{n_h1} identical")
    P(f"  pinned identical: {n_ok}/{len(pinned)}")

    P(f"\nE. build output / runtime (reported, not asserted): "
      f"{sorted(k for k in changed if bucket(k) == 'E-output/runtime')}")
    P()
    if failures:
        P(f"FALSIFIER 1 **FIRES** -- {len(failures)} problems:")
        for f in failures:
            P(f"  - {f}")
    else:
        P("FALSIFIER 1 PASS: the config diff is EXACTLY `cd4_deflate_k` plus the pre-declared CD4-frame "
          "consequences and the one merged preserve scalar; operator, alpha, cd4_target_f, h1_union_g and "
          "every effects.h1_union.* field are byte-identical.")
    (WORK / f"{out_tag}.log").write_text("\n".join(lines) + "\n")
    (WORK / f"{out_tag}.json").write_text(json.dumps(dict(
        added=added, removed=removed, changed=changed, failures=failures,
        moved_cd4={k: [inc.get(k), arm.get(k)] for k in EXPECTED_DIRECTION},
        preserve=[inc["effects.cd4.post_add_preserve.post_add_global_preserve_scale"],
                  arm["effects.cd4.post_add_preserve.post_add_global_preserve_scale"]],
        deflation={k: arm[EXPECTED_AUDIT_ADDED_PREFIX + k] for k in defl_keys},
        pinned_identical=n_ok, pinned_total=len(pinned)), indent=1, default=str))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
