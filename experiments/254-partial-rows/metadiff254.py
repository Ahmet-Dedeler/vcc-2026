"""254 falsifier 1 -- the rule-10 audit diff, the composed arm against its incumbent.

249's `metadiff249.py` with the pre-declared sets rewritten for THIS arm and **declared before
the diff was read**.  Two references:

  --against cd4svd5   (the mandate's incumbent) `submissions/249-cd4svd5-f034.h5ad.meta.json`, built
                      by 249's builder bytes `2d63b979...`: config ADDED exactly `cd4_partial_g`; audit
                      CHANGED at most `effects.final_sha256` (must) and the ONE merged preserve scalar
                      and the norm it is solved from, written twice (may -- the median over 300 rows
                      moves only if one of the 15 partial rows crosses it); audit ADDED only under
                      `effects.cd4_partial.`; `config.cd4_deflate_k` 5, `effects.cd4.fraction`, the
                      CD4 frame, the deflation block and every `effects.h1_union.*` field except the
                      preserve pair PINNED; implementation differs in the one file this arm edits.
  --against ka        the same-builder unflagged rebuild `submissions/254-ka-h1union-g028.h5ad.meta.json`
                      (= the h1union configuration through THIS builder): config ADDED exactly
                      {`cd4_deflate_k`, `cd4_partial_g`}; audit CHANGED = 249's 11 pre-declared keys;
                      ADDED under `effects.cd4.deflation.` (249's 5) and `effects.cd4_partial.`;
                      implementation differences NONE.

Every difference is printed, not only the expected ones (mandate).

    uv run python metadiff254.py <arm.meta.json> <reference.meta.json> --against cd4svd5|ka [--tag name]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
S = WORK / "stage"

# ---------------------------------------------------------------------------------------
# PRE-DECLARED BEFORE THE DIFF WAS READ (rule 10).
PRESERVE_PAIR = {
    "effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm",
    "effects.cd4.post_add_preserve.post_add_global_preserve_scale",
    "effects.h1_union.merged_preserve.candidate_pre_preserve_median_norm",
    "effects.h1_union.merged_preserve.post_add_global_preserve_scale",
}
CD4_FRAME_249 = {
    "effects.cd4.fraction",
    "effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction",
    "effects.cd4.residual_frame.f_per_unit_fraction",
    "effects.cd4.residual_frame.solved_fraction",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction",
    "effects.cd4.post_add_preserve.fraction",
}
DEFLATION_KEYS = {"k", "n_carried", "energy_removed_frac", "top_k_energy_frac_of_centred", "basis_u_k_sha256"}
PARTIAL_PREFIX = "effects.cd4_partial."
PARTIAL_EXPECTED = {
    "rows": 15, "rows_one_condition": 4, "rows_two_conditions": 11, "carrying_rows": 15, "target_g": 0.34,
    "sizing_rule": "median over CARRYING rows (coverage-aware)", "mode": "cd4_partial", "merged_add": True,
    "panel_residual_nonzero_entries": 0, "residual_frame.shipped_median_rule_would_solve": 0.0,
    "residual_frame.carrying_rows": 15, "residual_frame.coverage_frac": 0.05,
    "residual_frame.target_g_requested": 0.34, "residual_frame.whitened_base_median_row_norm": 4.229295604454299,
    "target_condition_count_histogram.0": 285, "target_condition_count_histogram.1": 4, "target_condition_count_histogram.2": 11,
}
PROFILES = {
    "cd4svd5": dict(
        config_added={"cd4_partial_g"},
        audit_changed_must={"effects.final_sha256"},
        audit_changed_may=set(PRESERVE_PAIR),
        audit_added_prefixes=(PARTIAL_PREFIX,),
        impl={"hashes.implementation.scripts/make_submission_pm_cd4.py"},
        pinned_extra=["config.cd4_deflate_k"] + sorted(CD4_FRAME_249)
        + [f"effects.cd4.deflation.{k}" for k in sorted(DEFLATION_KEYS)]
        + ["effects.cd4.post_add_preserve.baseline_output_median_norm"],
    ),
    "ka": dict(
        config_added={"cd4_deflate_k", "cd4_partial_g"},
        audit_changed_must={"effects.final_sha256"} | CD4_FRAME_249 | PRESERVE_PAIR,
        audit_changed_may=set(),
        audit_added_prefixes=("effects.cd4.deflation.", PARTIAL_PREFIX),
        impl=set(),
        pinned_extra=["effects.cd4.post_add_preserve.baseline_output_median_norm"],
    ),
}
PINNED_BASE = [
    "config.alpha", "config.cd4_target_f", "config.cd4_fraction", "config.cd4_mode", "config.cd4_conditions",
    "config.cd4_target_policy", "config.h1_union_g", "config.h1_union_sources", "config.h1_union_sizing_rule",
    "config.whitener_id", "config.whiten_gamma", "config.whiten_lambda", "config.emitter", "config.min_cpm",
    "config.keep_top", "config.keep_spread", "config.tail_spread", "config.overdispersion", "config.spread_map",
    "config.effect_scale", "config.emission_seed", "config.rc_key", "config.rc_callset_head", "config.rc_callset_cpm",
    "config.rc_callset_z_cap",
    "operator.sha256", "operator.metadata_sha256", "operator.columns_sha256", "operator.columns",
    "operator.all_panel_coordinates_excluded", "operator.axis_sha256",
    "effects.baseline_pre_sha256", "effects.candidate_pre_sha256", "effects.post_whiten_sha256",
    "effects.whitening.after_all_median_norm", "effects.whitening.before_all_median_norm",
    "effects.whitening.after_operator_median_norm", "effects.whitening.before_operator_median_norm",
    "effects.whitening.operator_preserve_scale", "effects.whitening.gamma", "effects.whitening.lambda",
    "effects.cd4.mode", "effects.cd4.fully_supported_targets",
    "effects.cd4.incomplete_targets_preserved_before_global_scale", "effects.cd4.supported_entries",
    "effects.cd4.panel_residual_nonzero_entries",
    "effects.cd4.target_condition_count_histogram.0", "effects.cd4.target_condition_count_histogram.1",
    "effects.cd4.target_condition_count_histogram.2", "effects.cd4.target_condition_count_histogram.3",
    "effects.cd4.residual_frame.target_f_requested", "effects.cd4.residual_frame.whitened_base_median_row_norm",
    "effects.cd4.residual_frame.cd4_fraction_constant", "effects.cd4.residual_frame.fraction_source",
    "frozen_bundle.sha256", "frozen_bundle.panel_id", "output.n_rows",
]
VOLATILE_PREFIX = ("output.", "peak_rss_gib", "runtime.")


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
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    against = sys.argv[sys.argv.index("--against") + 1]
    prof = PROFILES[against]
    out_tag = sys.argv[sys.argv.index("--tag") + 1] if "--tag" in sys.argv else f"metadiff254-{against}"
    arm_raw = json.loads(Path(args[0]).read_text())
    arm = flatten(arm_raw)
    inc = flatten(json.loads(Path(args[1]).read_text()))
    keys = sorted(set(arm) | set(inc))
    pinned = list(PINNED_BASE) + list(prof["pinned_extra"]) + sorted(
        k for k in inc if k.startswith("effects.h1_union.") and k not in PRESERVE_PAIR)
    pinned = [k for k in pinned if k in inc]     # (249's meta has no operator.axis_sha256? it does -- keep only what the ref carries)
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
    P(f"EVERY DIFFERENCE, not only the expected ones  (arm {args[0]}  vs  {args[1]}, profile {against})")
    P("=" * 100)
    for label, items in (("ADDED", added), ("REMOVED", removed), ("CHANGED", changed)):
        for key in items:
            tag = bucket(key)
            if label == "CHANGED":
                P(f"[{tag:20s}] {label:7s} {key}")
                P(f"{'':32s}  reference: {inc[key]!r}"[:200])
                P(f"{'':32s}  arm:       {arm[key]!r}"[:200])
            else:
                value = arm.get(key, inc.get(key))
                P(f"[{tag:20s}] {label:7s} {key} = {value!r}"[:200])
    P()
    failures = []
    cfg_added = {k.split("[")[0][len("config."):] for k in added if k.startswith("config.")}
    cfg_changed = {k[len("config."):] for k in changed if k.startswith("config.")}
    cfg_removed = {k[len("config."):] for k in removed if k.startswith("config.")}
    P(f"A. CONFIG added   {sorted(cfg_added)}\n   expected       {sorted(prof['config_added'])}")
    P(f"   CONFIG changed {sorted(cfg_changed)}   removed {sorted(cfg_removed)}")
    if cfg_added != prof["config_added"]:
        failures.append(f"config additions {sorted(cfg_added)} != pre-declared {sorted(prof['config_added'])}")
    if cfg_changed or cfg_removed:
        failures.append(f"config changed {sorted(cfg_changed)} removed {sorted(cfg_removed)}")
    if arm.get("config.cd4_partial_g") != 0.34:
        failures.append(f"config.cd4_partial_g is {arm.get('config.cd4_partial_g')!r}, not 0.34")
    if arm.get("config.cd4_deflate_k") != 5:
        failures.append(f"config.cd4_deflate_k is {arm.get('config.cd4_deflate_k')!r}, not 5")

    aud_changed = {k for k in changed if bucket(k) == "B-effects/operator"}
    aud_added = {k for k in added if bucket(k) == "B-effects/operator"}
    aud_removed = {k for k in removed if bucket(k) == "B-effects/operator"}
    allowed_changed = prof["audit_changed_must"] | prof["audit_changed_may"]
    unexpected_changed = aud_changed - allowed_changed
    missing_must = prof["audit_changed_must"] - aud_changed
    unexpected_added = {k for k in aud_added if not k.startswith(prof["audit_added_prefixes"])}
    partial_keys = {k[len(PARTIAL_PREFIX):] for k in aud_added if k.startswith(PARTIAL_PREFIX)}
    P(f"\nB. AUDIT changed {sorted(aud_changed)}\n   must {sorted(prof['audit_changed_must'])}\n   may  {sorted(prof['audit_changed_may'])}")
    P(f"   AUDIT added: {len(aud_added)} keys; under {prof['audit_added_prefixes']}: {len(aud_added) - len(unexpected_added)}; unexpected {sorted(unexpected_added)}")
    P(f"   `effects.cd4_partial.` keys: {sorted(partial_keys)}")
    if unexpected_changed or missing_must:
        failures.append(f"audit changes: unexpected {sorted(unexpected_changed)}, missing {sorted(missing_must)}")
    if unexpected_added or aud_removed:
        failures.append(f"unexpected audit added {sorted(unexpected_added)} removed {sorted(aud_removed)}")
    preserve_moved = bool(aud_changed & PRESERVE_PAIR)
    P(f"   the merged preserve pair {'MOVED' if preserve_moved else 'did NOT move'}: "
      f"{inc.get('effects.cd4.post_add_preserve.post_add_global_preserve_scale')!r} -> "
      f"{arm.get('effects.cd4.post_add_preserve.post_add_global_preserve_scale')!r}")
    same_scalar = (arm["effects.cd4.post_add_preserve.post_add_global_preserve_scale"]
                   == arm["effects.h1_union.merged_preserve.post_add_global_preserve_scale"]
                   == arm.get("effects.cd4_partial.merged_preserve.post_add_global_preserve_scale"))
    P(f"   {'OK  ' if same_scalar else 'FAIL'} the preserve scalar is ONE number written in all three blocks")
    if not same_scalar:
        failures.append("the three preserve scalars differ")
    P("\n   the partial block's pre-declared values:")
    for k, want in PARTIAL_EXPECTED.items():
        got = arm.get(PARTIAL_PREFIX + k, "<missing>")
        ok = got == want
        P(f"   {'OK  ' if ok else 'FAIL'} {PARTIAL_PREFIX + k:64s} {got!r} (expected {want!r})")
        if not ok:
            failures.append(f"{PARTIAL_PREFIX + k} = {got!r}, expected {want!r}")
    staged = json.loads((S / "stage-composed.json").read_text())
    sf = staged["audit"]["cd4_partial"]["residual_frame"]["solved_fraction"]
    got_sf = arm.get(PARTIAL_PREFIX + "residual_frame.solved_fraction")
    P(f"   {'OK  ' if got_sf == sf else 'FAIL'} the partial's solved fraction == the staged row's: {got_sf!r} vs {sf!r}")
    if got_sf != sf:
        failures.append("the partial's solved fraction differs from the staged row's")
    fs_ok = arm.get("effects.final_sha256") == staged["effects_array_sha256"]
    P(f"   {'OK  ' if fs_ok else 'FAIL'} effects.final_sha256 == the staged composed array {staged['effects_array_sha256'][:16]}...")
    if not fs_ok:
        failures.append("effects.final_sha256 != the staged composed array")

    frozen = [k for k in added + removed + changed if bucket(k) == "C-frozen inputs"]
    P(f"\nC. FROZEN INPUTS differences: {len(frozen)} (must be 0)  {frozen}")
    if frozen:
        failures.append(f"frozen input differences {frozen}")
    impl = {k for k in added + removed + changed if bucket(k) == "D-implementation"}
    P(f"\nD. IMPLEMENTATION differences: {sorted(impl)}\n   expected: {sorted(prof['impl'])}")
    if impl != prof["impl"]:
        failures.append(f"implementation differences {sorted(impl)} != pre-declared {sorted(prof['impl'])}")
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
    P(f"  ... and {n_h1} `effects.h1_union.*` fields (every one of the reference's except the preserve pair): "
      f"{sum(arm.get(k) == inc.get(k) for k in pinned if k.startswith('effects.h1_union.'))}/{n_h1} identical")
    P(f"  pinned identical: {n_ok}/{len(pinned)}")
    P(f"\nE. build output / runtime (reported, not asserted): {sorted(k for k in changed if bucket(k) == 'E-output/runtime')}")
    P()
    if failures:
        P(f"FALSIFIER 1 ({against}) **FIRES** -- {len(failures)} problems:")
        for f in failures:
            P(f"  - {f}")
    else:
        P(f"FALSIFIER 1 ({against}) PASS: the config diff is EXACTLY {sorted(prof['config_added'])} plus the pre-declared consequences.")
    (WORK / f"{out_tag}.log").write_text("\n".join(lines) + "\n")
    (WORK / f"{out_tag}.json").write_text(json.dumps(dict(
        against=against, added=added, removed=removed, changed=changed, failures=failures, preserve_moved=preserve_moved,
        preserve=[inc.get("effects.cd4.post_add_preserve.post_add_global_preserve_scale"),
                  arm.get("effects.cd4.post_add_preserve.post_add_global_preserve_scale")],
        partial={k: arm.get(PARTIAL_PREFIX + k) for k in partial_keys}, pinned_identical=n_ok, pinned_total=len(pinned)),
        indent=1, default=str))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
