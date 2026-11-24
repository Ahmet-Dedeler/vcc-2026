"""217 falsifier 1/2 — the rule-10 audit diff, `cd4f080` against the incumbent `cd4f034`.

`experiments/210-cd4-build/metadiff210.py` with the two pre-declared sets rewritten for
this arm and declared BEFORE the diff was read.  202's arm changed the EMITTER, so its section
A (the effects + operator audit) had to be EXACTLY zero differences.  This arm changes the CD4
DOSE, which is a profile change, so section A must differ in EXACTLY the deterministic
consequences of the dose and in nothing else -- above all `operator.sha256`,
`effects.whitening.after_all_median_norm`, `effects.post_whiten_sha256` (the whitened base,
computed BEFORE the CD4 add) and `effects.cd4.fully_supported_targets` must be untouched.

Every difference is printed, not only the expected ones (mandate).

    uv run python metadiff217.py <cd4f080.meta.json> <210-cd4f034.meta.json> [a15.meta.json]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# ---------------------------------------------------------------------------------------
# PRE-DECLARED BEFORE THE DIFF WAS READ (rule 10).
#
# CONFIG: the arm sets one environment variable, so exactly one config key may move.
EXPECTED_CONFIG_DIFF = {"cd4_target_f"}
#
# EFFECTS/OPERATOR: `VCC_CD4_TARGET_F` is consumed in exactly one place --
# `solve_cd4_fraction`, which divides it by the measured `f_per_unit_fraction` and hands the
# quotient to `add_post_whiten_residual`.  So the only audit fields that may move are that
# quotient (recorded three times under three names), the requested `f` echoed back, the two
# post-add preserve scalars it drives, and the sha of the FINAL effect matrix -- which must
# move, because a `final_sha256` that did not would mean the residual never applied.
EXPECTED_AUDIT_DIFF = {
    "effects.final_sha256",
    "effects.cd4.fraction",
    "effects.cd4.post_add_preserve.fraction",
    "effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm",
    "effects.cd4.post_add_preserve.post_add_global_preserve_scale",
    "effects.cd4.residual_frame.solved_fraction",
    "effects.cd4.residual_frame.target_f_requested",
}
# The scalars that must NOT move, named individually so a silent drift cannot hide in the
# aggregate count (AGENTS.md rule 10: the cheap form of a15's alpha falsifier).
PINNED = [
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
    "effects.cd4.fully_supported_targets",
    "effects.cd4.supported_entries",
    "effects.cd4.mode",
    "effects.cd4.panel_residual_nonzero_entries",
    "effects.cd4.residual_frame.cd4_fraction_constant",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction",
    "effects.cd4.residual_frame.f_per_unit_fraction",
    "effects.cd4.residual_frame.fraction_source",
    "effects.cd4.residual_frame.whitened_base_median_row_norm",
    "effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction",
    "effects.cd4.post_add_preserve.baseline_output_median_norm",
    "operator.sha256",
    "operator.columns_sha256",
    "operator.metadata.sha256",
    "operator.metadata_sha256",
]
# 207 section 8.2's required values, asserted against the BUILT meta rather than the command.
REQUIRED = {
    "config.cd4_target_f": 0.80,
    "config.alpha": 0.982,
    "config.emitter": "rc_callset",
    "config.whitener_id": "union_pruned",
    "effects.cd4.fully_supported_targets": 276,
    "effects.cd4.residual_frame.target_f_requested": 0.80,
    "effects.cd4.residual_frame.solved_fraction": 0.80 / 1.0556645369046553,
    "effects.cd4.residual_frame.fraction_source": "VCC_CD4_TARGET_F",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction": 0.05278322684523277,
    "effects.cd4.residual_frame.f_per_unit_fraction": 1.0556645369046553,
    "effects.whitening.after_all_median_norm": 4.229295604454299,
    "operator.sha256":
        "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
}
DRIFT_SCALARS = [
    "effects.whitening.after_all_median_norm",
    "effects.whitening.before_all_median_norm",
    "effects.whitening.after_operator_median_norm",
    "effects.whitening.before_operator_median_norm",
    "effects.whitening.operator_preserve_scale",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction",
    "effects.cd4.residual_frame.f_per_unit_fraction",
    "effects.cd4.residual_frame.whitened_base_median_row_norm",
    "effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction",
    "effects.cd4.post_add_preserve.baseline_output_median_norm",
]


def flatten(obj, prefix=""):
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.update(flatten(v, f"{prefix}{k}."))
    elif isinstance(obj, list):
        out[prefix[:-1]] = json.dumps(obj)
    else:
        out[prefix[:-1]] = obj
    return out


def diff(a: dict, b: dict) -> list[tuple[str, object, object]]:
    keys = sorted(set(a) | set(b))
    return [(k, a.get(k, "<absent>"), b.get(k, "<absent>"))
            for k in keys if a.get(k, "<absent>") != b.get(k, "<absent>")]


def main() -> None:
    built = json.loads(Path(sys.argv[1]).read_text())
    inc = json.loads(Path(sys.argv[2]).read_text())
    recorded = json.loads(Path(sys.argv[3]).read_text()) if len(sys.argv) > 3 else None
    fb, fi = flatten(built), flatten(inc)

    print("=== A: effects + operator audit, cd4f080 vs cd4f034 (rule 10) ===")
    bad = []
    fields = 0
    for section in ("effects", "operator"):
        left, right = flatten(built[section]), flatten(inc[section])
        fields += len(set(left) | set(right))
        for key, x, y in diff(left, right):
            bad.append((f"{section}.{key}", x, y))
    moved_audit = {k for k, _, _ in bad}
    print(f"  fields compared: {fields}   differences: {len(bad)}")
    for key, x, y in bad:
        flag = "expected" if key in EXPECTED_AUDIT_DIFF else "!! UNEXPECTED"
        print(f"    {flag:<14}{key:<58}{str(x)[:26]:<28} vs {str(y)[:26]}")
    missing = EXPECTED_AUDIT_DIFF - moved_audit
    for key in sorted(missing):
        print(f"    !! DID NOT MOVE {key}")
    gate_a = moved_audit == EXPECTED_AUDIT_DIFF
    print(f"  moved   : {sorted(moved_audit)}")
    print(f"  expected: {sorted(EXPECTED_AUDIT_DIFF)}")
    print(f"  A {'PASS' if gate_a else 'FAIL'}")

    print("\n  PINNED audit fields (must be EXACT):")
    pin_ok = True
    for key in PINNED:
        x, y = fb.get(key, "<absent>"), fi.get(key, "<absent>")
        mark = "EXACT" if x == y else "DIFFERS"
        pin_ok &= x == y
        print(f"    {key:<62}{str(x)[:22]:<24} vs {str(y)[:22]:<24} {mark}")
    print(f"  PINNED {'PASS' if pin_ok else 'FAIL'}")

    print("\n=== B: config diff (rule 10's intended-difference count) ===")
    cdiff = diff(flatten(built["config"]), flatten(inc["config"]))
    for key, x, y in cdiff:
        print(f"  {key:<26} cd4f080 {x!r:<22} cd4f034 {y!r}")
    moved = {k for k, _, _ in cdiff}
    gate_b = moved == EXPECTED_CONFIG_DIFF
    print(f"  config keys that moved: {sorted(moved)}")
    print(f"  pre-declared expected : {sorted(EXPECTED_CONFIG_DIFF)}")
    print(f"  B {'PASS' if gate_b else 'FAIL'} ({len(cdiff)} differences)")

    print("\n=== B2: 207 section 8.2's required meta values, against the BUILT meta ===")
    req_ok = True
    for key, want in REQUIRED.items():
        got = fb.get(key, "<absent>")
        ok = (abs(got - want) <= 1e-12 * max(1.0, abs(want))
              if isinstance(want, float) and isinstance(got, (int, float))
              else got == want)
        req_ok &= bool(ok)
        print(f"  {key:<62}{str(got)[:34]:<36}want {str(want)[:24]:<26}"
              f"{'OK' if ok else '!! WRONG'}")
    # the dose's algebra, checked rather than trusted
    f_req = fb.get("effects.cd4.residual_frame.target_f_requested")
    per_unit = fb.get("effects.cd4.residual_frame.f_per_unit_fraction")
    solved = fb.get("effects.cd4.residual_frame.solved_fraction")
    applied = fb.get("effects.cd4.fraction")
    preserved = fb.get("effects.cd4.post_add_preserve.fraction")
    algebra = (abs(solved - f_req / per_unit) < 1e-12 and solved == applied == preserved)
    req_ok &= algebra
    print(f"  solved_fraction == target_f / f_per_unit: {solved!r} == "
          f"{f_req}/{per_unit} = {f_req / per_unit!r}   "
          f"and cd4.fraction == post_add_preserve.fraction == solved   "
          f"{'OK' if algebra else '!! WRONG'}")
    inc_solved = fi.get("effects.cd4.residual_frame.solved_fraction")
    print(f"  dose ratio delivered: solved {solved:.9f} / incumbent {inc_solved:.9f} = "
          f"{solved / inc_solved:.6f}x   (requested 0.80 / 0.34 = {0.80 / 0.34:.6f}x)")
    print(f"  B2 {'PASS' if req_ok else 'FAIL'}")

    print("\n=== C: every other meta difference, stated rather than omitted ===")
    skip = ("config.", "effects.", "operator.")
    for key, x, y in diff(fb, fi):
        if key.startswith(skip):
            continue
        print(f"  {key:<56}{str(x)[:32]:<34} vs {str(y)[:32]}")

    gate_d = True
    if recorded is not None:
        print("\n=== D: the arm against a15's RECORDED audit (the arm64/x86_64 drift) ===")
        fr = flatten(recorded)
        print(f"  {'scalar':<62}{'cd4f080':>24}{'a15 recorded':>24}{'rel drift':>12}")
        worst = 0.0
        for key in DRIFT_SCALARS:
            left, right = fb.get(key), fr.get(key)
            if not isinstance(left, (int, float)) or not isinstance(right, (int, float)):
                continue
            rel = abs(left - right) / max(abs(right), 1e-300)
            worst = max(worst, rel)
            print(f"  {key:<62}{left:>24.15g}{right:>24.15g}{rel:>12.1e}")
        for key in ("effects.baseline_pre_sha256", "effects.candidate_pre_sha256",
                    "effects.post_whiten_sha256", "operator.sha256",
                    "operator.columns_sha256"):
            left, right = str(fb.get(key))[:16], str(fr.get(key))[:16]
            print(f"  {key:<62}{left:>24}{right:>24}"
                  f"{'  EXACT' if left == right else '  DIFFERS':>12}")
        print(f"  worst relative drift on the pinned audit scalars: {worst:.1e}")
        gate_d = worst < 1e-13
        print(f"  D {'PASS' if gate_d else 'FAIL'} (drift is last-bit noise)")

    ok = gate_a and gate_b and pin_ok and req_ok
    print(f"\nPREREG FALSIFIER 5 (rule-10 config diff) {'PASS' if ok else 'FAIL'}"
          f"   | drift row {'PASS' if gate_d else 'FAIL'}")


if __name__ == "__main__":
    main()
