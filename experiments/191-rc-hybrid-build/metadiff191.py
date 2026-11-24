"""191 falsifier 5 — the rule-10 audit diff, `rc_hybrid` against `a15_local`.

WHY `a15_local` AND NOT a15's RECORDED META.  170 §1 found a15 and a19 were built on
Linux/x86_64; this box is arm64, and `whitener.apply` is a pair of 17,942 x 17,942 GEMMs whose
BLAS reduction order is architecture-dependent.  So a15's recorded `effects.final_sha256`
cannot be byte-reproduced here and comparing against it would fail for a reason that has
nothing to do with the arm.  The correct local identity check is against a package built on
this machine from the same configuration, one variable apart.  a15's recorded audit is
compared separately, as a DRIFT row, by `--a15 <a15.meta.json>`.

WHAT IT PROVES.  `--emitter` is consumed only by `run()`'s validation, by the meta, and by
`write_submission_from_effects`, which is called AFTER `del donors, pm, cd4_banks, cd4,
operator, model`.  So the switch cannot reach the effect matrix, and the built package must
reproduce `a15_local`'s effects audit EXACTLY -- every sha, every whitening norm, the CD4
fraction, the operator sha.  Anything else means the patch reached the effect matrix.

    uv run python metadiff191.py <rc_hybrid.meta.json> <a15_local.meta.json> [a15.meta.json]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

# Pre-declared before the diff is read: the emitter and its four recorded settings, nothing
# else.  `resample_keep_top` is NOT here -- `rc_hybrid` does not take it, and a build that
# wrote it would mean two emission channels were live at once.
EXPECTED_CONFIG_DIFF = {
    "emitter", "rc_key", "rc_sign_policy", "rc_scorer_min_cpm", "rc_jitter_seed_base",
}
# The scalar a15's own alpha falsifier reads (AGENTS.md §11.1 rule 10).
DRIFT_SCALARS = [
    "effects.whitening.after_all_median_norm",
    "effects.whitening.before_all_median_norm",
    "effects.whitening.after_operator_median_norm",
    "effects.whitening.before_operator_median_norm",
    "effects.whitening.operator_preserve_scale",
    "effects.cd4.fraction",
    "effects.cd4.residual_frame.effective_f_at_cd4_fraction",
    "effects.cd4.residual_frame.f_per_unit_fraction",
    "effects.cd4.residual_frame.solved_fraction",
    "effects.cd4.residual_frame.whitened_base_median_row_norm",
    "effects.cd4.residual_frame.residual_median_row_norm_at_unit_fraction",
    "effects.cd4.post_add_preserve.baseline_output_median_norm",
    "effects.cd4.post_add_preserve.candidate_pre_preserve_median_norm",
    "effects.cd4.post_add_preserve.post_add_global_preserve_scale",
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

    print("=== A: effects + operator audit, rc_hybrid vs a15_local (rule 10) ===")
    fields = 0
    bad = []
    for section in ("effects", "operator"):
        left, right = flatten(built[section]), flatten(inc[section])
        fields += len(set(left) | set(right))
        for key, x, y in diff(left, right):
            bad.append((f"{section}.{key}", x, y))
    fb, fi = flatten(built["effects"]), flatten(inc["effects"])
    for key in sorted(fb):
        if key.endswith("sha256") or "median_norm" in key or key in (
            "whitening.operator_preserve_scale", "cd4.fraction"
        ):
            value, other = fb[key], fi.get(key, "<absent>")
            mark = "EXACT" if value == other else "DIFFERS"
            print(f"  effects.{key:<48}{str(value)[:34]:<36} vs {str(other)[:34]:<36} {mark}")
    for key in ("sha256", "columns_sha256"):
        value, other = built["operator"][key], inc["operator"][key]
        print(f"  operator.{key:<47}{value[:16]:<36} vs {other[:16]:<36} "
              f"{'EXACT' if value == other else 'DIFFERS'}")
    print(f"  fields compared: {fields}   differences: {len(bad)}")
    for key, x, y in bad:
        print(f"    !! {key}: {x!r} vs {y!r}")
    gate_a = not bad
    print(f"  A {'PASS (rule 10 clean)' if gate_a else 'FAIL'}")

    print("\n=== B: config diff (rule 10's intended-difference count) ===")
    cdiff = diff(flatten(built["config"]), flatten(inc["config"]))
    for key, x, y in cdiff:
        print(f"  {key:<26} rc_hybrid {x!r:<22} a15_local {y!r}")
    moved = {k for k, _, _ in cdiff}
    gate_b = moved == EXPECTED_CONFIG_DIFF
    print(f"  config keys that moved: {sorted(moved)}")
    print(f"  pre-declared expected : {sorted(EXPECTED_CONFIG_DIFF)}")
    print(f"  B {'PASS' if gate_b else 'FAIL'} ({len(cdiff)} differences)")

    print("\n=== C: every other meta difference, stated rather than omitted ===")
    skip = ("config.", "effects.", "operator.")
    for key, x, y in diff(flatten(built), flatten(inc)):
        if key.startswith(skip):
            continue
        print(f"  {key:<56}{str(x)[:32]:<34} vs {str(y)[:32]}")

    gate_d = True
    if recorded is not None:
        print("\n=== D: a15_local against a15's RECORDED audit (the arm64/x86_64 drift) ===")
        fr = flatten(recorded)
        fl = flatten(inc)
        print(f"  {'scalar':<62}{'a15_local':>24}{'a15 recorded':>24}{'rel drift':>12}")
        worst = 0.0
        for key in DRIFT_SCALARS:
            left, right = fl.get(key), fr.get(key)
            if not isinstance(left, (int, float)) or not isinstance(right, (int, float)):
                continue
            rel = abs(left - right) / max(abs(right), 1e-300)
            worst = max(worst, rel)
            print(f"  {key:<62}{left:>24.15g}{right:>24.15g}{rel:>12.1e}")
        for key in ("effects.final_sha256", "effects.post_whiten_sha256",
                    "effects.baseline_pre_sha256", "effects.candidate_pre_sha256",
                    "operator.sha256", "operator.columns_sha256"):
            left, right = str(fl.get(key))[:16], str(fr.get(key))[:16]
            print(f"  {key:<62}{left:>24}{right:>24}"
                  f"{'  EXACT' if left == right else '  DIFFERS':>12}")
        print(f"  worst relative drift on the 14 audit scalars: {worst:.1e}")
        print(f"  (170 §1: a15/a19 were built on Linux/x86_64; this box is arm64 and "
              f"`whitener.apply` is a pair of 17,942^2 GEMMs)")
        gate_d = worst < 1e-13
        print(f"  D {'PASS' if gate_d else 'FAIL'} (drift is last-bit noise)")

    print(f"\nFALSIFIER 5 {'PASS' if (gate_a and gate_b) else 'FAIL'}"
          f"   | drift row {'PASS' if gate_d else 'FAIL'}")


if __name__ == "__main__":
    main()
