"""179 gates 2 and 3 — the rule-10 audit diff between the built package and a23.

GATE 2 (rule 10, the bit-identity proof).
The `--emitter` / `--resample-keep-top` switches live entirely downstream of the effects
stage: `run()` consumes them only to validate, to write the meta and to hand them to
`write_submission_from_effects`, which is called AFTER `del donors, pm, cd4_banks, cd4,
operator, model`. So a build with the switches ON must reproduce the incumbent's effects
audit EXACTLY -- `baseline_pre_sha256`, `candidate_pre_sha256`, `post_whiten_sha256`,
`final_sha256`, every whitening norm, the CD4 fraction and the operator sha. 178 ran the same
check for its three riders and got 0 differences on 9 fields; anything else here means the
patch reached the effect matrix, which it must not.

GATE 3 (rule 10, the intended-difference count).
The finished artifact's `config` block must differ from a23's in EXACTLY two keys, `emitter`
and `resample_keep_top`, both of them added. Rule 10 exists because two near-misses on
2026-08-31 were invisible in the invocation and visible only here.

Differences OUTSIDE `config` are reported separately and explicitly, because some are
unavoidable and pretending otherwise would be the vacuous version of this check: the builder's
own source sha moved (it is the file that was patched), and the emitted artifact's size, nnz
and sha256 moved (that is the point of the arm).

    uv run python metadiff179.py <built.h5ad.meta.json> <a23.h5ad.meta.json>
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

EXPECTED_CONFIG_DIFF = {"emitter", "resample_keep_top"}


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

    print("=== GATE 2: effects + operator audit vs the incumbent (rule 10) ===")
    fields = 0
    bad = []
    for section in ("effects", "operator"):
        d = flatten(built[section]), flatten(inc[section])
        fields += len(set(d[0]) | set(d[1]))
        for key, x, y in diff(*d):
            bad.append((f"{section}.{key}", x, y))
    for key in sorted(flatten(built["effects"])):
        value = flatten(built["effects"])[key]
        other = flatten(inc["effects"]).get(key, "<absent>")
        mark = "EXACT" if value == other else "DIFFERS"
        if key.endswith("sha256") or "median_norm" in key or key in (
            "whitening.operator_preserve_scale", "cd4.fraction"
        ):
            print(f"  effects.{key:<45} {str(value)[:40]:<42} vs {str(other)[:40]:<42} {mark}")
    print(f"  operator.sha256{'':<38} {built['operator']['sha256'][:16]} vs "
          f"{inc['operator']['sha256'][:16]} "
          f"{'EXACT' if built['operator']['sha256'] == inc['operator']['sha256'] else 'DIFFERS'}")
    print(f"  operator.columns_sha256{'':<30} {built['operator']['columns_sha256'][:16]} vs "
          f"{inc['operator']['columns_sha256'][:16]} "
          f"{'EXACT' if built['operator']['columns_sha256'] == inc['operator']['columns_sha256'] else 'DIFFERS'}")
    print(f"  fields compared: {fields}   differences: {len(bad)}")
    for key, x, y in bad:
        print(f"    !! {key}: {x!r} vs {y!r}")
    gate2 = not bad
    print(f"  GATE 2 {'PASS (rule 10 clean)' if gate2 else 'FAIL'}")

    print("\n=== GATE 3: config diff (rule 10's intended-difference count) ===")
    cb, ci = built["config"], inc["config"]
    cdiff = diff(flatten(cb), flatten(ci))
    for key, x, y in cdiff:
        print(f"  {key:<24} built {x!r:<20} incumbent {y!r}")
    moved = {k for k, _, _ in cdiff}
    gate3 = moved == EXPECTED_CONFIG_DIFF
    print(f"  config keys that moved: {sorted(moved)}   expected {sorted(EXPECTED_CONFIG_DIFF)}")
    print(f"  GATE 3 {'PASS' if gate3 else 'FAIL'} ({len(cdiff)} differences)")

    print("\n=== outside `config`: every other meta difference, stated ===")
    skip = ("config.", "effects.", "operator.")
    other = [d for d in diff(flatten(built), flatten(inc))
             if not d[0].startswith(skip)]
    for key, x, y in other:
        sx, sy = str(x), str(y)
        print(f"  {key:<58} {sx[:34]:<36} vs {sy[:34]}")

    print(f"\nGATE 2 {'PASS' if gate2 else 'FAIL'} | GATE 3 {'PASS' if gate3 else 'FAIL'}")


if __name__ == "__main__":
    main()
