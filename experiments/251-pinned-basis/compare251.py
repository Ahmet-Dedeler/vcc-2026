"""251 known-answer row 2 -- the pinned staged effects against 249's re-fit staged effects.

Reads `work/stage/effects-pinned.npy` / `stage-pinned.json` (this experiment) and 249's
`work/stage/effects-svd5.npy` / `stage-svd5.json` (sha a28d5be2..., the armed `cd4svd5` arm) and
this experiment's own `k5` row, and reports: array sha256 of each, max |d| (absolute and relative
to the largest entry), the number of differing entries and which rows they sit on (CD4-carried
vs not), the solved fractions, the preserve scalars and the deflation audits side by side.

The pre-declared reading: the two are ONE projector in two floating-point associations
(`U(U^T R_c)` vs `(R_c V) V^T`), so byte identity is NOT expected; the row PASSES at numeric
identity if max |d| <= 1e-12 on the effect matrix, the fractions agree to 1e-14 and the energy
numbers to 1e-12; it FAILS on anything larger.  Both readings are printed.

    uv run python experiments/251-pinned-basis/compare251.py
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
S251 = HERE / "work" / "stage"
S249 = ROOT / "experiments/249-cd4-deflate/work/stage"


def sha(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def main() -> None:
    pin = np.load(S251 / "effects-pinned.npy")
    ref = np.load(S249 / "effects-svd5.npy")
    k5 = np.load(S251 / "effects-k5.npy") if (S251 / "effects-k5.npy").exists() else None
    off = np.load(S249 / "effects-off.npy")
    jp = json.loads((S251 / "stage-pinned.json").read_text())
    jr = json.loads((S249 / "stage-svd5.json").read_text())
    jk = json.loads((S251 / "stage-k5.json").read_text()) if (S251 / "stage-k5.json").exists() else None
    cd4 = np.load(S251 / "cd4-pinned.npz")
    carried = np.asarray(cd4["target_condition_count"]) > 0
    d = pin - ref
    amax = float(np.abs(d).max())
    rel = amax / float(np.abs(ref).max())
    ndiff = int((d != 0).sum())
    rows_diff = np.flatnonzero(np.abs(d).max(axis=1) > 0)
    rows_diff_carried = int(carried[rows_diff].sum())
    out = {
        "pinned_array_sha256": sha(pin), "pinned_json_sha256": jp["effects_array_sha256"],
        "refit_249_array_sha256": sha(ref), "refit_249_json_sha256": jr["effects_array_sha256"],
        "k5_this_builder_array_sha256": None if k5 is None else sha(k5),
        "k5_equals_249_bitwise": None if k5 is None else bool(np.array_equal(k5, ref)),
        "bitwise_equal": bool(np.array_equal(pin, ref)),
        "max_abs_diff": amax, "max_abs_entry": float(np.abs(ref).max()), "max_rel_to_largest_entry": rel,
        "n_entries_differing": ndiff, "n_entries": int(d.size),
        "n_rows_differing": int(rows_diff.size), "n_rows_differing_cd4_carried": rows_diff_carried,
        "n_rows_differing_not_carried": int(rows_diff.size - rows_diff_carried),
        "max_abs_diff_not_carried_rows": float(np.abs(d[~carried]).max()) if (~carried).any() else 0.0,
        "pinned_vs_off_max_abs_diff": float(np.abs(pin - off).max()),
        "fraction_pinned": jp["audit"]["cd4"]["fraction"], "fraction_refit": jr["audit"]["cd4"]["fraction"],
        "fraction_abs_diff": abs(jp["audit"]["cd4"]["fraction"] - jr["audit"]["cd4"]["fraction"]),
        "preserve_pinned": jp["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"],
        "preserve_refit": jr["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"],
        "union_fraction_pinned": jp["audit"]["h1_union"]["residual_frame"]["solved_fraction"],
        "union_fraction_refit": jr["audit"]["h1_union"]["residual_frame"]["solved_fraction"],
        "deflation_pinned": jp["audit"]["cd4"]["deflation"], "deflation_refit": jr["audit"]["cd4"]["deflation"],
        "frame_pinned": jp["audit"]["cd4"]["residual_frame"], "frame_refit": jr["audit"]["cd4"]["residual_frame"],
    }
    dp, dr = out["deflation_pinned"], out["deflation_refit"]
    d_energy = abs(dp["energy_removed_frac"] - dr["energy_removed_frac"])
    d_capt = abs(dp["captured_frac_of_centred"] - dr["top_k_energy_frac_of_centred"])
    out["energy_removed_abs_diff"] = d_energy
    out["captured_vs_topk_abs_diff"] = d_capt
    numeric = amax <= 1e-12 and out["fraction_abs_diff"] <= 1e-14 and d_energy <= 1e-12 and d_capt <= 1e-12 \
        and dp["n_carried"] == dr["n_carried"] == 291 \
        and abs(out["union_fraction_pinned"] - out["union_fraction_refit"]) == 0.0
    out["row2_byte_identity"] = out["bitwise_equal"]
    out["row2_numeric_identity"] = bool(numeric)
    (HERE / "work" / "compare251.json").write_text(json.dumps(out, indent=1))
    print("=" * 100)
    print("KNOWN-ANSWER ROW 2 -- pinned V5_abc vs 249's re-fit K=5 on A/B/C (staged effect matrices)")
    print("=" * 100)
    print(f"  pinned  array sha256 {out['pinned_array_sha256']}")
    print(f"  re-fit  array sha256 {out['refit_249_array_sha256']}  (249 stage-svd5.json {jr['effects_array_sha256'][:16]}...)")
    if k5 is not None:
        print(f"  k5 (this builder)    {out['k5_this_builder_array_sha256']}  == 249's bitwise: {out['k5_equals_249_bitwise']}")
    print(f"  bitwise equal: {out['bitwise_equal']}")
    print(f"  max |d| {amax:.3e} on entries up to {out['max_abs_entry']:.3f} (relative {rel:.3e}); "
          f"{ndiff:,} of {d.size:,} entries differ, on {rows_diff.size} rows "
          f"({rows_diff_carried} CD4-carried, {out['n_rows_differing_not_carried']} not carried; "
          f"max |d| on non-carried rows {out['max_abs_diff_not_carried_rows']:.3e})")
    print(f"  pinned vs OFF (abeb0749) max |d| {out['pinned_vs_off_max_abs_diff']:.3e} (the flag ACTS)")
    print(f"  cd4.fraction pinned {out['fraction_pinned']!r} vs re-fit {out['fraction_refit']!r} (|d| {out['fraction_abs_diff']:.1e})")
    print(f"  preserve {out['preserve_pinned']!r} vs {out['preserve_refit']!r}; union fraction "
          f"{out['union_fraction_pinned']!r} vs {out['union_fraction_refit']!r}")
    print(f"  deflation pinned: n_carried {dp['n_carried']}, energy_removed {dp['energy_removed_frac']!r}, "
          f"captured {dp['captured_frac_of_centred']!r}, own top-k {dp['own_top_k_energy_frac_of_centred']!r}, "
          f"basis {dp['basis']} file {dp['basis_sha256'][:16]} V {dp['basis_v_sha256'][:16]} U_k {dp['basis_u_k_sha256'][:16]} "
          f"axis verified {dp['gene_axis_verified']}")
    print(f"  deflation re-fit: n_carried {dr['n_carried']}, energy_removed {dr['energy_removed_frac']!r}, "
          f"top-k {dr['top_k_energy_frac_of_centred']!r}, U_k {dr['basis_u_k_sha256'][:16]}")
    print(f"  |d energy_removed| {d_energy:.1e}; |captured - top-k| {d_capt:.1e}")
    print(f"\nROW 2 byte identity: {'PASS' if out['row2_byte_identity'] else 'NO (expected: two associations of one projector)'}"
          f"   ROW 2 numeric identity (<= 1e-12 on the matrix, 1e-14 on the fraction): "
          f"{'PASS' if numeric else 'FAIL'}")
    raise SystemExit(0 if numeric else 1)


if __name__ == "__main__":
    main()
