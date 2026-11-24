"""417 builder window 2c: VCC_CARRIER_GAMMA -- the DE CARRIERS from a second whitening pass at a different gamma (416's
gXXX construction: apply_v9_whitening(candidate_pre, gamma) -> the same frozen CD4/H1 residuals -> each row rescaled to the
shipped row's norm). The bulk (pds/mse, through the ballast) keeps the gamma-0.65 matrix untouched; only the carriers
(fid/reach/jac/nmae) change. Requires VCC_PERT_SCALE_BALLAST. Unset is a no-op.
Apply only when no build has the builder pinned. Apply AFTER patch_next417.py and patch_bulkio417.py.
usage: uv run python patch_carrier417.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
"""
import difflib
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
P = Path(os.environ.get("PATCH_TARGET", ROOT / "scripts/make_submission_pm_cd4.py"))
src = P.read_text()
out = src


def sub(old, new, label):
    global out
    n = out.count(old)
    if n != 1:
        raise SystemExit(f"[{label}] expected exactly one match, found {n}")
    out = out.replace(old, new)


sub('''if BULK_DUMP is not None and BULK_LOAD is not None:
    raise SystemExit("VCC_BULK_DUMP and VCC_BULK_LOAD are exclusive")
''', '''if BULK_DUMP is not None and BULK_LOAD is not None:
    raise SystemExit("VCC_BULK_DUMP and VCC_BULK_LOAD are exclusive")
# `VCC_CARRIER_GAMMA` (417, from 416's H1 menu): the DE carriers from a second whitening pass at this gamma, the same frozen
# residuals, each row rescaled to the shipped row's norm. The bulk keeps gamma 0.65. Requires VCC_PERT_SCALE_BALLAST.
CARRIER_GAMMA = os.environ.get("VCC_CARRIER_GAMMA", "").strip() or None
if CARRIER_GAMMA is not None:
    CARRIER_GAMMA = float(CARRIER_GAMMA)
    if not 0.0 <= CARRIER_GAMMA <= 3.0:
        raise SystemExit(f"VCC_CARRIER_GAMMA must be in [0, 3] when set: {CARRIER_GAMMA!r}")
    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_CARRIER_GAMMA changes the carriers only and needs VCC_PERT_SCALE_BALLAST (welded it moves the bulk)")
''', "knob")

sub('''    union2_written = extra_audits.get("union2")
    cd4_partial_written = extra_audits.get("cd4_partial")

    if not np.all(np.isfinite(final)):
        raise ValueError("final effect matrix contains non-finite values")
    return final, {''', '''    union2_written = extra_audits.get("union2")
    cd4_partial_written = extra_audits.get("cd4_partial")

    # 417: VCC_CARRIER_GAMMA -- the carriers' own design (416's norm-matched gXXX). Deep-copied audits: the second residual
    # pass must not overwrite the shipped pass's audit blocks.
    carrier_final, carrier_audit = None, None
    if CARRIER_GAMMA is not None:
        import copy  # noqa: PLC0415
        if SIGN_REPAIR is not None:
            raise SystemExit("VCC_CARRIER_GAMMA with VCC_SIGN_REPAIR is not a priced combination")
        wc, wc_audit = apply_v9_whitening(candidate_pre, whitener, operator_columns, gamma=CARRIER_GAMMA)
        if cd4 is None:
            fc = wc
        else:
            fc, _, _ = apply_frozen_cd4_residual(
                wc, cd4, operator_columns, panel_positions,
                extra=h1_union, extra_audit=copy.deepcopy(h1_union_audit),
                extra2=union2, extra2_audit=copy.deepcopy(union2_audit),
                audits_out={}, partial=cd4_partial, partial_audit=copy.deepcopy(cd4_partial_audit),
            )
        n0 = np.linalg.norm(final, axis=1)
        n1 = np.linalg.norm(fc, axis=1)
        scale_rows = np.where(n1 > 0, n0 / np.maximum(n1, 1e-300), 0.0)
        carrier_final = fc * scale_rows[:, None]
        if not np.all(np.isfinite(carrier_final)):
            raise ValueError("carrier effect matrix contains non-finite values")
        row_cos = (carrier_final * final).sum(1) / np.maximum(n0 * np.linalg.norm(carrier_final, axis=1), 1e-300)
        carrier_audit = {"gamma": CARRIER_GAMMA, "whitening": wc_audit,
                         "norm_scale_q": np.quantile(scale_rows, [0, .5, 1]).tolist(),
                         "row_cos_vs_bulk_q": np.quantile(row_cos, [0, .1, .5, .9, 1]).tolist(),
                         "max_abs_diff_vs_bulk": float(np.max(np.abs(carrier_final - final))),
                         "sha256": array_sha256(carrier_final)}
        print(f"[carrier-gamma] gamma={CARRIER_GAMMA} row cos vs the shipped-gamma bulk q0/10/50/90/100="
              f"{np.round(carrier_audit['row_cos_vs_bulk_q'], 4).tolist()} norm scale q={np.round(carrier_audit['norm_scale_q'], 4).tolist()} "
              f"max|carrier-bulk|={carrier_audit['max_abs_diff_vs_bulk']:.3e}", flush=True)

    if not np.all(np.isfinite(final)):
        raise ValueError("final effect matrix contains non-finite values")
    return final, {
        **({} if carrier_final is None else {"carrier_gamma": carrier_audit, "_carrier_matrix": carrier_final}),''', "second_pass")

sub('''        cd4_partial_audit=cd4_partial_audit,
        weights=weights,
    )
    rank_scores = None''', '''        cd4_partial_audit=cd4_partial_audit,
        weights=weights,
    )
    carrier_matrix = effect_audit.pop("_carrier_matrix", None)   # 417: VCC_CARRIER_GAMMA's carriers (None when unset)
    rank_scores = None''', "pop")

sub('''    output_audit = write_submission_from_effects(
        spec,
        effects_unscaled if PERT_SCALE_BALLAST is not None else effects,''', '''    if CARRIER_GAMMA is not None:
        if carrier_matrix is None or PERT_SCALE is None or context_effects is not None:
            raise SystemExit("VCC_CARRIER_GAMMA needs the shared-matrix ballast path and its carrier matrix")
        effects_unscaled = carrier_matrix
    output_audit = write_submission_from_effects(
        spec,
        effects_unscaled if PERT_SCALE_BALLAST is not None else effects,''', "writer")

sub('''            **({} if BULK_LOAD is None else {"bulk_load": str(BULK_LOAD[0].name), "bulk_load_sha256": BULK_LOAD[1]}),
''', '''            **({} if BULK_LOAD is None else {"bulk_load": str(BULK_LOAD[0].name), "bulk_load_sha256": BULK_LOAD[1]}),
            **({} if CARRIER_GAMMA is None else {"carrier_gamma": CARRIER_GAMMA}),
''', "config")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
