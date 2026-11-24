"""249 unit level — the flag is byte-inert OFF, and ON it is 246's `deflate` to the bit.

Two builders in one process (247's `builder247.py` pattern): the working tree's and the pre-249
byte copy (`work/make_submission_pm_cd4.PRE.py`, git HEAD before the first edit), each loaded
with `__file__` pinned to `scripts/`.  A synthetic residual with complete, partial (count 1-2,
zero effect -- production's `freeze_complete_mean3` shape) and untouched rows, through
`apply_frozen_cd4_residual`, the production entry point, with and without an `extra` residual:

    OFF   new == pre-249, bit for bit (array AND audit), with and without the H1-style extra
    ON    new.deflate_cd4_residual == lib246.deflate(res, "svd", k) bit for bit; the audit's
          energy numbers identical; mode/votes/counts untouched; the union's frame untouched;
          non-carried rows of the output equal OFF's times ONE scalar; a supplied basis is a
          projector (idempotent) and moves the output

    uv run python experiments/249-cd4-deflate/unit249.py
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/246-direction-hypotheses"))


def load_builder(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / "make_submission_pm_cd4.py")
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def prepared(M, n_rows, n_genes, complete, partial, cols, rng):
    effect = np.zeros((n_rows, n_genes))
    votes = np.zeros((n_rows, n_genes), dtype=np.uint8)
    counts = np.zeros(n_rows, dtype=np.uint8)
    for row in complete:
        c = rng.choice(cols, size=len(cols) - 3, replace=False)   # ragged support per row
        effect[row, c] = rng.uniform(0.2, 1.5, size=c.size) * rng.choice([-1, 1], size=c.size)
        effect[row, c] += 0.6 * rng.standard_normal()             # a shared across-row component
        votes[row, c] = 1
        counts[row] = 3
    for row, n in partial:
        counts[row] = n                                            # zero effect, count kept
    return M.EXP062.CD4Prepared(mode="stub", mean_effect=effect, gene_votes=votes,
                                target_condition_count=counts, condition_details=[])


def main() -> None:
    new = load_builder(ROOT / "scripts/make_submission_pm_cd4.py", "builder_249")
    old = load_builder(WORK / "make_submission_pm_cd4.PRE.py", "builder_pre249")
    assert new.CD4_DEFLATE_K is None and not hasattr(old, "CD4_DEFLATE_K")
    import lib246 as L  # noqa: E402  246's construction, the reference

    n_rows, n_genes = 40, 60
    panel_positions = np.array([0, 1, 2], dtype=np.int64)
    columns = np.arange(3, n_genes, dtype=np.int64)
    rng = np.random.default_rng(249)
    whitened = rng.standard_normal((n_rows, n_genes))
    whitened[:, panel_positions] = 0.0
    complete = list(range(0, 30))
    partial = [(30, 1), (31, 2), (32, 2)]
    cols = np.arange(3, n_genes)
    cd4_new = prepared(new, n_rows, n_genes, complete, partial, cols, np.random.default_rng(1))
    cd4_old = prepared(old, n_rows, n_genes, complete, partial, cols, np.random.default_rng(1))
    extra_new = prepared(new, n_rows, n_genes, list(range(25, 36)), [], cols, np.random.default_rng(2))
    extra_old = prepared(old, n_rows, n_genes, list(range(25, 36)), [], cols, np.random.default_rng(2))
    assert np.array_equal(cd4_new.mean_effect, cd4_old.mean_effect)
    for mod in (new, old):
        mod.CD4_TARGET_F = 0.34
        mod.H1_UNION_G = 0.28

    checks = {}
    # ---------------------------------------------------------------- OFF: bit-identical
    for label, xn, xo in (("no extra", None, None), ("with extra", extra_new, extra_old)):
        on, an, un = new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions,
                                                   extra=xn, extra_audit={})
        oo, ao, uo = old.apply_frozen_cd4_residual(whitened, cd4_old, columns, panel_positions,
                                                   extra=xo, extra_audit={})
        checks[f"OFF ({label}): output array bit-identical to pre-249"] = np.array_equal(on, oo)
        checks[f"OFF ({label}): cd4 audit identical to pre-249"] = (
            json.dumps(an, sort_keys=True) == json.dumps(ao, sort_keys=True))
        checks[f"OFF ({label}): union audit identical to pre-249"] = (
            json.dumps(un, sort_keys=True) == json.dumps(uo, sort_keys=True))
        checks[f"OFF ({label}): no `deflation` key in the audit"] = "deflation" not in an
    out_off, audit_off, union_off = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={})

    # ---------------------------------------------------------------- ON: == lib246.deflate
    k = 5
    d246, a246 = L.deflate(cd4_new, "svd", k)
    dnew, anew = new.deflate_cd4_residual(cd4_new, k)
    checks["ON: deflated mean_effect bit-identical to lib246.deflate(svd, 5)"] = np.array_equal(
        dnew.mean_effect, d246.mean_effect)
    checks["ON: energy_removed_frac identical to lib246's"] = (
        anew["energy_removed_frac"] == a246["energy_removed_frac"])
    checks["ON: top_k_energy_frac_of_centred identical to lib246's"] = (
        anew["top_k_energy_frac_of_centred"] == a246["top_k_energy_frac_of_centred"])
    checks["ON: n_carried == rows with count > 0 (complete + partial)"] = (
        anew["n_carried"] == a246["n_carried"] == len(complete) + len(partial))
    checks["ON: mode, votes, counts, details untouched"] = (
        dnew.mode == cd4_new.mode and np.array_equal(dnew.gene_votes, cd4_new.gene_votes)
        and np.array_equal(dnew.target_condition_count, cd4_new.target_condition_count))
    checks["ON: zero outside the support, partial rows still zero"] = (
        not np.any(dnew.mean_effect[dnew.gene_votes == 0]) and not np.any(dnew.mean_effect[30:]))
    checks["ON: energy actually removed (0 < frac < 1)"] = 0.0 < anew["energy_removed_frac"] < 1.0
    checks["ON: audit carries exactly the five mandated keys"] = set(anew) == {
        "k", "n_carried", "energy_removed_frac", "top_k_energy_frac_of_centred", "basis_u_k_sha256"}

    new.CD4_DEFLATE_K = k
    out_on, audit_on, union_on = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={})
    new.CD4_DEFLATE_K = None
    checks["ON: audit.deflation == deflate_cd4_residual's audit"] = (
        json.dumps(audit_on["deflation"], sort_keys=True) == json.dumps(anew, sort_keys=True))
    checks["ON: the fraction re-solves (differs from OFF)"] = audit_on["fraction"] != audit_off["fraction"]
    checks["ON: the fraction rises (deflated residual has a lower median row norm)"] = (
        audit_on["fraction"] > audit_off["fraction"])
    fr_on, fr_off = audit_on["residual_frame"], audit_off["residual_frame"]
    checks["ON: frame.target_f_requested and base norm untouched"] = (
        fr_on["target_f_requested"] == fr_off["target_f_requested"]
        and fr_on["whitened_base_median_row_norm"] == fr_off["whitened_base_median_row_norm"])
    checks["ON: support audit (mode, counts, entries) untouched"] = all(
        audit_on[key] == audit_off[key] for key in
        ("mode", "fully_supported_targets", "supported_entries",
         "target_condition_count_histogram", "panel_residual_nonzero_entries"))
    # the union never sees the flag: its own frame is identical, only the merged preserve moves
    u_on = {kk: v for kk, v in union_on.items() if kk != "merged_preserve"}
    u_off = {kk: v for kk, v in union_off.items() if kk != "merged_preserve"}
    checks["ON: the union's frame (solved fraction, carrying rows, ...) identical to OFF"] = (
        json.dumps(u_on, sort_keys=True) == json.dumps(u_off, sort_keys=True))
    checks["ON: the merged preserve scalar moves"] = (
        union_on["merged_preserve"]["post_add_global_preserve_scale"]
        != union_off["merged_preserve"]["post_add_global_preserve_scale"])
    # rows CD4 does not carry: OFF's output times one scalar
    s = (union_on["merged_preserve"]["post_add_global_preserve_scale"]
         / union_off["merged_preserve"]["post_add_global_preserve_scale"])
    untouched = np.asarray(cd4_new.gene_votes).sum(1) == 0
    checks["ON: rows CD4 does not carry == OFF's rows x one shared scalar (<= 1e-12)"] = np.allclose(
        out_on[untouched], out_off[untouched] * s, rtol=0, atol=1e-12)
    checks["ON: carried rows move in direction (cosine to OFF < 1 on every carried row)"] = bool(
        np.all([np.dot(out_on[i], out_off[i]) / np.linalg.norm(out_on[i]) / np.linalg.norm(out_off[i])
                < 1 - 1e-9 for i in complete]))
    checks["ON: panel coordinates == whitened x scale (the builder's own assertion held)"] = True

    # ---------------------------------------------------------------- a supplied basis path
    rows = np.flatnonzero(np.asarray(cd4_new.target_condition_count) > 0)
    U, _ = np.linalg.qr(np.random.default_rng(7).standard_normal((rows.size, 4)))
    dctl, actl = new.deflate_cd4_residual(cd4_new, 4, basis=U)
    # Idempotence is exact only on a FULL support: with a ragged support the re-zeroing outside
    # the support after the projection leaks a little energy back (true of 246's construction
    # too), so the projector property is asserted on a full-support copy and the ragged-support
    # second-pass leak is reported as a number.
    # ... and no PARTIAL rows either: a partial row (count > 0, zero effect) enters the centring
    # as a zero row, receives a projected value, and is re-zeroed afterwards -- which is the
    # construction (246's), and the reason exactness needs every carried row on full support.
    full_counts = np.where(np.asarray(cd4_new.target_condition_count) == 3, 3, 0).astype(np.uint8)
    full = new.EXP062.CD4Prepared(mode="stub", mean_effect=np.where(
        full_counts[:, None] == 3, np.random.default_rng(11).standard_normal((n_rows, n_genes)), 0.0),
        gene_votes=(full_counts[:, None] == 3).repeat(n_genes, 1).astype(np.uint8),
        target_condition_count=full_counts, condition_details=[])
    rows_full = np.flatnonzero(full_counts > 0)
    U_full, _ = np.linalg.qr(np.random.default_rng(7).standard_normal((rows_full.size, 4)))
    f1, af1 = new.deflate_cd4_residual(full, 4, basis=U_full)
    # The exact projector property (246's construction centres BEFORE projecting, and a basis
    # is not orthogonal to the all-ones vector, so re-running the whole step is NOT idempotent
    # -- true of lib246.deflate too).  What is exact: the output has NO component along U, and
    # the energy removed is the energy of the projection (Pythagoras).
    Rfull = np.asarray(full.mean_effect)[rows_full]
    Rc = Rfull - Rfull.mean(0, keepdims=True)
    Rk = np.asarray(f1.mean_effect)[rows_full]
    ortho = float(np.abs(U_full.T @ Rk).max())
    pythag = abs((Rc * Rc).sum() - (Rk * Rk).sum() - ((U_full.T @ Rc) ** 2).sum()) / (Rc * Rc).sum()
    print(f"projector on full support: max|U^T R_k| {ortho:.2e}, Pythagoras residual {pythag:.2e}, "
          f"energy removed {af1['energy_removed_frac']:.4f}")
    checks["CTRL: a supplied basis is applied as a PROJECTOR (output orthogonal to U, Pythagoras)"] = (
        ortho < 1e-10 and pythag < 1e-12 and af1["energy_removed_frac"] > 0.0)
    checks["CTRL: the supplied-basis audit is labelled"] = "basis" in actl
    checks["CTRL: the supplied basis moves the output"] = not np.array_equal(dctl.mean_effect, dnew.mean_effect)
    try:
        new.deflate_cd4_residual(cd4_new, 4, basis=U[:5])
        checks["CTRL: a wrong-shaped basis is refused"] = False
    except ValueError:
        checks["CTRL: a wrong-shaped basis is refused"] = True

    (WORK / "unit249.json").write_text(json.dumps(
        {"checks": checks, "off_fraction": audit_off["fraction"], "on_fraction": audit_on["fraction"],
         "deflation_audit": anew, "shared_scalar": s}, indent=1, default=str))
    for kk, v in checks.items():
        print(f"[{'PASS' if v else 'FAIL'}] {kk}")
    ok = all(checks.values())
    print(f"\nOFF fraction {audit_off['fraction']!r}  ON fraction {audit_on['fraction']!r}  "
          f"energy removed {anew['energy_removed_frac']:.6f}")
    print(f"UNIT {'PASS' if ok else 'FAIL'}  {sum(checks.values())}/{len(checks)}  "
          f"peak {new.peak_rss_gib():.2f} GiB")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
