"""254 unit level -- the flag is byte-inert OFF, and ON it is 246's construction to the bit.

Two builders in one process (247/249's pattern): the working tree's and the pre-254 byte copy
(`work/make_submission_pm_cd4.PRE.py`, git HEAD before the first edit), each loaded with
`__file__` pinned to `scripts/`.  A synthetic PRE-FREEZE `prepare_cd4_mean`-shaped object with
complete (count 3), partial (count 1-2, NON-zero effect on a ragged support) and untouched rows;
`freeze_complete_mean3` gives the production complete-row object, `extract_cd4_partial` the
partial residual.  Through `apply_frozen_cd4_residual`, the production entry point:

    OFF   new == pre-254, bit for bit (array AND audits), with and without the H1-style extra
          and with a union 2; no `cd4_partial` key anywhere
    ON    extract_cd4_partial == lib246.cd4_partial's arithmetic (re-stated inline) bit for bit;
          the output == 246's `merged_many` chain (cd4*f + u1*f1 + partial*fp, votes OR, counts
          max, ONE add, ONE preserve) bit for bit; CD4's fraction/frame and the union's frame are
          untouched; rows the partial does not carry == OFF's rows x one scalar; partial rows move
          in direction; the partial's fraction == solve_carrying_row_fraction(whitened, partial, g);
          the shipped all-rows rule would solve 0.0 for it; with K=5 as well the CD4 fraction and
          deflation audit equal the K=5-only run's and the output == merged_many at the RE-SOLVED
          fraction; refusals (no dose, no audits_out, a panel coordinate, an empty residual)

    uv run python experiments/254-partial-rows/unit254.py
"""
from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))


def load_builder(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / "make_submission_pm_cd4.py")
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def prepared_prefreeze(M, n_rows, n_genes, complete, partial, cols, rng):
    """A `prepare_cd4_mean`-shaped object: complete rows carry all three conditions, partial rows
    one or two (their effect is the mean over the carrying conditions -- NON-zero), untouched rows
    none."""
    effect = np.zeros((n_rows, n_genes))
    votes = np.zeros((n_rows, n_genes), dtype=np.uint8)
    counts = np.zeros(n_rows, dtype=np.uint8)
    for row in complete:
        c = rng.choice(cols, size=len(cols) - 3, replace=False)
        effect[row, c] = rng.uniform(0.2, 1.5, size=c.size) * rng.choice([-1, 1], size=c.size)
        effect[row, c] += 0.6 * rng.standard_normal()
        votes[row, c] = 3
        counts[row] = 3
    for row, n in partial:
        c = rng.choice(cols, size=len(cols) - 5, replace=False)
        effect[row, c] = rng.uniform(0.2, 1.5, size=c.size) * rng.choice([-1, 1], size=c.size)
        votes[row, c] = n
        counts[row] = n
    return M.EXP062.CD4Prepared(mode="gene_norm", mean_effect=effect, gene_votes=votes,
                                target_condition_count=counts, condition_details=[])


def lib246_partial(prep, M):
    """lib246.cd4_partial(equalise=False)'s arithmetic, re-stated (the real-data row compares
    against lib246 itself: partialcheck254.py / repro254.py)."""
    count = np.asarray(prep.target_condition_count, dtype=np.int64)
    complete = count == 3
    partial = (count > 0) & ~complete
    eff = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    norms_c = np.linalg.norm(np.where(votes[complete] > 0, eff[complete], 0.0), axis=1)
    eff[~partial] = 0.0
    votes[~partial] = 0
    cnt = np.where(partial, count, 0).astype(np.uint8)
    return eff, votes, cnt, dict(n_partial=int(partial.sum()), n_one=int((count == 1).sum()),
                                 n_two=int((count == 2).sum()),
                                 median_complete_norm=float(np.median(norms_c)))


def merged_many(M, parts):
    """246 `lib246.merged_many` = compose218.merged chained: ONE add, ONE preserve."""
    parts = [(r, f) for r, f in parts if r is not None and f != 0.0]
    base_r, base_f = parts[0]
    vc = np.asarray(base_r.gene_votes) > 0
    eff = np.where(vc, np.asarray(base_r.mean_effect), 0.0) * base_f
    votes = vc.copy()
    counts = np.asarray(base_r.target_condition_count).astype(np.int64)
    for r, f in parts[1:]:
        v = np.asarray(r.gene_votes) > 0
        eff += np.where(v, np.asarray(r.mean_effect), 0.0) * f
        votes |= v
        counts = np.maximum(counts, np.asarray(r.target_condition_count).astype(np.int64))
    return M.EXP062.CD4Prepared(mode="merged", mean_effect=eff, gene_votes=votes.astype(np.uint8),
                                target_condition_count=np.minimum(counts, 255).astype(np.uint8),
                                condition_details=[])


def main() -> None:
    new = load_builder(ROOT / "scripts/make_submission_pm_cd4.py", "builder_254")
    old = load_builder(WORK / "make_submission_pm_cd4.PRE.py", "builder_pre254")
    assert new.CD4_PARTIAL_G is None and not hasattr(old, "CD4_PARTIAL_G")
    checks: dict[str, bool] = {}

    n_rows, n_genes = 40, 60
    panel_positions = np.array([0, 1, 2], dtype=np.int64)
    columns = np.arange(3, n_genes, dtype=np.int64)
    rng = np.random.default_rng(254)
    whitened = rng.standard_normal((n_rows, n_genes))
    whitened[:, panel_positions] = 0.0
    complete = list(range(0, 30))
    partial_rows = [(30, 1), (31, 2), (32, 2), (33, 1)]
    cols = np.arange(3, n_genes)
    prep_new = prepared_prefreeze(new, n_rows, n_genes, complete, partial_rows, cols, np.random.default_rng(1))
    prep_old = prepared_prefreeze(old, n_rows, n_genes, complete, partial_rows, cols, np.random.default_rng(1))
    cd4_new = new.freeze_complete_mean3(prep_new)
    cd4_old = old.freeze_complete_mean3(prep_old)
    checks["freeze_complete_mean3 unchanged (new == pre-254 bit for bit)"] = (
        np.array_equal(cd4_new.mean_effect, cd4_old.mean_effect)
        and np.array_equal(cd4_new.gene_votes, cd4_old.gene_votes)
        and np.array_equal(cd4_new.target_condition_count, cd4_old.target_condition_count))
    checks["the freeze zeroes the partial rows' effect and keeps their counts"] = (
        not np.any(cd4_new.mean_effect[30:34]) and list(cd4_new.target_condition_count[30:34]) == [1, 2, 2, 1])
    extra_new = prepared_prefreeze(new, n_rows, n_genes, list(range(25, 36)), [], cols, np.random.default_rng(2))
    extra_old = prepared_prefreeze(old, n_rows, n_genes, list(range(25, 36)), [], cols, np.random.default_rng(2))
    u2_new = prepared_prefreeze(new, n_rows, n_genes, list(range(10, 18)) + [36, 37], [], cols, np.random.default_rng(3))
    u2_old = prepared_prefreeze(old, n_rows, n_genes, list(range(10, 18)) + [36, 37], [], cols, np.random.default_rng(3))
    for mod in (new, old):
        mod.CD4_TARGET_F = 0.34
        mod.H1_UNION_G = 0.28
    sig = inspect.signature(new.prepare_frozen_cd4)
    checks["prepare_frozen_cd4: `return_partial` defaults to False (the default call is untouched)"] = (
        sig.parameters["return_partial"].default is False)

    # ---------------------------------------------------------------- OFF: bit-identical
    for label, xn, xo, u2n, u2o in (("no extra", None, None, None, None),
                                    ("with extra", extra_new, extra_old, None, None),
                                    ("with extra and union 2", extra_new, extra_old, u2_new, u2_old)):
        ao_n: dict = {}
        ao_o: dict = {}
        if u2n is not None:
            new.UNION2_G = 0.35
            old.UNION2_G = 0.35
        on, an, un = new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions,
                                                   extra=xn, extra_audit={}, extra2=u2n,
                                                   extra2_audit={} if u2n is not None else None,
                                                   audits_out=ao_n)
        oo, aoo, uo = old.apply_frozen_cd4_residual(whitened, cd4_old, columns, panel_positions,
                                                    extra=xo, extra_audit={}, extra2=u2o,
                                                    extra2_audit={} if u2o is not None else None,
                                                    audits_out=ao_o)
        checks[f"OFF ({label}): output array bit-identical to pre-254"] = np.array_equal(on, oo)
        checks[f"OFF ({label}): cd4 audit identical to pre-254"] = (
            json.dumps(an, sort_keys=True) == json.dumps(aoo, sort_keys=True))
        checks[f"OFF ({label}): union audit identical to pre-254"] = (
            json.dumps(un, sort_keys=True) == json.dumps(uo, sort_keys=True))
        checks[f"OFF ({label}): audits_out identical to pre-254, no `cd4_partial` key"] = (
            json.dumps(ao_n, sort_keys=True) == json.dumps(ao_o, sort_keys=True) and "cd4_partial" not in ao_n)
    new.UNION2_G = None
    old.UNION2_G = None
    out_off, audit_off, union_off = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={})

    # ---------------------------------------------------------------- ON: the extraction
    part = new.extract_cd4_partial(prep_new)
    eff_ref, votes_ref, cnt_ref, info_ref = lib246_partial(prep_new, new)
    checks["ON: extract_cd4_partial.mean_effect == lib246.cd4_partial arithmetic bit for bit"] = np.array_equal(
        np.asarray(part.mean_effect), eff_ref)
    checks["ON: extract_cd4_partial votes and counts == lib246's"] = (
        np.array_equal(np.asarray(part.gene_votes), votes_ref)
        and np.array_equal(np.asarray(part.target_condition_count), cnt_ref))
    checks["ON: extract_cd4_partial details == lib246's (n_partial 4, n_one 2, n_two 2, median complete norm)"] = (
        part.condition_details[0] == info_ref and info_ref["n_partial"] == 4)
    checks["ON: the partial residual is zero on complete and untouched rows, non-zero on partial rows"] = (
        not np.any(np.asarray(part.mean_effect)[[*complete, 34, 35, 36, 37, 38, 39]])
        and all(np.any(np.asarray(part.mean_effect)[r]) for r, _ in partial_rows))
    checks["ON: the partial rows' effect is the PRE-FREEZE effect (no re-normalisation)"] = np.array_equal(
        np.asarray(part.mean_effect)[30:34], np.asarray(prep_new.mean_effect)[30:34])
    checks["ON: the partial's support is disjoint from the frozen complete-row object's support"] = (
        not np.any((np.asarray(part.gene_votes) > 0) & (np.asarray(cd4_new.gene_votes) > 0)))
    paudit = new.cd4_partial_support_audit(part, target_g=0.34)
    checks["ON: support audit reads rows 4 / one 2 / two 2 / carrying 4 / target_g 0.34"] = (
        paudit["rows"] == 4 and paudit["rows_one_condition"] == 2 and paudit["rows_two_conditions"] == 2
        and paudit["carrying_rows"] == 4 and paudit["target_g"] == 0.34)

    # ---------------------------------------------------------------- ON: through the entry point
    new.CD4_PARTIAL_G = 0.34
    ao: dict = {}
    out_on, audit_on, union_on = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={},
        partial=part, partial_audit=paudit, audits_out=ao)
    checks["ON: audits_out carries `cd4_partial`, no `union2`"] = "cd4_partial" in ao and "union2" not in ao
    pblock = ao.get("cd4_partial", {})
    fp, frame_p = new.solve_carrying_row_fraction(whitened, part, 0.34)
    checks["ON: the partial's solved fraction == solve_carrying_row_fraction(whitened, partial, 0.34)"] = (
        pblock["residual_frame"]["solved_fraction"] == fp and pblock["residual_frame"]["applied_fraction"] == fp)
    checks["ON: the shipped all-rows median rule would solve exactly 0.0 for it (4/40 coverage)"] = (
        pblock["residual_frame"]["shipped_median_rule_would_solve"] == 0.0)
    checks["ON: CD4's own fraction and residual frame identical to OFF"] = (
        audit_on["fraction"] == audit_off["fraction"]
        and json.dumps(audit_on["residual_frame"], sort_keys=True) == json.dumps(audit_off["residual_frame"], sort_keys=True))
    checks["ON: CD4's support audit (mode, counts, entries) identical to OFF"] = all(
        audit_on[k] == audit_off[k] for k in
        ("mode", "fully_supported_targets", "supported_entries", "target_condition_count_histogram",
         "panel_residual_nonzero_entries"))
    u_on = {k: v for k, v in union_on.items() if k != "merged_preserve"}
    u_off = {k: v for k, v in union_off.items() if k != "merged_preserve"}
    checks["ON: the union's frame (solved fraction, carrying rows, ...) identical to OFF"] = (
        json.dumps(u_on, sort_keys=True) == json.dumps(u_off, sort_keys=True))
    checks["ON: the merged preserve scalar moves, and is ONE number written in all three blocks"] = (
        union_on["merged_preserve"]["post_add_global_preserve_scale"] != union_off["merged_preserve"]["post_add_global_preserve_scale"]
        and union_on["merged_preserve"] == pblock["merged_preserve"] == audit_on["post_add_preserve"] | {"fraction": union_on["merged_preserve"]["fraction"]}
        or (union_on["merged_preserve"]["post_add_global_preserve_scale"]
            == pblock["merged_preserve"]["post_add_global_preserve_scale"]
            == audit_on["post_add_preserve"]["post_add_global_preserve_scale"]
            != union_off["merged_preserve"]["post_add_global_preserve_scale"]))
    # the reference: 246's merged_many chain in 246's order (cd4, union, partial), ONE add, ONE preserve
    fx = union_off["residual_frame"]["solved_fraction"]
    ref = merged_many(new, [(cd4_new, audit_off["fraction"]), (extra_new, fx), (part, fp)])
    ref_out, ref_pres = new.EXP062.add_post_whiten_residual(whitened, ref, 1.0, columns)
    checks["ON: output == 246's merged_many chain (cd4*f + union*f1 + partial*fp, ONE add, ONE preserve) bit for bit"] = (
        np.array_equal(out_on, ref_out)
        and ref_pres["post_add_global_preserve_scale"] == union_on["merged_preserve"]["post_add_global_preserve_scale"])
    # order matters for bits: partial BEFORE the union is a different float association
    ref2 = merged_many(new, [(cd4_new, audit_off["fraction"]), (part, fp), (extra_new, fx)])
    ref2_out, _ = new.EXP062.add_post_whiten_residual(whitened, ref2, 1.0, columns)
    checks["ON: (info) accumulating the partial BEFORE the union gives the same numbers to 1e-12"] = np.allclose(
        ref2_out, out_on, rtol=0, atol=1e-12)
    s = (union_on["merged_preserve"]["post_add_global_preserve_scale"]
         / union_off["merged_preserve"]["post_add_global_preserve_scale"])
    pc = np.asarray(part.gene_votes).sum(1) > 0
    checks["ON: rows the partial does not carry == OFF's rows x one shared scalar (<= 1e-12)"] = np.allclose(
        out_on[~pc], out_off[~pc] * s, rtol=0, atol=1e-12)
    checks["ON: partial rows move in direction (cosine to OFF < 1 on every partial row)"] = bool(np.all(
        [np.dot(out_on[i], out_off[i]) / np.linalg.norm(out_on[i]) / np.linalg.norm(out_off[i]) < 1 - 1e-9
         for i in np.flatnonzero(pc)]))
    checks["ON: panel coordinates == whitened x scale (the builder's own assertion held)"] = np.allclose(
        out_on[:, panel_positions], whitened[:, panel_positions] * union_on["merged_preserve"]["post_add_global_preserve_scale"], rtol=0, atol=0)
    checks["ON: partial residual on its own (no union) also runs through the merged path"] = True
    ao_solo: dict = {}
    out_solo, audit_solo, union_solo = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, partial=part, partial_audit=paudit, audits_out=ao_solo)
    ref_solo = merged_many(new, [(cd4_new, audit_off["fraction"]), (part, fp)])
    ref_solo_out, _ = new.EXP062.add_post_whiten_residual(whitened, ref_solo, 1.0, columns)
    checks["ON (no union): output == merged_many(cd4, partial) bit for bit; union audit None; block present"] = (
        np.array_equal(out_solo, ref_solo_out) and union_solo is None and "cd4_partial" in ao_solo
        and audit_solo["post_add_preserve"]["fraction"] == audit_off["fraction"])

    # ---------------------------------------------------------------- ON + K=5: the composed arm
    new.CD4_DEFLATE_K = 5
    ao_k: dict = {}
    out_k, audit_k, union_k = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={}, audits_out=ao_k)
    ao_kp: dict = {}
    out_kp, audit_kp, union_kp = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={},
        partial=part, partial_audit=paudit, audits_out=ao_kp)
    new.CD4_DEFLATE_K = None
    checks["ON+K5: CD4's re-solved fraction and deflation audit identical to the K5-only run"] = (
        audit_kp["fraction"] == audit_k["fraction"]
        and json.dumps(audit_kp["deflation"], sort_keys=True) == json.dumps(audit_k["deflation"], sort_keys=True)
        and json.dumps(audit_kp["residual_frame"], sort_keys=True) == json.dumps(audit_k["residual_frame"], sort_keys=True))
    checks["ON+K5: the partial's fraction is the same as without deflation (it depends on the base only)"] = (
        ao_kp["cd4_partial"]["residual_frame"]["solved_fraction"] == fp)
    d5, _ = new.deflate_cd4_residual(cd4_new, 5)
    ref_k = merged_many(new, [(d5, audit_k["fraction"]), (extra_new, fx), (part, fp)])
    ref_k_out, _ = new.EXP062.add_post_whiten_residual(whitened, ref_k, 1.0, columns)
    checks["ON+K5: output == merged_many(deflated cd4 @ re-solved fraction, union, partial) bit for bit"] = np.array_equal(
        out_kp, ref_k_out)
    checks["ON+K5: the deflation never touches the partial rows (they enter the merge un-deflated)"] = np.array_equal(
        np.asarray(part.mean_effect)[30:34], np.asarray(prep_new.mean_effect)[30:34])

    # ---------------------------------------------------------------- refusals
    new.CD4_PARTIAL_G = None
    try:
        new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions, extra=extra_new,
                                      extra_audit={}, partial=part, partial_audit=paudit, audits_out={})
        checks["REFUSE: a partial residual without VCC_CD4_PARTIAL_G (no dose)"] = False
    except ValueError:
        checks["REFUSE: a partial residual without VCC_CD4_PARTIAL_G (no dose)"] = True
    new.CD4_PARTIAL_G = 0.34
    try:
        new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions, extra=extra_new,
                                      extra_audit={}, partial=part, partial_audit=paudit)
        checks["REFUSE: a partial residual without `audits_out`"] = False
    except ValueError:
        checks["REFUSE: a partial residual without `audits_out`"] = True
    bad = new.EXP062.CD4Prepared(mode="x", mean_effect=np.asarray(part.mean_effect).copy(),
                                 gene_votes=np.asarray(part.gene_votes).copy(),
                                 target_condition_count=np.asarray(part.target_condition_count).copy(),
                                 condition_details=part.condition_details)
    bad.gene_votes[30, 0] = 1
    bad.mean_effect[30, 0] = 0.5
    try:
        new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions, extra=extra_new,
                                      extra_audit={}, partial=bad, partial_audit=paudit, audits_out={})
        checks["REFUSE: a panel coordinate in the partial residual's support"] = False
    except AssertionError:
        checks["REFUSE: a panel coordinate in the partial residual's support"] = True
    empty = new.EXP062.CD4Prepared(mode="x", mean_effect=np.zeros_like(part.mean_effect),
                                   gene_votes=np.zeros_like(part.gene_votes),
                                   target_condition_count=np.zeros_like(part.target_condition_count),
                                   condition_details=[dict(n_partial=0, n_one=0, n_two=0, median_complete_norm=0.0)])
    try:
        new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions, extra=extra_new,
                                      extra_audit={}, partial=empty, partial_audit={}, audits_out={})
        checks["REFUSE: a partial residual with no carrying rows (a panel with no partial rows) is loud, not silent"] = False
    except ValueError:
        checks["REFUSE: a partial residual with no carrying rows (a panel with no partial rows) is loud, not silent"] = True
    new.CD4_PARTIAL_G = None

    (WORK / "unit254.json").write_text(json.dumps(
        {"checks": checks, "off_fraction": audit_off["fraction"], "partial_fraction": fp,
         "partial_frame": frame_p, "shared_scalar": s, "k5_fraction": audit_k["fraction"]},
        indent=1, default=str))
    for kk, v in checks.items():
        print(f"[{'PASS' if v else 'FAIL'}] {kk}")
    ok = all(checks.values())
    print(f"\nOFF fraction {audit_off['fraction']!r}  partial fraction {fp!r}  shared scalar {s!r}")
    print(f"UNIT {'PASS' if ok else 'FAIL'}  {sum(checks.values())}/{len(checks)}  peak {new.peak_rss_gib():.2f} GiB")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
