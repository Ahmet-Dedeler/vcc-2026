"""248 unit level -- the second union is byte-inert OFF, and ON it is 218's `merged()` chained.

Two builders in one process (247's `builder247.py` / 249's `unit249.py` pattern): the working
tree's and the pre-248 byte copy (`work/make_submission_pm_cd4.PRE.py`, git HEAD `2d63b979...`
before the first edit), each loaded with `__file__` pinned to `scripts/`.  Synthetic residuals
with complete, partial (count kept, zero effect) and untouched rows, through
`apply_frozen_cd4_residual`, the production entry point:

    OFF   new == pre-248 bit for bit (array, cd4 audit, union audit): no extras, union 1 only
    SWAP  union 2 alone (union 1 None), the SAME residual and dose as union 1 alone
          -> output bit-identical to the union-1-alone path (known-answer row 2 at unit level)
    ON    union 1 + union 2 == compose218.merged() chained twice, then ONE add + ONE preserve,
          bit for bit; CD4's fraction and union 1's frame untouched by union 2; union 2's frame
          == solve_carrying_row_fraction(whitened, u2, g2); rows neither union carries equal the
          union-1-only output times ONE scalar; rows only union 2 carries move in direction;
          `audits_out` carries the union2 block with the same keys as union 1's plus the two
          overlap counts; a missing `audits_out` with extra2 is refused; a missing dose is refused

    uv run python experiments/248-second-residual/unit248.py
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
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))


def load_builder(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / "make_submission_pm_cd4.py")
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def prepared(M, n_rows, n_genes, complete, partial, cols, rng, count=3):
    effect = np.zeros((n_rows, n_genes))
    votes = np.zeros((n_rows, n_genes), dtype=np.uint8)
    counts = np.zeros(n_rows, dtype=np.uint8)
    for row in complete:
        c = rng.choice(cols, size=len(cols) - 3, replace=False)
        effect[row, c] = rng.uniform(0.2, 1.5, size=c.size) * rng.choice([-1, 1], size=c.size)
        effect[row, c] += 0.6 * rng.standard_normal()
        votes[row, c] = 1
        counts[row] = count
    for row, n in partial:
        counts[row] = n
    return M.EXP062.CD4Prepared(mode="stub", mean_effect=effect, gene_votes=votes,
                                target_condition_count=counts, condition_details=[])


def same_audit(a, b) -> bool:
    return json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)


def main() -> None:
    new = load_builder(ROOT / "scripts/make_submission_pm_cd4.py", "builder_248")
    old = load_builder(WORK / "make_submission_pm_cd4.PRE.py", "builder_pre248")
    assert new.UNION2_SOURCES is None and new.UNION2_G is None and not hasattr(old, "UNION2_G")
    from compose218 import merged  # noqa: E402  218 section 7.2's construction, the reference

    n_rows, n_genes = 40, 60
    panel_positions = np.array([0, 1, 2], dtype=np.int64)
    columns = np.arange(3, n_genes, dtype=np.int64)
    rng = np.random.default_rng(248)
    whitened = rng.standard_normal((n_rows, n_genes))
    whitened[:, panel_positions] = 0.0
    cols = np.arange(3, n_genes)
    mk = {}
    for tag, mod in (("new", new), ("old", old)):
        mk[tag] = dict(
            cd4=prepared(mod, n_rows, n_genes, list(range(0, 30)), [(30, 1), (31, 2), (32, 2)], cols,
                         np.random.default_rng(1)),
            u1=prepared(mod, n_rows, n_genes, list(range(25, 36)), [], cols, np.random.default_rng(2), count=1),
            # union 2 (rows 33-39) overlaps union 1 on rows 33-35, lies entirely outside CD4 (rows 0-32), and is alone on 36-39
            u2=prepared(mod, n_rows, n_genes, list(range(33, 40)), [], cols, np.random.default_rng(3), count=1),
        )
    for key in ("cd4", "u1", "u2"):
        assert np.array_equal(mk["new"][key].mean_effect, mk["old"][key].mean_effect)
    for mod in (new, old):
        mod.CD4_TARGET_F = 0.34
        mod.H1_UNION_G = 0.28
    new.UNION2_G = None

    checks = {}
    # ---------------------------------------------------------------- OFF: bit-identical
    for label, has_u1 in (("no extra", False), ("union 1 only", True)):
        on, an, un = new.apply_frozen_cd4_residual(
            whitened, mk["new"]["cd4"], columns, panel_positions,
            extra=mk["new"]["u1"] if has_u1 else None, extra_audit={"tag": "u1"} if has_u1 else None)
        oo, ao, uo = old.apply_frozen_cd4_residual(
            whitened, mk["old"]["cd4"], columns, panel_positions,
            extra=mk["old"]["u1"] if has_u1 else None, extra_audit={"tag": "u1"} if has_u1 else None)
        checks[f"OFF ({label}): output array bit-identical to pre-248"] = np.array_equal(on, oo)
        checks[f"OFF ({label}): cd4 audit identical to pre-248"] = same_audit(an, ao)
        checks[f"OFF ({label}): union audit identical to pre-248"] = same_audit(un, uo)
    out_u1, audit_u1, union_u1 = new.apply_frozen_cd4_residual(
        whitened, mk["new"]["cd4"], columns, panel_positions, extra=mk["new"]["u1"], extra_audit={"tag": "u1"})

    # ---------------------------------------------------------------- SWAP: union 2 alone == union 1 alone
    new.UNION2_G = 0.28
    audits = {}
    out_sw, audit_sw, union_sw = new.apply_frozen_cd4_residual(
        whitened, mk["new"]["cd4"], columns, panel_positions,
        extra=None, extra2=mk["new"]["u1"], extra2_audit={"tag": "u1"}, audits_out=audits)
    checks["SWAP: union 2 alone (union 1 unset) == union 1 alone, output bit-identical"] = np.array_equal(out_sw, out_u1)
    checks["SWAP: cd4 audit identical (fraction, frame, preserve)"] = same_audit(audit_sw, audit_u1)
    checks["SWAP: union 1's return slot is None, the block is in audits_out['union2']"] = (
        union_sw is None and "union2" in audits)
    u2_frame = {k: v for k, v in audits["union2"]["residual_frame"].items()}
    u1_frame = dict(union_u1["residual_frame"])
    checks["SWAP: union 2's frame == union 1's frame, key for key"] = same_audit(u2_frame, u1_frame)
    checks["SWAP: union 2's merged preserve == union 1's"] = same_audit(
        audits["union2"]["merged_preserve"], union_u1["merged_preserve"])

    # ---------------------------------------------------------------- ON: union 1 + union 2 == merged() chained
    new.UNION2_G = 0.35
    audits = {}
    out_on, audit_on, union_on = new.apply_frozen_cd4_residual(
        whitened, mk["new"]["cd4"], columns, panel_positions,
        extra=mk["new"]["u1"], extra_audit={"tag": "u1"},
        extra2=mk["new"]["u2"], extra2_audit={"tag": "u2"}, audits_out=audits)
    cd4, u1, u2 = mk["new"]["cd4"], mk["new"]["u1"], mk["new"]["u2"]
    frac_cd4, _ = new.solve_cd4_fraction(whitened, cd4)
    f1, fr1 = new.solve_carrying_row_fraction(whitened, u1, 0.28)
    f2, fr2 = new.solve_carrying_row_fraction(whitened, u2, 0.35)
    chained = merged(merged(cd4, frac_cd4, u1, f1), 1.0, u2, f2)
    ref, ref_preserve = new.EXP062.add_post_whiten_residual(whitened, chained, 1.0, columns)
    checks["ON: output == compose218.merged() chained twice + ONE add + ONE preserve, bit for bit"] = np.array_equal(out_on, ref)
    checks["ON: the one preserve scalar == the chained reference's"] = (
        union_on["merged_preserve"]["post_add_global_preserve_scale"] == ref_preserve["post_add_global_preserve_scale"]
        == audits["union2"]["merged_preserve"]["post_add_global_preserve_scale"]
        == audit_on["post_add_preserve"]["post_add_global_preserve_scale"])
    checks["ON: CD4's solved fraction untouched by union 2"] = audit_on["fraction"] == audit_u1["fraction"] == frac_cd4
    frame1_on = {k: v for k, v in union_on["residual_frame"].items()}
    checks["ON: union 1's frame untouched by union 2 (solved fraction, carrying rows, ...)"] = same_audit(
        frame1_on, union_u1["residual_frame"])
    u2f = audits["union2"]["residual_frame"]
    checks["ON: union 2's solved fraction == solve_carrying_row_fraction(whitened, u2, 0.35)"] = (
        u2f["solved_fraction"] == f2 == u2f["applied_fraction"] and u2f["target_g_requested"] == 0.35)
    base_keys = set(union_u1["residual_frame"])
    checks["ON: union 2's frame has union 1's keys plus exactly the two overlap counts"] = (
        set(u2f) == base_keys | {"targets_shared_with_h1_union", "entries_shared_with_h1_union"})
    checks["ON: union 2's block has union 1's block keys"] = set(audits["union2"]) == set(union_on)
    checks["ON: overlap counts are right (3 rows shared with union 1, 7 outside CD4, 7 carried)"] = (
        u2f["targets_shared_with_h1_union"] == 3 and u2f["targets_outside_cd4_support"] == 7
        and u2f["carrying_rows"] == 7)
    s = (union_on["merged_preserve"]["post_add_global_preserve_scale"]
         / union_u1["merged_preserve"]["post_add_global_preserve_scale"])
    only_u2 = np.asarray(u2.gene_votes).sum(1) > 0
    untouched = ~only_u2
    checks["ON: rows union 2 does not carry == union-1-only output x ONE scalar (<= 1e-12)"] = np.allclose(
        out_on[untouched], out_u1[untouched] * s, rtol=0, atol=1e-12)
    checks["ON: rows union 2 carries move in direction (cosine to union-1-only < 1)"] = bool(np.all(
        [np.dot(out_on[i], out_u1[i]) / np.linalg.norm(out_on[i]) / np.linalg.norm(out_u1[i]) < 1 - 1e-9
         for i in np.flatnonzero(only_u2)]))
    checks["ON: panel coordinates == whitened x scale on every row"] = np.allclose(
        out_on[:, panel_positions], whitened[:, panel_positions] * union_on["merged_preserve"]["post_add_global_preserve_scale"], rtol=0, atol=0)
    checks["ON: the union-2 residual is not vacuous (fraction > 0, carrying rows > 0)"] = f2 > 0 and u2f["carrying_rows"] > 0
    # the same arm, built through the old builder with a hand-made union == this? no: the old
    # builder cannot express it (229 finding 1).  What it CAN do is union 1 = (u1 UNION u2) at
    # one dose, which must NOT equal the two-dose construction -- the whole reason 248 exists.
    single = new.union_residual([u1, u2])
    out_single, _, _ = old.apply_frozen_cd4_residual(
        whitened, mk["old"]["cd4"], columns, panel_positions, extra=single, extra_audit={})
    checks["ON: the one-union-one-dose form (229's ARM3 shape) is a DIFFERENT object"] = not np.allclose(
        out_single, out_on, rtol=0, atol=1e-9)

    # ---------------------------------------------------------------- refusals
    try:
        new.apply_frozen_cd4_residual(whitened, cd4, columns, panel_positions,
                                      extra=u1, extra_audit={}, extra2=u2, extra2_audit={})
        checks["REFUSE: extra2 without audits_out is refused"] = False
    except ValueError:
        checks["REFUSE: extra2 without audits_out is refused"] = True
    new.UNION2_G = None
    try:
        new.apply_frozen_cd4_residual(whitened, cd4, columns, panel_positions,
                                      extra=u1, extra_audit={}, extra2=u2, extra2_audit={}, audits_out={})
        checks["REFUSE: extra2 with UNION2_G unset is refused"] = False
    except ValueError:
        checks["REFUSE: extra2 with UNION2_G unset is refused"] = True
    bad = new.EXP062.CD4Prepared(mode="bad", mean_effect=u2.mean_effect.copy(), gene_votes=u2.gene_votes.copy(),
                                 target_condition_count=u2.target_condition_count.copy(), condition_details=[])
    bad.gene_votes[38, 0] = 1
    new.UNION2_G = 0.35
    try:
        new.apply_frozen_cd4_residual(whitened, cd4, columns, panel_positions,
                                      extra=u1, extra_audit={}, extra2=bad, extra2_audit={}, audits_out={})
        checks["REFUSE: a panel coordinate in union 2's support is refused"] = False
    except AssertionError as exc:
        checks["REFUSE: a panel coordinate in union 2's support is refused"] = "union2" in str(exc)

    (WORK / "unit248.json").write_text(json.dumps(
        {"checks": checks, "frac_cd4": frac_cd4, "f1": f1, "f2": f2, "shared_scalar": s,
         "union2_frame": u2f}, indent=1, default=str))
    for kk, v in checks.items():
        print(f"[{'PASS' if v else 'FAIL'}] {kk}")
    ok = all(checks.values())
    print(f"\nfrac_cd4 {frac_cd4!r}  f1 {f1!r}  f2 {f2!r}  preserve ratio {s!r}")
    print(f"UNIT {'PASS' if ok else 'FAIL'}  {sum(checks.values())}/{len(checks)}  peak {new.peak_rss_gib():.2f} GiB")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
