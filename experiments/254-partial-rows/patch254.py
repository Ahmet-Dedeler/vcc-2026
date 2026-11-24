"""254 — the ONE builder change, applied as exact-anchor replacements so the diff is auditable.

Adds `VCC_CD4_PARTIAL_G` (unset = byte-inert) to `scripts/make_submission_pm_cd4.py`:
  1. the env flag block (after 251's `VCC_CD4_DEFLATE_BASIS` validation);
  2. `extract_cd4_partial()` -- 246 `lib246.cd4_partial(targets, equalise=False)`'s arithmetic
     verbatim on the pre-freeze `prepare_cd4_mean` object -- and `cd4_partial_support_audit()`;
  3. `prepare_frozen_cd4(..., return_partial=True)` hands the partial residual back next to the
     frozen complete-row object (the default call is untouched);
  4. `apply_frozen_cd4_residual(..., partial=, partial_audit=)`: the partial residual is the
     LAST entry of the accumulation order (union 1, union 2, then the partial rows -- 246's
     `merged_many` order for the E arm), sized by `solve_carrying_row_fraction` at
     `CD4_PARTIAL_G` against the same base, ONE add + ONE preserve; its block comes back
     through `audits_out["cd4_partial"]`;
  5. `build_production_effects(..., cd4_partial=, cd4_partial_audit=)` writes `effects.cd4_partial`
     only when on; `run()` prepares it only when on and threads it through BOTH effect builds
     (the main one and `VCC_ALPHA=knee`'s per-context rebuild);
  6. `config.cd4_partial_g`, written only when on.
Every anchor must occur exactly once; the script refuses to run twice.

    uv run python experiments/254-partial-rows/patch254.py            # the working tree
    uv run python experiments/254-partial-rows/patch254.py <copy>     # any byte copy of HEAD
    uv run python experiments/254-partial-rows/patch254.py --check <copy-of-HEAD>
        # apply to the copy in memory and assert it equals the working tree byte for byte
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
WT = ROOT / "scripts/make_submission_pm_cd4.py"


def replace_once(text: str, old: str, new: str) -> str:
    n = text.count(old)
    assert n == 1, f"anchor found {n} times:\n{old[:300]}"
    return text.replace(old, new)


def patch(src: str) -> str:
    assert "VCC_CD4_PARTIAL_G" not in src, "already patched"

    # ------------------------------------------------------------------ 1. the env flag
    OLD1 = '''            "there, with its path and shas, first."
        )
# --- THE SOURCE REGISTRY (experiment 229).  Rule 11: a path AND a sha, never a description.
'''
    NEW1 = '''            "there, with its path and shas, first."
        )
# --- `VCC_CD4_PARTIAL_G` -- CD4's PARTIAL-condition rows as their OWN residual (246 section 3.7 / 254).
#
# THE ROWS.  `prepare_cd4_mean` averages a target over whichever of the three frozen CD4
# conditions carry it; `freeze_complete_mean3` (production's `complete mean3 only` rule) then
# ZEROES every row carried by one or two conditions and keeps only its count.  On A/B/C that is
# 15 of 300 targets (4 one-condition, 11 two-condition) -- CD4 rest/stim effects on 5% of the
# panel that no residual reaches, discarded.  246 section 3.7 measured them as a third residual
# on 162's PLOLO instrument (three cancer folds, carrying-row sizing at g = 0.34 on top of the
# shipped stack): macro cosine +0.587%, `pds` +0.00178 (3/3 folds each), target-shuffled
# control ~0, the row-matched self control +0.57% / -0.0001 (the cosine is "more residual on
# rows that had none", the `pds` is CD4-specific), leave-top-3-out retention -6.6%.  NOT
# build-worthy alone (price +0.0009 ... +0.0018); composed with the k = 5 deflation
# (`E_cd4_svd5+partial@0.34`) the pair reads `pds` +0.01405 (3/3), macro -0.699%, leave-top-3-out
# 61.5%, price +0.0032 / +0.0047 / +0.0100 -- additive to within 0.0002 of the parts.
#
# THE CONSTRUCTION (`experiments/246-direction-hypotheses/lib246.py` `cd4_partial(targets,
# equalise=False)`, reproduced verbatim in `extract_cd4_partial` below): from the SAME
# `prepare_cd4_mean` object the complete-row residual is frozen from, keep exactly the rows with
# `target_condition_count` in {1, 2} -- their effect is already the mean over the conditions
# that carry them, `gene_norm`, at the banks' own rescale (the per-condition multiplier is a
# property of the whole bank, so a partial row carries the same rescale as a complete one);
# zero every other row and its votes.  No re-normalisation (246's `equalise=True` variant read
# LESS `pds`).  Sized by the carrying-row rule (`solve_carrying_row_fraction`, 218 section 7.1
# -- the all-rows median rule would solve exactly 0.0 for a 15-row residual) at this dose
# against the same whitened base, and accumulated LAST into the one merged add (union 1, union 2,
# then the partial rows -- 246's `merged_many` order), then production's ONE global preserve.
#
# THE FLAG.  Unset (the default) is BYTE-INERT: `prepare_frozen_cd4` is called exactly as
# before, `apply_frozen_cd4_residual` accumulates exactly the residuals it accumulated before, and
# every prior build reproduces bit-for-bit (254's known-answer rows: the staged effects array
# `abeb0749...`, 249's arm `a28d5be2...` with `VCC_CD4_DEFLATE_K=5`, and the shipped package
# `8df4dadd...`).  Set to g in (0, 5) and the partial residual is prepared, sized and merged; the
# complete-row CD4 object, its deflation (249/251) and its solved fraction never see it -- the
# only CD4-block scalars it can move are the merged preserve pair.  `effects.cd4_partial` and
# `config.cd4_partial_g` are written only when on (the `h1_union_g` convention).
CD4_PARTIAL_G = os.environ.get("VCC_CD4_PARTIAL_G", "").strip()
if CD4_PARTIAL_G:
    try:
        CD4_PARTIAL_G = float(CD4_PARTIAL_G)
    except ValueError as exc:
        raise ValueError(
            f"VCC_CD4_PARTIAL_G must be a number in (0, 5); got {CD4_PARTIAL_G!r}"
        ) from exc
    if not 0.0 < CD4_PARTIAL_G < 5.0:
        raise ValueError(f"VCC_CD4_PARTIAL_G must be in (0, 5); got {CD4_PARTIAL_G}")
else:
    CD4_PARTIAL_G = None
# --- THE SOURCE REGISTRY (experiment 229).  Rule 11: a path AND a sha, never a description.
'''
    src = replace_once(src, OLD1, NEW1)

    # ------------------------------------------------------------------ 2/3. the functions
    OLD2 = '''def prepare_frozen_cd4(
    root: Path,
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
):
    banks = tuple(
        EXP062.Bank.load(name, root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    prepared = EXP062.prepare_cd4_mean(
        banks,
        targets,
        official_genes,
        set(map(str, panel)),
        CD4_MODE,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    return banks, freeze_complete_mean3(prepared)
'''
    NEW2 = '''def extract_cd4_partial(prepared):
    """246's `lib246.cd4_partial(targets, equalise=False)`, verbatim: the rows the freeze zeroes.

    From the pre-freeze `prepare_cd4_mean` object: rows with `target_condition_count` in
    {1, 2} keep their effect (already the mean over the conditions that carry them, `gene_norm`,
    the banks' own rescale) and their votes; every other row is zeroed.  The counts of the
    kept rows are kept (1 or 2), all others 0, so `target_condition_count > 0` is exactly
    "this row is a partial row".  `condition_details[0]` records the row census and the median
    norm of the complete rows (246's `median_complete_norm`), for the meta.
    """
    count = np.asarray(prepared.target_condition_count, dtype=np.int64)
    complete = count == len(CD4_NAMES)
    partial = (count > 0) & ~complete
    eff = np.asarray(prepared.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(prepared.gene_votes, dtype=np.uint8).copy()
    norms_c = np.linalg.norm(np.where(votes[complete] > 0, eff[complete], 0.0), axis=1)
    eff[~partial] = 0.0
    votes[~partial] = 0
    cnt = np.where(partial, count, 0).astype(np.uint8)
    return EXP062.CD4Prepared(
        mode="cd4_partial",
        mean_effect=eff,
        gene_votes=votes,
        target_condition_count=cnt,
        condition_details=[
            dict(
                n_partial=int(partial.sum()),
                n_one=int((count == 1).sum()),
                n_two=int((count == 2).sum()),
                median_complete_norm=float(np.median(norms_c)) if norms_c.size else 0.0,
            )
        ],
    )


def cd4_partial_support_audit(partial, *, target_g: float) -> dict[str, object]:
    """The partial residual's meta block (`effects.cd4_partial`), before sizing."""
    counts = np.asarray(partial.target_condition_count, dtype=np.uint8)
    votes = np.asarray(partial.gene_votes)
    info = dict(partial.condition_details[0]) if partial.condition_details else {}
    return {
        "mode": str(partial.mode),
        "conditions": list(CD4_NAMES),
        "rule": (
            "rows with target_condition_count in {1, 2} of the SAME prepare_cd4_mean object "
            "the complete-row residual is frozen from (production's `complete mean3 only` "
            "zeroes them); effect = mean over the conditions that carry the row, gene_norm, "
            "the banks' own rescale; no re-normalisation (246 section 3.7, "
            "lib246.cd4_partial equalise=False)"
        ),
        "rows": int(info.get("n_partial", int(np.sum(counts > 0)))),
        "rows_one_condition": int(info.get("n_one", int(np.sum(counts == 1)))),
        "rows_two_conditions": int(info.get("n_two", int(np.sum(counts == 2)))),
        "carrying_rows": int(np.sum(votes.sum(axis=1) > 0)),
        "target_condition_count_histogram": {
            str(value): int(np.sum(counts == value))
            for value in sorted(set(counts.tolist()))
        },
        "supported_entries": int(np.count_nonzero(votes)),
        "median_complete_row_norm": float(info.get("median_complete_norm", 0.0)),
        "sizing_rule": "median over CARRYING rows (coverage-aware)",
        "target_g": float(target_g),
    }


def prepare_frozen_cd4(
    root: Path,
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel: Sequence[str],
    *,
    chunk_rows: int,
    max_rss_gib: float,
    return_partial: bool = False,
):
    """`(banks, frozen)`; with `return_partial=True` (254) also the partial-row residual
    extracted from the SAME prepared object before the freeze: `(banks, frozen, partial)`.
    The default call performs exactly the operations it always performed."""
    banks = tuple(
        EXP062.Bank.load(name, root / name, require_cd4_scale=True)
        for name in CD4_NAMES
    )
    prepared = EXP062.prepare_cd4_mean(
        banks,
        targets,
        official_genes,
        set(map(str, panel)),
        CD4_MODE,
        chunk_rows=chunk_rows,
        max_rss_gib=max_rss_gib,
    )
    if not return_partial:
        return banks, freeze_complete_mean3(prepared)
    return banks, freeze_complete_mean3(prepared), extract_cd4_partial(prepared)
'''
    src = replace_once(src, OLD2, NEW2)

    # ------------------------------------------------------------------ 4. the merge
    OLD4a = '''    extra2=None,
    extra2_audit: dict[str, object] | None = None,
    audits_out: dict[str, object] | None = None,
) -> tuple[np.ndarray, dict[str, object], dict[str, object] | None]:
    """Apply frozen post-whitening CD4 and enforce zero panel support.
'''
    NEW4a = '''    extra2=None,
    extra2_audit: dict[str, object] | None = None,
    audits_out: dict[str, object] | None = None,
    partial=None,
    partial_audit: dict[str, object] | None = None,
) -> tuple[np.ndarray, dict[str, object], dict[str, object] | None]:
    """Apply frozen post-whitening CD4 and enforce zero panel support.

    `partial` is the optional residual of `VCC_CD4_PARTIAL_G` (254): CD4's partial-condition
    rows, sized by the same carrying-row rule at `CD4_PARTIAL_G` against the same base and
    accumulated LAST (after union 1 and union 2 -- 246's `merged_many` order), still ONE add
    and ONE preserve.  Its block comes back through `audits_out["cd4_partial"]`.  With
    `partial` None nothing here changes.
'''
    src = replace_once(src, OLD4a, NEW4a)

    OLD4b = '''    residuals = [
        (tag, res, g)
        for tag, res, g in (("h1_union", extra, H1_UNION_G), ("union2", extra2, UNION2_G))
        if res is not None
    ]
'''
    NEW4b = '''    residuals = [
        (tag, res, g)
        for tag, res, g in (
            ("h1_union", extra, H1_UNION_G),
            ("union2", extra2, UNION2_G),
            # 254: CD4's partial-condition rows, LAST in the accumulation order.
            ("cd4_partial", partial, CD4_PARTIAL_G),
        )
        if res is not None
    ]
'''
    src = replace_once(src, OLD4b, NEW4b)

    OLD4c = '''        audits_out["union2"] = union2_audit
    if extra is None:
        return output, cd4_audit, None
'''
    NEW4c = '''        audits_out["union2"] = union2_audit
    if partial is not None:
        # 254: the partial-row residual's block, the `h1_union` fields plus its row census,
        # handed back through the same out-dict.
        cd4_partial_written = {
            **(partial_audit or {}),
            "residual_frame": frames["cd4_partial"],
            "merged_add": True,
            "merged_preserve": preserve,
            "panel_residual_nonzero_entries": int(
                np.count_nonzero(np.asarray(partial.mean_effect)[:, panel_positions])
            ),
        }
        if audits_out is None:
            raise ValueError(
                "partial was given but no `audits_out` dict to return its audit through"
            )
        audits_out["cd4_partial"] = cd4_partial_written
    if extra is None:
        return output, cd4_audit, None
'''
    src = replace_once(src, OLD4c, NEW4c)

    # ------------------------------------------------------------------ 5. build_production_effects
    OLD5a = '''    union2=None,
    union2_audit: dict[str, object] | None = None,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Build one context-independent production effect matrix.

    `union2` / `union2_audit` are the second union residual and its support audit (248);
    None (the default) is the shipped path."""
'''
    NEW5a = '''    union2=None,
    union2_audit: dict[str, object] | None = None,
    cd4_partial=None,
    cd4_partial_audit: dict[str, object] | None = None,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Build one context-independent production effect matrix.

    `union2` / `union2_audit` are the second union residual and its support audit (248);
    `cd4_partial` / `cd4_partial_audit` the partial-condition CD4 residual and its audit (254);
    None (the default) is the shipped path."""
'''
    src = replace_once(src, OLD5a, NEW5a)

    OLD5b = '''            extra2=union2,
            extra2_audit=union2_audit,
            audits_out=extra_audits,
        )
        cd4_audit = {"enabled": True, **cd4_audit}
    union2_written = extra_audits.get("union2")
'''
    NEW5b = '''            extra2=union2,
            extra2_audit=union2_audit,
            audits_out=extra_audits,
            partial=cd4_partial,
            partial_audit=cd4_partial_audit,
        )
        cd4_audit = {"enabled": True, **cd4_audit}
    union2_written = extra_audits.get("union2")
    cd4_partial_written = extra_audits.get("cd4_partial")
'''
    src = replace_once(src, OLD5b, NEW5b)

    OLD5c = '''        # Emitted only when the second union is ON (248), the same convention.
        **({} if union2_written is None else {"union2": union2_written}),
    }
'''
    NEW5c = '''        # Emitted only when the second union is ON (248), the same convention.
        **({} if union2_written is None else {"union2": union2_written}),
        # Emitted only when the partial-row residual is ON (254), the same convention.
        **({} if cd4_partial_written is None else {"cd4_partial": cd4_partial_written}),
    }
'''
    src = replace_once(src, OLD5c, NEW5c)

    # ------------------------------------------------------------------ 5. run()
    OLD6a = '''    cd4_banks, cd4 = prepare_frozen_cd4(
        cd4_root,
        tuple(map(str, spec.perts)),
        tuple(map(str, spec.genes)),
        tuple(map(str, spec.perts)),
        chunk_rows=args.chunk_rows,
        max_rss_gib=args.max_rss_gib,
    )
    if H1_UNION_G is None:
'''
    NEW6a = '''    if CD4_PARTIAL_G is None:
        cd4_banks, cd4 = prepare_frozen_cd4(
            cd4_root,
            tuple(map(str, spec.perts)),
            tuple(map(str, spec.genes)),
            tuple(map(str, spec.perts)),
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
        )
        cd4_partial = None
        cd4_partial_audit = None
    else:
        # 254: the partial-condition rows, extracted from the SAME prepared object the
        # complete-row residual is frozen from; sized and merged inside
        # `apply_frozen_cd4_residual`.
        cd4_banks, cd4, cd4_partial = prepare_frozen_cd4(
            cd4_root,
            tuple(map(str, spec.perts)),
            tuple(map(str, spec.genes)),
            tuple(map(str, spec.perts)),
            chunk_rows=args.chunk_rows,
            max_rss_gib=args.max_rss_gib,
            return_partial=True,
        )
        cd4_partial_audit = cd4_partial_support_audit(cd4_partial, target_g=CD4_PARTIAL_G)
        print(
            f"[cd4-partial] g={CD4_PARTIAL_G:g}  {cd4_partial_audit['rows']} partial-condition "
            f"rows of {len(spec.perts)} ({cd4_partial_audit['rows_one_condition']} one-condition, "
            f"{cd4_partial_audit['rows_two_conditions']} two-condition), "
            f"{cd4_partial_audit['carrying_rows']} carrying, "
            f"{cd4_partial_audit['supported_entries']} supported entries; the complete-row "
            f"object is untouched",
            flush=True,
        )
        memory_checkpoint("cd4 partial residual", args.max_rss_gib)
    if H1_UNION_G is None:
'''
    src = replace_once(src, OLD6a, NEW6a)

    OLD6b = '''        union2=union2,
        union2_audit=union2_audit,
        weights=weights,
    )
    rank_scores = None
'''
    NEW6b = '''        union2=union2,
        union2_audit=union2_audit,
        cd4_partial=cd4_partial,
        cd4_partial_audit=cd4_partial_audit,
        weights=weights,
    )
    rank_scores = None
'''
    src = replace_once(src, OLD6b, NEW6b)

    OLD6c = '''                union2=union2,
                union2_audit=union2_audit,
                alpha=float(solved[context]["alpha"]),
'''
    NEW6c = '''                union2=union2,
                union2_audit=union2_audit,
                cd4_partial=cd4_partial,
                cd4_partial_audit=cd4_partial_audit,
                alpha=float(solved[context]["alpha"]),
'''
    src = replace_once(src, OLD6c, NEW6c)

    # ------------------------------------------------------------------ 6. the config key
    OLD7 = '''            # Written only when the deflation uses a PINNED basis (251), the same convention.
            **({} if CD4_DEFLATE_BASIS is None else {"cd4_deflate_basis": CD4_DEFLATE_BASIS}),
'''
    NEW7 = OLD7 + '''            # Written only when CD4's partial-condition rows ride as their own residual (254).
            **({} if CD4_PARTIAL_G is None else {"cd4_partial_g": CD4_PARTIAL_G}),
'''
    src = replace_once(src, OLD7, NEW7)
    return src


def main() -> None:
    if len(sys.argv) >= 3 and sys.argv[1] == "--check":
        copy = Path(sys.argv[2]).resolve()
        patched = patch(copy.read_text())
        same = patched == WT.read_text()
        print(f"--check: HEAD copy {copy} + patch254 == working tree: {same}")
        raise SystemExit(0 if same else 1)
    path = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else WT
    path.write_text(patch(path.read_text()))
    print(f"patched {path}")


if __name__ == "__main__":
    main()
