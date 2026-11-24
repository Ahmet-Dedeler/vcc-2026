"""248 -- the ONE builder change, as exact-anchor replacements (249's `patch249.py` pattern).

Adds a SECOND union residual to `scripts/make_submission_pm_cd4.py`:

    VCC_UNION2_SOURCES   ordered comma-separated registry names (unset = inert)
    VCC_UNION2_G         its carrying-row dose                  (unset = inert)

When set, CD4 (its own rule), union 1 (`VCC_UNION_SOURCES` / `VCC_H1_UNION_G`, carrying-row
rule) and union 2 (carrying-row rule) are each solved against the SAME whitened base, merged
into ONE residual, added in ONE step, followed by production's ONE global preserve -- 218
section 7.2's `merged()` chained over three terms, which is 225/229's `CUM` construction.  The
meta records `effects.union2` with the same fields as `effects.h1_union`, and the config gains
`union2_g`, `union2_sources`, `union2_sizing_rule`, written only when on.

Every replacement is anchored on a unique exact string; the script refuses to run if an
anchor is missing or repeated.  Applied to the working tree it must reproduce
`work/make_submission_pm_cd4.248.py` (the intended builder) byte for byte.

    uv run python experiments/248-second-residual/patch248.py [--check <path>] [--apply]
"""
from __future__ import annotations

import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
TARGET = ROOT / "scripts/make_submission_pm_cd4.py"

REPLACEMENTS: list[tuple[str, str, str]] = []


def R(label: str, old: str, new: str) -> None:
    REPLACEMENTS.append((label, old, new))


# ---------------------------------------------------------------------------------------
# 1. the env block, right after union 1's file list
R("env block",
  '''H1_UNION_FILES = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt")
MIN_CPM = 5.0
''',
  '''H1_UNION_FILES = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt")
# --- `VCC_UNION2_SOURCES` + `VCC_UNION2_G` -- a SECOND union residual (experiment 248).
#
# WHY.  225 priced the best arm on its table as a second residual RIDING ON TOP of the shipped
# `h1union`: `cd4@0.34 + h1pair@0.28`, then `+ dld1_lowmoi@0.35` (`CUM1`) or `+ 225's four
# @0.35` (`CUM4`), each union sized on ITS OWN carrying rows.  229 found that construction is
# not expressible with one union at one dose -- naming a source in `VCC_UNION_SOURCES` either
# replaces the H1 pair or dilutes DLD-1 into a six-source union at one shared dose (`ARM6`), and
# both read lower than the two-residual form.  This flag expresses it exactly.
#
# THE CONSTRUCTION (218 section 7.2's `merged()`, chained over THREE terms).  CD4 is sized by
# its own rule (`solve_cd4_dose`, untouched), union 1 by `solve_carrying_row_fraction` at
# `VCC_H1_UNION_G`, union 2 by `solve_carrying_row_fraction` at `VCC_UNION2_G` -- all three
# against the SAME whitened base, so every dose means the same thing -- then
#     merged = frac_cd4 * cd4 + frac_1 * union1 + frac_2 * union2
# is added in ONE step and production's ONE global preserve rescale follows.  The accumulation
# order is CD4, union 1, union 2, and each term is added with the same operations the union-1
# path already used, so with union 2 unset the floating-point path is the shipped one bit for
# bit (248 known-answer row 1: the staged effects `abeb0749...` and the package `8df4dadd...`),
# and union 2 alone with union 1 unset is the same arithmetic as union 1 alone (row 2: the H1
# pair named here with `VCC_H1_UNION_G` unset reproduces `abeb0749...` bit for bit).
#
# THE RULES.  Both variables or neither -- a dose without sources, or sources without a dose,
# is a knob that is set and cannot act (242 #17) and is refused at import.  Names must be in
# `UNION_SOURCE_REGISTRY`, may not repeat, and may not also be in union 1 while union 1 is on:
# a source in both unions would be dosed twice.  Union 2's banks are hash-verified before the
# build exactly as union 1's are (`verify_union_source_hashes`).  ORDER IS PART OF THE ARM
# (`union_residual` accumulates a float64 sum), and it is recorded in `config.union2_sources`.
#
# THE META.  `effects.union2` carries the same fields as `effects.h1_union` (sources with
# their four-file shas, carrying rows, `g_per_unit_fraction`, `solved_fraction`,
# `shipped_median_rule_would_solve`, the merged preserve) plus two overlap counts against
# union 1; `config.union2_g`, `config.union2_sources`, `config.union2_sizing_rule` are written
# only when on (the `h1_union_g` convention), so a default build's config block and audit are
# byte-identical to every prior build's.
UNION2_SOURCES_SELECTED = os.environ.get("VCC_UNION2_SOURCES", "").strip()
UNION2_G_ENV = os.environ.get("VCC_UNION2_G", "").strip()
UNION2_SOURCE_NAMES: tuple[str, ...] = tuple(
    part.strip() for part in UNION2_SOURCES_SELECTED.split(",") if part.strip()
)
if UNION2_SOURCES_SELECTED and not UNION2_SOURCE_NAMES:
    raise ValueError(f"VCC_UNION2_SOURCES names no source: {UNION2_SOURCES_SELECTED!r}")
if UNION2_SOURCE_NAMES and not UNION2_G_ENV:
    raise ValueError(
        "VCC_UNION2_SOURCES is set but VCC_UNION2_G is not, so the second union would have "
        "no dose; set VCC_UNION2_G or unset VCC_UNION2_SOURCES"
    )
if UNION2_G_ENV and not UNION2_SOURCE_NAMES:
    # A knob that is set and cannot act is the 242 #17 failure class; refuse it.
    raise ValueError(
        "VCC_UNION2_G is set but VCC_UNION2_SOURCES is not, so it would do nothing; "
        "unset it or name the second union's sources"
    )
UNION2_G = float(UNION2_G_ENV) if UNION2_G_ENV else None
if UNION2_G is not None and not 0.0 < UNION2_G < 5.0:
    raise ValueError(f"VCC_UNION2_G must be in (0, 5); got {UNION2_G}")
if len(set(UNION2_SOURCE_NAMES)) != len(UNION2_SOURCE_NAMES):
    raise ValueError(f"VCC_UNION2_SOURCES repeats a source: {UNION2_SOURCE_NAMES}")
for _name in UNION2_SOURCE_NAMES:
    if _name not in UNION_SOURCE_REGISTRY:
        raise ValueError(
            f"VCC_UNION2_SOURCES names {_name!r}, which is not in UNION_SOURCE_REGISTRY "
            f"({sorted(UNION_SOURCE_REGISTRY)}).  Add it there, with its shas, first."
        )
if UNION2_SOURCE_NAMES and H1_UNION_G is not None:
    _shared = sorted(set(UNION2_SOURCE_NAMES) & set(UNION_SOURCE_NAMES))
    if _shared:
        raise ValueError(
            f"VCC_UNION2_SOURCES shares {_shared} with union 1 ({list(UNION_SOURCE_NAMES)}) "
            "while VCC_H1_UNION_G is set: a source in both unions would be dosed twice"
        )
UNION2_SOURCES: tuple[tuple[str, Path], ...] | None = (
    tuple((name, UNION_SOURCE_REGISTRY[name]["path"]) for name in UNION2_SOURCE_NAMES)
    if UNION2_SOURCE_NAMES
    else None
)
MIN_CPM = 5.0
''')

# ---------------------------------------------------------------------------------------
# 2. the hash verifier takes the source list (default: union 1's), looked up at CALL time
R("verify_union_source_hashes signature",
  '''def verify_union_source_hashes() -> dict[str, dict[str, str]]:
    """Rule 11, ENFORCED: every selected bank's four files must hash to the registry's value.

    Three of the registered banks live under `experiments/*/work/`, a tree other experiments
    write to.  Recording a sha catches a changed file only if something compares it; this
    compares it, before a single byte of a 35-minute build is spent, and refuses to run
    otherwise.  Returns the measured table so the caller does not hash twice.
    """
    measured: dict[str, dict[str, str]] = {}
    for name, path in H1_UNION_SOURCES:
''',
  '''def verify_union_source_hashes(sources=None) -> dict[str, dict[str, str]]:
    """Rule 11, ENFORCED: every selected bank's four files must hash to the registry's value.

    Three of the registered banks live under `experiments/*/work/`, a tree other experiments
    write to.  Recording a sha catches a changed file only if something compares it; this
    compares it, before a single byte of a 35-minute build is spent, and refuses to run
    otherwise.  Returns the measured table so the caller does not hash twice.

    `sources` is a `(name, path)` list; None means union 1's (`H1_UNION_SOURCES`, read at
    call time so a harness that rebinds it still drives this).  Union 2 passes its own (248).
    """
    if sources is None:
        sources = H1_UNION_SOURCES
    measured: dict[str, dict[str, str]] = {}
    for name, path in sources:
''')

# ---------------------------------------------------------------------------------------
# 3. the union builder takes the source list too
R("prepare_h1_union signature",
  '''def prepare_h1_union(
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel: Sequence[str],
    *,
    max_rss_gib: float,
):
    """The union of the selected registry banks on the official axis, panel columns excluded.

    `chunk_rows` is pinned at 218's 256 rather than taken from `--chunk-rows`: the streamed
    column mean accumulates chunk by chunk, so the chunk size is a (16th-decimal) numerical
    choice and this residual has to reproduce the object 218 measured.
    """
    verify_union_source_hashes()
    panel_targets = set(map(str, panel))
    prepared = []
    details = []
    for name, path in H1_UNION_SOURCES:
''',
  '''def prepare_h1_union(
    targets: Sequence[str],
    official_genes: Sequence[str],
    panel: Sequence[str],
    *,
    max_rss_gib: float,
    sources=None,
):
    """The union of the selected registry banks on the official axis, panel columns excluded.

    `chunk_rows` is pinned at 218's 256 rather than taken from `--chunk-rows`: the streamed
    column mean accumulates chunk by chunk, so the chunk size is a (16th-decimal) numerical
    choice and this residual has to reproduce the object 218 measured.

    `sources` is a `(name, path)` list; None means union 1's (`H1_UNION_SOURCES`).  The second
    union (248) passes `UNION2_SOURCES`; the construction is identical.
    """
    if sources is None:
        sources = H1_UNION_SOURCES
    verify_union_source_hashes(sources)
    panel_targets = set(map(str, panel))
    prepared = []
    details = []
    for name, path in sources:
''')

# ---------------------------------------------------------------------------------------
# 4. the support audit takes the source list and the dose
R("h1_union_support_audit",
  '''def h1_union_support_audit(union, details) -> dict[str, object]:
    counts = np.asarray(union.target_condition_count, dtype=np.uint8)
    return {
''',
  '''def h1_union_support_audit(
    union, details, *, sources=None, target_g=None
) -> dict[str, object]:
    """The union's meta block.  `sources` / `target_g` default to union 1's (`H1_UNION_SOURCES`
    / `H1_UNION_G`, read at call time); union 2 (248) passes `UNION2_SOURCES` / `UNION2_G` and
    gets a block with exactly the same fields."""
    if sources is None:
        sources = H1_UNION_SOURCES
    if target_g is None:
        target_g = H1_UNION_G
    counts = np.asarray(union.target_condition_count, dtype=np.uint8)
    return {
''')
R("h1_union_support_audit sources loop",
  '''            for name, path in H1_UNION_SOURCES
        ],
        "source_details": details,
        "mode": str(union.mode),
        "sizing_rule": "median over CARRYING rows (coverage-aware)",
        "target_g": float(H1_UNION_G),
''',
  '''            for name, path in sources
        ],
        "source_details": details,
        "mode": str(union.mode),
        "sizing_rule": "median over CARRYING rows (coverage-aware)",
        "target_g": float(target_g),
''')
R("h1_union_support_audit source count",
  '''            str(k): int(np.sum(counts == k)) for k in range(0, len(H1_UNION_SOURCES) + 1)
''',
  '''            str(k): int(np.sum(counts == k)) for k in range(0, len(sources) + 1)
''')

# ---------------------------------------------------------------------------------------
# 5. the merge: three terms, one add, one preserve
R("apply_frozen_cd4_residual",
  '''def apply_frozen_cd4_residual(
    whitened: np.ndarray,
    cd4,
    operator_columns: np.ndarray,
    panel_positions: np.ndarray,
    *,
    extra=None,
    extra_audit: dict[str, object] | None = None,
) -> tuple[np.ndarray, dict[str, object], dict[str, object] | None]:
    """Apply frozen post-whitening CD4 and enforce zero panel support.

    `extra` is the optional second residual of `VCC_H1_UNION_G`.  When it is given, CD4 and the
    second source are solved against the SAME denominator (each by its own sizing rule), merged
    into one residual carrying `frac_cd4 * cd4 + frac_extra * extra`, and added in ONE step
    followed by production's single global preserve rescale -- 218 section 7.2's `merged()`.
    CD4's own solved fraction is bit-identical either way; the only thing the second residual
    can move in CD4's own audit is the shared preserve scalar, which is what "one global
    rescale" means.
    """
    panel_positions = np.asarray(panel_positions, dtype=np.int64)
''',
  '''def apply_frozen_cd4_residual(
    whitened: np.ndarray,
    cd4,
    operator_columns: np.ndarray,
    panel_positions: np.ndarray,
    *,
    extra=None,
    extra_audit: dict[str, object] | None = None,
    extra2=None,
    extra2_audit: dict[str, object] | None = None,
    audits_out: dict[str, object] | None = None,
) -> tuple[np.ndarray, dict[str, object], dict[str, object] | None]:
    """Apply frozen post-whitening CD4 and enforce zero panel support.

    `extra` is the optional second residual of `VCC_H1_UNION_G`.  When it is given, CD4 and the
    second source are solved against the SAME denominator (each by its own sizing rule), merged
    into one residual carrying `frac_cd4 * cd4 + frac_extra * extra`, and added in ONE step
    followed by production's single global preserve rescale -- 218 section 7.2's `merged()`.
    CD4's own solved fraction is bit-identical either way; the only thing the second residual
    can move in CD4's own audit is the shared preserve scalar, which is what "one global
    rescale" means.

    `extra2` is the optional THIRD residual of `VCC_UNION2_SOURCES` / `VCC_UNION2_G` (248),
    sized by the same carrying-row rule at `UNION2_G` against the same base and accumulated
    after `extra`, so the merged residual is `frac_cd4 * cd4 + frac_1 * extra + frac_2 *
    extra2` -- 218's `merged()` chained -- still ONE add and ONE preserve.  Its audit block is
    handed back through `audits_out["union2"]` (the return arity is unchanged so every existing
    caller is untouched).  With `extra2` None nothing here changes: the union-1 path performs
    exactly the operations it performed before, in the same order.
    """
    panel_positions = np.asarray(panel_positions, dtype=np.int64)
''')
R("apply_frozen_cd4_residual merge body",
  '''    fraction, frame = solve_cd4_dose(whitened, cd4)
    if extra is None:
        added = cd4
        applied = fraction
        extra_frame = None
        support = np.asarray(cd4.gene_votes) > 0
    else:
        if np.any(np.asarray(extra.gene_votes)[:, panel_positions] != 0):
            raise AssertionError("a VCC panel coordinate entered the second residual support")
        if np.any(np.asarray(extra.mean_effect)[:, panel_positions] != 0):
            raise AssertionError("a VCC panel coordinate entered the second residual")
        extra_fraction, extra_frame = solve_carrying_row_fraction(
            whitened, extra, float(H1_UNION_G)
        )
        cd4_votes = np.asarray(cd4.gene_votes) > 0
        extra_votes = np.asarray(extra.gene_votes) > 0
        merged_effect = np.where(cd4_votes, np.asarray(cd4.mean_effect), 0.0) * fraction
        merged_effect += (
            np.where(extra_votes, np.asarray(extra.mean_effect), 0.0) * extra_fraction
        )
        support = cd4_votes | extra_votes
        added = EXP062.CD4Prepared(
            mode="merged",
            mean_effect=merged_effect,
            gene_votes=support.astype(np.uint8),
            target_condition_count=np.maximum(
                np.asarray(cd4.target_condition_count),
                np.asarray(extra.target_condition_count),
            ).astype(np.uint8),
            condition_details=[],
        )
        applied = 1.0
        extra_frame["applied_fraction"] = float(extra_fraction)
        extra_frame["entries_outside_cd4_support"] = int(
            np.count_nonzero(extra_votes & ~cd4_votes)
        )
        extra_frame["targets_outside_cd4_support"] = int(
            np.sum(
                (np.asarray(extra.target_condition_count) > 0)
                & (np.asarray(cd4.target_condition_count) == 0)
            )
        )

    output, preserve = EXP062.add_post_whiten_residual(
''',
  '''    fraction, frame = solve_cd4_dose(whitened, cd4)
    # The extra residuals, in accumulation order: union 1 (`VCC_H1_UNION_G`), then union 2
    # (`VCC_UNION2_G`, 248).  Each present one is sized by the carrying-row rule at its own
    # dose against the same whitened base and added to the merged residual with the same
    # operations the union-1 path always used, so an absent union 2 changes nothing.
    residuals = [
        (tag, res, g)
        for tag, res, g in (("h1_union", extra, H1_UNION_G), ("union2", extra2, UNION2_G))
        if res is not None
    ]
    frames: dict[str, dict[str, object]] = {}
    if not residuals:
        added = cd4
        applied = fraction
        support = np.asarray(cd4.gene_votes) > 0
    else:
        cd4_votes = np.asarray(cd4.gene_votes) > 0
        merged_effect = np.where(cd4_votes, np.asarray(cd4.mean_effect), 0.0) * fraction
        support = cd4_votes
        counts = np.asarray(cd4.target_condition_count)
        for tag, res, g in residuals:
            if np.any(np.asarray(res.gene_votes)[:, panel_positions] != 0):
                raise AssertionError(
                    f"a VCC panel coordinate entered the {tag} residual support"
                )
            if np.any(np.asarray(res.mean_effect)[:, panel_positions] != 0):
                raise AssertionError(f"a VCC panel coordinate entered the {tag} residual")
            if g is None:
                raise ValueError(f"the {tag} residual was given without a dose")
            res_fraction, res_frame = solve_carrying_row_fraction(whitened, res, float(g))
            res_votes = np.asarray(res.gene_votes) > 0
            merged_effect += (
                np.where(res_votes, np.asarray(res.mean_effect), 0.0) * res_fraction
            )
            support = support | res_votes
            counts = np.maximum(counts, np.asarray(res.target_condition_count))
            res_frame["applied_fraction"] = float(res_fraction)
            res_frame["entries_outside_cd4_support"] = int(
                np.count_nonzero(res_votes & ~cd4_votes)
            )
            res_frame["targets_outside_cd4_support"] = int(
                np.sum(
                    (np.asarray(res.target_condition_count) > 0)
                    & (np.asarray(cd4.target_condition_count) == 0)
                )
            )
            if tag == "union2" and extra is not None:
                # Where the two unions meet: rows both carry receive both residuals at their
                # own doses (225's CUM frame), so the overlap is recorded, never hidden.
                res_frame["targets_shared_with_h1_union"] = int(
                    np.sum(
                        (np.asarray(res.target_condition_count) > 0)
                        & (np.asarray(extra.target_condition_count) > 0)
                    )
                )
                res_frame["entries_shared_with_h1_union"] = int(
                    np.count_nonzero(res_votes & (np.asarray(extra.gene_votes) > 0))
                )
            frames[tag] = res_frame
        added = EXP062.CD4Prepared(
            mode="merged",
            mean_effect=merged_effect,
            gene_votes=support.astype(np.uint8),
            target_condition_count=counts.astype(np.uint8),
            condition_details=[],
        )
        applied = 1.0

    output, preserve = EXP062.add_post_whiten_residual(
''')
R("apply_frozen_cd4_residual return",
  '''    if extra is None:
        reported_preserve: dict[str, object] = preserve
    else:
        # `add_post_whiten_residual` echoes the fraction it was handed, which is 1.0 for a
        # merged add.  Report CD4's OWN solved fraction here instead, so this key keeps the
        # meaning it has in every prior build's meta and a rule-10 diff cannot read the merge
        # as a change in CD4's dose.  The two norm scalars below it DO move -- one global
        # preserve over one merged add is exactly what 218 measured.
        reported_preserve = {**preserve, "fraction": fraction}
''',
  '''    if not residuals:
        reported_preserve: dict[str, object] = preserve
    else:
        # `add_post_whiten_residual` echoes the fraction it was handed, which is 1.0 for a
        # merged add.  Report CD4's OWN solved fraction here instead, so this key keeps the
        # meaning it has in every prior build's meta and a rule-10 diff cannot read the merge
        # as a change in CD4's dose.  The two norm scalars below it DO move -- one global
        # preserve over one merged add is exactly what 218 measured.
        reported_preserve = {**preserve, "fraction": fraction}
''')
R("apply_frozen_cd4_residual union audits",
  '''    if extra is None:
        return output, cd4_audit, None
    union_audit = {
        **(extra_audit or {}),
        "residual_frame": extra_frame,
        "merged_add": True,
        "merged_preserve": preserve,
        "panel_residual_nonzero_entries": int(
            np.count_nonzero(np.asarray(extra.mean_effect)[:, panel_positions])
        ),
    }
    return output, cd4_audit, union_audit
''',
  '''    if extra2 is not None:
        # 248: the second union's block, same fields as `h1_union`'s, handed back through the
        # optional out-dict so the return arity every existing caller relies on is unchanged.
        union2_audit = {
            **(extra2_audit or {}),
            "residual_frame": frames["union2"],
            "merged_add": True,
            "merged_preserve": preserve,
            "panel_residual_nonzero_entries": int(
                np.count_nonzero(np.asarray(extra2.mean_effect)[:, panel_positions])
            ),
        }
        if audits_out is None:
            raise ValueError(
                "extra2 was given but no `audits_out` dict to return its audit through"
            )
        audits_out["union2"] = union2_audit
    if extra is None:
        return output, cd4_audit, None
    union_audit = {
        **(extra_audit or {}),
        "residual_frame": frames["h1_union"],
        "merged_add": True,
        "merged_preserve": preserve,
        "panel_residual_nonzero_entries": int(
            np.count_nonzero(np.asarray(extra.mean_effect)[:, panel_positions])
        ),
    }
    return output, cd4_audit, union_audit
''')

# ---------------------------------------------------------------------------------------
# 6. build_production_effects threads union 2 through and writes `effects.union2`
R("build_production_effects signature",
  '''    h1_union=None,
    h1_union_audit: dict[str, object] | None = None,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Build one context-independent production effect matrix."""
''',
  '''    h1_union=None,
    h1_union_audit: dict[str, object] | None = None,
    union2=None,
    union2_audit: dict[str, object] | None = None,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Build one context-independent production effect matrix.

    `union2` / `union2_audit` are the second union residual and its support audit (248);
    None (the default) is the shipped path."""
''')
R("build_production_effects call",
  '''    union_audit: dict[str, object] | None = None
    if cd4 is None:
        final = whitened
        cd4_audit: dict[str, object] = {"enabled": False}
    else:
        final, cd4_audit, union_audit = apply_frozen_cd4_residual(
            whitened,
            cd4,
            operator_columns,
            panel_positions,
            extra=h1_union,
            extra_audit=h1_union_audit,
        )
        cd4_audit = {"enabled": True, **cd4_audit}
''',
  '''    union_audit: dict[str, object] | None = None
    extra_audits: dict[str, object] = {}
    if cd4 is None:
        final = whitened
        cd4_audit: dict[str, object] = {"enabled": False}
    else:
        final, cd4_audit, union_audit = apply_frozen_cd4_residual(
            whitened,
            cd4,
            operator_columns,
            panel_positions,
            extra=h1_union,
            extra_audit=h1_union_audit,
            extra2=union2,
            extra2_audit=union2_audit,
            audits_out=extra_audits,
        )
        cd4_audit = {"enabled": True, **cd4_audit}
    union2_written = extra_audits.get("union2")
''')
R("build_production_effects audit dict",
  '''        **({} if union_audit is None else {"h1_union": union_audit}),
    }


def _required_files(directory: Path, names: Sequence[str]) -> list[Path]:
''',
  '''        **({} if union_audit is None else {"h1_union": union_audit}),
        # Emitted only when the second union is ON (248), the same convention.
        **({} if union2_written is None else {"union2": union2_written}),
    }


def _required_files(directory: Path, names: Sequence[str]) -> list[Path]:
''')

# ---------------------------------------------------------------------------------------
# 7. run(): prepare union 2 after union 1, thread it into both effect builds, free it
R("run prepare union2",
  '''        memory_checkpoint("h1 union residual", args.max_rss_gib)
    effects, effect_audit = build_production_effects(
        model,
        tuple(map(str, spec.perts)),
        operator.whitener,
        operator.columns,
        panel_positions,
        pm=pm,
        cd4=cd4,
        h1_union=h1_union,
        h1_union_audit=h1_union_audit,
        weights=weights,
    )
''',
  '''        memory_checkpoint("h1 union residual", args.max_rss_gib)
    if UNION2_SOURCES is None:
        union2 = None
        union2_audit = None
    else:
        # 248: the second union, the same construction as union 1 over its own registry
        # banks, sized on its own carrying rows at `VCC_UNION2_G` inside
        # `apply_frozen_cd4_residual`.
        union2, union2_details = prepare_h1_union(
            tuple(map(str, spec.perts)),
            tuple(map(str, spec.genes)),
            tuple(map(str, spec.perts)),
            max_rss_gib=args.max_rss_gib,
            sources=UNION2_SOURCES,
        )
        union2_audit = h1_union_support_audit(
            union2, union2_details, sources=UNION2_SOURCES, target_g=UNION2_G
        )
        print(
            f"[union2] g={UNION2_G:g}  sources {[name for name, _ in UNION2_SOURCES]}  "
            f"{union2_audit['carrying_targets']}/{len(spec.perts)} targets carried "
            f"(by source count {union2_audit['targets_by_source_count']}), "
            f"{union2_audit['supported_entries']} supported entries",
            flush=True,
        )
        memory_checkpoint("union2 residual", args.max_rss_gib)
    effects, effect_audit = build_production_effects(
        model,
        tuple(map(str, spec.perts)),
        operator.whitener,
        operator.columns,
        panel_positions,
        pm=pm,
        cd4=cd4,
        h1_union=h1_union,
        h1_union_audit=h1_union_audit,
        union2=union2,
        union2_audit=union2_audit,
        weights=weights,
    )
''')
R("run knee build",
  '''                h1_union=h1_union,
                h1_union_audit=h1_union_audit,
                alpha=float(solved[context]["alpha"]),
''',
  '''                h1_union=h1_union,
                h1_union_audit=h1_union_audit,
                union2=union2,
                union2_audit=union2_audit,
                alpha=float(solved[context]["alpha"]),
''')
R("run del",
  '''    del donors, pm, cd4_banks, cd4, h1_union, operator, model
''',
  '''    del donors, pm, cd4_banks, cd4, h1_union, union2, operator, model
''')

# ---------------------------------------------------------------------------------------
# 8. the config block
R("config block",
  '''                    "h1_union_g": H1_UNION_G,
                    "h1_union_sources": [name for name, _ in H1_UNION_SOURCES],
                    "h1_union_sizing_rule": "median over carrying rows",
                }
            ),
            "min_cpm": min_cpm,
''',
  '''                    "h1_union_g": H1_UNION_G,
                    "h1_union_sources": [name for name, _ in H1_UNION_SOURCES],
                    "h1_union_sizing_rule": "median over carrying rows",
                }
            ),
            # Written only when the SECOND union is ON (248), the same convention.
            **(
                {}
                if UNION2_SOURCES is None
                else {
                    "union2_g": UNION2_G,
                    "union2_sources": [name for name, _ in UNION2_SOURCES],
                    "union2_sizing_rule": "median over carrying rows",
                }
            ),
            "min_cpm": min_cpm,
''')


def apply(text: str) -> str:
    for label, old, new in REPLACEMENTS:
        n = text.count(old)
        if n != 1:
            raise SystemExit(f"anchor {label!r} occurs {n} times, expected exactly 1")
        text = text.replace(old, new, 1)
    return text


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def main() -> None:
    src = Path(sys.argv[sys.argv.index("--check") + 1]) if "--check" in sys.argv else TARGET
    text = src.read_text()
    out = apply(text)
    WORK.mkdir(exist_ok=True)
    (WORK / "make_submission_pm_cd4.248.py").write_text(out)
    print(f"source {src}  sha256 {sha(text)}")
    print(f"patched -> {WORK / 'make_submission_pm_cd4.248.py'}  sha256 {sha(out)}  "
          f"({len(REPLACEMENTS)} replacements, +{out.count(chr(10)) - text.count(chr(10))} lines)")
    if "--apply" in sys.argv:
        TARGET.write_text(out)
        print(f"APPLIED to {TARGET}  sha256 {sha(TARGET.read_text())}")


if __name__ == "__main__":
    main()
