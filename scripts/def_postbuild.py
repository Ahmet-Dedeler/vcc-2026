"""Post-build falsifier for a D/E/F package: read the written meta, not the command typed.

AGENTS.md 11.1 rule 10 — "diff the finished artifact's config against the incumbent's and
require EXACTLY the intended number of differences" — plus the two D/E/F-specific traps
experiment 199 found, which no existing check catches:

* **The whitening operator is panel-dependent and fails SILENTLY on a new panel.** Its axis is
  `union(donor coverage on the official axis) − panel`
  (`scripts/build_union_whitener.py:269-273`). `load_production_operator` recomputes the
  column list from the NEW panel, then checks only `n_cols`, the evec shape and the file
  sha256 — all three of which pass with the A/B/C operator, because both panels give 17,942
  columns. Measured on a synthetic D/E/F panel: **98.8% of the 17,942 coordinates are paired
  with the wrong gene.** 235 made the builder assert the ordered-SYMBOL axis digest; this
  script asserts the same digest off the written meta, plus the operator FILE against the
  registry pin (247 — see below).
* **The CD4 banks are panel-restricted**, so `cd4_complete_targets` collapses on a new panel
  unless they are re-extracted. That failure is loud (`VCC_CD4_TARGET_F` raises on a
  zero-norm residual), but the near-miss — a partial re-extraction — is not.

* **The build can silently drop a shipped ARM** (240, from experiment 239 section 2.1). Rule 10's
  diff is only as good as its reference: it was pointed at `rc_hybrid`/`rc_callset`, both from
  before `cd4f034` and `h1union` landed, so a D/E/F build carrying NEITHER arm -- `-0.0199671`
  overall, more than our whole gap to first place -- diffed clean and printed PASS. The arms are
  now asserted by value against `SHIPPED_ARMS` below, the reference defaults to the current
  incumbent, and the reference is itself checked for carrying them.

* **The operator block was blind, and the row count was never compared** (242 #1, #7, #12; fixed
  by 247). The five operator rows read `n_cols`, the recomputed column COUNT and the recomputed
  column-INDEX digest -- none of which reads the operator file, and all of which are identical
  for every operator on disk -- plus a blacklist of exactly one sha. 242 demonstrated that
  `data/donors/whitener_union.npz` (the wrong operator) passed 5/5. The rows now assert the
  field that carries the answer: `operator.axis_sha256` (ordered symbols, 235's construction)
  against the axis recomputed from THIS bundle, `operator.sha256` against the registry pin the
  builder carries for `config.whitener_id` AND against the sha256 of the file at `operator.path`,
  the metadata sha against its pin, and `basis[0]` against the kept-row count. The row-count row
  used to count non-empty contexts while `rows` lived only in the message; it now compares
  `output.n_rows` and every context's `emitted_rows`.

* **CD4's dose inflates with missing coverage** (242 #2; fixed by 247). `solve_cd4_fraction`
  sizes on the all-rows median, which over-doses the carrying rows by
  `median_carrying / median_all` below ~92% coverage (1.196x on section 5's pessimistic draw,
  measured). On a non-A/B/C panel the build must carry `cd4_sizing == "carrying"` at
  `cd4_target_g == DEF_CD4_TARGET_G` -- the A/B/C DELIVERED carrying-row dose -- and on A/B/C it
  must NOT (the legacy path is what 240's row 1 reproduces byte-for-byte).

Everything asserted here is recomputed from the bundle and the donor tables, so the check is
independent of the build it is checking. The only thing read from the builder is its two pin
TABLES (`WHITENER_REGISTRY`, `LEGACY_OPERATOR_AXIS_SHA256`), parsed as literals out of the
source -- the builder is never imported.

    uv run python scripts/def_postbuild.py --meta submissions/<name>.h5ad.meta.json \
        --bundle data/challenge --incumbent submissions/222-h1union-g028.h5ad.meta.json
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import json
import os  # 554: DEF_DAY_ALLOW_LOW_DOSE
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts/make_submission_pm_cd4.py"
ABC_OPERATOR_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
ABC_OPERATOR_COLUMNS_SHA = "b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d"
# The A/B/C axis, ordered symbols (235): the value the recomputation must give on data/challenge.
ABC_AXIS_SHA = "edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f"
# Sorted A/B/C target_gene symbols, newline-joined. Pinned because `vcc datasets download
# controls` overwrites data/challenge/pert_counts.csv (DEF-RUNBOOK, the day page). A comparison
# against that path is true for every bundle that IS data/challenge.
ABC_TARGET_SET_SHA = "c9c4c9a69f76afd4507e9a7619fe9925c19ff34e5878c7854493683b23bea5ca"
PRODUCTION_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
DENSITY_CAP = 4_750_000_000
# Config keys a D/E/F build is ALLOWED to differ from the A/B/C incumbent in. Anything else
# is an unintended second variable (rule 10's two near-misses were both of that shape).
# `expect_contexts` is recorded in the config only when `--expect-contexts` was passed; the
# runbook passes it on purpose (the build must fail if the bundle is not the round you think
# it is), and the A/B/C incumbents never did, so it is an intended difference, not a variable.
ALLOWED_CONFIG_DIFFS = {"contexts", "targets", "expect_contexts"}
# ... and on a NEW panel, the two keys CD4's panel-invariant sizing writes (247, item C).
ALLOWED_CONFIG_DIFFS_NEW_PANEL = {"cd4_sizing", "cd4_target_g"}

# --- 247: CD4's dose on a NEW panel -----------------------------------------------------------
#
# `g_cd4 = 0.34 * (median_carrying / median_all)` measured on A/B/C's real CD4 residual:
# 0.34 * 1.0128397878381297. At this `g` the carrying-row rule reproduces the shipped
# `effects.cd4.fraction` 0.32207201067578145 to 5.6e-17 on A/B/C, and delivers the SAME
# carrying-row dose on any panel. `scripts/def_pipeline.py` sets it on every non-A/B/C bundle;
# this file requires it there, and requires its ABSENCE on A/B/C.
DEF_CD4_SIZING = "carrying"
DEF_CD4_TARGET_G = 0.3443655278649641
# The same floor `scripts/def_census.py` gates on (DEF-RUNBOOK section 5: the pessimistic draw
# has 167), re-asserted off the meta so a build that skipped the census stage still fails here.
CD4_MIN_COMPLETE_TARGETS = 100

# --- THE ARMS THE RUNBOOK SAYS WE SHIP (240) --------------------------------------------------
#
# The rule-10 diff below is only as good as what it is diffing AGAINST. Until 240 the reference
# defaulted to 191's `rc_hybrid` and the driver passed 202's `rc_callset` -- both from before the
# two direction arms landed -- so a D/E/F build that dropped both arms diffed CLEAN and this
# script printed PASS. That is experiment 239 section 2.1, and it was worth -0.0199671 overall on
# the only prize-eligible entry.
#
# So the arms are now asserted BY VALUE, here, in a tracked file, independently of whichever meta
# `--incumbent` points at -- and the reference itself is checked for carrying them, so re-pointing
# it at a stale incumbent fails instead of passing.
#
# config key          value                      arm         board delta   entry
SHIPPED_ARMS: dict[str, tuple[object, str, str, str]] = {
    "alpha":            (0.982, "VCC_ALPHA", "the a11->a15 amplitude series, closed by 211", ""),
    "cd4_target_f":     (0.34, "VCC_CD4_TARGET_F", "+0.0130557", "QSe3JlNUYZ2BzvlcsL14"),
    "h1_union_g":       (0.28, "VCC_H1_UNION_G", "+0.0069114", "iUziqE1gewfcPKzYUDaB"),
    "h1_union_sources": (["gse295214_h1", "h1_2025"], "VCC_UNION_SOURCES",
                         "the H1 pair 229 landed and 235 hardened; unset selects it", ""),
    "emitter":          ("rc_callset", "--emitter", "+0.0084448", "31oB9yWy8CSNBDEkY4iW"),
}
# The shipped package those values were measured on, and the default `--incumbent`.
SHIPPED_META = ROOT / "submissions/249-cd4svd5-f034.h5ad.meta.json"
SHIPPED_BOARD = "cd4svd5, board 0.2160944, rank 3/688, entry Xv0ayvuoHpsKYWv2iu42, 2026-09-05"
# 256 (BR249): the CD4 deflation is the one shipped arm that is PANEL-DEPENDENT -- the day-of
# re-fit k on A/B/C (what the incumbent was built with), the PINNED A/B/C basis anywhere else
# (250/251) -- so it is asserted below by panel, not in SHIPPED_ARMS.
SHIPPED_CD4_DEFLATE_K = 5
SHIPPED_CD4_DEFLATE_BASIS = "V5_abc"
SHIPPED_CD4_DEFLATE = "+0.0042634, entry Xv0ayvuoHpsKYWv2iu42, 2026-09-05"


def arm_rows(cfg: dict, allow_drift: set[str]) -> list[tuple[str, bool, str]]:
    """One check per shipped arm, by value, against `cfg`. Waived only by an explicit name."""
    rows = []
    for key, (want, env, delta, entry) in SHIPPED_ARMS.items():
        got = cfg.get(key, "<absent>")
        ok = got == want
        note = f"{got!r} vs the shipped {want!r} [{env}, {delta}{', entry ' + entry if entry else ''}]"
        if not ok and (key in allow_drift or env in allow_drift):
            rows.append((f"arm {key} DRIFT ALLOWED by --allow-arm-drift", True, note))
        else:
            rows.append((f"arm {key} == the shipped configuration", ok, note))
    return rows


def read_lines(p: Path) -> set[str]:
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def axis_sha256(symbols) -> str:
    """235's canonical operator-axis digest: ordered symbols, newline-joined UTF-8, no tail."""
    return hashlib.sha256("\n".join(map(str, symbols)).encode("utf-8")).hexdigest()


def builder_literal_or_empty(name: str, source: Path = BUILDER) -> dict:
    """`builder_literal`, but `{}` when the builder predates the table (251: the pinned-basis
    registry exists from `V5_abc` on; a build carrying `cd4_deflate_basis` against an older
    builder source then fails the 'is registered' row instead of crashing the check)."""
    try:
        return builder_literal(name, source)
    except SystemExit:
        return {}


def builder_literal(name: str, source: Path = BUILDER) -> dict:
    """A module-level dict LITERAL out of the builder's source, without importing it.

    The registry pin is what the runbook's step 7 rewrites on the day
    (`def_pin_panel_artifacts.py --apply`), so the operator a D/E/F build must carry is, by
    definition, the one pinned under `config.whitener_id` in this table.
    """
    tree = ast.parse(source.read_text())
    for node in tree.body:
        targets: list[str] = []
        if isinstance(node, ast.Assign):
            targets = [t.id for t in node.targets if isinstance(t, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target.id]
        if name in targets and node.value is not None:
            value = ast.literal_eval(node.value)
            if not isinstance(value, dict):
                raise SystemExit(f"{source}: {name} is not a dict literal")
            return value
    raise SystemExit(f"{source}: no module-level literal named {name}")


# 554: the K-class emitter can deliver its CLEAN4 block prior and coarse MDE (with the uncalled-gene hook riding on it) at
# dose 0 WITHOUT an error (553). Floors from 554's calibration (its RESULTS.md): healthy per-context medians are prior 1.0,
# MDE 0.951-1.0 (A / D the lowest), call-set misses 4-21.5; 429a's plain-MDE failure had A's MDE median 0.0.
DOSE_PRIOR_MIN = 0.9
DOSE_MDE_MIN = 0.9
DOSE_MISSED_MAX = 100.0
DOSE_KEYS = ("prior_lambda", "mde_lambda", "call_set_missed")


def dose_guard(meta: dict, rc_audit: Path | None, contexts: list[str]) -> dict:
    """Per-context medians of the emitter's delivered steering dose, from the builder's per-block VCC_RC_AUDIT_DUMP record
    (or, without one, the meta's POOLED medians: a global collapse only). Returns the record and the check rows (a row only
    when a floor trips or the dump does not belong to the meta, so a passing build's text is today's)."""
    rc = meta.get("output", {}).get("rc_callset")
    floors = {"prior_lambda": DOSE_PRIOR_MIN, "mde_lambda": DOSE_MDE_MIN, "call_set_missed": DOSE_MISSED_MAX}
    rec: dict = {"floors": floors, "source": None, "per_context": {}, "rows": []}
    if not isinstance(rc, dict):
        rec["source"] = "none (no rc_callset record: not a K-class build)"
        return rec
    allow = os.environ.get("DEF_DAY_ALLOW_LOW_DOSE", "") == "1"
    problems: list[str] = []
    if rc_audit is not None:
        rows = [json.loads(ln) for ln in rc_audit.read_text().splitlines() if ln.strip()]
        rec["source"] = f"rc_audit {rc_audit.name} ({len(rows)} blocks)"
        rec["rc_audit_sha256"] = sha256_file(rc_audit)
        # the dump must be THIS build's: same block count, and its pooled medians are the meta's (same values, np.median)
        same = len(rows) == rc.get("blocks")
        for k in DOSE_KEYS:
            if k in rc and same:
                v = [float(r[k]) for r in rows if k in r]
                same = bool(v) and float(np.median(v)) == float(rc[k]["median"])
        if not same:
            problems.append(f"--rc-audit {rc_audit.name} does not belong to this meta (block count or pooled medians differ)")
        for c in contexts:
            sub = [r for r in rows if str(r.get("context")) == str(c)]
            rec["per_context"][c] = {"blocks": len(sub), **{k: (float(np.median([float(r[k]) for r in sub if k in r]))
                                                               if any(k in r for r in sub) else None) for k in DOSE_KEYS}}
            if not sub:
                problems.append(f"{c}: no blocks in the dump")
    else:
        rec["source"] = "meta pooled medians (no --rc-audit: a one-context collapse is NOT visible here)"
        rec["per_context"]["pooled"] = {"blocks": rc.get("blocks"),
                                        **{k: (float(rc[k]["median"]) if k in rc else None) for k in DOSE_KEYS}}
    for c, v in rec["per_context"].items():
        m = {}
        for k in DOSE_KEYS:
            if v.get(k) is None:
                continue
            m[k] = (v[k] - floors[k]) if k != "call_set_missed" else (floors[k] - v[k])
            if m[k] < 0:
                problems.append(f"{c}: median {k} {v[k]:.4g} past its floor {floors[k]:g}")
        v["margin"] = m
    rec["allowed"] = allow if problems else None
    if problems:
        rec["rows"].append(("emitter steering dose per context: median prior_lambda >= "
                            f"{DOSE_PRIOR_MIN}, mde_lambda >= {DOSE_MDE_MIN}, call_set_missed <= {DOSE_MISSED_MAX:g} (554)",
                            allow, "; ".join(problems) + (" -- ALLOWED by DEF_DAY_ALLOW_LOW_DOSE=1" if allow else "")))
    return rec


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--meta", type=Path, required=True)
    ap.add_argument("--bundle", type=Path, default=ROOT / "data/challenge")
    ap.add_argument("--incumbent", type=Path, default=SHIPPED_META,
                    help=f"the CURRENT incumbent's meta ({SHIPPED_BOARD}). It is checked for "
                         "carrying the shipped arms itself, so pointing this at a stale "
                         "incumbent (191, 202) FAILS rather than passing a build that dropped them")
    ap.add_argument("--cd4-root", type=Path, default=CD4_ROOT)
    ap.add_argument("--builder", type=Path, default=BUILDER,
                    help="the builder SOURCE whose WHITENER_REGISTRY / LEGACY_OPERATOR_AXIS_SHA256 "
                         "literals are the operator pins (parsed, never imported)")
    ap.add_argument("--allow-config-diff", default="",
                    help="comma-separated extra config keys this build INTENDS to move")
    ap.add_argument("--require-config-diff", default="",
                    help="538: comma-separated config keys that MUST move against the incumbent. "
                         "Fails unless every listed key moved AND nothing outside the allow-list "
                         "(built-in allows, --allow-config-diff, the new-panel keys) moved. "
                         "Empty, the default, adds no check and does not change today's output")
    ap.add_argument("--allow-arm-drift", default="",
                    help="comma-separated arms this build DELIBERATELY does not ship, by config "
                         "key or env name (e.g. h1_union_g or VCC_H1_UNION_G; cd4_sizing waives "
                         "the panel-invariant CD4 rows). Each one waives a board-validated gain; "
                         "say why in the round's RESULTS.md")
    ap.add_argument("--rc-audit", type=Path, default=None,
                    help="554: the build's per-block VCC_RC_AUDIT_DUMP jsonl; per-context steering-dose floors. "
                         "Absent: the meta's pooled medians (a global collapse only); silent on PASS either way")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()

    meta = json.loads(args.meta.read_text())
    if not args.incumbent.is_file():
        raise SystemExit(
            f"the incumbent reference {args.incumbent} does not exist.\n"
            f"  It must be the meta of the CURRENTLY shipped package ({SHIPPED_BOARD}).\n"
            "  Its config is what rule 10 diffs against; without it this check is meaningless.\n"
            "  DEF-RUNBOOK section 1 lists the packages and section 2 step 10 the command.")
    inc = json.loads(args.incumbent.read_text())
    bundle = args.bundle if args.bundle.is_absolute() else ROOT / args.bundle
    manifest = json.loads((bundle / "manifest.json").read_text())
    panel = [str(t) for t in pd.read_csv(bundle / "pert_counts.csv")["target_gene"]]
    genes = [str(g) for g in pd.read_csv(bundle / "gene_names.csv")["gene_name"]]
    panel_set, gene_set = set(panel), set(genes)
    tindex = {g: i for i, g in enumerate(genes)}
    registry = builder_literal("WHITENER_REGISTRY", args.builder)
    legacy_axis = builder_literal("LEGACY_OPERATOR_AXIS_SHA256", args.builder)
    drift = {k.strip() for k in args.allow_arm_drift.split(",") if k.strip()}

    checks: list[tuple[str, bool, str]] = []

    def chk(name, ok, value=""):
        checks.append((name, bool(ok), str(value)))

    # --- contexts ------------------------------------------------------------------------
    chk("config.contexts == the manifest's",
        list(meta["config"]["contexts"]) == list(manifest["contexts"]),
        f"{meta['config']['contexts']} vs {manifest['contexts']}")
    chk("config.targets == n_constructs",
        meta["config"]["targets"] == manifest["n_constructs"], meta["config"]["targets"])

    # --- the operator: by registry pin, by file hash, and by recomputed SYMBOL axis --------
    union: set[str] = set()
    for d in PRODUCTION_DONORS:
        union |= read_lines(ROOT / f"data/donors/{d}/genes.txt") & gene_set
    expect_cols = np.asarray(sorted(tindex[g] for g in union - panel_set), dtype=np.int64)
    expect_cols_sha = hashlib.sha256(
        np.ascontiguousarray(expect_cols, dtype="<i8").tobytes()).hexdigest()
    expect_axis_sha = axis_sha256(genes[i] for i in expect_cols)
    op = meta["operator"]
    new_panel = list(manifest["contexts"]) != ["A", "B", "C"]
    wid = meta["config"].get("whitener_id")
    pin = registry.get(wid)
    chk("config.whitener_id is in the builder's WHITENER_REGISTRY", pin is not None,
        f"{wid!r} (registered: {sorted(registry)})")
    chk(f"operator.sha256 == WHITENER_REGISTRY[{wid!r}][0] (the FILE the builder pinned for this panel)",
        pin is not None and op["sha256"] == pin[0],
        f"{op['sha256'][:16]}… vs pin {pin[0][:16] + '…' if pin else None}")
    chk(f"operator.metadata_sha256 == WHITENER_REGISTRY[{wid!r}][1] (its metadata pair)",
        pin is not None and op.get("metadata_sha256") == pin[1],
        f"{str(op.get('metadata_sha256'))[:16]}… vs pin {pin[1][:16] + '…' if pin else None}")
    op_path = Path(op["path"])
    if op_path.is_file():
        on_disk = sha256_file(op_path)
        chk("sha256(the file at operator.path) == operator.sha256 (read the pin, hash the file)",
            on_disk == op["sha256"], f"{on_disk[:16]}… vs meta {op['sha256'][:16]}… ({op_path.name})")
    else:
        chk("the file at operator.path exists, so the pinned operator can be hashed",
            False, f"MISSING {op_path}")
    chk("operator.sha256 is NOT the A/B/C operator (it is panel-dependent; rebuild it)",
        (op["sha256"] != ABC_OPERATOR_SHA) if new_panel else True, op["sha256"])
    chk("operator.columns == the axis recomputed from THIS panel",
        op["columns"] == len(expect_cols), f"{op['columns']} vs {len(expect_cols)}")
    chk("operator.columns_sha256 == the axis recomputed from THIS panel",
        op["columns_sha256"] == expect_cols_sha,
        f"{op['columns_sha256'][:16]}… vs {expect_cols_sha[:16]}…")
    # THE ROW THAT CARRIES THE ANSWER (242 #1): the ordered-SYMBOL digest, 235's construction.
    recorded_axis = op.get("axis_sha256")
    axis_source = op.get("axis_sha256_checked_against")
    if recorded_axis:
        chk("operator.axis_sha256 == axis_sha256(ordered symbols) recomputed from THIS bundle",
            recorded_axis == expect_axis_sha,
            f"{recorded_axis[:16]}… vs recomputed {expect_axis_sha[:16]}…")
        real_sources = {f"metadata:{Path(str(op.get('metadata_path', ''))).name}",
                        "LEGACY_OPERATOR_AXIS_SHA256[operator sha256]"}
        chk("operator.axis_sha256_checked_against names a real source",
            axis_source in real_sources, f"{axis_source!r}")
    else:
        # A meta written before 235 carries no axis field. The builder's own fallback for such
        # operators is a table keyed by operator sha256; resolve it the same way, and FAIL --
        # never skip -- when neither the meta nor the table can say which panel it was fitted for.
        legacy = legacy_axis.get(op["sha256"])
        chk("operator.axis_sha256 ABSENT (meta predates 235): "
            "LEGACY_OPERATOR_AXIS_SHA256[operator.sha256] == axis recomputed from THIS bundle",
            legacy is not None and legacy == expect_axis_sha,
            f"legacy {legacy[:16] + '…' if legacy else 'NOT IN TABLE'} vs recomputed "
            f"{expect_axis_sha[:16]}… (no axis digest in the meta)")
    chk("operator.metadata.n_cols == operator.columns",
        op["metadata"].get("n_cols") == op["columns"], op["metadata"].get("n_cols"))
    md = op["metadata"]
    kept = md.get("rows_kept_per_donor")
    if isinstance(kept, dict) and kept:
        chk("operator.metadata.basis[0] == sum(rows_kept_per_donor) (the ROW count, 242 #12)",
            int(md["basis"][0]) == sum(int(v) for v in kept.values()),
            f"{md['basis'][0]} vs {sum(int(v) for v in kept.values())}")
    else:
        chk("operator.metadata.basis recorded (unpruned operator: no per-donor row table to cross-check)",
            isinstance(md.get("basis"), list) and md.get("basis_prune") is None, md.get("basis"))
    chk("all_panel_coordinates_excluded", op["all_panel_coordinates_excluded"], "")
    # The axis is union(donor coverage) - panel, a function of the TARGET set, not of the context labels.
    # abc_targets is the PINNED A/B/C target digest, not whatever data/challenge/pert_counts.csv currently
    # is: the Oct 22 download replaces that file (2026-09-28 review, 519). An unchanged symbol axis is
    # possible when the genes that differ are disjoint from the donor union; that is not a failed build.
    # The row that binds the operator to THIS panel is the symbol-axis equality above.
    abc_targets = hashlib.sha256("\n".join(sorted(panel_set)).encode()).hexdigest() == ABC_TARGET_SET_SHA
    if new_panel and abc_targets:
        chk("new contexts on A/B/C's own 300 targets: the recomputed SYMBOL axis is A/B/C's",
            expect_axis_sha == ABC_AXIS_SHA, f"symbols {expect_axis_sha[:16]}…")
    if not new_panel:
        chk("the recomputed SYMBOL axis is A/B/C's known value (this file's own known-answer row)",
            expect_axis_sha == ABC_AXIS_SHA, f"{expect_axis_sha[:16]}…")

    # --- panel_support, recomputed from the donor tables ----------------------------------
    k562 = read_lines(ROOT / "data/donors/k562_gwps/perts.txt")
    pm = read_lines(ROOT / "data/external/perturbmatch/perts.txt")
    cd4 = {n: read_lines(args.cd4_root / n / "perts.txt") for n in CD4_NAMES
           if (args.cd4_root / n / "perts.txt").exists()}
    expect_pm = sum(1 for t in panel if t in k562 and t in pm)
    expect_cd4 = (sum(1 for t in panel if all(t in s for s in cd4.values()))
                  if len(cd4) == len(CD4_NAMES) else -1)
    got_pm = int(meta["effects"]["perturbmatch"]["both_supported_modified_count"])
    got_cd4 = int(meta["effects"]["cd4"]["fully_supported_targets"])
    chk("perturbmatch replacement count == set-intersection prediction",
        got_pm == expect_pm, f"{got_pm} vs {expect_pm}")
    chk("cd4 complete targets == set-intersection prediction",
        got_cd4 == expect_cd4, f"{got_cd4} vs {expect_cd4}")
    chk("cd4 complete targets > 0 (the banks were re-extracted for this panel)",
        got_cd4 > 0, got_cd4)
    chk(f"cd4 complete targets >= {CD4_MIN_COMPLETE_TARGETS} (def_census.py's coverage floor, re-asserted)",
        got_cd4 >= CD4_MIN_COMPLETE_TARGETS,
        f"{got_cd4} = {got_cd4 / manifest['n_constructs']:.1%} of {manifest['n_constructs']}")
    fb = meta.get("frozen_bundle") or {}
    chk("frozen_bundle.panel_support is pinned and matches",
        (fb.get("panel_support") or {}).get("cd4_complete_targets") == got_cd4
        and (fb.get("panel_support") or {}).get("perturbmatch_replacement_targets") == got_pm,
        json.dumps(fb.get("panel_support")))

    # --- 247: CD4's sizing rule must be panel-invariant off A/B/C, and legacy on A/B/C ------
    cfg = meta["config"]
    rf = meta["effects"]["cd4"].get("residual_frame") or {}
    sizing = cfg.get("cd4_sizing")
    if "cd4_sizing" in drift or "VCC_CD4_SIZING" in drift:
        chk("cd4 sizing rows DRIFT ALLOWED by --allow-arm-drift", True,
            f"config.cd4_sizing={sizing!r}, cd4_target_g={cfg.get('cd4_target_g')!r}")
    elif new_panel:
        chk("config.cd4_sizing == 'carrying' on a non-A/B/C panel "
            "(the all-rows rule over-doses CD4 below ~92% coverage; 242 #2)",
            sizing == DEF_CD4_SIZING, f"{sizing!r}")
        chk(f"config.cd4_target_g == {DEF_CD4_TARGET_G!r} (A/B/C's DELIVERED carrying-row dose, 247)",
            cfg.get("cd4_target_g") == DEF_CD4_TARGET_G, f"{cfg.get('cd4_target_g')!r}")
    else:
        chk("config.cd4_sizing ABSENT on A/B/C (the legacy path 240's row 1 reproduces byte-for-byte)",
            sizing is None and "cd4_target_g" not in cfg,
            f"cd4_sizing={sizing!r}, cd4_target_g={cfg.get('cd4_target_g', '<absent>')!r}")
    chk("effects.cd4.residual_frame's sizing rule matches config.cd4_sizing",
        (rf.get("sizing_rule") == "carrying") == (sizing == "carrying"),
        f"frame {rf.get('sizing_rule', '<all-rows median>')!r} vs config {sizing!r}")

    # --- 251: a PINNED CD4 deflation basis, if (and only if) the build carries one ---------
    #
    # `VCC_CD4_DEFLATE_BASIS=<name>` writes `config.cd4_deflate_basis` and an
    # `effects.cd4.deflation` block with `mode == "pinned"`.  Rule 10 would refuse the new config
    # key as an unintended variable; it is accepted here ONLY when every pin it carries checks
    # out against the builder's `CD4_DEFLATION_BASIS_REGISTRY` (read as a literal, never
    # imported) and the file on disk, and the census diagnostic 250 section 6.2 names is echoed.
    # A build without the key sees none of these rows, so 222/249/244's metas are untouched.
    extra_allowed: set[str] = set()
    basis_name = cfg.get("cd4_deflate_basis")
    defl = meta["effects"]["cd4"].get("deflation") or {}
    if basis_name is not None:
        basis_registry = builder_literal_or_empty("CD4_DEFLATION_BASIS_REGISTRY", args.builder)
        entry = basis_registry.get(basis_name)
        chk("config.cd4_deflate_basis is in the builder's CD4_DEFLATION_BASIS_REGISTRY",
            entry is not None, f"{basis_name!r} (registered: {sorted(basis_registry)})")
        chk("config.cd4_deflate_k ABSENT while a pinned basis is on (one deflation variable, rule 10)",
            "cd4_deflate_k" not in cfg, f"cd4_deflate_k={cfg.get('cd4_deflate_k', '<absent>')!r}")
        chk("effects.cd4.deflation.mode == 'pinned' and .basis == config.cd4_deflate_basis",
            defl.get("mode") == "pinned" and defl.get("basis") == basis_name,
            f"mode={defl.get('mode')!r} basis={defl.get('basis')!r}")
        chk(f"effects.cd4.deflation.basis_sha256 == CD4_DEFLATION_BASIS_REGISTRY[{basis_name!r}]['sha256'] (the pin)",
            entry is not None and defl.get("basis_sha256") == entry.get("sha256"),
            f"{str(defl.get('basis_sha256'))[:16]}… vs pin {entry['sha256'][:16] + '…' if entry else None}")
        chk(f"effects.cd4.deflation.basis_v_sha256 == the registry's v_sha256 and .k == its k",
            entry is not None and defl.get("basis_v_sha256") == entry.get("v_sha256")
            and str(defl.get("k")) == str(entry.get("k")),
            f"V {str(defl.get('basis_v_sha256'))[:16]}… k={defl.get('k')!r} vs pin "
            f"{entry['v_sha256'][:16] + '…' if entry else None} k={entry.get('k') if entry else None!r}")
        bpath = (ROOT / str(entry["path"])) if entry else None
        if bpath is not None and bpath.is_file():
            on_disk_b = sha256_file(bpath)
            chk("sha256(the file at the basis pin's path) == the pin (read the pin, hash the file)",
                on_disk_b == entry["sha256"], f"{on_disk_b[:16]}… vs pin {entry['sha256'][:16]}… ({bpath.name})")
        else:
            chk("the pinned basis file exists at the registry's path", False, f"MISSING {bpath}")
        chk("effects.cd4.deflation.gene_axis_verified (the builder checked the basis's gene ORDER against this bundle)",
            defl.get("gene_axis_verified") is True, f"{defl.get('gene_axis_verified')!r}")
        captured = defl.get("captured_frac_of_centred")
        chk("effects.cd4.deflation.captured_frac_of_centred is a number in (0, 1) (250's census diagnostic, echoed below)",
            isinstance(captured, (int, float)) and 0.0 < float(captured) < 1.0, f"{captured!r}")
        extra_allowed.add("cd4_deflate_basis")
    elif defl.get("mode") == "pinned":
        chk("effects.cd4.deflation.mode == 'pinned' requires config.cd4_deflate_basis", False,
            "the audit says pinned but the config carries no basis name")

    # --- 256: the CD4 deflation, BY VALUE and BY PANEL (BR249: SHIPS) ---------------------
    if new_panel:
        if "VCC_CD4_DEFLATE_BASIS" in drift or "cd4_deflate_basis" in drift:
            chk("arm cd4_deflate_basis DRIFT ALLOWED by --allow-arm-drift", True, f"{basis_name!r}")
        else:
            chk(f"config.cd4_deflate_basis == {SHIPPED_CD4_DEFLATE_BASIS!r} on a non-A/B/C panel "
                f"(the pinned A/B/C basis SHIPS, {SHIPPED_CD4_DEFLATE})",
                basis_name == SHIPPED_CD4_DEFLATE_BASIS, f"{basis_name!r}")
        # The intended swap: the incumbent carries the re-fit k, the new panel the pin of the SAME k.
        if basis_name is not None and inc["config"].get("cd4_deflate_k") == SHIPPED_CD4_DEFLATE_K:
            extra_allowed.add("cd4_deflate_k")
    else:
        if "VCC_CD4_DEFLATE_K" in drift or "cd4_deflate_k" in drift:
            chk("arm cd4_deflate_k DRIFT ALLOWED by --allow-arm-drift", True,
                f"{cfg.get('cd4_deflate_k', '<absent>')!r}")
        else:
            chk(f"config.cd4_deflate_k == {SHIPPED_CD4_DEFLATE_K} on A/B/C (the re-fit the incumbent "
                f"was built with, {SHIPPED_CD4_DEFLATE})",
                cfg.get("cd4_deflate_k") == SHIPPED_CD4_DEFLATE_K, f"{cfg.get('cd4_deflate_k', '<absent>')!r}")
        chk("config.cd4_deflate_basis ABSENT on A/B/C (240's row 1 stays byte-for-byte)",
            basis_name is None, f"{basis_name!r}")

    # --- 240: the arms the runbook says we ship, BY VALUE ---------------------------------
    # Asserted against `SHIPPED_ARMS` in this file, not against the reference, so this cannot be
    # defeated by pointing `--incumbent` at an older package -- which is exactly how the
    # pre-240 check passed a build that had silently dropped both direction arms.
    for name, ok, value in arm_rows(meta["config"], drift):
        chk(name, ok, value)
    # And the reference itself must carry them, or it is not the current incumbent.
    stale = [k for k, (want, *_) in SHIPPED_ARMS.items() if inc["config"].get(k, "<absent>") != want]
    if inc["config"].get("cd4_deflate_k") != SHIPPED_CD4_DEFLATE_K:   # 256: the reference predates cd4svd5
        stale.append("cd4_deflate_k")
    chk(f"the --incumbent reference {args.incumbent.name} IS the current incumbent "
        f"({SHIPPED_BOARD})", not stale,
        "carries every shipped arm" if not stale else
        f"MISSING {stale} -- this reference predates them, so the diff below cannot see them")

    # --- rule 10: the config diff --------------------------------------------------------
    allowed = ALLOWED_CONFIG_DIFFS | {k.strip() for k in args.allow_config_diff.split(",") if k.strip()}
    if new_panel:
        allowed |= ALLOWED_CONFIG_DIFFS_NEW_PANEL
    allowed |= extra_allowed          # 251: `cd4_deflate_basis`, only when its rows above ran
    moved = sorted(k for k in set(meta["config"]) | set(inc["config"])
                   if meta["config"].get(k, "<absent>") != inc["config"].get(k, "<absent>"))
    unexpected = [k for k in moved if k not in allowed]
    chk(f"config differences vs {args.incumbent.name} are only {sorted(allowed)}",
        not unexpected, f"moved={moved}")
    # 538: empty --require-config-diff adds nothing, so the checks above and the exit code
    # stay today's. A listed key is not auto-allowed: --allow-config-diff keeps its meaning.
    required = [k.strip() for k in args.require_config_diff.split(",") if k.strip()]
    if required:
        missing = [k for k in required if k not in moved]
        outside = [k for k in moved if k not in allowed]
        chk("every --require-config-diff key moved and nothing outside the allow-list did",
            not missing and not outside,
            f"required={required} missing={missing} moved={moved} outside={outside}")

    # --- the emitted artifact ------------------------------------------------------------
    out = meta["output"]
    per_ctx = manifest["n_constructs"] * manifest["cells_per_pert"]
    rows = per_ctx * len(manifest["contexts"])
    chk("nnz under the density cap", out["nnz"] <= DENSITY_CAP,
        f"{out['nnz']:,} = {out['nnz'] / DENSITY_CAP:.1%} of cap")
    # 242 #7: the old row counted non-empty contexts and used `rows` only in its message.
    chk("output.n_rows == n_constructs x cells_per_pert x contexts",
        out.get("n_rows") == rows, f"{out.get('n_rows')} vs expected {rows:,}")
    ctx_out = out.get("contexts", {})
    chk("every manifest context emitted exactly n_constructs x cells_per_pert rows, and no other context",
        set(ctx_out) == set(manifest["contexts"])
        and all(ctx_out[c].get("emitted_rows") == per_ctx and ctx_out[c].get("stored_entries", 0) > 0
                for c in manifest["contexts"]),
        {c: ctx_out.get(c, {}).get("emitted_rows") for c in sorted(set(ctx_out) | set(manifest["contexts"]))})
    # 554: the K-class emitter's delivered steering dose (dose_guard). A row only when it trips or is ALLOWED.
    dose = dose_guard(meta, args.rc_audit, [str(c) for c in manifest["contexts"]])
    for _name, _ok, _value in dose["rows"]:
        chk(_name, _ok, _value)

    print(f"meta      {args.meta}")
    print(f"incumbent {args.incumbent}")
    print(f"bundle    {bundle}  panel {manifest.get('panel_id')}  contexts {manifest['contexts']}")
    print(f"operator  {op['sha256']}  columns {op['columns']}  whitener_id {wid!r}")
    print(f"axis      recomputed {expect_axis_sha}  meta {recorded_axis or '<absent>'}")
    print(f"whitening after_all_median_norm {meta['effects']['whitening']['after_all_median_norm']} "
          f"(A/B/C incumbent {inc['effects']['whitening']['after_all_median_norm']}; "
          f"reported, NOT asserted — it is a panel-dependent quantity)")
    if rf.get("sizing_rule") == "carrying":
        print(f"cd4       carrying rule g={cfg.get('cd4_target_g')} solved {rf.get('solved_fraction')} on "
              f"{rf.get('carrying_rows')} rows ({rf.get('coverage_frac', 0):.1%}); the all-rows rule "
              f"would solve {rf.get('legacy_all_rows_rule_would_solve')!r} "
              f"(over-dose x{rf.get('overdose_factor_all_rows_rule')!r})")
    else:
        print(f"cd4       all-rows rule f={cfg.get('cd4_target_f')} solved {rf.get('solved_fraction')} "
              f"({got_cd4} complete targets = {got_cd4 / manifest['n_constructs']:.1%} coverage)")
    if basis_name is not None:
        # 251: the census diagnostic, echoed off the meta (250 section 6.2: 0.10-0.14 on the four
        # disjoint panels, 0.198 on A/B/C; < 0.05 or > 0.19 is unmeasured territory; the day's own
        # top-k share 0.198-0.232 everywhere and NOT a gate).
        print(f"cd4       deflation PINNED {basis_name!r} sha256 {str(defl.get('basis_sha256'))[:16]}… k={defl.get('k')}: "
              f"{defl.get('n_carried')} carried rows, energy removed {defl.get('energy_removed_frac')!r}, "
              f"captured {defl.get('captured_frac_of_centred')!r} of the day's centred CD4 energy "
              f"(250: 0.10-0.14 disjoint, 0.198 A/B/C), own top-k share "
              f"{defl.get('own_top_k_energy_frac_of_centred')!r} (recorded, not gated)")
    print()
    failed = [c for c in checks if not c[1]]
    for name, okc, value in checks:
        print(f"  {'PASS' if okc else 'FAIL'}  {name}: {value}")
    print(f"\n{len(checks) - len(failed)}/{len(checks)} checks pass")
    if any(c[0].startswith(("arm ", "config.cd4_sizing", "config.cd4_target_g")) or "incumbent" in c[0]
           for c in failed):
        print("\n  *** THIS PACKAGE DOES NOT CARRY THE CONFIGURATION THE RUNBOOK SHIPS. ***\n"
              f"  The shipped one is {SHIPPED_BOARD}:\n"
              + "".join(f"    {k:<18} {want!r}   [{env}, {delta}]\n"
                        for k, (want, env, delta, _) in SHIPPED_ARMS.items())
              + f"    {'cd4_sizing':<18} {DEF_CD4_SIZING!r} at cd4_target_g {DEF_CD4_TARGET_G!r} on a "
              "NON-A/B/C panel only [VCC_CD4_SIZING/VCC_CD4_TARGET_G, 247]\n"
              "  Rebuild with the driver's defaults (scripts/def_pipeline.py), or, if a drift is\n"
              "  deliberate, name it in --allow-arm-drift and justify it in the round's RESULTS.md.\n"
              "  Do NOT submit this package.")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(
            {"meta": str(args.meta), "checks": [{"check": c, "pass": p, "value": v}
                                                for c, p, v in checks],
             "dose_guard": {k: v for k, v in dose.items() if k != "rows"}}, indent=2) + "\n")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
