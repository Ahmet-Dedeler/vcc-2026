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
  with the wrong gene.** The only thing that distinguishes a correct build is that the
  operator FILE changed, so this asserts `operator.sha256` is not the A/B/C one.
* **The CD4 banks are panel-restricted**, so `cd4_complete_targets` collapses on a new panel
  unless they are re-extracted. That failure is loud (`VCC_CD4_TARGET_F` raises on a
  zero-norm residual), but the near-miss — a partial re-extraction — is not.

* **The build can silently drop a shipped ARM** (240, from experiment 239 section 2.1). Rule 10's
  diff is only as good as its reference: it was pointed at `rc_hybrid`/`rc_callset`, both from
  before `cd4f034` and `h1union` landed, so a D/E/F build carrying NEITHER arm -- `-0.0199671`
  overall, more than our whole gap to first place -- diffed clean and printed PASS. The arms are
  now asserted by value against `SHIPPED_ARMS` below, the reference defaults to the current
  incumbent, and the reference is itself checked for carrying them.

Everything asserted here is recomputed from the bundle and the donor tables, so the check is
independent of the build it is checking.

    uv run python scripts/def_postbuild.py --meta submissions/<name>.h5ad.meta.json \
        --bundle data/challenge --incumbent submissions/222-h1union-g028.h5ad.meta.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ABC_OPERATOR_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
ABC_OPERATOR_COLUMNS_SHA = "b2cae0d381c3cd7065a3de33f1a6da5514577f2b3187b6dc84786ee789e27e4d"
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
SHIPPED_META = ROOT / "submissions/222-h1union-g028.h5ad.meta.json"
SHIPPED_BOARD = "h1union, board 0.2118310, rank 3/656, entry iUziqE1gewfcPKzYUDaB, 2026-09-04"


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
    ap.add_argument("--allow-config-diff", default="",
                    help="comma-separated extra config keys this build INTENDS to move")
    ap.add_argument("--allow-arm-drift", default="",
                    help="comma-separated arms this build DELIBERATELY does not ship, by config "
                         "key or env name (e.g. h1_union_g or VCC_H1_UNION_G). Each one waives a "
                         "board-validated gain; say why in the round's RESULTS.md")
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

    checks: list[tuple[str, bool, str]] = []

    def chk(name, ok, value=""):
        checks.append((name, bool(ok), str(value)))

    # --- contexts ------------------------------------------------------------------------
    chk("config.contexts == the manifest's",
        list(meta["config"]["contexts"]) == list(manifest["contexts"]),
        f"{meta['config']['contexts']} vs {manifest['contexts']}")
    chk("config.targets == n_constructs",
        meta["config"]["targets"] == manifest["n_constructs"], meta["config"]["targets"])

    # --- the operator, by sha and by recomputed axis --------------------------------------
    union: set[str] = set()
    for d in PRODUCTION_DONORS:
        union |= read_lines(ROOT / f"data/donors/{d}/genes.txt") & gene_set
    expect_cols = np.asarray(sorted(tindex[g] for g in union - panel_set), dtype=np.int64)
    expect_cols_sha = __import__("hashlib").sha256(
        np.ascontiguousarray(expect_cols, dtype="<i8").tobytes()).hexdigest()
    op = meta["operator"]
    new_panel = list(manifest["contexts"]) != ["A", "B", "C"]
    chk("operator.sha256 is NOT the A/B/C operator (it is panel-dependent; rebuild it)",
        (op["sha256"] != ABC_OPERATOR_SHA) if new_panel else True, op["sha256"])
    chk("operator.columns == the axis recomputed from THIS panel",
        op["columns"] == len(expect_cols), f"{op['columns']} vs {len(expect_cols)}")
    chk("operator.columns_sha256 == the axis recomputed from THIS panel",
        op["columns_sha256"] == expect_cols_sha,
        f"{op['columns_sha256'][:16]}… vs {expect_cols_sha[:16]}…")
    chk("operator.metadata.n_cols == operator.columns",
        op["metadata"].get("n_cols") == op["columns"], op["metadata"].get("n_cols"))
    chk("all_panel_coordinates_excluded", op["all_panel_coordinates_excluded"], "")
    if new_panel and op["columns_sha256"] == ABC_OPERATOR_COLUMNS_SHA:
        chk("the recomputed axis moved off A/B/C's", False, "IDENTICAL to A/B/C — impossible")

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
    fb = meta.get("frozen_bundle") or {}
    chk("frozen_bundle.panel_support is pinned and matches",
        (fb.get("panel_support") or {}).get("cd4_complete_targets") == got_cd4
        and (fb.get("panel_support") or {}).get("perturbmatch_replacement_targets") == got_pm,
        json.dumps(fb.get("panel_support")))

    # --- 240: the arms the runbook says we ship, BY VALUE ---------------------------------
    # Asserted against `SHIPPED_ARMS` in this file, not against the reference, so this cannot be
    # defeated by pointing `--incumbent` at an older package -- which is exactly how the
    # pre-240 check passed a build that had silently dropped both direction arms.
    drift = {k.strip() for k in args.allow_arm_drift.split(",") if k.strip()}
    for name, ok, value in arm_rows(meta["config"], drift):
        chk(name, ok, value)
    # And the reference itself must carry them, or it is not the current incumbent.
    stale = [k for k, (want, *_) in SHIPPED_ARMS.items() if inc["config"].get(k, "<absent>") != want]
    chk(f"the --incumbent reference {args.incumbent.name} IS the current incumbent "
        f"({SHIPPED_BOARD})", not stale,
        "carries every shipped arm" if not stale else
        f"MISSING {stale} -- this reference predates them, so the diff below cannot see them")

    # --- rule 10: the config diff --------------------------------------------------------
    allowed = ALLOWED_CONFIG_DIFFS | {k.strip() for k in args.allow_config_diff.split(",") if k.strip()}
    moved = sorted(k for k in set(meta["config"]) | set(inc["config"])
                   if meta["config"].get(k, "<absent>") != inc["config"].get(k, "<absent>"))
    unexpected = [k for k in moved if k not in allowed]
    chk(f"config differences vs {args.incumbent.name} are only {sorted(allowed)}",
        not unexpected, f"moved={moved}")

    # --- the emitted artifact ------------------------------------------------------------
    out = meta["output"]
    rows = manifest["n_constructs"] * manifest["cells_per_pert"] * len(manifest["contexts"])
    chk("nnz under the density cap", out["nnz"] <= DENSITY_CAP,
        f"{out['nnz']:,} = {out['nnz'] / DENSITY_CAP:.1%} of cap")
    chk("row count == n_constructs x cells_per_pert x contexts",
        sum(c["stored_entries"] > 0 for c in out.get("contexts", {}).values()) == len(manifest["contexts"]),
        f"expected {rows:,} rows")

    print(f"meta      {args.meta}")
    print(f"incumbent {args.incumbent}")
    print(f"bundle    {bundle}  panel {manifest.get('panel_id')}  contexts {manifest['contexts']}")
    print(f"operator  {op['sha256']}  columns {op['columns']}")
    print(f"whitening after_all_median_norm {meta['effects']['whitening']['after_all_median_norm']} "
          f"(A/B/C incumbent {inc['effects']['whitening']['after_all_median_norm']}; "
          f"reported, NOT asserted — it is a panel-dependent quantity)")
    print()
    failed = [c for c in checks if not c[1]]
    for name, okc, value in checks:
        print(f"  {'PASS' if okc else 'FAIL'}  {name}: {value}")
    print(f"\n{len(checks) - len(failed)}/{len(checks)} checks pass")
    if any(c[0].startswith("arm ") or "incumbent" in c[0] for c in failed):
        print("\n  *** THIS PACKAGE DOES NOT CARRY THE CONFIGURATION THE RUNBOOK SHIPS. ***\n"
              f"  The shipped one is {SHIPPED_BOARD}:\n"
              + "".join(f"    {k:<18} {want!r}   [{env}, {delta}]\n"
                        for k, (want, env, delta, _) in SHIPPED_ARMS.items())
              + "  Rebuild with the driver's defaults (scripts/def_pipeline.py), or, if a drift is\n"
              "  deliberate, name it in --allow-arm-drift and justify it in the round's RESULTS.md.\n"
              "  Do NOT submit this package.")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(
            {"meta": str(args.meta), "checks": [{"check": c, "pass": p, "value": v}
                                                for c, p, v in checks]}, indent=2) + "\n")
    if failed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
