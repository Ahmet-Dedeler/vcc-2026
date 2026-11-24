"""Donor coverage census for a new 300-target panel, in one command.

Transition item 3 (experiment 199). AGENTS.md §2.1 says the panel is selected for "adequate
representation in existing public datasets", which is why K562-gwps covers 269/300 of the
A/B/C panel and Orion 290–294. That is a reason to EXPECT coverage to hold on D/E/F. It is
not a reason to assume it, so this measures it.

What it reports, and why each number matters to the build:

* **Per-target carrier count over the three production donors** (`k562_gwps`, `orion_hct116`,
  `orion_hek293t`). This, not the marginal coverage, is what changes the model's character:
  `MultiDonorModel.effect` averages the donors that carry the target, and
  `donor_agreement` — the shrinkage key and the ordering key — is **all-zero on any target
  fewer than two donors carry** (`src/vcc2026/multi_donor.py`). A target with 0 carriers
  falls back to the panel-mean direction at the median donor magnitude
  (`fallback="mean"`), which is a real prediction, not a zero row.
* **`panel_support`, predicted exactly, in under a second.** The builder's two pinned
  numbers are pure set intersections:
  `perturbmatch_replacement_targets` = |panel ∩ k562_gwps ∩ perturbmatch| (the
  "both supported" branch of `apply_perturbmatch_k562_replacement`), and
  `cd4_complete_targets` = |panel ∩ cd4_rest ∩ cd4_stim8hr ∩ cd4_stim48hr| (a target's
  `condition_count` increments once per bank that carries it, and "complete" is all three).
  Verified on the A/B/C panel against the frozen table's 248 / 276. So on Oct 22 the
  outcome of the 4-minute, ~19 GiB `--repin-panel-support` effects build is known BEFORE it
  is started, and any surprise is caught for free.
* Every other donor table on disk, for context.

    uv run python scripts/def_census.py                       # current panel, known-answer row
    uv run python scripts/def_census.py --bundle <new bundle> --out work/census.json
    uv run python scripts/def_census.py --external            # + 166's data/external sweep
    uv run python scripts/def_census.py --bundle <new bundle> --cd4-root <re-extracted banks>
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DONOR_ROOT = ROOT / "data/donors"
PRODUCTION_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
PERTURBMATCH = ROOT / "data/external/perturbmatch"
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
# The A/B/C answers, as the known-answer row (AGENTS.md 11.1 rule 10's corollary).
KNOWN_ABC = {"panel_id": "vcc2026-val-1", "k562_gwps": 269,
             "perturbmatch_replacement_targets": 248, "cd4_complete_targets": 276}


def read_perts(path: Path) -> set[str]:
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


def census(bundle: Path, external: bool = False, cd4_root: Path = CD4_ROOT) -> dict:
    t0 = time.time()
    perts_csv = bundle / "pert_counts.csv"
    genes_csv = bundle / "gene_names.csv"
    manifest = json.loads((bundle / "manifest.json").read_text())
    # Both CSVs carry a HEADER ROW even though `vcc prep --genes` documents "headerless".
    panel = [str(x) for x in pd.read_csv(perts_csv)["target_gene"]]
    genes = [str(x) for x in pd.read_csv(genes_csv)["gene_name"]]
    panel_set = set(panel)
    out: dict = {
        "bundle": str(bundle),
        "panel_id": manifest.get("panel_id"),
        "contexts": manifest.get("contexts"),
        "n_panel": len(panel),
        "n_panel_unique": len(panel_set),
        "n_genes": len(genes),
        "panel_targets_on_gene_axis": len(panel_set & set(genes)),
    }

    # --- every donor table on disk ------------------------------------------------------
    tables: dict[str, set[str]] = {}
    donor_rows = []
    for d in sorted(DONOR_ROOT.iterdir()):
        if not (d / "perts.txt").exists():
            continue
        s = read_perts(d / "perts.txt")
        tables[d.name] = s
        donor_rows.append({"donor": d.name, "n_perts": len(s),
                           "panel_covered": len(panel_set & s),
                           "panel_fraction": round(len(panel_set & s) / len(panel_set), 4),
                           "production": d.name in PRODUCTION_DONORS})
    out["donors"] = sorted(donor_rows, key=lambda r: (-r["production"], -r["panel_covered"]))

    # --- the number that changes the model's character ----------------------------------
    prod = [tables.get(n, set()) for n in PRODUCTION_DONORS]
    carriers = {t: sum(t in s for s in prod) for t in panel}
    hist = {str(k): sum(1 for v in carriers.values() if v == k) for k in (3, 2, 1, 0)}
    out["production_donors"] = list(PRODUCTION_DONORS)
    out["carrier_histogram"] = hist
    out["consensus"] = {
        "targets_with_3_donors": hist["3"],
        "targets_with_2_donors": hist["2"],
        "targets_with_ge2_donors": hist["3"] + hist["2"],
        "targets_with_1_donor": hist["1"],
        "targets_with_0_donors": hist["0"],
        "agreement_key_defined_frac": round((hist["3"] + hist["2"]) / len(panel), 4),
        "agreement_key_all_zero_frac": round((hist["1"] + hist["0"]) / len(panel), 4),
        "note": ("donor_agreement is all-zero, and the shrinkage leaves alpha=1, on any "
                 "target fewer than two production donors carry; a 0-carrier target falls "
                 "back to the panel-mean direction (MultiDonorModel fallback='mean')"),
    }
    out["targets_below_2_donors"] = sorted(t for t, v in carriers.items() if v < 2)

    # --- panel_support, predicted without a build ---------------------------------------
    pm = read_perts(PERTURBMATCH / "perts.txt") if (PERTURBMATCH / "perts.txt").exists() else set()
    cd4 = {n: read_perts(cd4_root / n / "perts.txt") for n in CD4_NAMES
           if (cd4_root / n / "perts.txt").exists()}
    out["cd4_root"] = str(cd4_root)
    k562 = tables.get("k562_gwps", set())
    pm_replacement = sum(1 for t in panel if t in k562 and t in pm)
    cd4_complete = (sum(1 for t in panel if all(t in s for s in cd4.values()))
                    if len(cd4) == len(CD4_NAMES) else None)
    out["panel_support_predicted"] = {
        "perturbmatch_replacement_targets": pm_replacement,
        "cd4_complete_targets": cd4_complete,
        "perturbmatch_n_perts": len(pm),
        "cd4_per_bank_coverage": {n: len(panel_set & s) for n, s in cd4.items()},
        "derivation": ("pm = |panel n k562_gwps n perturbmatch| (both-supported branch of "
                       "apply_perturbmatch_k562_replacement); cd4 = |panel n all three CD4 "
                       "banks| (condition_count == len(CD4_NAMES) in 062's prepare_cd4_mean)"),
    }

    frozen = json.loads((ROOT / "scripts/frozen_inputs.json").read_text())
    out["frozen_inputs_panel_support"] = frozen.get("panel_support")
    out["frozen_inputs_panel_id"] = frozen.get("panel_id")

    # --- 166's external sweep, on THIS panel --------------------------------------------
    if external:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "census166", ROOT / "experiments/166-donor-acquisition/census.py")
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        mod.PANEL, mod.GENES = panel_set, set(genes)   # rebind to THIS panel
        rows = []
        for f in sorted((ROOT / "data/external").rglob("*.h5ad")):
            try:
                rows.append(mod.inspect(f))
            except Exception as e:
                rows.append({"file": str(f.relative_to(ROOT)), "error": f"{type(e).__name__}: {e}"})
        out["external_sweep"] = sorted(
            rows, key=lambda r: -(r.get("panel_overlap") or 0))

    out["seconds"] = round(time.time() - t0, 2)
    return out


def report(out: dict, verify: bool) -> bool:
    print(f"panel {out['panel_id']}  contexts {out['contexts']}  "
          f"{out['n_panel']} targets ({out['n_panel_unique']} unique), "
          f"{out['panel_targets_on_gene_axis']}/{out['n_panel_unique']} on the {out['n_genes']}-gene axis")
    print("\ndonor tables on disk")
    for r in out["donors"]:
        flag = "PROD" if r["production"] else "    "
        print(f"  {flag} {r['donor']:<18} {r['panel_covered']:>4}/{out['n_panel_unique']} "
              f"({r['panel_fraction']:6.1%})  of {r['n_perts']:>6} perts")
    c = out["consensus"]
    print(f"\nconsensus depth over {out['production_donors']}")
    for k in ("3", "2", "1", "0"):
        print(f"  {k} covering donors: {out['carrier_histogram'][k]:>4}")
    print(f"  agreement key DEFINED (>=2 donors): {c['targets_with_ge2_donors']} "
          f"({c['agreement_key_defined_frac']:.1%})   ALL-ZERO: "
          f"{c['targets_with_1_donor'] + c['targets_with_0_donors']} "
          f"({c['agreement_key_all_zero_frac']:.1%})")
    ps = out["panel_support_predicted"]
    print(f"\npanel_support predicted without a build:")
    print(f"  perturbmatch_replacement_targets = {ps['perturbmatch_replacement_targets']}")
    print(f"  cd4_complete_targets             = {ps['cd4_complete_targets']}")
    print(f"  frozen_inputs.json currently pins  {out['frozen_inputs_panel_support']} "
          f"for panel {out['frozen_inputs_panel_id']}")
    if "external_sweep" in out:
        print("\ndata/external sweep (166), top 15 by panel overlap")
        for r in out["external_sweep"][:15]:
            print(f"  {r.get('panel_overlap', '-'):>4}  {r.get('size_gb', 0) or 0:>7.2f}GB  {r['file']}")
    ok = True
    if verify:
        print("\n--- known-answer row: the A/B/C panel ---")
        got = {"panel_id": out["panel_id"],
               "k562_gwps": next(r["panel_covered"] for r in out["donors"] if r["donor"] == "k562_gwps"),
               **{k: ps[k] for k in ("perturbmatch_replacement_targets", "cd4_complete_targets")}}
        for k, expect in KNOWN_ABC.items():
            good = got[k] == expect
            ok &= good
            print(f"  {k:<34} {got[k]!s:<22} expected {expect}   {'OK' if good else 'MISMATCH'}")
        print(f"VERIFY {'PASS' if ok else 'FAIL'}")
    print(f"\n{out['seconds']}s")
    return ok


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bundle", type=Path, default=ROOT / "data/challenge")
    p.add_argument("--external", action="store_true", help="also run 166's data/external sweep on this panel")
    p.add_argument("--cd4-root", type=Path, default=CD4_ROOT,
                   help="CD4 compact banks to count against (re-extracted per panel; default = production A/B/C banks)")
    p.add_argument("--verify", action="store_true", help="assert the A/B/C known-answer row")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()
    bundle = args.bundle if args.bundle.is_absolute() else ROOT / args.bundle
    cd4_root = args.cd4_root if args.cd4_root.is_absolute() else ROOT / args.cd4_root
    out = census(bundle, external=args.external, cd4_root=cd4_root)
    ok = report(out, args.verify)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.out}")
    if args.verify and not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
