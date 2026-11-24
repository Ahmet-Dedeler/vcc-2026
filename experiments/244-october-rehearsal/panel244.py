"""244 step 1 -- the panel and its census, independently of the driver.

Which bundle, why, and the properties AGENTS.md 8.1g requires stated with every synthetic-panel
number: size, uniqueness, overlap with the real A/B/C panel, overlap with the common-essential
list, coverage by the three production donors, coverage by every registered union source, and
the operator axis digest this bundle implies (235's construction, via the driver's own
`bundle_axis`).  Known-answer rows are asserted against 238 (bundle238.json, census238),
235 (carry235.json coverage on the bundles it measured) and 247 (skip247: ef-bundle axis).

    uv run python experiments/244-october-rehearsal/panel244.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
import def_pipeline as P  # noqa: E402
import make_submission_pm_cd4 as M  # noqa: E402

BUNDLES = {
    "abc": ROOT / "data/challenge",
    "def199_pessimistic": ROOT / "experiments/199-def-runbook/work/def-bundle",
    "def199_realistic": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
    "ef238": ROOT / "experiments/238-essential-free-panel/work/ef-bundle",
}
ESSENTIAL = ROOT / "data/donors/k562_essential/perts.txt"
PRODUCTION = ("k562_gwps", "orion_hct116", "orion_hek293t")


def lines(p: Path) -> set[str]:
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def panel_of(b: Path) -> list[str]:
    return [str(t) for t in pd.read_csv(b / "pert_counts.csv")["target_gene"]]


def main() -> None:
    checks = []

    def chk(name, got, want):
        checks.append({"check": name, "got": got, "want": want, "verdict": "PASS" if got == want else "FAIL"})

    ess = lines(ESSENTIAL)
    chk("essential list sha256 (238 0.1)", sha(ESSENTIAL),
        "b1e5279534eb240f0f3487f0814d379e05b12799b21f6cb3818f5c7f55557e2d")
    chk("essential list size", len(ess), 2021)
    abc = panel_of(BUNDLES["abc"])
    out = {"bundles": {}}
    donors = {d: lines(ROOT / f"data/donors/{d}/perts.txt") for d in PRODUCTION}
    registry = {n: lines(Path(str(M.UNION_SOURCE_REGISTRY[n]["path"])) / "perts.txt")
                for n in sorted(M.UNION_SOURCE_REGISTRY)}
    carry235 = json.loads((ROOT / "experiments/235-october-hardening/work/carry235.json").read_text())["coverage"]
    for tag, b in BUNDLES.items():
        panel = panel_of(b)
        ps = set(panel)
        manifest = json.loads((b / "manifest.json").read_text())
        axis, ncols = P.bundle_axis(b)
        row = {
            "path": str(b.relative_to(ROOT)),
            "panel_id": manifest.get("panel_id"), "contexts": manifest["contexts"],
            "n": len(panel), "n_unique": len(ps),
            "overlap_abc": len(ps & set(abc)), "overlap_essential": len(ps & ess),
            "production_coverage": {d: len(ps & donors[d]) for d in PRODUCTION},
            "carriers_hist": {k: sum(1 for t in panel if sum(t in donors[d] for d in PRODUCTION) == k)
                              for k in (3, 2, 1, 0)},
            "registry_coverage": {n: len(ps & s) for n, s in registry.items()},
            "axis_sha256": axis, "axis_n_cols": ncols,
            "files_sha256": {f: sha(b / f) for f in ("pert_counts.csv", "gene_names.csv", "manifest.json")},
        }
        out["bundles"][tag] = row
        # KA: 235's coverage table on the bundles it measured
        if tag in carry235:
            for n, want in carry235[tag].items():
                chk(f"KA carry235 coverage {tag}/{n}", row["registry_coverage"][n], want)
    ef = out["bundles"]["ef238"]
    chk("ef-bundle: 300 targets", ef["n"], 300)
    chk("ef-bundle: 300 unique", ef["n_unique"], 300)
    chk("ef-bundle: disjoint from A/B/C", ef["overlap_abc"], 0)
    chk("ef-bundle: essential-free", ef["overlap_essential"], 0)
    chk("ef-bundle: contexts D,E,F", ef["contexts"], ["D", "E", "F"])
    b238 = json.loads((ROOT / "experiments/238-essential-free-panel/work/bundle238.json").read_text())["efp"]
    chk("KA 238 pert_counts.csv sha256", ef["files_sha256"]["pert_counts.csv"], b238["files"]["pert_counts.csv"])
    chk("KA 238 gene_names.csv sha256", ef["files_sha256"]["gene_names.csv"], b238["files"]["gene_names.csv"])
    chk("KA 238 k562_gwps 146", ef["production_coverage"]["k562_gwps"], 146)
    chk("KA 238 orion_hct116 293", ef["production_coverage"]["orion_hct116"], 293)
    chk("KA 238 orion_hek293t 296", ef["production_coverage"]["orion_hek293t"], 296)
    chk("KA 247 ef-bundle axis f060e926...", ef["axis_sha256"][:8], "f060e926")
    chk("KA 235 abc axis edf6ba88...", out["bundles"]["abc"]["axis_sha256"],
        "edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f")
    chk("KA 235 def199 axis c464a2ac...", out["bundles"]["def199_pessimistic"]["axis_sha256"][:8], "c464a2ac")
    chk("KA abc: 0 essential, 269/290/294", (out["bundles"]["abc"]["overlap_essential"],
        tuple(out["bundles"]["abc"]["production_coverage"].values())), (0, (269, 290, 294)))
    chk("KA 199 pessimistic: 30 essential, k562 150", (out["bundles"]["def199_pessimistic"]["overlap_essential"],
        out["bundles"]["def199_pessimistic"]["production_coverage"]["k562_gwps"]), (30, 150))
    out["checks"] = checks
    out["all_pass"] = all(c["verdict"] == "PASS" for c in checks)
    WORK.mkdir(exist_ok=True)
    (WORK / "panel244.json").write_text(json.dumps(out, indent=1, sort_keys=True))
    for tag, r in out["bundles"].items():
        print(f"[{tag}] n={r['n']} unique={r['n_unique']} abc={r['overlap_abc']} ess={r['overlap_essential']} "
              f"prod={r['production_coverage']} carriers={r['carriers_hist']} axis={r['axis_sha256'][:16]} ncols={r['axis_n_cols']}")
        print(f"   registry={r['registry_coverage']}")
    for c in checks:
        print(f"[{c['verdict']}] {c['check']}: got={c['got']!r} want={c['want']!r}")
    print(f"\nPANEL {'PASS' if out['all_pass'] else 'FAIL'} ({len(checks)} checks)")
    raise SystemExit(0 if out["all_pass"] else 1)


if __name__ == "__main__":
    main()
