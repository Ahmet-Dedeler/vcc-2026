"""244 -- CD4's over-dose factor on the rehearsal panel, measured from the STAGED residual.

Both sizing rules depend on the whitened base only through its median row norm, so 247's
construction (stub base with the shipped median) gives the base-INDEPENDENT quantity the mandate
asks for: residual_median_row_norm_over_carrying_rows / ...over_all_rows.  Known-answer rows:
A/B/C 1.0128397878381297 (247), 199's pessimistic 1.196293034774599 (247), ef-bundle on 238's
banks 1.1701602469240082 (247).  Then the same on the banks THIS rehearsal re-extracted
(work/run1/cd4_compact), which must equal 238's byte-for-byte on 15/15 files, and finally the
factor read off the built meta's residual_frame, which must equal it.

    uv run python experiments/244-october-rehearsal/cd4dose244.py [--meta <built meta>]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/247-oct22-guards"))
import cd4dose247 as D  # noqa: E402  (247's measure(); its main() is never called)

M = D.M
KA = {
    "abc": ("data/challenge", "data/external/cd4t/cloud-compact/cd4_compact", 1.0128397878381297, 276),
    "def199_pessimistic": ("experiments/199-def-runbook/work/def-bundle",
                           "experiments/199-def-runbook/work/cd4_compact_def199", 1.196293034774599, 167),
    "ef238_on_238_banks": ("experiments/238-essential-free-panel/work/ef-bundle",
                           "experiments/238-essential-free-panel/work/cd4/efp/cd4_compact", 1.1701602469240082, 182),
}
RUN1_BANKS = HERE / "work/run1/cd4_compact"
FILES = ("delta.npy", "b_ctrl.npy", "genes.txt", "perts.txt", "meta.json")


def sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--meta", type=Path, default=None)
    args = ap.parse_args()
    checks, rows = [], {}

    def chk(name, got, want, tol=0.0):
        ok = (got == want) if tol == 0 else abs(float(got) - float(want)) <= tol
        checks.append({"check": name, "got": got, "want": want, "verdict": "PASS" if ok else "FAIL"})

    for tag, (bundle, banks, want_ratio, want_complete) in KA.items():
        r = D.measure(tag, ROOT / bundle, ROOT / banks)
        rows[tag] = {k: v for k, v in r.items() if not k.startswith("_")}
        chk(f"KA 247 {tag}: overdose factor", r["overdose_factor_all_rows_rule"], want_ratio)
        chk(f"KA 247 {tag}: complete targets", r["fully_supported_targets"], want_complete)
        print(f"[{tag}] complete {r['fully_supported_targets']} carrying {r['carrying_rows']} "
              f"median_all {r['median_all_rows']!r} median_carrying {r['median_carrying_rows']!r} "
              f"overdose {r['overdose_factor_all_rows_rule']!r}")
    if RUN1_BANKS.exists():
        r = D.measure("ef238_on_run1_banks", ROOT / KA["ef238_on_238_banks"][0], RUN1_BANKS)
        rows["ef238_on_run1_banks"] = {k: v for k, v in r.items() if not k.startswith("_")}
        chk("run1 banks: overdose factor == 238 banks'", r["overdose_factor_all_rows_rule"],
            rows["ef238_on_238_banks"]["overdose_factor_all_rows_rule"])
        same = {f"{b}/{f}": sha(RUN1_BANKS / b / f) == sha(ROOT / KA["ef238_on_238_banks"][1] / b / f)
                for b in M.CD4_NAMES for f in FILES}
        data_same = {k: v for k, v in same.items() if not k.endswith("meta.json")}
        chk("run1 banks byte-identical to 238's efp banks (12 DATA files; meta.json is provenance)", sum(data_same.values()), 12)
        rows["run1_banks_meta_json_identical"] = {k: v for k, v in same.items() if k.endswith("meta.json")}
        rows["run1_banks_identical"] = same
        # the legacy rule at f=0.34 and the carrying rule at g_cd4 on the real residual, stub base
        g_cd4 = 0.3443655278649641
        cf, cfr = M.solve_carrying_row_fraction(r["_base"], r["_cd4"], g_cd4)
        rows["run1_stub_base_rules"] = {
            "legacy_solved_at_f034": r["legacy_rule"].get("solved_fraction"),
            "carrying_solved_at_g_cd4": cf,
            "carrying_frame": cfr,
            "legacy_over_carrying": (r["legacy_rule"].get("solved_fraction") / cf) if cf else None,
        }
        print(f"[run1 banks, stub base] legacy f=0.34 -> {r['legacy_rule'].get('solved_fraction')!r}; "
              f"carrying g={g_cd4} -> {cf!r}; ratio {rows['run1_stub_base_rules']['legacy_over_carrying']!r}")
    if args.meta and args.meta.exists():
        meta = json.loads(args.meta.read_text())
        rf = meta["effects"]["cd4"]["residual_frame"]
        rows["meta_residual_frame"] = rf
        ratio_meta = rf["residual_median_row_norm_over_carrying_rows"] / rf["residual_median_row_norm_over_all_rows"]
        chk("meta: carrying/all median ratio == staged overdose factor", ratio_meta,
            rows["ef238_on_238_banks"]["overdose_factor_all_rows_rule"], tol=1e-12)
        chk("meta: overdose_factor_all_rows_rule field == the ratio", rf.get("overdose_factor_all_rows_rule"),
            ratio_meta, tol=1e-12)
        chk("meta: sizing_rule carrying", rf.get("sizing_rule"), "carrying")
        chk("meta: target_g_requested", rf.get("target_g_requested"), 0.3443655278649641)
        chk("meta: carrying_rows", rf.get("carrying_rows"), 182)
        # the legacy rule solves at f = 0.34, the carrying rule at g_cd4 = 0.34 x 1.01284, so
        # legacy/carrying = overdose x (0.34 / g_cd4)  (247 section 3.3: "legacy delivered x1.1553 of A/B/C")
        chk("meta: legacy_would_solve / solved == overdose x (0.34 / g_cd4)",
            rf["legacy_all_rows_rule_would_solve"] / rf["solved_fraction"],
            ratio_meta * 0.34 / rf["target_g_requested"], tol=1e-12)
        rows["meta_legacy_over_carrying"] = rf["legacy_all_rows_rule_would_solve"] / rf["solved_fraction"]
        print(f"[meta] base median {rf['whitened_base_median_row_norm']!r} carrying {rf['carrying_rows']} "
              f"median_carrying {rf['residual_median_row_norm_over_carrying_rows']!r} "
              f"median_all {rf['residual_median_row_norm_over_all_rows']!r} ratio {ratio_meta!r} "
              f"solved {rf['solved_fraction']!r} legacy_would_solve {rf['legacy_all_rows_rule_would_solve']!r} "
              f"effects.cd4.fraction {meta['effects']['cd4']['fraction']!r}")
    out = {"rows": rows, "checks": checks, "peak_rss_gib": round(D.peak_gib(), 3),
           "all_pass": all(c["verdict"] == "PASS" for c in checks)}
    (HERE / "work/cd4dose244.json").write_text(json.dumps(out, indent=1, sort_keys=True, default=str))
    for c in checks:
        print(f"[{c['verdict']}] {c['check']}: got={c['got']!r} want={c['want']!r}")
    print(f"\nCD4 DOSE {'PASS' if out['all_pass'] else 'FAIL'}  {len(checks)} checks, peak {out['peak_rss_gib']} GiB")
    raise SystemExit(0 if out["all_pass"] else 1)


if __name__ == "__main__":
    main()
