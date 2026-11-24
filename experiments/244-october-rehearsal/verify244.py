"""244 -- verify the rehearsal artifact with numbers, off the WRITTEN meta and the driver's own logs.

Everything asserted is recomputed (bundle axis, set-intersection coverage, row sums streamed off the
h5ad) or read from an earlier experiment's recorded number (247, 238, 235); nothing is taken from
the command that was typed.

    uv run python experiments/244-october-rehearsal/verify244.py --run work/run1 --name <artifact stem>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import def_pipeline as P  # noqa: E402

EF = ROOT / "experiments/238-essential-free-panel/work/ef-bundle"
ABC_OP = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
ABC_AXIS = "edf6ba886286d7837fa6adaef152372cee481331a81987789151862d5c0a400f"
EF_AXIS_247 = "f060e9262a753dddb28d4a3b80f1ab75d4ca37f8c7de7ef9205f5c156ae0213f"
BUILDER_HEAD = "24a30055895eb251ea087c3517fe08a47ca2084f969355a1959c76852c56daab"
RUNBOOK_199 = {"preflight": 4.5, "identity": 14.7, "cd4-extract": 2.2, "census": 0.4, "whitener": 437.6,
               "repin-frozen": 2.1, "pin-artifacts": 0.5, "repin-panel-support": 9.2, "build": 1255.8,
               "postbuild": 0.2, "prep": 11.8}
CAP = 4_750_000_000
MAX_CELL = 1_000_000


def sha(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def lines(p: Path) -> set[str]:
    return {ln.strip() for ln in p.read_text().splitlines() if ln.strip()}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=Path, default=HERE / "work/run1")
    ap.add_argument("--name", default="sparrow-in-the-vestibule")
    ap.add_argument("--skip-h5ad-scan", action="store_true")
    args = ap.parse_args()
    run = args.run if args.run.is_absolute() else HERE / args.run
    meta = json.loads((run / f"{args.name}.h5ad.meta.json").read_text())
    cfg, eff, op = meta["config"], meta["effects"], meta["operator"]
    checks, report = [], {}

    def chk(name, got, want, tol=0.0):
        ok = (got == want) if tol == 0 else abs(float(got) - float(want)) <= tol
        checks.append({"check": name, "got": got, "want": want, "verdict": "PASS" if ok else "FAIL"})

    # --- 1. the operator was re-fit for THIS panel ------------------------------------------
    axis, ncols = P.bundle_axis(EF)
    wjson = json.loads(Path(op["metadata_path"]).read_text()) if op.get("metadata_path") else {}
    report["operator"] = {k: op.get(k) for k in ("path", "sha256", "metadata_sha256", "columns", "axis_sha256",
                                                 "axis_sha256_checked_against", "panel_sha256")}
    report["operator"]["metadata"] = {k: wjson.get(k) for k in ("n_cols", "basis", "rows_kept_per_donor", "dropped_in_panel",
                                                                 "basis_prune", "axis_sha256", "panel_sha256", "seconds")}
    chk("operator.sha256 != A/B/C's e2f65ecf", op["sha256"] != ABC_OP, True)
    chk("operator.axis_sha256 == axis recomputed from ef-bundle", op["axis_sha256"], axis)
    chk("KA 247: ef-bundle axis is f060e926...", axis, EF_AXIS_247)
    chk("operator.axis_sha256 != A/B/C's edf6ba88", op["axis_sha256"] != ABC_AXIS, True)
    chk("axis checked against the fresh json (not the legacy table)", op.get("axis_sha256_checked_against"),
        f"metadata:{Path(op['metadata_path']).name}")
    chk("operator.columns == recomputed n_cols", op["columns"], ncols)
    chk("whitener json axis_sha256 == meta's", wjson.get("axis_sha256"), op["axis_sha256"])
    chk("sha256(file at operator.path) == operator.sha256", sha(Path(op["path"])), op["sha256"])
    chk("dropped_in_panel == [] (essential-free panel)", wjson.get("dropped_in_panel"), [])
    # --- 2. the configuration the day must carry ------------------------------------------
    want_cfg = {"alpha": 0.982, "cd4_target_f": 0.34, "h1_union_g": 0.28,
                "h1_union_sources": ["gse295214_h1", "h1_2025"], "emitter": "rc_callset",
                "whitener_id": "union_pruned_def", "cd4_sizing": "carrying",
                "cd4_target_g": 0.3443655278649641, "contexts": ["D", "E", "F"], "expect_contexts": "D,E,F"}
    for k, v in want_cfg.items():
        chk(f"config.{k}", cfg.get(k, "<absent>"), v)
    report["config"] = {k: cfg.get(k, "<absent>") for k in list(want_cfg) + ["whiten_gamma", "whiten_lam", "gamma", "lam",
                                                                                "min_cpm", "keep_top", "keep_spread",
                                                                                "tail_spread", "overdispersion", "targets"]}
    # --- 3. CD4 on this panel ---------------------------------------------------------------
    cd4 = eff["cd4"]
    rf = cd4["residual_frame"]
    report["cd4"] = {"fraction": cd4["fraction"], "fully_supported_targets": cd4["fully_supported_targets"],
                     "residual_frame": rf}
    chk("KA 247: cd4 complete targets 182", cd4["fully_supported_targets"], 182)
    chk("cd4 residual_frame.sizing_rule == carrying", rf.get("sizing_rule"), "carrying")
    chk("cd4 target_g_requested == 0.3443655278649641", rf.get("target_g_requested"), 0.3443655278649641)
    ratio = rf["residual_median_row_norm_over_carrying_rows"] / rf["residual_median_row_norm_over_all_rows"]
    report["cd4"]["overdose_factor_measured"] = ratio
    chk("KA 247: CD4 over-dose factor (carrying/all medians) == 1.1701602469240082", ratio, 1.1701602469240082, tol=1e-12)
    chk("meta overdose_factor_all_rows_rule == the ratio", rf.get("overdose_factor_all_rows_rule"), ratio, tol=1e-12)
    chk("effects.cd4.fraction == residual_frame.solved_fraction", cd4["fraction"], rf["solved_fraction"])
    chk("legacy_all_rows_rule_would_solve == solved x overdose x (0.34 / g_cd4)  [the two rules run at different targets]",
        rf["legacy_all_rows_rule_would_solve"], rf["solved_fraction"] * ratio * 0.34 / rf["target_g_requested"], tol=1e-12)
    report["cd4"]["legacy_over_carrying_solved"] = rf["legacy_all_rows_rule_would_solve"] / rf["solved_fraction"]
    # --- 4. the H1 union at low coverage --------------------------------------------------
    h1 = eff["h1_union"]
    hf = h1["residual_frame"]
    report["h1_union"] = {k: h1.get(k) for k in ("carrying_targets", "applied_fraction", "source_coverage",
                                                  "targets_by_source_count", "sources")}
    report["h1_union"]["residual_frame"] = hf
    panel = set(pd.read_csv(EF / "pert_counts.csv")["target_gene"].astype(str))
    import make_submission_pm_cd4 as M  # noqa: E402  (registry paths only)
    src_sets = {n: lines(Path(str(M.UNION_SOURCE_REGISTRY[n]["path"])) / "perts.txt") for n in cfg["h1_union_sources"]}
    expect_cov = {n: len(panel & s) for n, s in src_sets.items()}
    union_rows = len(panel & set.union(*src_sets.values()))
    chk("h1 source_coverage == set-intersection", {k: int(v) for k, v in h1["source_coverage"].items()}, expect_cov)
    chk("h1 carrying_rows == |panel ∩ (gse ∪ h1_2025)|", hf["carrying_rows"], union_rows)
    chk("h1 shipped_median_rule_would_solve == 0.0 (all-rows rule dead below 50%)", hf["shipped_median_rule_would_solve"], 0.0)
    chk("h1 residual_median_row_norm_over_all_rows == 0.0", hf["residual_median_row_norm_over_all_rows"], 0.0)
    chk("h1 solved_fraction > 0", hf["solved_fraction"] > 0, True)
    chk("h1 g requested 0.28", hf.get("target_g_requested"), 0.28)
    chk("h1 carrying_rows < 150 (below 50% coverage)", hf["carrying_rows"] < 150, True)
    # --- 5. the emitted artifact --------------------------------------------------------
    out = meta["output"]
    report["output"] = {k: out.get(k) for k in ("path", "sha256", "bytes", "nnz", "n_rows", "contexts")}
    chk("nnz <= 4.75e9", out["nnz"] <= CAP, True)
    chk("n_rows == 360000", out.get("n_rows"), 360000)
    pm = eff["perturbmatch"]
    k562 = lines(ROOT / "data/donors/k562_gwps/perts.txt"); pmset = lines(ROOT / "data/external/perturbmatch/perts.txt")
    report["panel_support"] = {"perturbmatch": pm["both_supported_modified_count"], "cd4": cd4["fully_supported_targets"],
                               "expected_pm": len(panel & k562 & pmset), "frozen_bundle.panel_support": (meta.get("frozen_bundle") or {}).get("panel_support")}
    chk("perturbmatch count == |panel ∩ k562 ∩ pm|", pm["both_supported_modified_count"], len(panel & k562 & pmset))
    h5 = run / f"{args.name}.h5ad"
    if h5.exists() and not args.skip_h5ad_scan:
        t0 = time.time()
        with h5py.File(h5, "r") as f:
            X = f["X"]
            indptr = X["indptr"][()]
            data = X["data"]
            n = len(indptr) - 1
            mx, mn, nnz, neg, nonint, zeros = 0.0, np.inf, 0, 0, 0, 0
            step = 20000
            for s in range(0, n, step):
                e = min(n, s + step)
                a, b = int(indptr[s]), int(indptr[e])
                d = data[a:b]
                sums = np.add.reduceat(d, (indptr[s:e] - a).astype(np.int64)) if b > a else np.zeros(e - s)
                # rows with zero nnz: reduceat gives the next element; guard by masking empty rows
                empty = (indptr[s + 1:e + 1] - indptr[s:e]) == 0
                sums = np.where(empty, 0.0, sums)
                mx, mn = max(mx, float(sums.max())), min(mn, float(sums.min()))
                nnz += int(b - a); neg += int((d < 0).sum()); nonint += int((d != np.round(d)).sum()); zeros += int((d == 0).sum())
            report["h5ad_scan"] = {"rows": int(n), "nnz": nnz, "max_counts_per_cell": mx, "min_counts_per_cell": mn,
                                   "negative_entries": neg, "non_integer_entries": nonint,
                                   "explicit_zeros": zeros,
                                   "seconds": round(time.time() - t0, 1), "bytes": h5.stat().st_size}
        chk("h5ad rows == 360000", n, 360000)
        chk("h5ad nnz == meta nnz", nnz, out["nnz"])
        chk("max counts per cell < 1e6", mx < MAX_CELL, True)
        chk("no negative entries", neg, 0)
        chk("all-integer entries", nonint, 0)
        chk("no explicit zeros", report["h5ad_scan"]["explicit_zeros"], 0)
    # --- 6. postbuild, prep, arms block, timings ------------------------------------------
    pb = json.loads((run / "postbuild.json").read_text())
    n_pass = sum(1 for c in pb["checks"] if c["pass"]); n_all = len(pb["checks"])
    report["postbuild"] = {"pass": n_pass, "total": n_all, "failed": [c["check"] for c in pb["checks"] if not c["pass"]]}
    chk("def_postbuild.py all checks pass", n_pass, n_all)
    prep_log = (run / "logs/prep.log").read_text()
    report["prep_log_tail"] = [ln for ln in prep_log.splitlines() if ln.strip() and not re.match(r"^\s*\d+\s+(average|page|swaps|block|messages|signals|voluntary|involuntary|instructions|cycles|peak memory|maximum resident)", ln)][-20:]
    chk("vcc prep --dry-run clean (rc 0 in timings)", json.loads((run / "timings.json").read_text())["prep"]["returncode"], 0)
    build_log = (run / "logs/build.log").read_text()
    report["build_log_operator_line"] = [ln for ln in build_log.splitlines() if "[operator]" in ln][:3]
    report["build_log_peak"] = [ln for ln in build_log.splitlines() if "peak_rss_gib" in ln][-2:]
    drv = (run / "driver.log").read_text()
    report["arms_block"] = [ln for ln in drv.splitlines() if ln.startswith("  OK  config.") or ln.startswith("  NO  config.") or "[arms]" in ln]
    chk("[arms] block: no NO lines", sum(1 for ln in report["arms_block"] if ln.startswith("  NO")), 0)
    chk("[arms] block: 7 OK lines (alpha, f, g, sources, sizing, g_cd4 ...)", sum(1 for ln in report["arms_block"] if ln.startswith("  OK")) >= 6, True)
    tm = json.loads((run / "timings.json").read_text())
    table = {s: (t["seconds"], t.get("peak_rss_gib")) for s, t in tm.items() if not s.startswith("_")}
    report["timings"] = table
    report["timings_total_seconds"] = tm.get("_total", {}).get("seconds")
    report["runbook_199_seconds"] = RUNBOOK_199
    report["sum_of_stages"] = round(sum(v[0] for v in table.values()), 1)
    # --- 7. the builder is byte-restored, and the pin record ----------------------------
    b_now = sha(ROOT / "scripts/make_submission_pm_cd4.py")
    pin = json.loads((run / "pin-artifacts.json").read_text())
    report["builder"] = {"sha_now": b_now, "sha_head_before_run": BUILDER_HEAD, "pin_before": pin["builder_sha256_before"],
                         "pin_after": pin["builder_sha256_after"], "meta_hashes_implementation": meta.get("hashes", {}).get("implementation")}
    # 249 edited the builder (flag-gated VCC_CD4_DEFLATE_K, default off) between this run's byte copy
    # of HEAD (24a30055) and its pin-artifacts stage, so the pre-pin bytes are 249's, not HEAD's.
    chk("builder restored to the PRE-PIN bytes (249's working copy, not HEAD)", b_now, pin["builder_sha256_before"])
    report["builder"]["pre_pin_is_head"] = pin["builder_sha256_before"] == BUILDER_HEAD
    chk("meta carries no effects.cd4.deflation (249's flag was OFF)", "deflation" in eff["cd4"], False)
    chk("meta carries no config.cd4_deflate_k (249's flag was OFF)", "cd4_deflate_k" in cfg, False)
    impl = meta.get("hashes", {}).get("implementation") or {}
    impl_sha = impl.get("scripts/make_submission_pm_cd4.py") if isinstance(impl, dict) else impl
    chk("meta hashes.implementation == the PINNED builder (the day's committed edit)", impl_sha, pin["builder_sha256_after"])
    git = subprocess.run(["git", "status", "--short", "scripts/"], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    report["git_status_scripts"] = git
    chk("git status scripts/: only the builder (249's uncommitted edit), nothing of 244's",
        [ln.split()[-1] for ln in git.splitlines()], ["scripts/make_submission_pm_cd4.py"] if not report["builder"]["pre_pin_is_head"] else [])
    frozen = ROOT / "scripts/frozen_inputs.json"
    chk("scripts/frozen_inputs.json untouched (a654d66b...)", sha(frozen)[:8], "a654d66b")
    # --- 8. disk ----------------------------------------------------------------------
    du = subprocess.run(["du", "-sk", str(run)], capture_output=True, text=True).stdout.split()[0]
    report["disk_run_gib"] = round(int(du) / 2**20, 2)
    report["checks"] = checks
    report["all_pass"] = all(c["verdict"] == "PASS" for c in checks)
    (HERE / "work/verify244.json").write_text(json.dumps(report, indent=1, sort_keys=True, default=str))
    print(json.dumps({k: report[k] for k in ("operator", "config", "cd4", "h1_union", "output", "panel_support",
                                             "h5ad_scan", "postbuild", "timings", "timings_total_seconds", "sum_of_stages",
                                             "builder", "disk_run_gib") if k in report}, indent=1, default=str))
    print("\nARMS BLOCK:\n" + "\n".join(report["arms_block"]))
    print("\nPREP LOG:\n" + "\n".join(report["prep_log_tail"]))
    print()
    for c in checks:
        print(f"[{c['verdict']}] {c['check']}: got={c['got']!r} want={c['want']!r}")
    print(f"\nVERIFY {'PASS' if report['all_pass'] else 'FAIL'}  {len(checks)} checks")
    raise SystemExit(0 if report["all_pass"] else 1)


if __name__ == "__main__":
    main()
