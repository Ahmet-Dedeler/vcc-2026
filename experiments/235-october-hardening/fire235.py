"""235 known-answer row 5 — THE ASSERTION MUST FIRE, and the pre-change builder must not.

An assertion never observed to fire is not a guard.  This runs the REAL production
`load_production_operator` against the REAL three-donor model and the REAL A/B/C operator
(`data/donors/whitener_union_pruned.npz`, sha `e2f65ecf...`), four ways:

  PASS         the A/B/C panel                       -> loads, prints its axis digest
  FIRE         `DEF-RUNBOOK` section 5's pessimistic synthetic D/E/F panel
               (`experiments/199-def-runbook/work/def-bundle`, 300 targets, ZERO overlap with
               A/B/C)                                 -> must raise
  FIRE_BADAXIS 199's `def-bundle-badaxis` (same 300 targets, two symbols TRANSPOSED in
               `gene_names.csv`), whose integer column digest is IDENTICAL to `def-bundle`'s
                                                      -> must raise
  FIRE_NOAXIS  `data/donors/whitener_union_full.npz`, a 45,531-row operator that is not in
               `LEGACY_OPERATOR_AXIS_SHA256` and whose json predates `axis_sha256`
                                                      -> must raise (unvalidatable, refused)

  BLIND        the PRE-235 builder (byte copy at `work/make_submission_pm_cd4.PRE235.py`)
               on the SAME disjoint panel               -> must SUCCEED, silently

The last row is the point: it measures what the change bought.  Everything the old builder
checks -- `n_cols`, the basis shape, the eigenvector shape, the file sha256 -- passes on a panel
whose axis pairs 98.8% of the coordinates with the wrong gene.

Serial, one process, ~12-15 GiB peak (the three-donor model is 8-16 GiB in the production build
itself).  Run it alone.

    /usr/bin/time -l uv run python experiments/235-october-hardening/fire235.py
"""

from __future__ import annotations

import gc
import importlib.util
import json
import resource
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import make_submission_pm_cd4 as M  # noqa: E402

OPERATOR = ROOT / "data/donors/whitener_union_pruned.npz"
OPERATOR_META = ROOT / "data/donors/whitener_union_pruned.json"
UNREGISTERED = ROOT / "data/donors/whitener_union_full.npz"
UNREGISTERED_META = ROOT / "data/donors/whitener_union_full.json"
ABC = ROOT / "data/challenge"
DEF_PESS = ROOT / "experiments/199-def-runbook/work/def-bundle"
DEF_BAD = ROOT / "experiments/199-def-runbook/work/def-bundle-badaxis"
PRE = WORK / "make_submission_pm_cd4.PRE235.py"


def peak_gib() -> float:
    rss = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    return rss / 1024**3 if sys.platform == "darwin" else rss / 1024**2


def load_builder(path: Path, name: str):
    """Import a builder byte copy with `ROOT` still resolving to the repository root."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / path.name)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def panel_of(bundle: Path) -> tuple[str, ...]:
    return tuple(
        map(str, pd.read_csv(bundle / "pert_counts.csv")["target_gene"].to_numpy())
    )


def genes_of(bundle: Path) -> np.ndarray:
    return pd.read_csv(bundle / "gene_names.csv")["gene_name"].to_numpy()


def build_model(module, genes: np.ndarray):
    donors = {
        name: module.Donor.load(ROOT / f"data/donors/{name}")
        for name in ("k562_gwps", "orion_hct116", "orion_hek293t")
    }
    return module.MultiDonorModel(donors, genes, gene_norm=False, eb=False)


def attempt(label: str, fn, *, must_raise: bool, results: dict) -> None:
    t0 = time.time()
    row: dict[str, object] = {"must_raise": must_raise}
    try:
        value = fn()
        row["raised"] = False
        row["value"] = value
    except Exception as exc:  # noqa: BLE001 -- the raise IS the measurement
        row["raised"] = True
        row["exception_type"] = type(exc).__name__
        row["message"] = str(exc)
        row["traceback_tail"] = traceback.format_exc().strip().splitlines()[-1]
    row["seconds"] = round(time.time() - t0, 2)
    row["verdict"] = "PASS" if row["raised"] == must_raise else "FAIL"
    results[label] = row
    print(f"[{label}] {row['verdict']}  raised={row['raised']}  {row['seconds']}s", flush=True)
    if row["raised"]:
        print("--- verbatim ---")
        print(row["message"])
        print("--- end ---", flush=True)


def main() -> None:
    WORK.mkdir(parents=True, exist_ok=True)
    results: dict[str, object] = {
        "operator": str(OPERATOR.relative_to(ROOT)),
        "operator_sha256": M.file_sha256(OPERATOR),
        "builder_sha256": M.file_sha256(ROOT / "scripts/make_submission_pm_cd4.py"),
        "pre_builder_sha256": M.file_sha256(PRE),
    }
    abc_genes = genes_of(ABC)
    abc_panel, def_panel = panel_of(ABC), panel_of(DEF_PESS)
    results["panel_overlap_abc_vs_def"] = len(set(abc_panel) & set(def_panel))

    print(f"[setup] building the three-donor model on the A/B/C axis "
          f"({len(abc_genes)} genes)", flush=True)
    model = build_model(M, abc_genes)
    print(f"[setup] model ready, peak {peak_gib():.2f} GiB", flush=True)

    def _load(panel, operator=OPERATOR, meta=OPERATOR_META, wid="union_pruned", sha="auto"):
        expected = M.WHITENER_REGISTRY["union_pruned"][0] if sha == "auto" else sha
        op = M.load_production_operator(
            model, panel, operator, meta, expected_sha256=expected, whitener_id=wid
        )
        out = {
            "columns": int(len(op.columns)),
            "axis_sha256": op.axis_sha256,
            "axis_sha256_source": op.axis_sha256_source,
        }
        del op
        gc.collect()
        return out

    attempt("PASS_abc", lambda: _load(abc_panel), must_raise=False, results=results)
    attempt("FIRE_def_pessimistic", lambda: _load(def_panel), must_raise=True, results=results)
    attempt(
        "FIRE_unregistered_operator",
        lambda: _load(abc_panel, UNREGISTERED, UNREGISTERED_META, None, None),
        must_raise=True,
        results=results,
    )

    # --- what the PRE-235 builder does with the same inputs -------------------------------
    pre = load_builder(PRE, "pre235_builder")
    results["pre_builder_has_axis_guard"] = hasattr(pre, "operator_axis_sha256")

    def _load_pre(panel):
        op = pre.load_production_operator(
            model,
            panel,
            OPERATOR,
            OPERATOR_META,
            expected_sha256=pre.WHITENER_REGISTRY["union_pruned"][0],
        )
        out = {
            "columns": int(len(op.columns)),
            "metadata_n_cols": op.metadata.get("n_cols"),
            "evecs_shape": list(op.whitener.evecs.shape),
            "note": "the pre-235 builder accepted the WRONG panel with no error",
        }
        del op
        gc.collect()
        return out

    attempt("BLIND_pre235_def_pessimistic", lambda: _load_pre(def_panel),
            must_raise=False, results=results)

    del model
    gc.collect()

    # --- the transposed gene axis: identical column indices, different symbols -------------
    bad_genes = genes_of(DEF_BAD)
    print(f"[setup] rebuilding the model on def-bundle-badaxis's gene axis, "
          f"peak {peak_gib():.2f} GiB", flush=True)
    model = build_model(M, bad_genes)
    attempt("FIRE_badaxis", lambda: _load(panel_of(DEF_BAD)), must_raise=True, results=results)

    results["peak_rss_gib"] = round(peak_gib(), 3)
    rows = [k for k in results if isinstance(results[k], dict) and "verdict" in results[k]]
    results["all_pass"] = all(results[k]["verdict"] == "PASS" for k in rows)
    (WORK / "fire235.json").write_text(json.dumps(results, indent=1, sort_keys=True))
    print(json.dumps({k: results[k] for k in rows}, indent=1, sort_keys=True))
    print(f"\nROW 5 {'PASS' if results['all_pass'] else 'FAIL'}  "
          f"peak {results['peak_rss_gib']} GiB")
    raise SystemExit(0 if results["all_pass"] else 1)


if __name__ == "__main__":
    main()
