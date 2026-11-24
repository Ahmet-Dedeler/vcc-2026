"""233 Phase 1.1 — the Oct-22 operator trap, demonstrated deliberately so the assertion that
catches it is PROVEN to fire, and priced on the six members.

AGENTS.md §8.3 / DEF-RUNBOOK §0.1.  The whitening operator's gene axis is
`union(donor coverage) - panel`.  Point the A/B/C operator at a new panel and every check the
production builder runs PASSES while the coordinates are paired with the wrong genes.  212
priced that at -0.0429 overall on `pds` alone.  This script does four things:

  1. NEGATIVE CONTROL — run every production-side check (file sha256, `donors`, `n_cols`,
     eigenvector shape) with the A/B/C operator against BOTH synthetic panels, and show they
     all pass.  Report the fraction of axis positions that are paired with a DIFFERENT gene.
  2. POSITIVE CONTROL — the same operator against its OWN panel: the axis hash matches, so
     the assertion is not one that always fires.
  3. THE ASSERTION — `panels233.load_frozen_operator_233` refuses the mismatched pair, with
     the axis hash in the message.
  4. THE PRICE — one PLOLO fold on the pessimistic panel, scored twice: once with the
     correctly re-fitted operator and once with the A/B/C one, all six members.

    uv run python trap233.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import panels233 as P  # noqa: E402
import plolo  # noqa: E402
from plolo import E062, OFFICIAL  # noqa: E402

PROD_OP = ROOT / "data/donors/whitener_union_pruned.npz"
PROD_SHA = "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123"
PROD_DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")


def wcols_for(panel: str, donors) -> np.ndarray:
    banks = [E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in donors]
    union: set[str] = set()
    for b in banks:
        union |= set(b.genes)
    ps = set(P.read_panel(panel))
    return np.asarray([i for i, g in enumerate(OFFICIAL) if g in union and g not in ps],
                      dtype=np.int64)


def production_checks(op_path: Path, meta: dict, wcols: np.ndarray) -> dict:
    """Exactly what the production builder verifies about a frozen operator."""
    with np.load(op_path) as z:
        shape = tuple(z["evecs"].shape)
        n_evals = int(z["evals"].shape[0])
    return dict(
        file_sha256_matches_pin=E062.sha256_file(op_path) == PROD_SHA,
        donors_match=list(meta.get("donors", [])) == list(PROD_DONORS),
        n_cols_matches_recomputed_axis=int(meta["n_cols"]) == int(len(wcols)),
        eigenvector_shape_ok=shape == (len(wcols), len(wcols)),
        n_evals_ok=n_evals == len(wcols),
        dropped_in_panel=meta.get("dropped_in_panel"),
    )


def main() -> None:
    meta = json.loads(PROD_OP.with_suffix(".json").read_text())
    out: dict[str, object] = {"operator": str(PROD_OP.relative_to(ROOT)),
                              "n_cols_recorded": int(meta["n_cols"])}
    w_abc = wcols_for("abc", PROD_DONORS)
    print(f"A/B/C production operator: n_cols {meta['n_cols']}, "
          f"recomputed A/B/C axis {len(w_abc)}\n")

    for panel in ("abc", "real", "pess"):
        w = wcols_for(panel, PROD_DONORS)
        checks = production_checks(PROD_OP, meta, w)
        same_len = len(w) == len(w_abc)
        mis = (float(np.mean(w[:len(w_abc)] != w_abc[:len(w)])) if same_len
               else float("nan"))
        sha_p, _ = P.axis_sha256(set(P.read_panel(panel)),
                                 [E062.Bank.load(n, plolo.DONOR_PATHS[n])
                                  for n in PROD_DONORS])
        sha_a, _ = P.axis_sha256(set(P.read_panel("abc")),
                                 [E062.Bank.load(n, plolo.DONOR_PATHS[n])
                                  for n in PROD_DONORS])
        rec = dict(n_cols_this_panel=int(len(w)),
                   production_checks=checks,
                   all_production_checks_pass=all(v is True for k, v in checks.items()
                                                  if isinstance(v, bool)),
                   axis_hash_matches_abc_operator=(sha_p == sha_a),
                   fraction_of_positions_paired_with_a_DIFFERENT_gene=mis,
                   genes_only_in_abc_axis=int(len(set(w_abc.tolist()) - set(w.tolist()))),
                   genes_only_in_this_axis=int(len(set(w.tolist()) - set(w_abc.tolist()))))
        out[panel] = rec
        print(f"[{panel}] n_cols {len(w)}  production checks all pass: "
              f"{rec['all_production_checks_pass']}  axis hash matches the A/B/C operator: "
              f"{rec['axis_hash_matches_abc_operator']}")
        print(f"      positions paired with a DIFFERENT gene: {100 * mis:.2f}%   "
              f"axis symmetric difference {rec['genes_only_in_abc_axis']}"
              f"+{rec['genes_only_in_this_axis']} genes")

    # ---- 3. the assertion, on a 2-donor fold operator that exists on disk
    print("\nTHE ASSERTION (`panels233.load_frozen_operator_233`):")
    for panel, expect in (("abc", "load"), ("pess", "refuse")):
        op = P.op_path("abc", "orion_hct116")
        if not op.exists():
            print("      (abc__held_orion_hct116 not fitted yet — skipped)")
            break
        banks = [E062.Bank.load(n, plolo.DONOR_PATHS[n])
                 for n in P.FOLD_SOURCES["orion_hct116"]]
        try:
            P.load_frozen_operator_233(op, banks, OFFICIAL, set(P.read_panel(panel)),
                                       expected_sha256=E062.sha256_file(op))
            got = "load"
            msg = "loaded"
        except ValueError as e:
            got = "refuse"
            msg = str(e)[:160]
        ok = got == expect
        out[f"assertion|{panel}"] = dict(expected=expect, got=got, message=msg, pass_=ok)
        print(f"   asked for panel {panel:5} -> {got:6} (expected {expect:6}) "
              f"{'PASS' if ok else 'FAIL'}\n      {msg}")

    (P.WORK / "trap233.json").write_text(json.dumps(out, indent=2))
    print("\n-> work/trap233.json")


if __name__ == "__main__":
    main()
