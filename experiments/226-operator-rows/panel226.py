"""226 falsifier 5 — the PANEL-IDENTITY check, and a deliberate demonstration that it fires.

`docs/DEF-RUNBOOK.md` §0.1: the whitening operator's gene axis is `union(donor coverage) - panel`,
the `.npz` stores **no gene names at all** (only `evals` and `evecs`), and pointing an A/B/C
operator at a new 300-gene panel passes every check the builder runs while pairing 98.8% of the
coordinates with the wrong gene.  212 priced that trap at -0.0429 overall on `pds` alone.

Any operator this experiment refits inherits that trap, so the guard has to be proven, not
asserted.  This script:

  1. recomputes the operator axis for the A/B/C bundle          -> the axis the shipped file has;
  2. recomputes it for 199's synthetic D/E/F bundle             -> the axis a new panel produces;
  3. measures the positional mispairing between the two         -> the size of the silent failure;
  4. runs `rows226`'s guard on the mismatched pair              -> it must RAISE;
  5. runs the same guard on the matched pair                    -> it must PASS.

Cheap: donor coverage comes from each donor's `genes.txt`, so no `delta.npy` is ever loaded.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
BUNDLES = {
    "abc (production)": ROOT / "data/challenge",
    "def199 (synthetic new panel)": ROOT / "experiments/199-def-runbook/work/def-bundle",
    "def199r (second synthetic)": ROOT / "experiments/199-def-runbook/work/def-bundle-r",
}


def operator_axis(bundle: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`sorted(union(donor coverage) - panel)` as positions on the bundle's official axis."""
    target_genes = pd.read_csv(bundle / "gene_names.csv")["gene_name"].to_numpy()
    panel = pd.read_csv(bundle / "pert_counts.csv")["target_gene"].to_numpy()
    tindex = {g: i for i, g in enumerate(target_genes)}
    union: set[str] = set()
    for n in DONORS:
        genes = (ROOT / "data/donors" / n / "genes.txt").read_text().split("\n")
        union |= {g for g in genes if g in tindex}
    union -= set(panel)
    wcols = np.array(sorted(tindex[g] for g in union), dtype=np.int64)
    return wcols, target_genes[wcols], panel


def guard(fold_axis: np.ndarray, operator_axis_: np.ndarray) -> None:
    """rows226's guard, verbatim in meaning: the operator's axis must BE the fold's axis."""
    if len(fold_axis) != len(operator_axis_) or not np.array_equal(fold_axis, operator_axis_):
        raise SystemExit(
            f"PANEL-IDENTITY ABORT: fold axis ({len(fold_axis)}) != operator axis "
            f"({len(operator_axis_)}) -- see falsifier 5")


def main() -> None:
    out: dict = {"axes": {}, "pairings": {}, "guard": {}}
    axes = {}
    for name, b in BUNDLES.items():
        if not b.exists():
            print(f"  {name}: MISSING ({b}) -- skipped", flush=True)
            continue
        w, syms, panel = operator_axis(b)
        axes[name] = (w, syms)
        out["axes"][name] = {"n_axis": int(len(w)), "n_panel": int(len(panel)),
                             "bundle": str(b.relative_to(ROOT))}
        print(f"  {name:32s} axis {len(w):6d} genes | panel {len(panel)}", flush=True)

    ref = "abc (production)"
    wref, sref = axes[ref]
    print(f"\nshipped operator was fitted on the {ref} axis: {len(wref)} columns")
    print("the .npz stores evals+evecs and NO gene names, so nothing downstream can "
          "detect a swap by itself.\n")

    for name, (w, syms) in axes.items():
        if name == ref:
            continue
        same_len = len(w) == len(wref)
        k = min(len(w), len(wref))
        # what the operator would actually do: coordinate i of the prediction, gathered on
        # the NEW axis, is multiplied by the eigenvector row fitted for the OLD axis's gene i.
        mispaired = int((sref[:k] != syms[:k]).sum())
        frac = mispaired / k
        out["pairings"][name] = {
            "same_length": bool(same_len), "n_compared": int(k),
            "mispaired": mispaired, "mispaired_frac": round(frac, 6),
            "identical_axis": bool(same_len and np.array_equal(w, wref)),
        }
        print(f"  {ref}  vs  {name}")
        print(f"    same length ..... {same_len}  ({len(wref)} vs {len(w)})")
        print(f"    mispaired ....... {mispaired}/{k} = {100*frac:.2f}% of coordinates")
        print(f"    identical axis .. {out['pairings'][name]['identical_axis']}\n")

    # ---- 4. the guard MUST fire on the mismatched pair -------------------------------
    victim = next((n for n in axes if n != ref), None)
    fired = False
    if victim is not None:
        try:
            guard(axes[victim][0], wref)
        except SystemExit as exc:
            fired = True
            print(f"  GUARD FIRED as required: {exc}")
    out["guard"]["fires_on_mismatch"] = fired
    if not fired:
        raise SystemExit("FALSIFIER 5 FAILED: the guard did not fire on a mismatched axis")

    # ---- 5. and MUST pass on the matched pair ---------------------------------------
    guard(wref, wref)
    out["guard"]["passes_on_match"] = True
    print("  GUARD PASSES on the matched pair.\n")

    # ---- what the fold harness actually hands us ------------------------------------
    sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
    import plolo  # noqa: E402
    fold = plolo.build_fold("kolf2")
    fw = np.asarray(fold["wcols"])
    out["fold_kolf2"] = {"wcols": int(len(fw)), "matches_abc_axis": bool(np.array_equal(fw, wref))}
    print(f"  kolf2 fold wcols {len(fw)} | matches the A/B/C operator axis: "
          f"{np.array_equal(fw, wref)}")
    guard(fw, wref)

    (HERE / "work/panel226.json").write_text(json.dumps(out, indent=1))
    print("\nwrote work/panel226.json — FALSIFIER 5 PASS")


if __name__ == "__main__":
    main()
