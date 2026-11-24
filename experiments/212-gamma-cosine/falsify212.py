"""212 Phase 4 — falsifier 3: the PERMUTED-OPERATOR control, swept in gamma.

209's random control was the maximiser over all 36 of its arms. That is the outcome to watch
for, and the only way to watch for it is to put a deliberately meaningless operator through
the identical sweep and see where it lands.

THE CONTROL, AND WHY THIS ONE
-----------------------------
The operator is `W = (B'B/n)^(-gamma/2)` and it acts on the effect matrix's OPERATOR COLUMNS
(`wcols`, the union of donor coverage minus the panel). Permuting the *rows* of the donor
basis `B` leaves `B'B` unchanged, so a row permutation is not a control at all. The
structure the operator carries is entirely in **which gene each coordinate is paired with**,
so the control is a common random permutation of that pairing:

    W_perm = P W P'    <=>    apply the same operator through a permuted `wcols`.

This has the IDENTICAL eigenvalue spectrum, the identical condition number (7,943x), the
identical amplification at every gamma, and identical everything except that 100% of the
coordinates are attached to the wrong gene. It is precisely AGENTS.md 8.3's Oct-22 trap --
the failure that passes every check the builder runs -- reproduced deliberately.

WHAT EACH OUTCOME MEANS
-----------------------
* If the control's **cosine** also rises as gamma falls, the cosine rise is NOT evidence for
  the operator's structure: it is the trivial statement that less whitening leaves more bulk
  cosine. Expected, and it is why the cosine alone cannot select gamma.
* If the control's **pds** does NOT collapse, the pds channel is not reading the operator and
  the whole screen is void.
* If the control's composite BEATS the real operator at any gamma, this experiment is broken
  and reports nothing. That is 209's outcome and it is the one to watch for.

A second control is run alongside: `spectrum_shuffle`, the same eigenvectors with the
eigenvalues randomly permuted among them. That keeps the gene pairing and the multiset of
eigenvalues and destroys only the pairing of a direction with its variance.

    uv run python falsify212.py
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

import score204  # noqa: E402  # noqa: F401  (kept for the identical import order)
from plolo212 import DEV_KA, GRID, score_at, sha256_file  # noqa: E402
from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for  # noqa: E402
from score204 import DONORS, fold_sources, weighted_baseline  # noqa: E402

SEED = 212
CONTROLS = ("real", "gene_perm", "spectrum_shuffle")
# a coarser grid than the primary sweep -- the control only has to be located, not resolved.
CGRID = [0.00, 0.30, 0.40, 0.50, 0.575, 0.65, 0.75, 0.80]


class Permuted(Operator):
    """The shipped operator with its gene pairing or its spectrum deliberately scrambled."""

    def __init__(self, path, kind, rng):
        super().__init__(path)
        self.kind = kind
        n = self.evals.size
        if kind == "gene_perm":
            self.perm = rng.permutation(n)
        elif kind == "spectrum_shuffle":
            self.evals = self.evals[rng.permutation(n)]
            self.perm = None
        else:
            self.perm = None

    def apply(self, M, wcols, gamma, lam):
        cols = wcols if self.perm is None else np.asarray(wcols)[self.perm]
        # `before`/`after` medians are taken on the SAME columns the operator writes, so the
        # renormalisation is the production one; only the pairing moved.
        return super().apply(M, cols, gamma, lam)


def main() -> None:
    out_path = WORK / "falsify212.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["gamma"], r["control"]) for r in rows}
    t0 = time.time()
    for held in FOLDS:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        srcs = fold_sources(held, fold)
        cd4 = cd4_for(fold["targets"])
        op_path = Path(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                       else PROD_OPS / "whitener_union_pruned.npz")
        sha = sha256_file(op_path)
        pre, aud = weighted_baseline(srcs, fold["targets"], "equal",
                                     {d: 1.0 for d in DONORS})
        dev = float(np.max(np.abs(pre - np.asarray(fold["pre"], dtype=np.float64))))
        assert dev == DEV_KA[held], f"`equal` drifted on {held}: {dev!r}"
        for kind in CONTROLS:
            op = Permuted(op_path, kind, np.random.default_rng(SEED))
            for gamma in CGRID:
                if (held, gamma, kind) in done:
                    continue
                e = score_at(fold, pre, op, cd4, fold["scored"], truth, gamma)
                rows.append(dict(fold=held, gamma=gamma, control=kind,
                                 operator_sha256=sha, seed=SEED,
                                 n_targets=len(fold["targets"]), **e))
                out_path.write_text(json.dumps(rows, indent=1))
                print(f"[{held}/{kind}/g{gamma:.3f}] cos {e['a23_cos']:.7f} "
                      f"a15cos {e['a15_cos']:.7f} pds {e['a23_pds']:.6f} "
                      f"rr {e['rr']:.6f} ({time.time()-t0:.0f}s)", flush=True)
            del op
            gc.collect()
        del srcs, cd4, pre
        gc.collect()

    # --------------------------------------------------------------- identity assertions
    folds = sorted(FOLDS)
    prim = {(r["fold"], r["gamma"]): r
            for r in json.loads((WORK / "plolo212.json").read_text())}
    print("\n===== FALSIFIER 4: `real` must reproduce the primary sweep EXACTLY =====")
    ok = True
    for r in rows:
        if r["control"] != "real":
            continue
        p = prim.get((r["fold"], r["gamma"]))
        if p is None:
            continue
        for k in ("a23_cos", "a15_cos", "a23_pds", "a15_pds", "a15_sgn_w"):
            same = r[k] == p[k]
            ok &= same
            if not same:
                print(f"  MISMATCH {r['fold']} g{r['gamma']} {k}: {r[k]!r} vs {p[k]!r}")
    print(f"  ==> {'IDENTICAL on every checked scalar' if ok else 'FAILED'}; the operator "
          f"sha is asserted at load, which is the gate 170 warns the numeric identity "
          f"check cannot be (orthonormality reproduces the baseline at production "
          f"settings whatever basis is used).")

    print("\n===== the control, macro over four folds =====")
    print("gamma |" + "".join(f"  {k:>21}" for k in CONTROLS))
    print("      |" + "".join("      cos       pds  " for _ in CONTROLS))
    for g in CGRID:
        line = f"{g:5.3f} |"
        for kind in CONTROLS:
            sel = [r for r in rows if r["control"] == kind and r["gamma"] == g
                   and r["fold"] in folds]
            if len(sel) != len(folds):
                line += "          --        --  "
                continue
            line += (f"  {np.mean([s['a23_cos'] for s in sel]):.6f}  "
                     f"{np.mean([s['a23_pds'] for s in sel]):.6f}  ")
        print(line)


if __name__ == "__main__":
    main()
