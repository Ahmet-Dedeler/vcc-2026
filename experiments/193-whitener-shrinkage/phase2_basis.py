"""193 Phase 2 — the BASIS.  How many donor effect vectors should `C` be fitted on?

Production's shipped operator (`data/donors/whitener_union_pruned.npz`, sha `e2f65ecf`) was
fitted on **29,394** rows: `k562_gwps` whole (9,787), the two Orion lines subsampled to a
12,000-row cap, then the common-essential perturbations dropped (experiment 091).  Two other
operators on the SAME 17,942-gene union axis already exist on disk:

    whitener_union.npz       sha 5608e807   33,787 rows  (12,000 cap, no prune)
    whitener_union_full.npz  sha 965d3cef   45,531 rows  (no cap, no prune -- every donor row)

`union -> union_full` is a PURE SIZE increase (9,787 + 17,767 + 17,977, nothing dropped);
`union -> union_pruned` is 091's composition change.  So the two questions the mandate asks --
"is `C` fitted on too little data" and "does pruning still pay" -- are both expressible without
fitting anything new, on the ONE fold where an operator swap is leakage-free.

LEAKAGE (162's rule, restated because it is the whole reason only one fold is used)
-----------------------------------------------------------------------------------
All three of these operators contain rows from all three production donors.  Scoring
`orion_hct116`, `orion_hek293t` or `k562_gwps` through them would put the held library's own
perturbations in the basis.  `kolf2` is not a production donor and contributed no row to any of
them, so it is the only fold on which these three operators can be compared at all.  This is
exactly the constraint 162 imposed on its own `a2-prune` arm.

KNOWN-ANSWER ROW (rule 10)
--------------------------
162 measured the prune at the v11 configuration (alpha 1.0796, fixed CD4 fraction 0.05):
`union -> union_pruned` = **+0.001360** raw `pds` on kolf2, challenge frame.  This script
recomputes it from the operator files, so a mis-wired operator swap cannot pass silently.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))

import plolo                                                        # noqa: E402
from plolo import build_fold, deliver, endpoints, solve_fraction    # noqa: E402
from run_arms import CONFIGS                                        # noqa: E402
from screen193 import SPAN_PDS, cd4_cached, frames_for, sha256_file  # noqa: E402

OPERATORS = {
    "union_pruned(PROD)": (ROOT / "data/donors/whitener_union_pruned.npz", 29394),
    "union":              (ROOT / "data/donors/whitener_union.npz", 33787),
    "union_full":         (ROOT / "data/donors/whitener_union_full.npz", 45531),
}


def score(fold, cfg, op_path, cd4, frames, truth, wcols):
    with np.load(op_path) as z:
        evals = np.asarray(z["evals"], float)
        evecs = np.asarray(z["evecs"], float)
    assert evecs.shape[0] == len(wcols) == len(evals)
    pre = np.asarray(fold["pre"], dtype=np.float64) * cfg["alpha"]
    sub = pre[:, wcols]
    before = float(np.median(np.linalg.norm(sub, axis=1))) or 1.0
    w = (evals + cfg["lam"]) ** (-cfg["gamma"] / 2.0)
    t = ((sub @ evecs) * w) @ evecs.T
    after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
    whitened = pre.copy()
    whitened[:, wcols] = t * (before / after)
    del t, evecs, sub
    if "cd4_f" in cfg:
        frac, _ = solve_fraction(whitened, cd4, cfg["cd4_f"])
    else:
        frac = cfg["cd4_frac"]
    final, _ = plolo.E062.add_post_whiten_residual(whitened, cd4, frac, wcols)
    final[:, fold["panel_positions"]] = 0.0
    del whitened
    per = {}
    for fname, b, gate in frames:
        mat = deliver(final, b, gate, cfg["mass"])
        e = endpoints(mat[:, fold["scored"]], truth)
        per[fname] = dict(pds=e["pds"], cos=e["cos"], sgn_w=e["sgn_w"], sse=e["sse"])
        del mat
    ei = endpoints(final[:, fold["scored"]], truth)
    return dict(chal_pds=float(np.mean([v["pds"] for v in per.values()])),
                intended_pds=ei["pds"], per_frame=per, max_ev=float(evals.max()),
                cd4_fraction=float(frac))


def main() -> None:
    held = "kolf2"
    fold = build_fold(held)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    wcols = np.asarray(fold["wcols"])
    cd4 = cd4_cached(fold["targets"])
    frames = frames_for(fold)
    out = {"fold": held, "n_targets": len(fold["targets"]),
           "scored_genes": int(len(fold["scored"])), "arms": {}}

    print("KNOWN-ANSWER ROW: v11 config, union -> union_pruned, kolf2 challenge frame")
    v11 = CONFIGS["v11"]
    ka = {}
    for name in ("union", "union_pruned(PROD)"):
        p, n = OPERATORS[name]
        r = score(fold, v11, p, cd4, frames, truth, wcols)
        ka[name] = r
        print(f"   {name:20s} sha {sha256_file(p)[:8]}  chal_pds {r['chal_pds']:.6f}")
    d = ka["union_pruned(PROD)"]["chal_pds"] - ka["union"]["chal_pds"]
    print(f"   prune effect = {d:+.6f}   162 expected +0.001360   |diff| {abs(d-0.001360):.2e}")
    ok = abs(d - 0.001360) < 5e-6
    print("   GATE 3 " + ("PASSED" if ok else "FAILED"))
    out["gate3"] = dict(delta=d, expected=0.001360, pass_=bool(ok))
    if not ok:
        (HERE / "work/phase2-basis.json").write_text(json.dumps(out, indent=1))
        raise SystemExit("GATE 3 FAILED — operator swap path is not reproducing 162")

    print("\na15 config, basis sweep at production's f(ev):")
    a15 = CONFIGS["a15"]
    for name, (p, n) in OPERATORS.items():
        r = score(fold, a15, p, cd4, frames, truth, wcols)
        r.update(path=str(p.relative_to(ROOT)), sha256=sha256_file(p), n_basis=n)
        out["arms"][name] = r
        print(f"   {name:20s} n={n:6d} sha {r['sha256'][:8]} maxev {r['max_ev']:.4f} "
              f"chal_pds {r['chal_pds']:.6f}  intended {r['intended_pds']:.6f}")
    base = out["arms"]["union_pruned(PROD)"]["chal_pds"]
    print()
    for name, r in out["arms"].items():
        dd = r["chal_pds"] - base
        r["delta_vs_prod"] = dd
        r["delta_scaled"] = dd / SPAN_PDS
        print(f"   {name:20s} d_pds {dd:+.6f}  scaled {dd/SPAN_PDS:+.6f}")
    (HERE / "work/phase2-basis.json").write_text(json.dumps(out, indent=1))
    print("\nwrote work/phase2-basis.json")


if __name__ == "__main__":
    main()
