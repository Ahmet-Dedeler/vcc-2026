"""195 Phase 0 — known-answer rows.  Abort the experiment if these do not reproduce.

Target (AGENTS.md 4.2 / experiments/019-whiten-interact/RESULTS.md):

    4:1, gamma = 0     -> 0.7378   (incumbent)
    4:1, gamma = 0.55  -> 0.7873   (the "+ whitening" row)
    1:1, gamma = 0.65  -> 0.7978   (equal donor weights, what shipped)
    1:1, gamma = 0.55  -> 0.7963

Rule 11: the fitted operator for every fold is written to disk with its sha256, so that every
later arm can name a path and a hash rather than a description.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from harness import (LODO, apply_operator, build_fold, common_perts, fit_operator,  # noqa: E402
                     fold_basis_rows, load_banks, pds_cosine, sha256_arr)

HERE = Path(__file__).resolve().parent
OPS = HERE / "work" / "operators"

KNOWN = {("4:1", 0.0): 0.7378, ("4:1", 0.55): 0.7873,
         ("1:1", 0.65): 0.7978, ("1:1", 0.55): 0.7963}
WEIGHTS = ["1:0", "8:1", "4:1", "2:1", "1:1", "1:2"]
GAMMAS = [0.0, 0.4, 0.5, 0.55, 0.65, 0.75, 0.9, 1.1]


def main() -> None:
    OPS.mkdir(parents=True, exist_ok=True)
    banks = load_banks()
    common = common_perts()
    brows = fold_basis_rows(banks)
    print(f"panel perturbations {len(common)}   folds {LODO}", flush=True)

    folds, opmeta = {}, {}
    for held in LODO:
        t0 = time.time()
        f = build_fold(held, banks, common, brows[held])
        ev, V = fit_operator(f["basis"])
        del f["basis"]
        p = OPS / f"op.{held}.npz"
        np.savez(p, evals=ev, evecs=V)
        sha = sha256_arr(ev, V)
        opmeta[held] = {"path": str(p.relative_to(HERE.parents[1])), "sha256": sha,
                        "max_eigenvalue": float(ev.max()), "min_eigenvalue": float(ev.min()),
                        "n_basis": int(BASIS_N := sum(len(v) for v in [brows[held][n] for n in f["others"]])),
                        "ref": f["ref"], "sc": f["sc"]}
        f["ev"], f["V"] = ev, V
        folds[held] = f
        print(f"fold held={held:14} basis={BASIS_N}x{V.shape[0]}  max_ev={ev.max():.4f}  "
              f"sha={sha[:16]}  {time.time()-t0:.1f}s", flush=True)

    rows = {}
    for wspec in WEIGHTS:
        ws = [float(x) for x in wspec.split(":")]
        for g in GAMMAS:
            per = {}
            for h, f in folds.items():
                w = dict(zip(f["others"], ws))
                tot = sum(w.values())
                P = sum((w[n] / tot) * f["parts"][n] for n in f["others"])
                per[h] = pds_cosine(apply_operator(P, f["ev"], f["V"], g), f["T"])
            macro = float(np.mean(list(per.values())))
            rows[f"{wspec}|{g}"] = {**per, "macro": macro}
            tag = ""
            if (wspec, g) in KNOWN:
                d = macro - KNOWN[(wspec, g)]
                tag = f"   KNOWN-ANSWER {KNOWN[(wspec, g)]:.4f}  delta {d:+.5f}"
            fs = "  ".join(f"{k.replace('orion_', '')[:9]}={v:.4f}" for k, v in per.items())
            print(f"  w={wspec:5} gamma={g:<5g} macro {macro:.4f}  [{fs}]{tag}", flush=True)

    ok = True
    for (wspec, g), want in KNOWN.items():
        got = rows[f"{wspec}|{g}"]["macro"]
        if abs(got - want) > 5e-4:
            ok = False
            print(f"MISMATCH w={wspec} gamma={g}: got {got:.4f} want {want:.4f}")
    (HERE / "work" / "phase0.json").write_text(json.dumps(
        {"rows": rows, "operators": opmeta, "known_answer_pass": ok}, indent=2))
    print("\nPhase 0", "PASS" if ok else "FAIL")
    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
