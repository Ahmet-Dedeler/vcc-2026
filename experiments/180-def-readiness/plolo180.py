"""180 — `pds` for the knee alphas, on 162's PLOLO, in the challenge frames.

Rule 13 names PLOLO as THE `pds` instrument (board-validated: LOO sigma 0.00149 raw, sign
6/6 over 12 board points), and 178 measured h1-test understating `pds` costs ~6x on
un-whitening arms. So the knee's `pds` is priced here and not on h1-test.

It is scored in **delivered** space (`profile_from_effect` + the scorer's read-back), because
156 measured `identity/intended` beating `identity/delivered` by +0.0032/+0.0038 raw pds --
an effect-level `pds` measures an object the board never sees.

WHY THIS IS NOT VACUOUS. `pds` is `pds_cosine`, and a cosine is invariant to a positive
per-row scale, so "a per-context alpha cannot move `pds`" is close to an identity -- close,
because a23's minimum-norm mass fix makes the renderer non-homogeneous, and 180's own
`knee180.py` measures the delivered row direction moving by up to 1 - cos = 1.3e-4. A screen
that can only ever return 0.000000 proves nothing, so two arms with a DELIBERATELY LARGE
amplitude change (alpha 0.90 and 1.10 at a23's exact configuration) are run alongside: if
those read 0.000000 too, this instrument is blind to amplitude and its verdict on the knee is
worthless. 162's own recorded a11 (alpha 0.838) vs a14 (0.982) rows say it is not blind --
they differ by -0.0004 to -0.0016 raw pds across four folds.

KNOWN-ANSWER ROW: `a23` re-run here must reproduce 162's own recorded `arms.json` value for
the same (fold, space) bit-for-bit.

    uv run python plolo180.py [fold ...]
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
EXP162 = REPO / "experiments/162-retrodicting-validator"
WORK = HERE / "work"
sys.path.insert(0, str(EXP162))

A23 = dict(alpha=0.982, gamma=0.65, lam=0.0, cd4_f=0.10, mass="minnorm", op="union_pruned")


def main(folds) -> int:
    os.chdir(EXP162)
    import plolo
    import run_arms
    from plolo import FOLDS, PROD_OPS, Operator, build_fold, cd4_for, deliver, endpoints

    knee = json.loads((WORK / "knee180.json").read_text())
    configs = {"a23": dict(A23)}
    for context in ("A", "B", "C"):
        configs[f"knee_{context}"] = dict(
            A23, alpha=float(knee["per_context"][context]["alpha"])
        )
    # the non-vacuity pair: same configuration, deliberately large amplitude moves
    configs["amp_0.90"] = dict(A23, alpha=0.90)
    configs["amp_1.10"] = dict(A23, alpha=1.10)
    print("[plolo] arms: " + "  ".join(f"{k}={v['alpha']:.6f}" for k, v in configs.items()),
          flush=True)

    incumbent = {
        (r["fold"], r["arm"], r["space"]): r["pds"]
        for r in json.loads((EXP162 / "work/arms.json").read_text())
    }

    rows = []
    for held in folds:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        frames = run_arms.frames_for(fold)
        path = (PROD_OPS / "whitener_union_pruned.npz"
                if FOLDS[held]["op_kind"] == "production" else fold["op_path"])
        operator = Operator(path)
        print(f"[{held}] operator {path.name} maxeig {operator.evals.max():.4f}", flush=True)
        for arm, cfg in configs.items():
            final, caud = run_arms.effects_for(fold, cfg, operator, cd4)
            for fname, b, gate in frames:
                if not fname.startswith("chal"):
                    continue
                matrix = deliver(final, b, gate, cfg["mass"])
                e = endpoints(matrix[:, fold["scored"]], truth)
                space = f"{fname}/delivered"
                rows.append(dict(fold=held, arm=arm, alpha=cfg["alpha"], space=space,
                                 pds=float(e["pds"]), cos=float(e["cos"]),
                                 sse=float(e["sse"])))
                mark = ""
                if arm == "a23" and (held, arm, space) in incumbent:
                    gap = abs(e["pds"] - incumbent[(held, arm, space)])
                    mark = f"  | 162 recorded {incumbent[(held, arm, space)]:.6f} d {gap:.1e}"
                    if gap > 1e-12:
                        raise SystemExit(
                            f"KNOWN-ANSWER FAIL: a23 {held}/{space} {e['pds']} != "
                            f"{incumbent[(held, arm, space)]}"
                        )
                print(f"  {held:<14}{arm:<10}{space:<18} pds {e['pds']:.6f} "
                      f"cos {e['cos']:.5f}{mark}", flush=True)
            (WORK / "plolo180.json").write_text(json.dumps(rows, indent=1))
        del operator

    print("\n[plolo] delta pds against a23, per fold and challenge frame "
          "(162's LOO sigma on this instrument is 0.00149 raw):")
    for held in folds:
        base = {r["space"]: r["pds"] for r in rows if r["fold"] == held and r["arm"] == "a23"}
        for arm in configs:
            if arm == "a23":
                continue
            deltas = {r["space"]: r["pds"] - base[r["space"]]
                      for r in rows if r["fold"] == held and r["arm"] == arm}
            print(f"  {held:<14}{arm:<10}" + "  ".join(
                f"{k.replace('/delivered', '')} {v:+.6f}" for k, v in sorted(deltas.items())))
    # the arm that ships: each context reads its OWN knee alpha in its OWN frame
    print("\n[plolo] the arm as it would ship -- context c at alpha_c, in frame chal c:")
    for held in folds:
        base = {r["space"]: r["pds"] for r in rows if r["fold"] == held and r["arm"] == "a23"}
        total = 0.0
        parts = []
        for context in ("A", "B", "C"):
            space = f"chal{context}/delivered"
            got = [r["pds"] for r in rows
                   if r["fold"] == held and r["arm"] == f"knee_{context}"
                   and r["space"] == space][0]
            delta = got - base[space]
            total += delta
            parts.append(f"{context} {delta:+.6f}")
        print(f"  {held:<14}" + "  ".join(parts) + f"  | mean {total / 3:+.6f} raw pds "
              f"= {total / 3 / 0.455 / 6:+.6f} overall")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or ["kolf2"]))
