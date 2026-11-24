"""213 — 204's paired target-cluster bootstrap on the law sweep.

The unit is a (fold, target) pair; the resample is over TARGETS, paired across laws, exactly as
204 and 212 did it, so what comes out is an interval on the DIFFERENCE and not on the scatter
of four folds.  204 calibrated this instrument at +-4.3e-5 absolute / 0.20% relative per arm on
these same folds; 212 reproduced 3.90e-5 / 0.178%.  Anything inside that is not a finding.

    uv run python boot213.py [--src work/plolo213_main.json]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]
REF = "logc20"
DRAWS = 4000
SEED = 20260904


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", default=str(WORK / "plolo213_main.json"))
    ap.add_argument("--out", default=None)
    ap.add_argument("--render", default="a15", choices=["a15", "a23"])
    args = ap.parse_args()
    rows = json.loads(Path(args.src).read_text())
    by = {(r["fold"], r["law"]): r for r in rows}
    tag = f"_cos_per_target_{args.render}"
    ptag = f"_pds_per_target_{args.render}"
    # the rank laws have no a23 render (no analytic density for the minnorm mass fix)
    laws = [l for l in dict.fromkeys(r["law"] for r in rows)
            if all((f, l) in by and tag in by[(f, l)] for f in FOLDS)]

    rng = np.random.default_rng(SEED)
    out = {}
    for f in FOLDS:
        n = len(by[(f, REF)][tag])
        out.setdefault("_n", {})[f] = n
    idx = {f: rng.integers(0, out["_n"][f], size=(DRAWS, out["_n"][f])) for f in FOLDS}

    def macro(law, key):
        return float(np.mean([np.mean(by[(f, law)][key]) for f in FOLDS]))

    def boot(law, key):
        """DRAWS macro means of `law`, and of the reference, on the SAME resample."""
        a = np.zeros(DRAWS)
        b = np.zeros(DRAWS)
        for f in FOLDS:
            va = np.asarray(by[(f, law)][key])
            vb = np.asarray(by[(f, REF)][key])
            a += va[idx[f]].mean(axis=1) / len(FOLDS)
            b += vb[idx[f]].mean(axis=1) / len(FOLDS)
        return a, b

    res = []
    for law in laws:
        row = dict(law=law, render=args.render,
                   cos=macro(law, tag), pds=macro(law, ptag),
                   d_cos=macro(law, tag) - macro(REF, tag),
                   d_pds=macro(law, ptag) - macro(REF, ptag),
                   ratio_cos=macro(law, tag) / macro(REF, tag))
        for key, name in ((tag, "cos"), (ptag, "pds")):
            a, b = boot(law, key)
            d = a - b
            lo, hi = np.percentile(d, [2.5, 97.5])
            row[f"ci_{name}"] = [float(lo), float(hi)]
            row[f"halfwidth_{name}"] = float((hi - lo) / 2.0)
            row[f"p_gt0_{name}"] = float((d > 0).mean())
        res.append(row)
        print(f"  {law:<10} cos {row['cos']:.7f}  ratio {row['ratio_cos']:.4f}  "
              f"dcos {row['d_cos']:+.7f} [{row['ci_cos'][0]:+.6f}, {row['ci_cos'][1]:+.6f}] "
              f"P>0 {row['p_gt0_cos']:.3f} | dpds {row['d_pds']:+.6f} "
              f"[{row['ci_pds'][0]:+.5f}, {row['ci_pds'][1]:+.5f}]", flush=True)
    path = Path(args.out) if args.out else WORK / f"boot213_{args.render}.json"
    path.write_text(json.dumps(res, indent=1))
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
