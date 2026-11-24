"""196 — the four DE members for the ROUNDING RULE, on h1-test, three pools.

Same three pools as `six196.py` (native / wald89 / call-set-matched), same known-answer
rows, but the axis is `floor` against largest-remainder at production's depth rather than
the depth.  The `floor` arms are bit-identical to the depth sweep's `m = 1` arms
(`ham196.py` asserts it), so 185/191's a15 and rc_hybrid rows are the gate.

    uv run python sixham196.py
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/153-de-bench"))
from debench import Surface                                            # noqa: E402
from six196 import INV, WaldSurface, row, subset_rows, truncated_call   # noqa: E402

TAG = os.environ.get("TAG", "ham196")
ARMS = ["quasi_floor", "quasi_ham", "rc_hybrid_floor", "rc_hybrid_ham"]
KA = {"quasi_floor": dict(reach=0.18283, n_pred=6034, jac=0.06824, nmae=0.96659),
      "rc_hybrid_floor": dict(reach=0.21681, n_pred=6043, jac=0.06827, nmae=0.96670)}


def main() -> None:
    path = {a: WORK / f"{TAG}.{a}.npz" for a in ARMS if (WORK / f"{TAG}.{a}.npz").exists()}
    arms = [a for a in ARMS if a in path]
    meta = np.load(ROOT / "experiments/185-rank-controlled-emission/work/ctrl_meta.npz",
                   allow_pickle=True)
    cpm = meta["cpm"][meta["keep"]]

    S = {a: Surface(path[a]) for a in arms}
    print(f"h1-test rounding rule, {len(S[arms[0]].perts)} perturbations, "
          f"{S[arms[0]].genes.size} gated genes; median n_real "
          f"{np.median(S[arms[0]].real_sig.sum(1)):.0f}\n")

    print("KNOWN-ANSWER ROWS (rule 10) -- 185/191's a15 and rc_hybrid rows")
    bad = []
    for a, want in KA.items():
        if a not in S:
            continue
        got = row(S[a])
        for k, w in want.items():
            g = got[k]
            hit = abs(g - w) <= (0.5 if k == "n_pred" else 2e-5)
            print(f"  {a:<18} {k:<7} {g:>10.5f}  expect {w:<10.5f} {'PASS' if hit else 'FAIL'}")
            if not hit:
                bad.append(f"{a}.{k}")
    if bad and not os.environ.get("NOABORT"):
        raise SystemExit(f"ABORT (rule 10): {bad}")
    print("  all known-answer rows PASS\n")

    ok = np.ones(len(S[arms[0]].perts), bool)
    for a in arms:
        ok &= np.isfinite(S[a].pred_p).any(1)
    if not ok.all():
        print(f"PAIRED COHORT: dropping {int((~ok).sum())} perturbations "
              f"{[str(x) for x in S[arms[0]].perts[~ok]]}\n")
        for a in arms:
            subset_rows(S[a], ok)

    res = {}
    for pool, mk in (("native", lambda a: S[a]),
                     ("wald89", lambda a: subset_rows(WaldSurface(path[a], cpm), ok))):
        print(f"=== POOL: {pool} ===")
        print(f"{'arm':<18}{'n_pred':>8}{'reach':>9}{'fid':>9}{'jac':>9}{'nmae':>9}"
              f"{'prec':>9}{'d overall':>11}")
        for a in arms:
            r = row(mk(a))
            e = a.rsplit("_", 1)[0]
            if a.endswith("_floor"):
                res[(pool, e, "base")] = r
            b = res[(pool, e, "base")]
            d = ((r["reach"] - b["reach"]) * INV["reach"] + (r["fid"] - b["fid"]) * INV["fid"]
                 + (r["jac"] - b["jac"]) * INV["jac"] - (r["nmae"] - b["nmae"]) * INV["nmae"]) / 6
            res[(pool, a, "v")] = dict(r, d_overall=d)
            print(f"{a:<18}{r['n_pred']:>8.0f}{r['reach']:>9.5f}{r['fid']:>9.5f}"
                  f"{r['jac']:>9.5f}{r['nmae']:>9.5f}{r['prec']:>9.5f}{d:>+11.5f}")
        print()

    print("=== POOL: native, CALL SET MATCHED to the `floor` arm's own per-row n_pred ===")
    print(f"{'arm':<18}{'n_pred':>8}{'reach':>9}{'fid':>9}{'jac':>9}{'nmae':>9}"
          f"{'prec':>9}{'d overall':>11}")
    for e in ("quasi", "rc_hybrid"):
        b0 = f"{e}_floor"
        if b0 not in S:
            continue
        sizes = (S[b0].pred_sig & S[b0].adj).sum(1)
        base = None
        for a in (b0, f"{e}_ham"):
            r = row(S[a], call=truncated_call(S[a], sizes))
            base = base or r
            d = ((r["reach"] - base["reach"]) * INV["reach"]
                 + (r["fid"] - base["fid"]) * INV["fid"]
                 + (r["jac"] - base["jac"]) * INV["jac"]
                 - (r["nmae"] - base["nmae"]) * INV["nmae"]) / 6
            res[("matched", a, "v")] = dict(r, d_overall=d)
            print(f"{a:<18}{r['n_pred']:>8.0f}{r['reach']:>9.5f}{r['fid']:>9.5f}"
                  f"{r['jac']:>9.5f}{r['nmae']:>9.5f}{r['prec']:>9.5f}{d:>+11.5f}")
    json.dump({f"{k[0]}|{k[1]}": v for k, v in res.items()},
              open(WORK / "sixham196.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
