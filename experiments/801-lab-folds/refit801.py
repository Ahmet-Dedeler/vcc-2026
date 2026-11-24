"""801 -- how much does a fold's operator refit move a price? The same efp fold dumped twice: with 238's operator (fitted on
the Mac; the KA dump) and with an EC2 refit by the same command (a different sha256: fits are not bit-reproducible across
machines). Prices the named arms on both dumps and reports the absolute differences (cos_agg, literal pds, pds lam 0.5,
K_mean) and the difference of the incumbent-vs-arm r. Information only: nothing is asserted, nothing is tuned.

    uv run python experiments/801-lab-folds/refit801.py <held>      -> work/refit801-<held>.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import json  # noqa: E402

import engine801 as E  # noqa: E402

ARMS = ("cap3", "cap6", "pct", "incumbent", "pct+S_P2", "cap3+S_P2", "cap3+P2g")


def main(held):
    A = E.DumpFold(f"efp/{held}", HERE / "work/dumps" / f"dump708-efp-{held}.npz")
    B = E.DumpFold(f"efp/{held}", HERE / "work/dumps" / f"dump708-efprefit-{held}.npz")
    assert A.targets == B.targets
    out = dict(held=held, n=A.n, kv_max_rel=float(abs(A.kv - B.kv).max() / abs(A.kv).max()),
               fin_max_abs=float(abs(A.fin - B.fin).max()), fin_scale=float(abs(A.fin).max()), arms={})
    sA, sB = {}, {}
    for nm in ARMS:
        sA[nm] = E.score_rows_of(A, *E.render(A, E.named(A, nm)))[0]
        sB[nm] = E.score_rows_of(B, *E.render(B, E.named(B, nm)))[0]
        out["arms"][nm] = {k: abs(sA[nm][k] - sB[nm][k]) for k in ("cos_agg", "pds", "pds_lam0.5", "K_mean")}
    for nm in ARMS:
        rA = sA["incumbent"]["cos_agg"] / sA[nm]["cos_agg"] - 1
        rB = sB["incumbent"]["cos_agg"] / sB[nm]["cos_agg"] - 1
        out["arms"][nm]["r_incumbent_vs_arm_mac_op"] = rA
        out["arms"][nm]["r_incumbent_vs_arm_refit_op"] = rB
    (HERE / "work" / f"refit801-{held}.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main(sys.argv[1])
