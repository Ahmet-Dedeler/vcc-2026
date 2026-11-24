"""534 Part B -- descriptive extras for RESULTS-B.md (light; reads the priced jsons, changes no reading rule):
  * each step's decomposition per panel (levelled k = 1): d_rel (relative aggregate cosine), board-emulated dpds, dmse raw;
  * the exact milestone designs on A/B/C's own panel (supplementary) and on efp (the known-answer panel), at k = 1 and at
    the calibrated transfer, so the board-unit projection has an A/B/C counterpart in the same instrument;
  * the per-panel calibrated cumulative.

    uv run python extras534.py        -> work/partB/extras534.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402

import numpy as np  # noqa: E402

import designs534 as D  # noqa: E402
import price534 as P5  # noqa: E402

MAIN = [s[0] for s in D.MAIN_STEPS]


def decomp(P):
    out = {}
    for st in D.STEPS:
        a, r = D.M[st[1]], D.M[st[2]]
        pr = {f: P.price(f, a, 1.0, r, True) for f in P.folds}
        out[st[0]] = dict(ov=float(np.mean([x["ov"] for x in pr.values()])),
                          d_rel=float(np.mean([x["d_rel"] for x in pr.values()])),
                          dpds_bd=float(np.mean([x["dpds_bd"] for x in pr.values()])),
                          dpds_lit=float(np.mean([x["dpds_lit"] for x in pr.values()])),
                          dmse=float(np.mean([x["dmse"] for x in pr.values()])),
                          per_fold={f: x["ov"] for f, x in pr.items()})
    return out


def main():
    out = dict(decomposition={}, exact={})
    for p in ("abc", "efp", "S1", "S2", "S3", "S4"):
        P = P5.Panel(p)
        out["decomposition"][p] = decomp(P)
        k1 = P5.point_table({p: P}, level=True)[p]
        cal = P5.point_table({p: P}, level=True, cal=True)[p]
        out["exact"][p] = dict(k1={s: k1[s]["mean"] for s in k1}, calibrated={s: cal[s]["mean"] for s in cal},
                               cum_k1=float(sum(k1[s]["mean"] for s in MAIN)), cum_calibrated=float(sum(cal[s]["mean"] for s in MAIN)),
                               n_pos={s: k1[s]["n_pos"] for s in k1}, ret=P.ret)
    (P5.HERE / "extras534.json").write_text(json.dumps(out, indent=1))
    for p, e in out["exact"].items():
        print(f"{p:4s} cum k1 {e['cum_k1']:+.5f}  cum calibrated {e['cum_calibrated']:+.5f}  "
              + " ".join(f"{s} {e['calibrated'][s]:+.5f}" for s in MAIN))


if __name__ == "__main__":
    main()
