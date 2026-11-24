"""534 Part B -- the per-step values FROM THE RECORD (A/B/C, real, efp, efr), quoted, not recomputed.

Two kinds of record value, each labelled with its source:
  * RESULTS text: 438 sections 3.1 / 3.4 / 3.5 and section 1's CTRL row, 522 section 7, FINDINGS section 0's board chain. Where
    the record never priced the exact milestone design (h2c1 at q 1.75, the LODO key at q 1.75, the percentile cap in
    438's levelled route), the closest recorded row is quoted and named (e.g. the q 1.5 -> 1.75 step on the CD4 x1.5
    stack, the LODO key step at q 1.5).
  * 438's STORED rows through 438's OWN pricer (price438.Panel, read-only, 438's jsons): the exact milestone steps whose
    two designs are both in 438's grid (M0->M1, M1->M2, M2->M3). This reads the record; it re-derives nothing (the
    identical pricer on the identical stored statistics; M1->M2 equals 438 section 3.1's printed H1 x2 half).

    uv run python record534.py            -> work/partB/record534.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import designs534 as D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))

# RESULTS-text values (levelled k = 1 unless noted), quoted verbatim from the files named
TEXT = {
    "M0->M1": dict(abc=None, real=None, efp=None, efr=None, board=0.0069612,
                   src="no October-like re-price of the 417a step in the record (438 prices everything against the 417a "
                       "analogue); board: 417a 0.2435334 - 410b 0.2365722 (FINDINGS section 0)"),
    "M1->M2": dict(abc=0.00238, real=0.00511, efp=0.00255, efr=0.00031, board=None,
                   src="438 section 3.1 'H1 x2 half' (exact design); board: 417e +0.0051 is the whole stack incl. CD4 x1.5"),
    "M2->M3": dict(abc=0.00847, real=0.00492, efp=0.00412, efr=0.00454, board=0.0046474,
                   src="438 section 3.4 'q 1.5, levelled k = 1' on the CD4 x1.5 stack (analogue; on PP: +0.00844 / +0.00481 / "
                       "+0.00373 / +0.00419); board: 433k 0.2677594 - 433h 0.2631120"),
    "M3->M4": dict(abc=0.01311 - 0.00847, real=0.00746 - 0.00492, efp=0.00607 - 0.00412, efr=0.00758 - 0.00454, board=0.0022834,
                   src="438 section 3.4, the stack's 'q 1.75' row minus its 'q 1.5' row (analogue: CD4 x1.5, K3; Part A's "
                       "'a further +0.0020 ... +0.0030'); board: 446b 0.2700428 - 433k 0.2677594"),
    "M4->M5": dict(abc=0.00329, real=0.00343, efp=0.00349, efr=0.00187, board=0.0022338,
                   src="438 section 3.5 'LODO V1 clip (0.5, 3), q 1.5 (vs K3 q 1.5)' on the stack (analogue at q 1.5; at q 1.25: "
                       "+0.00274 / +0.00318 / +0.00306 / +0.00169); board: 446bj2 0.2722766 - 446b 0.2700428"),
    "M5->M6": dict(abc=0.0021, real=0.0095, efp=0.0087, efr=0.0029, board=None,
                   src="522 section 7 (price522: vs cap 3, k 1, 522's unlevelled route, retention-scaled; at k 1.31: +0.0026 / "
                       "+0.0123 / +0.0114 / +0.0037); board: 495 (cap 6) reads 10-03, pending"),
    "C-a": dict(abc=0.00170, real=-0.00160, efp=-0.00111, efr=-0.00236, board=0.0031484,
                src="438 section 3.1 'CD4 x1.5 on top of H1 x2' at q 1.25 (analogue at q 1.25); board: 457dr's CD4 half +0.0031"),
    "C-b": dict(abc=-0.02225, real=-0.02351, efp=-0.01489, efr=-0.00476, board=None,
                src="438 section 1 / 3.2 CTRL: the LODO key (0.3, 3) permuted, on PP at q 1.25, UNLEVELLED (analogue)"),
}


def stored():
    """Exact milestone steps from 438's stored rows via 438's own pricer (levelled k = 1, cancer-fold mean)."""
    import price438 as P8
    out = {}
    panels = {"abc": P8.Panel("abc"), "real": P8.Panel("real"), "efp": P8.Panel("efp"),
              "efr": P8.Panel("efr", full_tag="efr_s1", energy="fold")}
    for sid in ("M0->M1", "M1->M2", "M2->M3"):
        st = [s for s in D.STEPS if s[0] == sid][0]
        a, r = D.M[st[1]], D.M[st[2]]
        out[sid] = {}
        for t, P in panels.items():
            v = {f: P.price(f, a, 1.0, r, True)["ov"] for f in ("k562_gwps", "orion_hct116", "orion_hek293t")}
            out[sid][t] = dict(mean=float(np.mean(list(v.values()))), n_pos=int(sum(x > 0 for x in v.values())), per_fold=v)
    return out


def main():
    rec = dict(text=TEXT, stored=stored())
    for sid, v in rec["stored"].items():
        print(sid, " ".join(f"{t} {x['mean']:+.5f} ({x['n_pos']}/3)" for t, x in v.items()))
    s = rec["stored"]["M1->M2"]
    for t in ("abc", "real", "efp", "efr"):                    # the stored route reproduces the printed H1 x2 half
        assert abs(round(s[t]["mean"], 5) - TEXT["M1->M2"][t]) < 6e-6, (t, s[t]["mean"], TEXT["M1->M2"][t])
    (HERE / "record534.json").write_text(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
