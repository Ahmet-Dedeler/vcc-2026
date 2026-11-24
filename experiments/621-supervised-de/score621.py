#!/usr/bin/env python3
"""621 — score the supervised DE heads (and same-axis baselines) on d4cd5f's 636 instrument, paired vs CLEAN4 LOLO.

Rows per surface: the model's direction head (P_down), the family-balanced donor consensus sign, the gene-level
programme (family-balanced), and a confidence switch (consensus where |Zbal| >= 3, programme elsewhere); the
membership head's AUROC / precision@n_real. Writes work/score621.json and prints a table.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
from inst636 import Inst, clean4_lolo  # noqa: E402


def main() -> None:
    I = Inst()
    C4 = {s: clean4_lolo(I, s)[0] for s in I.S}
    P, M, CONS, PROG = {}, {}, {}, {}
    for s in I.S:
        z = np.load(HERE / f"work/pred-{s}.npz", allow_pickle=True)
        assert list(z["targets"].astype(str)) == list(I.S[s]["T"]) and list(z["genes"].astype(str)) == list(I.S[s]["G"])
        P[s], M[s] = z["P_down"], z["member"]
        CONS[s], PROG[s] = z["P_down_consensus"], z["P_down_programme"]
    rows = []
    rows += I.score("CLEAN4 (LOLO)", C4)
    rows += I.score("621 direction head (HGB)", P, member=M, ref=C4)
    rows += I.score("donor consensus sign (family-balanced)", CONS, ref=C4)
    rows += I.score("programme (family-balanced, gene-level)", PROG, ref=C4)
    (HERE / "work/score621.json").write_text(json.dumps(rows, indent=1, default=float))
    keys = ["surface", "source", "targets", "cov", "acc", "acc_w", "d_vs_ref", "d_ci", "auroc", "prec_at_nreal"]
    for r in rows:
        print("  ".join(f"{k}={r[k]:.4f}" if isinstance(r.get(k), float) else f"{k}={r.get(k)}" for k in keys if k in r))


if __name__ == "__main__":
    main()
