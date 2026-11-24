#!/usr/bin/env python3
"""511's tie rules vs the previous function (39d38196^).

The scalar-q path never calls this function. On the per-target path, a lower-cut tie
block moves: stercile used to leave it to np.digitize; it now goes LOW. A NaN used to
be filled with the median and included in the quantile; it is now excluded from the
cuts and assigned MID. Exits 0.
"""
from __future__ import annotations

import ast
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
SRC = (ROOT / "scripts/make_submission_pm_cd4.py").read_text()


def extract(text: str, name: str) -> str:
    tree = ast.parse(text)
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == name:
            return ast.get_source_segment(text, node)
    raise SystemExit(f"{name} not found")


def load_old() -> str:
    return subprocess.check_output(
        ["git", "show", "39d38196^:scripts/make_submission_pm_cd4.py"], cwd=ROOT, text=True)


def fn(source: str):
    ns: dict = {"np": np}
    exec(source, ns)
    return ns["sharpen_q_by_tercile"]


def bundle(norms: list[float]):
    # stercile reads the carrier norm on non-panel genes. One non-panel column, so the
    # norm equals the absolute value we put there. No panel genes in `genes`.
    n = len(norms)
    carriers = np.array(norms, dtype=np.float64).reshape(n, 1)
    bulk = carriers.copy()
    genes = np.array(["G0"])
    perts = np.array([f"T{i}" for i in range(n)])
    return bulk, carriers, genes, perts


def main() -> int:
    new = fn(extract(SRC, "sharpen_q_by_tercile"))
    old = fn(extract(load_old(), "sharpen_q_by_tercile"))
    pt = {"mode": "stercile", "q": [2.5, 2.0, 1.75]}
    # The lower cut lands on a repeated 1.0. np.digitize (right=False) puts that edge in
    # MID. 511 sends a lower-cut tie block of size >= 2 to LOW.
    vals = [1.0, 1.0, 1.0, 2.0, 3.0, 4.0]
    b, c, g, p = bundle(vals)
    q_old, a_old = old(b, c, g, p, pt)
    q_new, a_new = new(b, c, g, p, pt)
    print("stercile tie block")
    print("  old q", q_old.tolist(), "counts", a_old["counts"])
    print("  new q", q_new.tolist(), "counts", a_new["counts"])
    assert q_old.tolist() != q_new.tolist()
    assert q_new[:3].tolist() == [2.5, 2.5, 2.5]   # the tied weakest mass is LOW
    assert q_old[:3].tolist() == [2.0, 2.0, 2.0]   # the previous function left the edge in MID

    # No tie, no NaN: the two functions agree. Audit keys are the same set.
    vals = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    b, c, g, p = bundle(vals)
    q_old, a_old = old(b, c, g, p, pt)
    q_new, a_new = new(b, c, g, p, pt)
    print("stercile no tie", q_old.tolist())
    assert q_old.tolist() == q_new.tolist()
    assert set(a_old) == set(a_new) == {"mode", "q", "edges", "counts", "floor_to_low", "tie_to_mid", "nonfinite"}

    # A NaN moves the quantile under the old rule (it is filled, then included) and is MID under the new one.
    vals = [1.0, 2.0, 3.0, 4.0, 100.0, float("nan")]
    b, c, g, p = bundle(vals)
    q_old, a_old = old(b, c, g, p, pt)
    q_new, a_new = new(b, c, g, p, pt)
    print("stercile with NaN")
    print("  old edges", [round(e, 4) for e in a_old["edges"]], "q", q_old.tolist())
    print("  new edges", [round(e, 4) for e in a_new["edges"]], "q", q_new.tolist(), "nonfinite", a_new["nonfinite"])
    assert a_new["nonfinite"] == 1
    assert float(q_new[-1]) == 2.0   # the NaN sits in MID and is not a cut point
    assert [round(e, 6) for e in a_old["edges"]] != [round(e, 6) for e in a_new["edges"]]

    # The scalar path does not call the function.
    assert "if BULK_SHARPEN_PT is not None else None" in SRC
    print("SHARPEN: unset/scalar path does not call sharpen_q_by_tercile; the ON path reassigns lower-cut ties")
    return 0


if __name__ == "__main__":
    sys.exit(main())
