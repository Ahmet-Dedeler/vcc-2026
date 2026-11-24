#!/usr/bin/env python3
"""Execute the 2026-09-28 abc_targets block in scripts/def_postbuild.py on a temp ROOT.

The block is the source lines themselves (not a paraphrase). ROOT in that exec is a temp
tree, so data/challenge/pert_counts.csv is whatever this script writes. The ABC_* constants
are the module's. Exits 0 when the current source shows the fail-closed October cases below.
"""
from __future__ import annotations

import hashlib
import importlib.util
import sys
import tempfile
import textwrap
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
SRC_LINES = (REPO / "scripts/def_postbuild.py").read_text().splitlines()
BLOCK = textwrap.dedent("\n".join(SRC_LINES[317:329]))  # lines 318-329
BIND = textwrap.dedent("\n".join(SRC_LINES[283:287]))  # lines 284-287, the symbol-axis binding row

spec = importlib.util.spec_from_file_location("def_postbuild", REPO / "scripts/def_postbuild.py")
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

import pandas as pd

ABC_AXIS = mod.ABC_AXIS_SHA
ABC_COLS = mod.ABC_OPERATOR_COLUMNS_SHA
PIN = hashlib.sha256(
    "\n".join(sorted(str(t) for t in pd.read_csv(REPO / "data/challenge/pert_counts.csv")["target_gene"])).encode()
).hexdigest()
MOVED = "0" * 64  # a digest that is not A/B/C's
assert PIN != ABC_AXIS and MOVED != ABC_AXIS


def write_perts(root: Path, names: list[str]) -> None:
    d = root / "data" / "challenge"
    d.mkdir(parents=True, exist_ok=True)
    (d / "pert_counts.csv").write_text("target_gene\n" + "\n".join(names) + "\n")


def run_block(root: Path, panel: list[str], new_panel: bool, expect_axis: str, columns_sha: str):
    checks: list[tuple[str, bool, str]] = []

    def chk(name, ok, value=""):
        checks.append((name, bool(ok), str(value)))

    ns = dict(
        panel_set=set(panel), new_panel=new_panel, expect_axis_sha=expect_axis,
        op={"columns_sha256": columns_sha}, chk=chk, pd=pd, ROOT=root,
        ABC_OPERATOR_COLUMNS_SHA=ABC_COLS, ABC_AXIS_SHA=ABC_AXIS,
    )
    exec(BLOCK, ns)
    return checks, ns["abc_targets"]


def bind(recorded: str, expect: str) -> bool:
    checks: list[tuple[str, bool, str]] = []

    def chk(name, ok, value=""):
        checks.append((name, bool(ok), str(value)))

    exec(BIND, dict(chk=chk, recorded_axis=recorded, expect_axis_sha=expect, op={}))
    return checks[0][1]


def fixed_abc(panel: list[str]) -> bool:
    """The pin in patches/def_postbuild_abc_targets.patch. Not applied to the source."""
    digest = hashlib.sha256("\n".join(sorted(panel)).encode()).hexdigest()
    return digest == PIN


def fails(checks) -> list[str]:
    return [n for n, ok, _ in checks if not ok]


def main() -> int:
    print("BLOCK source line 318:")
    print(" ", SRC_LINES[317])
    abc = [str(t) for t in pd.read_csv(REPO / "data/challenge/pert_counts.csv")["target_gene"]]
    new = ["NEW%03d" % i for i in range(300)]
    # partial: drop one ABC target that we treat as outside the donor union, add one new gene
    partial = abc[1:] + ["NOT_IN_UNION"]
    bad = 0

    def expect(label, cond, detail):
        nonlocal bad
        print(f"{'OK' if cond else 'UNEXPECTED'}  {label}: {detail}")
        if not cond:
            bad += 1

    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        # 1. Oct 22: --bundle data/challenge after the download. The file IS the new panel.
        write_perts(root, new)
        checks, abc_t = run_block(root, new, True, MOVED, MOVED)
        expect("oct22 tautology abc_targets", abc_t is True, f"abc_targets={abc_t}")
        expect("oct22 false-FAIL on a moved axis", any("own 300 targets" in n for n in fails(checks)),
               fails(checks))
        expect("oct22 pin would not call them A/B/C targets", fixed_abc(new) is False, fixed_abc(new))

        # 2. Before the download: data/challenge still A/B/C, bundle has new targets, axis moved.
        write_perts(root, abc)
        checks, abc_t = run_block(root, new, True, MOVED, MOVED)
        expect("pre-download new targets", abc_t is False and not fails(checks),
               f"abc_targets={abc_t} fails={fails(checks)}")

        # 3. New contexts, A/B/C's own targets, axis is A/B/C's (497dr's rehearsal). File still A/B/C.
        checks, abc_t = run_block(root, abc, True, ABC_AXIS, ABC_COLS)
        expect("same targets, axis is A/B/C", abc_t is True and not fails(checks),
               f"abc_targets={abc_t} fails={fails(checks)}")

        # 4. A/B/C itself.
        checks, abc_t = run_block(root, abc, False, ABC_AXIS, ABC_COLS)
        expect("A/B/C known-answer row passes", any(ok and "known-answer" in n for n, ok, _ in checks),
               [(n, ok) for n, ok, _ in checks])

        # 5. Side bundle keeps A/B/C targets, but data/challenge was already overwritten.
        write_perts(root, new)
        checks, abc_t = run_block(root, abc, True, ABC_AXIS, ABC_COLS)
        expect("preserved A/B/C bundle after download is not abc_targets", abc_t is False, abc_t)
        expect("and the unchanged axis is declared impossible",
               any("SYMBOL axis" in n for n in fails(checks)), fails(checks))
        expect("pin still recognises the real A/B/C set", fixed_abc(abc) is True, PIN)

        # 6. Partial overlap whose swapped genes miss the donor union: axis digest stays A/B/C's.
        #    Set algebra, no donor tables: axis = union - panel, so the digests match iff
        #    (symdiff ∩ union) is empty.
        union = {"IN_UNION"}
        abc_set, part_set = set(abc), set(partial)
        sym = abc_set.symmetric_difference(part_set)
        axis_same = (union - abc_set) == (union - part_set)
        expect("partial symdiff misses a union that lacks those genes",
               axis_same and sym.isdisjoint(union) and not abc_set == part_set,
               f"symdiff={sorted(sym)[:4]}… n={len(sym)}")
        write_perts(root, abc)  # file is still the official list; panel is the partial set
        checks, abc_t = run_block(root, partial, True, ABC_AXIS, ABC_COLS)
        expect("partial overlap false-FAIL (axis unchanged)",
               abc_t is False and any("impossible" in v for _, ok, v in checks if not ok),
               fails(checks))

        # 7. Same partial panel but the axis DID move: this block adds no failure.
        checks, abc_t = run_block(root, partial, True, MOVED, MOVED)
        expect("partial overlap, axis moved: block adds no row", not fails(checks), fails(checks))

        # 8. Fail-open attempt: targets differ, recomputed axis moved, recorded axis is A/B/C's.
        write_perts(root, abc)
        checks, abc_t = run_block(root, new, True, MOVED, MOVED)
        bound = bind(ABC_AXIS, MOVED)
        expect("wrong recorded axis does not get a passing moved-off row", not fails(checks), fails(checks))
        expect("line 285 still fails that operator", bound is False, bound)
        expect("line 285 passes when the recorded axis is this bundle's", bind(MOVED, MOVED) is True, True)

    print("PIN", PIN)
    print("FAILS" if bad else "ALL DEMONSTRATED", bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
