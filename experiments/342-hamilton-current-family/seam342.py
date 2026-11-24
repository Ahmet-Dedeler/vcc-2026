"""Explicit AST seam on production `emit_counts`: one statement, in one branch.

The production function is parsed from `src/vcc2026/donor_transfer.py` at import time,
and inside its `mode == "quasi"` branch the single statement

    totals = np.floor(n * L * profile_p).astype(np.int64)

is replaced by

    totals = __vcc342_round_totals__(n * L * profile_p)

Nothing else is edited: the other three `totals = np.floor(...)` sites (`resample_mix`,
`mixed`, `even`) stay in place, `_quasi_block`, `_repair_totals` and
`_stochastic_round` are the production objects, and the production module itself is
never modified (no monkeypatch, no global state). `attestation()` records exactly what
differs so a freeze can hash-bind it; `emit_counts_with(rounder)` binds a per-block
rounder into a fresh function object whose globals carry that one name.

If the production source changes so that the statement is not found exactly once, this
module refuses to build rather than guessing.
"""
from __future__ import annotations

import ast
import builtins
import hashlib
import inspect
import types
from pathlib import Path

import numpy as np
from vcc2026 import donor_transfer as production

HOOK = "__vcc342_round_totals__"
FUNCTION = "emit_counts"
BRANCH = "mode == 'quasi'"
ORIGINAL = "totals = np.floor(n * L * profile_p).astype(np.int64)"
SEAMED = f"totals = {HOOK}(n * L * profile_p)"
EXPECTED_UNTOUCHED_FLOOR_SITES = 3


def production_path() -> Path:
    return Path(inspect.getsourcefile(production)).resolve()


def _unbound(_x):
    raise RuntimeError("seamed emit_counts called without a bound rounder")


def _function(tree: ast.Module) -> ast.FunctionDef:
    found = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == FUNCTION]
    if len(found) != 1:
        raise ValueError(f"expected exactly one {FUNCTION} definition, found {len(found)}")
    return found[0]


def _floor_sites(fn: ast.FunctionDef) -> int:
    count = 0
    for node in ast.walk(fn):
        if (isinstance(node, ast.Assign) and len(node.targets) == 1
                and ast.unparse(node.targets[0]) == "totals"
                and ast.unparse(node.value).startswith("np.floor(")):
            count += 1
    return count


class Seam:
    """The seamed production function plus its attestation."""

    def __init__(self):
        path = production_path()
        source_bytes = path.read_bytes()
        tree = ast.parse(source_bytes.decode())
        fn = _function(tree)
        original_text = ast.unparse(fn)
        floor_sites_before = _floor_sites(fn)
        branches = [n for n in fn.body if isinstance(n, ast.If) and ast.unparse(n.test) == BRANCH]
        if len(branches) != 1:
            raise ValueError(f"expected exactly one `{BRANCH}` branch, found {len(branches)}")
        branch = branches[0]
        replaced = 0
        for i, stmt in enumerate(branch.body):
            if isinstance(stmt, ast.Assign) and ast.unparse(stmt) == ORIGINAL:
                new = ast.parse(SEAMED).body[0]
                ast.copy_location(new, stmt)
                branch.body[i] = new
                replaced += 1
        if replaced != 1:
            raise ValueError(f"expected exactly one rounding statement in the quasi branch, "
                             f"found {replaced}")
        module = ast.Module(body=[fn], type_ignores=[])
        ast.fix_missing_locations(module)
        seamed_text = ast.unparse(fn)
        before, after = original_text.splitlines(), seamed_text.splitlines()
        if len(before) != len(after):
            raise ValueError("the seam changed the line count of emit_counts")
        # unparse indents statements; the attestation records the statements themselves
        differences = [[a.strip(), b.strip()] for a, b in zip(before, after) if a != b]
        if differences != [[ORIGINAL, SEAMED]]:
            raise ValueError(f"unexpected seam differences: {differences}")
        floor_sites_after = _floor_sites(fn)
        if floor_sites_after != EXPECTED_UNTOUCHED_FLOOR_SITES or \
                floor_sites_before != EXPECTED_UNTOUCHED_FLOOR_SITES + 1:
            raise ValueError("the other rounding sites are not what this seam expects")
        namespace = {
            "__builtins__": builtins,
            "np": np,
            "_quasi_block": production._quasi_block,
            "_repair_totals": production._repair_totals,
            "_stochastic_round": production._stochastic_round,
            HOOK: _unbound,
        }
        code = compile(module, filename=f"{path}#vcc342-seam", mode="exec")
        exec(code, namespace)
        self._base = namespace[FUNCTION]
        repo = Path(__file__).resolve().parents[2]
        try:
            shown = str(path.relative_to(repo))
        except ValueError:
            shown = str(path)
        self.attestation = {
            "production_path": shown,
            "production_sha256": hashlib.sha256(source_bytes).hexdigest(),
            "function": FUNCTION,
            "branch": BRANCH,
            "hook": HOOK,
            "original_statement": ORIGINAL,
            "seamed_statement": SEAMED,
            "replacements": replaced,
            "line_differences": differences,
            "n_line_differences": len(differences),
            "floor_sites_before": floor_sites_before,
            "untouched_floor_sites_after": floor_sites_after,
            "borrowed_production_objects": ["_quasi_block", "_repair_totals",
                                            "_stochastic_round"],
            "original_function_sha256": hashlib.sha256(original_text.encode()).hexdigest(),
            "seamed_function_sha256": hashlib.sha256(seamed_text.encode()).hexdigest(),
        }

    def emit_counts_with(self, rounder):
        """A fresh function object with `rounder` bound as the ONLY totals hook."""
        if not callable(rounder):
            raise TypeError("rounder must be callable")
        namespace = dict(self._base.__globals__)
        namespace[HOOK] = rounder
        return types.FunctionType(self._base.__code__, namespace, self._base.__name__,
                                  self._base.__defaults__, self._base.__closure__)


def attestation() -> dict:
    return Seam().attestation
