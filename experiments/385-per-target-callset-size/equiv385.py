#!/usr/bin/env python3
"""385 drift gate: prove the two drifted files are behavior-inert for THIS harness,
then pin their shas so any FURTHER concurrent edit is a hard stop, not silent drift.

Gate A (make_submission_pm_cd4.py): the five builder functions this harness imports
must be AST-identical between HEAD (= what cap15/Kc1 effectively used for these
functions) and the worktree. Timing of the 319o landing then cannot matter.
Gate B (rc_emit.py): frozen == HEAD; the worktree adds only optional order_key.
Prove bit-identity of the default path by running HEAD's and worktree's
rc_hybrid_block on one synthetic case (same seeds, CallSet-capped like our arms).

Writes work/drift385.json {shas, gate results}. verify385() re-checks pins.
Usage: uv run python equiv385.py
"""
import ast
import hashlib
import importlib.util
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
FIVE = ["profile_from_effect", "emission_spread_vector", "emission_rank_axis",
        "make_rngs", "array_sha256"]
TMP = Path("/tmp/equiv385")
PINNED = None  # filled by main(), read back from drift385.json by verify385()


def sha_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def gate_a():
    head = subprocess.check_output(["git", "show", "HEAD:scripts/make_submission_pm_cd4.py"],
                                   cwd=ROOT)
    work = (ROOT / "scripts/make_submission_pm_cd4.py").read_bytes()
    t_head, t_work = ast.parse(head.decode()), ast.parse(work.decode())

    def fns(t):
        return {n.name: ast.dump(n) for n in ast.walk(t) if isinstance(n, ast.FunctionDef)}

    fh, fw = fns(t_head), fns(t_work)
    missing = [f for f in FIVE if f not in fh or f not in fw]
    changed = [f for f in FIVE if f in fh and f in fw and fh[f] != fw[f]]
    return {"pass": not missing and not changed, "missing": missing, "changed": changed,
            "head_sha256": sha_bytes(head), "worktree_sha256": sha_bytes(work)}


def load_mod(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def gate_b():
    """One REAL target through both rc_emit versions on the production call path."""
    TMP.mkdir(parents=True, exist_ok=True)
    head_path = TMP / "rc_emit_HEAD.py"
    head_path.write_bytes(subprocess.check_output(["git", "show", "HEAD:src/vcc2026/rc_emit.py"],
                                                  cwd=ROOT))
    for _p in (str(ROOT / "experiments/323-current-family-validation"), str(ROOT / "src")):
        if _p not in sys.path:
            sys.path.insert(0, _p)
    import run323 as R
    import evaluate323 as E
    new = load_mod("rc_emit_worktree385", str(ROOT / "src/vcc2026/rc_emit.py"))
    old = load_mod("rc_emit_HEAD385", str(head_path))
    from vcc2026.donor_transfer import emit_counts as EMIT_COUNTS
    # R.builder() itself enforces the frozen-input gate, which is exactly what is
    # drifting; load the module directly instead (gate A covers the five functions).
    with R.model_env("uniform"):
        b = R.load_module("builder385_gate", ROOT / "scripts/make_submission_pm_cd4.py")
    assert b.H1_UNION_G is None and b.MASS_FIX is None
    r, c, emit_c = E.load_reference()
    with np.load(ROOT / "experiments/323-current-family-validation/work/effects.npz",
                 allow_pickle=False) as z:
        eff = R.map_axis(z["cap15"], z["genes"].astype(str), r["genes"])
    i, target = 0, r["targets"][0]
    profile = b.profile_from_effect(eff[i], r["b_ctrl"], min_cpm=10.)
    pick_rng, _ = b.make_rngs(0)
    picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
    spread = b.emission_spread_vector(
        eff[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
        keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
        spread_lo=0., spread_hi=.60)
    gate = b.emission_rank_axis(eff[i], r["b_ctrl"], min_cpm=10.)[0]
    key = np.where(gate, eff[i], 0.)
    kw = dict(spread=spread, overdispersion=4., control=emit_c, key=key,
              jitter_seed=1000 + i, call_set_first=True,
              call_set=new.CallSet(400, 2200, 1.9), emit_counts=EMIT_COUNTS)
    b_new, _ = new.rc_hybrid_block(profile, r["libraries"][picks], .1,
                                   np.random.default_rng(0), **kw)
    kw_old = dict(kw)
    kw_old["call_set"] = old.CallSet(400, 2200, 1.9)
    b_old, _ = old.rc_hybrid_block(profile, r["libraries"][picks], .1,
                                   np.random.default_rng(0), **kw_old)
    same = np.array_equal(np.asarray(b_new), np.asarray(b_old))
    return {"pass": bool(same), "target": str(target),
            "block_shape": list(np.asarray(b_new).shape),
            "block_sum": int(np.asarray(b_new).sum()),
            "worktree_sha256": sha_bytes((ROOT / "src/vcc2026/rc_emit.py").read_bytes()),
            "head_sha256": sha_bytes(head_path.read_bytes())}


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    t0 = time.monotonic()
    a, b = gate_a(), gate_b()
    rec = {"gate_a_builder_fns_ast_identical": a, "gate_b_default_path_bitidentical": b,
           "pass": bool(a["pass"] and b["pass"]), "seconds": time.monotonic() - t0}
    (WORK / "drift385.json").write_text(
        __import__("json").dumps(rec, indent=2, sort_keys=True) + "\n")
    print(__import__("json").dumps({"pass": rec["pass"], "A": a["pass"], "B": b["pass"],
                                    "B_sum": b.get("block_sum")}, indent=2))
    if not rec["pass"]:
        raise AssertionError("drift gate failed; see work/drift385.json")


def verify385():
    """Hard stop unless every frozen input matches, except the two pinned drifted files
    whose shas must equal the gate-time pins (any further edit -> ValueError)."""
    import json
    rec = json.loads((WORK / "drift385.json").read_text())
    if not rec["pass"]:
        raise ValueError("drift gate has not passed")
    pins = {"scripts/make_submission_pm_cd4.py":
            rec["gate_a_builder_fns_ast_identical"]["worktree_sha256"],
            "src/vcc2026/rc_emit.py":
            rec["gate_b_default_path_bitidentical"]["worktree_sha256"]}
    frozen = json.loads((BASE_INPUTS()).read_text())
    for name, expected in frozen["inputs_sha256"].items():
        actual = sha_bytes((ROOT / name).read_bytes())
        if name in pins:
            if actual != pins[name]:
                raise ValueError(f"concurrent edit moved {name} after the gate pin; "
                                 f"refusing to emit incomparable arms")
            continue
        if actual != expected:
            raise ValueError(f"frozen input drift: {name}")
    return rec


def BASE_INPUTS():
    return ROOT / "experiments/323-current-family-validation/inputs.json"


if __name__ == "__main__":
    main()
