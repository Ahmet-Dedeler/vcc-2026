#!/usr/bin/env python3
"""366 tests: the patched builder copy, tiny synthetic frame + known-answer row.

Runs with `uv run python experiments/366-ctx-a-prior-patch/test366.py` from the repo
root. Imports work/patched/make_submission_pm_cd4.py BY PATH (spec_from_file_location;
`__file__` pointed at the real builder so its ROOT resolves to the repo). No h5ad, no
build, peak RSS < 1 GiB.
"""
from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PATCHED = HERE / "work/patched/make_submission_pm_cd4.py"
REAL = ROOT / "scripts/make_submission_pm_cd4.py"
BASE = ROOT / "experiments/323-current-family-validation/work"

PASS: list[str] = []


def check(name: str, cond: bool) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if cond:
        PASS.append(name)
    else:
        raise AssertionError(f"test failed: {name}")


def load_patched(env: dict[str, str], modname: str):
    """Exec the patched copy under `env` (restored after). Raises whatever import raises."""
    saved = dict(os.environ)
    try:
        os.environ.update(env)
        for var in ("VCC_CTX_A_PRIOR", "VCC_ALPHA"):
            if var not in env and var in os.environ:
                del os.environ[var]
        spec = importlib.util.spec_from_file_location(modname, PATCHED)
        mod = importlib.util.module_from_spec(spec)
        mod.__file__ = str(REAL)  # ROOT must resolve to the repo, not work/patched/
        sys.modules[modname] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        os.environ.clear()
        os.environ.update(saved)


def write_npz(path: Path, **kw) -> Path:
    np.savez(path, **kw)
    return path


def frame(tmp: Path):
    """3 contexts x 5 targets x 50 genes synthetic frame; genes G000-G002 are 'panel'."""
    rng = np.random.default_rng(366)
    genes = np.array([f"G{i:03d}" for i in range(50)])
    design = rng.normal(0.0, 0.5, (5, 50))
    b_ctrl = {c: np.abs(rng.normal(0.5, 0.3, 50)) + 0.01 for c in ("A", "B", "C")}
    panel = np.array([0, 1, 2], dtype=np.int64)
    return genes, design, b_ctrl, panel


def expected_block(row, ctrl, hk, bk):
    head = np.argsort(-np.abs(row), kind="stable")[:hk]
    mask = np.ones(row.size, bool)
    mask[head] = False
    rest = np.argsort(-ctrl, kind="stable")
    return head, rest[mask[rest]][:bk]


def test_known_answer() -> None:
    print("known-answer: 323 cap15 members + 361 block-size recomputation (own code)")
    m = __import__("json").loads((BASE / "metrics-cap15.json").read_text())["raw"]
    expect = {"pds": 0.88964, "mse": 0.97495, "jac": 0.07476,
              "nmae": 0.95084, "fid": 0.49755, "reach": 0.25435}
    for k, v in expect.items():
        check(f"metrics-cap15.json {k}=={v} (got {m[k]:.8f})", abs(m[k] - v) < 5e-6)
    with np.load(BASE / "effects.npz", allow_pickle=False) as e:
        cap15, egenes = e["cap15"], e["genes"].astype(str)
    with np.load(BASE / "reference.npz", allow_pickle=True) as r:
        rgenes, keep, ref_lin = r["genes"].astype(str), r["keep"].astype(bool), r["ref_lin"]
    full = np.full(len(egenes), -1.0)
    pos = {g: i for i, g in enumerate(egenes)}
    kept_genes = rgenes[keep]
    miss = sorted({str(g) for g in kept_genes} - set(pos))
    check(f"universe misses == 361's 2 absent genes (got {miss})",
          miss == ["HSPA14-1", "TBCE-1"])
    full[[pos[g] for g in kept_genes if g in pos]] = np.asarray(
        [v for g, v in zip(kept_genes, np.asarray(ref_lin, dtype=np.float64)) if g in pos])
    sizes = set()
    for t in range(cap15.shape[0]):
        _, block = expected_block(cap15[t], full, 400, 2200)
        sizes.add(block.size)
    check(f"block size exactly 2200 on all 85 targets (got {sorted(sizes)})", sizes == {2200})


def test_unset_inert(mod) -> None:
    print("unset => inert")
    check("CTX_A_PRIOR_PATH is None", mod.CTX_A_PRIOR_PATH is None)
    check("CTX_A_PRIOR is None", mod.CTX_A_PRIOR is None)
    genes, design, b_ctrl, panel = frame(Path("."))
    out, audit = mod.build_ctx_a_context_effects(
        design, None, b_ctrl["A"], genes, panel, ("A", "B", "C"))
    check("helper returns (None, None) when off", out is None and audit is None)


def test_sign_block(mod, tmp: Path) -> None:
    print("sign_block on synthetic frame (head_k=5, block_k=10)")
    genes, design, b_ctrl, panel = frame(tmp)
    prior_genes = genes[3:]  # never touches panel G000-G002
    rng = np.random.default_rng(7)
    vec = rng.choice([-1.0, 0.0, 1.0], size=47, p=[0.45, 0.1, 0.45])
    path = write_npz(tmp / "sign.npz", genes=np.asarray(prior_genes),
                     mode=np.asarray("sign_block"), vector=vec,
                     tau=np.asarray(0.0), alpha=np.asarray(1.0),
                     scale_rule=np.asarray("gate_median"), context=np.asarray("A"))
    mod2 = load_patched({"VCC_CTX_A_PRIOR": str(path)}, "builder366_sign")
    prior = mod2.CTX_A_PRIOR
    check("mode/context parsed", prior["mode"] == "sign_block" and prior["context"] == "A")
    out, audit = mod2.apply_ctx_a_prior(design, prior, b_ctrl["A"], genes, panel,
                                        head_k=5, block_k=10)
    sgn = np.zeros(50)
    sgn[3:] = np.where(vec != 0.0, np.sign(vec), 0.0)
    exp = design.copy()
    n_use = 0
    for t in range(5):
        head, block = expected_block(design[t], b_ctrl["A"], 5, 10)
        check(f"row {t}: head entries bitwise identical",
              bool(np.array_equal(out[t, head], design[t, head])))
        use = block[sgn[block] != 0.0]
        exp[t, use] = sgn[use] * np.abs(design[t, use])
        n_use += use.size
    check("block entries == sign*|design| exactly", bool(np.array_equal(out, exp)))
    check(f"audit touched == {n_use}", audit["touched_entries"] == n_use)
    ctx, _ = mod2.build_ctx_a_context_effects(design, prior, b_ctrl["A"], genes, panel,
                                              ("A", "B", "C"), head_k=5, block_k=10)
    check("B rows are the SAME object (untouched)", ctx["B"] is design)
    check("C rows bitwise identical", bool(np.array_equal(ctx["C"], design)))
    check("A rows differ", not bool(np.array_equal(ctx["A"], design)))


def test_sign_block_pdown(mod, tmp: Path) -> None:
    print("sign_block with p_down/tau")
    genes, design, b_ctrl, panel = frame(tmp)
    prior_genes = genes[3:]
    rng = np.random.default_rng(11)
    down = rng.uniform(0, 1, size=47)
    tau = 0.1
    path = write_npz(tmp / "pdown.npz", genes=np.asarray(prior_genes),
                     mode=np.asarray("sign_block"), vector=np.zeros(47),
                     p_down=down, tau=np.asarray(tau), context=np.asarray("A"))
    mod2 = load_patched({"VCC_CTX_A_PRIOR": str(path)}, "builder366_pdown")
    out, _ = mod2.apply_ctx_a_prior(design, mod2.CTX_A_PRIOR, b_ctrl["A"], genes, panel,
                                    head_k=5, block_k=10)
    sgn = np.zeros(50)
    sgn[3:] = np.where(np.abs(down - 0.5) > tau, np.sign(down - 0.5), -1.0)
    exp = design.copy()
    for t in range(5):
        _, block = expected_block(design[t], b_ctrl["A"], 5, 10)
        use = block[sgn[block] != 0.0]
        exp[t, use] = sgn[use] * np.abs(design[t, use])
    check("p_down tau-rule exact (else down)", bool(np.array_equal(out, exp)))


def test_add_vector(mod, tmp: Path) -> None:
    print("add_vector on synthetic frame")
    genes, design, b_ctrl, panel = frame(tmp)
    prior_genes = genes[3:]
    rng = np.random.default_rng(13)
    vec = rng.normal(0, 1, size=47)
    vec[np.abs(vec) < 0.2] = 0.0
    alpha = 0.7
    path = write_npz(tmp / "add.npz", genes=np.asarray(prior_genes),
                     mode=np.asarray("add_vector"), vector=vec,
                     alpha=np.asarray(alpha), context=np.asarray("A"))
    mod2 = load_patched({"VCC_CTX_A_PRIOR": str(path)}, "builder366_add")
    out, audit = mod2.apply_ctx_a_prior(design, mod2.CTX_A_PRIOR, b_ctrl["A"], genes, panel,
                                        head_k=5, block_k=10)
    full = np.zeros(50)
    full[3:] = vec
    univ = b_ctrl["A"] > 0
    exp = design.copy()
    for t in range(5):
        s = float(np.median(np.abs(design[t, univ]))) / float(np.median(np.abs(full[univ])))
        exp[t] += alpha * s * full
    check("add_vector == alpha*s*vector exactly", bool(np.array_equal(out, exp)))
    check("scale audit present", audit["scale_median"] > 0)
    ctx, _ = mod2.build_ctx_a_context_effects(design, mod2.CTX_A_PRIOR, b_ctrl["A"], genes,
                                              panel, ("A", "B", "C"), head_k=5, block_k=10)
    check("B/C bitwise identical", bool(np.array_equal(ctx["B"], design))
          and bool(np.array_equal(ctx["C"], design)))


def test_malformed(mod, tmp: Path) -> None:
    print("malformed artifacts refused at load")
    genes = np.array([f"G{i:03d}" for i in range(50)])
    good = dict(genes=genes[3:], mode=np.asarray("sign_block"), vector=np.ones(47))
    cases = {
        "missing vector": {k: v for k, v in good.items() if k != "vector"},
        "NaN vector": {**good, "vector": np.full(47, np.nan)},
        "bad mode": {**good, "mode": np.asarray("replace")},
        "bad context": {**good, "context": np.asarray("Z")},
        "extra key": {**good, "bogus": np.asarray(1)},
        "p_down range": {**good, "p_down": np.full(47, 1.5)},
        "tau range": {**good, "tau": np.asarray(0.9)},
        "alpha zero": {**good, "mode": np.asarray("add_vector"), "alpha": np.asarray(0.0)},
        "shape mismatch": {**good, "vector": np.ones(10)},
    }
    for name, kw in cases.items():
        p = write_npz(tmp / f"bad-{len(PASS)}.npz", **kw)
        try:
            mod._load_ctx_a_prior(str(p))
            check(f"{name} refused", False)
        except ValueError:
            check(f"{name} refused", True)
    try:
        mod._load_ctx_a_prior(str(tmp / "does-not-exist.npz"))
        check("missing path refused", False)
    except ValueError:
        check("missing path refused", True)


def test_import_refusals(tmp: Path) -> None:
    print("import-time refusals")
    genes = np.array([f"G{i:03d}" for i in range(50)])
    path = write_npz(tmp / "ok.npz", genes=genes[3:], mode=np.asarray("sign_block"),
                     vector=np.ones(47), context=np.asarray("A"))
    try:
        load_patched({"VCC_CTX_A_PRIOR": str(path), "VCC_ALPHA": "knee"}, "builder366_knee")
        check("knee combo refused", False)
    except ValueError:
        check("knee combo refused", True)
    bad = write_npz(tmp / "badmode.npz", genes=genes[3:], mode=np.asarray("nope"),
                    vector=np.ones(47))
    try:
        load_patched({"VCC_CTX_A_PRIOR": str(bad)}, "builder366_badmode")
        check("bad mode refused at import", False)
    except ValueError:
        check("bad mode refused at import", True)


def test_apply_refusals(mod, tmp: Path) -> None:
    print("apply-time refusals")
    genes, design, b_ctrl, panel = frame(tmp)
    vec = np.ones(47)
    p = write_npz(tmp / "panel.npz", genes=np.array(["G000"] + [f"G{i:03d}" for i in range(3, 49)]),
                  mode=np.asarray("sign_block"), vector=vec, context=np.asarray("A"))
    prior = mod._load_ctx_a_prior(str(p))
    try:
        mod.apply_ctx_a_prior(design, prior, b_ctrl["A"], genes, panel, head_k=5, block_k=10)
        check("panel support refused", False)
    except ValueError:
        check("panel support refused", True)
    p2 = write_npz(tmp / "miss.npz", genes=np.array(["NOT_A_GENE"] * 47),
                   mode=np.asarray("sign_block"), vector=vec, context=np.asarray("A"))
    prior2 = mod._load_ctx_a_prior(str(p2))
    try:
        mod.apply_ctx_a_prior(design, prior2, b_ctrl["A"], genes, panel, head_k=5, block_k=10)
        check("gene miss refused", False)
    except ValueError:
        check("gene miss refused", True)
    p3 = write_npz(tmp / "zero.npz", genes=genes[3:], mode=np.asarray("add_vector"),
                   vector=np.zeros(47), context=np.asarray("A"))
    prior3 = mod._load_ctx_a_prior(str(p3))
    try:
        mod.apply_ctx_a_prior(design, prior3, b_ctrl["A"], genes, panel, head_k=5, block_k=10)
        check("vacuous add_vector refused", False)
    except ValueError:
        check("vacuous add_vector refused", True)
    p4 = write_npz(tmp / "ctxb.npz", genes=genes[3:], mode=np.asarray("sign_block"),
                   vector=np.ones(47), context=np.asarray("B"))
    prior4 = mod._load_ctx_a_prior(str(p4))
    try:
        mod.build_ctx_a_context_effects(design, prior4, b_ctrl["A"], genes, panel, ("A", "C"))
        check("wrong-context refused", False)
    except ValueError:
        check("wrong-context refused", True)


def test_scripts_untouched() -> None:
    print("scripts/ untouched")
    check("real builder has no CTX_A_PRIOR", "CTX_A_PRIOR" not in REAL.read_text())
    check("patched copy exists", PATCHED.is_file())
    src = PATCHED.read_text()
    check("patched copy carries all 4 anchors",
          all(s in src for s in ("_load_ctx_a_prior", "apply_ctx_a_prior",
                                 "build_ctx_a_context_effects", "ctx_a_prior_mode")))


def main() -> None:
    assert PATCHED.is_file(), "run patch366.py first (work/patched copy missing)"
    with tempfile.TemporaryDirectory(prefix="t366-") as td:
        tmp = Path(td)
        test_known_answer()
        mod = load_patched({}, "builder366_base")
        test_unset_inert(mod)
        test_sign_block(mod, tmp)
        test_sign_block_pdown(mod, tmp)
        test_add_vector(mod, tmp)
        test_malformed(mod, tmp)
        test_import_refusals(tmp)
        test_apply_refusals(mod, tmp)
        test_scripts_untouched()
    print(f"\n{len(PASS)} checks passed, 0 failed")


if __name__ == "__main__":
    main()
