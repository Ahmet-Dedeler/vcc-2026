#!/usr/bin/env python3
"""373 tests: 366's 45 re-run on the combined copy + new `blend_rows` tests.

Runs with `uv run python experiments/373-blend-rows-patch/test373.py` from the repo
root. Part 1 imports 366's test module BY PATH and overrides its PATCHED to the
combined copy (no write into 366's dir). Part 2 exercises `blend_rows` on a tiny
synthetic frame with independently recomputed expectations. No h5ad, no build,
peak RSS < 1 GiB.
"""
from __future__ import annotations

import hashlib
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
T366 = ROOT / "experiments/366-ctx-a-prior-patch/test366.py"

PASS: list[str] = []


def check(name: str, cond: bool) -> None:
    print(f"  {'PASS' if cond else 'FAIL'}  {name}")
    if cond:
        PASS.append(name)
    else:
        raise AssertionError(f"test failed: {name}")


def load_patched(env: dict[str, str], modname: str):
    """Exec the COMBINED patched copy under `env` (restored after)."""
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


def bframe():
    """5 rows (T000-T004) x 50 genes (G000-G049); gene positions 0,1,2 are 'panel'."""
    rng = np.random.default_rng(373)
    genes = np.array([f"G{i:03d}" for i in range(50)])
    panel_targets = [f"T{i:03d}" for i in range(5)]
    design = rng.normal(0.0, 0.5, (5, 50))
    b_ctrl = np.abs(rng.normal(0.5, 0.3, 50)) + 0.01
    panel = np.array([0, 1, 2], dtype=np.int64)
    return genes, panel_targets, design, b_ctrl, panel


def blend_artifact(tmp: Path, name: str, beta=0.35, targets=("T001", "T003"),
                   genes=None, seed=21, **extra):
    genes_all = np.array([f"G{i:03d}" for i in range(50)])
    cov = genes_all[10:30] if genes is None else np.asarray(genes)
    rng = np.random.default_rng(seed)
    mat = rng.normal(0.0, 1.0, (len(targets), len(cov)))
    kw = dict(targets=np.asarray(list(targets)), genes=np.asarray(cov),
              matrix=mat, beta=np.asarray(beta), mode=np.asarray("blend_rows"),
              context=np.asarray("A"), norm_rule=np.asarray("fixed_norm_on_covered"))
    kw.update(extra)
    return write_npz(tmp / name, **kw), cov, mat


def expected_blend(design, art_targets, panel_targets, cov_genes, all_genes, mat, beta):
    pos = {g: i for i, g in enumerate([str(x) for x in all_genes])}
    at = np.array([pos[str(g)] for g in cov_genes])
    row_of = {t: i for i, t in enumerate(panel_targets)}
    exp = design.copy()
    for k, t in enumerate(art_targets):
        i = row_of[t]
        r, m = design[i, at], mat[k]
        rn, mn = float(np.linalg.norm(r)), float(np.linalg.norm(m))
        ms = m * (rn / mn)
        bl = (1.0 - beta) * r + beta * ms
        exp[i, at] = bl * (rn / float(np.linalg.norm(bl)))
    return exp, at


def test_blend_exact(mod, tmp: Path) -> None:
    print("blend_rows exact vs independent recomputation")
    genes, pt, design, b_ctrl, panel = bframe()
    path, cov, mat = blend_artifact(tmp, "blend.npz")
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    mod2 = load_patched({"VCC_CTX_A_PRIOR": str(path)}, "builder373_blend")
    prior = mod2.CTX_A_PRIOR
    check("loader keys (targets/matrix/beta/norm_rule)", prior["mode"] == "blend_rows"
          and prior["beta"] == 0.35 and prior["norm_rule"] == "fixed_norm_on_covered"
          and prior["targets"].tolist() == ["T001", "T003"] and prior["matrix"].shape == (2, 20))
    check("loader sha", prior["sha256"] == sha)
    out, audit = mod2.apply_ctx_a_blend_rows(design, prior, genes, panel, pt)
    exp, at = expected_blend(design, ["T001", "T003"], pt, cov, genes, mat, 0.35)
    check("full matrix bitwise identical to independent recomputation", bool(np.array_equal(out, exp)))
    check("audit rows_touched==2, covered==20", audit["rows_touched"] == 2 and audit["covered_genes"] == 20)
    check("audit beta/sha/context", audit["beta"] == 0.35 and audit["artifact_sha256"] == sha
          and audit["context"] == "A" and audit["mode"] == "blend_rows")
    cos = []
    for k, t in enumerate(["T001", "T003"]):
        r, m = design[pt.index(t)][at], mat[k]
        cos.append(float(r @ m / (np.linalg.norm(r) * np.linalg.norm(m))))
    check(f"audit mean_cosine=={np.mean(cos):.6f}", abs(audit["mean_cosine"] - float(np.mean(cos))) < 1e-15)
    worst = 0.0
    for t in ["T001", "T003"]:
        r = design[pt.index(t)][at]
        worst = max(worst, abs(float(np.linalg.norm(out[pt.index(t)][at])) / float(np.linalg.norm(r)) - 1.0))
    check(f"fixed-norm identity worst rel err {worst:.2e} < 1e-12", worst < 1e-12)
    for t in ["T000", "T002", "T004"]:
        check(f"untouched row {t} bitwise", bool(np.array_equal(out[pt.index(t)], design[pt.index(t)])))
    check("genes outside C bitwise on touched rows",
          bool(np.array_equal(np.delete(out[1], at), np.delete(design[1], at)))
          and bool(np.array_equal(np.delete(out[3], at), np.delete(design[3], at))))
    ctx, _ = mod2.build_ctx_a_context_effects(design, prior, b_ctrl, genes, panel,
                                              ("A", "B", "C"), panel_targets=pt)
    check("B is the SAME object (untouched)", ctx["B"] is design)
    check("C bitwise identical", bool(np.array_equal(ctx["C"], design)))
    check("A differs", not bool(np.array_equal(ctx["A"], design)))
    check("A equals direct applicator output", bool(np.array_equal(ctx["A"], out)))


def test_beta_edges(mod, tmp: Path) -> None:
    print("beta edges: 0 => identity (hand-built prior), 1 => pure direction")
    genes, pt, design, b_ctrl, panel = bframe()
    _, cov, mat = blend_artifact(tmp, "edge.npz")
    pos = {g: i for i, g in enumerate([str(x) for x in genes])}
    at = np.array([pos[str(g)] for g in cov])
    zero_prior = {"targets": np.asarray(["T001", "T003"]), "genes": np.asarray(cov),
                  "matrix": mat, "beta": 0.0, "mode": "blend_rows",
                  "context": "A", "norm_rule": "fixed_norm_on_covered", "sha256": "x"}
    out0, _ = mod.apply_ctx_a_blend_rows(design, zero_prior, genes, panel, pt)
    check("beta=0 hand-built prior bitwise identical", bool(np.array_equal(out0, design)))
    one_prior = dict(zero_prior, beta=1.0)
    out1, _ = mod.apply_ctx_a_blend_rows(design, one_prior, genes, panel, pt)
    for k, t in enumerate(["T001", "T003"]):
        r, m = design[pt.index(t)][at], mat[k]
        got, want_dir = out1[pt.index(t)][at], m / float(np.linalg.norm(m))
        c = float(got @ want_dir / (np.linalg.norm(got) * np.linalg.norm(want_dir)))
        check(f"beta=1 row {t}: direction==artifact cos={c:.12f}", abs(c - 1.0) < 1e-12)
        check(f"beta=1 row {t}: norm preserved", abs(float(np.linalg.norm(got)) / float(np.linalg.norm(r)) - 1.0) < 1e-12)


def test_blend_load_refusals(mod, tmp: Path) -> None:
    print("blend_rows malformed artifacts refused at load")
    genes_all = np.array([f"G{i:03d}" for i in range(50)])
    cov = genes_all[10:30]
    rng = np.random.default_rng(5)
    mat = rng.normal(0, 1, (2, 20))
    good = dict(targets=np.asarray(["T001", "T003"]), genes=np.asarray(cov), matrix=mat,
                beta=np.asarray(0.35), mode=np.asarray("blend_rows"), context=np.asarray("A"),
                norm_rule=np.asarray("fixed_norm_on_covered"))
    cases = {
        "p_down here": {**good, "p_down": np.zeros(20)},
        "vector here": {**good, "vector": np.zeros(20)},
        "tau here": {**good, "tau": np.asarray(0.1)},
        "alpha here": {**good, "alpha": np.asarray(1.0)},
        "scale_rule here": {**good, "scale_rule": np.asarray("gate_median")},
        "NaN matrix": {**good, "matrix": np.where(np.arange(40).reshape(2, 20) == 0, np.nan, mat)},
        "inf matrix": {**good, "matrix": np.where(np.arange(40).reshape(2, 20) == 0, np.inf, mat)},
        "beta 0": {**good, "beta": np.asarray(0.0)},
        "beta 1.5": {**good, "beta": np.asarray(1.5)},
        "beta NaN": {**good, "beta": np.asarray(np.nan)},
        "dup targets": {**good, "targets": np.asarray(["T001", "T001"])},
        "dup genes": {**good, "genes": np.concatenate([cov[:1], cov[:1], cov[2:]]),
                      "matrix": rng.normal(0, 1, (2, 21))},
        "shape mismatch": {**good, "matrix": np.zeros((2, 19))},
        "missing beta": {k: v for k, v in good.items() if k != "beta"},
        "extra key": {**good, "bogus": np.asarray(1)},
        "bad norm_rule": {**good, "norm_rule": np.asarray("gate_median")},
        "empty targets": {**good, "targets": np.asarray([], dtype="<U1"),
                          "matrix": np.zeros((0, 20))},
        "object targets": {**good, "targets": np.asarray(["T001", "T003"], dtype=object)},
        "sign_block+matrix": dict(genes=cov, mode=np.asarray("sign_block"),
                                  vector=np.ones(20), matrix=np.zeros((2, 20))),
        "add_vector+beta": dict(genes=cov, mode=np.asarray("add_vector"),
                                vector=np.ones(20), beta=np.asarray(0.5)),
    }
    for name, kw in cases.items():
        p = write_npz(tmp / f"blend-bad-{len(PASS)}.npz", **kw)
        try:
            mod._load_ctx_a_prior(str(p))
            check(f"{name} refused", False)
        except (ValueError, TypeError):
            check(f"{name} refused", True)


def test_blend_apply_refusals(mod, tmp: Path) -> None:
    print("blend_rows apply-time refusals")
    genes, pt, design, b_ctrl, panel = bframe()
    _, cov, mat = blend_artifact(tmp, "ref.npz")
    base = dict(targets=np.asarray(["T001", "T003"]), genes=np.asarray(cov), matrix=mat,
                beta=0.35, mode="blend_rows", context="A",
                norm_rule="fixed_norm_on_covered", sha256="x")

    def run(prior, plabels=pt):
        return mod.apply_ctx_a_blend_rows(design, prior, genes, panel, plabels)

    p = dict(base, targets=np.asarray(["T001", "NOPE"]))
    try:
        run(p)
        check("off-panel target refused", False)
    except ValueError:
        check("off-panel target refused", True)
    p = dict(base, genes=np.asarray(["NOT_A_GENE"] * 20), matrix=mat)
    try:
        run(p)
        check("gene miss refused", False)
    except ValueError:
        check("gene miss refused", True)
    p = dict(base, genes=np.asarray(["G000"] + [f"G{i:03d}" for i in range(10, 29)]),
             matrix=np.random.default_rng(0).normal(0, 1, (2, 20)))
    try:
        run(p)
        check("panel-coordinate genes refused", False)
    except ValueError:
        check("panel-coordinate genes refused", True)
    m0 = mat.copy()
    m0[1] = 0.0
    p = dict(base, matrix=m0)
    try:
        run(p)
        check("zero artifact row refused", False)
    except ValueError:
        check("zero artifact row refused", True)
    try:
        run(base, plabels=None)
        check("panel_targets=None refused", False)
    except ValueError:
        check("panel_targets=None refused", True)
    p = dict(base, context="B")
    try:
        mod.build_ctx_a_context_effects(design, p, b_ctrl, genes, panel, ("A", "C"),
                                        panel_targets=pt)
        check("wrong-context refused", False)
    except ValueError:
        check("wrong-context refused", True)


def test_old_path_with_panel_kw(mod, tmp: Path) -> None:
    print("old modes unaffected by the new kwarg")
    genes = np.array([f"G{i:03d}" for i in range(50)])
    rng = np.random.default_rng(366)
    design = rng.normal(0.0, 0.5, (5, 50))
    b_ctrl = np.abs(rng.normal(0.5, 0.3, 50)) + 0.01
    panel = np.array([0, 1, 2], dtype=np.int64)
    pt = [f"T{i:03d}" for i in range(5)]
    vec = rng.choice([-1.0, 1.0], size=47)
    p = write_npz(tmp / "oldmode.npz", genes=genes[3:], mode=np.asarray("sign_block"),
                  vector=vec, context=np.asarray("A"))
    prior = mod._load_ctx_a_prior(str(p))
    a, _ = mod.apply_ctx_a_prior(design, prior, b_ctrl, genes, panel, head_k=5, block_k=10)
    ctx, _ = mod.build_ctx_a_context_effects(design, prior, b_ctrl, genes, panel,
                                             ("A", "B", "C"), head_k=5, block_k=10,
                                             panel_targets=pt)
    check("sign_block identical with panel_targets kwarg", bool(np.array_equal(ctx["A"], a)))
    check("loader has no blend keys on old mode", "beta" not in prior and "targets" not in prior)


def test_config_expression() -> None:
    print("E8 config expression key sets (evaluated from the shipped line)")
    src = PATCHED.read_text()
    line = next(l.strip() for l in src.splitlines() if "ctx_a_prior_beta" in l)
    expr = line[line.index("{"): line.rindex("}") + 1]
    fake_old = {"mode": "sign_block", "sha256": "s", "context": "A", "alpha": 1.0, "tau": 0.0}
    fake_new = {"mode": "blend_rows", "sha256": "s", "context": "A", "beta": 0.35,
                "targets": np.asarray(["T001", "T003"])}
    for stub, want in ((fake_old, {"ctx_a_prior_mode", "ctx_a_prior_sha256", "ctx_a_prior_context",
                                   "ctx_a_prior_alpha", "ctx_a_prior_tau"}),
                       (fake_new, {"ctx_a_prior_mode", "ctx_a_prior_sha256", "ctx_a_prior_context",
                                   "ctx_a_prior_beta", "ctx_a_prior_targets"})):
        got = eval(expr, {"CTX_A_PRIOR": stub, "np": np})  # noqa: S307 — local test stub
        check(f"{stub['mode']} config keys == {sorted(want)}", set(got) == want)
    check("blend targets count==2", eval(expr, {"CTX_A_PRIOR": fake_new, "np": np})["ctx_a_prior_targets"] == 2)  # noqa: S307


def test_scripts_untouched373() -> None:
    print("scripts/ untouched by 373")
    txt = REAL.read_text()
    check("real builder has no CTX_A_PRIOR", "CTX_A_PRIOR" not in txt)
    check("real builder has no blend_rows", "blend_rows" not in txt)
    src = PATCHED.read_text()
    check("combined copy carries 366+373 anchors",
          all(s in src for s in ("_load_ctx_a_prior", "apply_ctx_a_prior", "apply_ctx_a_blend_rows",
                                 "build_ctx_a_context_effects", "ctx_a_prior_mode",
                                 "ctx_a_prior_beta", "panel_targets")))


def test_defaults(mod, tmp: Path) -> None:
    print("loader defaults (context/norm_rule absent)")
    cov = np.array([f"G{i:03d}" for i in range(10, 30)])
    mat = np.random.default_rng(9).normal(0, 1, (1, 20))
    p = write_npz(tmp / "defaults.npz", targets=np.asarray(["T002"]), genes=np.asarray(cov),
                  matrix=mat, beta=np.asarray(1.0), mode=np.asarray("blend_rows"))
    prior = mod._load_ctx_a_prior(str(p))
    check("default context A + norm_rule fixed", prior["context"] == "A"
          and prior["norm_rule"] == "fixed_norm_on_covered")


def main() -> None:
    assert PATCHED.is_file(), "run patch373.py first (work/patched copy missing)"
    # Part 1: 366's suite, unmodified, against the combined copy.
    print("=== part 1: 366's 45 tests on the combined copy ===")
    spec = importlib.util.spec_from_file_location("test366_rerun", T366)
    t366 = importlib.util.module_from_spec(spec)
    sys.modules["test366_rerun"] = t366
    spec.loader.exec_module(t366)
    t366.PATCHED = PATCHED
    t366.main()
    n366 = len(t366.PASS)
    # Part 2: blend_rows.
    print("=== part 2: blend_rows ===")
    with tempfile.TemporaryDirectory(prefix="t373-") as td:
        tmp = Path(td)
        mod = load_patched({}, "builder373_base")
        test_blend_exact(mod, tmp)
        test_beta_edges(mod, tmp)
        test_blend_load_refusals(mod, tmp)
        test_blend_apply_refusals(mod, tmp)
        test_old_path_with_panel_kw(mod, tmp)
        test_config_expression()
        test_scripts_untouched373()
        test_defaults(mod, tmp)
    print(f"\n{n366 + len(PASS)} checks passed ({n366} inherited from 366), 0 failed")


if __name__ == "__main__":
    main()
