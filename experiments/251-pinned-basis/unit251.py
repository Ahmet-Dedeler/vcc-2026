"""251 unit level -- the flag is byte-inert OFF, and ON it is 250's `deflate_gene_basis` to the bit.

Two builders in one process (249's `unit249.py` pattern): the working tree's and the pre-251
byte copy (`work/make_submission_pm_cd4.PRE.py`, git HEAD 89cb22a), each loaded with `__file__`
pinned to `scripts/`.  Sections:

  OFF     `apply_frozen_cd4_residual` output + both audits bit-identical to pre-251, with and
          without an H1-style extra; no `deflation` key.
  ON      synthetic residual (complete / partial / untouched rows): the builder's
          `deflate_cd4_residual_pinned` with the residual's OWN gene basis == the re-fit
          `deflate_cd4_residual(k)` to 1e-12 (the two arithmetics of one projector), ==
          `lib250.deflate_gene_basis` bit for bit; audits; through `apply_frozen_cd4_residual`
          the non-carried rows move by ONE scalar and the panel coordinates stay zero even when
          the basis has mass on them; a foreign basis obeys Pythagoras on full-support rows;
          wrong shapes refused.
  LOADER  the real artifact loads and == 250's `basis250-production.npz[V5_abc]` bit for bit;
          every refusal fires: unregistered name, wrong file sha, tampered file (not
          orthonormal / V sha), wrong k, missing file, wrong gene axis.
  REAL    the real A/B/C residual (249's `cd4-off.npz`) on a one-hot stub base whose median row
          norm is EXACTLY the shipped 4.229295604454299: the re-fit path (K=5) re-solves the
          production fraction 0.34104860341941284 bit for bit (known answer); the pinned path
          (`V5_abc`) agrees to <= 1e-14 on the fraction and <= 1e-12 on the output array -- the
          unit-level preview of known-answer row 2.

    uv run python experiments/251-pinned-basis/unit251.py
"""
from __future__ import annotations

import importlib.util
import json
import os
import shutil
import sys
import tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
for key in ("VCC_CD4_DEFLATE_K", "VCC_CD4_DEFLATE_BASIS", "VCC_CD4_SIZING", "VCC_CD4_TARGET_G",
            "VCC_UNION_SOURCES", "VCC_UNION2_SOURCES", "VCC_UNION2_G"):
    os.environ.pop(key, None)
os.environ["VCC_ALPHA"] = "0.982"
os.environ["VCC_CD4_TARGET_F"] = "0.34"
os.environ["VCC_H1_UNION_G"] = "0.28"
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(HERE))                       # lib250 (byte copy of 250's) for the reference


def load_builder(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / "make_submission_pm_cd4.py")
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def prepared(M, n_rows, n_genes, complete, partial, cols, rng):
    effect = np.zeros((n_rows, n_genes))
    votes = np.zeros((n_rows, n_genes), dtype=np.uint8)
    counts = np.zeros(n_rows, dtype=np.uint8)
    for row in complete:
        c = rng.choice(cols, size=len(cols) - 3, replace=False)   # ragged support per row
        effect[row, c] = rng.uniform(0.2, 1.5, size=c.size) * rng.choice([-1, 1], size=c.size)
        effect[row, c] += 0.6 * rng.standard_normal()             # a shared across-row component
        votes[row, c] = 1
        counts[row] = 3
    for row, n in partial:
        counts[row] = n                                            # zero effect, count kept
    return M.EXP062.CD4Prepared(mode="stub", mean_effect=effect, gene_votes=votes,
                                target_condition_count=counts, condition_details=[])


def gene_basis(res, k):
    """250 `lib250.gene_basis`, verbatim (V_k = R_c^T U_k / s_k)."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    sup = np.asarray(res.gene_votes) > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(np.asarray(res.target_condition_count) > 0)
    R = eff[rows]
    Rc = R - R.mean(axis=0, keepdims=True)
    Rc[~sup[rows]] = 0.0
    w, U = np.linalg.eigh(Rc @ Rc.T)
    idx = np.argsort(w)[::-1][:k]
    Uk, sk = U[:, idx], np.sqrt(np.maximum(w[idx], 1e-300))
    return (Rc.T @ Uk) / sk[None, :], float(w[idx].sum() / max(w.sum(), 1e-300))


def stub_info(name="stub", k=3):
    return {"mode": "pinned", "basis": name, "basis_path": f"<{name}>", "basis_sha256": "0" * 64,
            "basis_v_sha256": "1" * 64, "basis_u_k_sha256": "2" * 64, "k": k, "n_genes": None,
            "orthonormality_err": 0.0, "gene_axis_verified": False, "provenance": {"experiment": "unit251"}}


def main() -> None:
    import lib250 as L  # noqa: E402  (250's construction, byte copy; imports the 233 instrument)
    new = load_builder(ROOT / "scripts/make_submission_pm_cd4.py", "builder_251")
    old = load_builder(WORK / "make_submission_pm_cd4.PRE.py", "builder_pre251")
    assert new.CD4_DEFLATE_K is None and new.CD4_DEFLATE_BASIS is None
    assert not hasattr(old, "CD4_DEFLATE_BASIS") and old.CD4_DEFLATE_K is None
    checks: dict[str, bool] = {}
    notes: dict[str, object] = {}

    n_rows, n_genes = 40, 60
    panel_positions = np.array([0, 1, 2], dtype=np.int64)
    columns = np.arange(3, n_genes, dtype=np.int64)
    rng = np.random.default_rng(251)
    whitened = rng.standard_normal((n_rows, n_genes))
    whitened[:, panel_positions] = 0.0
    complete = list(range(0, 30))
    partial = [(30, 1), (31, 2), (32, 2)]
    cols = np.arange(3, n_genes)
    cd4_new = prepared(new, n_rows, n_genes, complete, partial, cols, np.random.default_rng(1))
    cd4_old = prepared(old, n_rows, n_genes, complete, partial, cols, np.random.default_rng(1))
    extra_new = prepared(new, n_rows, n_genes, list(range(25, 36)), [], cols, np.random.default_rng(2))
    extra_old = prepared(old, n_rows, n_genes, list(range(25, 36)), [], cols, np.random.default_rng(2))
    assert np.array_equal(cd4_new.mean_effect, cd4_old.mean_effect)

    # ---------------------------------------------------------------- OFF: bit-identical
    for label, xn, xo in (("no extra", None, None), ("with extra", extra_new, extra_old)):
        on, an, un = new.apply_frozen_cd4_residual(whitened, cd4_new, columns, panel_positions,
                                                   extra=xn, extra_audit={})
        oo, ao, uo = old.apply_frozen_cd4_residual(whitened, cd4_old, columns, panel_positions,
                                                   extra=xo, extra_audit={})
        checks[f"OFF ({label}): output array bit-identical to pre-251"] = np.array_equal(on, oo)
        checks[f"OFF ({label}): cd4 audit identical to pre-251"] = (
            json.dumps(an, sort_keys=True) == json.dumps(ao, sort_keys=True))
        checks[f"OFF ({label}): union audit identical to pre-251"] = (
            json.dumps(un, sort_keys=True) == json.dumps(uo, sort_keys=True))
        checks[f"OFF ({label}): no `deflation` key in the audit"] = "deflation" not in an
    out_off, audit_off, union_off = new.apply_frozen_cd4_residual(
        whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={})

    # ---------------------------------------------------------------- ON, synthetic: one projector, two arithmetics
    k = 3
    V_own, topk = gene_basis(cd4_new, k)
    d_refit, a_refit = new.deflate_cd4_residual(cd4_new, k)
    d_pin, a_pin = new.deflate_cd4_residual_pinned(cd4_new, V_own, basis_info=stub_info(k=k))
    d_250, a_250 = L.deflate_gene_basis(cd4_new, V_own, "own")
    dmax = float(np.abs(d_pin.mean_effect - d_refit.mean_effect).max())
    notes["synthetic pinned(own V) vs refit max_abs_diff"] = dmax
    checks["ON: pinned(own basis) == re-fit(k) to 1e-12 (same projector, different association)"] = dmax < 1e-12
    checks["ON: pinned energy_removed_frac == re-fit's to 1e-12"] = (
        abs(a_pin["energy_removed_frac"] - a_refit["energy_removed_frac"]) < 1e-12)
    checks["ON: captured_frac_of_centred == the re-fit's top_k share to 1e-12"] = (
        abs(a_pin["captured_frac_of_centred"] - a_refit["top_k_energy_frac_of_centred"]) < 1e-12
        and abs(a_pin["own_top_k_energy_frac_of_centred"] - topk) < 1e-12)
    checks["ON: builder pinned deflation == lib250.deflate_gene_basis BIT FOR BIT"] = (
        np.array_equal(d_pin.mean_effect, d_250.mean_effect)
        and a_pin["energy_removed_frac"] == a_250["energy_removed_frac"]
        and a_pin["captured_frac_of_centred"] == a_250["captured_frac_of_centred"]
        and a_pin["n_carried"] == a_250["n_carried"] == len(complete) + len(partial))
    checks["ON: mode, votes, counts, details untouched"] = (
        d_pin.mode == cd4_new.mode and np.array_equal(d_pin.gene_votes, cd4_new.gene_votes)
        and np.array_equal(d_pin.target_condition_count, cd4_new.target_condition_count))
    checks["ON: zero outside the support, partial rows still zero"] = (
        not np.any(d_pin.mean_effect[d_pin.gene_votes == 0]) and not np.any(d_pin.mean_effect[30:]))
    checks["ON: audit carries the mandated keys (mode pinned, basis name/path/shas, energies, the diagnostic)"] = (
        {"mode", "basis", "basis_path", "basis_sha256", "basis_v_sha256", "basis_u_k_sha256", "k",
         "n_carried", "energy_removed_frac", "captured_frac_of_centred",
         "own_top_k_energy_frac_of_centred", "centred_energy_frac_of_raw", "gene_axis_verified"}
        <= set(a_pin) and a_pin["mode"] == "pinned")
    # a foreign basis on full-support rows: Pythagoras (energy removed == captured x centred/raw)
    cd4_full = prepared(new, n_rows, n_genes, complete, [], cols, np.random.default_rng(7))
    cd4_full.gene_votes[:30, 3:] = 1                               # full support on the carried rows
    cd4_full.mean_effect[:30, 3:] += 0.05 * np.random.default_rng(8).standard_normal((30, n_genes - 3))
    other = prepared(new, n_rows, n_genes, list(range(0, 20)), [], cols, np.random.default_rng(3))
    V_for, _ = gene_basis(other, 4)
    d_for, a_for = new.deflate_cd4_residual_pinned(cd4_full, V_for, basis_info=stub_info("foreign", 4))
    pyth = a_for["captured_frac_of_centred"] * a_for["centred_energy_frac_of_raw"]
    expect_removed = 1.0 - a_for["centred_energy_frac_of_raw"] * (1.0 - a_for["captured_frac_of_centred"])
    notes["foreign basis: energy_removed vs 1 - (1 - captured) x centred/raw"] = (
        a_for["energy_removed_frac"], expect_removed)
    checks["ON: foreign basis on full-support rows obeys Pythagoras to 1e-12 (a projector, not a norm cut)"] = (
        abs(a_for["energy_removed_frac"] - expect_removed) < 1e-12 and 0.0 < pyth < 1.0)
    Rk = d_for.mean_effect[:30]
    checks["ON: foreign basis: output rows orthogonal to V (max |R_k V| < 1e-10)"] = (
        float(np.abs(Rk @ V_for).max()) < 1e-10)
    for bad, label in ((V_own[:10], "wrong n_genes"), (V_own[:, :0], "k = 0"), (V_own.ravel(), "1-D")):
        try:
            new.deflate_cd4_residual_pinned(cd4_new, bad, basis_info=stub_info())
            checks[f"ON: a {label} basis is refused"] = False
        except ValueError:
            checks[f"ON: a {label} basis is refused"] = True
    empty = prepared(new, n_rows, n_genes, [], [], cols, np.random.default_rng(1))
    try:
        new.deflate_cd4_residual_pinned(empty, V_own, basis_info=stub_info())
        checks["ON: a residual with no carried rows is refused"] = False
    except ValueError:
        checks["ON: a residual with no carried rows is refused"] = True

    # ---------------------------------------------------------------- ON through the production entry point
    # a basis with mass on a PANEL gene (index 0): the re-zeroing must keep panel coordinates at 0
    Q, _ = np.linalg.qr(np.random.default_rng(9).standard_normal((n_genes, 3)))
    assert abs(Q[0]).max() > 0
    new.CD4_DEFLATE_BASIS = "stub"
    original_loader = new.load_cd4_deflation_basis
    try:
        new.load_cd4_deflation_basis = lambda name, official_genes=None: (V_own, stub_info(k=k))
        out_on, audit_on, union_on = new.apply_frozen_cd4_residual(
            whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={})
        new.load_cd4_deflation_basis = lambda name, official_genes=None: (Q, stub_info("panelmass", 3))
        out_q, audit_q, _ = new.apply_frozen_cd4_residual(
            whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={})
    finally:
        new.load_cd4_deflation_basis = original_loader
        new.CD4_DEFLATE_BASIS = None
    checks["ON (entry point): audit has deflation.mode == 'pinned' and basis == the name"] = (
        audit_on.get("deflation", {}).get("mode") == "pinned" and audit_on["deflation"]["basis"] == "stub")
    checks["ON (entry point): the fraction re-solves on the deflated object (moves vs OFF)"] = (
        audit_on["fraction"] != audit_off["fraction"])
    carried = np.zeros(n_rows, bool)
    carried[[r for r in complete] + [r for r, _ in partial]] = True
    touched = carried | (np.asarray(extra_new.target_condition_count) > 0)
    s = audit_on["post_add_preserve"]["post_add_global_preserve_scale"] / \
        audit_off["post_add_preserve"]["post_add_global_preserve_scale"]
    d_nc = float(np.abs(out_on[~touched] - s * out_off[~touched]).max())
    notes["entry point: non-carried rows vs OFF x preserve ratio max_abs_diff"] = d_nc
    checks["ON (entry point): rows CD4 does not carry == OFF's rows x ONE scalar (<= 1e-12)"] = d_nc < 1e-12
    checks["ON (entry point): carried rows move in direction (not a rescale)"] = (
        float(np.abs(out_on[:30] - s * out_off[:30]).max()) > 1e-3)
    checks["ON (entry point): the union's frame is untouched by the CD4 flag"] = (
        json.dumps(union_on["residual_frame"], sort_keys=True)
        == json.dumps(union_off["residual_frame"], sort_keys=True))
    checks["ON (entry point): a basis with mass on a panel gene leaves the panel coordinates at whitened x scale"] = (
        np.array_equal(out_q[:, panel_positions],
                       whitened[:, panel_positions] * audit_q["post_add_preserve"]["post_add_global_preserve_scale"]))
    # the union2 / K paths untouched: K=5 through the same entry point still runs the re-fit branch
    new.CD4_DEFLATE_K = k
    try:
        out_k, audit_k, _ = new.apply_frozen_cd4_residual(
            whitened, cd4_new, columns, panel_positions, extra=extra_new, extra_audit={})
    finally:
        new.CD4_DEFLATE_K = None
    dk = float(np.abs(out_k - out_on).max())
    notes["entry point: K=3 re-fit vs pinned(own V) output max_abs_diff"] = dk
    checks["ON (entry point): pinned(own basis) == the K re-fit through the whole entry point to 1e-12"] = dk < 1e-12
    checks["ON (entry point): the re-fit audit is 249's (no `mode`, U_k sha present)"] = (
        "mode" not in audit_k["deflation"] and "basis_u_k_sha256" in audit_k["deflation"])

    # ---------------------------------------------------------------- LOADER: the real artifact, every refusal
    genes = tuple(ln.strip() for ln in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:] if ln.strip())
    V5, info = new.load_cd4_deflation_basis("V5_abc", official_genes=genes)
    ref = np.load(ROOT / "experiments/250-deflation-october/work/basis250-production.npz")["V5_abc"]
    checks["LOADER: V5_abc loads; V == 250's basis250-production.npz[V5_abc] bit for bit"] = np.array_equal(V5, ref)
    checks["LOADER: info carries the pins (file sha, V sha, U_k sha f06e0262..., k=5, 18,533 genes, axis verified)"] = (
        info["basis_sha256"] == new.CD4_DEFLATION_BASIS_REGISTRY["V5_abc"]["sha256"]
        and info["basis_v_sha256"] == new.CD4_DEFLATION_BASIS_REGISTRY["V5_abc"]["v_sha256"]
        and info["basis_u_k_sha256"].startswith("f06e0262") and info["k"] == 5
        and info["n_genes"] == 18533 and info["gene_axis_verified"] is True)
    checks["LOADER: a second call without genes returns the verified cache"] = (
        new.load_cd4_deflation_basis("V5_abc")[1]["gene_axis_verified"] is True)
    wrong_axis = genes[1:] + genes[:1]
    try:
        new.load_cd4_deflation_basis("V5_abc", official_genes=wrong_axis)
        checks["LOADER: a wrong gene ORDER is refused"] = False
    except ValueError as exc:
        checks["LOADER: a wrong gene ORDER is refused"] = "gene ORDER" in str(exc) or "official axis" in str(exc)
    try:
        new.load_cd4_deflation_basis("nope")
        checks["LOADER: an unregistered name is refused"] = False
    except ValueError:
        checks["LOADER: an unregistered name is refused"] = True
    reg = new.CD4_DEFLATION_BASIS_REGISTRY
    saved = dict(reg["V5_abc"])
    tmp = Path(tempfile.mkdtemp(prefix="unit251-"))
    try:
        # wrong file sha in the registry
        new._CD4_DEFLATION_BASIS_CACHE.clear()
        reg["V5_abc"] = {**saved, "sha256": "0" * 64}
        try:
            new.load_cd4_deflation_basis("V5_abc")
            checks["LOADER: a file sha that is not the pin is refused"] = False
        except ValueError as exc:
            checks["LOADER: a file sha that is not the pin is refused"] = "does not match" in str(exc)
        # a tampered COPY (V scaled by 2: not orthonormal, V sha moves) pinned by its own file sha
        src = ROOT / saved["path"]
        with np.load(src, allow_pickle=False) as z:
            payload = {key: z[key] for key in z.files}
        payload["V"] = payload["V"] * 2.0
        tampered = tmp / "tampered.npz"
        np.savez(tampered, **payload)
        reg["V5_abc"] = {**saved, "path": str(tampered), "sha256": new.file_sha256(tampered)}
        new._CD4_DEFLATION_BASIS_CACHE.clear()
        # ROOT / absolute path == absolute path, so the registry entry resolves to the temp file
        try:
            new.load_cd4_deflation_basis("V5_abc")
            checks["LOADER: a tampered V (not orthonormal) is refused even when the file sha is re-pinned"] = False
        except ValueError as exc:
            checks["LOADER: a tampered V (not orthonormal) is refused even when the file sha is re-pinned"] = (
                "orthonormal" in str(exc))
        # wrong k in the registry
        reg["V5_abc"] = {**saved, "k": "6"}
        new._CD4_DEFLATION_BASIS_CACHE.clear()
        try:
            new.load_cd4_deflation_basis("V5_abc")
            checks["LOADER: a k that is not the registry's is refused"] = False
        except ValueError as exc:
            checks["LOADER: a k that is not the registry's is refused"] = "k 5 != the registry's 6" in str(exc)
        # a V sha pin that is not the file's V
        reg["V5_abc"] = {**saved, "v_sha256": "1" * 64}
        new._CD4_DEFLATION_BASIS_CACHE.clear()
        try:
            new.load_cd4_deflation_basis("V5_abc")
            checks["LOADER: a V sha that is not the pin is refused"] = False
        except ValueError as exc:
            checks["LOADER: a V sha that is not the pin is refused"] = "V hashes to" in str(exc)
        # a U_k provenance sha that is not the pin
        reg["V5_abc"] = {**saved, "u_k_sha256": "2" * 64}
        new._CD4_DEFLATION_BASIS_CACHE.clear()
        try:
            new.load_cd4_deflation_basis("V5_abc")
            checks["LOADER: a provenance U_k sha that is not the pin is refused"] = False
        except ValueError as exc:
            checks["LOADER: a provenance U_k sha that is not the pin is refused"] = "provenance U_k sha" in str(exc)
        # missing file
        reg["V5_abc"] = {**saved, "path": str(tmp / "missing.npz")}
        new._CD4_DEFLATION_BASIS_CACHE.clear()
        try:
            new.load_cd4_deflation_basis("V5_abc")
            checks["LOADER: a missing file is refused"] = False
        except FileNotFoundError:
            checks["LOADER: a missing file is refused"] = True
    finally:
        reg["V5_abc"] = saved
        new._CD4_DEFLATION_BASIS_CACHE.clear()
        shutil.rmtree(tmp, ignore_errors=True)
    V5, info = new.load_cd4_deflation_basis("V5_abc", official_genes=genes)
    checks["LOADER: the registry entry is restored and loads again"] = np.array_equal(V5, ref)

    # ---------------------------------------------------------------- REAL: 249's A/B/C residual on a one-hot stub base
    z = np.load(ROOT / "experiments/249-cd4-deflate/work/stage/cd4-off.npz")
    res = new.EXP062.CD4Prepared(mode="gene_norm", mean_effect=np.asarray(z["mean_effect"], np.float64),
                                 gene_votes=np.asarray(z["gene_votes"], np.uint8),
                                 target_condition_count=np.asarray(z["target_condition_count"], np.uint8),
                                 condition_details=[])
    targets = tuple(str(t) for t in z["targets"])
    gene_index = {g: i for i, g in enumerate(genes)}
    pp = np.asarray([gene_index[t] for t in targets], dtype=np.int64)
    n_g = len(genes)
    # the operator columns: every gene that is not a panel target (a superset of the real 17,942;
    # only the medians over `wcols` enter, and the stub rows are one-hot on non-panel genes)
    wcols = np.asarray(sorted(set(range(n_g)) - set(pp.tolist())), dtype=np.int64)
    shipped_median = 4.229295604454299
    base = np.zeros((len(targets), n_g))
    hot = [c for c in wcols if c not in set(pp.tolist())][:len(targets)]
    for i, c in enumerate(hot):
        base[i, c] = shipped_median                                 # every row norm == the median, exactly
    assert float(np.median(np.linalg.norm(base, axis=1))) == shipped_median
    new.CD4_DEFLATE_K = 5
    try:
        out_refit, aud_refit, _ = new.apply_frozen_cd4_residual(base, res, wcols, pp)
    finally:
        new.CD4_DEFLATE_K = None
    new.CD4_DEFLATE_BASIS = "V5_abc"
    try:
        out_pin, aud_pin, _ = new.apply_frozen_cd4_residual(base, res, wcols, pp)
    finally:
        new.CD4_DEFLATE_BASIS = None
    fr_refit, fr_pin = aud_refit["fraction"], aud_pin["fraction"]
    d_fr = abs(fr_pin - fr_refit)
    d_out = float(np.abs(out_pin - out_refit).max())
    notes["REAL: re-fit fraction (shipped 0.34104860341941284)"] = fr_refit
    notes["REAL: pinned fraction"] = fr_pin
    notes["REAL: |d fraction|"] = d_fr
    notes["REAL: output max_abs_diff pinned vs re-fit"] = d_out
    notes["REAL: pinned audit"] = {k_: v for k_, v in aud_pin["deflation"].items() if k_ != "provenance"}
    notes["REAL: re-fit audit"] = aud_refit["deflation"]
    checks["REAL: the re-fit path (K=5) re-solves the SHIPPED fraction 0.34104860341941284 bit for bit (known answer)"] = (
        fr_refit == 0.34104860341941284)
    checks["REAL: the re-fit audit is the shipped cd4svd5 audit (291 rows, 0.19866046506100588, f06e0262...)"] = (
        aud_refit["deflation"]["n_carried"] == 291
        and aud_refit["deflation"]["energy_removed_frac"] == 0.19866046506100588
        and aud_refit["deflation"]["basis_u_k_sha256"].startswith("f06e0262"))
    checks["REAL: pinned(V5_abc) fraction == the re-fit's to 1e-14"] = d_fr < 1e-14
    checks["REAL: pinned(V5_abc) output == the re-fit's to 1e-12 (numeric identity; NOT expected bitwise)"] = d_out < 1e-12
    checks["REAL: pinned energy_removed_frac == shipped 0.19866046506100588 to 1e-12; captured == 0.1982359547102955 to 1e-12"] = (
        abs(aud_pin["deflation"]["energy_removed_frac"] - 0.19866046506100588) < 1e-12
        and abs(aud_pin["deflation"]["captured_frac_of_centred"] - 0.1982359547102955) < 1e-12)
    checks["REAL: pinned audit records the pins and the axis check (basis V5_abc, file sha ba7ade3f..., verified axis)"] = (
        aud_pin["deflation"]["basis"] == "V5_abc" and aud_pin["deflation"]["basis_sha256"].startswith("ba7ade3f")
        and aud_pin["deflation"]["gene_axis_verified"] is True and aud_pin["deflation"]["n_carried"] == 291)
    notes["REAL: bitwise equal"] = bool(np.array_equal(out_pin, out_refit))

    n_fail = sum(not v for v in checks.values())
    for name, ok in checks.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}")
    print("\nnotes:")
    for key, val in notes.items():
        print(f"  {key}: {val}")
    (WORK / "unit251.json").write_text(json.dumps({"checks": checks, "notes": notes,
                                                  "peak_rss_gib": new.peak_rss_gib()}, indent=1, default=str))
    print(f"\nUNIT {len(checks) - n_fail}/{len(checks)} {'PASS' if not n_fail else 'FAIL'}  "
          f"peak RSS {new.peak_rss_gib():.2f} GiB")
    raise SystemExit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
