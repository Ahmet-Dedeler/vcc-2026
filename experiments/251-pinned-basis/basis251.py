"""251 -- export the PINNED CD4 deflation basis as a registry artifact (rule 11: a path and a sha).

Source: 249's staged A/B/C CD4 residual BEFORE deflation (`experiments/249-cd4-deflate/work/
stage/cd4-off.npz`: production complete-mean3, 300 rows, 291 carried, 276 complete) -- the object
whose k = 5 deflation audit the shipped `cd4svd5` meta records (`basis_u_k_sha256 f06e0262...`).

Construction (250 `lib250.gene_basis`, verbatim; 246/249's centring):
    rows = target_condition_count > 0 (291);  R = mean_effect[rows] masked to gene_votes > 0
    R_c  = R - mean over rows (per gene), re-zeroed outside the support
    G    = R_c R_c^T;  eigh;  U_k = top-k eigenvectors;  s_k = sqrt(eigenvalues)
    V_k  = R_c^T U_k / s_k          (genes x k, orthonormal columns)

KNOWN-ANSWER ROWS (abort on failure):
  KA-1  the BUILDER's `deflate_cd4_residual(res, 5)` on the source reproduces 249's shipped audit:
        n_carried 291, energy_removed_frac 0.19866046506100588,
        top_k_energy_frac_of_centred 0.1982359547102955, basis_u_k_sha256 f06e0262...
  KA-2  this script's U_k hashes to the same f06e0262...; its V_5 == 250's
        `work/basis250-production.npz[V5_abc]` BIT FOR BIT (the object 250's pinned rows used)
  KA-3  the gene axis is the official one (`data/challenge/gene_names.csv`, 18,533 genes,
        sha256 25bfa667...) and the source's target order is `data/challenge/pert_counts.csv`'s

Writes `data/donors/cd4_deflation_basis/V5_abc.npz` (+ a `.json` sidecar with every sha) and prints
the registry entry to paste into `scripts/make_submission_pm_cd4.py`.  Refuses to overwrite.

    uv run python experiments/251-pinned-basis/basis251.py [--force]
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
SOURCE = ROOT / "experiments/249-cd4-deflate/work/stage/cd4-off.npz"
REF250 = ROOT / "experiments/250-deflation-october/work/basis250-production.npz"
GENES_CSV = ROOT / "data/challenge/gene_names.csv"
PERTS_CSV = ROOT / "data/challenge/pert_counts.csv"
OUT_DIR = ROOT / "data/donors/cd4_deflation_basis"
OUT = OUT_DIR / "V5_abc.npz"
SIDECAR = OUT_DIR / "V5_abc.json"
K = 5
SHIPPED = dict(n_carried=291, energy_removed_frac=0.19866046506100588,
               top_k_energy_frac_of_centred=0.1982359547102955,
               basis_u_k_sha256="f06e0262cec2affb6c0594d2cebf8285888fd36ed73b518d514c7f0129b74e4c")
GENES_CSV_SHA = "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201"


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def array_sha256(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def load_builder(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    module.__file__ = str(ROOT / "scripts" / "make_submission_pm_cd4.py")
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true", help="overwrite an existing artifact")
    ap.add_argument("--builder", type=Path, default=ROOT / "scripts/make_submission_pm_cd4.py")
    args = ap.parse_args()
    for key in ("VCC_CD4_DEFLATE_K", "VCC_CD4_DEFLATE_BASIS"):
        os.environ.pop(key, None)
    MOD = load_builder(args.builder, "basis251_builder")      # sets the BLAS thread defaults (4)
    threads = {k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")}
    checks: dict[str, bool] = {}

    # ---------------------------------------------------------------- the source object
    src_sha = file_sha256(SOURCE)
    z = np.load(SOURCE)
    mean_effect = np.asarray(z["mean_effect"], dtype=np.float64)
    gene_votes = np.asarray(z["gene_votes"], dtype=np.uint8)
    counts = np.asarray(z["target_condition_count"], dtype=np.uint8)
    targets = tuple(str(t) for t in z["targets"])
    res = MOD.EXP062.CD4Prepared(mode="production_abc_off", mean_effect=mean_effect,
                                 gene_votes=gene_votes, target_condition_count=counts,
                                 condition_details=[])
    genes_lines = GENES_CSV.read_text().splitlines()
    assert genes_lines[0].strip() == "gene_name"
    genes = tuple(ln.strip() for ln in genes_lines[1:] if ln.strip())
    perts_lines = PERTS_CSV.read_text().splitlines()
    assert perts_lines[0].strip() == "target_gene"
    perts = tuple(ln.split(",")[0].strip() for ln in perts_lines[1:] if ln.strip())
    checks["KA-3 gene_names.csv sha256 == 25bfa667..."] = file_sha256(GENES_CSV) == GENES_CSV_SHA
    checks["KA-3 official axis has 18,533 genes == the source's columns"] = (
        len(genes) == 18533 == mean_effect.shape[1])
    checks["KA-3 source target order == data/challenge/pert_counts.csv"] = targets == perts
    print(f"[source] {SOURCE.relative_to(ROOT)} sha256 {src_sha}; mean_effect {mean_effect.shape} "
          f"sha256 {array_sha256(mean_effect)[:16]}...; carried {int((counts > 0).sum())}, "
          f"complete {int((counts == 3).sum())}", flush=True)

    # ---------------------------------------------------------------- KA-1: the builder's own deflation
    _d, audit = MOD.deflate_cd4_residual(res, K)
    for key, want in SHIPPED.items():
        got = audit[key]
        ok = (got == want) if isinstance(want, str) else (got == want if isinstance(want, int) else abs(got - want) < 1e-15)
        checks[f"KA-1 builder deflate_cd4_residual audit {key} == shipped {want}"] = bool(ok)
        print(f"[KA-1] {key}: builder {got!r} vs shipped {want!r} -> {'PASS' if ok else 'FAIL'}", flush=True)

    # ---------------------------------------------------------------- the basis (lib250.gene_basis, verbatim)
    eff = mean_effect.copy()
    sup = gene_votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(counts > 0)
    R = eff[rows]
    Rc = R - R.mean(axis=0, keepdims=True)
    Rc[~sup[rows]] = 0.0
    mean_rows = R.mean(axis=0)
    G = Rc @ Rc.T
    w, U = np.linalg.eigh(G)
    order = np.argsort(w)[::-1]
    idx = order[:K]
    Uk = np.ascontiguousarray(U[:, idx])
    sk = np.sqrt(np.maximum(w[idx], 1e-300))
    Vk = np.ascontiguousarray((Rc.T @ Uk) / sk[None, :])
    tot = float(max(w.sum(), 1e-300))
    u_sha = array_sha256(Uk)
    v_sha = array_sha256(Vk)
    gram_err = float(np.abs(Vk.T @ Vk - np.eye(K)).max())
    checks["KA-2 U_k sha256 == f06e0262... (249's shipped basis_u_k_sha256)"] = u_sha == SHIPPED["basis_u_k_sha256"]
    ref = np.load(REF250)["V5_abc"]
    checks["KA-2 V_5 == 250's basis250-production.npz[V5_abc] bit for bit"] = bool(np.array_equal(Vk, ref))
    checks["V_5 columns orthonormal (max |V^T V - I| < 1e-12)"] = gram_err < 1e-12
    # the two arithmetics on the source's own rows: same operator, different association
    Rk_refit = Rc - Uk @ (Uk.T @ Rc)
    Rk_refit[~sup[rows]] = 0.0
    P = Rc @ Vk
    Rk_pin = Rc - P @ Vk.T
    Rk_pin[~sup[rows]] = 0.0
    d_arith = float(np.abs(Rk_refit - Rk_pin).max())
    e0 = float((R * R).sum())
    ec = float((Rc * Rc).sum())
    captured = float((P * P).sum() / ec)
    removed_pin = float(1.0 - (Rk_pin * Rk_pin).sum() / e0)
    removed_refit = float(1.0 - (Rk_refit * Rk_refit).sum() / e0)
    print(f"[basis] U_k sha256 {u_sha}\n[basis] V_5 sha256 {v_sha}  shape {Vk.shape}  nonzero genes "
          f"{int((np.abs(Vk).sum(axis=1) > 0).sum())}  orthonormality err {gram_err:.2e}\n"
          f"[basis] top-{K} share of the centred energy {w[idx].sum() / tot!r}; captured by V_5 on its own rows "
          f"{captured!r}; energy removed pinned {removed_pin!r} vs re-fit {removed_refit!r}\n"
          f"[basis] pinned vs re-fit arithmetic on the SAME rows: max|d| {d_arith:.3e} on entries up to "
          f"{float(np.abs(Rk_refit).max()):.3f} (bitwise equal: {np.array_equal(Rk_refit, Rk_pin)})",
          flush=True)

    for name, ok in checks.items():
        print(f"  {'PASS' if ok else 'FAIL'}  {name}", flush=True)
    if not all(checks.values()):
        raise SystemExit("[basis251] KNOWN-ANSWER FAILURE -- nothing written")

    # ---------------------------------------------------------------- write
    if OUT.exists() and not args.force:
        raise SystemExit(f"{OUT} exists; pass --force to overwrite (and re-pin its sha everywhere)")
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    axis_sha = hashlib.sha256("\n".join(genes).encode("utf-8")).hexdigest()
    provenance = {
        "name": "V5_abc",
        "experiment": "251-pinned-basis",
        "date": dt.date.today().isoformat(),
        "k": K,
        "fitted_on": "A/B/C (data/challenge, panel vcc2026-val-1): 249's staged production CD4 residual, "
                     "300 rows, 291 carried (target_condition_count > 0), 276 complete",
        "source": {
            "path": str(SOURCE.relative_to(ROOT)),
            "sha256": src_sha,
            "mean_effect_sha256": array_sha256(mean_effect),
            "gene_votes_sha256": hashlib.sha256(np.ascontiguousarray(gene_votes).tobytes()).hexdigest(),
            "target_condition_count_sha256": hashlib.sha256(np.ascontiguousarray(counts).tobytes()).hexdigest(),
            "what": "the CD4 residual as scripts/make_submission_pm_cd4.py holds it BEFORE deflation "
                    "(prepare_frozen_cd4 -> freeze_complete_mean3), staged by 249 (stage249.py, builder "
                    "2d63b979...) at VCC_ALPHA=0.982 VCC_CD4_TARGET_F=0.34 VCC_H1_UNION_G=0.28",
            "n_rows": int(mean_effect.shape[0]), "n_carried": int(rows.size),
            "n_complete": int((counts == 3).sum()),
        },
        "construction": "rows = target_condition_count > 0; R = mean_effect[rows] masked to gene_votes > 0; "
                        "R_c = R - mean over rows (per gene), re-zeroed outside the support; G = R_c R_c^T; "
                        "np.linalg.eigh; U_k = top-k eigenvectors; s_k = sqrt(eigenvalues); V_k = R_c^T U_k / s_k "
                        "(experiments/250-deflation-october/lib250.py gene_basis, verbatim; U_k == "
                        "scripts/make_submission_pm_cd4.py deflate_cd4_residual's U_k, 249)",
        "application": "on ANY panel's CD4 residual (scripts/make_submission_pm_cd4.py "
                       "deflate_cd4_residual_pinned == lib250.deflate_gene_basis): centre over ITS carried rows "
                       "(per gene), R_k = R_c - (R_c V) V^T, re-zero outside ITS support; never re-fit (250 "
                       "sections 4.3-4.5)",
        "centring": "per-gene mean over the carried rows of the residual the basis is APPLIED to; the mean "
                    "of the A/B/C source rows is stored as `mean_abc` for provenance only and is never used",
        "basis_u_k_sha256": u_sha,
        "v_sha256": v_sha,
        "top_k_energy_frac_of_centred": float(w[idx].sum() / tot),
        "energy_removed_frac_on_source": removed_refit,
        "captured_frac_of_centred_on_source": captured,
        "eigenvalues_top15": [float(x) for x in w[order[:15]]],
        "spectrum_frac_top15": [float(x / tot) for x in w[order[:15]]],
        "orthonormality_err": gram_err,
        "pinned_vs_refit_max_abs_diff_on_source_rows": d_arith,
        "gene_axis": {"n": len(genes), "axis_sha256": axis_sha, "gene_names_csv": str(GENES_CSV.relative_to(ROOT)),
                      "gene_names_csv_sha256": GENES_CSV_SHA},
        "known_answers": {
            "builder_audit_equals_shipped_cd4svd5": True,
            "V5_equals_basis250_production_bitwise": True,
            "reference_250": str(REF250.relative_to(ROOT)),
        },
        "threads": threads,
        "numpy": np.__version__,
        "python": sys.version.split()[0],
    }
    np.savez(OUT, V=Vk, U=Uk, s=sk, eigenvalues_top15=np.asarray(w[order[:15]], dtype=np.float64),
             k=np.int64(K), genes=np.asarray(genes), carried_targets=np.asarray([targets[i] for i in rows]),
             mean_abc=np.asarray(mean_rows, dtype=np.float64), provenance=json.dumps(provenance))
    out_sha = file_sha256(OUT)
    with np.load(OUT, allow_pickle=False) as back:
        assert np.array_equal(back["V"], Vk) and int(back["k"]) == K
        assert tuple(str(g) for g in back["genes"]) == genes
        assert json.loads(str(back["provenance"]))["basis_u_k_sha256"] == u_sha
    sidecar = {"artifact": str(OUT.relative_to(ROOT)), "sha256": out_sha, "bytes": OUT.stat().st_size,
               "v_sha256": v_sha, "u_k_sha256": u_sha, "k": K, **{"provenance": provenance}}
    SIDECAR.write_text(json.dumps(sidecar, indent=2) + "\n")
    (WORK / "basis251.json").write_text(json.dumps({"checks": checks, **sidecar}, indent=2) + "\n")
    print(f"\n[artifact] {OUT.relative_to(ROOT)}  sha256 {out_sha}  {OUT.stat().st_size:,} B\n"
          f"[sidecar]  {SIDECAR.relative_to(ROOT)}\n\nregistry entry:\n"
          f'    "V5_abc": {{\n        "path": "{OUT.relative_to(ROOT)}",\n        "sha256": "{out_sha}",\n'
          f'        "v_sha256": "{v_sha}",\n        "u_k_sha256": "{u_sha}",\n        "k": "{K}",\n    }},',
          flush=True)


if __name__ == "__main__":
    main()
