#!/usr/bin/env python3
"""D/E/F bundle shapes against preflight, the census gate, same-line, and the LODO axis assert.

No challenge matrix is read. Control files are tiny h5py CSR fixtures. Exits 0.
The classifications this locks:

  fail closed   bad CSV header, duplicate targets, context count != 3, obs label
                mismatch, var order != gene_names.csv, ntc_id present but the wrong
                category count, programme-source with no rows, LODO panel length != 300
                (lib438.bind), a gene axis that is not the CD4 V5 basis axis
  silent        ntc_id column absent (preflight skips it; the builder falls back to
                random_draw_cp), 1..5 targets with zero production donors (census
                passes; those rows are the panel-mean fallback), no registered
                same-line bank (artifact written with zero contexts)
  tautology     lib438 `assert official == OFFICIAL` when the bundle IS data/challenge:
                both sides are that file. A reordered axis of length 18533 is not
                what LODO refuses. The CD4 basis check is the one that raises.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[3]


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pre = load("def_preflight_526", ROOT / "scripts/def_preflight.py")
census = load("def_census_526", ROOT / "scripts/def_census.py")


def write_csr(path: Path, genes: list[str], ctx: str, ntc: list[str] | None, n_cells: int = 4) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    indptr = np.zeros(n_cells + 1, dtype=np.int32)
    indptr[1:] = 1
    with h5py.File(path, "w") as h:
        x = h.create_group("X")
        x.attrs["encoding-type"] = "csr_matrix"
        x.attrs["shape"] = (n_cells, len(genes))
        x.create_dataset("data", data=np.ones(n_cells, dtype=np.float32))
        x.create_dataset("indices", data=np.zeros(n_cells, dtype=np.int32))
        x.create_dataset("indptr", data=indptr)
        var = h.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=np.array(genes, dtype=h5py.string_dtype(encoding="utf-8")))
        obs = h.create_group("obs")

        def cat(group, name, values):
            g = group.create_group(name)
            cats = list(dict.fromkeys(values))
            codes = np.array([cats.index(v) for v in values], dtype=np.int32)
            g.create_dataset("categories", data=np.array(cats, dtype=h5py.string_dtype(encoding="utf-8")))
            g.create_dataset("codes", data=codes)

        cat(obs, "context", [ctx] * n_cells)
        cat(obs, "target_gene", ["non-targeting"] * n_cells)
        if ntc is not None:
            # repeat the category list across the few rows
            cat(obs, "ntc_id", [ntc[i % len(ntc)] for i in range(n_cells)])


def manifest(contexts, n_genes, ntc_ids, n_constructs=300, cells=400) -> dict:
    return {
        "season": "2026", "partition": "final", "panel_id": "synthetic",
        "contexts": contexts, "pert_col": "target_gene", "context_col": "context",
        "control_label": "non-targeting", "n_genes": n_genes, "n_constructs": n_constructs,
        "cells_per_pert": cells,
        "per_context": {c: {"n_perturbations": n_constructs, "control_cells": 4, "n_ntc_ids": ntc_ids} for c in contexts},
    }


def write_bundle(root: Path, *, contexts=("D", "E", "F"), swap_var=False, header=True,
                 dup_target=False, ntc: list[str] | None = None, ntc_ids=2,
                 pert_cols=("target_gene", "n_cells")) -> None:
    genes = [f"G{i:05d}" for i in range(18533)]
    var_genes = list(reversed(genes)) if swap_var else genes
    panel = genes[:300]
    if dup_target:
        panel = [genes[0]] * 300
    root.mkdir(parents=True, exist_ok=True)
    if header:
        (root / "gene_names.csv").write_text("gene_name\n" + "\n".join(genes) + "\n")
    else:
        (root / "gene_names.csv").write_text("\n".join(genes) + "\n")
    lines = [",".join(pert_cols)]
    for i, t in enumerate(panel):
        if pert_cols[0] == "target_gene":
            lines.append(f"{t},400")
        else:
            lines.append(f"400,{t}")
    (root / "pert_counts.csv").write_text("\n".join(lines) + "\n")
    (root / "manifest.json").write_text(json.dumps(manifest(list(contexts), 18533, ntc_ids)))
    for c in contexts:
        write_csr(root / f"context_{c}.h5ad", var_genes, c, ntc)


def fails(bundle: Path) -> list[str]:
    out = pre.preflight(bundle)
    bad = [c["check"] for c in out["checks"] if not c["pass"]]
    return bad


def main() -> int:
    import tempfile
    td = Path(tempfile.mkdtemp())
    ntc = ["N0", "N1"]

    good = td / "good"
    write_bundle(good, ntc=ntc)
    bad = fails(good)
    print(f"matching D/E/F fixture: {len(bad)} failed checks")
    assert bad == [], bad

    missing_ntc = td / "no-ntc"
    write_bundle(missing_ntc, ntc=None)
    bad = fails(missing_ntc)
    print(f"ntc_id column absent: {len(bad)} failed checks (preflight does not ask)")
    assert bad == []
    builder = (ROOT / "scripts/make_submission_pm_cd4.py").read_text()
    assert 'source = f"random_{cells_per_pert}_cell_draw' in builder

    bad = fails(td / "header") if False else None
    write_bundle(td / "header", header=False, ntc=ntc)
    bad = fails(td / "header")
    print("headerless gene_names:", bad[:2])
    assert any("HEADER" in b for b in bad)

    write_bundle(td / "dup", dup_target=True, ntc=ntc)
    bad = fails(td / "dup")
    print("duplicate targets:", [b for b in bad if "unique" in b or "target" in b][:3])
    assert any("unique" in b for b in bad)

    write_bundle(td / "two", contexts=("D", "E"), ntc=ntc)
    bad = fails(td / "two")
    print("two contexts:", [b for b in bad if "3 contexts" in b])
    assert any("3 contexts" in b for b in bad)

    write_bundle(td / "order", swap_var=True, ntc=ntc)
    bad = fails(td / "order")
    print("var order swapped:", [b for b in bad if "var index" in b])
    assert any("var index" in b for b in bad)

    write_bundle(td / "ntc-count", ntc=ntc, ntc_ids=46)
    bad = fails(td / "ntc-count")
    print("ntc categories 2, manifest 46:", [b for b in bad if "ntc" in b])
    assert any("n_ntc_ids" in b for b in bad)

    # obs label: rewrite one context after a good bundle
    write_bundle(td / "label", ntc=ntc)
    write_csr(td / "label" / "context_D.h5ad", [f"G{i:05d}" for i in range(18533)], "NOT_D", ntc)
    bad = fails(td / "label")
    print("obs context label:", [b for b in bad if "obs" in b][:2])
    assert any("NOT_D" in b or "context" in b for b in bad)

    # Census: 5 zero-donor targets pass, 6 fail. Real gate().
    def census_out(n_zero: int) -> dict:
        return {
            "donors": [{"donor": d, "panel_covered": 200} for d in census.PRODUCTION_DONORS],
            "production_donors": list(census.PRODUCTION_DONORS),
            "panel_support_predicted": {"cd4_banks_present": list(census.CD4_NAMES),
                                        "cd4_union_targets": 200, "cd4_complete_targets": 150},
            "consensus": {"targets_with_0_donors": n_zero},
        }
    ok5, _ = census.gate(census_out(5), census.GATE_DEFAULTS)
    ok6, rows6 = census.gate(census_out(6), census.GATE_DEFAULTS)
    print(f"census zero-donor targets: 5 -> {'PASS' if ok5 else 'FAIL'}, 6 -> {'PASS' if ok6 else 'FAIL'}")
    assert ok5 and not ok6
    assert census.GATE_DEFAULTS["max_zero_carrier"] == 5

    # Same-line: an unregistered call writes an artifact with no contexts. Mode off, not an exception.
    mini = td / "mini"
    mini.mkdir()
    (mini / "gene_names.csv").write_text("gene_name\nG0\nG1\n")
    (mini / "pert_counts.csv").write_text("target_gene\nG0\n")
    out_npz = td / "sameline.npz"
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts/def_sameline.py"),
         "--bundle", str(mini), "--force", "D=NOT_A_LINE,E=ALSO_NO,F=NOPE", "--out", str(out_npz)],
        cwd=ROOT, capture_output=True, text=True,
    )
    print("sameline unregistered rc", r.returncode)
    print(r.stdout.strip())
    assert r.returncode == 0, r.stderr
    nctx = subprocess.check_output(
        [sys.executable, str(ROOT / "scripts/def_day_helpers.py"), "nctx", str(out_npz)], text=True).strip()
    print("nctx", nctx)
    assert nctx == "0"
    # programme source that is not a bank dies
    r = subprocess.run(
        [sys.executable, str(ROOT / "scripts/def_sameline.py"),
         "--bundle", str(mini), "--force", "D=NOT_A_LINE", "--programme-source", "no-such-bank",
         "--out", str(td / "sameline2.npz")],
        cwd=ROOT, capture_output=True, text=True,
    )
    print("unknown programme-source rc", r.returncode, r.stderr.strip()[-200:])
    assert r.returncode != 0

    # In-mix banks: the same-line registry's primary banks overlap the production donor mix.
    # Nothing in def_sameline or def_day refuses that overlap. w is whatever the recipe passed.
    sam = (ROOT / "scripts/def_sameline.py").read_text()
    assert '"K562_HAEMATOPOIETIC_AND_LYMPHOID_TISSUE": ("k562_gwps"' in sam
    assert '"HCT116_LARGE_INTESTINE": ("orion_hct116"' in sam
    assert "PRODUCTION_DONORS" not in sam
    print("in-mix banks K562 and HCT116 are registry primaries; no production-mix check")

    # LODO tautology vs a different bundle path, and the CD4 axis raise.
    lib = (ROOT / "experiments/438-october-stack-repricing/lib438.py").read_text()
    assert "assert official == OFFICIAL" in lib
    plolo = (ROOT / "experiments/162-retrodicting-validator/plolo.py").read_text()
    assert 'split(",")[0]' in plolo
    assert 'read_csv_col(bundle / "pert_counts.csv", "target_gene")' in lib or \
        'read_csv_col(bundle / "pert_counts.csv", "target_gene")' in lib
    # column 0 vs the named column
    swapped = td / "swapped.csv"
    swapped.write_text("n_cells,target_gene\n400,TP53\n400,MYC\n")
    col0 = [ln.split(",")[0] for ln in swapped.read_text().splitlines()[1:]]
    named = [ln.split(",")[1] for ln in swapped.read_text().splitlines()[1:]]
    print("plolo column 0", col0, "named target_gene", named)
    assert col0 != named

    # same file read twice always compares equal, including after a permutation
    axis = td / "gene_names.csv"
    axis.write_text("gene_name\nB\nA\nC\n")

    def read_axis(p: Path) -> tuple:
        return tuple(ln.strip() for ln in p.read_text().splitlines()[1:] if ln.strip())

    official_pin = read_axis(axis)
    assert read_axis(axis) == official_pin          # bundle path is the pin's path
    other = td / "other_genes.csv"
    other.write_text("gene_name\nA\nB\nC\n")
    try:
        assert read_axis(other) == official_pin
        raised = False
    except AssertionError:
        raised = True
    print("LODO-style assert: same path passes; a second file with a different order raises:", raised)
    assert raised

    # The real CD4 check, exec'd from the builder.
    start = builder.index("def _verify_cd4_deflation_basis_axis")
    chunk = builder[start:builder.index("\ndef deflate_cd4_residual_pinned")]
    ns: dict = {}
    exec(chunk, ns)
    try:
        ns["_verify_cd4_deflation_basis_axis"]("V5_abc", ("A", "B", "C"), ("C", "B", "A"))
        cd4_raised = False
    except ValueError as exc:
        cd4_raised = True
        print("CD4 basis axis:", str(exc)[:140])
    assert cd4_raised

    # bind() refuses a panel that is not 300 long. The source assert is the October stop
    # (LODO is decided on). A length-300 panel with different symbols passes this assert.
    assert "assert len(panel) == 300 and len(set(panel)) == 300" in lib
    lodo = (ROOT / "scripts/def_lodo_key.py").read_text()
    assert "assert not (set(targets) & {L.OFFICIAL[i] for i in wcols})" in lodo
    print("Q1 CLASSIFICATION LOCKED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
