"""One-command bundle preflight: check a control bundle against the §6.1 format contract.

Transition item 4 (experiment 199). The writer ASSERTS the gene axis and `cells_per_pert`
deep inside a 13-minute, ~19 GiB build; an assertion that fires there costs the build. This
reads the same facts off the bundle in seconds, before anything is started, and prints one
PASS/FAIL table.

Every check is against a value READ FROM THE BUNDLE, never a constant of A/B/C's. Two traps
this exists to catch:

* **Both challenge CSVs carry a header row**, even though `vcc prep --genes` documents
  "headerless". `pd.read_csv` with a header is the correct reader for both, and `vcc prep`
  skips the header itself (its `--expected-gene-dim 18533` default passes on the real file,
  which would be 18,534 symbols if it did not). A hand-rolled `readlines()` is off by one
  and silently shifts the whole gene axis.
* **Gene ORDER must be read from the file.** `.var` index equality is order-sensitive and
  the writer compares per context. A D/E/F bundle is allowed to ship a different axis.

    uv run python scripts/def_preflight.py                          # A/B/C, known-answer row
    uv run python scripts/def_preflight.py --bundle <new bundle> --out work/preflight.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "scripts/frozen_inputs.json"
DENSITY_CAP = 4_750_000_000        # AGENTS.md §2, enforced by `vcc prep --max-nnz`
MAX_COUNTS_PER_CELL = 1_000_000
# 554: the K-class emitter's library classes are the control quartiles; above q75/q25 ~3.5 its steering dose (CLEAN4
# block prior, MDE, uncalled-gene hook) collapses to 0 WITHOUT an error (553, H1 sweep), by ~7 the placement crashes.
# A/B/C 1.57-1.95, H1 1.58, R7b as built 7.24. DEF_DAY_ALLOW_WIDE_LIBS=1 is the explicit override.
LIB_SPREAD_MAX = 3.0
# 641: the depth row. Arc downsampled every released context to one median depth (A/B/C control medians 20,109 /
# 19,946 / 20,034; the full data's median is 29,482 UMI/cell, Arc's 2026-10-01 post, Fig. 2). The recipe's
# depth-tuned constants (the DE call-set size above all) were chosen on A/B/C at that depth. DEF_DAY_ALLOW_DEPTH=1 is
# the explicit override (experiments/641-depth-row).
LIB_MEDIAN_REF = 20030.0
LIB_MEDIAN_TOL = 0.15
# AGENTS.md §2's measured A/B/C numbers, as the known-answer row.
KNOWN_ABC = {
    "n_genes": 18533, "n_constructs": 300, "cells_per_pert": 400, "control_cells": 18400,
    "n_ntc_ids": 46, "submission_rows": 360000,
    "A": {"nnz": 109910478, "median_lib": 20109.0, "min_lib": 3275.0, "max_lib": 52420.0},
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def _decode(a):
    if a.dtype.kind in "SO":
        return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in a])
    return np.asarray(a)


def _strings(node):
    """Read a string column under any of anndata's encodings.

    The 2026 control files write `var/_index` and `obs/_index` as
    `nullable-string-array` groups (`values` + `mask`), not as plain datasets, and
    `obs/context` as a `categorical` group. A reader that assumes one of these shapes
    silently fails on the others.
    """
    if isinstance(node, h5py.Group):
        if "categories" in node:                       # categorical
            cats = _decode(node["categories"][:])
            return cats[np.clip(node["codes"][:], 0, len(cats) - 1)]
        if "values" in node:                           # nullable-string-array
            return _decode(node["values"][:])
        raise KeyError(f"unrecognised string encoding: {list(node.keys())}")
    return _decode(node[:])


def _cat(g, name):
    node = g[name]
    if isinstance(node, h5py.Group) and "categories" in node:
        return _decode(node["categories"][:]), node["codes"][:]
    return None, _strings(node)


def inspect_control(path: Path, ctx: str, genes: list[str], manifest: dict) -> dict:
    """Read one control file's shape, dtypes, gene axis and count statistics, streaming X."""
    rec: dict = {"file": path.name, "bytes": path.stat().st_size, "sha256": sha256(path)}
    checks: list[tuple[str, bool, str]] = []
    with h5py.File(path, "r") as h:
        X = h["X"]
        enc = X.attrs.get("encoding-type", b"")
        enc = enc.decode() if isinstance(enc, bytes) else str(enc)
        shape = tuple(int(v) for v in X.attrs["shape"])
        rec.update(encoding=enc, shape=list(shape),
                   data_dtype=str(X["data"].dtype), indices_dtype=str(X["indices"].dtype),
                   indptr_dtype=str(X["indptr"].dtype), nnz=int(X["data"].shape[0]))
        per = manifest["per_context"][ctx]
        checks.append(("X is csr_matrix", enc == "csr_matrix", enc))
        checks.append(("X.data is float32", str(X["data"].dtype) == "float32", str(X["data"].dtype)))
        checks.append((f"shape == ({per['control_cells']}, {manifest['n_genes']})",
                       shape == (per["control_cells"], manifest["n_genes"]), str(shape)))
        checks.append(("indptr length == n_cells + 1",
                       int(X["indptr"].shape[0]) == shape[0] + 1, str(int(X["indptr"].shape[0]))))

        # gene axis, order-sensitive, read from the file
        v = h["var"]
        key = v.attrs.get("_index", "_index")
        key = key.decode() if isinstance(key, bytes) else key
        vgenes = _strings(v[key])
        same = list(map(str, vgenes)) == genes
        rec["var_index_equals_gene_names_csv_in_order"] = bool(same)
        checks.append(("var index == gene_names.csv, in order", same,
                       "identical" if same else f"{sum(a != b for a, b in zip(map(str, vgenes), genes))} positions differ"))

        # obs
        obs = h["obs"]
        cols = [c for c in obs.keys() if not c.startswith("_")]
        rec["obs_columns"] = cols
        ccol, pcol = manifest["context_col"], manifest["pert_col"]
        checks.append((f"obs has '{ccol}' and '{pcol}'", ccol in cols and pcol in cols, str(cols)))
        if ccol in cols:
            cats, codes = _cat(obs, ccol)
            labels = sorted(set(map(str, cats if cats is not None else codes)))
            rec["context_labels"] = labels
            checks.append((f"obs['{ccol}'] == ['{ctx}']", labels == [ctx], str(labels)))
        if pcol in cols:
            cats, codes = _cat(obs, pcol)
            labels = sorted(set(map(str, cats if cats is not None else codes)))
            rec["target_gene_labels"] = labels
            checks.append((f"all control rows are '{manifest['control_label']}'",
                           labels == [manifest["control_label"]], str(labels)))
        if "ntc_id" in cols:
            cats, _ = _cat(obs, "ntc_id")
            n = int(len(cats)) if cats is not None else None
            rec["n_ntc_ids"] = n
            checks.append((f"n_ntc_ids == {per['n_ntc_ids']}", n == per["n_ntc_ids"], str(n)))
        elif per.get("n_ntc_ids") is not None:
            # 526 finding 2: a manifest that declares NTC guides but a control file with no ntc_id column used to pass,
            # and the builder then took the 096 random C_p proxy instead of the guide-balanced draw
            # (make_submission_pm_cd4.py ~5220): a silent dispersion change. Fail closed; DEF_DAY_ALLOW_NO_NTC_ID=1 is the
            # explicit override for a release that genuinely drops the column.
            allow = os.environ.get("DEF_DAY_ALLOW_NO_NTC_ID", "") == "1"
            rec["n_ntc_ids"] = None
            checks.append((f"obs has 'ntc_id' (manifest declares n_ntc_ids {per['n_ntc_ids']})", allow,
                           "absent" + (" -- allowed by DEF_DAY_ALLOW_NO_NTC_ID=1" if allow else "")))

        # count statistics, streamed: never materialise 110M nonzeros
        indptr = X["indptr"][:].astype(np.int64)
        data = X["data"]
        lib = np.zeros(shape[0], dtype=np.float64)
        integral, nonneg, finite = True, True, True
        step = 1 << 23
        for start in range(0, int(data.shape[0]), step):
            block = data[start:start + step].astype(np.float64)
            integral &= bool(np.all(block == np.floor(block)))
            nonneg &= bool(np.all(block >= 0))
            finite &= bool(np.all(np.isfinite(block)))
            lo = np.searchsorted(indptr, start, "right") - 1
            hi = np.searchsorted(indptr, start + len(block), "left")
            for row in range(max(lo, 0), min(hi + 1, shape[0])):
                a, b = indptr[row], indptr[row + 1]
                a2, b2 = max(a, start), min(b, start + len(block))
                if b2 > a2:
                    lib[row] += block[a2 - start:b2 - start].sum()
        rec.update(
            nnz_per_cell=round(float(rec["nnz"]) / shape[0], 1),
            median_library=float(np.median(lib)), min_library=float(lib.min()),
            max_library=float(lib.max()), mean_library=round(float(lib.mean()), 1),
        )
        # 554: the K-class emitter's depth cliff (LIB_SPREAD_MAX). Always recorded; a row only when it trips, so a passing
        # bundle's report text is today's.
        _q25, _q75 = (float(v) for v in np.quantile(lib, [0.25, 0.75]))
        _spread = _q75 / _q25 if _q25 > 0 else float("inf")
        rec.update(library_q25=_q25, library_q75=_q75, library_spread_q75_q25=round(_spread, 4),
                   library_spread_max=LIB_SPREAD_MAX)
        if _spread > LIB_SPREAD_MAX:
            _allow = os.environ.get("DEF_DAY_ALLOW_WIDE_LIBS", "") == "1"
            rec["library_spread_allowed"] = _allow
            checks.append((f"control library q75/q25 <= {LIB_SPREAD_MAX} (554: above ~3.5 the emitter's steering dose "
                           "collapses silently)", _allow,
                           f"{_spread:.3f} (q25 {_q25:.0f}, q75 {_q75:.0f})"
                           + (" -- ALLOWED by DEF_DAY_ALLOW_WIDE_LIBS=1" if _allow else "")))
        # 641: the depth row (LIB_MEDIAN_REF / LIB_MEDIAN_TOL). Always recorded; a row only when it trips, so a passing
        # bundle's report text is today's.
        _ratio = float(np.median(lib)) / LIB_MEDIAN_REF
        rec.update(library_median_ref=LIB_MEDIAN_REF, library_median_ratio=round(_ratio, 4),
                   library_median_tol=LIB_MEDIAN_TOL)
        if abs(_ratio - 1.0) > LIB_MEDIAN_TOL:
            _allow_d = os.environ.get("DEF_DAY_ALLOW_DEPTH", "") == "1"
            rec["library_median_allowed"] = _allow_d
            checks.append((f"control library median within {LIB_MEDIAN_TOL:.0%} of A/B/C's {LIB_MEDIAN_REF:,.0f} "
                           "(641: the recipe's depth-tuned constants were chosen at that depth)", _allow_d,
                           f"{float(np.median(lib)):,.0f} (x{_ratio:.3f})"
                           + (" -- ALLOWED by DEF_DAY_ALLOW_DEPTH=1" if _allow_d else "")))
        checks.append(("X values are integral", integral, "yes" if integral else "NON-INTEGER"))
        checks.append(("X values are non-negative and finite", nonneg and finite, "yes"))
        checks.append((f"max library <= {MAX_COUNTS_PER_CELL}", lib.max() <= MAX_COUNTS_PER_CELL,
                       f"{lib.max():.0f}"))
    rec["checks"] = [{"check": c, "pass": bool(p), "value": v} for c, p, v in checks]
    return rec


def preflight(bundle: Path) -> dict:
    t0 = time.time()
    out: dict = {"bundle": str(bundle), "checks": []}

    def check(name, ok, value=""):
        out["checks"].append({"check": name, "pass": bool(ok), "value": str(value)})
        return bool(ok)

    mpath = bundle / "manifest.json"
    check("manifest.json present", mpath.exists(), str(mpath))
    manifest = json.loads(mpath.read_text())
    contexts = list(manifest["contexts"])
    out["manifest"] = manifest
    out["contexts"] = contexts
    out["files"] = {"manifest.json": {"sha256": sha256(mpath), "bytes": mpath.stat().st_size}}

    for key in ("season", "partition", "panel_id", "contexts", "pert_col", "context_col",
                "control_label", "n_genes", "n_constructs", "per_context", "cells_per_pert"):
        check(f"manifest has '{key}'", key in manifest, manifest.get(key, "ABSENT"))
    check("manifest declares 3 contexts", len(contexts) == 3, contexts)
    # 531 finding 14: the challenge's format is 400 cells per perturbation (AGENTS §2); a manifest saying otherwise would pass
    # here and fail at `vcc prep` after the whole build. DEF_DAY_ALLOW_CELLS_PER_PERT=<n> accepts a declared change on purpose.
    _cpp_ok = str(manifest.get("cells_per_pert")) in ("400", os.environ.get("DEF_DAY_ALLOW_CELLS_PER_PERT", "400"))
    check("manifest cells_per_pert == 400 (or DEF_DAY_ALLOW_CELLS_PER_PERT)", _cpp_ok, manifest.get("cells_per_pert", "ABSENT"))
    check("per_context covers every context",
          set(manifest["per_context"]) == set(contexts), list(manifest["per_context"]))
    for ctx in contexts:
        per = manifest["per_context"][ctx]
        check(f"per_context[{ctx}].n_perturbations == n_constructs",
              per["n_perturbations"] == manifest["n_constructs"], per["n_perturbations"])

    # --- gene_names.csv --------------------------------------------------------------
    gpath = bundle / "gene_names.csv"
    raw = gpath.read_text().splitlines()
    gdf = pd.read_csv(gpath)
    genes = [str(x) for x in gdf[gdf.columns[0]]]
    out["files"]["gene_names.csv"] = {"sha256": sha256(gpath), "bytes": gpath.stat().st_size,
                                      "header": raw[0], "n_rows_including_header": len(raw)}
    check("gene_names.csv has a HEADER row (vcc prep documents 'headerless')",
          gdf.columns[0] == "gene_name", gdf.columns[0])
    check(f"gene count == manifest n_genes ({manifest['n_genes']})",
          len(genes) == manifest["n_genes"], len(genes))
    check("gene symbols unique", len(set(genes)) == len(genes),
          f"{len(genes) - len(set(genes))} duplicates")
    check("no whitespace-padded gene symbols", all(g == g.strip() for g in genes), "")
    out["gene_axis"] = {"n": len(genes), "first3": genes[:3], "last3": genes[-3:]}

    # --- pert_counts.csv -------------------------------------------------------------
    ppath = bundle / "pert_counts.csv"
    pdf = pd.read_csv(ppath)
    panel = [str(x) for x in pdf[manifest["pert_col"]]]
    out["files"]["pert_counts.csv"] = {"sha256": sha256(ppath), "bytes": ppath.stat().st_size,
                                       "columns": list(pdf.columns)}
    check(f"pert_counts.csv has a '{manifest['pert_col']}' column",
          manifest["pert_col"] in pdf.columns, list(pdf.columns))
    check(f"target count == manifest n_constructs ({manifest['n_constructs']})",
          len(panel) == manifest["n_constructs"], len(panel))
    check("targets unique", len(set(panel)) == len(panel), f"{len(panel) - len(set(panel))} duplicates")
    check("no control label in the panel", manifest["control_label"] not in set(panel), "")
    on_axis = len(set(panel) & set(genes))
    check("every target is on the gene axis", on_axis == len(set(panel)),
          f"{on_axis}/{len(set(panel))}")
    if "n_cells" in pdf.columns:
        check("n_cells column is constant == cells_per_pert",
              set(pdf["n_cells"]) == {manifest["cells_per_pert"]}, sorted(set(pdf["n_cells"])))
    out["panel"] = {"n": len(panel), "first3": panel[:3], "last3": panel[-3:],
                    "has_n_cells_column": "n_cells" in pdf.columns}

    # --- the submission the writer will have to produce ------------------------------
    rows = manifest["n_constructs"] * manifest["cells_per_pert"] * len(contexts)
    out["submission_contract"] = {
        "rows": rows, "genes": manifest["n_genes"],
        "dense_entries": rows * manifest["n_genes"],
        "density_cap": DENSITY_CAP,
        "dense_over_cap_factor": round(rows * manifest["n_genes"] / DENSITY_CAP, 3),
        "max_nnz_per_cell_at_cap": round(DENSITY_CAP / rows, 1),
        "max_counts_per_cell": MAX_COUNTS_PER_CELL,
    }
    check("dense is over the density cap, so CSR is mandatory",
          rows * manifest["n_genes"] > DENSITY_CAP,
          f"{rows * manifest['n_genes'] / DENSITY_CAP:.2f}x")

    # --- control files ---------------------------------------------------------------
    out["control_files"] = {}
    for ctx in contexts:
        path = bundle / f"context_{ctx}.h5ad"
        if not check(f"context_{ctx}.h5ad present", path.exists(), str(path)):
            continue
        rec = inspect_control(path, ctx, genes, manifest)
        out["control_files"][ctx] = rec
        out["files"][path.name] = {"sha256": rec["sha256"], "bytes": rec["bytes"]}
        for c in rec["checks"]:
            out["checks"].append({"check": f"[{ctx}] {c['check']}", "pass": c["pass"],
                                  "value": c["value"]})

    # --- what the frozen table says about this bundle --------------------------------
    frozen = json.loads(FROZEN.read_text())
    fmatch = {k: frozen["files"].get(f"bundle/{k}") == v["sha256"] for k, v in out["files"].items()}
    out["frozen_inputs"] = {
        "path": str(FROZEN.relative_to(ROOT)), "sha256": sha256(FROZEN),
        "pinned_panel_id": frozen.get("panel_id"), "pinned_contexts": frozen.get("contexts"),
        "pinned_n_constructs": frozen.get("n_constructs"),
        "pinned_cells_per_pert": frozen.get("cells_per_pert"),
        "pinned_panel_support": frozen.get("panel_support"),
        "file_hash_matches": fmatch,
        "needs_repin_frozen_inputs": not all(fmatch.values()) or frozen.get("contexts") != contexts,
        "needs_repin_panel_support": (frozen.get("panel_id") != manifest.get("panel_id")
                                      or frozen.get("panel_support") is None),
    }
    out["seconds"] = round(time.time() - t0, 2)
    return out


def report(out: dict, verify: bool) -> bool:
    m = out["manifest"]
    print(f"bundle {out['bundle']}")
    print(f"  panel_id {m.get('panel_id')}  season {m.get('season')}  partition {m.get('partition')}  "
          f"contexts {out['contexts']}")
    print(f"  n_genes {m['n_genes']}  n_constructs {m['n_constructs']}  "
          f"cells_per_pert {m['cells_per_pert']}  control_label {m['control_label']!r}")
    sc = out["submission_contract"]
    print(f"  submission: {sc['rows']:,} rows x {sc['genes']:,} genes; dense is "
          f"{sc['dense_over_cap_factor']}x the {sc['density_cap']:,} cap; "
          f"max {sc['max_nnz_per_cell_at_cap']:,} stored entries/cell")
    print(f"  gene axis: {out['gene_axis']['first3']} ... {out['gene_axis']['last3']}")
    print(f"  panel:     {out['panel']['first3']} ... {out['panel']['last3']}")
    print("\ncontrol files")
    for ctx, r in out["control_files"].items():
        print(f"  {ctx}  {r['shape'][0]:,} x {r['shape'][1]:,}  {r['encoding']}/{r['data_dtype']}  "
              f"nnz {r['nnz']:,} ({r['nnz_per_cell']:,}/cell)  library med {r['median_library']:,.0f} "
              f"[{r['min_library']:,.0f}, {r['max_library']:,.0f}]  ntc {r.get('n_ntc_ids')}")
        print(f"     sha256 {r['sha256']}")
    fi = out["frozen_inputs"]
    print(f"\nfrozen_inputs.json ({fi['sha256'][:16]}...) pins panel {fi['pinned_panel_id']} "
          f"contexts {fi['pinned_contexts']} support {fi['pinned_panel_support']}")
    print(f"  --repin-frozen-inputs needed: {fi['needs_repin_frozen_inputs']}")
    print(f"  --repin-panel-support needed: {fi['needs_repin_panel_support']}")

    failed = [c for c in out["checks"] if not c["pass"]]
    print(f"\n{len(out['checks']) - len(failed)}/{len(out['checks'])} checks pass")
    for c in failed:
        print(f"  FAIL  {c['check']}: {c['value']}")
    ok = not failed

    if verify:
        print("\n--- known-answer row: AGENTS.md §2's measured A/B/C numbers ---")
        rowsA = out["control_files"].get("A", {})
        got = {"n_genes": m["n_genes"], "n_constructs": m["n_constructs"],
               "cells_per_pert": m["cells_per_pert"],
               "control_cells": rowsA.get("shape", [None])[0],
               "n_ntc_ids": rowsA.get("n_ntc_ids"),
               "submission_rows": out["submission_contract"]["rows"]}
        for k, expect in KNOWN_ABC.items():
            if k == "A":
                continue
            good = got[k] == expect
            ok &= good
            print(f"  {k:<16} {got[k]!s:<12} expected {expect}   {'OK' if good else 'MISMATCH'}")
        for k, expect in KNOWN_ABC["A"].items():
            gotv = rowsA.get({"nnz": "nnz", "median_lib": "median_library",
                              "min_lib": "min_library", "max_lib": "max_library"}[k])
            good = gotv == expect
            ok &= good
            print(f"  A.{k:<14} {gotv!s:<12} expected {expect}   {'OK' if good else 'MISMATCH'}")
        print(f"VERIFY {'PASS' if ok else 'FAIL'}")
    print(f"\n{out['seconds']}s")
    return ok


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bundle", type=Path, default=ROOT / "data/challenge")
    p.add_argument("--verify", action="store_true", help="assert AGENTS.md §2's A/B/C numbers")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()
    bundle = args.bundle if args.bundle.is_absolute() else ROOT / args.bundle
    out = preflight(bundle)
    ok = report(out, args.verify)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(out, indent=2) + "\n")
        print(f"wrote {args.out}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
