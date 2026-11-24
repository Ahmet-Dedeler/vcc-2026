"""216 Phase 1 — is the published Cell atlas the same object our banks were built from?

Read-only.  Two network surfaces, both public and both carrying nothing about our method:
the Crossref record for `10.1016/j.cell.2026.08.002`, and the S3 bucket named in the
paper's own data-sharing readme.  Nothing is posted anywhere.

What this answers, concretely:

  1. is the object our banks were built from (`GWCD4i.DE_stats.h5ad`, sha `c355f535...`,
     16,786,240,107 B, pinned at `experiments/199-def-runbook/work/`) still the object the
     publication points at?  Byte-range spot checks against the live bucket, not a size
     comparison;
  2. what ELSE does the release provide, when was each piece published, and does anything
     at all land inside the 2026-08-28 window the two leaders' cosine doubled in;
  3. panel coverage: ours is 276 fully-supported / 291 union of 300.  What is the ceiling
     the published release could support, and what stops the missing nine;
  4. the gene axis.  Our banks reach 10,113 of the official 18,533 genes because the DE
     table has 10,282 columns.  The pseudobulk object has 18,129.  How much of the effect
     energy the scorer reads sits on genes our CD4 residual can never touch.

    uv run python atlas216.py
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
DL = ROOT / "data/external/cd4_cell2026"
BUCKET = "https://genome-scale-tcell-perturb-seq.s3.amazonaws.com"
SOURCE_KEY = "marson2025_data/GWCD4i.DE_stats.h5ad"
LOCAL_SOURCE = ROOT / "experiments/199-def-runbook/work/GWCD4i.DE_stats.h5ad"
PINNED_SHA = "c355f535ff32cf7ba1edc49cf9c6039fe84f2c9ebe4d005515cba75790cfbb62"
PINNED_BYTES = 16_786_240_107
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"

PANEL = tuple(line.split(",")[0]
              for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:]
              if line.strip())
OFFICIAL = tuple(line.strip()
                 for line in (ROOT / "data/challenge/gene_names.csv").read_text().splitlines()[1:]
                 if line.strip())


def http_range(key: str, start: int, length: int) -> bytes:
    url = f"{BUCKET}/{key}"
    out = subprocess.run(
        ["curl", "-sS", "-r", f"{start}-{start + length - 1}", url],
        capture_output=True, check=True)
    return out.stdout


def main() -> None:
    rep: dict[str, object] = {}

    # ------------------------------------------------------------------ 1. same object?
    assert LOCAL_SOURCE.exists(), f"pinned CD4 source missing at {LOCAL_SOURCE}"
    size = LOCAL_SOURCE.stat().st_size
    inv = json.loads((WORK / "s3_inventory.json").read_text())
    remote = {r["key"]: r for r in inv}
    rspec = remote[SOURCE_KEY]

    rng = np.random.default_rng(216)
    probes = [0, 4096, size // 3, size // 2, (2 * size) // 3, size - 65536]
    probes += [int(x) for x in rng.integers(0, size - 65536, 6)]
    checks = []
    with open(LOCAL_SOURCE, "rb") as fh:
        for off in sorted(set(int(p) for p in probes)):
            n = min(65536, size - off)
            fh.seek(off)
            local = fh.read(n)
            got = http_range(SOURCE_KEY, off, n)
            checks.append(dict(offset=off, n=n, match=bool(got == local),
                               local_md5=hashlib.md5(local).hexdigest()[:12],
                               remote_md5=hashlib.md5(got).hexdigest()[:12]))
            print(f"[range] {off:>14,}  {'MATCH' if checks[-1]['match'] else 'DIFFER'}",
                  flush=True)
    rep["source_identity"] = dict(
        local_path=str(LOCAL_SOURCE.relative_to(ROOT)), local_bytes=size,
        pinned_sha256=PINNED_SHA, pinned_bytes=PINNED_BYTES,
        remote_bytes=rspec["size"], remote_last_modified=rspec["mtime"],
        remote_etag=rspec["etag"],
        size_matches=bool(size == rspec["size"] == PINNED_BYTES),
        n_byte_ranges=len(checks), all_ranges_match=bool(all(c["match"] for c in checks)),
        ranges=checks)

    # ----------------------------------------------------- 2. what else, and when
    rep["release_inventory"] = [
        dict(key=r["key"], gib=round(r["size"] / 2**30, 2), last_modified=r["mtime"])
        for r in inv if r["size"] > 1e9]
    windows = {}
    for r in inv:
        windows.setdefault(r["mtime"][:10], []).append(r["key"])
    rep["objects_by_date"] = {d: len(v) for d, v in sorted(windows.items())}
    rep["objects_in_publication_window_0826_0903"] = sorted(
        k for d, v in windows.items() if "2026-08-26" <= d <= "2026-09-03" for k in v)

    # the one artifact that DID change in the window: the DE_stats supplementary table
    old = pd.read_csv(ROOT / "data/external/cd4t/DE_stats.suppl_table.csv")
    new = pd.read_csv(DL / "suppl_tables_2026-08-24/DE_stats.suppl_table.csv")
    rep["de_stats_suppl_table_change"] = dict(
        old_shape=list(old.shape), new_shape=list(new.shape),
        rows_identical=bool(len(old) == len(new)),
        new_columns=[c for c in new.columns if c not in old.columns],
        dropped_columns=[c for c in old.columns if c not in new.columns],
        same_perturbation_set=bool(
            set(zip(old.target_contrast, old.culture_condition))
            == set(zip(new.target_contrast, new.culture_condition))))

    # are those "new" columns actually new DATA, or were they already in the h5ad we hold?
    import h5py
    with h5py.File(LOCAL_SOURCE, "r") as h:
        obs_cols = sorted(h["obs"].keys())
        var_cols = sorted(h["var"].keys())
        layers = sorted(h["layers"].keys())
        varm = sorted(h["varm"].keys())
        n_obs, n_var = h["layers/log_fc"].shape
        gene_names = h["var/gene_name"]
        gene_names = (gene_names["categories"][:][gene_names["codes"][:]]
                      if isinstance(gene_names, h5py.Group) else gene_names[:])
        gene_names = np.asarray([g.decode() if isinstance(g, bytes) else str(g)
                                 for g in gene_names])
    rep["pinned_source_structure"] = dict(
        n_obs=int(n_obs), n_var=int(n_var), obs_columns=obs_cols, var_columns=var_cols,
        layers=layers, varm=varm)
    rep["suppl_new_columns_already_in_pinned_h5ad"] = sorted(
        c for c in rep["de_stats_suppl_table_change"]["new_columns"] if c in obs_cols)
    rep["suppl_new_columns_NOT_in_pinned_h5ad"] = sorted(
        c for c in rep["de_stats_suppl_table_change"]["new_columns"] if c not in obs_cols)

    # ----------------------------------------------------- 3. panel target coverage
    man = json.loads((CD4_ROOT / "manifest.json").read_text())
    ours_union = 300 - len(man["panel_union_missing"])
    trm = pd.read_csv(CD4_ROOT / "target_row_map.csv")
    per_cond = trm[trm.selected].groupby("condition").target_gene.nunique().to_dict()
    complete = (trm[trm.selected].groupby("target_gene").condition.nunique() == 3).sum()

    de = new  # the published per-perturbation table, 33,983 rows
    de_panel = de[de.target_contrast_gene_name.isin(PANEL)]
    lib = pd.read_csv(DL / "suppl_tables_2026-08-24/sgrna_library_metadata.suppl_table.csv")
    lib_panel = set(lib.target_gene_name.dropna()) & set(PANEL)
    lib_designed = set(lib.designed_target_gene_name.dropna()) & set(PANEL)
    missing = man["panel_union_missing"]
    rep["panel_coverage"] = dict(
        panel=300,
        ours_union=int(ours_union), ours_complete_all3=int(complete),
        ours_per_condition=per_cond,
        published_de_table_panel_targets=int(de_panel.target_contrast_gene_name.nunique()),
        published_de_table_panel_complete_all3=int(
            (de_panel.groupby("target_contrast_gene_name").culture_condition.nunique()
             == 3).sum()),
        guide_library_panel_targets_validated=len(lib_panel),
        guide_library_panel_targets_designed=len(lib_designed),
        never_covered=missing,
        never_covered_in_guide_library=sorted(g for g in missing if g in lib_designed),
        never_covered_in_de_table=sorted(
            g for g in missing if g in set(de.target_contrast_gene_name)))

    # why do the per-condition sets differ?  the published table's own eligibility flags
    flags = {}
    for g in missing:
        sub = de[de.target_contrast_gene_name == g]
        flags[g] = dict(rows=int(len(sub)),
                        conditions=sorted(sub.culture_condition.unique().tolist()))
    rep["panel_coverage"]["never_covered_rows_in_published_table"] = flags

    partial = sorted(set(PANEL) - set(trm[trm.selected].groupby("target_gene").filter(
        lambda d: d.condition.nunique() == 3).target_gene.unique()) - set(missing))
    rep["panel_coverage"]["covered_but_not_all_three_conditions"] = partial
    pub_flags = {}
    for g in partial:
        sub = de[de.target_contrast_gene_name == g]
        pub_flags[g] = sorted(sub.culture_condition.unique().tolist())
    rep["panel_coverage"]["partial_targets_conditions_in_published_table"] = pub_flags

    # ----------------------------------------------------- 4. the gene axis
    gm = pd.read_csv(CD4_ROOT / "gene_map.csv")
    ours_genes = set(gm.iloc[:, -1].dropna()) if "official" not in gm.columns else set()
    rep["gene_axis"] = dict(
        de_table_n_var=int(n_var),
        de_table_genes_on_official_axis=int(len(set(gene_names) & set(OFFICIAL))),
        our_bank_mapped_official_genes=int(man["mapped_official_genes"]),
        pseudobulk_n_var_per_readme=18129,
        official_axis=len(OFFICIAL),
        gene_map_rows=int(len(gm)), gene_map_columns=list(gm.columns))
    rep["gene_axis"]["_ours_genes_placeholder"] = len(ours_genes)

    (WORK / "atlas216.json").write_text(json.dumps(rep, indent=2, default=str))
    print(json.dumps({k: v for k, v in rep.items()
                      if k in ("panel_coverage", "gene_axis",
                               "objects_by_date",
                               "objects_in_publication_window_0826_0903",
                               "suppl_new_columns_already_in_pinned_h5ad",
                               "suppl_new_columns_NOT_in_pinned_h5ad")},
                     indent=2, default=str))


if __name__ == "__main__":
    main()
