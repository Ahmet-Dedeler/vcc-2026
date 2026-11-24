#!/usr/bin/env python3
"""Metadata-only feasibility audit for four IGVF TF Perturb-seq lineages.

This script intentionally downloads only small public metadata artifacts. Every
HTTP response is capped at 5 MB, and every output is written below this
experiment directory.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import re
import statistics
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
EVIDENCE = ROOT / "evidence"
RAW = EVIDENCE / "raw"
MAX_BYTES = 5_000_000

REPO = "IGVF/tf_perturb_seq"
COMMIT = "589f5ff606e7f9d09ad55e938b10727d206d5a60"
RAW_BASE = f"https://raw.githubusercontent.com/{REPO}/{COMMIT}"
API_BASE = f"https://api.github.com/repos/{REPO}"

PANEL_PATH = ROOT.parents[1] / "data/challenge/pert_counts.csv"
PRIOR_CARDIOMYOCYTE_PROVENANCE = (
    ROOT.parent / "081-igvf-cardiomyocyte" / "PROVENANCE.json"
)
PRIOR_CARDIOMYOCYTE_RESULTS = (
    ROOT.parent / "081-igvf-cardiomyocyte" / "results.json"
)

DATASETS = {
    "hues8_esc_current_qc": {
        "short": "Huangfu_ESC",
        "lineage": "HUES8 embryonic stem cell",
        "run": "qc_barcode_filter_none",
        "cells": 130_789,
        "large_uri": (
            "gs://igvf-pertub-seq-pipeline-data/"
            "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/"
            "2026_04_13/outs/qc_barcode_filter_none/inference_mudata.h5mu"
        ),
        "published_size_gb": 19.4,
    },
    "hues8_esc_scout_run": {
        "short": "Huangfu_ESC",
        "lineage": "HUES8 embryonic stem cell",
        "run": "sceptre_v1",
        "cells": 189_599,
        "large_uri": (
            "gs://igvf-pertub-seq-pipeline-data/"
            "Huangfu_HUES8-embryonic-stemcell-differentiation_TF-Perturb-seq/"
            "2026_04_13/outs/sceptre_v1/inference_mudata.h5mu"
        ),
        "published_size_gb": 17.52,
        "energy_distance_path": (
            "datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_"
            "TF-Perturb-seq/sceptre_v1/energy_distance/pval_edist_full.csv"
        ),
    },
    "hues8_endoderm": {
        "short": "Huangfu_DE",
        "lineage": "HUES8 definitive endoderm",
        "run": "muddy_penguin",
        "cells": 269_491,
        "large_uri": (
            "gs://igvf-pertub-seq-pipeline-data/"
            "Huangfu_HUES8-definitive-endoderm-differentiation_TF-Perturb-seq/"
            "2026_04_09/outs/muddy_penguin/inference_mudata.h5mu"
        ),
        "published_size_gb": 18.02,
        "energy_distance_path": (
            "datasets/Huangfu_HUES8-definitive-endoderm-differentiation_"
            "TF-Perturb-seq/muddy_penguin/energy_distance/pval_edist_full.csv"
        ),
    },
    "wtc11_hepatocyte": {
        "short": "Gersbach_Hep",
        "lineage": "WTC11-derived hepatocyte",
        "run": "sara_synapse_syn74842722",
        "cells": 1_053_820,
        "large_uri": "https://www.synapse.org/Synapse:syn74728027",
        "published_size_gb": 30.31,
    },
    "wtc11_cardiomyocyte": {
        "short": "Hon_CM",
        "lineage": "WTC11-derived cardiomyocyte",
        "run": "2026_04_19_no_spacer",
        "cells": 1_004_091,
        "large_uri": "https://www.synapse.org/Synapse:syn74522725",
        "large_folder_uri": "https://www.synapse.org/Synapse:syn74520421",
        "published_size_gb": 16.65,
        "energy_distance_path": (
            "datasets/Hon_WTC11-cardiomyocyte-differentiation_"
            "TF-Perturb-seq/2026_04_19_no_spacer/energy_distance/"
            "pval_edist_full.csv"
        ),
        "annotation_table_path": (
            "datasets/Hon_WTC11-cardiomyocyte-differentiation_"
            "TF-Perturb-seq/2026_04_19_no_spacer/energy_distance/"
            "annotation_table.csv"
        ),
    },
}

PER_DATASET_BASE = (
    "docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/"
    "pipeline_qc/results/per_dataset_distributions"
)

SMALL_SOURCE_PATHS = {
    "production_h5mu.tsv": "CRISPRi_lineage_atlas/manifests/production_h5mu.tsv",
    "pipeline_runs_catalog.tsv": (
        "docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/"
        "pipeline_qc/pipeline_runs_catalog.tsv"
    ),
    "pipeline_manifest.tsv": (
        "docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/"
        "pipeline_qc/manifest.tsv"
    ),
    "library_counts.tsv": (
        "docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/"
        "summary/library_counts.tsv"
    ),
    "tf_notebook.ipynb": (
        "docs/jamborees/2026_UTSW/working_groups/wg1_data_qc/"
        "TF_expr/TF_expr_and_DEGs.ipynb"
    ),
}

SYNAPSE_IDS = {
    "hues8_esc_folder": "syn74835010",
    "hues8_endoderm_folder": "syn74834952",
    "hues8_esc_cnmf": "syn74893846",
    "hues8_endoderm_cnmf": "syn74893844",
    "hepatocyte_project": "syn70518849",
    "hepatocyte_h5mu": "syn74728027",
    "hepatocyte_processed_folder": "syn74842722",
    "cardiomyocyte_h5mu": "syn74522725",
    "cardiomyocyte_run_folder": "syn74520421",
    "cardiomyocyte_params": "syn74918982",
}

IGVF_CARDIOMYOCYTE_ANALYSIS_SET = "IGVFDS6332VCTO"
IGVF_CARDIOMYOCYTE_FILES = {
    "filtered_global_de": "IGVFFI0830FXFI",
    "unfiltered_global_de": "IGVFFI5989UAVX",
    "local_de": "IGVFFI7298IERA",
    "full_h5ad": "IGVFFI7637STPX",
    "shared_guide_library": "IGVFFI8270UPKB",
}

HUES8_MEASUREMENT_SETS = {
    "esc": [
        "IGVFDS0746MYRH",
        "IGVFDS0956UUQN",
        "IGVFDS2908DIHX",
        "IGVFDS2940LYGK",
        "IGVFDS3348KRMQ",
        "IGVFDS5934EZKT",
        "IGVFDS6247FKDR",
        "IGVFDS8623ONYE",
    ],
    "endoderm": [
        "IGVFDS0788TUIL",
        "IGVFDS1313VGHL",
        "IGVFDS1403KZYV",
        "IGVFDS1437XFJK",
        "IGVFDS2129VHBD",
        "IGVFDS2520ZEYH",
        "IGVFDS3434HGEI",
        "IGVFDS4079VXPX",
    ],
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fetch(url: str, *, method: str = "GET", max_bytes: int = MAX_BYTES) -> dict[str, Any]:
    request = urllib.request.Request(
        url,
        method=method,
        headers={
            "User-Agent": "vcc-experiment-131-metadata-audit",
            "Accept": "application/json,text/plain,*/*",
        },
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            declared = response.headers.get("Content-Length")
            if declared is not None and int(declared) > max_bytes:
                raise RuntimeError(f"refusing {declared}-byte response from {url}")
            body = b"" if method == "HEAD" else response.read(max_bytes + 1)
            if len(body) > max_bytes:
                raise RuntimeError(f"response exceeded {max_bytes} bytes: {url}")
            return {
                "url": url,
                "status": response.status,
                "headers": dict(response.headers.items()),
                "body": body,
            }
    except urllib.error.HTTPError as error:
        body = error.read(min(max_bytes, 64_000))
        return {
            "url": url,
            "status": error.code,
            "headers": dict(error.headers.items()),
            "body": body,
        }


def save_source(name: str, path: str) -> dict[str, Any]:
    result = fetch(f"{RAW_BASE}/{path}")
    if result["status"] != 200:
        raise RuntimeError(f"source fetch failed: {path}: HTTP {result['status']}")
    destination = RAW / name
    destination.write_bytes(result["body"])
    return {
        "name": name,
        "repo_path": path,
        "raw_url": result["url"],
        "bytes": len(result["body"]),
        "sha256": sha256(result["body"]),
    }


def read_tsv(data: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(data.decode()), delimiter="\t"))


def load_panel() -> tuple[list[str], str]:
    data = PANEL_PATH.read_bytes()
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    panel = [row["target_gene"].strip() for row in rows]
    if len(panel) != 300 or len(set(panel)) != 300:
        raise RuntimeError("VCC panel metadata is not exactly 300 unique symbols")
    return panel, sha256(data)


def parse_guide_var(data: bytes) -> dict[str, Any]:
    rows = read_tsv(gzip.decompress(data))
    id_to_symbols: dict[str, set[str]] = defaultdict(set)
    id_to_types: dict[str, set[str]] = defaultdict(set)
    unique_guides: set[str] = set()
    type_guides: Counter[str] = Counter()
    for row in rows:
        guide_id = row.get("guide_id", "").strip()
        target_id = row.get("intended_target_name", "").strip()
        symbol = row.get("gene_name", "").strip()
        guide_type = row.get("type", "").strip()
        if guide_id:
            unique_guides.add(guide_id)
        if guide_type:
            type_guides[guide_type] += 1
        if target_id:
            if symbol:
                id_to_symbols[target_id].add(symbol)
            if guide_type:
                id_to_types[target_id].add(guide_type)
    return {
        "rows": rows,
        "id_to_symbols": id_to_symbols,
        "id_to_types": id_to_types,
        "guide_rows": len(rows),
        "unique_guide_ids": len(unique_guides),
        "guide_type_counts": dict(sorted(type_guides.items())),
    }


def summarize_cells_per_target(
    data: bytes,
    guide_metadata: dict[str, Any],
    panel: set[str],
) -> dict[str, Any]:
    rows = read_tsv(data)
    id_to_symbols = guide_metadata["id_to_symbols"]
    id_to_types = guide_metadata["id_to_types"]

    non_ntc = [row for row in rows if row["is_ntc"].lower() != "true"]
    targeting = [
        row for row in non_ntc if "targeting" in id_to_types.get(row["target"], set())
    ]
    targeting_symbols = set().union(
        *(id_to_symbols.get(row["target"], set()) for row in targeting)
    )
    overlap = sorted(targeting_symbols & panel)
    cells = [float(row["n_cells"]) for row in targeting]
    cells_by_symbol: dict[str, float] = defaultdict(float)
    for row in targeting:
        for symbol in id_to_symbols.get(row["target"], set()):
            cells_by_symbol[symbol] += float(row["n_cells"])
    panel_cells = [cells_by_symbol[symbol] for symbol in overlap]

    thresholds = {}
    for threshold in (1, 10, 20, 30, 50, 100):
        retained = [row for row in targeting if float(row["n_cells"]) >= threshold]
        symbols = set().union(
            *(id_to_symbols.get(row["target"], set()) for row in retained)
        )
        thresholds[str(threshold)] = {
            "target_ids": len(retained),
            "target_symbols": len(symbols),
            "vcc_overlap": len(symbols & panel),
        }

    return {
        "all_rows_including_ntc": len(rows),
        "non_ntc_target_ids": len(non_ntc),
        "targeting_target_ids": len(targeting),
        "targeting_gene_symbols": len(targeting_symbols),
        "vcc_overlap_count": len(overlap),
        "vcc_overlap_symbols": overlap,
        "targeting_cells_min": min(cells),
        "targeting_cells_median": statistics.median(cells),
        "targeting_cells_max": max(cells),
        "vcc_target_cells_min": min(panel_cells),
        "vcc_target_cells_median": statistics.median(panel_cells),
        "vcc_target_cells_max": max(panel_cells),
        "vcc_targets_below_20_cells": sorted(
            symbol for symbol in overlap if cells_by_symbol[symbol] < 20
        ),
        "vcc_targets_below_100_cells": sorted(
            symbol for symbol in overlap if cells_by_symbol[symbol] < 100
        ),
        "thresholds": thresholds,
    }


def summarize_energy_distance(
    data: bytes,
    guide_metadata: dict[str, Any],
    panel: set[str],
) -> dict[str, Any]:
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    id_to_symbols = guide_metadata["id_to_symbols"]
    types = Counter(row["type"] for row in rows)
    targeting_symbols: set[str] = set()
    targeting_target_ids: set[str] = set()
    for row in rows:
        if row["type"] != "targeting":
            continue
        ensembl_id = row[""].split("|", 1)[0]
        targeting_target_ids.add(row[""])
        targeting_symbols.update(id_to_symbols.get(ensembl_id, set()))
    overlap = sorted(targeting_symbols & panel)
    return {
        "target_units_total": len(rows),
        "type_counts": dict(sorted(types.items())),
        "targeting_locus_units": len(targeting_target_ids),
        "targeting_gene_symbols": len(targeting_symbols),
        "vcc_overlap_count": len(overlap),
        "vcc_overlap_symbols": overlap,
        "note": (
            "Energy-distance units preserve target locus coordinates; counts are "
            "therefore larger than Ensembl-collapsed cells_per_target counts."
        ),
    }


def extract_notebook_evidence(data: bytes) -> dict[str, Any]:
    notebook = json.loads(data)
    gene_shapes: dict[str, list[int]] = {}
    calibrated_shapes: dict[str, list[int]] = {}

    def text_value(value: Any) -> str:
        if isinstance(value, list):
            return "".join(str(item) for item in value)
        return str(value or "")

    for cell in notebook["cells"]:
        source = "".join(cell.get("source", [])).strip()
        for prefix in ("HP", "DE", "ES"):
            if source == f"{prefix}_gene_var.shape":
                text = "".join(
                    text_value(output.get("data", {}).get("text/plain", ""))
                    for output in cell.get("outputs", [])
                )
                match = re.search(r"\((\d+),\s*(\d+)\)", text)
                if match:
                    gene_shapes[prefix] = [int(match.group(1)), int(match.group(2))]
        output_text = "\n".join(
            text_value(output.get("text", ""))
            for output in cell.get("outputs", [])
            if output.get("output_type") == "stream"
        )
        for prefix in ("HP", "DE", "ES"):
            match = re.search(
                rf"=== {prefix} ===.*?table shape: \((\d+),\s*(\d+)\)",
                output_text,
                flags=re.DOTALL,
            )
            if match:
                calibrated_shapes[prefix] = [int(match.group(1)), int(match.group(2))]
    return {
        "gene_var_shapes": gene_shapes,
        "calibrated_no_ribo_mt_table_shapes": calibrated_shapes,
        "interpretation": {
            "HP": "WTC11-derived hepatocyte",
            "DE": "HUES8 definitive endoderm",
            "ES": "HUES8 embryonic stem cell sceptre_v1",
        },
    }


def github_metadata(path: str) -> dict[str, Any]:
    response = fetch(f"{API_BASE}/contents/{path}?ref={COMMIT}")
    if response["status"] != 200:
        return {"path": path, "status": response["status"]}
    payload = json.loads(response["body"])
    return {
        "path": path,
        "status": 200,
        "bytes": payload.get("size"),
        "git_blob_sha": payload.get("sha"),
        "html_url": payload.get("html_url"),
        "download_url": payload.get("download_url"),
    }


def access_probe(
    url: str,
    *,
    method: str = "GET",
    max_bytes: int = 64_000,
) -> dict[str, Any]:
    response = fetch(url, method=method, max_bytes=max_bytes)
    body = response["body"].decode(errors="replace")
    parsed: Any
    try:
        parsed = json.loads(body) if body else None
    except json.JSONDecodeError:
        parsed = body[:2_000]
    selected_headers = {
        key.lower(): value
        for key, value in response["headers"].items()
        if key.lower()
        in {
            "content-length",
            "content-type",
            "etag",
            "last-modified",
            "x-goog-generation",
            "x-goog-hash",
            "x-goog-storage-class",
        }
    }
    return {
        "url": url,
        "method": method,
        "status": response["status"],
        "headers": selected_headers,
        "body": parsed,
    }


def main() -> None:
    EVIDENCE.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)

    source_manifest = []
    for name, path in SMALL_SOURCE_PATHS.items():
        source_manifest.append(save_source(name, path))

    per_dataset_data: dict[str, dict[str, bytes]] = {}
    for short in ("Huangfu_ESC", "Huangfu_DE", "Gersbach_Hep"):
        per_dataset_data[short] = {}
        for filename in ("cells_per_target.tsv", "guide_var.tsv.gz"):
            local_name = f"{short}__{filename}"
            path = f"{PER_DATASET_BASE}/{short}/{filename}"
            source_manifest.append(save_source(local_name, path))
            per_dataset_data[short][filename] = (RAW / local_name).read_bytes()

    for key in ("hues8_esc_scout_run", "hues8_endoderm"):
        path = DATASETS[key]["energy_distance_path"]
        name = f"{key}__pval_edist_full.csv"
        source_manifest.append(save_source(name, path))
        DATASETS[key]["energy_distance_local"] = name

    panel_list, panel_sha = load_panel()
    panel = set(panel_list)
    summaries: dict[str, Any] = {}
    guide_metadata_by_short: dict[str, dict[str, Any]] = {}

    for short, files in per_dataset_data.items():
        guide_metadata_by_short[short] = parse_guide_var(files["guide_var.tsv.gz"])

    for key, dataset in DATASETS.items():
        short = dataset["short"]
        current_summary = summarize_cells_per_target(
            per_dataset_data[short]["cells_per_target.tsv"],
            guide_metadata_by_short[short],
            panel,
        )
        summaries[key] = {
            **dataset,
            "current_post_qc_target_summary": current_summary,
        }
        if "energy_distance_local" in dataset:
            summaries[key]["run_specific_energy_distance_summary"] = (
                summarize_energy_distance(
                    (RAW / dataset["energy_distance_local"]).read_bytes(),
                    guide_metadata_by_short[short],
                    panel,
                )
            )

    notebook_evidence = extract_notebook_evidence(
        (RAW / "tf_notebook.ipynb").read_bytes()
    )
    response_axis = {
        "hues8_esc_scout_run": notebook_evidence["gene_var_shapes"]["ES"][0],
        "hues8_endoderm": notebook_evidence["gene_var_shapes"]["DE"][0],
        "wtc11_hepatocyte": notebook_evidence["gene_var_shapes"]["HP"][0],
    }
    for key, count in response_axis.items():
        summaries[key]["response_axis_genes"] = count

    public_artifact_paths = [
        (
            "datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_"
            "TF-Perturb-seq/sceptre_v1/energy_distance/pval_edist_full.csv"
        ),
        (
            "datasets/Huangfu_HUES8-embryonic-stemcell-differentiation_"
            "TF-Perturb-seq/sceptre_v1/energy_distance/target_by_target_matrix.csv"
        ),
        (
            "datasets/Huangfu_HUES8-definitive-endoderm-differentiation_"
            "TF-Perturb-seq/muddy_penguin/energy_distance/pval_edist_full.csv"
        ),
        (
            "datasets/Huangfu_HUES8-definitive-endoderm-differentiation_"
            "TF-Perturb-seq/muddy_penguin/energy_distance/target_by_target_matrix.csv"
        ),
        (
            "datasets/Gersbach_WTC11-hepatocyte-differentiation_TF-Perturb-seq/"
            "sara_synapse_syn74842722/energy_distance/config1_2.json"
        ),
    ]
    compact_artifacts = [github_metadata(path) for path in public_artifact_paths]

    access_probes: dict[str, Any] = {
        "synapse": {},
        "gcs": {},
        "igvf": {},
    }
    for label, synapse_id in SYNAPSE_IDS.items():
        access_probes["synapse"][label] = {
            "entity": access_probe(
                f"https://repo-prod.prod.sagebase.org/repo/v1/entity/{synapse_id}"
            ),
            "permissions": access_probe(
                f"https://repo-prod.prod.sagebase.org/repo/v1/entity/"
                f"{synapse_id}/permissions"
            ),
        }

    for key in ("hues8_esc_scout_run", "hues8_esc_current_qc", "hues8_endoderm"):
        gs_uri = DATASETS[key]["large_uri"]
        bucket_and_path = gs_uri.removeprefix("gs://")
        bucket, object_path = bucket_and_path.split("/", 1)
        url = f"https://storage.googleapis.com/{bucket}/{object_path}"
        access_probes["gcs"][key] = access_probe(url, method="HEAD")

    for lineage, accessions in HUES8_MEASUREMENT_SETS.items():
        access_probes["igvf"][lineage] = {}
        for accession in accessions:
            access_probes["igvf"][lineage][accession] = access_probe(
                f"https://api.data.igvf.org/measurement-sets/"
                f"{accession}/?format=json"
            )
    access_probes["igvf"]["construct_library"] = access_probe(
        "https://api.data.igvf.org/construct-library-sets/"
        "IGVFDS3299AXST/?format=json",
        max_bytes=1_000_000,
    )

    overlap_rows = []
    for key, summary in summaries.items():
        current = summary["current_post_qc_target_summary"]
        for symbol in current["vcc_overlap_symbols"]:
            overlap_rows.append(
                {
                    "dataset_key": key,
                    "run": summary["run"],
                    "target_gene": symbol,
                    "basis": "current committed cells_per_target + guide_var",
                }
            )
    with (ROOT / "vcc_overlap.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["dataset_key", "run", "target_gene", "basis"],
        )
        writer.writeheader()
        writer.writerows(overlap_rows)

    result = {
        "audit_date": "2026-08-30",
        "scope": "metadata-only; no multi-gigabyte assets; no expression values",
        "upstream_repo": {
            "repo": REPO,
            "commit": COMMIT,
            "url": f"https://github.com/{REPO}/tree/{COMMIT}",
        },
        "vcc_panel": {
            "path": str(PANEL_PATH),
            "sha256": panel_sha,
            "count": len(panel_list),
        },
        "shared_library": {
            "construct_library_accession": "IGVFDS3299AXST",
            "guide_file_accession": "IGVFFI8270UPKB",
            "published_counts": {
                "guides": 14_151,
                "unique_intended_target_loci": 2_061,
                "unique_gene_names": 2_062,
                "targeting_loci": 1_950,
                "targeting_gene_names": 1_951,
            },
        },
        "datasets": summaries,
        "notebook_evidence": notebook_evidence,
        "compact_artifacts": compact_artifacts,
        "compact_effect_artifact_conclusion": {
            "public_gene_level_effect_artifact_found": False,
            "public_target_level_artifacts": (
                "HUES8 ESC and endoderm publish small energy-distance phenotype "
                "summaries, but these do not contain signed gene-level response vectors."
            ),
            "internal_gene_level_tables": (
                "The committed notebook references calibrated_all_results.tsv files "
                "for ES, DE, and HP and records broad table shapes, but those files "
                "are absent from the public Git tree and their containing Synapse "
                "entities are not anonymously readable."
            ),
        },
        "source_manifest": source_manifest,
    }

    (ROOT / "summary.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    (ROOT / "access_probes.json").write_text(
        json.dumps(access_probes, indent=2, sort_keys=True) + "\n"
    )
    (ROOT / "source_manifest.json").write_text(
        json.dumps(source_manifest, indent=2, sort_keys=True) + "\n"
    )

    print(json.dumps(
        {
            "datasets": {
                key: {
                    "run": value["run"],
                    "targeting_gene_symbols": value[
                        "current_post_qc_target_summary"
                    ]["targeting_gene_symbols"],
                    "vcc_overlap": value["current_post_qc_target_summary"][
                        "vcc_overlap_count"
                    ],
                    "response_axis_genes": value.get("response_axis_genes"),
                }
                for key, value in summaries.items()
            },
            "outputs": [
                "summary.json",
                "access_probes.json",
                "source_manifest.json",
                "vcc_overlap.csv",
            ],
        },
        indent=2,
    ))


if __name__ == "__main__":
    main()
