#!/usr/bin/env python3
"""Metadata-only, fail-closed source audit for GSE337988.

This module intentionally has no large-file download path. It may download the
three small metadata assets declared in SOURCE_MANIFEST.json and may issue a
bounded HTTP range request that is sufficient to read the first CSV row of the
two candidate LFC matrices.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import tempfile
import urllib.request
import zlib
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MANIFEST_PATH = HERE / "SOURCE_MANIFEST.json"
REPORT_PATH = HERE / "METADATA_AUDIT.json"
METADATA_ROOT = HERE / "metadata"
TARGET_LABELS_PATH = METADATA_ROOT / "lfc_target_labels.txt"
LARGE_CANDIDATES = frozenset({"low1_lfc_matrix", "low2_lfc_matrix"})
USER_AGENT = "vcc-experiment-116-source-audit/1"


class SourceAuditError(RuntimeError):
    """A source contract, provenance, or authorization check failed."""


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def atomic_write(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
    )
    try:
        with os.fdopen(descriptor, "wb") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary_name, path)
    except BaseException:
        Path(temporary_name).unlink(missing_ok=True)
        raise


def load_manifest(path: Path = MANIFEST_PATH) -> dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise SourceAuditError("unsupported manifest schema")
    if manifest.get("experiment") != "116-dld1-lowmoi":
        raise SourceAuditError("manifest experiment changed")
    if manifest.get("dataset") != "GSE337988":
        raise SourceAuditError("manifest dataset must be GSE337988")
    if manifest.get("scope") != "source_only_metadata":
        raise SourceAuditError("manifest scope must remain source_only_metadata")

    authorization = manifest.get("authorization")
    if not isinstance(authorization, dict):
        raise SourceAuditError("authorization block is missing")
    if set(authorization.get("eligible_large_assets", [])) != LARGE_CANDIDATES:
        raise SourceAuditError(
            "large-asset eligibility must contain only Low1 and Low2"
        )
    if authorization.get("authorized_large_assets") != []:
        raise SourceAuditError("no large asset is currently authorized")
    for forbidden in (
        "download_geo_raw_archive",
        "download_any_other_processed_asset",
        "open_vcc_held_perturbation_truth",
        "model_fitting",
        "evaluation",
        "cell_eval2",
        "count_generation",
        "challenge_submission",
        "commit",
        "push",
    ):
        if authorization.get(forbidden) is not False:
            raise SourceAuditError(f"{forbidden} must remain false")

    assets = manifest.get("assets")
    if not isinstance(assets, dict) or not assets:
        raise SourceAuditError("manifest assets are missing")
    for name, asset in assets.items():
        group = asset.get("group")
        if group not in {"metadata", "large_candidate", "denied_large"}:
            raise SourceAuditError(f"{name}: invalid asset group")
        size = asset.get("size_bytes")
        if not isinstance(size, int) or size <= 0:
            raise SourceAuditError(f"{name}: invalid size")
        if not isinstance(asset.get("url"), str):
            raise SourceAuditError(f"{name}: URL is missing")
        if group == "metadata":
            relative = asset.get("local_path")
            if not isinstance(relative, str):
                raise SourceAuditError(f"{name}: metadata local_path is missing")
            destination = (HERE / relative).resolve()
            if not destination.is_relative_to(METADATA_ROOT.resolve()):
                raise SourceAuditError(f"{name}: metadata path escapes experiment")
            digest = asset.get("sha256")
            if not isinstance(digest, str) or len(digest) != 64:
                raise SourceAuditError(f"{name}: metadata SHA-256 is missing")
        else:
            if asset.get("local_path") is not None:
                raise SourceAuditError(f"{name}: large asset may not have a local path")
            if asset.get("download_authorized") is not False:
                raise SourceAuditError(f"{name}: large download must remain locked")

    if assets["geo_raw_archive"]["group"] != "denied_large":
        raise SourceAuditError("raw archive must remain denied")
    if {
        name for name, asset in assets.items() if asset["group"] == "large_candidate"
    } != LARGE_CANDIDATES:
        raise SourceAuditError("only Low1 and Low2 may be large candidates")
    return manifest


def is_large_download_authorized(
    manifest: dict[str, Any],
    asset_name: str,
) -> bool:
    """Return true only for an allowlisted, fully hash-bound, explicitly enabled asset."""
    if asset_name not in LARGE_CANDIDATES:
        return False
    asset = manifest["assets"].get(asset_name, {})
    digest = asset.get("sha256")
    return (
        asset.get("group") == "large_candidate"
        and asset.get("download_authorized") is True
        and isinstance(digest, str)
        and re.fullmatch(r"[0-9a-f]{64}", digest) is not None
        and asset_name in manifest["authorization"].get("authorized_large_assets", [])
    )


def verify_small_asset(path: Path, asset: dict[str, Any]) -> None:
    if path.is_symlink() or not path.is_file():
        raise SourceAuditError(f"missing or symlinked metadata asset: {path}")
    if path.stat().st_size != asset["size_bytes"]:
        raise SourceAuditError(f"{path.name}: byte count mismatch")
    if sha256_file(path) != asset["sha256"]:
        raise SourceAuditError(f"{path.name}: SHA-256 mismatch")


def acquire_small_metadata(manifest: dict[str, Any]) -> list[Path]:
    acquired: list[Path] = []
    for name, asset in manifest["assets"].items():
        if asset["group"] != "metadata":
            continue
        destination = (HERE / asset["local_path"]).resolve()
        if destination.exists():
            verify_small_asset(destination, asset)
            acquired.append(destination)
            continue
        request = urllib.request.Request(
            asset["url"],
            headers={"User-Agent": USER_AGENT},
        )
        with urllib.request.urlopen(request, timeout=60) as response:
            declared = response.headers.get("Content-Length")
            if declared is not None and int(declared) != asset["size_bytes"]:
                raise SourceAuditError(f"{name}: remote byte count drift")
            payload = response.read(asset["size_bytes"] + 1)
        if len(payload) != asset["size_bytes"]:
            raise SourceAuditError(f"{name}: downloaded byte count mismatch")
        if sha256_bytes(payload) != asset["sha256"]:
            raise SourceAuditError(f"{name}: downloaded SHA-256 mismatch")
        atomic_write(destination, payload)
        acquired.append(destination)
    return acquired


def head_binding(asset: dict[str, Any]) -> dict[str, Any]:
    request = urllib.request.Request(
        asset["url"],
        method="HEAD",
        headers={"User-Agent": USER_AGENT},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        headers = response.headers
        result = {
            "status": response.status,
            "content_length": int(headers["Content-Length"]),
            "last_modified": headers.get("Last-Modified"),
            "etag": headers.get("ETag"),
            "accept_ranges": headers.get("Accept-Ranges"),
        }
    if result["status"] != 200:
        raise SourceAuditError(f"{asset['filename']}: HEAD status changed")
    if result["content_length"] != asset["size_bytes"]:
        raise SourceAuditError(f"{asset['filename']}: remote byte count drift")
    if result["last_modified"] != asset["last_modified_http"]:
        raise SourceAuditError(f"{asset['filename']}: Last-Modified drift")
    return result


def parse_gzip_csv_header(prefix: bytes) -> tuple[bytes, list[str]]:
    decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
    try:
        plain = decompressor.decompress(prefix)
    except zlib.error as error:
        raise SourceAuditError("invalid gzip prefix") from error
    newline = plain.find(b"\n")
    if newline < 0:
        raise SourceAuditError("bounded prefix does not contain the CSV header")
    header = plain[: newline + 1]
    try:
        columns = next(csv.reader([header.decode("utf-8").rstrip("\r\n")]))
    except (UnicodeDecodeError, csv.Error) as error:
        raise SourceAuditError("invalid UTF-8 CSV header") from error
    if not columns or columns[0] != "gene":
        raise SourceAuditError("LFC matrix must start with a gene column")
    labels = columns[1:]
    if not labels or len(labels) != len(set(labels)):
        raise SourceAuditError("LFC target labels are empty or duplicated")
    return header, labels


def bounded_header_probe(asset: dict[str, Any]) -> tuple[dict[str, Any], list[str]]:
    cap = asset["bounded_header_probe_bytes"]
    request = urllib.request.Request(
        asset["url"],
        headers={
            "Range": f"bytes=0-{cap - 1}",
            "User-Agent": USER_AGENT,
            "Accept-Encoding": "identity",
        },
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        content_range = response.headers.get("Content-Range")
        payload = response.read(cap + 1)
        status = response.status
    if status != 206:
        raise SourceAuditError(f"{asset['filename']}: server ignored bounded range")
    if len(payload) != cap:
        raise SourceAuditError(f"{asset['filename']}: unexpected range byte count")
    expected_range = f"bytes 0-{cap - 1}/{asset['size_bytes']}"
    if content_range != expected_range:
        raise SourceAuditError(f"{asset['filename']}: Content-Range drift")
    header, labels = parse_gzip_csv_header(payload)
    labels_payload = ("\n".join(labels) + "\n").encode("utf-8")
    result = {
        "range_status": status,
        "content_range": content_range,
        "compressed_bytes_read": len(payload),
        "header_line_bytes": len(header),
        "header_sha256": sha256_bytes(header),
        "target_label_count": len(labels),
        "target_labels_sha256": sha256_bytes(labels_payload),
    }
    for key in (
        "header_line_bytes",
        "header_sha256",
        "target_label_count",
        "target_labels_sha256",
    ):
        if result[key] != asset[key]:
            raise SourceAuditError(
                f"{asset['filename']}: bounded header drift for {key}"
            )
    return result, labels


def parse_family_soft(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    fields: dict[str, Any] = {}
    for key in (
        "Series_title",
        "Series_geo_accession",
        "Series_status",
        "Series_submission_date",
        "Series_last_update_date",
        "Series_summary",
        "Series_overall_design",
    ):
        match = re.search(rf"^!{key} = (.*)$", text, re.MULTILINE)
        if match is None:
            raise SourceAuditError(f"SOFT metadata is missing {key}")
        fields[key] = match.group(1)
    fields["sample_count"] = text.count("^SAMPLE = ")
    return fields


def parse_geo_filelist(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    archives = [row for row in rows if row["#Archive/File"] == "Archive"]
    if len(archives) != 1 or archives[0]["Name"] != "GSE337988_RAW.tar":
        raise SourceAuditError("GEO inventory has an unexpected archive record")
    files = [row for row in rows if row["#Archive/File"] == "File"]
    return {
        "raw_archive_size_bytes": int(archives[0]["Size"]),
        "raw_member_count": len(files),
        "raw_member_bytes": sum(int(row["Size"]) for row in files),
    }


def parse_guide_library(path: Path) -> dict[str, Any]:
    guide_labels: set[str] = set()
    promoter_labels: set[str] = set()
    row_count = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle):
            if len(row) != 3 or row[2] != "crispr":
                raise SourceAuditError("unexpected CRISPR library row")
            sequence, label, _ = row
            if not sequence or not label or label in guide_labels:
                raise SourceAuditError("empty or duplicate CRISPR guide label")
            guide_labels.add(label)
            match = re.fullmatch(r"(.+)_([1-4])", label)
            promoter_labels.add(match.group(1) if match else label)
            row_count += 1
    return {
        "row_count": row_count,
        "unique_guide_labels": len(guide_labels),
        "promoter_labels": sorted(promoter_labels),
    }


def read_vcc_panel(manifest: dict[str, Any]) -> list[str]:
    axis = manifest["local_axes"]["vcc_2026_panel"]
    path = ROOT / axis["path"]
    if path.is_symlink() or not path.is_file():
        raise SourceAuditError(f"missing or symlinked VCC panel: {path}")
    if sha256_file(path) != axis["sha256"]:
        raise SourceAuditError("VCC panel SHA-256 drift")
    with path.open(encoding="utf-8", newline="") as handle:
        labels = [row[axis["column"]].strip() for row in csv.DictReader(handle)]
    if len(labels) != axis["row_count"] or len(labels) != len(set(labels)):
        raise SourceAuditError("VCC panel row-count or uniqueness drift")
    return labels


def panel_overlap(
    panel: list[str],
    target_labels: list[str],
) -> tuple[list[str], dict[str, list[str]]]:
    matches = {
        gene: sorted(
            label
            for label in target_labels
            if label == gene or label.startswith(f"{gene}_")
        )
        for gene in panel
    }
    overlap = sorted(gene for gene, labels in matches.items() if labels)
    multi = {gene: labels for gene, labels in matches.items() if len(labels) > 1}
    return overlap, multi


def require_expected(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    for key, expected_value in expected.items():
        if actual.get(key) != expected_value:
            raise SourceAuditError(
                f"audit drift for {key}: {actual.get(key)!r} != {expected_value!r}"
            )


def build_report(*, live_probe: bool) -> dict[str, Any]:
    manifest = load_manifest()
    metadata_assets = {
        name: asset
        for name, asset in manifest["assets"].items()
        if asset["group"] == "metadata"
    }
    for asset in metadata_assets.values():
        verify_small_asset(HERE / asset["local_path"], asset)

    soft = parse_family_soft(HERE / metadata_assets["geo_family_soft"]["local_path"])
    inventory = parse_geo_filelist(
        HERE / metadata_assets["geo_file_inventory"]["local_path"]
    )
    library = parse_guide_library(
        HERE / metadata_assets["scaled_crispr_library"]["local_path"]
    )

    probes: dict[str, Any] = {}
    if live_probe:
        labels_by_asset: dict[str, list[str]] = {}
        for name in sorted(LARGE_CANDIDATES):
            asset = manifest["assets"][name]
            head = head_binding(asset)
            header, labels = bounded_header_probe(asset)
            probes[name] = {"head": head, "bounded_header": header}
            labels_by_asset[name] = labels
        if labels_by_asset["low1_lfc_matrix"] != labels_by_asset["low2_lfc_matrix"]:
            raise SourceAuditError("Low1 and Low2 target-label headers differ")
        target_labels = labels_by_asset["low1_lfc_matrix"]
        atomic_write(
            TARGET_LABELS_PATH,
            ("\n".join(target_labels) + "\n").encode("utf-8"),
        )
        probes["geo_raw_archive"] = {
            "head": head_binding(manifest["assets"]["geo_raw_archive"])
        }
    else:
        derived = manifest["derived_metadata"]["target_labels"]
        verify_small_asset(TARGET_LABELS_PATH, derived)
        target_labels = TARGET_LABELS_PATH.read_text(encoding="utf-8").splitlines()

    panel = read_vcc_panel(manifest)
    overlap, multi = panel_overlap(panel, target_labels)
    target_gene_prefixes = {
        label.split("_", 1)[0]
        for label in target_labels
        if not label.startswith("NTC_")
    }
    summary = {
        "geo_sample_count": soft["sample_count"],
        "raw_archive_size_bytes": inventory["raw_archive_size_bytes"],
        "guide_library_rows": library["row_count"],
        "guide_library_unique_labels": library["unique_guide_labels"],
        "guide_library_promoter_labels": len(library["promoter_labels"]),
        "lfc_target_labels": len(target_labels),
        "lfc_headers_identical": True,
        "lfc_unique_target_labels": len(set(target_labels)),
        "lfc_non_ntc_target_labels": sum(
            not label.startswith("NTC_") for label in target_labels
        ),
        "lfc_unique_non_ntc_gene_prefixes": len(target_gene_prefixes),
        "vcc_panel_gene_overlap": len(overlap),
        "vcc_panel_promoter_label_overlap": sum(
            1
            for label in target_labels
            if any(label == gene or label.startswith(f"{gene}_") for gene in panel)
        ),
        "vcc_panel_overlap_genes": overlap,
        "vcc_panel_multi_promoter_labels": multi,
    }
    require_expected(summary, manifest["expected_audit"])

    return {
        "schema_version": 1,
        "experiment": manifest["experiment"],
        "dataset": manifest["dataset"],
        "audited_at": manifest["recorded_at"],
        "status": "source_metadata_verified_large_assets_locked",
        "source_manifest_sha256": sha256_file(MANIFEST_PATH),
        "geo_series": soft,
        "inventory": inventory,
        "remote_probes": probes,
        "summary": summary,
        "authorization": {
            "authorized_large_assets": [
                name
                for name in sorted(LARGE_CANDIDATES)
                if is_large_download_authorized(manifest, name)
            ],
            "raw_archive_authorized": False,
            "reason": manifest["large_asset_authorization_contract"]["reason"],
        },
        "interpretation": {
            "held_truth_opened": False,
            "large_file_downloaded": False,
            "model_fit_or_evaluation_performed": False,
            "overlap_unit": "VCC gene symbols with one or more exact-prefix official LFC target labels",
            "overlap_is_post_effect_qc": False,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser(
        "acquire-metadata",
        help="download only the three small hash-pinned metadata files",
    )
    audit_parser = subparsers.add_parser(
        "audit",
        help="verify metadata and recompute the target overlap",
    )
    audit_parser.add_argument(
        "--offline",
        action="store_true",
        help="use the previously frozen target-label file instead of HTTP ranges",
    )
    audit_parser.add_argument(
        "--write-report",
        action="store_true",
        help="write METADATA_AUDIT.json inside this experiment",
    )
    args = parser.parse_args()

    manifest = load_manifest()
    if args.command == "acquire-metadata":
        for path in acquire_small_metadata(manifest):
            print(path)
        return
    report = build_report(live_probe=not args.offline)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.write_report:
        atomic_write(REPORT_PATH, rendered.encode("utf-8"))
        print(REPORT_PATH)
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
