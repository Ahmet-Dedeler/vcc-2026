#!/usr/bin/env python3
"""Reproduce the metadata-only GSE295214 inventory and overlap audit."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import re
import zipfile
from pathlib import Path
from typing import Any, Iterable
from xml.etree import ElementTree

from acquire import AcquisitionError, DATA_ROOT, HERE, ROOT, load_manifest, verify_asset


REPORT_PATH = HERE / "METADATA_AUDIT.json"
TARGET_PATTERN = re.compile(r"^(?P<gene>.+?)_[+-]_.+?_f_(?P<copy>[0-9]+)_left$")
CELL_REFERENCE = re.compile(r"^([A-Z]+)([0-9]+)$")
XML_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
XML_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
XML_PACKAGE_REL = "http://schemas.openxmlformats.org/package/2006/relationships"


def sha256_file(path: Path, block_bytes: int = 8 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(block_bytes):
            digest.update(block)
    return digest.hexdigest()


def column_index(reference: str) -> int:
    match = CELL_REFERENCE.fullmatch(reference)
    if match is None:
        raise AcquisitionError(f"invalid XLSX cell reference: {reference}")
    result = 0
    for character in match.group(1):
        result = result * 26 + ord(character) - ord("A") + 1
    return result - 1


def first_sheet_rows(path: Path) -> Iterable[list[str]]:
    """Read values from the first XLSX worksheet using only the standard library."""
    with zipfile.ZipFile(path) as archive:
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        sheets = workbook.find(f"{{{XML_MAIN}}}sheets")
        if sheets is None or len(sheets) == 0:
            raise AcquisitionError(f"{path.name}: workbook has no sheets")
        relationship_id = sheets[0].attrib[f"{{{XML_REL}}}id"]
        relationships = ElementTree.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        target = None
        for relationship in relationships:
            if relationship.attrib.get("Id") == relationship_id:
                target = relationship.attrib["Target"]
                break
        if target is None:
            raise AcquisitionError(f"{path.name}: first-sheet relationship is missing")
        sheet_path = target.lstrip("/")
        if not sheet_path.startswith("xl/"):
            sheet_path = f"xl/{sheet_path}"

        shared: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared_root.findall(f"{{{XML_MAIN}}}si"):
                shared.append(
                    "".join(node.text or "" for node in item.iter(f"{{{XML_MAIN}}}t"))
                )

        sheet = ElementTree.fromstring(archive.read(sheet_path))
        sheet_data = sheet.find(f"{{{XML_MAIN}}}sheetData")
        if sheet_data is None:
            raise AcquisitionError(f"{path.name}: sheetData is missing")
        for row in sheet_data:
            values: dict[int, str] = {}
            for cell in row.findall(f"{{{XML_MAIN}}}c"):
                index = column_index(cell.attrib["r"])
                cell_type = cell.attrib.get("t")
                value_node = cell.find(f"{{{XML_MAIN}}}v")
                if cell_type == "inlineStr":
                    inline = cell.find(f"{{{XML_MAIN}}}is")
                    value = (
                        ""
                        if inline is None
                        else "".join(
                            node.text or "" for node in inline.iter(f"{{{XML_MAIN}}}t")
                        )
                    )
                elif value_node is None:
                    value = ""
                elif cell_type == "s":
                    value = shared[int(value_node.text or "0")]
                else:
                    value = value_node.text or ""
                values[index] = value
            if values:
                width = max(values) + 1
                yield [values.get(index, "") for index in range(width)]


def parse_library_workbook(path: Path) -> dict[str, Any]:
    rows = iter(first_sheet_rows(path))
    try:
        header = next(rows)
    except StopIteration as error:
        raise AcquisitionError("library workbook is empty") from error
    if header[:2] != ["element", "sequence"]:
        raise AcquisitionError(f"unexpected library header: {header[:2]}")
    elements: set[str] = set()
    sequences: set[str] = set()
    target_sequences: set[str] = set()
    controls: set[str] = set()
    targets: set[str] = set()
    row_count = 0
    for row in rows:
        if len(row) < 2 or not row[0]:
            continue
        element, sequence = row[:2]
        row_count += 1
        if element in elements:
            raise AcquisitionError(f"duplicate oligo element: {element}")
        elements.add(element)
        sequences.add(sequence)
        if element.startswith("non-targeting_"):
            controls.add(element.split("_f_", 1)[0])
            continue
        match = TARGET_PATTERN.fullmatch(element)
        if match is None:
            raise AcquisitionError(f"unrecognized oligo element: {element}")
        targets.add(match.group("gene"))
        target_sequences.add(sequence)
    return {
        "row_count": row_count,
        "unique_oligo_sequences": len(sequences),
        "unique_targeting_sequences": len(target_sequences),
        "unique_negative_controls": len(controls),
        "literal_target_labels": sorted(targets),
    }


def parse_geo_filelist(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    archives = [row for row in rows if row["#Archive/File"] == "Archive"]
    files = [row for row in rows if row["#Archive/File"] == "File"]
    if len(archives) != 1 or archives[0]["Name"] != "GSE295214_RAW.tar":
        raise AcquisitionError("GEO file list has an unexpected archive record")
    if any(not row["Name"].endswith(".h5") for row in files):
        raise AcquisitionError("GEO file list contains a non-H5 member")
    return {
        "archive_size_bytes": int(archives[0]["Size"]),
        "h5_member_count": len(files),
        "h5_member_bytes": sum(int(row["Size"]) for row in files),
        "h5_members": [
            {"name": row["Name"], "size_bytes": int(row["Size"])} for row in files
        ],
    }


def parse_family_soft(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        text = handle.read()
    fields: dict[str, str] = {}
    for key in (
        "Series_title",
        "Series_geo_accession",
        "Series_status",
        "Series_submission_date",
        "Series_last_update_date",
        "Series_web_link",
        "Series_overall_design",
    ):
        match = re.search(rf"^!{key} = (.*)$", text, re.MULTILINE)
        if match is None:
            raise AcquisitionError(f"SOFT metadata is missing {key}")
        fields[key] = match.group(1)
    fields["sample_count"] = str(text.count("^SAMPLE = "))
    fields["sample_h5_url_count"] = str(
        len(
            re.findall(
                r"^!Sample_supplementary_file_1 = .*[.]h5$",
                text,
                re.MULTILINE,
            )
        )
    )
    return fields


def read_axis(path: Path, column: str | None = None) -> list[str]:
    if column is None:
        values = [
            line.strip() for line in path.read_text().splitlines() if line.strip()
        ]
    else:
        with path.open(encoding="utf-8", newline="") as handle:
            values = [row[column].strip() for row in csv.DictReader(handle)]
    if not values or len(values) != len(set(values)):
        raise AcquisitionError(f"{path}: axis is empty or contains duplicates")
    return values


def require_expected(actual: dict[str, Any], expected: dict[str, Any]) -> None:
    for key, expected_value in expected.items():
        if actual.get(key) != expected_value:
            raise AcquisitionError(
                f"audit drift for {key}: {actual.get(key)!r} != {expected_value!r}"
            )


def build_report() -> dict[str, Any]:
    manifest = load_manifest()
    for asset in manifest["assets"].values():
        if asset["group"] in {"metadata", "derived_metadata"}:
            verify_asset(ROOT / asset["local_path"], asset)
    for axis in manifest["local_axes"].values():
        path = ROOT / axis["path"]
        if path.is_symlink() or not path.is_file():
            raise AcquisitionError(f"missing or symlinked local axis: {path}")
        if sha256_file(path) != axis["sha256"]:
            raise AcquisitionError(f"local axis hash drift: {path}")

    inventory = parse_geo_filelist(
        DATA_ROOT / manifest["assets"]["geo_file_inventory"]["filename"]
    )
    soft = parse_family_soft(
        DATA_ROOT / manifest["assets"]["geo_family_soft"]["filename"]
    )
    library = parse_library_workbook(
        DATA_ROOT / manifest["assets"]["supplementary_table_6"]["filename"]
    )
    panel = set(
        read_axis(
            ROOT / manifest["local_axes"]["vcc_2026_panel"]["path"], "target_gene"
        )
    )
    h1 = set(
        read_axis(ROOT / manifest["local_axes"]["vcc_2025_h1_donor_targets"]["path"])
    )
    retained = set(
        read_axis(DATA_ROOT / manifest["assets"]["post_qc_target_axis"]["filename"])
    )
    library_targets = set(library.pop("literal_target_labels"))

    library_panel = sorted(library_targets & panel)
    retained_panel = sorted(retained & panel)
    retained_panel_new = sorted((retained & panel) - h1)
    retained_h1 = sorted(retained & h1)
    failed_panel = sorted((library_targets & panel) - retained)

    summary = {
        "geo_h5_member_count": inventory["h5_member_count"],
        "geo_h5_member_bytes": inventory["h5_member_bytes"],
        "library_workbook_rows": library["row_count"],
        "library_unique_oligo_sequences": library["unique_oligo_sequences"],
        "library_unique_negative_controls": library["unique_negative_controls"],
        "library_literal_target_labels": len(library_targets),
        "post_qc_targets": len(retained),
        "vcc_panel_library_overlap": len(library_panel),
        "vcc_panel_post_qc_overlap": len(retained_panel),
        "vcc_panel_post_qc_new_vs_vcc2025_h1": len(retained_panel_new),
        "vcc2025_h1_post_qc_overlap": len(retained_h1),
        "vcc_panel_library_targets_removed_by_qc": failed_panel,
    }
    require_expected(summary, manifest["expected_audit"])
    if (
        inventory["archive_size_bytes"]
        != manifest["assets"]["geo_raw_archive"]["size_bytes"]
    ):
        raise AcquisitionError("raw archive size drift")
    if int(soft["sample_count"]) != manifest["official_sources"]["geo"]["sample_count"]:
        raise AcquisitionError("SOFT sample count drift")
    if int(soft["sample_h5_url_count"]) != inventory["h5_member_count"]:
        raise AcquisitionError("SOFT and file-list sample inventories disagree")

    return {
        "schema_version": 1,
        "dataset": "GSE295214",
        "status": "metadata_audit_complete_large_download_locked",
        "manifest_sha256": sha256_file(HERE / "SOURCE_MANIFEST.json"),
        "summary": summary,
        "overlaps": {
            "vcc_panel_library": library_panel,
            "vcc_panel_post_qc": retained_panel,
            "vcc_panel_post_qc_new_vs_vcc2025_h1": retained_panel_new,
            "vcc2025_h1_post_qc": retained_h1,
        },
        "inventory": inventory,
        "geo_series": soft,
        "interpretation": {
            "large_source_recommendation": "acquire",
            "reason": (
                "64 VCC-panel targets survive the authors' final H1 QC, including "
                "57 targets absent from the existing VCC-2025 H1 donor axis."
            ),
            "full_download_performed": False,
            "held_truth_opened": False,
        },
    }


def main() -> None:
    report = build_report()
    REPORT_PATH.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
