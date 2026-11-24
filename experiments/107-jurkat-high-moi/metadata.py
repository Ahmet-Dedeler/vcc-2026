#!/usr/bin/env python3
"""Deterministic metadata crosswalks for experiment 107.

This module reads only metadata: the pinned guide FASTA, workbook tables,
challenge CSV axes, and compressed feature tables inside the source archive.
It never opens a Matrix Market coordinate stream.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
import tarfile
import zipfile
from collections import Counter, defaultdict
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import TypeAlias
from xml.etree import ElementTree

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

GUIDE_FASTA = HERE / "downloads" / "GSE249595_guides_reference.fa.gz"
WORKBOOK = HERE / "downloads" / "MOESM3_ESM.xlsx"
RAW_ARCHIVE = HERE / "downloads" / "GSE249595_RAW.tar"
CHALLENGE_TARGETS = ROOT / "data" / "challenge" / "pert_counts.csv"
CHALLENGE_GENES = ROOT / "data" / "challenge" / "gene_names.csv"

PINNED_SHA256 = {
    GUIDE_FASTA: "99e99d6c14907a2cc60c4c11ed3d93e9da5cc113740489b7e8d1e8c0abe00ca6",
    WORKBOOK: "307f346279b4bda2e9c19aca39c692a144197c0117dfd9d595724cc353e2f7ff",
    CHALLENGE_TARGETS: (
        "f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276"
    ),
    CHALLENGE_GENES: (
        "25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201"
    ),
}

EXPECTED_GUIDES = 83_401
EXPECTED_NON_TARGETING = 250
EXPECTED_SOURCE_TARGET_ROWS = 18_466
EXPECTED_SOURCE_TARGET_SYMBOLS = 18_465
EXPECTED_RESPONSE_GENES = 374
EXPECTED_HASH_LABELS = 7
EXPECTED_PANEL_TARGETS = 300
NO_UNIQUE_TARGET = -1

GUIDE_SUPPORTED = "supported"
GUIDE_NON_TARGETING = "non_targeting"
GUIDE_COMPOSITE = "composite"
GUIDE_UNMAPPED = "unmapped"
GUIDE_AMBIGUOUS = "ambiguous"

_GUIDE_ID = re.compile(r"^(?P<label>.+)_(?P<ordinal>[1-9][0-9]*)$")
_HASH_LABEL = re.compile(r"^(?P<condition>activated|untreated)(?P<number>[0-9]*)_1$")
_FEATURE_MEMBER = re.compile(
    r"^(?P<accession>GSM[0-9]+)_channel(?P<channel>[0-9]+)_"
    r"(?P<modality>transcriptome|guides|labels)_features\.tsv\.gz$"
)
_CELL_REFERENCE = re.compile(r"^(?P<column>[A-Z]+)(?P<row>[1-9][0-9]*)$")

_SPREADSHEET_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
_RELATIONSHIP_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
_PACKAGE_RELATIONSHIP_NS = (
    "http://schemas.openxmlformats.org/package/2006/relationships"
)

CellValue: TypeAlias = str | int | float | bool | None


class MetadataError(RuntimeError):
    """A source metadata invariant failed."""


@dataclass(frozen=True)
class TargetRecord:
    """One populated row from supplementary Table S1."""

    symbol: str
    entrez_id: int
    workbook_row: int


@dataclass(frozen=True)
class ResponseRecord:
    """One integer-Entrez gene row from supplementary Table S2."""

    symbol: str
    entrez_id: int
    group: str
    workbook_row: int


@dataclass(frozen=True)
class HashOligo:
    """One hash oligo and its condition from supplementary Table S6."""

    oligo: str
    condition: str
    workbook_row: int


@dataclass(frozen=True)
class WorkbookMetadata:
    """The three workbook tables used by experiment 107."""

    targets: tuple[TargetRecord, ...]
    responses: tuple[ResponseRecord, ...]
    hash_oligos: tuple[HashOligo, ...]
    ignored_response_rows: tuple[tuple[int, tuple[CellValue, ...]], ...]


@dataclass(frozen=True)
class FeatureAxis:
    """A 10x feature axis, without any matrix coordinates."""

    identifiers: tuple[str, ...]
    symbols: tuple[str, ...]
    feature_types: tuple[str, ...]


@dataclass(frozen=True)
class ArchiveFeatureAxes:
    """Guide, hash-label, and transcriptome feature axes for one channel."""

    channel: int
    accession: str
    guides: FeatureAxis
    labels: FeatureAxis
    transcriptome: FeatureAxis


@dataclass(frozen=True)
class DuplicateSymbol:
    """A symbol associated with more than one integer Entrez identifier."""

    symbol: str
    entrez_ids: tuple[int, ...]
    workbook_rows: tuple[int, ...]


@dataclass(frozen=True)
class GuideCrosswalk:
    """Canonical guide axis and its target-level design mapping."""

    guide_ids: tuple[str, ...]
    guide_labels: tuple[str, ...]
    guide_status: tuple[str, ...]
    target_symbols: tuple[str, ...]
    target_entrez_ids: tuple[tuple[int, ...], ...]
    guide_to_target: np.ndarray
    non_targeting_guides: np.ndarray
    duplicate_target_symbols: tuple[DuplicateSymbol, ...]
    composite_labels: tuple[str, ...]
    unmapped_labels: tuple[str, ...]
    ambiguous_labels: tuple[str, ...]
    panel_targets: tuple[str, ...]
    panel_covered_targets: tuple[str, ...]
    panel_missing_targets: tuple[str, ...]

    @property
    def status_counts(self) -> dict[str, int]:
        """Return deterministic guide-record counts by classification."""

        return dict(sorted(Counter(self.guide_status).items()))


@dataclass(frozen=True)
class HashLabelRecord:
    """One archive hash feature mapped to a workbook-backed condition."""

    feature_index: int
    label: str
    condition: str


@dataclass(frozen=True)
class HashCrosswalk:
    """The seven hash labels and the condition inventory that validates them."""

    records: tuple[HashLabelRecord, ...]
    workbook_oligos: tuple[HashOligo, ...]

    @property
    def label_to_condition(self) -> dict[str, str]:
        return {record.label: record.condition for record in self.records}


@dataclass(frozen=True)
class ResponseGeneCrosswalk:
    """One targeted response gene mapped onto source and VCC feature axes."""

    response_index: int
    symbol: str
    entrez_id: int
    group: str
    source_feature_indices: tuple[int, ...]
    vcc_gene_indices: tuple[int, ...]
    is_panel_target: bool

    @property
    def vcc_gene_index(self) -> int | None:
        """Return the unique VCC index, or ``None`` if absent/ambiguous."""

        if len(self.vcc_gene_indices) == 1:
            return self.vcc_gene_indices[0]
        return None


@dataclass(frozen=True)
class ResponseCrosswalk:
    """The 374 response genes plus explicit source/VCC ambiguity reports."""

    records: tuple[ResponseGeneCrosswalk, ...]
    duplicate_response_symbols: tuple[str, ...]
    duplicate_response_entrez_ids: tuple[int, ...]
    missing_source_symbols: tuple[str, ...]
    ambiguous_source_symbols: tuple[str, ...]
    missing_vcc_symbols: tuple[str, ...]
    ambiguous_vcc_symbols: tuple[str, ...]
    visible_vcc_indices: tuple[int, ...]


@dataclass(frozen=True)
class MetadataCrosswalks:
    """All metadata products needed by later source-processing stages."""

    archive_axes: ArchiveFeatureAxes
    guide: GuideCrosswalk
    hashes: HashCrosswalk
    response: ResponseCrosswalk


def sha256_file(path: Path, block_bytes: int = 1 << 20) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(block_bytes), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_pinned_file(path: Path) -> None:
    """Require a regular file with the preregistered SHA-256."""

    expected = PINNED_SHA256.get(path)
    if expected is None:
        raise MetadataError(f"no pinned SHA-256 is registered for {path}")
    if not path.is_file() or path.is_symlink():
        raise MetadataError(f"missing or unsafe pinned metadata file: {path}")
    observed = sha256_file(path)
    if observed != expected:
        raise MetadataError(
            f"pinned metadata hash mismatch for {path.name}: "
            f"expected {expected}, observed {observed}"
        )


def _column_index(reference: str) -> int:
    match = _CELL_REFERENCE.fullmatch(reference)
    if match is None:
        raise MetadataError(f"malformed spreadsheet cell reference: {reference!r}")
    index = 0
    for letter in match.group("column"):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def _safe_xlsx_target(target: str) -> str:
    path = PurePosixPath("xl") / target
    if path.is_absolute() or ".." in path.parts:
        raise MetadataError(f"unsafe workbook relationship target: {target!r}")
    return str(path)


def _shared_strings(archive: zipfile.ZipFile) -> tuple[str, ...]:
    try:
        root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    except KeyError:
        return ()
    except ElementTree.ParseError as error:
        raise MetadataError("malformed workbook shared strings") from error
    namespace = {"x": _SPREADSHEET_NS}
    return tuple(
        "".join(node.text or "" for node in item.findall(".//x:t", namespace))
        for item in root.findall("x:si", namespace)
    )


def _sheet_path(archive: zipfile.ZipFile, sheet_name: str) -> str:
    try:
        workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
        relationships = ElementTree.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
    except (KeyError, ElementTree.ParseError) as error:
        raise MetadataError("malformed workbook structure") from error

    relation_targets = {
        relation.attrib["Id"]: relation.attrib["Target"]
        for relation in relationships.findall(
            f"{{{_PACKAGE_RELATIONSHIP_NS}}}Relationship"
        )
        if "Id" in relation.attrib and "Target" in relation.attrib
    }
    relationship_key = f"{{{_RELATIONSHIP_NS}}}id"
    for sheet in workbook.findall(
        f".//{{{_SPREADSHEET_NS}}}sheet",
    ):
        if sheet.attrib.get("name") != sheet_name:
            continue
        relation_id = sheet.attrib.get(relationship_key)
        if relation_id not in relation_targets:
            raise MetadataError(f"worksheet {sheet_name!r} has no relationship")
        return _safe_xlsx_target(relation_targets[relation_id])
    raise MetadataError(f"workbook is missing worksheet {sheet_name!r}")


def _cell_value(
    cell: ElementTree.Element,
    shared_strings: Sequence[str],
) -> CellValue:
    cell_type = cell.attrib.get("t")
    value = cell.find(f"{{{_SPREADSHEET_NS}}}v")
    raw = None if value is None else value.text

    if cell_type == "inlineStr":
        return "".join(
            node.text or ""
            for node in cell.findall(
                f".//{{{_SPREADSHEET_NS}}}t",
            )
        )
    if raw is None:
        return None
    if cell_type == "s":
        try:
            return shared_strings[int(raw)]
        except (ValueError, IndexError) as error:
            raise MetadataError("invalid workbook shared-string index") from error
    if cell_type == "b":
        if raw not in {"0", "1"}:
            raise MetadataError(f"invalid workbook boolean value: {raw!r}")
        return raw == "1"
    if cell_type in {"str", "e"}:
        return raw
    if cell_type not in {None, "n"}:
        raise MetadataError(f"unsupported workbook cell type: {cell_type!r}")
    try:
        numeric = float(raw)
    except ValueError as error:
        raise MetadataError(f"invalid workbook numeric value: {raw!r}") from error
    if not np.isfinite(numeric):
        raise MetadataError(f"non-finite workbook numeric value: {raw!r}")
    if numeric.is_integer():
        return int(numeric)
    return numeric


def read_xlsx_sheet(
    path: Path,
    sheet_name: str,
) -> tuple[tuple[int, tuple[CellValue, ...]], ...]:
    """Read one XLSX worksheet with only Python's standard library."""

    try:
        with zipfile.ZipFile(path) as archive:
            strings = _shared_strings(archive)
            sheet_path = _sheet_path(archive, sheet_name)
            root = ElementTree.fromstring(archive.read(sheet_path))
    except (OSError, zipfile.BadZipFile, KeyError, ElementTree.ParseError) as error:
        raise MetadataError(
            f"cannot read worksheet {sheet_name!r} from {path}"
        ) from error

    rows: list[tuple[int, tuple[CellValue, ...]]] = []
    for row in root.findall(
        f".//{{{_SPREADSHEET_NS}}}sheetData/{{{_SPREADSHEET_NS}}}row"
    ):
        try:
            row_number = int(row.attrib["r"])
        except (KeyError, ValueError) as error:
            raise MetadataError("worksheet row has no valid row number") from error
        values: dict[int, CellValue] = {}
        for cell in row.findall(f"{{{_SPREADSHEET_NS}}}c"):
            reference = cell.attrib.get("r")
            if reference is None:
                raise MetadataError("worksheet cell has no reference")
            column = _column_index(reference)
            if column in values:
                raise MetadataError(f"duplicate worksheet cell: {reference}")
            values[column] = _cell_value(cell, strings)
        width = max(values, default=-1) + 1
        rows.append(
            (
                row_number,
                tuple(values.get(column) for column in range(width)),
            )
        )
    return tuple(rows)


def _find_header(
    rows: Sequence[tuple[int, tuple[CellValue, ...]]],
    expected: tuple[str, ...],
) -> int:
    matches = []
    for index, (_, values) in enumerate(rows):
        leading = tuple(values[: len(expected)])
        if leading == expected:
            matches.append(index)
    if len(matches) != 1:
        raise MetadataError(
            f"expected exactly one workbook header {expected!r}, found {len(matches)}"
        )
    return matches[0]


def _integer_entrez(value: CellValue) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int) and value > 0:
        return value
    return None


def read_workbook_metadata(path: Path) -> WorkbookMetadata:
    """Read and validate supplementary Tables S1, S2, and S6."""

    table_1 = read_xlsx_sheet(path, "Table_S1")
    table_2 = read_xlsx_sheet(path, "Table_S2")
    table_6 = read_xlsx_sheet(path, "Table_S6")

    target_header = _find_header(table_1, ("symbol", "entrez"))
    targets: list[TargetRecord] = []
    for row_number, values in table_1[target_header + 1 :]:
        padded = values + (None,) * max(0, 2 - len(values))
        if padded[0] is None and padded[1] is None:
            continue
        symbol, raw_entrez = padded[:2]
        entrez_id = _integer_entrez(raw_entrez)
        if not isinstance(symbol, str) or not symbol or entrez_id is None:
            raise MetadataError(f"Table_S1 row {row_number} is malformed")
        targets.append(TargetRecord(symbol, entrez_id, row_number))

    response_header = _find_header(table_2, ("symbol", "entrez", "group"))
    responses: list[ResponseRecord] = []
    ignored: list[tuple[int, tuple[CellValue, ...]]] = []
    for row_number, values in table_2[response_header + 1 :]:
        padded = values + (None,) * max(0, 3 - len(values))
        if all(value is None for value in padded):
            continue
        symbol, raw_entrez, group = padded[:3]
        entrez_id = _integer_entrez(raw_entrez)
        if entrez_id is None:
            ignored.append((row_number, values))
            continue
        if (
            not isinstance(symbol, str)
            or not symbol
            or not isinstance(group, str)
            or not group
        ):
            raise MetadataError(f"Table_S2 row {row_number} is malformed")
        responses.append(ResponseRecord(symbol, entrez_id, group, row_number))

    hash_header = _find_header(table_6, ("oligo", "condition"))
    hash_oligos: list[HashOligo] = []
    for row_number, values in table_6[hash_header + 1 :]:
        padded = values + (None,) * max(0, 2 - len(values))
        if padded[0] is None and padded[1] is None:
            continue
        oligo, condition = padded[:2]
        if (
            not isinstance(oligo, str)
            or not re.fullmatch(r"[ACGT]+", oligo)
            or condition not in {"untreated", "activated"}
        ):
            raise MetadataError(f"Table_S6 row {row_number} is malformed")
        hash_oligos.append(HashOligo(oligo, condition, row_number))

    return WorkbookMetadata(
        targets=tuple(targets),
        responses=tuple(responses),
        hash_oligos=tuple(hash_oligos),
        ignored_response_rows=tuple(ignored),
    )


def read_fasta_ids(path: Path) -> tuple[str, ...]:
    """Return FASTA record identifiers in source order."""

    opener = gzip.open if path.suffix == ".gz" else open
    identifiers: list[str] = []
    current_identifier: str | None = None
    current_sequence_length = 0
    try:
        with opener(path, "rt", encoding="utf-8", newline=None) as handle:
            for line_number, raw in enumerate(handle, start=1):
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if current_identifier is not None and current_sequence_length == 0:
                        raise MetadataError(
                            f"FASTA record {current_identifier!r} has no sequence"
                        )
                    identifier = line[1:].strip()
                    if not identifier or any(
                        character.isspace() for character in identifier
                    ):
                        raise MetadataError(
                            f"FASTA line {line_number} has an invalid identifier"
                        )
                    identifiers.append(identifier)
                    current_identifier = identifier
                    current_sequence_length = 0
                    continue
                if current_identifier is None:
                    raise MetadataError(
                        "FASTA sequence appears before the first header"
                    )
                if re.fullmatch(r"[ACGTN]+", line.upper()) is None:
                    raise MetadataError(
                        f"FASTA record {current_identifier!r} has invalid sequence text"
                    )
                current_sequence_length += len(line)
    except (OSError, UnicodeError) as error:
        raise MetadataError(f"cannot read FASTA {path}") from error

    if current_identifier is not None and current_sequence_length == 0:
        raise MetadataError(f"FASTA record {current_identifier!r} has no sequence")
    if not identifiers:
        raise MetadataError("FASTA contains no records")
    duplicates = _duplicates(identifiers)
    if duplicates:
        raise MetadataError(f"FASTA contains duplicate identifiers: {duplicates!r}")
    return tuple(identifiers)


def _duplicates(values: Iterable[str | int]) -> tuple[str | int, ...]:
    counts = Counter(values)
    return tuple(sorted(value for value, count in counts.items() if count > 1))


def read_single_column_csv(path: Path, column: str) -> tuple[str, ...]:
    """Read a named one-column CSV axis and require unique nonempty values."""

    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != [column]:
                raise MetadataError(
                    f"{path.name} must contain exactly the {column!r} column"
                )
            values = tuple(row[column] for row in reader)
    except (OSError, csv.Error, KeyError) as error:
        raise MetadataError(f"cannot read {path}") from error
    if any(not value for value in values):
        raise MetadataError(f"{path.name} contains an empty {column}")
    duplicates = _duplicates(values)
    if duplicates:
        raise MetadataError(f"{path.name} contains duplicate values: {duplicates!r}")
    return values


def _read_feature_member(
    archive: tarfile.TarFile,
    member: tarfile.TarInfo,
) -> FeatureAxis:
    if not member.isfile():
        raise MetadataError(f"archive feature member is not regular: {member.name}")
    extracted = archive.extractfile(member)
    if extracted is None:
        raise MetadataError(f"cannot open archive feature member: {member.name}")
    identifiers: list[str] = []
    symbols: list[str] = []
    feature_types: list[str] = []
    try:
        with extracted, gzip.GzipFile(fileobj=extracted, mode="rb") as stream:
            for line_number, raw in enumerate(stream, start=1):
                fields = raw.decode("utf-8").rstrip("\r\n").split("\t")
                if len(fields) != 3 or any(not field for field in fields):
                    raise MetadataError(
                        f"{member.name} row {line_number} is not a three-field feature"
                    )
                identifiers.append(fields[0])
                symbols.append(fields[1])
                feature_types.append(fields[2])
    except (OSError, EOFError, UnicodeError) as error:
        raise MetadataError(
            f"cannot read archive feature member: {member.name}"
        ) from error
    duplicate_identifiers = _duplicates(identifiers)
    if duplicate_identifiers:
        raise MetadataError(
            f"{member.name} has duplicate feature identifiers: "
            f"{duplicate_identifiers!r}"
        )
    return FeatureAxis(
        identifiers=tuple(identifiers),
        symbols=tuple(symbols),
        feature_types=tuple(feature_types),
    )


def read_archive_feature_axes(
    path: Path,
    *,
    channel: int = 1,
) -> ArchiveFeatureAxes:
    """Read only the three compressed feature tables for one archive channel."""

    if not isinstance(channel, int) or isinstance(channel, bool) or channel <= 0:
        raise MetadataError("channel must be a positive integer")
    if not path.is_file() or path.is_symlink():
        raise MetadataError(f"missing or unsafe raw archive: {path}")

    try:
        with tarfile.open(path, mode="r:") as archive:
            matches: dict[str, tuple[str, tarfile.TarInfo]] = {}
            for member in archive.getmembers():
                parsed = _FEATURE_MEMBER.fullmatch(member.name)
                if parsed is None or int(parsed.group("channel")) != channel:
                    continue
                modality = parsed.group("modality")
                if modality in matches:
                    raise MetadataError(
                        f"duplicate {modality} feature member for channel {channel}"
                    )
                matches[modality] = (parsed.group("accession"), member)
            if set(matches) != {"transcriptome", "guides", "labels"}:
                raise MetadataError(
                    f"channel {channel} does not have exactly three feature tables"
                )
            accessions = {value[0] for value in matches.values()}
            if len(accessions) != 1:
                raise MetadataError(
                    f"channel {channel} feature members disagree on accession"
                )
            axes = {
                modality: _read_feature_member(archive, member)
                for modality, (_, member) in matches.items()
            }
    except (OSError, tarfile.TarError) as error:
        raise MetadataError(f"cannot read raw archive {path}") from error

    return ArchiveFeatureAxes(
        channel=channel,
        accession=accessions.pop(),
        guides=axes["guides"],
        labels=axes["labels"],
        transcriptome=axes["transcriptome"],
    )


def _read_guide_label(guide_id: str) -> str:
    match = _GUIDE_ID.fullmatch(guide_id)
    if match is None:
        raise MetadataError(f"guide ID has no terminal positive ordinal: {guide_id!r}")
    return match.group("label")


def build_guide_crosswalk(
    guide_ids: Sequence[str],
    target_records: Sequence[TargetRecord],
    panel_targets: Sequence[str],
    *,
    archive_guide_ids: Sequence[str] | None = None,
) -> GuideCrosswalk:
    """Build the canonical guide-to-target mapping.

    Composite labels, labels absent from Table S1, and symbols associated with
    multiple Entrez IDs receive ``NO_UNIQUE_TARGET``.
    """

    canonical_guides = tuple(guide_ids)
    if not canonical_guides or any(
        not isinstance(value, str) for value in canonical_guides
    ):
        raise MetadataError("guide IDs must be nonempty strings")
    duplicates = _duplicates(canonical_guides)
    if duplicates:
        raise MetadataError(f"guide IDs are not unique: {duplicates!r}")
    if archive_guide_ids is not None and tuple(archive_guide_ids) != canonical_guides:
        raise MetadataError("FASTA and archive guide feature axes differ")

    source_rows: defaultdict[str, list[TargetRecord]] = defaultdict(list)
    target_symbols: list[str] = []
    for record in target_records:
        if not isinstance(record, TargetRecord):
            raise MetadataError("target_records must contain TargetRecord values")
        if record.symbol not in source_rows:
            target_symbols.append(record.symbol)
        source_rows[record.symbol].append(record)
    if not target_symbols:
        raise MetadataError("source target table is empty")

    target_index = {symbol: index for index, symbol in enumerate(target_symbols)}
    target_entrez_ids = tuple(
        tuple(sorted({record.entrez_id for record in source_rows[symbol]}))
        for symbol in target_symbols
    )
    duplicate_target_symbols = tuple(
        DuplicateSymbol(
            symbol=symbol,
            entrez_ids=tuple(record.entrez_id for record in source_rows[symbol]),
            workbook_rows=tuple(record.workbook_row for record in source_rows[symbol]),
        )
        for symbol in target_symbols
        if len({record.entrez_id for record in source_rows[symbol]}) != 1
    )

    labels: list[str] = []
    statuses: list[str] = []
    mapping = np.full(len(canonical_guides), NO_UNIQUE_TARGET, dtype=np.int32)
    non_targeting = np.zeros(len(canonical_guides), dtype=np.bool_)
    composite_labels: set[str] = set()
    unmapped_labels: set[str] = set()
    ambiguous_labels: set[str] = set()

    for guide_index, guide_id in enumerate(canonical_guides):
        label = _read_guide_label(guide_id)
        labels.append(label)
        if label == "non-targeting":
            non_targeting[guide_index] = True
            statuses.append(GUIDE_NON_TARGETING)
        elif "~" in label:
            composite_labels.add(label)
            statuses.append(GUIDE_COMPOSITE)
        elif label not in source_rows:
            unmapped_labels.add(label)
            statuses.append(GUIDE_UNMAPPED)
        elif len({record.entrez_id for record in source_rows[label]}) != 1:
            ambiguous_labels.add(label)
            statuses.append(GUIDE_AMBIGUOUS)
        else:
            mapping[guide_index] = target_index[label]
            statuses.append(GUIDE_SUPPORTED)

    canonical_panel = tuple(panel_targets)
    if any(not isinstance(value, str) or not value for value in canonical_panel):
        raise MetadataError("panel targets must be nonempty strings")
    panel_duplicates = _duplicates(canonical_panel)
    if panel_duplicates:
        raise MetadataError(f"panel targets are not unique: {panel_duplicates!r}")
    supported_symbols = {
        target_symbols[index] for index in mapping if index != NO_UNIQUE_TARGET
    }
    covered = tuple(symbol for symbol in canonical_panel if symbol in supported_symbols)
    missing = tuple(
        symbol for symbol in canonical_panel if symbol not in supported_symbols
    )

    mapping.flags.writeable = False
    non_targeting.flags.writeable = False
    return GuideCrosswalk(
        guide_ids=canonical_guides,
        guide_labels=tuple(labels),
        guide_status=tuple(statuses),
        target_symbols=tuple(target_symbols),
        target_entrez_ids=target_entrez_ids,
        guide_to_target=mapping,
        non_targeting_guides=non_targeting,
        duplicate_target_symbols=duplicate_target_symbols,
        composite_labels=tuple(sorted(composite_labels)),
        unmapped_labels=tuple(sorted(unmapped_labels)),
        ambiguous_labels=tuple(sorted(ambiguous_labels)),
        panel_targets=canonical_panel,
        panel_covered_targets=covered,
        panel_missing_targets=missing,
    )


def build_hash_crosswalk(
    hash_labels: Sequence[str],
    hash_oligos: Sequence[HashOligo],
) -> HashCrosswalk:
    """Map archive labels to conditions and validate Table S6 counts."""

    labels = tuple(hash_labels)
    if not labels or any(not isinstance(label, str) or not label for label in labels):
        raise MetadataError("hash labels must be nonempty strings")
    duplicates = _duplicates(labels)
    if duplicates:
        raise MetadataError(f"hash labels are not unique: {duplicates!r}")

    records: list[HashLabelRecord] = []
    label_condition_counts: Counter[str] = Counter()
    condition_ordinals: defaultdict[str, list[int]] = defaultdict(list)
    for feature_index, label in enumerate(labels):
        match = _HASH_LABEL.fullmatch(label)
        if match is None:
            raise MetadataError(f"unsupported hash label: {label!r}")
        condition = match.group("condition")
        number = match.group("number")
        ordinal = 1 if condition == "untreated" and number == "" else int(number or 0)
        if ordinal <= 0:
            raise MetadataError(f"hash label has no positive ordinal: {label!r}")
        condition_ordinals[condition].append(ordinal)
        label_condition_counts[condition] += 1
        records.append(HashLabelRecord(feature_index, label, condition))

    for condition, ordinals in condition_ordinals.items():
        if sorted(ordinals) != list(range(1, len(ordinals) + 1)):
            raise MetadataError(
                f"{condition} hash labels do not form a contiguous one-based axis"
            )

    workbook_counts: Counter[str] = Counter()
    workbook_oligos = tuple(hash_oligos)
    seen_oligos: set[str] = set()
    for record in workbook_oligos:
        if not isinstance(record, HashOligo):
            raise MetadataError("hash_oligos must contain HashOligo values")
        if record.oligo in seen_oligos:
            raise MetadataError(f"duplicate workbook hash oligo: {record.oligo}")
        seen_oligos.add(record.oligo)
        workbook_counts[record.condition] += 1
    if label_condition_counts != workbook_counts:
        raise MetadataError(
            "archive hash-label condition counts disagree with supplementary Table S6"
        )
    return HashCrosswalk(records=tuple(records), workbook_oligos=workbook_oligos)


def build_response_crosswalk(
    response_records: Sequence[ResponseRecord],
    source_feature_symbols: Sequence[str],
    vcc_genes: Sequence[str],
    panel_targets: Sequence[str],
) -> ResponseCrosswalk:
    """Map Table S2 response symbols to source feature rows and VCC indices."""

    responses = tuple(response_records)
    if not responses:
        raise MetadataError("response table is empty")
    for record in responses:
        if not isinstance(record, ResponseRecord):
            raise MetadataError("response_records must contain ResponseRecord values")

    response_symbols = tuple(record.symbol for record in responses)
    response_entrez = tuple(record.entrez_id for record in responses)
    duplicate_response_symbols = tuple(
        str(value) for value in _duplicates(response_symbols)
    )
    duplicate_response_entrez_ids = tuple(
        int(value) for value in _duplicates(response_entrez)
    )

    source_lookup: defaultdict[str, list[int]] = defaultdict(list)
    for index, symbol in enumerate(source_feature_symbols):
        if not isinstance(symbol, str) or not symbol:
            raise MetadataError("source feature symbols must be nonempty strings")
        source_lookup[symbol].append(index)

    vcc_lookup: defaultdict[str, list[int]] = defaultdict(list)
    for index, symbol in enumerate(vcc_genes):
        if not isinstance(symbol, str) or not symbol:
            raise MetadataError("VCC genes must be nonempty strings")
        vcc_lookup[symbol].append(index)

    panel = tuple(panel_targets)
    if _duplicates(panel):
        raise MetadataError("panel targets must be unique")
    panel_set = set(panel)

    records: list[ResponseGeneCrosswalk] = []
    for response_index, record in enumerate(responses):
        records.append(
            ResponseGeneCrosswalk(
                response_index=response_index,
                symbol=record.symbol,
                entrez_id=record.entrez_id,
                group=record.group,
                source_feature_indices=tuple(source_lookup.get(record.symbol, ())),
                vcc_gene_indices=tuple(vcc_lookup.get(record.symbol, ())),
                is_panel_target=record.symbol in panel_set,
            )
        )

    missing_source = tuple(
        record.symbol for record in records if not record.source_feature_indices
    )
    ambiguous_source = tuple(
        record.symbol for record in records if len(record.source_feature_indices) > 1
    )
    missing_vcc = tuple(
        record.symbol for record in records if not record.vcc_gene_indices
    )
    ambiguous_vcc = tuple(
        record.symbol for record in records if len(record.vcc_gene_indices) > 1
    )
    visible_vcc_indices = tuple(
        record.vcc_gene_indices[0]
        for record in records
        if len(record.vcc_gene_indices) == 1 and not record.is_panel_target
    )
    return ResponseCrosswalk(
        records=tuple(records),
        duplicate_response_symbols=duplicate_response_symbols,
        duplicate_response_entrez_ids=duplicate_response_entrez_ids,
        missing_source_symbols=missing_source,
        ambiguous_source_symbols=ambiguous_source,
        missing_vcc_symbols=missing_vcc,
        ambiguous_vcc_symbols=ambiguous_vcc,
        visible_vcc_indices=visible_vcc_indices,
    )


def validate_pinned_crosswalks(crosswalks: MetadataCrosswalks) -> None:
    """Hard-fail if any preregistered source metadata count has drifted."""

    guide = crosswalks.guide
    hashes = crosswalks.hashes
    response = crosswalks.response
    axes = crosswalks.archive_axes

    if len(guide.guide_ids) != EXPECTED_GUIDES:
        raise MetadataError("guide axis is not exactly 83,401 records")
    if len(guide.target_symbols) != EXPECTED_SOURCE_TARGET_SYMBOLS:
        raise MetadataError("source target axis is not exactly 18,465 symbols")
    if np.count_nonzero(guide.non_targeting_guides) != EXPECTED_NON_TARGETING:
        raise MetadataError("guide axis does not contain exactly 250 NTC records")
    if guide.status_counts != {
        GUIDE_AMBIGUOUS: 5,
        GUIDE_COMPOSITE: 201,
        GUIDE_NON_TARGETING: 250,
        GUIDE_SUPPORTED: 82_945,
    }:
        raise MetadataError(f"guide classification drifted: {guide.status_counts}")
    if guide.unmapped_labels:
        raise MetadataError(
            f"unexpected unmapped guide labels: {guide.unmapped_labels}"
        )
    if guide.ambiguous_labels != ("MEMO1",):
        raise MetadataError(
            f"unexpected ambiguous guide labels: {guide.ambiguous_labels}"
        )
    if guide.panel_missing_targets != ("ZNF714",):
        raise MetadataError(
            f"panel coverage drifted; missing {guide.panel_missing_targets}"
        )
    if len(guide.panel_covered_targets) != 299:
        raise MetadataError("guide library does not cover exactly 299 panel targets")

    if len(hashes.records) != EXPECTED_HASH_LABELS:
        raise MetadataError("hash-label axis is not exactly seven records")
    if Counter(record.condition for record in hashes.records) != {
        "activated": 6,
        "untreated": 1,
    }:
        raise MetadataError("hash-label condition mapping drifted")

    if len(response.records) != EXPECTED_RESPONSE_GENES:
        raise MetadataError("response axis is not exactly 374 genes")
    if response.duplicate_response_symbols or response.duplicate_response_entrez_ids:
        raise MetadataError("response rows are not unique by symbol and Entrez ID")
    if (
        len(
            [record for record in response.records if record.vcc_gene_index is not None]
        )
        != 358
    ):
        raise MetadataError("response-to-VCC coverage is not exactly 358 genes")
    if len(response.visible_vcc_indices) != 349:
        raise MetadataError("PDS-visible response axis is not exactly 349 genes")
    if len(set(response.visible_vcc_indices)) != 349:
        raise MetadataError("PDS-visible VCC response indices are not unique")

    if len(axes.guides.identifiers) != EXPECTED_GUIDES:
        raise MetadataError("archive guide feature count drifted")
    if len(axes.labels.identifiers) != EXPECTED_HASH_LABELS:
        raise MetadataError("archive label feature count drifted")
    if len(axes.transcriptome.identifiers) != 20_606:
        raise MetadataError("archive transcriptome feature count drifted")


def load_pinned_crosswalks(
    *,
    guide_fasta: Path = GUIDE_FASTA,
    workbook_path: Path = WORKBOOK,
    raw_archive: Path = RAW_ARCHIVE,
    challenge_targets: Path = CHALLENGE_TARGETS,
    challenge_genes: Path = CHALLENGE_GENES,
    channel: int = 1,
) -> MetadataCrosswalks:
    """Load all pinned metadata and enforce the frozen experiment contract."""

    for path in (guide_fasta, workbook_path, challenge_targets, challenge_genes):
        verify_pinned_file(path)
    workbook = read_workbook_metadata(workbook_path)
    if len(workbook.targets) != EXPECTED_SOURCE_TARGET_ROWS:
        raise MetadataError("Table S1 does not contain exactly 18,466 rows")
    if len(workbook.responses) != EXPECTED_RESPONSE_GENES:
        raise MetadataError("Table S2 does not contain exactly 374 integer-Entrez rows")
    if len(workbook.hash_oligos) != EXPECTED_HASH_LABELS:
        raise MetadataError("Table S6 does not contain exactly seven hash oligos")

    axes = read_archive_feature_axes(raw_archive, channel=channel)
    guide_ids = read_fasta_ids(guide_fasta)
    panel = read_single_column_csv(challenge_targets, "target_gene")
    genes = read_single_column_csv(challenge_genes, "gene_name")
    if len(panel) != EXPECTED_PANEL_TARGETS:
        raise MetadataError("challenge panel does not contain exactly 300 targets")

    guide = build_guide_crosswalk(
        guide_ids,
        workbook.targets,
        panel,
        archive_guide_ids=axes.guides.identifiers,
    )
    hashes = build_hash_crosswalk(axes.labels.identifiers, workbook.hash_oligos)
    response = build_response_crosswalk(
        workbook.responses,
        axes.transcriptome.symbols,
        genes,
        panel,
    )
    crosswalks = MetadataCrosswalks(
        archive_axes=axes,
        guide=guide,
        hashes=hashes,
        response=response,
    )
    validate_pinned_crosswalks(crosswalks)
    return crosswalks
