"""Bounded official metadata retrieval. Does not request expression data."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
PAIRS = ROOT / "experiments/324-residual-label-feasibility/work/panel__id_fallback__same_method__pairs.csv"
MAX_BYTES = 30 << 20


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def fields_soft(text, expected, prefix="Series"):
    if f"^{prefix.upper()} = {expected}" not in text:
        raise ValueError(f"requested {expected} missing from returned SOFT header")
    fields = {}
    for line in text.splitlines():
        if line.startswith(f"!{prefix}_") and " = " in line:
            key, value = line.split(" = ", 1)
            fields.setdefault(key[1:], []).append(value)
    if fields.get(f"{prefix}_geo_accession") != [expected]:
        raise ValueError("returned accession does not match request")
    return fields


def fields_encode(text, expected):
    value = json.loads(text)
    if value.get("accession") != expected:
        raise ValueError("ENCODE accession mismatch")
    return value


def source_url(accession):
    if re.fullmatch(r"GSE\d+", accession):
        return ("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"
                f"acc={accession}&targ=self&form=text&view=full")
    if re.fullmatch(r"ENCSR[A-Z0-9]{6}", accession):
        return f"https://www.encodeproject.org/experiments/{accession}/?format=json"
    raise ValueError("unregistered accession type")


def publication_and_samples(summary):
    """Publication grouping and sample intersection are separate facts."""
    if summary["kind"] == "GEO":
        f = summary["fields"]
        return (set(f.get("Series_pubmed_id", [])), set(f.get("Series_sample_id", [])))
    f = summary["fields"]
    refs = f.get("references", [])
    publications = set()
    for ref in refs:
        if isinstance(ref, str):
            publications.add(ref)
        elif isinstance(ref, dict):
            publications.add(ref.get("@id") or ref.get("uuid") or ref.get("title", "UNKNOWN"))
    # ENCODE replicate objects are not GEO samples; inspect them separately.
    return publications, set()


def freeze():
    out = HERE / "inputs.json"
    if out.exists():
        raise FileExistsError(out)
    files = [PAIRS] + [HERE / n for n in
                       ("PREREGISTRATION.md", "fetch325.py", "test325.py")]
    write_json(out, {"inputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in files}})


def verify():
    manifest = json.loads((HERE / "inputs.json").read_text())
    for path, expected in manifest["inputs_sha256"].items():
        if sha(ROOT / path) != expected:
            raise ValueError(f"input drift: {path}")
    local_files = [str(p.relative_to(ROOT)) for p in HERE.glob("*.py")]
    local_files += [str((HERE / n).relative_to(ROOT)) for n in
                    ("inputs.json", "PREREGISTRATION.md")]
    for path in local_files:
        subprocess.run(["git", "ls-files", "--error-unmatch", path], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    if subprocess.check_output(["git", "status", "--porcelain", "--", *local_files],
                               cwd=ROOT, text=True):
        raise ValueError("Commit frozen retrieval code before fetching")


def fetch():
    verify()
    with PAIRS.open(newline="") as handle:
        pairs = list(csv.DictReader(handle))
    if len(pairs) != 10:
        raise AssertionError("Expected exactly 324's ten pairs")
    accessions = sorted({r[k] for r in pairs for k in ("accession_a", "accession_b")})
    rawdir = WORK / "sources"
    rawdir.mkdir(parents=True, exist_ok=True)
    if (WORK / "series-summary.json").exists():
        raise FileExistsError("Source retrieval already complete")
    records = {}
    started = time.monotonic()
    for accession in accessions:
        suffix = ".soft" if accession.startswith("GSE") else ".json"
        path = rawdir / f"{accession}{suffix}"
        receipt = rawdir / f"{accession}.receipt.json"
        if time.monotonic()-started > 2400:
            raise TimeoutError("Leave time for the report inside the 45-minute bound")
        url = source_url(accession)
        if receipt.exists():
            rec = json.loads(receipt.read_text())
            if rec.get("status") != "OK":
                records[accession] = rec
                continue
            if sha(path) != rec["sha256"]:
                raise ValueError(f"cached source changed: {accession}")
        else:
            process = subprocess.run(
                ["curl", "--fail", "--location", "--silent", "--show-error",
                 "--max-time", "25", "--max-filesize", str(2 << 20), url],
                capture_output=True, timeout=30)
            rec = {"accession": accession, "url": url,
                   "queried_at_utc": datetime.now(timezone.utc).isoformat()}
            if process.returncode:
                rec.update(status="HTTP_ERROR", error=process.stderr.decode(errors="replace"))
                write_json(receipt, rec)
                records[accession] = rec
                print(f"{accession}: HTTP_ERROR", flush=True)
                continue
            body = process.stdout
            if len(body) > (2 << 20) or sum(p.stat().st_size for p in rawdir.iterdir())+len(body) > MAX_BYTES:
                raise ValueError("metadata byte budget exceeded")
            path.write_bytes(body)
            rec.update(status="OK", bytes=len(body), sha256=sha(path))
            write_json(receipt, rec)
        text = path.read_text()
        kind = "GEO" if accession.startswith("GSE") else "ENCODE"
        fields = fields_soft(text, accession) if kind == "GEO" else fields_encode(text, accession)
        records[accession] = {**rec, "kind": kind, "fields": fields}
        print(f"{accession}: metadata ID verified", flush=True)

    comparisons = []
    for row in pairs:
        a, b = records[row["accession_a"]], records[row["accession_b"]]
        item = dict(row)
        if a.get("status") != "OK" or b.get("status") != "OK":
            item.update(status="UNRESOLVED_SOURCE_UNAVAILABLE")
        else:
            pubs_a, samples_a = publication_and_samples(a)
            pubs_b, samples_b = publication_and_samples(b)
            item.update(shared_publication_identifiers=sorted(pubs_a & pubs_b),
                        shared_geo_sample_ids=sorted(samples_a & samples_b),
                        publication_identifier_types_match=a["kind"] == b["kind"],
                        status="METADATA_ONLY_REQUIRES_UNIT_INTERPRETATION")
        comparisons.append(item)
    write_json(WORK / "series-summary.json", {
        "input_manifest_sha256": sha(HERE / "inputs.json"), "records": records,
        "pairs": comparisons, "seconds": time.monotonic()-started,
        "response_data_read": False, "reserve_access": False})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "fetch"))
    args = parser.parse_args()
    globals()[args.mode]()


if __name__ == "__main__":
    main()
