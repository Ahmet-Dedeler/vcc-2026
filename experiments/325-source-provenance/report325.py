"""Summarize retrieved metadata without requesting or measuring response data."""

import csv
import json

import fetch325 as F


def summarize(record):
    if record["kind"] == "GEO":
        f = record["fields"]
        return {
            "accession": record["accession"],
            "kind": "GEO",
            "title": f.get("Series_title", []),
            "design": f.get("Series_overall_design", []),
            "publications": ["PMID:" + p for p in f.get("Series_pubmed_id", [])],
            "samples": f.get("Series_sample_id", []),
            "relations": f.get("Series_relation", []),
        }
    f = record["fields"]
    biosamples = [r["library"]["biosample"] for r in f["replicates"]]
    return {
        "accession": record["accession"],
        "kind": "ENCODE",
        "description": f.get("description"),
        "publications": sorted({
            identifier for ref in f.get("references", []) if isinstance(ref, dict)
            for identifier in ref.get("identifiers", [])
        }),
        "lab": f["lab"]["name"],
        "biosamples": sorted({b["accession"] for b in biosamples}),
        "modifications": sorted({
            m["accession"] for b in biosamples for m in b.get("genetic_modifications", [])
        }),
        "reagents": sorted({
            r["identifier"] for b in biosamples
            for m in b.get("genetic_modifications", []) for r in m.get("reagents", [])
        }),
        "controls": sorted({c["accession"] for c in f.get("possible_controls", [])}),
    }


JUDGMENTS = {
    ("EZH2", "CVCL_0004"): (
        "UNRESOLVED_GROUP_MODALITY",
        "Different PMIDs. GSE162623 explicitly describes EZH2 KO; GSE193486 "
        "describes silencer-region KOs and drug treatments. The exact GPSAdb "
        "group-to-sample mapping is needed; a series title cannot adjudicate it."),
    ("EZH2", "CVCL_Y019"): (
        "UNRESOLVED_GROUP_MODALITY",
        "Different PMIDs. GSE133183 includes EZH2KO, while GSE285461 includes "
        "mixed CRISPR regulatory-element experiments. Exact groups and constructs unresolved."),
    ("EZH2", "CVCL_0031"): (
        "UNRESOLVED_MATERIAL_INDEPENDENCE",
        "Different PMIDs and shEZH2 designs; exact samples, controls, constructs "
        "and material independence not verified."),
    ("EZH2", "NAME:NHDF"): (
        "UNRESOLVED_CONTEXT_AND_MATERIAL",
        "Both series describe shEZH2, but one samples days 4/8 and the other "
        "describes stable fibroblast lines. NHDF is name-only; donor and construct "
        "identity and a publication for GSE98768 remain unresolved."),
    ("KDM5A", "CVCL_0063"): (
        "UNRESOLVED_GROUP_AND_CONSTRUCT",
        "Different parent series. GSE277160 describes many ZFP knockdowns; the "
        "specific KDM5A samples, controls, constructs and independence remain unresolved."),
    ("KHSRP", "CVCL_0004"): (
        "EXPLICIT_ENCODE_MATERIAL_REUSE",
        "GSE145984 explicitly says its KD/control samples were originally generated "
        "for ENCODE by the Graveley lab. Not an independent biological study replicate; "
        "exact GSM-to-ENCBS matching has not been verified."),
    ("KHSRP", "CVCL_0027"): (
        "EXPLICIT_ENCODE_MATERIAL_REUSE",
        "Same explicit reuse statement for HepG2. Publication-level crossing exists "
        "as a grouping, but does not create independent biological perturbation material."),
    ("MTA1", "CVCL_0291"): (
        "SHARED_STUDY_DIFFERENT_STATE",
        "Both SubSeries of GSE123959. One design is nonsynchronized and the other "
        "synchronized HCT116; not interchangeable same-state independent-study repeats."),
    ("SRSF5", "CVCL_0027"): (
        "DISTINCT_CONSTRUCT_SHARED_STUDY",
        "Same ENCODE lab/shared publications, but distinct TRC constructs, biosample "
        "accessions and control-experiment accessions. A within-study construct contrast, "
        "not independent publication replication; underlying material independence unresolved."),
    ("TAF15", "CVCL_0004"): (
        "SAME_CONSTRUCT_SHARED_STUDY",
        "Same ENCODE lab/shared publications and the same TRC construct, despite "
        "different modification IDs, biosample accessions and control experiments."),
}


def main():
    out = F.WORK / "evidence.json"
    if out.exists():
        raise FileExistsError(out)
    manifest = json.loads((F.HERE / "inputs.json").read_text())
    for path, expected in manifest["inputs_sha256"].items():
        if F.sha(F.ROOT / path) != expected:
            raise ValueError(f"frozen retrieval input drift: {path}")
    data = json.loads((F.WORK / "series-summary.json").read_text())
    records = {}
    receipts = {}
    for accession, rec in data["records"].items():
        if rec["status"] != "OK":
            raise ValueError(f"unavailable source: {accession}")
        suffix = ".soft" if rec["kind"] == "GEO" else ".json"
        raw = F.WORK / "sources" / (accession + suffix)
        if F.sha(raw) != rec["sha256"] or raw.stat().st_size != rec["bytes"]:
            raise ValueError(f"source receipt mismatch: {accession}")
        if rec["kind"] == "GEO":
            parsed = F.fields_soft(raw.read_text(), accession)
        else:
            parsed = F.fields_encode(raw.read_text(), accession)
        if parsed != rec["fields"]:
            raise ValueError(f"parsed source mismatch: {accession}")
        records[accession] = summarize(rec)
        receipts[accession] = {k: rec[k] for k in
                              ("url", "queried_at_utc", "sha256", "bytes")}

    khsrp_design = " ".join(records["GSE145984"]["design"]).lower()
    assert "originally generated for the encode project" in khsrp_design
    assert "graveley" in khsrp_design
    assert "SubSeries of: GSE123959" in records["GSE123957"]["relations"]
    assert "SubSeries of: GSE123959" in records["GSE123958"]["relations"]
    assert records["ENCSR447UCG"]["reagents"] == ["trc:TRCN0000000141"]
    assert records["ENCSR781YNI"]["reagents"] == ["trc:TRCN0000284799"]
    assert records["ENCSR031RZS"]["reagents"] == records["ENCSR611ZAL"]["reagents"]
    assert records["ENCSR031RZS"]["reagents"] == ["trc:TRCN0000020143"]

    pairs = []
    for p in data["pairs"]:
        status, conclusion = JUDGMENTS[p["target"], p["context"]]
        a, b = records[p["accession_a"]], records[p["accession_b"]]
        pairs.append({
            **{k: p[k] for k in ("target", "context", "method", "accession_a", "accession_b")},
            "status": status, "conclusion": conclusion,
            "shared_publications": ";".join(sorted(set(a["publications"]) & set(b["publications"]))),
            "fully_independent_same_state_repeat_verified": False,
        })
    assert len(pairs) == len(JUDGMENTS) == 10
    result = {
        "status": "COMPLETE_METADATA_ONLY_NO_LABEL_PILOT",
        "sources_verified": len(records), "pairs": pairs,
        "records": records, "receipts": receipts,
        "source_bytes": sum(r["bytes"] for r in receipts.values()),
        "followup_metadata_objects": 0,
        "expression_downloaded": False, "reserve_access": False,
        "limitation": "Unresolved is not evidence of nonindependence or incorrect GPSAdb annotation.",
        "series_summary_sha256": F.sha(F.WORK / "series-summary.json"),
        "input_manifest_sha256": F.sha(F.HERE / "inputs.json"),
    }
    F.write_json(out, result)
    with (F.WORK / "pairs.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pairs[0]))
        writer.writeheader()
        writer.writerows(pairs)
    print(json.dumps({k: result[k] for k in
                      ("status", "sources_verified", "source_bytes")}, indent=2))


if __name__ == "__main__":
    main()
