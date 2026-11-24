"""Check released description/embedding identity; never repairs an embedding."""

import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile
import numpy as np
import acquire333 as A
import convert_recovery333 as R
import convert_canonical333 as C

MANIFEST = A.HERE / "identity-inputs.json"
FILES = [Path(__file__).resolve(), A.HERE / "IDENTITY-AUDIT-SCOPE.md",
         A.HERE / "canonical-inputs.json", A.WORK / "canonical-results.json"]
TEXT_MEMBERS = ("NCBI_summary_of_genes.json", "NCBI_UniProt_summary_of_genes.json")


def freeze():
    C.verify()
    if MANIFEST.exists():
        raise FileExistsError("identity audit already frozen")
    A.write(MANIFEST, {"sha256": {str(p.relative_to(A.ROOT)): A.digest(p) for p in FILES},
                       "description_members": TEXT_MEMBERS})


def verify():
    C.verify()
    manifest = json.loads(MANIFEST.read_text())
    for name, expected in manifest["sha256"].items():
        assert A.digest(A.ROOT / name) == expected, name
    assert set(manifest["sha256"]) == {str(p.relative_to(A.ROOT)) for p in FILES}
    own = [str(p.relative_to(A.ROOT)) for p in FILES] + [str(MANIFEST.relative_to(A.ROOT))]
    for name in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", name], cwd=A.ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=A.ROOT)


def audit():
    verify()
    output = A.WORK / "identity-audit.json"
    if output.exists():
        raise FileExistsError("identity audit exists")
    results = json.loads((A.WORK / "canonical-results.json").read_text())
    npz = A.WORK / "canonical-embeddings.npz"
    assert A.digest(npz) == results["embedding_sha256"]
    with np.load(npz, allow_pickle=False) as z:
        genes, vectors = z["genes"], z["vectors"]
    groups = defaultdict(list)
    for gene, row in zip(genes, vectors):
        groups[hashlib.sha256(row.tobytes()).hexdigest()].append(str(gene))
    duplicates = [g for g in groups.values() if len(g) > 1]
    with (A.ROOT / "data/challenge/pert_counts.csv").open() as f:
        targets = {r["target_gene"] for r in csv.DictReader(f)}
    text_data = {}
    with zipfile.ZipFile(A.ZIP) as archive:
        for name in TEXT_MEMBERS:
            matches = [i for i in archive.infolist() if Path(i.filename).name == name]
            assert len(matches) == 1 and matches[0].file_size < 30_000_000
            raw = archive.read(matches[0])
            text = json.loads(raw)
            assert isinstance(text, dict) and all(isinstance(k, str) and isinstance(v, str)
                                                  for k, v in text.items())
            text_data[name] = (text, hashlib.sha256(raw).hexdigest())
        A.MEMBER = R.ACTUAL_MEMBER
        member = A.select_member(archive.infolist())
        with archive.open(member) as stream:
            original = A.RestrictedNumpyUnpickler(stream).load()
    comparisons = []
    for group in duplicates:
        record = {"genes": group, "panel_targets": sorted(targets & set(group)), "text": {}}
        for name, (text, _) in text_data.items():
            hashes = {g: hashlib.sha256(text[g].encode()).hexdigest()
                      for g in group if g in text and text[g].strip()}
            record["text"][name] = {"nonempty_primary_texts": len(hashes),
                                     "distinct_text_hashes": len(set(hashes.values())),
                                     "hashes": hashes}
        if record["panel_targets"]:
            arrays = [np.asarray(original[g]) for g in group]
            record["raw_vectors_equal"] = all(np.array_equal(arrays[0], v) for v in arrays[1:])
            record["same_python_object"] = all(original[group[0]] is original[g] for g in group[1:])
            record["raw_dtypes"] = [str(v.dtype) for v in arrays]
        comparisons.append(record)
    flagged = sorted(targets & {g for group in duplicates for g in group})
    summaries = {}
    for name, (text, text_sha) in text_data.items():
        differing = [r for r in comparisons if r["text"][name]["distinct_text_hashes"] > 1]
        summaries[name] = {
            "member_sha256": text_sha, "primary_keys": len(text),
            "canonical_nonempty_texts": sum(g in text and bool(text[g].strip()) for g in genes),
            "panel_nonempty_texts": sum(g in text and bool(text[g].strip()) for g in targets),
            "duplicate_vector_groups_with_distinct_texts": len(differing),
            "panel_targets_in_those_groups": sorted({g for r in differing for g in r["panel_targets"]})}
    payload = {"status": "QUARANTINED_IDENTITY_UNRESOLVED",
               "duplicate_vector_groups": len(duplicates),
               "canonical_genes_in_duplicate_groups": sum(map(len, duplicates)),
               "panel_targets_in_duplicate_groups": flagged,
               "source_texts": summaries,
               "target_groups_checked_raw": sum(bool(r["panel_targets"]) for r in comparisons),
               "target_groups_raw_vectors_equal": sum(bool(r.get("raw_vectors_equal")) for r in comparisons),
               "comparisons": comparisons, "inputs_sha256": A.digest(MANIFEST),
               "corrected_vectors_written": False, "biological_outcomes_read": False,
               "limitation": "Distinct released texts sharing exact vectors raise an identity concern; "
                             "this does not identify which vector is wrong or repair any alias."}
    A.write(output, payload)
    print(json.dumps({k: v for k, v in payload.items() if k != "comparisons"}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "audit", "verify"))
    globals()[parser.parse_args().mode]()
