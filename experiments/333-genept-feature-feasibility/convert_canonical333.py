"""Exact-symbol-only feature conversion; never normalize an alias key."""

import argparse
import csv
import json
from pathlib import Path
import resource
import subprocess
import time
import zipfile
import numpy as np
import acquire333 as A
import convert_recovery333 as R

MANIFEST = A.HERE / "canonical-inputs.json"
FILES = [Path(__file__).resolve(), A.HERE / "test_canonical333.py",
         A.HERE / "RECOVERY-002.md", A.WORK / "name-audit.json",
         A.HERE / "recovery-inputs.json", A.HERE / "NAME-GATE-FAILURE.md"]


def freeze():
    R.verify()
    if MANIFEST.exists():
        raise FileExistsError("canonical recovery frozen")
    A.write(MANIFEST, {"sha256": {str(p.relative_to(A.ROOT)): A.digest(p) for p in FILES},
                       "identifier_policy": "canonical-exact-case-sensitive-only"})


def verify():
    R.verify()
    m = json.loads(MANIFEST.read_text())
    assert m["identifier_policy"] == "canonical-exact-case-sensitive-only"
    assert set(m["sha256"]) == {str(p.relative_to(A.ROOT)) for p in FILES}
    for name, expected in m["sha256"].items():
        assert A.digest(A.ROOT / name) == expected, name
    own = [str(p.relative_to(A.ROOT)) for p in FILES] + [str(MANIFEST.relative_to(A.ROOT))]
    for name in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", name], cwd=A.ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=A.ROOT)


def exact_subset(mapping, canonical):
    if not isinstance(mapping, dict):
        raise ValueError("dictionary required")
    if len(set(canonical)) != len(canonical):
        raise ValueError("canonical identifiers not unique")
    if any(not isinstance(g, str) or not g or g != g.strip() for g in canonical):
        raise ValueError("invalid canonical symbol")
    names = [g for g in canonical if g in mapping]
    if not names:
        raise ValueError("no exact symbols")
    vectors, dim = [], None
    for gene in names:
        value = np.asarray(mapping[gene])
        if value.ndim != 1 or not len(value) or value.dtype.kind not in "fi":
            raise ValueError("invalid selected vector: "+gene)
        if not np.isfinite(value).all() or not np.any(value):
            raise ValueError("nonfinite/zero selected vector: "+gene)
        if dim is None:
            dim = len(value)
        if len(value) != dim:
            raise ValueError("inconsistent selected dimensions")
        vectors.append(value.astype(np.float32))
    matrix = np.stack(vectors)
    if not np.isfinite(matrix).all():
        raise ArithmeticError("float32 overflow")
    return np.array(names), matrix


def convert():
    verify()
    output = A.WORK / "canonical-embeddings.npz"
    result = A.WORK / "canonical-results.json"
    if output.exists() or result.exists():
        raise FileExistsError("canonical conversion exists")
    tick = time.monotonic()
    download = json.loads((A.WORK / "download.json").read_text())
    assert A.digest(A.ZIP) == download["sha256"]
    with (A.ROOT / "data/challenge/gene_names.csv").open() as f:
        genes = [r["gene_name"] for r in csv.DictReader(f)]
    with (A.ROOT / "data/challenge/pert_counts.csv").open() as f:
        targets = [r["target_gene"] for r in csv.DictReader(f)]
    assert len(genes) == 18533 and len(targets) == 300
    canonical = list(dict.fromkeys(genes + targets))
    with zipfile.ZipFile(A.ZIP) as archive:
        A.MEMBER = R.ACTUAL_MEMBER
        member = A.select_member(archive.infolist())
        with archive.open(member) as stream:
            mapping = A.RestrictedNumpyUnpickler(stream).load()
    selected, matrix = exact_subset(mapping, canonical)
    source_keys = len(mapping)
    del mapping
    available = set(selected)
    missing_genes = [g for g in genes if g not in available]
    missing_targets = [g for g in targets if g not in available]
    np.savez_compressed(output, genes=selected, vectors=matrix)
    passed = len(missing_targets) <= 15 and (18533-len(missing_genes))/18533 >= .7
    payload = {
        "status": "FEATURE_FEASIBLE_ONLY" if passed else "FEATURE_FEASIBILITY_NO_GO",
        "original_all_uppercase_gate": "FAILED_AND_PRESERVED",
        "identifier_policy": "canonical-exact-case-sensitive-only",
        "source_dictionary_keys": source_keys, "selected_genes": len(selected),
        "dimensions": matrix.shape[1], "target_coverage": 300-len(missing_targets),
        "response_coverage": 18533-len(missing_genes), "missing_targets": missing_targets,
        "missing_response_genes": missing_genes,
        "embedding_sha256": A.digest(output), "archive_sha256": download["sha256"],
        "canonical_inputs_sha256": A.digest(MANIFEST), "seconds": time.monotonic()-tick,
        "max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "alias_resolution": False, "biological_outcomes_read": False,
        "paid_api": False, "submission": False}
    A.write(result, payload)
    verify()
    print(json.dumps({k: v for k, v in payload.items() if not k.startswith("missing_")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "convert", "verify"))
    globals()[parser.parse_args().mode]()
