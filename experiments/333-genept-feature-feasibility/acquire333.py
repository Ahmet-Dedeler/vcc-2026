"""Public, hash-bound download; conversion is a separate guarded command."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import pickle
import subprocess
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
ZIP = WORK / "GenePT_emebdding_v2.zip"
METADATA_URL = "https://zenodo.org/api/records/10833191"
SIZE = 574395233
MD5 = "3f6ce4317e3a0091978ae5cb8fbf05a3"
MEMBER = "GenePT_gene_protein_embedding_model_3_text.pickle"


def digest(path, algorithm="sha256"):
    h = hashlib.new(algorithm)
    with path.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+".tmp")
    temp.write_text(json.dumps(data, indent=2, sort_keys=True, allow_nan=False)+"\n")
    temp.replace(path)


def inputs():
    return [HERE / "PREREGISTRATION.md", HERE / "acquire333.py", HERE / "test333.py",
            ROOT / "data/challenge/gene_names.csv", ROOT / "data/challenge/pert_counts.csv",
            ROOT / "uv.lock"]


def freeze():
    if (HERE / "inputs.json").exists():
        raise FileExistsError("already frozen")
    write(HERE / "inputs.json", {"sha256": {str(p.relative_to(ROOT)): digest(p) for p in inputs()},
                                "record": 10833191, "expected_bytes": SIZE, "expected_md5": MD5})


def verify():
    m = json.loads((HERE / "inputs.json").read_text())
    assert set(m["sha256"]) == {str(p.relative_to(ROOT)) for p in inputs()}
    for p, expected in m["sha256"].items():
        assert digest(ROOT / p) == expected, p
    own = [str(p.relative_to(ROOT)) for p in inputs() if p.parent == HERE]
    own.append(str((HERE / "inputs.json").relative_to(ROOT)))
    for p in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", p], check=True,
                       stdout=subprocess.DEVNULL, cwd=ROOT)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT)


def download():
    verify()
    if ZIP.exists():
        raise FileExistsError("archive exists; inspect rather than overwrite")
    raw = subprocess.check_output(["curl", "--fail", "--silent", "--show-error",
                                   "--max-time", "30", METADATA_URL], text=True)
    meta = json.loads(raw)
    asset, = [f for f in meta["files"] if f["key"] == ZIP.name]
    assert asset["size"] == SIZE and asset["checksum"] == "md5:"+MD5
    assert meta["metadata"]["license"]["id"] == "cc-by-4.0"
    url = asset["links"]["self"]
    expected = METADATA_URL + "/files/" + ZIP.name + "/content"
    assert url == expected
    write(WORK / "official-metadata.json", meta)
    subprocess.run(["curl", "--fail", "--location", "--show-error", "--max-time", "1200",
                    "--output", str(ZIP), url], check=True)
    assert ZIP.stat().st_size == SIZE
    assert digest(ZIP, "md5") == MD5
    write(WORK / "download.json", {"bytes": SIZE, "md5": MD5, "sha256": digest(ZIP),
                                  "record": 10833191, "url": url,
                                  "utc": datetime.now(timezone.utc).isoformat(),
                                  "paid_api": False, "biological_outcomes_read": False})
    verify()


class RestrictedNumpyUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) in {
                ("numpy", "ndarray"), ("numpy", "dtype"),
                ("numpy.core.multiarray", "_reconstruct"),
                ("numpy._core.multiarray", "_reconstruct"),
                ("numpy.core.multiarray", "scalar"),
                ("numpy._core.multiarray", "scalar"),
                ("numpy.core.numeric", "_frombuffer"),
                ("numpy._core.numeric", "_frombuffer")}:
            return super().find_class(module, name)
        raise pickle.UnpicklingError(f"unapproved pickle global: {module}.{name}")


def select_member(infos):
    for info in infos:
        p = Path(info.filename)
        if p.is_absolute() or ".." in p.parts:
            raise ValueError("unsafe archive path")
    found = [i for i in infos if Path(i.filename).name == MEMBER]
    if len(found) != 1:
        raise ValueError("expected one exact embedding member")
    if found[0].file_size > 2*1024**3:
        raise ValueError("member exceeds declared resource bound")
    return found[0]


def validate_mapping(mapping):
    if not isinstance(mapping, dict) or not mapping:
        raise ValueError("expected nonempty dictionary")
    genes = sorted(mapping)
    if any(not isinstance(g, str) or not g or g != g.upper() for g in genes):
        raise ValueError("symbols must be unique nonempty uppercase strings")
    vectors = []
    dimension = None
    for g in genes:
        v = np.asarray(mapping[g])
        if v.ndim != 1 or v.dtype.kind not in "fi" or not len(v) or not np.isfinite(v).all():
            raise ValueError("invalid vector: "+g)
        if dimension is None:
            dimension = len(v)
        if len(v) != dimension or not np.any(v):
            raise ValueError("dimension/zero-vector error: "+g)
        vectors.append(v.astype(np.float32))
    matrix = np.stack(vectors)
    if not np.isfinite(matrix).all():
        raise ValueError("float32 conversion overflow")
    return np.array(genes), matrix


def convert():
    verify()
    if (WORK / "results.json").exists() or (WORK / "embeddings.npz").exists():
        raise FileExistsError("conversion exists")
    receipt = json.loads((WORK / "download.json").read_text())
    assert ZIP.stat().st_size == SIZE and digest(ZIP) == receipt["sha256"]
    with zipfile.ZipFile(ZIP) as archive:
        infos = archive.infolist()
        member = select_member(infos)
        inventory = [{"path": i.filename, "bytes": i.file_size,
                      "compressed_bytes": i.compress_size, "crc32": i.CRC} for i in infos]
        write(WORK / "archive.json", {"files": inventory, "selected": member.filename})
        with archive.open(member) as stream:
            mapping = RestrictedNumpyUnpickler(stream).load()
    genes, matrix = validate_mapping(mapping)
    del mapping
    with (ROOT / "data/challenge/gene_names.csv").open() as f:
        response_genes = [r["gene_name"] for r in csv.DictReader(f)]
    with (ROOT / "data/challenge/pert_counts.csv").open() as f:
        targets = [r["target_gene"] for r in csv.DictReader(f)]
    assert len(response_genes) == 18533 and len(targets) == 300
    available = set(genes)
    target_missing = sorted(set(targets)-available)
    response_missing = sorted(set(response_genes)-available)
    np.savez_compressed(WORK / "embeddings.npz", genes=genes, vectors=matrix)
    passed = len(target_missing) <= 15 and (18533-len(response_missing))/18533 >= .7
    write(WORK / "results.json", {
        "status": "FEATURE_FEASIBLE_ONLY" if passed else "FEATURE_FEASIBILITY_NO_GO",
        "n_genes": len(genes), "dimensions": matrix.shape[1],
        "target_coverage": 300-len(target_missing), "response_coverage": 18533-len(response_missing),
        "missing_targets": target_missing, "missing_response_genes": response_missing,
        "embedding_sha256": digest(WORK / "embeddings.npz"),
        "archive_sha256": receipt["sha256"], "inputs_sha256": digest(HERE / "inputs.json"),
        "biological_outcomes_read": False, "paid_api": False, "submission": False})
    verify()
    print(f"Feature schema {matrix.shape}; targets {300-len(target_missing)}/300; "
          f"genes {18533-len(response_missing)}/18533; feasible={passed}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "download", "convert", "verify"))
    globals()[parser.parse_args().mode]()
