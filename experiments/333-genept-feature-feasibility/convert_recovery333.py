"""Filename-only engineering correction, preserving the original freeze."""

import argparse
import json
from pathlib import Path
import subprocess
import acquire333 as A

ACTUAL_MEMBER = "GenePT_gene_protein_embedding_model_3_text.pickle."
MANIFEST = A.HERE / "recovery-inputs.json"
FILES = [Path(__file__).resolve(), A.HERE / "RECOVERY-001.md",
         A.HERE / "inputs.json", A.WORK / "download.json"]
BASE_VERIFY = A.verify


def freeze():
    BASE_VERIFY()
    if MANIFEST.exists():
        raise FileExistsError("recovery already frozen")
    receipt = json.loads((A.WORK / "download.json").read_text())
    assert receipt["sha256"] == "6193575dcbd7bbf214c8ca3eb518cb3d13272443f98ecd6c26402411a7ca745e"
    A.write(MANIFEST, {"sha256": {str(p.relative_to(A.ROOT)): A.digest(p) for p in FILES},
                       "replacement_basename": ACTUAL_MEMBER})


def verify():
    BASE_VERIFY()
    m = json.loads(MANIFEST.read_text())
    assert m["replacement_basename"] == ACTUAL_MEMBER
    assert set(m["sha256"]) == {str(p.relative_to(A.ROOT)) for p in FILES}
    for name, expected in m["sha256"].items():
        assert A.digest(A.ROOT / name) == expected, name
    own = [str(p.relative_to(A.ROOT)) for p in FILES] + [str(MANIFEST.relative_to(A.ROOT))]
    for name in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", name],
                       cwd=A.ROOT, check=True, stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=A.ROOT)


def convert():
    verify()
    A.MEMBER = ACTUAL_MEMBER
    A.verify = verify
    A.convert()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("freeze", "convert", "verify"))
    globals()[parser.parse_args().mode]()
