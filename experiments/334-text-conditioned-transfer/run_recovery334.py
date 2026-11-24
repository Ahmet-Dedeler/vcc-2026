"""Explicit pre-outcome text-source recovery; original experiment stays frozen."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile
import numpy as np
from threadpoolctl import threadpool_limits
import core334 as C
import run334 as R

FALLBACK_SHA = "319897d7e2f8a0306bccdf2d022982b3c8321dbec8273739094021d57e414c31"
ORIGINAL_VERIFY = R.verify
RECOVERY_WORK = R.HERE/"recovery-work"


def sources(primary, fallback, keys):
    docs, provenance = [], []
    for key in keys:
        text = C.description(primary, key)
        source = "NCBI" if text else "missing"
        if not text:
            text = C.description(fallback, key)
            if text:
                source = "NCBI_UniProt_fallback"
        docs.append(text)
        provenance.append(source)
    return docs, provenance


def paths():
    return [R.HERE/n for n in ("inputs.json", "COVERAGE-FAILURE.md",
                              "RECOVERY-001.md", "run_recovery334.py",
                              "test_recovery334.py")]


def freeze():
    p = R.HERE/"recovery-inputs.json"
    if p.exists():
        raise FileExistsError("already frozen")
    ORIGINAL_VERIFY()
    R.write(p, {"sha256": {str(p.relative_to(R.ROOT)): R.sha(p) for p in paths()},
                "fallback_member_sha256": FALLBACK_SHA,
                "biological_fit_already_run": False})


def verify():
    ORIGINAL_VERIFY()
    record = json.loads((R.HERE/"recovery-inputs.json").read_text())
    assert set(record["sha256"]) == {str(p.relative_to(R.ROOT)) for p in paths()}
    for name, expected in record["sha256"].items():
        assert R.sha(R.ROOT/name) == expected, name
    own = [str(p.relative_to(R.ROOT)) for p in paths()]
    own += [str((R.HERE/"recovery-inputs.json").relative_to(R.ROOT))]
    for name in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", name], cwd=R.ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=R.ROOT)


def representations(targets, genes, train, seed):
    keys = targets+genes
    fit = np.r_[train, np.arange(len(targets), len(keys))]
    mappings = []
    with zipfile.ZipFile(R.ZIP) as z:
        for suffix, expected in (("NCBI_summary_of_genes.json", R.TEXT_SHA),
                                 ("NCBI_UniProt_summary_of_genes.json", FALLBACK_SHA)):
            member, = [n for n in z.namelist() if Path(n).name == suffix]
            raw = z.read(member)
            assert hashlib.sha256(raw).hexdigest() == expected
            mappings.append(json.loads(raw))
    docs, provenance = sources(*mappings, keys)
    have = np.array([bool(t) for t in docs])
    nt = len(targets)
    assert have[:nt].sum() >= .95*nt and have[nt:].sum() >= .9*len(genes)
    text, tfidf, svd = C.text_representation(docs, fit, seed)
    with np.load(R.STRING, allow_pickle=False) as z:
        index = {g:i for i,g in enumerate(z["genes"])}
        val = np.zeros((len(keys), z["unit"].shape[1]))
        available = np.zeros(len(keys), bool)
        for i,g in enumerate(keys):
            if g in index and z["available"][index[g]]:
                val[i] = z["unit"][index[g]]
                available[i] = True
    string, pca = C.string_representation(val, available, fit, seed)
    pt = C.derangement(nt, C.SEED+10)
    pg = C.derangement(len(genes), C.SEED+11)
    blocks = {"source_only": (np.zeros_like(text[:nt]), np.zeros_like(text[nt:])),
              "string": (string[:nt], string[nt:]),
              "text": (text[:nt], text[nt:]),
              "text_shuffled": (text[:nt][pt], text[nt:][pg]),
              "labels_shuffled": (text[:nt], text[nt:])}
    meta = {"source_policy": "RECOVERY-001.md",
            "text_targets": int(have[:nt].sum()), "text_genes": int(have[nt:].sum()),
            "string_targets": int(available[:nt].sum()),
            "string_genes": int(available[nt:].sum()),
            "vocabulary_size": len(tfidf.vocabulary_), "fit_document_rows": fit.tolist(),
            "primary_member_sha256": R.TEXT_SHA, "fallback_member_sha256": FALLBACK_SHA,
            "description_sources": dict(zip(keys, provenance)),
            "description_sha256": {k:hashlib.sha256(d.encode()).hexdigest()
                                   for k,d in zip(keys, docs)},
            "recovery_inputs_sha256": R.sha(R.HERE/"recovery-inputs.json")}
    return blocks, (tfidf, svd, pca), meta


def run():
    verify()
    R.WORK = RECOVERY_WORK
    R.verify = verify
    R.representations = representations
    R.run()
    R.write(RECOVERY_WORK/"recovery-attestation.json",
            {"recovery_inputs_sha256": R.sha(R.HERE/"recovery-inputs.json"),
             "original_inputs_sha256": R.sha(R.HERE/"inputs.json"),
             "result_sha256": R.sha(RECOVERY_WORK/"results.json"),
             "execution_sha256": R.sha(RECOVERY_WORK/"execution.json"),
             "source_policy": "RECOVERY-001.md",
             "models_or_gates_changed": False})


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("action", choices=("freeze", "verify", "run"))
    with threadpool_limits(limits=6):
        globals()[p.parse_args().action]()
