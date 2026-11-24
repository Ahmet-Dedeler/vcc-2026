"""Frozen local response pilot. No vector pickle, network or production writes."""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import pickle
import resource
import subprocess
import time
import zipfile

import numpy as np
from threadpoolctl import threadpool_limits
import core334 as C

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BANKS = ("data/donors/k562_gwps", "data/donors/orion_hct116",
         "data/donors/orion_hek293t", "data/donors/rpe1",
         "data/donors/nadig_jurkat",
         "experiments/166-donor-acquisition/donors/nadig_hepg2")
LINES = tuple(Path(p).name for p in BANKS)
ZIP = ROOT / "experiments/333-genept-feature-feasibility/work/GenePT_emebdding_v2.zip"
TEXT_SHA = "3db721b7cfbc35795428c6a7ab4e8a9c59fb5f5ccb8c64af50dd3c97c2d4de2e"
STRING = ROOT / "experiments/100-sign-confidence/work/string-pca.npz"
HIST = ROOT / "experiments/224-learned-transfer/work/surface-S5.npz"


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for block in iter(lambda: f.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write(p, obj):
    p.parent.mkdir(parents=True, exist_ok=True)
    q = p.with_suffix(p.suffix+".tmp")
    q.write_text(json.dumps(obj, indent=2, sort_keys=True, allow_nan=False)+"\n")
    q.replace(p)


def names(p):
    return p.read_text().splitlines()


def paths():
    own = [HERE/n for n in ("PREREGISTRATION.md", "core334.py", "run334.py",
                            "report334.py", "test334.py", "selection.json")]
    data = [ROOT/b/n for b in BANKS for n in
            ("perts.txt", "genes.txt", "delta.npy", "b_ctrl.npy", "meta.json")]
    return own + data + [ZIP, STRING, HIST, ROOT/"uv.lock",
                         ROOT/"data/challenge/pert_counts.csv",
                         ROOT/"data/challenge/gene_names.csv",
                         ROOT/"experiments/323-current-family-validation/inputs.json"]


def freeze():
    if (HERE/"inputs.json").exists() or (HERE/"selection.json").exists():
        raise FileExistsError("already selected/frozen")
    ps = [set(names(ROOT/b/"perts.txt")) for b in BANKS]
    gs = [set(names(ROOT/b/"genes.txt")) for b in BANKS]
    panel = {r["target_gene"] for r in csv.DictReader(
        (ROOT/"data/challenge/pert_counts.csv").open())}
    panel.update(json.loads((ROOT/"experiments/323-current-family-validation/inputs.json")
                            .read_text())["h1_targets"])
    common = set.intersection(*ps)
    target_pool = sorted(common-panel)
    canonical = {r["gene_name"] for r in csv.DictReader(
        (ROOT/"data/challenge/gene_names.csv").open())}
    gene_pool = sorted((set.intersection(*gs) & canonical)-common-panel)
    assert len(target_pool) >= C.N_TARGETS and len(gene_pool) >= C.N_GENES
    targets = [target_pool[i] for i in C.order(target_pool, "334-target|")[:C.N_TARGETS]]
    genes = [gene_pool[i] for i in C.order(gene_pool, "334-gene|")[:C.N_GENES]]
    meta = [json.loads((ROOT/b/"meta.json").read_text()) for b in BANKS]
    assert all(m["target_sum"] == 50000. for m in meta)
    depth = [m["cells_per_pert"]["median"] for m in meta]
    write(HERE/"selection.json", {"targets": targets, "genes": genes,
          "target_pool": len(target_pool), "gene_pool": len(gene_pool),
          "lines": LINES, "median_cells_per_perturbation": depth,
          "lower_cell_lines": np.argsort(depth)[:3].tolist(),
          "higher_cell_lines": np.argsort(depth)[3:].tolist(),
          "targets_excluded": sorted(panel), "no_response_values_read": True})
    write(HERE/"inputs.json", {"sha256": {str(p.relative_to(ROOT)): sha(p) for p in paths()}})


def verify():
    m = json.loads((HERE/"inputs.json").read_text())
    assert set(m["sha256"]) == {str(p.relative_to(ROOT)) for p in paths()}
    for p, value in m["sha256"].items():
        if sha(ROOT/p) != value:
            raise ValueError("input drift: "+p)
    own = [str(p.relative_to(ROOT)) for p in paths() if p.parent == HERE]
    own += [str((HERE/"inputs.json").relative_to(ROOT))]
    for p in own:
        subprocess.run(["git", "ls-files", "--error-unmatch", p], cwd=ROOT,
                       check=True, stdout=subprocess.DEVNULL)
    assert not subprocess.check_output(["git", "status", "--porcelain", "--", *own], cwd=ROOT)


def load(bank, targets, genes):
    p = ROOT/bank
    pi = {t:i for i,t in enumerate(names(p/"perts.txt"))}
    gi = {g:i for i,g in enumerate(names(p/"genes.txt"))}
    ti, ix = [pi[t] for t in targets], [gi[g] for g in genes]
    raw = np.load(p/"delta.npy", mmap_mode="r")
    out = np.array(raw[np.ix_(ti, ix)], dtype=np.float64)
    basal = np.array(np.load(p/"b_ctrl.npy")[ix], dtype=np.float64)
    assert np.isfinite(out).all() and np.isfinite(basal).all()
    # Direct scalar-index fixture is not a replacement for historical check.
    np.testing.assert_array_equal(out[:4, :4],
        [[raw[pi[t], gi[g]] for g in genes[:4]] for t in targets[:4]])
    return out, basal


def known():
    receipt = []
    with np.load(HIST, allow_pickle=False) as z:
        members, targets, genes = list(z["members"]), list(z["targets"]), list(z["genes"])
        cached = z["C"]
        for bank in BANKS[:5]:
            line = Path(bank).name
            raw, _ = load(bank, targets, genes)
            recreated = (raw-raw.mean(0, keepdims=True)).astype(np.float32)
            original = cached[members.index(line)]
            difference = float(np.abs(recreated-original).max())
            assert difference == 0, (line, difference)
            receipt.append({"line": line, "maximum_historical_difference": difference})
    return receipt


def representations(targets, genes, train, seed):
    all_names = targets+genes
    fit = np.r_[train, np.arange(len(targets), len(all_names))]
    with zipfile.ZipFile(ZIP) as z:
        member, = [n for n in z.namelist() if Path(n).name == "NCBI_summary_of_genes.json"]
        raw = z.read(member)
        assert hashlib.sha256(raw).hexdigest() == TEXT_SHA
        mapping = json.loads(raw)
    docs = [C.description(mapping, g) for g in all_names]
    have = np.array([bool(d) for d in docs])
    tc, gc = int(have[:len(targets)].sum()), int(have[len(targets):].sum())
    assert tc >= .95*len(targets) and gc >= .9*len(genes), (tc, gc)
    text, tfidf, svd = C.text_representation(docs, fit, seed)
    with np.load(STRING, allow_pickle=False) as z:
        index = {g:i for i,g in enumerate(z["genes"])}
        val = np.zeros((len(all_names), z["unit"].shape[1]))
        available = np.zeros(len(all_names), bool)
        for i, g in enumerate(all_names):
            if g in index and z["available"][index[g]]:
                val[i] = z["unit"][index[g]]
                available[i] = True
    string, pca = C.string_representation(val, available, fit, seed)
    perm_t = C.derangement(len(targets), C.SEED+10)
    perm_g = C.derangement(len(genes), C.SEED+11)
    nt = len(targets)
    blocks = {"source_only": (np.zeros_like(text[:nt]), np.zeros_like(text[nt:])),
              "string": (string[:nt], string[nt:]),
              "text": (text[:nt], text[nt:]),
              "text_shuffled": (text[:nt][perm_t], text[nt:][perm_g]),
              "labels_shuffled": (text[:nt], text[nt:])}
    meta = {"text_targets": tc, "text_genes": gc,
            "string_targets": int(available[:nt].sum()),
            "string_genes": int(available[nt:].sum()),
            "vocabulary_size": len(tfidf.vocabulary_), "fit_document_rows": fit.tolist(),
            "text_identity_failures_or_missing": [g for g,d in zip(all_names, docs) if not d],
            "text_sha256": TEXT_SHA}
    return blocks, (tfidf, svd, pca), meta


def run():
    verify()
    if (WORK/"execution.json").exists() or list(WORK.glob("*/prediction-*.npz")):
        raise FileExistsError("do not overwrite a prior run")
    start = time.monotonic()
    work_start = datetime.now(timezone.utc).isoformat()
    selection = json.loads((HERE/"selection.json").read_text())
    targets, genes = selection["targets"], selection["genes"]
    checks = known()
    loaded = [load(b, targets, genes) for b in BANKS]
    raw, basal = np.array([a for a,b in loaded]), np.array([b for a,b in loaded])
    direction, rms = C.rms_direction(raw)
    del loaded, raw
    work = []
    permutation = C.order(targets, "334-split|")
    tblocks = [np.array(permutation[::2]), np.array(permutation[1::2])]
    static = {}
    for block in (0, 1):
        static[block] = representations(targets, genes, tblocks[1-block], C.SEED+block)
    write(WORK/"features.json", {str(b): static[b][2] for b in static})
    for outer, line in enumerate(LINES):
        directory = WORK/line
        directory.mkdir(parents=True, exist_ok=True)
        baseline = direction[[i for i in range(len(LINES)) if i != outer]].mean(0)
        evaluation = directory/"evaluation.npz"
        np.savez_compressed(evaluation, truth=direction[outer], baseline=baseline,
                            targets=targets, genes=genes)
        receipts = []
        for block in (0, 1):
            if time.monotonic()-start > 5400:
                raise TimeoutError("registered 90-minute budget exhausted")
            split_start = time.monotonic()
            train, test = tblocks[1-block], tblocks[block]
            inner = C.source_blocks(direction, basal, rms, outer, train)
            query = C.query_blocks(direction, basal, rms, outer, test)
            gain = C.diagonal(inner)
            d2, r2 = direction.copy(), rms.copy()
            d2[outer], r2[outer] = 1e6, 1e6
            poisoned = C.source_blocks(d2, basal, r2, outer, train)
            for (_, x, y), (_, xp, yp) in zip(inner, poisoned):
                np.testing.assert_array_equal(x, xp)
                np.testing.assert_array_equal(y, yp)
            for x, xp in zip(query, C.query_blocks(d2, basal, r2, outer, test)):
                np.testing.assert_array_equal(x, xp)
            np.testing.assert_array_equal(gain, C.diagonal(poisoned))
            np.testing.assert_allclose(np.mean([q[..., 0] for q in query], 0),
                                       baseline[test], atol=1e-14, rtol=1e-14)
            ix = C.sample_indices(5, len(train), len(genes), C.N_TRAIN,
                                  C.SEED+100*outer+block)
            assert not set(train)&set(test)
            assert all(r != outer for r,_,_ in inner)
            x = np.empty((C.N_TRAIN, len(C.SOURCE_FEATURES)))
            y, ys = np.empty(C.N_TRAIN), np.empty(C.N_TRAIN)
            labperm = C.derangement(len(train), C.SEED+block+20)
            for j, (_, block_x, block_y) in enumerate(inner):
                mask = ix[:, 0] == j
                t, g = ix[mask, 1], ix[mask, 2]
                x[mask] = block_x[t, g]
                y[mask] = block_y[t, g]
                ys[mask] = block_y[labperm[t], g]
            training_path = directory/f"training-{block}.npz"
            np.savez_compressed(training_path, x=x, y=y, y_shuffled=ys,
                                indices=ix, train_indices=train, test_indices=test)
            vectors, transforms, _ = static[block]
            models = {}
            outputs = {"source_mean": baseline[test], "zero": np.zeros_like(baseline[test]),
                       "diagonal": baseline[test]*gain}
            for arm in C.FITTED:
                tv, gv = vectors[arm]
                f = C.pair_features(tv, gv, train[ix[:, 1]], ix[:, 2])
                fit_x = np.column_stack((x, f))
                response = ys if arm == "labels_shuffled" else y
                fitted = C.model(C.SEED+100*outer+block)
                fitted.fit(fit_x, response-x[:, 0])
                models[arm] = fitted
                outputs[arm] = np.zeros_like(baseline[test])
                ti = np.repeat(test, len(genes))
                gi = np.tile(np.arange(len(genes)), len(test))
                pairs = C.pair_features(tv, gv, ti, gi)
                for q in query:
                    flat = q.reshape(-1, len(C.SOURCE_FEATURES))
                    # Chunked inference keeps 49-feature pair matrix bounded.
                    predicted = np.empty(len(flat))
                    for s in range(0, len(flat), 8192):
                        v = np.column_stack((flat[s:s+8192], pairs[s:s+8192]))
                        predicted[s:s+8192] = fitted.predict(v)+flat[s:s+8192, 0]
                    outputs[arm] += predicted.reshape(len(test), len(genes))/len(query)
            for a in outputs.values():
                assert a.shape == (len(test), len(genes)) and np.isfinite(a).all()
            path = directory/f"prediction-{block}.npz"
            np.savez_compressed(path, test_indices=test, **outputs)
            model_path = directory/f"models-{block}.pickle"
            with model_path.open("wb") as f:
                pickle.dump({"models": models, "transforms": transforms, "vectors": vectors,
                             "gain": gain}, f, protocol=5)
            record = {"outer": line, "block": block, "seconds": time.monotonic()-split_start,
                      "poison_isolation_pass": True,
                      "training_sha256": sha(training_path), "prediction_sha256": sha(path),
                      "models_sha256": sha(model_path),
                      "evaluation_sha256": sha(evaluation),
                      "inputs_sha256": sha(HERE/"inputs.json"),
                      "selection_sha256": sha(HERE/"selection.json")}
            write(directory/f"fit-{block}.json", record)
            receipts.append(record)
            print(f"{line} block {block} complete; {record['seconds']:.1f}s", flush=True)
        work.extend(receipts)
    verify()
    write(WORK/"execution.json", {"start_utc": work_start,
          "finish_utc": datetime.now(timezone.utc).isoformat(),
          "seconds": time.monotonic()-start, "historical_checks": checks,
          "max_rss_bytes_macos": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
          "fits": work, "inputs_sha256": sha(HERE/"inputs.json")})
    from report334 import report
    report()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("freeze", "verify", "run"))
    args = parser.parse_args()
    with threadpool_limits(limits=6):
        globals()[args.action]()
