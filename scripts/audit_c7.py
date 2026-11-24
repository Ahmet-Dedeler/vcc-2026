"""Independent arithmetic/provenance audit of C7, saved product only.

No import of the C7 analysis kernel/driver, and no source HDF5 access.
This checks the frozen result; it does not choose another statistic.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]
EXP = ROOT / "experiments/303-second-overexpression-library"
WORK = EXP / "work"


def digest(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def correlation(x, y):
    a, b = rankdata(x), rankdata(y)
    return float(np.corrcoef(a, b)[0, 1])


def main():
    recorded = json.loads((WORK / "stats303.json").read_text())
    access = json.loads((WORK / "access303.json").read_text())
    release = json.loads((WORK / "RELEASE.json").read_text())
    assert digest(WORK / "nak303.npz") == access["product_sha256"]
    assert digest(WORK / "RELEASE.json") == access["release_sha256"]
    assert digest(EXP / "PREREGISTRATION.md") == recorded["provenance"]["final_preregistration_sha256"]
    for rel, expected in recorded["provenance"]["input_sha256"].items():
        assert digest(ROOT / rel) == expected, rel
    cohort = list(csv.DictReader((WORK / "COHORT.tsv").open(), delimiter="\t"))
    symbols = [r["symbol"] for r in cohort]
    with np.load(WORK / "nak303.npz", allow_pickle=False) as z:
        axis = z["axis"].tolist()
        op = {s: i for i, s in enumerate(z["symbols"].tolist())}
        O = z["delta"][[op[s] for s in symbols]].astype(float)
        assert np.all(z["mask"])
    bank = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
    columns = {g: i for i, g in enumerate((bank / "genes.txt").read_text().splitlines())}
    rows = {p: i for i, p in enumerate((bank / "perts.txt").read_text().splitlines())}
    K = np.asarray(np.load(bank / "delta.npy", mmap_mode="r")[
        [rows[s] for s in symbols]][:, [columns[g] for g in axis]], float)
    # Independent explicit fixed-axis implementation.
    keep = np.array([g not in set(symbols) for g in axis])
    assert keep.sum() == 15211
    O = (O - O.mean(0))[:, keep]
    K = (K - K.mean(0))[:, keep]
    C = (K @ O.T) / (np.linalg.norm(K, axis=1)[:, None] *
                      np.linalg.norm(O, axis=1)[None, :])
    direct = np.array([np.dot(k, o) / (np.linalg.norm(k) * np.linalg.norm(o))
                       for k, o in zip(K, O)])
    with np.load(WORK / "analysis303.npz", allow_pickle=False) as a:
        error = float(np.max(np.abs(C - a["cosine_matrix"])))
        assert error < 1e-12
        assert np.max(np.abs(direct - a["S"])) < 1e-12
        rng = np.random.default_rng(20260907)
        perms = np.stack([rng.permutation(100) for _ in range(2000)])
        assert np.array_equal(perms, a["permutations"])
        null = C[np.arange(100), perms]
        assert np.max(np.abs(null - a["null_S"])) < 1e-12
    tail = lambda s: np.mean(s < -.05, axis=-1) - np.mean(s > .05, axis=-1)
    p_med = (1 + np.sum(np.median(null, axis=1) <= np.median(direct))) / 2001
    p_tail = (1 + np.sum(tail(null) >= tail(direct))) / 2001
    p = {"C7a": float(max(p_med, p_tail))}
    heldpath = ROOT / "experiments/279-znf-state-dependence/work/held279.npz"
    with np.load(heldpath, allow_pickle=True) as z:
        gp = {g: i for i, g in enumerate(z["genes"].tolist())}
        q = z["keys"].tolist().index("gse295214_h1")
        used = (z["pairs_cross"] != q).all(axis=1)
        assert int(used.sum()) == 338
        a = z["Ccross"][used]
        count = np.isfinite(a).sum(0)
        mean = np.divide(np.nansum(a, axis=0), count,
                         out=np.full(len(count), np.nan), where=count > 0)
        conservation = np.array([mean[gp[s]] if s in gp else np.nan for s in symbols])
    with np.load(ROOT / "experiments/290-remainder-response-features/work/persrc290.npz",
                 allow_pickle=True) as z:
        gp = {g: i for i, g in enumerate(z["genes"].tolist())}
        values = z["pc20"][z["keys"].tolist().index("gse295214_h1")]
        geometry = np.array([values[gp[s]] if s in gp else np.nan for s in symbols])
    effects = {}
    for name, x, y, yn in [("C7b", conservation, -direct, -null),
                           ("C7c", geometry, -direct, -null),
                           ("C7d", geometry, np.abs(direct), np.abs(null))]:
        take = np.isfinite(x)
        assert take.sum() == 97
        r = correlation(x[take], y[take])
        rn = np.array([correlation(x[take], v[take]) for v in yn])
        p[name] = float((1 + np.sum(rn >= r)) / 2001)
        effects[name] = r
        assert abs(r - recorded["cells"][name]["observed"]) < 1e-12
    adjusted = {}
    largest = 0.
    for i, name in enumerate(sorted(p, key=p.get)):
        largest = max(largest, (4 - i) * p[name])
        adjusted[name] = min(1., largest)
    for name in p:
        assert abs(p[name] - recorded["cells"][name]["p_raw"]) < 1e-12
        assert abs(adjusted[name] - recorded["cells"][name]["p_holm"]) < 1e-12
    assert np.median(direct) < 0 and tail(direct) > 0
    assert access["cpm_gate_pass"] is True and recorded["validity"]["validity"] == "VALID"
    result = dict(
        status="VERIFIED", source_reopened=False, experiment="303", claim="C7",
        freeze_manifest_sha256=access["manifest_sha256"],
        source_product_sha256=access["product_sha256"],
        readout_sha256=digest(WORK / "stats303.json"),
        report_sha256=digest(EXP / "RESULTS.md"),
        audit_code_sha256=digest(Path(__file__)),
        max_cosine_matrix_difference=error,
        median=float(np.median(direct)), n_negative=int((direct < 0).sum()),
        n_below_minus005=int((direct < -.05).sum()),
        n_above_plus005=int((direct > .05).sum()),
        association_effects=effects, p_raw=p, p_holm=adjusted,
        verdicts={k: v["verdict"] for k, v in recorded["cells"].items()},
        limitations="Checks arithmetic/provenance, not independent biological replication.")
    (WORK / "COORDINATOR-VERIFICATION.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
