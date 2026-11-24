"""Independent, saved-statistics-only checks for the Round 5 coordinator.

Never opens either reserve source HDF5 file. Reconstructs the important
304/305 numbers without importing those workers' analysis helpers.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(4 << 20), b""):
            h.update(block)
    return h.hexdigest()


def rho(x, y):
    use = np.isfinite(x) & np.isfinite(y)
    a, b = rankdata(x[use]), rankdata(y[use])
    if use.sum() < 3 or np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def profile(counts, totals=None):
    counts = np.asarray(counts, dtype=np.float64)
    if totals is None:
        totals = counts.sum(axis=-1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.log1p(50000.0 * counts / totals)


def cos_matrix(a, b, own):
    """Literal pairwise removal of both members' own genes; n <= 175."""
    dot = a @ b.T
    aa = np.square(a).sum(axis=1)
    bb = np.square(b).sum(axis=1)
    out = np.empty(dot.shape)
    for i in range(len(a)):
        for j in range(len(b)):
            take = np.array(sorted({own[i], own[j]} - {-1}), dtype=int)
            num = dot[i, j] - (a[i, take] * b[j, take]).sum()
            na = aa[i] - np.square(a[i, take]).sum()
            nb = bb[j] - np.square(b[j, take]).sum()
            out[i, j] = num / np.sqrt(na * nb) if na > 0 and nb > 0 else np.nan
    return out


def audit304():
    e = ROOT / "experiments/302-orf-vs-crispri-symmetry/work"
    z = np.load(e / "sums302.npz", allow_pickle=False)
    symbols, axis = z["scored"].tolist(), z["axis"].tolist()
    primary = np.isin(symbols, z["primary"])
    positions = {s: i for i, s in enumerate(axis)}
    own = np.array([positions.get(s, -1) for s in symbols])
    n = len(symbols)
    counts, totals = z["acc_sym"][:n], z["lib_sym"][:n]
    a = profile(counts.sum(axis=1), totals.sum(axis=1)[:, None])
    a -= a.mean(axis=0)
    bank = ROOT / "experiments/115-gse295214-h1/work/compact-donor"
    bp = {s: i for i, s in enumerate((bank / "perts.txt").read_text().splitlines())}
    b = np.asarray(np.load(bank / "delta.npy", mmap_mode="r")[
        [bp[s] for s in symbols]][:, z["axis_h1"]], dtype=np.float64)
    b -= b.mean(axis=0)
    c = cos_matrix(a, b, own)
    sg = c.diagonal()
    ref = {r["symbol"]: r for r in csv.DictReader((e / "symmetry302.csv").open())}
    error = np.max(np.abs(sg - [float(ref[s]["S_g"]) for s in symbols]))
    assert error <= 5e-7, error
    observed_median = float(np.median(sg[primary]))
    expected = json.loads((e / "stats302.json").read_text())["C6a"]["median"]
    assert abs(observed_median - expected) < 1e-12
    gen = np.array([float(ref[s]["gen_pc20_h1"]) if ref[s]["gen_pc20_h1"]
                    else np.nan for s in symbols])
    ix = np.flatnonzero(primary)
    observed = rho(gen[ix], np.abs(sg[ix]))
    rng = np.random.default_rng(20260907)
    # A different seed/draw count is a reproducibility check, not a new claim.
    null = np.array([rho(gen[ix], np.abs(c[rng.permutation(ix), ix]))
                     for _ in range(2000)])
    ladder = {}
    for k in (1, 2, 4, 8):
        cors = []
        for draw in z["draws"]:
            x = profile(counts[:, draw[:8][:k]].sum(1),
                        totals[:, draw[:8][:k]].sum(1)[:, None])
            y = profile(counts[:, draw[8:][:k]].sum(1),
                        totals[:, draw[8:][:k]].sum(1)[:, None])
            x -= x.mean(0)
            y -= y.mean(0)
            for i, col in enumerate(own):
                if col >= 0:
                    x[i, col] = y[i, col] = 0
            den = np.linalg.norm(x, axis=1) * np.linalg.norm(y, axis=1)
            # Mirrors the historical fallback; disclose zero-depth rows below.
            cors.append(np.divide((x * y).sum(1), den,
                                  out=np.zeros(n), where=den > 0))
        ladder[str(k)] = float(np.median(np.mean(cors, axis=0)[primary]))
    return {
        "inputs": {str(p.relative_to(ROOT)): sha(p) for p in
                   [e / "sums302.npz", e / "symmetry302.csv"]},
        "csv_max_abs_error": float(error), "median_S": observed_median,
        "geometry_rho": observed,
        "uniform_permutation_draws": 2000,
        "uniform_permutation_p_plus1": float((1 + (null >= observed).sum()) / 2001),
        "uniform_null_mean": float(null.mean()),
        "uniform_null_interval95": np.quantile(null, [.025, .975]).tolist(),
        "same_saved_shard_ladder_r": ladder,
        "zero_shard_groups_all175": int((z["cells_sym"][:n] == 0).sum()),
        "limitation": "Post-hoc saved-object audit, not independent biological "
                      "validation. The shard ladder is conditional on the original "
                      "normalization, centering and nested sorted-shard design.",
    }


def holm(p):
    order = np.argsort(p, kind="stable")
    ans = np.empty(len(p))
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (len(p) - rank) * p[idx])
        ans[idx] = min(1.0, running)
    return ans


def audit305():
    e = ROOT / "experiments/305-t7-expression-test/work"
    rows = list(csv.DictReader((e / "expression305.csv").open()))
    splits = {}
    for split in ("train", "validation", "test"):
        z = np.load(e / f"sums305_{split}.npz", allow_pickle=False)
        meta = json.loads((e / f"sums305_{split}.json").read_text())
        halves = {(str(c), int(h)): a.astype(np.float64)
                  for c, h, a in zip(z["names"], z["halves"], z["sums"])}
        total = {c: halves[c, 0] + halves[c, 1] for c in meta["constructs"]}
        ctrls = [sum(total[c] for c in meta["ntc_constructs"]
                     if meta["ntc_parity"][c] == h) for h in (0, 1)]
        splits[split] = (z["genes"].tolist(), meta, halves, total, ctrls)
    verified = []
    max_error = 0.
    for r in rows:
        genes, meta, halves, total, ctrl = splits[r["split"]]
        c, g = r["construct"], r["candidate_gene"]
        j = genes.index(g)
        ca, cb = (a.copy() for a in ctrl)
        if c in meta["ntc_constructs"]:
            if meta["ntc_parity"][c] == 0:
                ca -= total[c]
            else:
                cb -= total[c]
        pooled = ca + cb
        e_pool = float(profile(total[c])[j] - profile(pooled)[j])
        e_halves = [float(profile(halves[c, h])[j] - profile(ref)[j])
                    for h, ref in enumerate((ca, cb))]
        keep = [cc for cc in meta["constructs"] if cc != c and
                meta["target_of_construct"][cc] not in ("non-targeting", g)]
        null = [float(profile(total[cc])[j] - profile(pooled)[j]) for cc in keep]
        p = (1 + sum(v <= e_pool for v in null)) / (1 + len(null))
        fc = float((total[c][j] / total[c].sum()) / (pooled[j] / pooled.sum()))
        for key, val in [("e_pool", e_pool), ("e_half0", e_halves[0]),
                         ("e_half1", e_halves[1]), ("p_raw", p), ("fold_change", fc)]:
            max_error = max(max_error, abs(float(r[key]) - val))
        assert max_error < 1e-9, (r["pair_id"], max_error)
        verified.append(dict(pair_id=r["pair_id"], group=r["group"], split=r["split"],
                             testable=r["testable"] == "True", e_pool=e_pool,
                             e_halves=e_halves, fold_change=fc, p_raw=p))
    primary = [r for r in verified if r["group"] != "ntc"]
    assert len(primary) == 11
    adj = holm([r["p_raw"] if r["testable"] else 1.0 for r in primary])
    for r, p in zip(primary, adj):
        r["p_holm_full11"] = float(p)
    ntc = {r["pair_id"] for r in verified if r["group"] == "ntc"}
    ntc_equal = all(
        np.array_equal(splits["train"][3][c], splits[s][3][c])
        for s in ("validation", "test")
        for c in splits["train"][1]["ntc_constructs"])
    return {
        "inputs": {str(p.relative_to(ROOT)): sha(p) for p in
                   [e / f"sums305_{s}.npz" for s in splits] + [e / "expression305.csv"]},
        "rows_verified": len(verified), "max_abs_error": max_error,
        "primary_family_n": 11, "primary_rows": primary,
        "primary_family_significant": int((adj <= .05).sum()),
        "unique_ntc_pairs": len(ntc), "ntc_evaluation_rows": len(verified) - len(primary),
        "ntc_sums_identical_across_splits": ntc_equal,
        "interpretation": "Empirical construct-rank diagnostics, no full-family "
                          "discovery. Repeated NTC rows are not independent controls.",
    }


def main():
    result = {"304": audit304(), "305": audit305()}
    out = ROOT / "work/round5-coordinator-review.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
