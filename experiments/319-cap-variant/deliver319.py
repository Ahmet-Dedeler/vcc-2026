"""319 addendum -- the CAP arm's DELIVERED magnitude profile on the REAL 300-target panel.

  PYTHONDONTWRITEBYTECODE=1 uv run python experiments/319-cap-variant/deliver319.py

318's Correction 5 is the reason this exists: every per-target norm this family has quoted as
"delivered" was measured on 313's 20-target synthetic FIXTURE, and the real panel reads 6x narrower.
So the cap arm's own delivered profile is measured here where it acts -- on the emitted counts of
the 12.6 GB prediction -- and never quoted from `test319.py`.

Recipe, `experiments/318-per-target-magnitude/profile318.py` verbatim (the delivered log-effect
every DE member of the scorer reads):

    p(t,c)  = the 400 emitted count vectors of target t in context c, summed
    effect  = log1p(50000 * p / p.sum())  -  b_ctrl(c)      # b_ctrl from the challenge's controls
    m(t,c)  = || effect_arm || / || effect_reference ||

KNOWN-ANSWER ROW: the shipped-vs-parity pair is recomputed first and must reproduce 318's published
summary (`work/profile318.json`: median 1.011353, sd(log) 0.032885, whole matrix 1.0085497).

Read-only outside `experiments/319-cap-variant/`.
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import h5py
import numpy as np

REPO = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "work"
E318 = REPO / "experiments/318-per-target-magnitude/work/profile318.json"
TS = 50_000.0

BUILDS = {
    "incumbent": REPO / "submissions/249-cd4svd5-f034.h5ad",
    "shipped": REPO / "experiments/312-normrel-arm/work/build/the-cat-sat-on-the-modem.h5ad",
    "parity": REPO / "experiments/316-parity-arm/work/build/a-duck-reviews-the-diff.h5ad",
    "cap": REPO / "experiments/319-cap-variant/work/build/the-parrot-approved-my-expenses.h5ad",
}
CONTROLS = {c: REPO / f"data/challenge/context_{c}.h5ad" for c in ("A", "B", "C")}


def rss() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3) if sys.platform == "darwin" else r / (1024 ** 2)


def labels(h, path):
    o = h[path]
    if hasattr(o, "keys") and "categories" in o:
        cat = np.asarray([x.decode() if isinstance(x, bytes) else str(x)
                          for x in o["categories"][:]])
        return cat[o["codes"][:]]
    v = o["values"][:] if hasattr(o, "keys") else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def control_profile(path: Path, n_genes: int) -> np.ndarray:
    with h5py.File(path, "r") as h:
        X = h["X"]
        shape = tuple(int(x) for x in X.attrs["shape"])
        indptr = X["indptr"][:]
        tot = np.zeros(n_genes, dtype=np.float64)
        step = 4000
        for s in range(0, shape[0], step):
            e = min(s + step, shape[0])
            lo, hi = int(indptr[s]), int(indptr[e])
            tot += np.bincount(X["indices"][lo:hi].astype(np.int64),
                               weights=X["data"][lo:hi].astype(np.float64),
                               minlength=n_genes)
    return np.log1p(TS * tot / tot.sum())


def block_effects(path: Path, b_ctrl: dict[str, np.ndarray]):
    with h5py.File(path, "r") as h:
        X = h["X"]
        n_cells, n_genes = (int(x) for x in X.attrs["shape"])
        indptr = X["indptr"][:]
        ctx = labels(h, "obs/context")
        tgt = labels(h, "obs/target_gene")
        key = np.char.add(np.char.add(ctx, "|"), tgt)
        bounds = np.flatnonzero(np.r_[True, key[1:] != key[:-1], True])
        out = {}
        data, indices = X["data"], X["indices"]
        for bi in range(bounds.size - 1):
            s, e = int(bounds[bi]), int(bounds[bi + 1])
            lo, hi = int(indptr[s]), int(indptr[e])
            p = np.bincount(indices[lo:hi].astype(np.int64),
                            weights=data[lo:hi].astype(np.float64), minlength=n_genes)
            eff = np.log1p(TS * p / p.sum()) - b_ctrl[ctx[s]]
            out[(ctx[s], tgt[s])] = eff.astype(np.float32)
    return out, n_genes


def summarise(eff_a: dict, eff_b: dict, non_panel: np.ndarray) -> tuple[dict, list]:
    rows = []
    for kk in sorted(eff_b):
        a = eff_a[kk].astype(np.float64)
        p = eff_b[kk].astype(np.float64)
        na, nb = float(np.linalg.norm(a)), float(np.linalg.norm(p))
        naw = float(np.linalg.norm(a[non_panel]))
        nbw = float(np.linalg.norm(p[non_panel]))
        rows.append({"context": kk[0], "target": kk[1], "norm_a": na, "norm_b": nb,
                     "m_full": na / nb, "m_nonpanel": naw / nbw,
                     "row_cosine": float(a @ p / max(na * nb, 1e-300))})
    m = np.array([r["m_full"] for r in rows])
    whole = float(np.sqrt(sum(r["norm_a"] ** 2 for r in rows)
                          / sum(r["norm_b"] ** 2 for r in rows)))
    return {
        "n_blocks": len(rows),
        "m_full": {"median": float(np.median(m)), "mean": float(m.mean()),
                   "min": float(m.min()), "max": float(m.max()),
                   "sd_log": float(np.std(np.log(m))),
                   "q05": float(np.quantile(m, 0.05)), "q95": float(np.quantile(m, 0.95)),
                   "frac_gt1": float((m > 1).mean()),
                   "frac_within_1e6": float((np.abs(m - 1) < 1e-6).mean())},
        "whole_matrix_norm_ratio": whole,
        "row_cosine": {"median": float(np.median([r["row_cosine"] for r in rows])),
                       "min": float(min(r["row_cosine"] for r in rows))},
    }, rows


def main() -> None:
    t0 = time.time()
    with h5py.File(BUILDS["cap"], "r") as h:
        genes = labels(h, "var/_index")
    n_genes = genes.size
    panel_targets = set(np.loadtxt(REPO / "data/challenge/pert_counts.csv", dtype=str,
                                   delimiter=",", skiprows=1, usecols=0).tolist())
    panel_pos = np.asarray([i for i, g in enumerate(genes) if g in panel_targets], dtype=np.int64)
    non_panel = np.setdiff1d(np.arange(n_genes), panel_pos)
    print(f"[319] {n_genes} genes; {len(panel_targets)} panel targets, {panel_pos.size} on the "
          f"gene axis", flush=True)

    b = {c: control_profile(p, n_genes) for c, p in CONTROLS.items()}
    eff = {}
    for name, path in BUILDS.items():
        eff[name], _ = block_effects(path, b)
        print(f"[319] {name}: {len(eff[name])} blocks ({time.time() - t0:.0f}s, "
              f"RSS {rss():.2f} GiB)", flush=True)

    KA = []
    ref = json.loads(E318.read_text())["summary"]
    s_par, _ = summarise(eff["shipped"], eff["parity"], non_panel)
    for k, got, want in (("median", s_par["m_full"]["median"], ref["m_full"]["median"]),
                         ("sd_log", s_par["m_full"]["sd_log"], ref["m_full"]["sd_log"]),
                         ("whole", s_par["whole_matrix_norm_ratio"],
                          ref["whole_matrix_norm_ratio"])):
        d = abs(float(got) - float(want))
        KA.append(dict(label=f"KA-D 318's shipped/parity {k}", got=float(got),
                       want=float(want), abs_diff=d, tol=1e-9, ok=bool(d <= 1e-9)))
        print(f"  KA {'PASS' if d <= 1e-9 else 'FAIL'} 318's shipped/parity {k}: "
              f"{got!r} vs {want!r} |d| {d:.3e}", flush=True)

    out = {"known_answers": KA, "pairs": {}}
    for label, a, bname in (("shipped/incumbent", "shipped", "incumbent"),
                            ("cap/incumbent", "cap", "incumbent"),
                            ("cap/shipped", "cap", "shipped"),
                            ("shipped/parity (318's pair)", "shipped", "parity")):
        s, rows = summarise(eff[a], eff[bname], non_panel)
        out["pairs"][label] = s
        mf = s["m_full"]
        print(f"[319] {label:28s} m median {mf['median']:.6f} range "
              f"[{mf['min']:.4f}, {mf['max']:.4f}] sd(log) {mf['sd_log']:.6f} "
              f"whole {s['whole_matrix_norm_ratio']:.7f} row cos median "
              f"{s['row_cosine']['median']:.6f} min {s['row_cosine']['min']:.6f}", flush=True)
        if label == "cap/incumbent":
            out["cap_rows"] = rows
    out["peak_rss_gib"] = rss()
    out["seconds"] = round(time.time() - t0, 1)
    (OUT / "deliver319.json").write_text(json.dumps(out, indent=1))
    print(f"[319] KA rows {sum(1 for k in KA if k['ok'])}/{len(KA)} pass; "
          f"peak {rss():.2f} GiB, {out['seconds']} s")
    if not all(k["ok"] for k in KA):
        raise SystemExit("KNOWN-ANSWER FAILURE")


if __name__ == "__main__":
    main()
