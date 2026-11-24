#!/usr/bin/env python3
"""517 — Arc programme ESTIMATOR variants as the hook source, priced on H1 (price517's validated analytic path).
progest517 (Xaira, cross-line within a lab) found that for a small, Arc-like programme (275 rows, 36 % essential), the
essential rows do not transfer to another line's non-essential targets' uncalled genes (0.48–0.50), the non-essential rows
alone are as good or better than the mix, and a mean z-scored delta beats the P_down sign count (+0.01 … +0.03). H1 is Arc's
own line, so it cannot show the cross-line benefit; this checks that the variants cost nothing on H1 (same line).
Variants (h1-test's rows out, 467's rule): pdown_all (production, KA == price517's Arc), pdown_ne (non-essential rows),
meanz_all, meanz_ne (mean over rows of delta / per-gene sd). Hook on the c4 base, gated (le1 on H1 = 467's strong 57;
weak -> CLEAN4 tau .05), at Arc tau .05's coverage (63 % of gated genes, matched per variant by |z| quantile).
usage: uv run python experiments/517-flex-programme/variant517.py      (< 2 GB)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import price517 as P  # noqa: E402
RC = P.RC
ESS = {l.strip() for l in open(ROOT / "data/donors/k562_essential/perts.txt") if l.strip()}


def main():
    zr = np.load(RC.REF, allow_pickle=False)
    genes_full = zr["genes"]; keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str)
    zr.close()
    h1t = set(T)
    d = ROOT / "data/donors/h1_2025_full"
    delta = np.load(d / "delta.npy"); perts = [l.strip() for l in open(d / "perts.txt")]
    bg = [l.strip() for l in open(d / "genes.txt")]
    pos = {g: j for j, g in enumerate(bg)}
    col = np.array([pos.get(str(g), -1) for g in np.asarray(genes_full)[keep]])
    assert (col >= 0).all()
    rows = np.array([i for i, p in enumerate(perts) if p not in h1t])
    ne = np.array([perts[i] not in ESS for i in rows])
    X = np.asarray(delta[rows][:, col], dtype=np.float64)
    print(f"Arc rows {rows.size} (h1-test's out): non-essential {int(ne.sum())}, essential {int((~ne).sum())}", flush=True)

    def pdown(M):
        nz = (M != 0).sum(0); dn = (M < 0).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    sd = X.std(0); sd[sd == 0] = 1
    Z = {"pdown_all": 0.5 - pdown(X), "pdown_ne": 0.5 - pdown(X[ne]),
         "meanz_all": (X / sd).mean(0), "meanz_ne": (X[ne] / sd).mean(0)}
    SRC, _ = P.load_sources(genes_full, keep)
    assert np.allclose(Z["pdown_all"], 0.5 - SRC["arc"], equal_nan=True), "KA: pdown_all != price517's Arc"
    print("KA: pdown_all == price517's Arc programme", flush=True)
    n_cov = int((np.isfinite(SRC["arc"]) & (np.abs(SRC["arc"] - 0.5) >= 0.05)).sum())
    c4 = RC.prior_signs_full(genes_full)[keep]
    hc = np.where(np.isfinite(SRC["c4"]) & (np.abs(SRC["c4"] - 0.5) >= 0.05), np.where(SRC["c4"] > 0.5, -1.0, 1.0), 0.0)
    strong, _ = P.strong_set(T)
    out = {"n_cov": n_cov}
    for seed in (0, 1, 2):
        ctx = P.Ctx("c4", seed, lambda i, t: c4)
        rng = np.random.default_rng(5176 + seed)
        for name, z in Z.items():
            a = np.abs(np.nan_to_num(z))
            thr = np.sort(a)[::-1][n_cov - 1]
            hs = np.where(a >= thr, np.sign(z), 0.0)
            for gated in (True, False):
                fn = (lambda i, t, hs=hs: hs if strong[i] else hc) if gated else (lambda i, t, hs=hs: hs)
                lfc, rec = ctx.hook(fn)
                sc, _ = ctx.score(lfc, rng, n_boot=500)
                acc = ctx.gate_accuracy(fn, "comp")
                key = f"{name}_{'g' if gated else 'u'}"
                out.setdefault(key, []).append({**{m: sc[m] for m in P.SP}, "ov": sc["overall_equiv"], "ci95": sc["ci95"],
                                                "ne_ov": sc["ne_overall_equiv"], "ne_nmae": sc["ne_nmae"], "acc": acc["acc"]})
                print(f"s{seed} {key:14s} nmae {sc['nmae']:+.5f} reach {sc['reach']:+.5f} ov {sc['overall_equiv']:+.5f} "
                      f"ne_nmae {sc['ne_nmae']:+.5f} | complement-gate acc {acc['acc']:.4f} (n {acc['n']})", flush=True)
        del ctx
    (HERE / "work/variant517.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
