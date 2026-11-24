"""458d — Arc's H1 programme from ALL of Arc's public H1 cells, not the 400-cell harness bank.

462c measured that the `h1_2025` bank's programme (255 perturbations x 400 cells vs 18,399 control rows, 15,571 unique)
is ~1/3 its control sample's noise by variance on expressed genes (shared-control split-half r 0.79, disjoint 0.39). The
raw VCC 2025 files hold ~1,100 cells per perturbation and 38,176 control cells. Summing every cell cuts both noises.

What it does (streaming, row chunks, < 3 GB):
  1. per perturbation and for the controls: summed counts over all cells, and over two seeded random cell halves;
  2. delta = log1p(5e4 p / sum p) - the same for the controls (the bank's frame): a denoised bank
     `data/donors/h1_2025_full/{delta.npy, perts.txt, genes.txt, meta.json}` for `def_sameline.py --programme-source`;
  3. reliability exactly as 462c (target halves vs the full control; vs disjoint control halves);
  4. against the harness bank: r of the programmes, and sign flips on each A/B/C block;
  5. the A/B/C artifact in artifact458's format (`work/arch1d_programme_abc.npz`) for a possible 458d arm.
Nothing is submitted or built here.   usage: uv run python experiments/458-arc-lab-programme/denoise458d.py
"""
from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "scripts"))
import def_sameline as S  # noqa: E402

S.BANK_DIRS.setdefault("h1_2025", ROOT / "data/donors/h1_2025")
FILES = [ROOT / "data/external/vcc2025/train/adata_Training.h5ad", ROOT / "data/external/vcc2025/validation/adata_Validation.h5ad",
         ROOT / "data/external/vcc2025/test/adata_Test.h5ad"]
COUNTS = [ROOT / "data/external/vcc2025/train/pert_counts_Training.csv",
          ROOT / "data/external/vcc2025/validation/pert_counts_Validation.csv",
          ROOT / "data/external/vcc2025/test/pert_counts_Test.csv"]
OUT_BANK = ROOT / "data/donors/h1_2025_full"
TS, CHUNK = 5e4, 4000
t0 = time.time()
log = lambda *x: print("[458d]", *x, f"({time.time() - t0:.0f}s)", flush=True)
rng = np.random.default_rng(4580)


def cat(n):
    return n["categories"][:].astype(str)[n["codes"][:]] if isinstance(n, h5py.Group) else n[:].astype(str)


def main():
    genes = None
    sums: dict[str, np.ndarray] = {}          # key -> (3, G): full, half1, half2
    ncell: dict[str, int] = {}
    seen_ctrl: set[str] = set()
    for path in FILES:
        with h5py.File(path, "r") as f:
            g = f["var/_index"]
            gv = (g[:] if not isinstance(g, h5py.Group) else g["values"][:]).astype(str)
            if genes is None:
                genes = gv
            assert np.array_equal(genes, gv), f"{path.name}: gene axis differs"
            tv = cat(f["obs/target_gene"])
            ix = f["obs/_index"]
            ids = (ix[:] if not isinstance(ix, h5py.Group) else ix["values"][:]).astype(str)
            # controls shared across files are counted once
            is_ctrl = tv == "non-targeting"
            dup = is_ctrl & np.isin(ids, list(seen_ctrl)) if seen_ctrl else np.zeros(tv.size, bool)
            seen_ctrl |= set(ids[is_ctrl])
            half = rng.integers(0, 2, tv.size)                     # a random half for every cell
            keys = np.where(is_ctrl, "__ctrl__", tv)
            ukeys, kidx = np.unique(keys, return_inverse=True)
            for k in ukeys:
                sums.setdefault(k, np.zeros((3, len(genes))))
            kk, kc = np.unique(keys[~dup], return_counts=True)
            for k, c in zip(kk, kc):
                ncell[k] = ncell.get(k, 0) + int(c)
            n = tv.size
            ip = f["X/indptr"]
            log(f"{path.name}: {n} cells, {int(is_ctrl.sum())} controls ({int(dup.sum())} already counted), "
                f"{len(ukeys) - 1} perturbations")
            for s in range(0, n, CHUNK):
                e = min(n, s + CHUNK)
                a, b = int(ip[s]), int(ip[e])
                M = sp.csr_matrix((f["X/data"][a:b], f["X/indices"][a:b], ip[s:e + 1] - a), shape=(e - s, len(genes)))
                keep = ~dup[s:e]
                for hv in (0, 1):
                    sel = keep & (half[s:e] == hv)
                    if not sel.any():
                        continue
                    G = sp.csr_matrix((np.ones(int(sel.sum())), (kidx[s:e][sel], np.flatnonzero(sel))),
                                      shape=(len(ukeys), e - s))
                    part = np.asarray((G @ M).todense())
                    for j in np.flatnonzero(np.asarray(G.sum(1)).ravel() > 0):
                        sums[ukeys[j]][0] += part[j]
                        sums[ukeys[j]][1 + hv] += part[j]
    log(f"summed: {len(sums) - 1} perturbations, controls {len(seen_ctrl)} unique cells")
    # KA (parse): every perturbation's cell count equals the official pert_counts CSV of its split
    want = {}
    for c in COUNTS:
        for r_ in csv.DictReader(open(c)):
            want[r_["target_gene"]] = int(r_["n_cells"])
    bad = {k: (ncell.get(k), v) for k, v in want.items() if ncell.get(k) != v}
    assert not bad and set(want) == set(sums) - {"__ctrl__"}, f"KA cell counts: {list(bad.items())[:5]}"
    log(f"KA: all {len(want)} perturbations' cell counts equal the official pert_counts CSVs")
    lp = lambda c: np.log1p(TS * c / c.sum())
    ctrl = sums.pop("__ctrl__")
    perts = sorted(sums)
    cF, c1, c2 = lp(ctrl[0]), lp(ctrl[1]), lp(ctrl[2])
    delta = np.stack([lp(sums[p][0]) - cF for p in perts]).astype(np.float32)
    # the denoised bank
    OUT_BANK.mkdir(parents=True, exist_ok=True)
    np.save(OUT_BANK / "delta.npy", delta)
    (OUT_BANK / "perts.txt").write_text("\n".join(perts) + "\n")
    (OUT_BANK / "genes.txt").write_text("\n".join(genes) + "\n")
    (OUT_BANK / "meta.json").write_text(json.dumps({
        "source": "VCC2025 H1 hESC, ALL cells of adata_Training + adata_Validation + adata_Test (458d)", "target_sum": TS,
        "n_perturbations": len(perts), "n_genes": int(len(genes)), "n_control_cells": len(seen_ctrl),
        "cells_per_pert": "all (no downsampling)", "builder": "experiments/458-arc-lab-programme/denoise458d.py"}, indent=1))
    log(f"wrote {OUT_BANK} delta {delta.shape} sha256 {S.file_sha256(OUT_BANK / 'delta.npy')[:16]}")
    # reliability, 462c's way, on the A/B/C-panel-excluded perturbations
    official = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
    panel = {r["target_gene"] for r in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))}
    keep = [i for i, p in enumerate(perts) if p not in panel]

    def pdn(M):
        nz = (M != 0).sum(0); dn = (M < 0).sum(0)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(nz > 0, dn / np.maximum(nz, 1), np.nan)
    T1C = np.stack([lp(sums[perts[i]][1]) - cF for i in keep]); T2C = np.stack([lp(sums[perts[i]][2]) - cF for i in keep])
    T1C1 = np.stack([lp(sums[perts[i]][1]) - c1 for i in keep]); T2C2 = np.stack([lp(sums[perts[i]][2]) - c2 for i in keep])
    cpm = 1e6 * ctrl[0] / ctrl[0].sum(); expressed = cpm >= 5

    def r(a, b, m):
        m = m & np.isfinite(a) & np.isfinite(b); return float(np.corrcoef(a[m] - .5, b[m] - .5)[0, 1])
    res = {"n_perts": len(perts), "n_perts_nonpanel": len(keep), "n_control_cells": len(seen_ctrl),
           "cells_per_pert_median": float(np.median([ncell[p] for p in perts])),
           "shared_control_r": r(pdn(T1C), pdn(T2C), expressed),
           "disjoint_control_r": r(pdn(T1C1), pdn(T2C2), expressed)}
    res["control_noise_share_of_reliability"] = 1 - res["disjoint_control_r"] / res["shared_control_r"]
    # vs the harness bank's programme, on the official axis and the A/B/C blocks
    psign_d, pd_d, n_d, _ = S.programme(delta, perts, list(genes), official, panel)
    d0, p0, g0, _ = S.load_bank("h1_2025"); psign_0, pd_0, n_0, _ = S.programme(d0, p0, g0, official, panel)
    m = np.isfinite(pd_d) & np.isfinite(pd_0)
    res["vs_bank"] = dict(rows_full=n_d, rows_bank=n_0, r=float(np.corrcoef(pd_d[m] - .5, pd_0[m] - .5)[0, 1]),
                          sign_flip_all=float(np.mean(psign_d[m] != psign_0[m])))
    blocks = {}
    c4 = np.load(ROOT / "experiments/391-consensus-block-signs/work/clean4_prior_challenge.npz", allow_pickle=True)
    c4s = np.asarray(c4["psign"], float)
    for c in "ABC":
        b = np.load(ROOT / f"experiments/438-october-stack-repricing/work/abc_frames/bctrl-{c}.npy"); e = np.expm1(b)
        order = np.argsort(-(e / e.sum())); blk = np.asarray([i for i in order if official[i] not in panel][:2200])
        blocks[c] = dict(flip_vs_bank=float(np.mean(psign_d[blk] != psign_0[blk])),
                         agree_clean4_bank=float(np.mean(psign_0[blk] == np.sign(c4s[blk]))),
                         agree_clean4_full=float(np.mean(psign_d[blk] == np.sign(c4s[blk]))),
                         frac_down_full=float(np.mean(psign_d[blk] < 0)))
    res["abc_blocks"] = blocks
    # the A/B/C artifact (artifact458's format, no rows, w 0)
    perts_abc = [r_["target_gene"] for r_ in csv.DictReader(open(ROOT / "data/challenge/pert_counts.csv"))]
    out = {"genes": np.asarray(official), "perts": np.asarray(perts_abc), "contexts": np.asarray(["A", "B", "C"])}
    for c in "ABC":
        out[f"programme_{c}"] = psign_d; out[f"programme_ok_{c}"] = np.asarray(True); out[f"w_{c}"] = np.asarray(0.0)
        out[f"pdown_{c}"] = pd_d; out[f"targets_{c}"] = np.zeros(0, dtype=np.int64)
        out[f"rows_{c}"] = np.zeros((0, len(official)), dtype=np.float32)
        out[f"signs_{c}"] = np.zeros((0, len(official)), dtype=np.int8)
    prov = {"experiment": "458d", "bank": "h1_2025_full", "bank_delta_sha256": S.file_sha256(OUT_BANK / "delta.npy"),
            "programme_rows": n_d, "programme_coverage": float(np.isfinite(pd_d).mean()),
            "programme_frac_down": float((psign_d < 0).mean()),
            "note": "Arc's H1 programme from all public H1 cells (458d); no rows"}
    out["provenance"] = np.asarray(json.dumps(prov))
    dst = HERE / "work/arch1d_programme_abc.npz"
    np.savez_compressed(dst, **out)
    res["artifact"] = dict(path=str(dst.relative_to(ROOT)), sha256=hashlib.sha256(dst.read_bytes()).hexdigest(), **prov)
    (HERE / "work/denoise458d.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
