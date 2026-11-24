"""208 Phase 4 — the falsifiers, on the real challenge artifact.

FALSIFIER 1 (and the loudest thing in this experiment).  The previous two arms -- `rc_hybrid`
and `rc_callset` -- were PSEUDOBULK-PRESERVING: their per-gene panel totals were digest-
identical to a15's on 900/900 blocks, which is why the board reproduced `pds` and `mse` to
sixteen significant figures and why four of six members were pinned by construction.

**The Hamilton rounding rule is NOT in that class.**  It moves the delivered per-gene panel
totals themselves.  `pds`, `mse` and `nmae` are therefore NOT pinned for any arm carrying it,
and this experiment prices all six members without the exemption.  This script measures
exactly how far the totals move: how many blocks break the digest (expected 900/900), the
worst per-gene move (expected 1 count), and how many genes move per block.

The flat-1.00 spread IS in the pseudobulk-preserving class -- `emit_counts` computes `totals`
from `(n, L, profile)` before `spread` is read -- and `emit208.py` / `mse208.py` assert that
per emitted block rather than claiming it.

FALSIFIER 3.  nnz against the 4.75e9 cap (`rc_callset` sits at 33.22%) and the 1,000,000
counts/cell cap, both scaled off `mse208`'s measured per-block ratios on the real contexts.

    uv run python falsify208.py
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
WORK = HERE / "work"
WORK.mkdir(exist_ok=True)
E196 = ROOT / "experiments/196-emitted-depth"
sys.path.insert(0, str(E196))
import recon196  # noqa: E402
from ham196 import hamilton  # noqa: E402

PKG = ROOT / "submissions/202-rc_callset.h5ad"
PKG_SHA = "2a894aa543db7f574ebc024dbb3ebd5e36f923827c4e9a01f263092fe49692e0"
NNZ_CAP, CELL_CAP = 4_750_000_000, 1_000_000
RC_CALLSET_NNZ = 1_578_026_567          # 202 falsifier 6, the built package
n = 400


def main() -> None:
    t0 = time.time()
    out = {}
    print(f"artifact {PKG.name}  {PKG.stat().st_size:,} bytes")

    # ---- read the 900 delivered blocks' per-gene panel totals off the artifact ----------
    f = h5py.File(PKG, "r")
    tg = f["obs"]["target_gene"]
    cx = f["obs"]["context"]

    def decode(node):
        if isinstance(node, h5py.Group):          # categorical
            cats = np.asarray(node["categories"]).astype(str)
            return cats[np.asarray(node["codes"])]
        return np.asarray(node).astype(str)

    targets = decode(tg)
    contexts = decode(cx)
    indptr = f["X"]["indptr"][:]
    data, indices = f["X"]["data"], f["X"]["indices"]
    G = int(f["X"].attrs["shape"][1]) if "shape" in f["X"].attrs else 18533
    print(f"  {targets.size:,} cells, {G} genes, {int(indptr[-1]):,} stored entries "
          f"[{time.time()-t0:.0f}s]")

    keys, blocks = [], {}
    order = np.lexsort((targets, contexts))
    lab = np.char.add(np.char.add(contexts, "|"), targets)
    uniq, first = np.unique(lab[order], return_index=True)
    tot_art, nnz_art, maxcell_art = {}, {}, {}
    for k, u in enumerate(uniq):
        rows = order[first[k]:first[k] + n] if k + 1 == uniq.size else \
            order[first[k]:first[k + 1]]
        assert rows.size == n, (u, rows.size)
        rows = np.sort(rows)
        assert rows[-1] - rows[0] == n - 1, f"{u}: block is not contiguous"
        lo, hi = int(indptr[rows[0]]), int(indptr[rows[-1] + 1])
        col = indices[lo:hi].astype(np.int64)
        val = data[lo:hi].astype(np.float64)
        tot_art[u] = np.bincount(col, weights=val, minlength=G).astype(np.int64)
        nnz_art[u] = hi - lo
        ptr = np.asarray(indptr[rows[0]:rows[-1] + 2]) - lo
        maxcell_art[u] = int(max(val[ptr[r]:ptr[r + 1]].sum() for r in range(n)))
        if (k + 1) % 150 == 0:
            print(f"  {k+1}/{uniq.size} blocks read [{time.time()-t0:.0f}s]", flush=True)
    f.close()
    keys = list(uniq)
    blocks = tot_art

    # ---- gate A: the artifact's totals must be a15's, on 900/900 -----------------------
    z = np.load(E196 / "work/chal-totals-a15.npz", allow_pickle=True)
    a15_blocks = np.asarray(z["blocks"]).astype(str)
    a15_tot = z["totals"]
    a15_genes = np.asarray(z["genes"]).astype(str)
    idx = {b: i for i, b in enumerate(a15_blocks)}
    same, digests = 0, {}
    for u in keys:
        i = idx[u]
        d_art = hashlib.sha256(blocks[u].tobytes()).hexdigest()[:16]
        d_a15 = hashlib.sha256(a15_tot[i].astype(np.int64).tobytes()).hexdigest()[:16]
        digests[u] = (d_art, d_a15)
        same += d_art == d_a15
    print(f"\nGATE A — `rc_callset` artifact vs a15's recorded per-gene panel totals: "
          f"identical on {same}/{len(keys)} blocks")
    print(f"  first block {keys[0]}  digest {digests[keys[0]][0]}")
    out["gate_a_identical"] = same

    # ---- falsifier 1: Hamilton is NOT pseudobulk-preserving ----------------------------
    print("\nFALSIFIER 1 — the Hamilton rounding rule against the SAME 900 blocks")
    rowstat = []
    for ctx in "ABC":
        D = recon196.load_chal(ctx)
        assert np.array_equal(D["genes"], a15_genes)
        for i, t in enumerate(D["targets"]):
            u = f"{ctx}|{t}"
            R = recon196.Recon(D["sums"][i], D["r"], D["gate"], n=n, seed=1000 + int(i))
            assert np.array_equal(R.totals(1.0), blocks[u]), \
                f"{u}: the reconstruction does not re-deliver the artifact's own totals"
            Th = hamilton(R.x)
            dT = Th - blocks[u]
            rowstat.append(dict(
                block=u, digest_same=int(hashlib.sha256(Th.tobytes()).hexdigest()[:16]
                                         == digests[u][0]),
                max_abs=int(np.abs(dT).max()), moved=int((dT != 0).sum()),
                moved_gated=int((dT[D["gate"]] != 0).sum()),
                up=int((dT > 0).sum()), down=int((dT < 0).sum()),
                S_floor=int(blocks[u].sum()), S_ham=int(Th.sum())))
        print(f"  ctx {ctx} done [{time.time()-t0:.0f}s]", flush=True)
    ds = sum(r["digest_same"] for r in rowstat)
    print(f"\n  per-gene panel-total digest identical to `rc_callset` on {ds}/{len(rowstat)} "
          f"blocks   <-- NOT a pseudobulk-preserving arm")
    for k in ("max_abs", "moved", "moved_gated", "up", "down"):
        v = np.array([r[k] for r in rowstat])
        print(f"    {k:<12} min {v.min():>7,}  median {np.median(v):>9,.0f}  "
              f"max {v.max():>7,}")
    sf = np.array([r["S_floor"] for r in rowstat], float)
    sh = np.array([r["S_ham"] for r in rowstat], float)
    print(f"    panel total S  median {np.median(sf):,.0f} -> {np.median(sh):,.0f} "
          f"({100*np.median(sh/sf-1):+.4f}%)")
    out["falsifier1"] = dict(digest_same=ds, n_blocks=len(rowstat),
                             max_abs=int(max(r["max_abs"] for r in rowstat)),
                             moved_median=float(np.median([r["moved"] for r in rowstat])),
                             moved_gated_median=float(np.median(
                                 [r["moved_gated"] for r in rowstat])),
                             S_ratio_median=float(np.median(sh / sf)))

    # ---- falsifier 3: the two format caps ---------------------------------------------
    print("\nFALSIFIER 3 — the format contract")
    art_nnz = int(sum(nnz_art.values()))
    art_max = int(max(maxcell_art.values()))
    print(f"  artifact as built: nnz {art_nnz:,} = {100*art_nnz/NNZ_CAP:.2f}% of the "
          f"{NNZ_CAP:,.0f} cap   (202 recorded {RC_CALLSET_NNZ:,} = "
          f"{100*RC_CALLSET_NNZ/NNZ_CAP:.2f}%)")
    print(f"  artifact max counts/cell {art_max:,} = {100*art_max/CELL_CAP:.2f}% of the "
          f"{CELL_CAP:,.0f} cap")
    m = json.loads((WORK / "mse208.json").read_text()) if (WORK / "mse208.json").exists() \
        else None
    if m:
        base = m["rows"]["rc_callset|floor|tier"]
        print(f"\n  scaled off mse208's own emitted blocks on the real contexts "
              f"(MAXT={m['MAXT']} blocks x 3 contexts):")
        print(f"{'arm':<24}{'nnz/blk':>12}{'ratio':>9}{'pkg nnz':>17}{'% cap':>8}"
              f"{'max cell':>10}{'% cap':>8}")
        for k, r in m["rows"].items():
            if not k.startswith("rc_callset"):
                continue
            ratio = r["nnz"] / base["nnz"]
            pkg = art_nnz * ratio
            print(f"{k:<24}{r['nnz']:>12,.0f}{ratio:>9.4f}{pkg:>17,.0f}"
                  f"{100*pkg/NNZ_CAP:>8.2f}{r['max_cell']:>10,.0f}"
                  f"{100*r['max_cell']/CELL_CAP:>8.2f}")
            out.setdefault("format", {})[k] = dict(nnz_ratio=ratio, pkg_nnz=pkg,
                                                   pct_cap=100 * pkg / NNZ_CAP,
                                                   max_cell=r["max_cell"])
    out.update(artifact=dict(path=str(PKG), sha256_expected=PKG_SHA, nnz=art_nnz,
                             max_cell=art_max, blocks=len(keys)))
    json.dump(dict(out, rowstat=rowstat), open(WORK / "falsify208.json", "w"), default=float)
    print(f"\ndone [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
