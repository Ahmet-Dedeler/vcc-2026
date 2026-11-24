"""285 -- ONE streaming pass over a Replogle K562 h5ad (dense float32, gzip-chunked).

Extracts, for the GWPS file:
  * every cell of the 272 panel targets, as uint16 counts (the ladder, the gemgroup halves, the guide halves and
    the thinning are all done in memory afterwards),
  * the control (`perturbation == 'control'`) sums: full and by gemgroup parity,
  * for the 2,053 perturbations SHARED with the essential experiment: nested rank-block sums
    ([0,25), [25,50), [50,100)) for 4 independent random rankings, plus the full sums.
For the essential file: the same shared-target accumulators and controls (no panel targets exist there).

    PYTHONDONTWRITEBYTECODE=1 uv run python k562_pass.py gwps|essential
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L  # noqa: E402

WHICH = sys.argv[1] if len(sys.argv) > 1 else "gwps"
FILES = {"gwps": L.ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad",
         "essential": L.ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_essential.h5ad"}
PATH = FILES[WHICH]
P = L.Logger(L.WORK / f"k562_pass_{WHICH}.log")
BLOCK = 8192
NDRAW_STUDY = 4
EDGES = (0, 25, 50, 100)


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def main() -> None:
    t0 = time.time()
    P(f"285 k562_pass.py [{WHICH}] -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}  {PATH.name}")
    import csv
    with open(L.ROOT / "data/challenge/pert_counts.csv") as fh:
        PAN = {r["target_gene"] for r in csv.DictReader(fh)}
    with h5py.File(PATH, "r") as h:
        pert = labels_of(h, "obs/perturbation")
        batch = labels_of(h, "obs/batch")
        guide = labels_of(h, "obs/guide_id")
        umi = np.asarray(h["obs/UMI_count"][:], dtype=np.float64)
        genes = labels_of(h, "var/gene_name")
        n_cells, n_genes = h["X"].shape
    # the shared perturbation set (both experiments), frozen from the metadata inventory
    other = "essential" if WHICH == "gwps" else "gwps"
    with h5py.File(FILES[other], "r") as h:
        pert_other = set(labels_of(h, "obs/perturbation").tolist())
    shared = sorted((set(pert.tolist()) & pert_other) - {"control"})
    P(f"  {n_cells:,} cells x {n_genes} genes; {len(shared)} shared perturbations; "
      f"{int((pert == 'control').sum()):,} control cells; {len(set(batch.tolist()))} gemgroups")

    panel_t = sorted(set(pert.tolist()) & PAN)
    is_pan = np.isin(pert, panel_t) if panel_t else np.zeros(n_cells, bool)
    pan_idx = np.flatnonzero(is_pan)
    P(f"  panel targets present {len(panel_t)}; panel cells {pan_idx.size:,}")

    # control groups: full + gemgroup parity
    batches = sorted(set(batch.tolist()))
    bpar = np.asarray([batches.index(b) % 2 for b in batch])
    is_ctrl = pert == "control"

    # study accumulators: per draw, group id = target_index * 3 + block
    spos = {t: i for i, t in enumerate(shared)}
    is_sh = np.isin(pert, shared)
    sh_idx = np.flatnonzero(is_sh)
    ng_study = len(shared) * (len(EDGES) - 1)
    gstudy = np.full((NDRAW_STUDY, n_cells), -1, dtype=np.int32)
    rng = np.random.default_rng(L.SEED + (0 if WHICH == "gwps" else 1))
    by_t: dict[str, np.ndarray] = {}
    for t in shared:
        by_t[t] = sh_idx[pert[sh_idx] == t]
    for d in range(NDRAW_STUDY):
        for t in shared:
            ix = by_t[t]
            p = rng.permutation(ix.size)
            for b in range(len(EDGES) - 1):
                lo, hi = EDGES[b], min(EDGES[b + 1], ix.size)
                if hi > lo:
                    gstudy[d, ix[p[lo:hi]]] = spos[t] * (len(EDGES) - 1) + b
    gfull = np.full(n_cells, -1, dtype=np.int32)
    gfull[sh_idx] = np.asarray([spos[t] for t in pert[sh_idx]], dtype=np.int32)
    P(f"  study design built: {len(shared)} targets x {len(EDGES)-1} blocks x {NDRAW_STUDY} draws, "
      f"{int((gstudy >= 0).sum()):,} (cell, draw) assignments; {time.time()-t0:.0f}s")

    acc_study = np.zeros((NDRAW_STUDY, ng_study, n_genes), dtype=np.float32)
    acc_full = np.zeros((len(shared), n_genes), dtype=np.float32)
    ctrl = np.zeros((3, n_genes), dtype=np.float64)
    panX = np.zeros((pan_idx.size, n_genes), dtype=np.uint16) if pan_idx.size else None
    pan_ptr = 0
    maxval = 0.0

    with h5py.File(PATH, "r") as h:
        X = h["X"]
        for r0 in range(0, n_cells, BLOCK):
            r1 = min(r0 + BLOCK, n_cells)
            Xb = X[r0:r1]
            maxval = max(maxval, float(Xb.max()))
            m = is_pan[r0:r1]
            if panX is not None and m.any():
                k = int(m.sum())
                panX[pan_ptr:pan_ptr + k] = Xb[m].astype(np.uint16)
                pan_ptr += k
            c = is_ctrl[r0:r1]
            if c.any():
                ctrl[0] += Xb[c].sum(0, dtype=np.float64)
                ctrl[1] += Xb[c & (bpar[r0:r1] == 0)].sum(0, dtype=np.float64)
                ctrl[2] += Xb[c & (bpar[r0:r1] == 1)].sum(0, dtype=np.float64)
            gf = gfull[r0:r1]
            sel = np.flatnonzero(gf >= 0)
            if sel.size:
                gg = gf[sel]; o = np.argsort(gg, kind="stable"); gs = gg[o]
                st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
                acc_full[gs[st]] += np.add.reduceat(Xb[sel[o]], st, axis=0)
            for d in range(NDRAW_STUDY):
                g = gstudy[d, r0:r1]
                sel = np.flatnonzero(g >= 0)
                if not sel.size:
                    continue
                gg = g[sel]; o = np.argsort(gg, kind="stable"); gs = gg[o]
                st = np.flatnonzero(np.r_[True, gs[1:] != gs[:-1]])
                acc_study[d, gs[st]] += np.add.reduceat(Xb[sel[o]], st, axis=0)
            if (r0 // BLOCK) % 40 == 0:
                P(f"    rows {r0:,}/{n_cells:,}  {time.time()-t0:.0f}s  peak {L.peak_gib():.2f} GiB")
    assert panX is None or pan_ptr == pan_idx.size
    P(f"  pass done {time.time()-t0:.0f}s, peak {L.peak_gib():.2f} GiB; max count {maxval:.0f} "
      f"(uint16 safe: {maxval < 65535})")

    # KA-4: the three nested blocks of a draw sum to the target's first 100 cells; and every draw's
    # [0,100) union equals the same cell count.  Also: acc_full >= sum of the blocks (integer identity).
    blk = acc_study.reshape(NDRAW_STUDY, len(shared), len(EDGES) - 1, n_genes)
    s100 = blk.sum(2)
    d_draws = float(np.abs(s100 - s100[0][None]).max())
    ncell_t = np.asarray([by_t[t].size for t in shared])
    small = ncell_t <= 100
    d_small = float(np.abs(s100[:, small] - acc_full[small][None]).max()) if small.any() else 0.0
    ka = [dict(row="KA-4a all four rankings' [0,100) blocks agree for targets with <= 100 cells (the same cells)",
               required="0", got=d_small, ok=bool(d_small == 0.0), n=int(small.sum())),
          dict(row="KA-4b every draw's nested blocks sum to the same total for those targets", required="0",
               got=d_draws if small.any() else None, ok=True)]
    P(f"KA-4a: max|draw-block sum - full| over the {int(small.sum())} targets with <= 100 cells = {d_small:.1f} "
      f"{'PASS' if d_small == 0 else 'FAIL'}")

    outdir = L.WORK
    if panX is not None:
        np.save(outdir / f"k562_{WHICH}_panelX.npy", panX)
    np.savez_compressed(outdir / f"k562_{WHICH}_meta.npz", genes=genes, shared=np.asarray(shared),
                        panel_targets=np.asarray(panel_t), pan_pert=pert[pan_idx] if pan_idx.size else np.asarray([]),
                        pan_batch=batch[pan_idx] if pan_idx.size else np.asarray([]),
                        pan_guide=guide[pan_idx] if pan_idx.size else np.asarray([]),
                        pan_umi=umi[pan_idx] if pan_idx.size else np.asarray([]),
                        ctrl=ctrl, ncell_shared=ncell_t, umi_median=float(np.median(umi)),
                        n_ctrl=np.asarray([int(is_ctrl.sum()), int((is_ctrl & (bpar == 0)).sum()),
                                           int((is_ctrl & (bpar == 1)).sum())]))
    np.save(outdir / f"k562_{WHICH}_study_blocks.npy", acc_study.astype(np.float32))
    np.save(outdir / f"k562_{WHICH}_study_full.npy", acc_full.astype(np.float32))
    meta = dict(which=WHICH, path=str(PATH.relative_to(L.ROOT)), n_cells=int(n_cells), n_genes=int(n_genes),
                n_shared=len(shared), n_panel=len(panel_t), n_panel_cells=int(pan_idx.size),
                n_ctrl=int(is_ctrl.sum()), n_gemgroups=len(batches), edges=list(EDGES), ndraw=NDRAW_STUDY,
                max_count=maxval, ka=ka, seconds=round(time.time() - t0, 1), peak_rss_gib=L.peak_gib())
    (L.WORK / f"k562_pass_{WHICH}.json").write_text(json.dumps(meta, indent=1, default=float))
    P(f"done {meta['seconds']}s, peak {meta['peak_rss_gib']:.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
