"""299 -- INPUT STAGE (no fold number, no arm, no operator): the per-row split-half reliability and magnitude of
EVERY bank row of the three PLOLO donor lines, from 267's batch-parity half banks.

For a bank row t of line j (a perturbation the line's production donor table carries):

    A_j(t), B_j(t)  the row built from the line's batch-A / batch-B cells, each half bank centred on ITS OWN
                    column mean over its own rows (271 halfrows271 / 286 stage286, verbatim arithmetic), on the
                    official-axis columns BOTH halves map onto;
    r_j(t) = cos(A_j(t), B_j(t))                      the split-half reliability at half depth,
    R_j(t) = 2 r / (1 + r)                            Spearman-Brown to the full depth (271 est271.spearman_brown),
    m_j(t) = || (A_j(t) + B_j(t)) / 2 ||              the magnitude of the same object (the `mag` control's key).

The cosine is taken over ALL shared mapped columns (panel columns INCLUDED) because that is 271's definition --
271's rows are `E062.prepare_main_source` output, which fills every mapped position -- so the panel rows reproduce
271's reliabilities.  `r_nopanel` (the same cosine with the 300 panel columns dropped, the operator's own support)
is reported alongside as a sensitivity; the SELECTION uses `r`.

    PYTHONDONTWRITEBYTECODE=1 uv run python rel299.py k562_gwps orion_hct116 orion_hek293t
"""
from __future__ import annotations

import gc
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/299-reliable-training-rows"
WORK = HERE / "work"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
E272 = ROOT / "experiments/272-noise-corrected-operator"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E272))
sys.path.insert(0, str(E271))
import plolo  # noqa: E402
import op272 as O  # noqa: E402
import est271 as R  # noqa: E402

plolo.WORK = WORK / "plolo"
E062 = plolo.E062
OFFICIAL, PANEL, PANEL_SET = plolo.OFFICIAL, list(plolo.PANEL), plolo.PANEL_SET
BLOCK = 512
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_self() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


def main(line: str) -> None:
    t0 = time.time()
    (WORK / "rel").mkdir(parents=True, exist_ok=True)
    P(f"[{line}] rel299.py sha256 {sha256_self()[:16]}; op272 {O.sha256_self()[:16]}; est271 {R.sha256_self()[:16]}")
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos], dtype=np.int64)

    # ---------------------------------------------------------------- 267's half banks (op272.write_half_bank)
    sm = O.Sums(line)
    sides = {}
    for side in ("A", "B"):
        bd = WORK / "half" / line / f"bank-batch{side}"
        if not (bd / "delta.npy").exists():
            meta = O.write_half_bank(sm, side, bd)
        else:
            meta = json.loads((bd / "meta.json").read_text())
        hb = E062.Bank.load(f"{line}-batch{side}", bd)
        hmap = E062.map_bank_to_axis(hb, OFFICIAL)
        mean_all = E062.stream_column_mean(hb, hmap.all_source_columns, chunk_rows=2048, max_rss_gib=9.0)
        allpos = {int(s): i for i, s in enumerate(hmap.all_source_columns.tolist())}
        umean = np.asarray([mean_all[allpos[int(s)]] for s in hmap.unique_source_columns], dtype=np.float64)
        sides[side] = dict(bank=hb, meta=meta, umean=umean, src=hmap.unique_source_columns, pos=hmap.unique_target_positions)
        P(f"  half {side}: {meta['n_perts']} rows x {meta['n_genes']} genes, {hmap.unique_target_positions.size} mapped "
          f"official columns, median {meta['cells_per_pert']['median']:.0f} cells/row, {time.time()-t0:.0f}s, rss {rss():.2f} GiB")
    del sm
    gc.collect()

    cols = np.intersect1d(sides["A"]["pos"], sides["B"]["pos"])                 # official positions, panel INCLUDED
    keep_np = ~np.isin(cols, panel_positions)                                   # the no-panel sensitivity mask
    for side in ("A", "B"):
        take = {int(p): i for i, p in enumerate(sides[side]["pos"].tolist())}
        sel = np.asarray([take[int(p)] for p in cols], dtype=np.int64)
        sides[side]["sel"] = sel
        sides[side]["s_src"] = sides[side]["src"][sel]
        sides[side]["s_mu"] = sides[side]["umean"][sel]
    P(f"  shared mapped columns {cols.size} ({int((~keep_np).sum())} of them panel columns)")

    names = [n for n in sides["A"]["bank"].perts if n in sides["B"]["bank"].pert_pos]
    names.sort()
    nrows = len(names)
    rowA = np.asarray([sides["A"]["bank"].pert_pos[n] for n in names], dtype=np.int64)
    rowB = np.asarray([sides["B"]["bank"].pert_pos[n] for n in names], dtype=np.int64)
    is_panel = np.asarray([n in PANEL_SET for n in names], dtype=bool)
    P(f"  rows in both halves: {nrows} ({int(is_panel.sum())} panel targets)")

    delt = {s: np.load(sides[s]["bank"].delta_path, mmap_mode="r") for s in ("A", "B")}
    r = np.full(nrows, np.nan)
    r_np = np.full(nrows, np.nan)
    mag = np.full(nrows, np.nan)
    nrmA = np.full(nrows, np.nan)
    nrmB = np.full(nrows, np.nan)
    for lo in range(0, nrows, BLOCK):
        hi = min(lo + BLOCK, nrows)
        blk = {}
        for s, rr in (("A", rowA), ("B", rowB)):
            idx = rr[lo:hi]
            order = np.argsort(idx)
            raw = np.asarray(delt[s][np.sort(idx)][:, sides[s]["s_src"]], dtype=np.float64)
            out = np.empty_like(raw)
            out[order] = raw
            blk[s] = out - sides[s]["s_mu"]
            del raw, out
        a, b = blk["A"], blk["B"]
        r[lo:hi] = O.row_cos(a, b)
        r_np[lo:hi] = O.row_cos(a[:, keep_np], b[:, keep_np])
        h = 0.5 * (a + b)
        mag[lo:hi] = np.linalg.norm(h, axis=1)
        nrmA[lo:hi] = np.linalg.norm(a, axis=1)
        nrmB[lo:hi] = np.linalg.norm(b, axis=1)
        del blk, a, b, h
    Rsb = R.spearman_brown(r)
    P(f"  reliabilities: median r {np.nanmedian(r):+.4f}, share r<=0 {np.nanmean(r <= 0):.3f}, "
      f"median R(SB) {np.nanmedian(Rsb):+.4f}, undefined {int((~np.isfinite(r)).sum())}; {time.time()-t0:.0f}s")

    # ---------------------------------------------------------------- KA against 271's panel-row reliabilities
    z271 = np.load(E271 / f"work/half/rows-{line}.npz", allow_pickle=True)
    t271 = [str(t) for t in z271["targets"]]
    ok271 = np.asarray(z271["avail_A"], bool) & np.asarray(z271["avail_B"], bool)
    r271_full = O.row_cos(z271["A"].astype(np.float64), z271["B"].astype(np.float64))          # 271's own definition
    r271_cols = O.row_cos(z271["A"].astype(np.float64)[:, cols], z271["B"].astype(np.float64)[:, cols])
    pos_of = {n: i for i, n in enumerate(names)}
    shared = [(pos_of[n], t271.index(n)) for n in names if n in t271 and ok271[t271.index(n)]]
    mine = np.asarray([r[i] for i, _ in shared])
    theirs_full = np.asarray([r271_full[k] for _, k in shared])
    theirs_cols = np.asarray([r271_cols[k] for _, k in shared])
    d_full = float(np.nanmax(np.abs(mine - theirs_full)))
    d_cols = float(np.nanmax(np.abs(mine - theirs_cols)))
    P(f"  KA-R vs 271 on {len(shared)} shared panel rows: max|d| vs 271's full-axis row_cos {d_full:.3e}; "
      f"vs 271's rows restricted to `cols` {d_cols:.3e}")

    np.savez_compressed(WORK / "rel" / f"{line}.npz", names=np.asarray(names), r=r, r_sb=Rsb, r_nopanel=r_np,
                        mag=mag, norm_A=nrmA, norm_B=nrmB, is_panel=is_panel, cols=cols)
    meta = dict(line=line, rel299_sha256=sha256_self(), op272_sha256=O.sha256_self(), est271_sha256=R.sha256_self(),
                n_rows=nrows, n_panel_rows=int(is_panel.sum()), n_cols=int(cols.size), n_panel_cols=int((~keep_np).sum()),
                halves={s: sides[s]["meta"] for s in ("A", "B")},
                r=dict(median=float(np.nanmedian(r)), mean=float(np.nanmean(r)), share_le0=float(np.nanmean(r <= 0)),
                       p10=float(np.nanpercentile(r, 10)), p90=float(np.nanpercentile(r, 90)), n_undefined=int((~np.isfinite(r)).sum())),
                r_nonpanel_rows=dict(median=float(np.nanmedian(r[~is_panel])), n=int((~is_panel).sum())),
                spearman_r_mag=float(np.corrcoef(np.argsort(np.argsort(np.nan_to_num(r, nan=-9))),
                                                 np.argsort(np.argsort(np.nan_to_num(mag, nan=-9))))[0, 1]),
                ka=dict(shared_panel_rows=len(shared), max_abs_diff_vs_271_full_axis=d_full,
                        max_abs_diff_vs_271_on_cols=d_cols),
                seconds=round(time.time() - t0, 1), peak_rss_gib=rss())
    (WORK / "rel" / f"{line}.json").write_text(json.dumps(meta, indent=1, default=str))
    (WORK / "rel" / f"{line}.log").write_text("\n".join(LOG) + "\n")
    for s in ("A", "B"):
        (sides[s]["bank"].path / "delta.npy").unlink(missing_ok=True)     # reproducible from 267's sums in ~1 min
    P(f"[{line}] done {(time.time()-t0)/60:.1f} min, peak rss {rss():.2f} GiB")


if __name__ == "__main__":
    for ln in sys.argv[1:]:
        main(ln)
