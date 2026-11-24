"""271 -- item 1 inputs: each donor line's rows built from its two BATCH-PARITY halves, in the production row's units.

Reads 267's exact integer pseudobulk sums (`experiments/267-.../work/half/<line>/{sums.npz, acc_*.npy}`, read-only),
re-writes the batch-A / batch-B half-banks in the donor-table format under OUR work/half/<line>/ (267 deleted the
half `delta.npy` files after use; the arithmetic is truth267.write_half_bank verbatim), and prepares, on the FULL
300-target panel, through the instrument's own path (`E062.prepare_main_source`: half-bank column mean over every
perturbation of the half, x EFFECT_SCALE / raw_median):

    real   the production donor row (the real donor table; the K562 row here is UN-blended -- the fold's PM blend is
           applied downstream by lib263.stage, exactly as the cache was built)
    union  267's bank-union (== the donor table on the mapped columns, 267 KA-1), through the same path: the KNOWN-
           ANSWER row -- it must reproduce `real` to float32 precision (the fold cache's `pre` is the float32 mean)
    A, B   the half rows, converted to the real row's units: effects_half x (raw_median_half / raw_median_real)
           (so A, B and real share one scale and (A + B)/2 approximates real)

Writes work/half/rows-<line>.npz (float32 rows on the official axis, availability, per-target cell counts per half and
in full, the raw medians, the KA numbers) and work/half/rows-<line>.json; the half delta.npy files are removed after
use (reproducible from the sums in ~1 min).  One line per process:
    PYTHONDONTWRITEBYTECODE=1 uv run python halfrows271.py k562_gwps
"""
from __future__ import annotations

import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/271-donor-reliability-estimation"
WORK = HERE / "work"
HALF267 = ROOT / "experiments/267-replicate-validated-oracle/work/half"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
E062, OFFICIAL, PANEL = plolo.E062, plolo.OFFICIAL, plolo.PANEL
TS = 5e4
MIN_HALF = 5
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


# ------------------------------------------------------------------------------ 267's sums (truth267.Sums, verbatim)
class Sums:
    def __init__(self, line: str):
        self.line = line
        d = HALF267 / line
        z = np.load(d / "sums.npz", allow_pickle=True)
        self.keep = np.asarray(z["keep"], dtype=np.int64)
        self.keep_names = [str(x) for x in z["keep_names"]]
        self.pert_names = [str(x) for x in z["pert_names"]]
        self.tot = np.asarray(z["tot"], dtype=np.int64)
        self.ncell = np.asarray(z["ncell"], dtype=np.int64)
        self.ctrl_sum = np.asarray(z["ctrl_sum"], dtype=np.float64)
        self.ctrl_n = np.asarray(z["ctrl_n"], dtype=np.int64)
        self.acc = {h: np.load(d / f"acc_{h}.npy", mmap_mode="r") for h in ("batchA", "batchB")}
        self.k562 = line == "k562_gwps"
        self.on_all_genes = self.k562

    def half(self, side: str):
        A, B = self.acc["batchA"], self.acc["batchB"]
        cs, cn, tot, nc = self.ctrl_sum, self.ctrl_n, self.tot, self.ncell
        if side == "A":
            return (lambda lo, hi: np.asarray(A[lo:hi], dtype=np.float64)), tot[:, 0], nc[:, 0], cs[0], cn[0]
        return (lambda lo, hi: np.asarray(B[lo:hi], dtype=np.float64)), tot[:, 1], nc[:, 1], cs[1], cn[1]


def b_ctrl_of(ctrl_sum_all: np.ndarray) -> np.ndarray:
    return np.log1p(TS * ctrl_sum_all / ctrl_sum_all.sum()).astype(np.float32)


def write_half_bank(sm: Sums, side: str, out: Path) -> dict:
    """truth267.write_half_bank for the batch split (min 5 cells in the half), verbatim arithmetic."""
    S, tot, nc, cs, cn = sm.half(side)
    b_all = b_ctrl_of(cs)
    keep = sm.keep
    b_keep = b_all[keep]
    bank_perts = set((ROOT / f"data/donors/{sm.line}/perts.txt").read_text().split("\n"))
    rows = [i for i, nm in enumerate(sm.pert_names) if nm in bank_perts and nc[i] >= MIN_HALF and tot[i] > 0]
    rows = sorted(rows, key=lambda i: sm.pert_names[i])
    names = [sm.pert_names[i] for i in rows]
    out.mkdir(parents=True, exist_ok=True)
    delta = np.lib.format.open_memmap(out / "delta.npy", mode="w+", dtype=np.float32, shape=(len(rows), len(keep)))
    rows_a = np.asarray(rows, dtype=np.int64)
    step = 2048
    for lo in range(0, len(rows), step):
        r = rows_a[lo:lo + step]
        blk = S(0, len(sm.pert_names))[r] if len(sm.pert_names) <= 12000 else S(int(r.min()), int(r.max()) + 1)[r - int(r.min())]
        if sm.on_all_genes:
            blk = blk[:, keep]
        prof = np.log1p(TS * blk / tot[r][:, None]).astype(np.float32)
        delta[lo:lo + len(r)] = prof - b_keep
    delta.flush()
    del delta
    np.save(out / "b_ctrl.npy", b_keep.astype(np.float32))
    (out / "perts.txt").write_text("\n".join(names))
    (out / "genes.txt").write_text("\n".join(sm.keep_names))
    meta = dict(line=sm.line, split="batch", side=side, n_perts=len(names), min_cells=MIN_HALF, n_control_cells=int(cn),
                n_genes=len(keep), cells_per_pert=dict(min=int(nc[rows_a].min()), median=float(np.median(nc[rows_a])), max=int(nc[rows_a].max())))
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    return meta


def main(line: str) -> None:
    t0 = time.time()
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos], dtype=np.int64)
    targets = list(PANEL)
    P(f"[{line}] 300 panel targets; rss {rss():.2f} GiB")

    # the production row (un-blended), on the full panel
    real_bank = E062.Bank.load(line, plolo.DONOR_PATHS[line])
    real = E062.prepare_main_source(real_bank, targets, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=40.0)
    P(f"  real: {int(real.target_available.sum())}/300 targets, raw_median {real.raw_median_norm:.6f}, mapped {real.mapped_gene_count}, "
      f"duplicates {len(real.duplicate_symbols)}, rss {rss():.2f} GiB, {time.time() - t0:.0f}s")

    # KA: 267's union bank (== the donor table on the mapped columns) through the same path
    ub = HALF267 / line / "bank-union"
    union = E062.prepare_main_source(E062.Bank.load(f"{line}-union", ub), targets, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=40.0)
    assert np.array_equal(union.target_available, real.target_available)
    av = real.target_available
    d_eff = float(np.abs(union.effects[av] - real.effects[av]).max())
    rel = d_eff / float(np.abs(real.effects[av]).max())
    ratio_rm = union.raw_median_norm / real.raw_median_norm
    # direction agreement per row (what a scalar could hide) and the float32 rounding scale
    cosr = np.array([float(np.dot(union.effects[t], real.effects[t]) / (np.linalg.norm(union.effects[t]) * np.linalg.norm(real.effects[t]))) for t in np.flatnonzero(av)])
    f32 = float(np.abs(union.effects[av].astype(np.float32).astype(np.float64) - real.effects[av].astype(np.float32).astype(np.float64)).max())
    # the real bank's `raw_median` is the centred median row norm over ALL its columns that map onto the axis,
    # duplicates included (E062.centered_statistics on all_source_columns); the union half-bank carries the unique
    # columns only, so with duplicate symbols the two rows differ by ONE global scalar per donor -- the direction is
    # identical.  The known answer is therefore the row after that scalar, at float32 precision, plus the scalar.
    scaled = union.effects[av] * ratio_rm
    d_scaled = float(np.abs(scaled - real.effects[av]).max())
    ulp = float(np.max(np.abs(scaled.astype(np.float32).astype(np.float64) - real.effects[av].astype(np.float32).astype(np.float64))
                       / np.maximum(np.spacing(np.abs(real.effects[av]).astype(np.float32)).astype(np.float64), 1e-45)))
    ka = dict(label="KA-A union-of-halves bank through prepare_main_source == the production donor row (after the duplicate-column raw-median scalar)",
              max_abs_diff_unscaled=d_eff, rel_to_max_unscaled=rel, raw_median_ratio=float(ratio_rm), n_duplicate_symbols=len(real.duplicate_symbols),
              min_row_cos=float(cosr.min()), max_abs_diff_scaled=d_scaled, max_diff_float32_ulps=ulp,
              ok=bool(ulp <= 1.0 + 1e-9 and cosr.min() > 1 - 1e-12))
    P(f"  {ka['label']}: {'PASS' if ka['ok'] else 'FAIL'} unscaled max|d| {d_eff:.3e} (rel {rel:.1e}); raw_median ratio {ratio_rm:.9f} "
      f"({len(real.duplicate_symbols)} duplicate symbols); after the scalar max|d| {d_scaled:.2e} = {ulp:.2f} float32 ulp; min row cos {cosr.min():.9f}")

    # the halves
    sm = Sums(line)
    pert_idx = {nm: i for i, nm in enumerate(sm.pert_names)}
    out = dict(targets=np.asarray(targets), avail=av, real=real.effects.astype(np.float32), union=union.effects.astype(np.float32),
               raw_median_real=real.raw_median_norm, raw_median_union=union.raw_median_norm)
    metas = {}
    for side in ("A", "B"):
        bd = WORK / "half" / line / f"bank-batch{side}"
        meta = write_half_bank(sm, side, bd)
        hb = E062.Bank.load(f"{line}-batch{side}", bd)
        prep = E062.prepare_main_source(hb, targets, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=40.0)
        rows = prep.effects * (prep.raw_median_norm / real.raw_median_norm)      # the real row's units
        nc = sm.half(side)[2]
        ncells = np.asarray([nc[pert_idx[t]] if t in pert_idx else 0 for t in targets], dtype=np.int64)
        out[f"{side}"] = rows.astype(np.float32)
        out[f"avail_{side}"] = np.asarray(prep.target_available, bool)
        out[f"ncells_{side}"] = ncells
        out[f"raw_median_{side}"] = prep.raw_median_norm
        metas[side] = dict(meta, raw_median=prep.raw_median_norm, n_targets_available=int(prep.target_available.sum()))
        P(f"  batch{side}: {int(prep.target_available.sum())}/300 targets (>= {MIN_HALF} cells), half-bank {meta['n_perts']} perts, controls "
          f"{meta['n_control_cells']:,}, raw_median {prep.raw_median_norm:.6f} (x{prep.raw_median_norm / real.raw_median_norm:.4f} of real), "
          f"median cells/target {float(np.median(ncells[av])):.0f}, rss {rss():.2f} GiB, {time.time() - t0:.0f}s")
        (bd / "delta.npy").unlink()
    ncf = sm.ncell[:, 0] + sm.ncell[:, 1]
    out["ncells_full"] = np.asarray([ncf[pert_idx[t]] if t in pert_idx else 0 for t in targets], dtype=np.int64)
    # the reliability in the un-whitened residual space on the line's own mapped columns minus the panel (267's item-3
    # quantity for THIS line), for the known-answer check against 267's rel|batch|cos_resid
    A, B = out["A"].astype(np.float64), out["B"].astype(np.float64)
    ok = out["avail_A"] & out["avail_B"]
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (A * B).sum(1) / (na * nb)
    r = np.where(ok & (na > 0) & (nb > 0), r, np.nan)
    out["r_half_raw"] = r
    tr = np.load(HALF267 / f"truth-{line}.npz", allow_pickle=True)
    t267 = [str(t) for t in tr["targets"]]
    r267 = np.asarray(tr["rel|batch|cos_resid"], dtype=np.float64)
    pos = {t: i for i, t in enumerate(targets)}
    mine = np.array([r[pos[t]] for t in t267])
    d267 = float(np.nanmax(np.abs(mine - r267)))
    ka2 = dict(label="KA-B split-half reliability of the half rows == 267 rel|batch|cos_resid (float32 truths)", max_abs_diff=d267, n=len(t267),
               ok=bool(d267 < 1e-5))
    P(f"  {ka2['label']}: {'PASS' if ka2['ok'] else 'FAIL'} max|d| {d267:.2e} over {len(t267)} targets; median r {np.nanmedian(r[av]):+.4f}, "
      f"share <= 0 {float(np.mean(r[av] <= 0)):.3f}")
    # (A + B)/2 vs the real row: how close the half average is to the full-depth row (the Wiener filter's assumption)
    half_mean = 0.5 * (A + B)
    R = out["real"].astype(np.float64)
    with np.errstate(invalid="ignore", divide="ignore"):
        c_hm = (half_mean * R).sum(1) / (np.linalg.norm(half_mean, axis=1) * np.linalg.norm(R, axis=1))
    out["cos_halfmean_real"] = np.where(av, c_hm, np.nan)
    P(f"  cos((A+B)/2, real) median {np.nanmedian(c_hm[av]):.5f} min {np.nanmin(c_hm[av]):.5f}; ||(A+B)/2|| / ||real|| median "
      f"{np.nanmedian(np.linalg.norm(half_mean[av], axis=1) / np.linalg.norm(R[av], axis=1)):.4f}")
    np.savez_compressed(WORK / "half" / f"rows-{line}.npz", **out)
    js = dict(line=line, ka=[ka, ka2], halves=metas, n_available=int(av.sum()), raw_median_real=real.raw_median_norm,
              raw_median_union=union.raw_median_norm, seconds=time.time() - t0, peak_rss_gib=rss(),
              inputs=dict(sums=str((HALF267 / line / "sums.npz").relative_to(ROOT)), union_bank=str(ub.relative_to(ROOT)),
                          real_bank=str(plolo.DONOR_PATHS[line].relative_to(ROOT))))
    (WORK / "half" / f"rows-{line}.json").write_text(json.dumps(js, indent=1, default=str))
    (WORK / "half" / f"rows-{line}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{line}] done {(time.time() - t0) / 60:.1f} min, peak rss {rss():.2f} GiB")


if __name__ == "__main__":
    for line in sys.argv[1:]:
        main(line)
