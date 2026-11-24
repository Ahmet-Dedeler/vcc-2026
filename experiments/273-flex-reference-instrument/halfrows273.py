"""273 item 4(b) inputs -- each production donor's rows for h1-test's 85 targets and its two BATCH-PARITY half rows, in
the production row's units.  271's construction verbatim (267's exact integer pseudobulk sums -> half-banks in the
donor-table format via op272.write_half_bank == halfrows271.write_half_bank -> E062.prepare_main_source), the ONLY
difference being the target list (h1-test's 85 targets instead of the 300 panel) and the panel zeroing (none: on the H1
context the A/B/C panel genes are ordinary scored genes and the builder's `cand_h1` carries them raw).

KNOWN-ANSWER ROW: on the targets shared with the A/B/C panel that the donor carries, my `real`, `A`, `B` rows must equal
271's `work/half/rows-<line>.npz` rows bit for bit (float32) on every column outside the A/B/C panel positions, and the
raw medians must be identical -- `prepare_main_source`'s centring and scaling do not depend on the target list.

    PYTHONDONTWRITEBYTECODE=1 uv run python halfrows273.py <line>        (one line per process)
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib273 as X  # noqa: E402

plolo, E062, O = X.plolo, X.E062, X.O
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def main(line: str) -> None:
    t0 = time.time()
    st = X.staged170()
    targets = [str(t) for t in st["h1_targets"]]
    no_panel = np.zeros(0, dtype=np.int64)
    P(f"[{line}] {len(targets)} h1-test targets; no panel zeroing (H1 context); rss {X.peak_gib():.2f} GiB")
    real_bank = E062.Bank.load(line, plolo.DONOR_PATHS[line])
    real = E062.prepare_main_source(real_bank, targets, X.OFFICIAL, no_panel, chunk_rows=256, max_rss_gib=40.0)
    P(f"  real: {int(real.target_available.sum())}/{len(targets)} targets, raw_median {real.raw_median_norm:.6f}, mapped {real.mapped_gene_count}, "
      f"{len(real.duplicate_symbols)} duplicate symbols; {time.time() - t0:.0f}s, rss {X.peak_gib():.2f} GiB")
    sm = O.Sums(line)
    pert_idx = {nm: i for i, nm in enumerate(sm.pert_names)}
    out = dict(targets=np.asarray(targets), avail=np.asarray(real.target_available, bool), real=real.effects.astype(np.float64),
               raw_median_real=real.raw_median_norm, mean_on_axis=np.asarray(real.mean_on_axis, dtype=np.float64))
    metas = {}
    for side in ("A", "B"):
        bd = X.WORK / "half" / line / f"bank-batch{side}"
        meta = O.write_half_bank(sm, side, bd)
        hb = E062.Bank.load(f"{line}-batch{side}", bd)
        prep = E062.prepare_main_source(hb, targets, X.OFFICIAL, no_panel, chunk_rows=256, max_rss_gib=40.0)
        rows = prep.effects * (prep.raw_median_norm / real.raw_median_norm)
        nc = sm.half(side)[2]
        out[side] = rows.astype(np.float64)
        out[f"avail_{side}"] = np.asarray(prep.target_available, bool)
        out[f"ncells_{side}"] = np.asarray([nc[pert_idx[t]] if t in pert_idx else 0 for t in targets], dtype=np.int64)
        out[f"raw_median_{side}"] = prep.raw_median_norm
        metas[side] = dict(meta, raw_median=prep.raw_median_norm, n_targets_available=int(prep.target_available.sum()))
        P(f"  batch{side}: {int(prep.target_available.sum())}/{len(targets)} targets (>= {O.MIN_HALF} cells), half-bank {meta['n_perts']} perts, controls "
          f"{meta['n_control_cells']:,}, raw_median {prep.raw_median_norm:.6f} (x{prep.raw_median_norm / real.raw_median_norm:.4f} of real), median cells/target "
          f"{float(np.median(out[f'ncells_{side}'][out['avail']])):.0f}; {time.time() - t0:.0f}s, rss {X.peak_gib():.2f} GiB")
        (bd / "delta.npy").unlink()
    ncf = sm.ncell[:, 0] + sm.ncell[:, 1]
    out["ncells_full"] = np.asarray([ncf[pert_idx[t]] if t in pert_idx else 0 for t in targets], dtype=np.int64)

    # ---- KA against 271's rows on the targets shared with the A/B/C panel
    z = np.load(X.E271 / f"work/half/rows-{line}.npz", allow_pickle=True)
    t271 = [str(t) for t in z["targets"]]
    shared = [t for t in targets if t in t271 and real.target_available[targets.index(t)] and z["avail"][t271.index(t)]]
    keep = np.ones(len(X.OFFICIAL), dtype=bool)
    keep[X.PANEL_POS] = False
    worst = 0.0
    for key in ("real", "A", "B"):
        mine = out[key][[targets.index(t) for t in shared]].astype(np.float32)[:, keep]
        theirs = np.asarray(z[key])[[t271.index(t) for t in shared]][:, keep]
        worst = max(worst, float(np.abs(mine - theirs).max()) if shared else 0.0)
    d_rm = max(abs(out["raw_median_real"] - float(z["raw_median_real"])), abs(out["raw_median_A"] - float(z["raw_median_A"])),
               abs(out["raw_median_B"] - float(z["raw_median_B"])))
    ka = dict(label=f"KA-H rows for the {len(shared)} h1 targets shared with the A/B/C panel == 271's rows (float32, off-panel columns) + raw medians",
              shared_targets=shared, max_abs_diff=worst, max_abs_diff_raw_medians=d_rm, ok=bool(worst == 0.0 and d_rm == 0.0))
    P(f"  {ka['label']}: {'PASS' if ka['ok'] else 'FAIL'} max|d| {worst:.2e}, raw medians |d| {d_rm:.2e}")

    # ---- the reliability of the donor's rows on THIS target set (on the mapped columns minus the 85 target genes)
    ok = out["avail_A"] & out["avail_B"]
    mapped = np.asarray(real.mean_on_axis) != 0
    mapped[st["h1_positions"]] = False
    r = X.row_cos(out["A"][:, mapped], out["B"][:, mapped])
    r = np.where(ok, r, np.nan)
    out["r_half_raw"] = r
    P(f"  split-half reliability of the donor rows on the 85 h1 targets (batch halves, raw): median {np.nanmedian(r):+.4f}, share <= 0 "
      f"{float(np.nanmean(r <= 0)):.3f}, n {int(np.isfinite(r).sum())}; SB median {2 * np.nanmedian(r) / (1 + np.nanmedian(r)):+.4f}")
    np.savez_compressed(X.WORK / "half" / f"rows273-{line}.npz", **out)
    js = dict(line=line, ka=ka, halves=metas, n_available=int(out["avail"].sum()), raw_median_real=real.raw_median_norm,
              reliability=dict(median=float(np.nanmedian(r)), share_le0=float(np.nanmean(r <= 0)), n=int(np.isfinite(r).sum())),
              seconds=round(time.time() - t0, 1), peak_rss_gib=X.peak_gib(),
              inputs=dict(sums=str((O.HALF267 / line / "sums.npz").relative_to(X.ROOT)), real_bank=str(plolo.DONOR_PATHS[line].relative_to(X.ROOT)),
                          rows271=str((X.E271 / f"work/half/rows-{line}.npz").relative_to(X.ROOT))))
    (X.WORK / "half" / f"rows273-{line}.json").write_text(json.dumps(js, indent=1, default=str))
    (X.WORK / "half" / f"rows273-{line}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{line}] done {(time.time() - t0) / 60:.1f} min, peak rss {X.peak_gib():.2f} GiB")
    if not ka["ok"]:
        raise SystemExit(1)


if __name__ == "__main__":
    for ln in sys.argv[1:]:
        main(ln)
