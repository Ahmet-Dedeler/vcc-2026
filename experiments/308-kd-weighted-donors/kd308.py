"""308 -- THE FEATURE: on-target knockdown |kd|(t, d) of every panel target t in every donor screen d.

Provenance, stated so it can be checked (this is the leakage argument):
  * the ONLY input is `data/donors/<d>/{delta.npy, genes.txt, perts.txt}` -- donor screen d's own pseudobulk
    perturbation table, `delta[p, g] = log1p(TS * pert_frac_g) - log1p(TS * ctrl_frac_g)`, TS = 50,000
    (`scripts/build_donor_table.py:10`, `scripts/build_orion_donor.py:112-116`), i.e. exactly 293's `|kd|` scale;
  * kd(t, d) = delta[row of t, column of gene t] -- t's OWN gene coordinate in d's OWN measured response to
    perturbing t, in d's own units.  No held line, no truth, no fold, no operator, no panel-mean, no centring
    touches it: it is a function of ONE donor bank and the target symbol.
  * the production prepared rows ZERO every panel column (`E062.prepare_main_source`: `effects[:, panel_positions] =
    0.0`), so this coordinate is NOT in the incumbent's base and cannot be recovered from it -- it has to be read
    from the bank, which is what this script does.

Outputs: work/kd308_feature.npz, work/kd308_feature.json, work/kd308.log
    PYTHONDONTWRITEBYTECODE=1 uv run python kd308.py
"""
from __future__ import annotations

import csv
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

E062 = plolo.E062
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
PANEL, OFFICIAL = plolo.PANEL, plolo.OFFICIAL
CHUNK = 512
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def own_column_map(bank) -> dict[str, int]:
    """gene symbol -> its LAST column in the bank (E062.map_bank_to_axis's `last_source` convention)."""
    pos: dict[str, int] = {}
    for col, g in enumerate(bank.genes):
        pos[g] = col
    return pos


def read_diagonal(bank, pairs: list[tuple[int, int]]) -> np.ndarray:
    """delta[r, c] for a list of (row, col), read in row-sorted chunks off the memmap."""
    delta = np.load(bank.delta_path, mmap_mode="r")
    out = np.full(len(pairs), np.nan)
    order = np.argsort([r for r, _ in pairs], kind="stable")
    for lo in range(0, len(order), CHUNK):
        idx = order[lo:lo + CHUNK]
        rows = np.array([pairs[i][0] for i in idx])
        cols = np.array([pairs[i][1] for i in idx])
        r0, r1 = int(rows.min()), int(rows.max()) + 1
        block = np.asarray(delta[r0:r1], dtype=np.float64)
        out[idx] = block[rows - r0, cols]
        del block
    return out


def main() -> None:
    t0 = time.time()
    P(f"308 feature: on-target knockdown per (target, donor).  panel {len(PANEL)}, official axis {len(OFFICIAL)}")
    out: dict = dict(donors={}, panel=list(PANEL))
    kd = np.full((len(PANEL), len(DONORS)), np.nan)          # raw delta at the own-gene column
    measured = np.zeros((len(PANEL), len(DONORS)), bool)      # gene column exists AND the target is perturbed
    perturbed = np.zeros((len(PANEL), len(DONORS)), bool)     # the target is a perturbation in the bank
    on_axis = np.zeros((len(PANEL), len(DONORS)), bool)       # the target gene has a column in the bank

    for j, name in enumerate(DONORS):
        bank = E062.Bank.load(name, plolo.DONOR_PATHS[name])
        gpos = own_column_map(bank)
        meta = json.loads((plolo.DONOR_PATHS[name] / "meta.json").read_text())
        # --- the panel targets
        pairs, rows_idx = [], []
        for i, t in enumerate(PANEL):
            r, c = bank.pert_pos.get(t), gpos.get(t)
            perturbed[i, j] = r is not None
            on_axis[i, j] = c is not None
            if r is not None and c is not None:
                pairs.append((int(r), int(c)))
                rows_idx.append(i)
        vals = read_diagonal(bank, pairs)
        kd[np.asarray(rows_idx, dtype=np.int64), j] = vals
        measured[np.asarray(rows_idx, dtype=np.int64), j] = True
        # --- the library-wide distribution: every perturbation whose own gene is on the bank's axis
        lib_pairs = []
        for p, pert in enumerate(bank.perts):
            c = gpos.get(pert)
            if c is not None:
                lib_pairs.append((int(p), int(c)))
        lib = read_diagonal(bank, lib_pairs)
        lib_abs = np.abs(lib)
        med_lib = float(np.median(lib_abs))
        pk = kd[:, j]
        P(f"[{name}] bank {bank.shape} perts {len(bank.perts):,} genes {len(bank.genes):,}; "
          f"panel targets perturbed {int(perturbed[:, j].sum())}/300, own gene on the bank axis "
          f"{int(on_axis[:, j].sum())}/300, BOTH (feature defined) {int(measured[:, j].sum())}/300")
        P(f"[{name}] library-wide |kd| over {len(lib_pairs):,} self-measured perturbations: median {med_lib:.4f}, "
          f"q10 {np.quantile(lib_abs, .10):.4f} q25 {np.quantile(lib_abs, .25):.4f} q75 {np.quantile(lib_abs, .75):.4f} "
          f"q90 {np.quantile(lib_abs, .90):.4f}; |kd| < 1.0: {float((lib_abs < 1.0).mean())*100:.1f} %, "
          f"< 0.5: {float((lib_abs < 0.5).mean())*100:.1f} %, positive lfc {float((lib > 0).mean())*100:.1f} %")
        good = np.isfinite(pk)
        P(f"[{name}] PANEL |kd|: median {np.median(np.abs(pk[good])):.4f}, "
          f"q10 {np.quantile(np.abs(pk[good]), .10):.4f} q90 {np.quantile(np.abs(pk[good]), .90):.4f}; "
          f"< 1.0: {float((np.abs(pk[good]) < 1.0).mean())*100:.1f} %, < 0.5: {float((np.abs(pk[good]) < 0.5).mean())*100:.1f} %, "
          f"lfc > 0 (no repression at all): {float((pk[good] > 0).mean())*100:.1f} %")
        out["donors"][name] = dict(
            n_perts=len(bank.perts), n_genes=len(bank.genes), meta_cells_per_pert=meta.get("cells_per_pert"),
            n_self_measured=len(lib_pairs), lib_median_abs=med_lib,
            lib_q=dict(q10=float(np.quantile(lib_abs, .10)), q25=float(np.quantile(lib_abs, .25)),
                       q75=float(np.quantile(lib_abs, .75)), q90=float(np.quantile(lib_abs, .90))),
            lib_frac_below_1=float((lib_abs < 1.0).mean()), lib_frac_below_05=float((lib_abs < 0.5).mean()),
            lib_frac_positive=float((lib > 0).mean()),
            panel_perturbed=int(perturbed[:, j].sum()), panel_on_axis=int(on_axis[:, j].sum()),
            panel_measured=int(measured[:, j].sum()),
            panel_median_abs=float(np.median(np.abs(pk[good]))),
            panel_frac_below_1=float((np.abs(pk[good]) < 1.0).mean()),
            panel_frac_below_05=float((np.abs(pk[good]) < 0.5).mean()),
            panel_frac_positive=float((pk[good] > 0).mean()))
        del bank

    # ------------------------------------------------------------------ cuts and the disagreeing set
    A = np.abs(kd)
    cut_abs = 1.0                                              # 293's absolute cut
    cut_rel = np.array([out["donors"][d]["lib_median_abs"] for d in DONORS])   # 296's C2, per library
    P("")
    P(f"cuts: absolute (293) |kd| >= {cut_abs}; library-relative (296 C2) per donor = "
      + ", ".join(f"{d} {c:.4f}" for d, c in zip(DONORS, cut_rel)))
    for label, cut in (("abs1.0", np.full(3, cut_abs)), ("rel", cut_rel)):
        works = (A >= cut[None, :]) & measured
        fails = (A < cut[None, :]) & measured
        nboth = measured.sum(1)
        n_work = works.sum(1)
        n_fail = fails.sum(1)
        dis = (n_work >= 1) & (n_fail >= 1)                     # at least one working and one failing donor
        allfail = (nboth >= 1) & (n_work == 0)
        P(f"[{label}] per-donor pass rate: " + ", ".join(
            f"{d} {int(works[:, j].sum())}/{int(measured[:, j].sum())} ({100*works[:, j].sum()/max(measured[:, j].sum(),1):.0f} %)"
            for j, d in enumerate(DONORS)))
        P(f"[{label}] targets with >= 2 donors measured: {int((nboth >= 2).sum())}; "
          f"DISAGREEING (>= 1 works, >= 1 fails, among measured): {int(dis.sum())} of 300; "
          f"all measured donors FAIL: {int(allfail.sum())}; all measured donors WORK: "
          f"{int(((nboth >= 1) & (n_fail == 0)).sum())}")
        out[f"cut_{label}"] = dict(cut=[float(c) for c in cut], n_disagree=int(dis.sum()),
                                   n_allfail=int(allfail.sum()), n_allwork=int(((nboth >= 1) & (n_fail == 0)).sum()),
                                   per_donor_pass={d: int(works[:, j].sum()) for j, d in enumerate(DONORS)})

    # ------------------------------------------------------------------ depth covariate (278, read-only)
    cells = np.full((len(PANEL), len(DONORS)), np.nan)
    with (ROOT / "experiments/278-depth-inventory/work/depth278.csv").open() as fh:
        for row in csv.DictReader(fh):
            if row["source"] in DONORS and row["target"] in set(PANEL):
                i, j = PANEL.index(row["target"]), DONORS.index(row["source"])
                try:
                    cells[i, j] = float(row["cells"])
                except (TypeError, ValueError):
                    pass
    for j, d in enumerate(DONORS):
        m = np.isfinite(cells[:, j]) & measured[:, j]
        if m.sum():
            P(f"[{d}] depth covariate (278 cells/target) on {int(m.sum())} measured targets: median "
              f"{np.median(cells[m, j]):.0f}, q10 {np.quantile(cells[m, j], .1):.0f}, q90 {np.quantile(cells[m, j], .9):.0f}; "
              f"Spearman(cells, |kd|) = {spearman(cells[m, j], A[m, j]):+.3f}")

    np.savez_compressed(WORK / "kd308_feature.npz", panel=np.asarray(PANEL), donors=np.asarray(DONORS),
                        kd=kd, abs_kd=A, measured=measured, perturbed=perturbed, on_axis=on_axis,
                        cells=cells, cut_rel=cut_rel, cut_abs=np.float64(cut_abs))
    out["peak_gib"] = peak_gib()
    (WORK / "kd308_feature.json").write_text(json.dumps(out, indent=2))
    P(f"\nwrote work/kd308_feature.npz + .json in {time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB")
    (WORK / "kd308.log").write_text("\n".join(LINES) + "\n")


def spearman(a, b) -> float:
    from scipy.stats import spearmanr
    return float(spearmanr(a, b).statistic)


if __name__ == "__main__":
    main()
