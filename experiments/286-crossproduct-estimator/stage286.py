"""286 -- Stage B inputs, ONE DONOR LINE PER PROCESS: the cross-product moments S-hat, N-hat from 267's batch-parity
half banks, fitted on ALL the line's bank rows EXCEPT the 300 panel targets (held targets never enter the fit).

    PYTHONDONTWRITEBYTECODE=1 uv run python stage286.py k562_gwps

Views.  267's half banks are rebuilt from its exact integer pseudobulk sums with `op272.write_half_bank` (which is
`halfrows271.write_half_bank`, which is `truth267.write_half_bank`, verbatim), then every bank row is converted to the
PRODUCTION row's units exactly as `E062.prepare_main_source` would:

    r_i = (delta_half[i, unique_source] - mean_half[unique_source]) * EFFECT_SCALE / raw_median_REAL

(the half bank's own column mean over ITS rows -- "each view centred on its own mean over the bank's rows" -- and the
REAL bank's raw median, which is halfrows271's `prep.effects * (raw_median_half / raw_median_real)` written out).
Panel COLUMNS are dropped (production zeroes them); panel TARGET ROWS are dropped (they are the held targets).

Outputs (work/moments/<line>.npz + .json):
    S, N          float32 (m, m) on `cols` (official-axis positions, panel columns excluded)
    S_shuf, N_shuf the Cs control: view a of row t paired with view b of row pi(t), a fixed derangement
    cols, n_rows, the KA-M moment identities, the leading-eigenvalue ratio input, timings, sha256 of est286.py
Also work/prod/<line>.npz: the line's production un-whitened rows for ALL 300 panel targets (E062.prepare_main_source
on the real bank), the availability and the raw median -- the object every arm filters.
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
HERE = ROOT / "experiments/286-crossproduct-estimator"
WORK = HERE / "work"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
E272 = ROOT / "experiments/272-noise-corrected-operator"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(E272))
sys.path.insert(0, str(HERE))
import plolo  # noqa: E402
import op272 as O  # noqa: E402
import est286 as E  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)
E062 = plolo.E062
OFFICIAL, PANEL, PANEL_SET = plolo.OFFICIAL, list(plolo.PANEL), plolo.PANEL_SET
EFFECT_SCALE = plolo.EFFECT_SCALE
BLOCK = 512
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def moments_inplace(G_U: np.ndarray, G_D: np.ndarray, n: int, *, block: int = 1024):
    """est286.moments + est286.moment_ka, in place: S overwrites G_U, N overwrites G_D, and every Frobenius norm is
    accumulated over row blocks so C = S + 2N is never materialised.  Algebraically identical (KA-Mi)."""
    m = G_U.shape[0]
    fC2 = fS2 = fN2 = fEx2 = fLi2 = 0.0
    trC = trS = trN = 0.0
    for lo in range(0, m, block):
        hi = min(lo + block, m)
        s = (G_U[lo:hi] - 0.25 * G_D[lo:hi]) / n
        nn = G_D[lo:hi] / (4.0 * n)
        c = s + 2.0 * nn
        fC2 += float((c * c).sum()); fS2 += float((s * s).sum()); fN2 += float((nn * nn).sum())
        fEx2 += float(((c - s - 2.0 * nn) ** 2).sum()); fLi2 += float(((c - s - nn) ** 2).sum())
        d = np.arange(lo, hi)
        trC += float(c[np.arange(hi - lo), d].sum()); trS += float(s[np.arange(hi - lo), d].sum())
        trN += float(nn[np.arange(hi - lo), d].sum())
        G_U[lo:hi] = s
        G_D[lo:hi] = nn
        del s, nn, c
    fC = fC2 ** 0.5
    ka = dict(frob_C=fC, frob_S=fS2 ** 0.5, frob_N=fN2 ** 0.5, ka_exact_C_eq_S_plus_2N=(fEx2 ** 0.5) / fC,
              literal_C_minus_S_minus_N=(fLi2 ** 0.5) / fC, trace_C=trC, trace_S=trS, trace_N=trN,
              signal_share_of_trace=trS / trC if trC else float("nan"))
    return G_U, G_D, ka


def main(line: str) -> None:
    t0 = time.time()
    P(f"[{line}] est286.py sha256 {E.sha256_self()}")
    (WORK / "moments").mkdir(parents=True, exist_ok=True)
    (WORK / "prod").mkdir(parents=True, exist_ok=True)

    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos], dtype=np.int64)

    # ---------------------------------------------------------------- the production rows (all 300 panel targets)
    real_bank = E062.Bank.load(line, plolo.DONOR_PATHS[line])
    prep = E062.prepare_main_source(real_bank, PANEL, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=40.0)
    raw_median_real = float(prep.raw_median_norm)
    np.savez_compressed(WORK / "prod" / f"{line}.npz", targets=np.asarray(PANEL), rows=prep.effects.astype(np.float32),
                        avail=np.asarray(prep.target_available, bool), raw_median=raw_median_real,
                        mapped=int(prep.mapped_gene_count))
    P(f"  production rows: {int(prep.target_available.sum())}/300 targets, raw_median {raw_median_real:.6f}, "
      f"mapped {prep.mapped_gene_count}, duplicates {len(prep.duplicate_symbols)}, rss {rss():.2f} GiB")
    # KA against 271's saved `real` rows (the same object through halfrows271's path)
    z271 = np.load(E271 / f"work/half/rows-{line}.npz", allow_pickle=True)
    t271 = [str(t) for t in z271["targets"]]
    idx271 = np.asarray([t271.index(t) for t in PANEL], dtype=np.int64)
    d271 = float(np.abs(z271["real"][idx271].astype(np.float64) - prep.effects.astype(np.float32).astype(np.float64)).max())
    P(f"  KA-P production rows == 271's `real` rows: max|d| {d271:.3e} ({'PASS' if d271 == 0.0 else 'CHECK'})")
    amap = E062.map_bank_to_axis(real_bank, OFFICIAL)
    del prep
    gc.collect()

    # ---------------------------------------------------------------- rebuild the two half banks (267's sums)
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
        all_position = {int(s): i for i, s in enumerate(hmap.all_source_columns.tolist())}
        unique_mean = np.asarray([mean_all[all_position[int(s)]] for s in hmap.unique_source_columns], dtype=np.float64)
        sides[side] = dict(bank=hb, meta=meta, unique_mean=unique_mean, src=hmap.unique_source_columns,
                           pos=hmap.unique_target_positions)
        P(f"  half {side}: {meta['n_perts']} bank rows x {meta['n_genes']} genes, {meta['n_control_cells']:,} control "
          f"cells, median {meta['cells_per_pert']['median']:.0f} cells/row, mapped {hmap.unique_target_positions.size} "
          f"official columns, rss {rss():.2f} GiB, {time.time()-t0:.0f}s")

    # the fitted columns: the official positions BOTH half banks map onto, minus the panel columns (production zeroes
    # those).  The half banks carry 267's `keep` gene set, which is not the real donor table's gene set: the fitted
    # subspace is therefore a SUBSET of the production row's support, and the share of the row's energy inside it is
    # reported (columns outside `cols` pass through every arm unchanged).
    pos = np.intersect1d(sides["A"]["pos"], sides["B"]["pos"])
    pos = pos[~np.isin(pos, panel_positions)]
    cols = pos
    m = int(cols.size)
    for side in ("A", "B"):
        take = {int(p): i for i, p in enumerate(sides[side]["pos"].tolist())}
        sel = np.asarray([take[int(p)] for p in cols], dtype=np.int64)
        sides[side]["sel"] = sel
    in_map = np.isin(cols, amap.unique_target_positions)
    P(f"  fitted columns: {m} (half-A {sides['A']['pos'].size}, half-B {sides['B']['pos'].size}, real bank maps "
      f"{amap.unique_target_positions.size}); {int(in_map.sum())}/{m} of them are in the real bank's map")
    # the fitted rows: bank perturbations present in BOTH halves and NOT a panel target
    names = [n for n in sides["A"]["bank"].perts if n in sides["B"]["bank"].pert_pos and n not in PANEL_SET]
    names.sort()
    n_rows = len(names)
    rowA = np.asarray([sides["A"]["bank"].pert_pos[n] for n in names], dtype=np.int64)
    rowB = np.asarray([sides["B"]["bank"].pert_pos[n] for n in names], dtype=np.int64)
    n_panel_in_bank = sum(1 for n in sides["A"]["bank"].perts if n in PANEL_SET)
    P(f"  fit: n = {n_rows} bank rows (both halves, panel targets excluded: {n_panel_in_bank} panel rows in the half "
      f"bank), m = {m} columns (mapped {pos.size} minus {pos.size - m} panel columns)")

    scale = EFFECT_SCALE / raw_median_real
    src = {s: sides[s]["src"][sides[s]["sel"]] for s in ("A", "B")}
    mu = {s: sides[s]["unique_mean"][sides[s]["sel"]] for s in ("A", "B")}
    delt = {s: np.load(sides[s]["bank"].delta_path, mmap_mode="r") for s in ("A", "B")}
    row_of = {"A": rowA, "B": rowB}

    def reader(side: str):
        def f(lo, hi):
            if hi is None or isinstance(hi, np.ndarray):     # permuted gather (the Cs control)
                r = row_of[side][hi]
            else:
                if lo >= n_rows:
                    return None
                r = row_of[side][lo:min(hi, n_rows)]
            order = np.argsort(r)
            blk = np.asarray(delt[side][np.sort(r)][:, src[side]], dtype=np.float64)
            out = np.empty_like(blk)
            out[order] = blk
            return (out - mu[side]) * scale
        return f

    rA, rB = reader("A"), reader("B")

    # KA-H: the same arithmetic on the PANEL targets must reproduce 271's half rows on `cols` (271 built them through
    # E062.prepare_main_source on the half bank and rescaled by raw_median_half / raw_median_real).
    ka_h = {}
    for side, key in (("A", "A"), ("B", "B")):
        pn = [n for n in PANEL if n in sides[side]["bank"].pert_pos]
        r = np.asarray([sides[side]["bank"].pert_pos[n] for n in pn], dtype=np.int64)
        blk = (np.asarray(delt[side][np.sort(r)][:, src[side]], dtype=np.float64) - mu[side]) * scale
        blk = blk[np.argsort(np.argsort(r))]
        ref = z271[key][np.asarray([t271.index(n) for n in pn], dtype=np.int64)][:, cols].astype(np.float64)
        ka_h[side] = float(np.abs(blk - ref).max()) / float(np.abs(ref).max())
        del blk, ref
    P(f"  KA-H half rows == 271's rows-{line}.npz on `cols`: rel max|d| A {ka_h['A']:.2e}, B {ka_h['B']:.2e}")
    prow = z271["real"][idx271].astype(np.float64)
    e_in = float((prow[:, cols] ** 2).sum() / (prow ** 2).sum())
    P(f"  production-row energy inside `cols`: {100*e_in:.2f} %")
    del prow

    out: dict = dict(line=line, est286_sha256=E.sha256_self(), n_rows=n_rows, m=m, raw_median_real=raw_median_real,
                     ka_half_rows_rel=ka_h, prod_energy_share_in_cols=e_in,
                     n_panel_rows_in_half_bank=n_panel_in_bank, mapped=int(pos.size), n_panel_cols=int(pos.size - m),
                     halves={s: sides[s]["meta"] for s in ("A", "B")}, ka_prod_vs_271=d271, block=BLOCK)

    for tag, perm in (("", None), ("_shuf", E.derangement(n_rows, E.DERANGE_SEED))):
        t1 = time.time()
        G_U, G_D, n = E.accumulate(rA, rB, m, perm=perm, block=BLOCK)
        assert n == n_rows
        # est286.moments / est286.moment_ka arithmetic done IN PLACE with blockwise Frobenius norms: at m = 17.8k the
        # library form materialises S, N and C plus temporaries (~23 GiB); this form holds two m x m float64 buffers
        # (~5.3 GiB).  Algebraically identical -- KA-Mi below re-runs the library form on k562_gwps and compares.
        S, N, ka = moments_inplace(G_U, G_D, n)
        del G_U, G_D
        gc.collect()
        lead = float(np.linalg.eigvalsh(S[:64, :64]).max()) if m >= 64 else float("nan")   # cheap sanity only
        out[f"moments{tag}"] = dict(ka, seconds=round(time.time() - t1, 1), lead64=lead)
        np.save(WORK / "moments" / f"{line}_S{tag}.npy", S.astype(np.float32))
        np.save(WORK / "moments" / f"{line}_N{tag}.npy", N.astype(np.float32))
        P(f"  moments{tag or ' (real pairing)'}: n {n}, KA-M exact ||C-S-2N||/||C|| {ka['ka_exact_C_eq_S_plus_2N']:.3e}; "
          f"literal ||C-S-N||/||C|| {ka['literal_C_minus_S_minus_N']:.4f} (= ||N||/||C||); tr S / tr C "
          f"{ka['signal_share_of_trace']:+.4f}; {time.time()-t1:.0f}s, rss {rss():.2f} GiB")
        del S, N
        gc.collect()

    np.save(WORK / "moments" / f"{line}_cols.npy", cols)
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = rss()
    out["sha256"] = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
                     for p in sorted((WORK / "moments").glob(f"{line}_*.npy"))}
    (WORK / "moments" / f"{line}.json").write_text(json.dumps(out, indent=1, default=str))
    (WORK / "moments" / f"{line}.log").write_text("\n".join(LOG) + "\n")
    for s in ("A", "B"):
        (sides[s]["bank"].path / "delta.npy").unlink(missing_ok=True)      # reproducible from 267's sums in ~1 min
    P(f"[{line}] done {(time.time()-t0)/60:.1f} min, peak rss {rss():.2f} GiB")


if __name__ == "__main__":
    for ln in sys.argv[1:]:
        main(ln)
