"""267 -- items 2 and 3: HALF-TRUTHS IN THE INSTRUMENT'S FORMAT, the known-answer rows, the reliability arrays.

For each held line and split (batch, rand; plus the union for the KA and the guide split for the two-construct
targets) a half-bank is written in the donor-table format (delta.npy = log1p(TS * S / tot) - b_ctrl in float32, exactly
build_orion_donor / build_donor_table's arithmetic, on the bank columns that map onto the official axis), and the
half-truth is produced by the instrument's OWN path: `E062.prepare_main_source` (per-column mean over EVERY
perturbation of the half-bank, i.e. the residual centring of the truth) followed by `E062.raw_truth`, exactly as
`plolo.build_fold` builds the fold truth; the panel positions are zeroed there and the scored columns are taken here.

KA-1  the union of the halves reproduces the donor table: delta rows bit-identical (float32) on the mapped columns,
      b_ctrl bit-identical, the same perturbation set at min_cells 20.
KA-2  the union truth through the instrument's path reproduces the fold cache's truth (float32) on the scored columns.
KA-3  the explicit formula used for the guide split (delta - the batch half's mean) reproduces the E062 path on the
      batch halves (so the guide-split truths are in the same format).

Reliability (item 3): cos(truth_A, truth_B) per target in the instrument's residual space; also the un-centred
('full') cosine and the panel-centred cosine (159's convention) for reference.

The `orion_hek293t` half-B truths (batch B, rand B, guide B) are written ONLY to work/reserve/ with a sha256.
    PYTHONDONTWRITEBYTECODE=1 uv run python truth267.py k562_gwps
"""
from __future__ import annotations

import hashlib
import json
import resource
import shutil
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/267-replicate-validated-oracle"
HALF = HERE / "work/half"
RESERVE = HERE / "work/reserve"
RESERVE.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = HERE / "work/plolo"
E062, OFFICIAL, PANEL = plolo.E062, plolo.OFFICIAL, plolo.PANEL
TS = 5e4
MIN_CELLS_BANK = 20
MIN_HALF = 5
RESERVED = "orion_hek293t"
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(8 << 20), b""):
            h.update(blk)
    return h.hexdigest()


def sha256_array(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


# --------------------------------------------------------------------------------------------- the sums on disk
class Sums:
    def __init__(self, line: str):
        self.line = line
        d = HALF / line
        z = np.load(d / "sums.npz", allow_pickle=True)
        self.keep = np.asarray(z["keep"], dtype=np.int64)
        self.keep_names = [str(x) for x in z["keep_names"]]
        self.official_positions = np.asarray(z["official_positions"], dtype=np.int64)
        self.pert_names = [str(x) for x in z["pert_names"]]
        self.tot = np.asarray(z["tot"], dtype=np.int64)
        self.ncell = np.asarray(z["ncell"], dtype=np.int64)
        self.ctrl_sum = np.asarray(z["ctrl_sum"], dtype=np.float64)     # (4, n_all) batchA, batchB, randA, guideA
        self.ctrl_n = np.asarray(z["ctrl_n"], dtype=np.int64)
        self.acc = {h: np.load(d / f"acc_{h}.npy", mmap_mode="r") for h in ("batchA", "batchB", "randA")}
        self.k562 = line == "k562_gwps"
        if self.k562:
            self.acc["guideA"] = np.load(d / "acc_guideA.npy", mmap_mode="r")
            self.ci = int(z["control_index"])
            self.guide_targets = [str(x) for x in z["guide_split_targets"]]
        else:
            self.guide_codes = np.asarray(z["guide_codes"], dtype=np.int64)
            self.guide_accA = np.asarray(z["guide_accA"])
            self.guide_totA = np.asarray(z["guide_totA"], dtype=np.int64)
            self.guide_ncellA = np.asarray(z["guide_ncellA"], dtype=np.int64)
            self.guide_ncellB = np.asarray(z["guide_ncellB"], dtype=np.int64)
        # the K562 accumulators are on ALL 8,248 genes; the Orion ones on the mapped columns only
        self.on_all_genes = self.k562
        self.n_all = self.ctrl_sum.shape[1]

    def half(self, split: str, side: str):
        """Returns (S rows x mapped columns as a callable over row blocks, tot, ncell, ctrl_sum_all_tokens, ctrl_n)."""
        A, B, R = self.acc["batchA"], self.acc["batchB"], self.acc["randA"]
        cs, cn, tot, nc = self.ctrl_sum, self.ctrl_n, self.tot, self.ncell
        if split == "batch":
            if side == "A":
                return (lambda lo, hi: np.asarray(A[lo:hi], dtype=np.float64)), tot[:, 0], nc[:, 0], cs[0], cn[0]
            return (lambda lo, hi: np.asarray(B[lo:hi], dtype=np.float64)), tot[:, 1], nc[:, 1], cs[1], cn[1]
        if split == "rand":
            if side == "A":
                return (lambda lo, hi: np.asarray(R[lo:hi], dtype=np.float64)), tot[:, 2], nc[:, 2], cs[2], cn[2]
            return ((lambda lo, hi: np.asarray(A[lo:hi], dtype=np.float64) + np.asarray(B[lo:hi], dtype=np.float64)
                     - np.asarray(R[lo:hi], dtype=np.float64)),
                    tot[:, 0] + tot[:, 1] - tot[:, 2], nc[:, 0] + nc[:, 1] - nc[:, 2], cs[0] + cs[1] - cs[2], cn[0] + cn[1] - cn[2])
        if split == "union":
            return ((lambda lo, hi: np.asarray(A[lo:hi], dtype=np.float64) + np.asarray(B[lo:hi], dtype=np.float64)),
                    tot[:, 0] + tot[:, 1], nc[:, 0] + nc[:, 1], cs[0] + cs[1], cn[0] + cn[1])
        raise ValueError(split)


def b_ctrl_of(ctrl_sum_all: np.ndarray) -> np.ndarray:
    """build_*_donor's control profile on ALL tokens, float32."""
    return np.log1p(TS * ctrl_sum_all / ctrl_sum_all.sum()).astype(np.float32)


def write_half_bank(sm: Sums, split: str, side: str, out: Path) -> dict:
    """The half-bank in the donor-table format; rows = bank perturbations with >= MIN_HALF cells in the half
    (>= MIN_CELLS_BANK for the union, the bank's own rule)."""
    S, tot, nc, cs, cn = sm.half(split, side)
    b_all = b_ctrl_of(cs)
    keep = sm.keep
    b_keep = b_all[keep] if sm.on_all_genes else b_all[keep]     # ctrl_sum is always on ALL tokens
    min_cells = MIN_CELLS_BANK if split == "union" else MIN_HALF
    bank_perts = set((ROOT / f"data/donors/{sm.line}/perts.txt").read_text().split("\n"))
    rows = [i for i, nm in enumerate(sm.pert_names) if nm in bank_perts and nc[i] >= min_cells and tot[i] > 0]
    rows = sorted(rows, key=lambda i: sm.pert_names[i])            # the donor tables are written in sorted order
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
    meta = dict(line=sm.line, split=split, side=side, n_perts=len(names), min_cells=min_cells, n_control_cells=int(cn),
                n_genes=len(keep), cells_per_pert=dict(min=int(nc[rows_a].min()), median=float(np.median(nc[rows_a])),
                                                       max=int(nc[rows_a].max())),
                n_bank_perts_below_min=int(sum(1 for i, nm in enumerate(sm.pert_names) if nm in bank_perts and nc[i] < min_cells)))
    (out / "meta.json").write_text(json.dumps(meta, indent=1))
    return meta


def instrument_truth(bank_dir: Path, name: str, targets, panel_positions):
    """The instrument's own path: E062.prepare_main_source + raw_truth (panel positions zeroed inside)."""
    bank = E062.Bank.load(name, bank_dir)
    prep = E062.prepare_main_source(bank, list(targets), OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=40.0)
    truth = E062.raw_truth(prep)
    return truth, np.asarray(prep.target_available, bool), bank


def explicit_truth(bank_dir: Path, targets, panel_positions, mean_override=None):
    """The same arithmetic written out: delta[t] - column mean over every row of the bank, on the official axis, panel
    zeroed.  `mean_override` (a float64 vector on the bank's columns) replaces the bank's own mean (guide split)."""
    delta = np.load(bank_dir / "delta.npy", mmap_mode="r")
    names = (bank_dir / "perts.txt").read_text().split("\n")
    genes = (bank_dir / "genes.txt").read_text().split("\n")
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    cols = np.asarray([pos[g] for g in genes], dtype=np.int64)
    if mean_override is None:
        mean = np.zeros(delta.shape[1])
        for lo in range(0, delta.shape[0], 256):
            mean += np.asarray(delta[lo:lo + 256], dtype=np.float64).sum(0)
        mean /= delta.shape[0]
    else:
        mean = np.asarray(mean_override, dtype=np.float64)
    idx = {n: i for i, n in enumerate(names)}
    T = np.zeros((len(targets), len(OFFICIAL)))
    avail = np.zeros(len(targets), bool)
    for k, t in enumerate(targets):
        i = idx.get(t)
        if i is None:
            continue
        T[k, cols] = np.asarray(delta[i], dtype=np.float64) - mean
        avail[k] = True
    T[:, panel_positions] = 0.0
    return T, avail, mean, cols


def main(line: str) -> None:
    t0 = time.time()
    sm = Sums(line)
    fold = np.load(HERE / f"work/plolo/fold-{line}.npz", allow_pickle=True)
    targets = [str(t) for t in fold["targets"]]
    scored = np.asarray(fold["scored"], dtype=np.int64)
    truth_cache = np.asarray(fold["truth"], dtype=np.float32)
    official_pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([official_pos[g] for g in PANEL if g in official_pos], dtype=np.int64)
    n = len(targets)
    P(f"[{line}] {n} fold targets, {len(scored):,} scored columns, {len(sm.keep):,} mapped bank columns, rss {rss():.2f} GiB")
    out = dict(line=line, targets=targets, ka=[], halves={}, n_targets=n, n_scored=int(len(scored)))
    tidx = {t: k for k, t in enumerate(targets)}

    # ---------------------------------------------------------------- KA-1: the union reproduces the donor table
    ub = HALF / line / "bank-union"
    meta_u = write_half_bank(sm, "union", "U", ub)
    bank = E062.Bank.load(line, plolo.DONOR_PATHS[line])
    d_bank = np.load(bank.delta_path, mmap_mode="r")
    d_union = np.load(ub / "delta.npy", mmap_mode="r")
    names_u = (ub / "perts.txt").read_text().split("\n")
    same_set = names_u == list(bank.perts)
    max_d = 0.0
    n_diff_rows = 0
    for lo in range(0, len(names_u), 1024):
        rows = [bank.pert_pos[nm] for nm in names_u[lo:lo + 1024]]
        a = np.asarray(d_bank[rows][:, sm.keep], dtype=np.float32)
        b = np.asarray(d_union[lo:lo + 1024], dtype=np.float32)
        dd = np.abs(a - b).max(1)
        max_d = max(max_d, float(dd.max()))
        n_diff_rows += int((dd > 0).sum())
    b_bank = np.asarray(np.load(bank.b_ctrl_path, mmap_mode="r"), dtype=np.float32)[sm.keep]
    b_u = np.load(ub / "b_ctrl.npy")
    ka1 = dict(label="KA-1 union of halves == donor table (float32 delta on the mapped columns, b_ctrl, pert set)",
               same_pert_set=bool(same_set), n_union=len(names_u), n_bank=len(bank.perts), max_abs_diff_delta=max_d,
               rows_with_any_diff=n_diff_rows, max_abs_diff_bctrl=float(np.abs(b_bank - b_u).max()),
               n_control_cells_union=meta_u["n_control_cells"], ok=bool(same_set and max_d == 0.0 and np.array_equal(b_bank, b_u)))
    out["ka"].append(ka1)
    P(f"  {ka1['label']}: {'PASS' if ka1['ok'] else 'FAIL'} max|d| {max_d:.3e} rows differing {n_diff_rows}, "
      f"b_ctrl max|d| {ka1['max_abs_diff_bctrl']:.3e}, pert set {'same' if same_set else 'DIFFERENT'} ({len(names_u)} vs {len(bank.perts)})")

    # ---------------------------------------------------------------- KA-2: the union truth reproduces the cache
    T_u, av_u, _ = instrument_truth(ub, f"{line}-union", targets, panel_positions)
    assert av_u.all()
    Tu32 = T_u[:, scored].astype(np.float32)
    dd = np.abs(Tu32.astype(np.float64) - truth_cache[:, scored].astype(np.float64))
    rel = dd.max() / max(float(np.abs(truth_cache[:, scored]).max()), 1e-300)
    ulp = float(np.max(dd / np.maximum(np.spacing(np.abs(truth_cache[:, scored]).astype(np.float32)), 1e-45)))
    ka2 = dict(label="KA-2 union truth through E062.prepare_main_source/raw_truth == fold cache truth (float32, scored columns)",
               max_abs_diff=float(dd.max()), max_rel_to_truth_max=float(rel), max_diff_in_float32_ulps=ulp,
               n_entries_differing=int((dd > 0).sum()), n_entries=int(dd.size),
               ok=bool(ulp <= 1.0 + 1e-9))
    out["ka"].append(ka2)
    P(f"  {ka2['label']}: {'PASS' if ka2['ok'] else 'FAIL'} max|d| {dd.max():.3e} ({ulp:.2f} float32 ulp; {ka2['n_entries_differing']}/{dd.size} entries differ)")

    # ---------------------------------------------------------------- the halves through the instrument's path
    halves = {}
    means = {}
    reserved = {}
    for split in ("batch", "rand"):
        for side in ("A", "B"):
            bd = HALF / line / f"bank-{split}{side}"
            meta = write_half_bank(sm, split, side, bd)
            T, av, _ = instrument_truth(bd, f"{line}-{split}{side}", targets, panel_positions)
            Te, ave, mean_e, cols = explicit_truth(bd, targets, panel_positions)
            d3 = float(np.abs(T[:, scored] - Te[:, scored]).max())
            if split == "batch":
                out["ka"].append(dict(label=f"KA-3 explicit formula == E062 path ({split}{side})", max_abs_diff=d3, ok=bool(d3 < 1e-9)))
                P(f"  KA-3 explicit == E062 path ({split}{side}): {'PASS' if d3 < 1e-9 else 'FAIL'} max|d| {d3:.2e}")
            means[(split, side)] = (mean_e, cols)
            # the uncentred ('full') rows on the scored columns, for the reference cosines
            delta = np.load(bd / "delta.npy", mmap_mode="r")
            names = (bd / "perts.txt").read_text().split("\n")
            idx = {nm: i for i, nm in enumerate(names)}
            F = np.zeros((n, len(OFFICIAL)))
            for k, t in enumerate(targets):
                if t in idx:
                    F[k, cols] = np.asarray(delta[idx[t]], dtype=np.float64)
            F[:, panel_positions] = 0.0
            nc = sm.half(split, side)[2]
            ncells = np.asarray([nc[sm.pert_names.index(t)] if t in sm.pert_names else 0 for t in targets], dtype=np.int64)
            rec = dict(truth=T[:, scored].astype(np.float32), full=F[:, scored].astype(np.float32), avail=av, ncells=ncells,
                       meta=meta)
            key = f"{split}{side}"
            if line == RESERVED and side == "B":
                reserved[key] = rec
            else:
                halves[key] = rec
            P(f"  {key}: {int(av.sum())}/{n} targets available (>= {MIN_HALF} cells), perts in half-bank {meta['n_perts']}, "
              f"controls {meta['n_control_cells']:,}, median cells/target {float(np.median(ncells[av])):.0f}, rss {rss():.2f} GiB")
            shutil.rmtree(bd / "__pycache__", ignore_errors=True)
            (bd / "delta.npy").unlink()                     # 0.3-1.3 GB each; reproducible from the sums
    # ---------------------------------------------------------------- the guide split (two-construct targets)
    if sm.k562:
        g_targets = [t for t in sm.guide_targets if t in tidx]
        S_all, totU, ncU, csU, cnU = sm.half("union", "U")
        gA = sm.acc["guideA"]
        rows = [sm.pert_names.index(t) for t in g_targets]
        SA = np.asarray(gA[rows], dtype=np.float64)
        SU = S_all(0, len(sm.pert_names))[rows] if False else np.stack([S_all(r, r + 1)[0] for r in rows])
        totA = sm.tot[rows, 3]
        ncA = sm.ncell[rows, 3]
        csA = sm.ctrl_sum[3]
        cnA = int(sm.ctrl_n[3])
    else:
        gmap = {int(c): i for i, c in enumerate(sm.guide_codes)}
        g_targets = [t for t in targets if sm.pert_names.index(t) in gmap]
        rows = [sm.pert_names.index(t) for t in g_targets]
        gi = [gmap[r] for r in rows]
        S_all = sm.half("union", "U")[0]
        SA = np.asarray(sm.guide_accA[gi], dtype=np.float64)
        SU = np.stack([S_all(r, r + 1)[0] for r in rows])
        totA = sm.guide_totA[gi]
        ncA = sm.guide_ncellA[gi]
        csA = sm.ctrl_sum[3]
        cnA = int(sm.ctrl_n[3])
        totU, ncU, csU, cnU = sm.half("union", "U")[1:]
    if g_targets:
        if sm.on_all_genes:
            SA = SA[:, sm.keep]
            SU = SU[:, sm.keep]
        SB = SU - SA
        totB = totU[rows] - totA
        ncB = ncU[rows] - ncA
        csB = csU - csA
        bA = b_ctrl_of(csA)[sm.keep]
        bB = b_ctrl_of(csB)[sm.keep]
        dA = np.log1p(TS * SA / totA[:, None]).astype(np.float32) - bA
        dB = np.log1p(TS * SB / totB[:, None]).astype(np.float32) - bB
        for side, d, ncs in (("A", dA, ncA), ("B", dB, ncB)):
            mean_e, cols = means[("batch", side)]
            T = np.zeros((n, len(OFFICIAL)))
            F = np.zeros((n, len(OFFICIAL)))
            av = np.zeros(n, bool)
            ncells = np.zeros(n, np.int64)
            for k, t in enumerate(g_targets):
                T[tidx[t], cols] = d[k].astype(np.float64) - mean_e
                F[tidx[t], cols] = d[k].astype(np.float64)
                av[tidx[t]] = ncs[k] >= MIN_HALF
                ncells[tidx[t]] = ncs[k]
            T[:, panel_positions] = 0.0
            F[:, panel_positions] = 0.0
            rec = dict(truth=T[:, scored].astype(np.float32), full=F[:, scored].astype(np.float32), avail=av, ncells=ncells,
                       meta=dict(line=line, split="guide", side=side, targets=g_targets, n_control_cells=int(cnA if side == "A" else cnU - cnA),
                                 centring="the batch-" + side + " half-bank's mean (same depth regime); controls = the guide-" + side + " control set"))
            key = f"guide{side}"
            if line == RESERVED and side == "B":
                reserved[key] = rec
            else:
                halves[key] = rec
        P(f"  guide split: {len(g_targets)} two-construct targets {g_targets}; cells A {ncA.tolist()} B {ncB.tolist()}; "
          f"controls A {cnA:,} B {int(cnU - cnA):,}")
    out["guide_targets"] = g_targets

    # ---------------------------------------------------------------- reliability arrays (item 3)
    def cosr(A, B):
        na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            c = (A * B).sum(1) / (na * nb)
        return np.where((na > 0) & (nb > 0), c, np.nan)

    rel = {}
    allh = dict(halves)
    allh.update(reserved)
    for split in ("batch", "rand", "guide"):
        if f"{split}A" not in allh or f"{split}B" not in allh:
            continue
        A, B = allh[f"{split}A"], allh[f"{split}B"]
        ok = A["avail"] & B["avail"]
        TA, TB = A["truth"].astype(np.float64), B["truth"].astype(np.float64)
        FA, FB = A["full"].astype(np.float64), B["full"].astype(np.float64)
        # panel-centred (159's convention): minus the mean over the AVAILABLE fold targets of that side
        PA = FA - FA[ok].mean(0)
        PB = FB - FB[ok].mean(0)
        rel[split] = dict(cos_resid=np.where(ok, cosr(TA, TB), np.nan), cos_full=np.where(ok, cosr(FA, FB), np.nan),
                          cos_panel=np.where(ok, cosr(PA, PB), np.nan), ok=ok,
                          ncells_A=A["ncells"], ncells_B=B["ncells"],
                          cos_truthA_cache=cosr(TA, truth_cache[:, scored].astype(np.float64)),
                          cos_truthB_cache=cosr(TB, truth_cache[:, scored].astype(np.float64)))
        c = rel[split]["cos_resid"][ok]
        P(f"  reliability {split}: n {int(ok.sum())}  cos(truth_A, truth_B) residual space: mean {np.mean(c):+.4f} "
          f"p10 {np.percentile(c, 10):+.4f} median {np.median(c):+.4f} p90 {np.percentile(c, 90):+.4f}  share<=0 {np.mean(c <= 0):.3f} | "
          f"full {np.nanmedian(rel[split]['cos_full']):+.4f} panel-centred {np.nanmedian(rel[split]['cos_panel']):+.4f}")

    # ---------------------------------------------------------------- save
    save = dict(targets=np.asarray(targets), scored=scored, truth_cache_scored=truth_cache[:, scored])
    for k, r in halves.items():
        for kk in ("truth", "full", "avail", "ncells"):
            save[f"{k}|{kk}"] = r[kk]
    for split, r in rel.items():
        for kk, v in r.items():
            save[f"rel|{split}|{kk}"] = v
    np.savez_compressed(HALF / f"truth-{line}.npz", **save)
    out["halves"] = {k: r["meta"] for k, r in halves.items()}
    if reserved:
        rs = dict(targets=np.asarray(targets), scored=scored)
        for k, r in reserved.items():
            for kk in ("truth", "full", "avail", "ncells"):
                rs[f"{k}|{kk}"] = r[kk]
        rp = RESERVE / f"truthB-{line}.npz"
        np.savez_compressed(rp, **rs)
        sha = sha256_file(rp)
        (RESERVE / "SHA256SUMS.txt").write_text(f"{sha}  {rp.name}\n")
        out["reserve"] = dict(file=str(rp.relative_to(ROOT)), sha256=sha, keys=sorted(reserved),
                              content_sha256={k: sha256_array(r["truth"]) for k, r in reserved.items()},
                              note="the orion_hek293t half-B truths; used ONCE, in sweep267.py's item-4 scoring pass, and in item 3's reliability above")
        out["halves"].update({k: r["meta"] for k, r in reserved.items()})
        P(f"  RESERVE written: {rp.relative_to(ROOT)} sha256 {sha}")
    out["seconds"] = time.time() - t0
    out["peak_rss_gib"] = rss()
    (HALF / f"truth-{line}.json").write_text(json.dumps(out, indent=1, default=str))
    (HALF / f"truth-{line}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{line}] done {out['seconds'] / 60:.1f} min, peak rss {rss():.2f} GiB")


if __name__ == "__main__":
    for line in sys.argv[1:]:
        main(line)
