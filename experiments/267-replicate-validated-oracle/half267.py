"""267 -- item 2, the heavy step: ONE streaming pass per held line over the raw cells, accumulating EXACT integer
pseudobulk sums per (perturbation, half) for every perturbation the bank carries (the centring mean of the
instrument is over ALL of a bank's perturbations, so every perturbation is split, not only the panel targets).

Halves (assignments fixed by inv267.py, read from work/inv/<line>-cells.npz, never re-drawn here):
  batch split : A = even batch / gemgroup number, B = odd (controls split by the same rule -> own controls per half)
  random split: A = cells with half_rand == 0 (seed 20260905), B = the rest; B is derived as (A_batch + B_batch) - A_rand
  guide split (K562 only, the 8 two-construct panel targets): A = the construct that sorts first; controls: the
                non-targeting constructs that sort even -> A, odd -> B.
Accumulators live on the bank's columns that map onto the official axis (E062.map_bank_to_axis's
unique_source_columns, the columns the truth is read from); per-(perturbation, half) TOTALS over ALL tokens are kept
separately because the pseudobulk normalises by the total over every gene (build_orion_donor / build_donor_table).
Writes work/half/<line>/sums.npz (+ .npy memmaps for the big int32 arrays).  Peak RSS ~6 GiB on HEK293T.

    PYTHONDONTWRITEBYTECODE=1 uv run python half267.py orion_hct116
"""
from __future__ import annotations

import json
import re
import resource
import sys
import time
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/267-replicate-validated-oracle"
INV = HERE / "work/inv"
HALF = HERE / "work/half"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402  (E062, OFFICIAL, DONOR_PATHS; read-only)

plolo.WORK = HERE / "work/plolo"
E062, OFFICIAL = plolo.E062, plolo.OFFICIAL
CHUNK_ROWS = 2500
K562_BLOCK = 3886 * 4


def ctrl_guide_is_A(g: str) -> bool:
    """Control guide-split: parity of the first non-targeting id ('non-targeting_03016|...' -> 3016 even -> A).  The
    donor builder's control regex (^neg, ^ntc, ^control, ^safe) also pools a few GENES (e.g. NEGR1) into the control
    set; their constructs carry no id, so they are assigned by the parity of the construct's byte sum."""
    head = g.split("|")[0]
    if head.startswith("non-targeting_"):
        return int(head.rsplit("_", 1)[1]) % 2 == 0
    return sum(head.encode()) % 2 == 0


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def bank_and_map(name: str):
    bank = E062.Bank.load(name, plolo.DONOR_PATHS[name])
    amap = E062.map_bank_to_axis(bank, OFFICIAL)
    keep = np.asarray(amap.unique_source_columns, dtype=np.int64)          # bank columns, in official-position order
    return bank, amap, keep


# --------------------------------------------------------------------------------------------------- Orion
def orion(name: str) -> None:
    """One file at a time (each parquet is a single row group, ~2 GB in arrow); the arrays are sliced in place and the
    group sums are ONE keyed np.bincount per slice (key = local group * n_keep + column) -- no sparse matrix."""
    line = {"orion_hct116": "HCT116", "orion_hek293t": "HEK293T"}[name]
    out = HALF / name
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    bank, amap, keep = bank_and_map(name)
    n_tok = len(bank.genes)
    n_keep = len(keep)
    col_map = np.full(n_tok, -1, dtype=np.int64)
    col_map[keep] = np.arange(n_keep)
    inv = np.load(INV / f"{name}-cells.npz", allow_pickle=True)
    pert_code = inv["pert_code"]
    pert_names = [str(x) for x in inv["pert_names"]]
    ctrl = inv["ctrl"]
    half_batch = inv["half_batch"]
    half_rand = inv["half_rand"]
    file_order = [str(x) for x in inv["file_order"]]
    file_kept = inv["file_kept"]
    P = len(pert_names)
    print(f"[{name}] {len(pert_code):,} cells, {P:,} gene_target values, {n_keep:,}/{n_tok:,} mapped columns, "
          f"bank {bank.shape}", flush=True)
    acc = {h: np.lib.format.open_memmap(out / f"acc_{h}.npy", mode="w+", dtype=np.int32, shape=(P, n_keep))
           for h in ("batchA", "batchB", "randA")}
    for a in acc.values():
        a[:] = 0
    tot = np.zeros((P, 3), dtype=np.int64)          # totals over ALL tokens: batchA, batchB, randA
    ncell = np.zeros((P, 3), dtype=np.int64)
    ctrl_sum = np.zeros((4, n_tok), dtype=np.float64)   # batchA, batchB, randA, guideA on ALL tokens (exact integers)
    ctrl_n = np.zeros(4, dtype=np.int64)
    acc_gA: dict[int, np.ndarray] = {}
    tot_gA: dict[int, int] = {}
    ncell_gA: dict[int, int] = {}
    ncell_gB: dict[int, int] = {}
    files = [ROOT / "data/external/orion/data" / f for f in file_order]
    pos = 0
    n_slices = 0
    pool = pa.system_memory_pool()          # malloc/free: returns memory to the OS (mimalloc retained ~1 GB per file)
    pa.set_memory_pool(pool)
    for i, f in enumerate(files, 1):
        bno = int(re.search(r"Batch(\d+)", f.name)[1])
        hb = bno % 2
        pf = pq.ParquetFile(f, memory_map=False)
        got = 0
        for rb in pf.iter_batches(batch_size=CHUNK_ROWS, columns=["gene_token_id", "gene_expression", "gene_target",
                                                                   "guide_target", "pass_guide_filter", "total_counts"]):
            tok = rb.column("gene_token_id")
            off = np.asarray(tok.offsets, dtype=np.int64)
            vals = np.asarray(tok.values)
            expr = np.asarray(rb.column("gene_expression").values)
            pgf = rb.column("pass_guide_filter").to_numpy()
            file_tot = rb.column("total_counts").to_numpy()
            gt = np.asarray(rb.column("gene_target").to_numpy(zero_copy_only=False), dtype=object)
            guide_all = np.asarray(rb.column("guide_target").to_numpy(zero_copy_only=False), dtype=object)
            n_rows = len(pgf)
            keep_rows = pgf == 1
            kept = np.flatnonzero(keep_rows)
            nk = len(kept)
            tot_all = np.diff(np.concatenate([[0.0], np.cumsum(expr)])[off])        # exact: integer-valued doubles
            assert np.array_equal(tot_all, file_tot), f"{f.name}: sum(gene_expression) != total_counts"
            pc_k = pert_code[pos:pos + nk]
            assert all(pert_names[c] == g for c, g in zip(pc_k[:20], gt[kept][:20])), "inventory order mismatch"
            assert (half_batch[pos:pos + nk] == hb).all()
            is_c_k = ctrl[pos:pos + nk]
            hr_k = half_rand[pos:pos + nk]
            guide_k = guide_all[kept]
            g_p1 = np.fromiter(("_P1-" in g for g in guide_k), bool, nk)
            g_p2 = np.fromiter(("_P2-" in g for g in guide_k), bool, nk)
            g_ce = np.fromiter((ctrl_guide_is_A(g) if c else False for g, c in zip(guide_k, is_c_k)), bool, nk)
            is_c_row = np.zeros(n_rows, bool); is_c_row[kept] = is_c_k
            hr_row = np.full(n_rows, -1, np.int8); hr_row[kept] = hr_k
            pc_row = np.full(n_rows, -1, np.int64); pc_row[kept] = pc_k
            p1_row = np.zeros(n_rows, bool); p1_row[kept] = g_p1 & ~is_c_k
            ge_row = np.zeros(n_rows, bool); ge_row[kept] = g_ce
            row_id = np.repeat(np.arange(n_rows), np.diff(off))
            keep_e = keep_rows[row_id]
            cm = keep_e & is_c_row[row_id]
            if cm.any():
                ctrl_sum[hb] += np.bincount(vals[cm], weights=expr[cm], minlength=n_tok)
                rm = cm & (hr_row[row_id] == 0)
                ctrl_sum[2] += np.bincount(vals[rm], weights=expr[rm], minlength=n_tok)
                gm = cm & ge_row[row_id]
                ctrl_sum[3] += np.bincount(vals[gm], weights=expr[gm], minlength=n_tok)
            c = col_map[vals]
            pm = keep_e & ~is_c_row[row_id] & (c >= 0)
            if pm.any():
                rid = row_id[pm]
                codes = pc_row[rid]
                uniq, local = np.unique(codes, return_inverse=True)
                U = len(uniq)
                key = (local * 2 + hr_row[rid]) * n_keep + c[pm]
                w = expr[pm]
                blk = np.bincount(key, weights=w, minlength=2 * U * n_keep).reshape(U, 2, n_keep)
                assert blk.max() < 2**30
                bi = blk.astype(np.int32)
                del blk
                acc["batchA" if hb == 0 else "batchB"][uniq] += bi[:, 0] + bi[:, 1]
                acc["randA"][uniq] += bi[:, 0]
                del bi
                p1 = p1_row[rid]
                if p1.any():
                    u1, l1 = np.unique(codes[p1], return_inverse=True)
                    k1 = l1 * n_keep + c[pm][p1]
                    b1 = np.bincount(k1, weights=w[p1], minlength=len(u1) * n_keep).reshape(len(u1), n_keep).astype(np.int32)
                    for k, u in enumerate(u1):
                        u = int(u)
                        if u in acc_gA:
                            acc_gA[u] += b1[k]
                        else:
                            acc_gA[u] = b1[k].copy()
                    del b1
                del rid, codes, key, w
            # per-cell totals and counts
            pk = kept[~is_c_k]
            if pk.size:
                codes = pc_row[pk]
                uniq, local = np.unique(codes, return_inverse=True)
                hr0 = hr_row[pk] == 0
                tot[uniq, hb] += np.rint(np.bincount(local, weights=tot_all[pk], minlength=len(uniq))).astype(np.int64)
                tot[uniq, 2] += np.rint(np.bincount(local, weights=tot_all[pk] * hr0, minlength=len(uniq))).astype(np.int64)
                ncell[uniq, hb] += np.bincount(local, minlength=len(uniq))
                ncell[uniq, 2] += np.bincount(local, weights=hr0, minlength=len(uniq)).astype(np.int64)
            p1k = kept[g_p1 & ~is_c_k]
            if p1k.size:
                u1, l1, cnt1 = np.unique(pc_row[p1k], return_inverse=True, return_counts=True)
                ts1 = np.bincount(l1, weights=tot_all[p1k], minlength=len(u1))
                for u, cnt, tsum in zip(u1, cnt1, ts1):
                    u = int(u)
                    ncell_gA[u] = ncell_gA.get(u, 0) + int(cnt)
                    tot_gA[u] = tot_gA.get(u, 0) + int(round(tsum))
            p2k = kept[g_p2 & ~is_c_k]
            if p2k.size:
                for u, cnt in zip(*np.unique(pc_row[p2k], return_counts=True)):
                    ncell_gB[int(u)] = ncell_gB.get(int(u), 0) + int(cnt)
            ctrl_n[hb] += int(is_c_k.sum())
            ctrl_n[2] += int((is_c_k & (hr_k == 0)).sum())
            ctrl_n[3] += int(g_ce.sum())
            pos += nk
            got += nk
            n_slices += 1
            del rb, tok, off, vals, expr, pgf, file_tot, gt, guide_all, row_id, keep_e, cm, c, pm
        assert got == file_kept[i - 1], f"{f.name}: kept {got} != inventory {file_kept[i-1]}"
        pool.release_unused()
        if i % 10 == 0 or i == len(files):
            print(f"  [{name}] {i}/{len(files)} files, {pos:,} cells, {time.time() - t0:.0f}s, rss {rss():.2f} GiB", flush=True)
    assert pos == len(pert_code)
    for a in acc.values():
        a.flush()
    g_codes = np.asarray(sorted(u for u in acc_gA if u in ncell_gB), dtype=np.int64)   # perts with BOTH constructs
    np.savez(out / "sums.npz", tot=tot, ncell=ncell, ctrl_sum=ctrl_sum, ctrl_n=ctrl_n, keep=keep,
             keep_names=np.asarray([bank.genes[c] for c in keep], dtype=object),
             official_positions=np.asarray(amap.unique_target_positions, dtype=np.int64),
             pert_names=np.asarray(pert_names, dtype=object), halves=np.asarray(["batchA", "batchB", "randA", "guideA"], dtype=object),
             guide_codes=g_codes, guide_accA=np.stack([acc_gA[int(u)] for u in g_codes]) if len(g_codes) else np.zeros((0, n_keep), np.int32),
             guide_totA=np.asarray([tot_gA[int(u)] for u in g_codes], dtype=np.int64),
             guide_ncellA=np.asarray([ncell_gA[int(u)] for u in g_codes], dtype=np.int64),
             guide_ncellB=np.asarray([ncell_gB[int(u)] for u in g_codes], dtype=np.int64))
    json.dump(dict(line=name, n_cells=int(pos), n_pert=int(P), n_keep=int(n_keep), n_tok=int(n_tok), slices=n_slices,
                   ctrl_n=ctrl_n.tolist(), n_two_construct_perts=int(len(g_codes)), minutes=(time.time() - t0) / 60, peak_rss_gib=rss()),
              open(out / "half267.json", "w"), indent=1)
    print(f"[{name}] done in {(time.time() - t0) / 60:.1f} min, peak rss {rss():.2f} GiB", flush=True)


# --------------------------------------------------------------------------------------------------- K562
def k562() -> None:
    import h5py
    name = "k562_gwps"
    out = HALF / name
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    bank, amap, keep = bank_and_map(name)
    inv = np.load(INV / f"{name}-cells.npz", allow_pickle=True)
    pert_code = inv["pert_code"]
    pert_names = [str(x) for x in inv["pert_names"]]
    guide_code = inv["guide_code"]
    guide_names = [str(x) for x in inv["guide_names"]]
    ctrl = inv["ctrl"]
    half_batch = inv["half_batch"]
    half_rand = inv["half_rand"]
    P, G = len(pert_names), len(bank.genes)
    ci = pert_names.index("control")
    # guide split: the 8 two-construct panel targets (A = the construct sorting first); controls by parity of the
    # sorted non-targeting construct index
    targets_inv = json.loads((INV / f"{name}.json").read_text())["targets_with_ge2_guides"]
    gsplit = {}
    for t in targets_inv:
        m = pert_code == pert_names.index(t)
        gs = sorted(set(guide_code[m].tolist()), key=lambda k: guide_names[k])
        gsplit[t] = gs
    guide_is_A = np.zeros(len(guide_names), dtype=bool)
    for t, gs in gsplit.items():
        guide_is_A[gs[0]] = True
    ctrl_guides = sorted(set(guide_code[ctrl].tolist()), key=lambda k: guide_names[k])
    for j, k in enumerate(ctrl_guides):
        guide_is_A[k] = (j % 2 == 0)
    half_guide = np.where(guide_is_A[guide_code], 0, 1).astype(np.int8)
    P_h = h5py.File(ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad", "r")
    X = P_h["X"]
    n, G2 = X.shape
    assert G2 == G and n == len(pert_code)
    acc = np.zeros((4, P, G), dtype=np.int64)          # batchA, batchB, randA, guideA on ALL 8,248 genes (exact)
    ncell = np.zeros((P, 4), dtype=np.int64)
    tot = np.zeros((P, 4), dtype=np.int64)
    ncounts = P_h["obs/ncounts"][:]
    for lo in range(0, n, K562_BLOCK):
        hi = min(lo + K562_BLOCK, n)
        blk = np.asarray(X[lo:hi], dtype=np.float32)
        rows_tot = blk.sum(1, dtype=np.float64)
        assert np.allclose(rows_tot, ncounts[lo:hi])
        pc = pert_code[lo:hi]
        uniq, local = np.unique(pc, return_inverse=True)
        for h, bits in ((0, half_batch[lo:hi] == 0), (1, half_batch[lo:hi] == 1), (2, half_rand[lo:hi] == 0),
                        (3, half_guide[lo:hi] == 0)):
            sel = np.flatnonzero(bits)
            if not sel.size:
                continue
            ind = sp.csr_matrix((np.ones(sel.size, np.float32), (local[sel], sel)), shape=(len(uniq), hi - lo))
            gs = ind @ blk                                    # float32, exact integers within a block
            assert float(gs.max()) < 2**24, "float32 exactness bound exceeded"
            acc[h][uniq] += gs.astype(np.int64)
            ncell[uniq, h] += np.bincount(local[sel], minlength=len(uniq))
            tot[uniq, h] += np.rint(np.bincount(local[sel], weights=rows_tot[sel], minlength=len(uniq))).astype(np.int64)
        if (lo // K562_BLOCK) % 16 == 0:
            print(f"  [{name}] {hi:,}/{n:,} rows, {time.time() - t0:.0f}s, rss {rss():.2f} GiB", flush=True)
    P_h.close()
    ctrl_sum = acc[:, ci, :].astype(np.float64)
    ctrl_n = ncell[ci].copy()
    np.save(out / "acc_batchA.npy", acc[0]); np.save(out / "acc_batchB.npy", acc[1])
    np.save(out / "acc_randA.npy", acc[2]); np.save(out / "acc_guideA.npy", acc[3])
    np.savez(out / "sums.npz", tot=tot, ncell=ncell, ctrl_sum=ctrl_sum, ctrl_n=ctrl_n, keep=keep,
             keep_names=np.asarray([bank.genes[c] for c in keep], dtype=object),
             official_positions=np.asarray(amap.unique_target_positions, dtype=np.int64),
             pert_names=np.asarray(pert_names, dtype=object), control_index=ci,
             halves=np.asarray(["batchA", "batchB", "randA", "guideA"], dtype=object),
             guide_split_targets=np.asarray(list(gsplit), dtype=object),
             guide_split_constructs=np.asarray([[guide_names[k] for k in gsplit[t]] for t in gsplit], dtype=object),
             n_ctrl_guides=len(ctrl_guides))
    json.dump(dict(line=name, n_cells=int(n), n_pert=int(P), n_genes=int(G), n_keep=int(len(keep)), ctrl_n=ctrl_n.tolist(),
                   guide_split_targets=list(gsplit), minutes=(time.time() - t0) / 60, peak_rss_gib=rss()),
              open(out / "half267.json", "w"), indent=1)
    print(f"[{name}] done in {(time.time() - t0) / 60:.1f} min, peak rss {rss():.2f} GiB", flush=True)


if __name__ == "__main__":
    for w in sys.argv[1:]:
        if w == "k562_gwps":
            k562()
        else:
            orion(w)
