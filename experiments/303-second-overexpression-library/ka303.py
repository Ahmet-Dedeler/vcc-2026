"""303 — the LOADER known-answer rows.  Reads NO reserve value.

KA-L1  the bank's `delta.npy` rows reproduce through lib303's transform for a
       fixed seeded 20-target sample of `Pan2025_hESC` (free rows).
KA-L2  the bank's `b_ctrl.npy` reproduces from the first Pan2025 control row.
KA-L3  the bank's `measured_mask.npy` rows reproduce for the same 20 targets.
KA-L4  the bank's `rows.tsv` source-row map reproduces from `obs` alone.
KA-C   279's `Ccross` rows involving `gse295214_h1` reproduce at 0.0e+00 through
       279's own `load_unit` (302's KA-C, recomputed in THIS experiment's work dir).
KA-W   importing 276/279 wrote nothing outside `experiments/303-*/work`.

A failure of any LOADER row is an ABORT: no release may be issued.

    PYTHONDONTWRITEBYTECODE=1 uv run python ka303.py
"""
from __future__ import annotations

import csv
import resource
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib303 as L  # noqa: E402

KA_SEED = 20260907
N_SAMPLE = 20
DELTA_TOL = 1e-6      # float32 storage of values up to ~1e1 -> ~1 ulp at magnitude 8


def _snapshot(dirs) -> dict:
    out = {}
    for d in dirs:
        if d.exists():
            for p in sorted(d.rglob("*")):
                if p.is_file():
                    st = p.stat()
                    out[str(p)] = (st.st_size, st.st_mtime_ns)
    return out


def main() -> None:
    t0 = time.time()
    L.WORK.mkdir(exist_ok=True)
    watch = [L.ROOT / "experiments/276-conservation-map/work",
             L.ROOT / "experiments/279-znf-state-dependence/work",
             L.ROOT / "experiments/290-remainder-response-features/work",
             L.ROOT / "experiments/302-orf-vs-crispri-symmetry/work",
             L.H1BANK]
    before = _snapshot(watch)
    res: dict = {}

    m = L.read_db_meta()
    canon = L.canonical_genes()
    bank = L.read_bank()
    proj = L.build_projection(m.var_symbol, canon)
    bnk_rows = np.flatnonzero(m.dataset == L.BANK_DATASET).astype(np.int64)

    # ---------------------------------------------------------------- KA-L4
    tgt_to_src = {str(s): int(r) for s, r in zip(m.pert_symbol[bnk_rows], bnk_rows)}
    stated = {}
    with open(L.H1BANK / "rows.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stated[r["target_gene"]] = int(r["source_row_index"])
    same = (set(stated) == set(tgt_to_src)) and all(stated[k] == tgt_to_src[k] for k in stated)
    order_ok = [r["target_gene"] for r in csv.DictReader(
        open(L.H1BANK / "rows.tsv"), delimiter="\t")] == bank["perts"]
    res["KA-L4"] = dict(n=len(stated), map_matches=bool(same), row_order_matches=bool(order_ok),
                        pass_=bool(same and order_ok))
    print(f"KA-L4  rows.tsv source map from obs alone: {len(stated)} rows, "
          f"map {same}, order {order_ok} -> {'PASS' if res['KA-L4']['pass_'] else 'FAIL'}")

    # ---------------------------------------------------------------- KA-L1/3
    rng = np.random.default_rng(KA_SEED)
    pick = np.sort(rng.choice(len(bank["perts"]), N_SAMPLE, replace=False))
    targets = [bank["perts"][int(i)] for i in pick]
    src_rows = np.asarray([tgt_to_src[t] for t in targets], np.int64)
    order = np.argsort(src_rows, kind="stable")
    src_sorted = src_rows[order]

    blk = L.read_source_rows(L.DB, src_sorted, L.FREE_DATASETS, m.dataset)
    delta, measured, diag = L.bank_transform(
        blk["X"], blk["layers/control"], blk["layers/measure_mask"], proj,
        [m.obs_id[i] for i in src_sorted], strict_sums=True)
    got_delta = np.empty_like(delta)
    got_mask = np.empty_like(measured)
    got_delta[order] = delta
    got_mask[order] = measured

    ref_delta = np.asarray(np.load(L.H1BANK / "delta.npy", mmap_mode="r")[pick])
    ref_mask = np.asarray(np.load(L.H1BANK / "measured_mask.npy", mmap_mode="r")[pick])
    d = float(np.max(np.abs(got_delta.astype(np.float64) - ref_delta.astype(np.float64))))
    bit = bool(np.array_equal(got_delta, ref_delta))
    mask_same = bool(np.array_equal(got_mask, ref_mask))
    res["KA-L1"] = dict(n_targets=N_SAMPLE, seed=KA_SEED, targets=targets,
                        max_abs_diff=d, bit_identical=bit, tol=DELTA_TOL,
                        cpm_sum_drift=diag["drift"], pass_=bool(d <= DELTA_TOL))
    res["KA-L3"] = dict(mask_identical=mask_same, pass_=mask_same)
    print(f"KA-L1  delta.npy for {N_SAMPLE} seeded Pan2025 targets: max|diff| {d:.1e} "
          f"(bit-identical {bit}, tol {DELTA_TOL:.0e}) -> {'PASS' if res['KA-L1']['pass_'] else 'FAIL'}")
    print(f"KA-L3  measured_mask.npy rows identical: {mask_same} -> "
          f"{'PASS' if mask_same else 'FAIL'}")

    # ---------------------------------------------------------------- KA-L2
    first = int(bnk_rows.min())
    c1 = L.read_source_rows(L.DB, np.asarray([first], np.int64), L.FREE_DATASETS, m.dataset,
                            what=("layers/control",))["layers/control"][0]
    b = L.control_profile(c1, proj)
    ref_b = np.asarray(np.load(L.H1BANK / "b_ctrl.npy"))
    db = float(np.max(np.abs(b.astype(np.float64) - ref_b.astype(np.float64))))
    res["KA-L2"] = dict(source_row=first, max_abs_diff=db,
                        bit_identical=bool(np.array_equal(b, ref_b)),
                        pass_=bool(db <= DELTA_TOL))
    print(f"KA-L2  b_ctrl.npy from source row {first}: max|diff| {db:.1e} -> "
          f"{'PASS' if res['KA-L2']['pass_'] else 'FAIL'}")

    # -------------------------------------------- Pan2025 control-row identity
    # 115 asserts every Pan2025 control row equals the first.  Re-check on the
    # sample: this property is NOT guaranteed for R5 and is an access-time gate.
    same_ctrl = bool(np.array_equal(blk["layers/control"],
                                    np.broadcast_to(c1, blk["layers/control"].shape)))
    res["KA-L5"] = dict(pan2025_control_rows_identical_on_sample=same_ctrl, pass_=same_ctrl)
    print(f"KA-L5  Pan2025 control rows identical to row {first} on the sample: {same_ctrl} -> "
          f"{'PASS' if same_ctrl else 'FAIL'}")
    del blk

    # ---------------------------------------------------------------- KA-C
    sys.path.insert(0, str(L.ROOT / "experiments/276-conservation-map"))
    sys.path.insert(0, str(L.ROOT / "experiments/279-znf-state-dependence"))
    import build279 as B  # noqa: E402
    import inv279 as I  # noqa: E402
    SRC = I.CANDIDATES

    z = np.load(L.HELD279, allow_pickle=True)
    keys = np.asarray(z["keys"]).astype(str).tolist()
    axis279 = np.asarray(z["axis"]).astype(str).tolist()
    genes279 = np.asarray(z["genes"]).astype(str).tolist()
    pairs = np.asarray(z["pairs_cross"])
    Ccross = np.asarray(z["Ccross"])
    h1 = keys.index("gse295214_h1")
    sel = [p for p in range(len(pairs)) if h1 in tuple(pairs[p])]
    idxH, UH, _ = B.load_unit(SRC["gse295214_h1"][0], axis279, genes279)
    worst, n_cmp, per_pair = 0.0, 0, {}
    for p in sel:
        i, j = int(pairs[p][0]), int(pairs[p][1])
        other = keys[j] if i == h1 else keys[i]
        idxO, UO, _ = B.load_unit(SRC[other][0], axis279, genes279)
        g, c = (B.pair_cos(idxH, UH, idxO, UO) if i == h1 else B.pair_cos(idxO, UO, idxH, UH))
        ref = Ccross[p]
        got = np.full(ref.shape, np.nan, np.float32)
        got[g] = c
        both = np.isfinite(ref) | np.isfinite(got)
        dd = float(np.nanmax(np.abs(np.nan_to_num(ref[both], nan=0.0)
                                    - np.nan_to_num(got[both], nan=0.0)))) if both.any() else 0.0
        same_n = int(np.isfinite(ref).sum()) == int(np.isfinite(got).sum()) == g.size
        worst = max(worst, dd); n_cmp += int(g.size)
        per_pair[f"{keys[i]}|{keys[j]}"] = dict(n=int(g.size), max_abs_diff=dd, same_n=same_n)
        del idxO, UO
    ok = bool(worst == 0.0 and all(v["same_n"] for v in per_pair.values()))
    res["KA-C"] = dict(n_pairs=len(sel), n_cosines=n_cmp, max_abs_diff=worst,
                       pass_=ok, per_pair=per_pair)
    print(f"KA-C   279 Ccross, {len(sel)} pairs / {n_cmp} (pair,gene) cosines: "
          f"max|diff| {worst:.1e} -> {'PASS' if ok else 'FAIL'}")

    # ---------------------------------------------------------------- KA-W
    after = _snapshot(watch)
    changed = sorted(k for k in set(before) | set(after) if before.get(k) != after.get(k))
    res["KA-W"] = dict(n_watched=len(before), changed=changed, pass_=not changed)
    print(f"KA-W   files changed outside 303 during imports/reads: {len(changed)} -> "
          f"{'PASS' if not changed else 'FAIL'}")

    res["summary"] = dict(
        all_pass=all(v.get("pass_") for k, v in res.items() if k.startswith("KA")),
        peak_rss_gib=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, 2),
        seconds=round(time.time() - t0, 1))
    L.dump_json(L.WORK / "ka303.json", res)
    print(f"\nALL LOADER KA ROWS PASS: {res['summary']['all_pass']}  "
          f"(peak RSS {res['summary']['peak_rss_gib']:.2f} GiB, {res['summary']['seconds']}s)")
    if not res["summary"]["all_pass"]:
        raise SystemExit("KA FAILURE — abort before any reserve access")


if __name__ == "__main__":
    sys.exit(main())
