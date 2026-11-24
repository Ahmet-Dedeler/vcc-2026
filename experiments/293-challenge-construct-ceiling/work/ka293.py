"""293 -- KNOWN-ANSWER ROWS (PREREGISTRATION.md sec 2).  Exit 1 on any failure.

KA-1  273's lane-parity 200/200 median on h1-test = +0.24291977562969289 at 0.0e+00 with this code.
KA-2  273's 5-target construct row from h1-test (median +0.013512255772155394 and all five per-target
      values, and the construct cell counts 37/363, 191/209, 150/250, 351/49, 338/62) at 0.0e+00.
KA-3  integer identity: lane-half A + lane-half B == the full row, exactly, per target and for the counts.
KA-4  the construct parser on adata_Test.h5ad's obs reproduces 273's five multi-construct targets and their
      guide ids; every h1-test barcode is in the test source with a matching lane suffix.

    PYTHONDONTWRITEBYTECODE=1 uv run python ka293.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib293 as X  # noqa: E402
from lib293 import L  # noqa: E402

P = L.Logger(X.WORK / "ka293.log")
FAIL: list[str] = []

REF_BATCH_MEDIAN = 0.24291977562969289
REF_GUIDE = {
    "DPF2": (-0.03155149708284103, 37, 363),
    "GRB10": (0.013512255772155394, 191, 209),
    "MED15": (-0.09916872748565872, 150, 250),
    "PCCB": (0.027145874176245796, 351, 49),
    "SOX4": (0.35394153008771090, 338, 62),
}
REF_GUIDE_MEDIAN = 0.013512255772155394
REF_GUIDE_FULL = {"DPF2": -0.018979442389958825, "GRB10": 0.030497869069015183,
                  "MED15": 0.008770527769598165, "PCCB": -0.0076139630659650525,
                  "SOX4": 0.32593939920179843}


def main() -> None:
    t0 = time.time()
    P(f"293 ka293.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")
    P(f"  h1-test sha256 {L.sha256_file(X.H1TEST)[:16]}")

    with h5py.File(X.H1TEST, "r") as h:
        lab = X.labels_of(h, "obs/target_gene")
        barcodes = X.labels_of(h, "obs/_index")
        genes = X.labels_of(h, "var/_index")
        n_cells, n_genes = (int(v) for v in h["X"].attrs["shape"])
    targets = np.asarray(sorted(set(lab.tolist()) - {X.CTRL}))
    tpos = {t: i for i, t in enumerate(targets)}
    assert targets.size == 85 and n_cells == 52400

    SRC = X.SPLITS["test"]
    with h5py.File(SRC, "r") as h:
        src_idx = X.labels_of(h, "obs/_index")
        src_batch = X.labels_of(h, "obs/batch")
        src_guide = X.labels_of(h, "obs/guide_id")
        src_tgt = X.labels_of(h, "obs/target_gene")
    pos = {b: i for i, b in enumerate(src_idx)}
    ka4_all_found = all(b in pos for b in barcodes)
    hit = np.asarray([pos[b] for b in barcodes], dtype=np.int64)
    suffix = np.asarray([b.split("-", 1)[1] for b in barcodes])
    ka4_suffix = bool(np.array_equal(suffix, src_batch[hit]))
    guide = src_guide[hit]
    batches = sorted(set(suffix.tolist()))
    bidx = {b: i for i, b in enumerate(batches)}
    parity = np.asarray([bidx[b] % 2 for b in suffix])

    # --- KA-4: the parser on the FULL test split reproduces 273's five multi-construct targets
    cons_src = np.asarray([X.construct_of(g) for g in src_guide])
    by_t: dict[str, set] = {}
    for t, c in zip(src_tgt, cons_src):
        by_t.setdefault(t, set()).add(c)
    multi_src = sorted([t for t, cs in by_t.items() if len(cs) > 1 and t != X.CTRL])
    multi_h1 = {}
    for t in targets:
        gs = sorted(set(guide[lab == t].tolist()))
        if len(gs) > 1:
            multi_h1[t] = gs
    ka4 = dict(label="KA-4 construct parser on adata_Test.h5ad obs == 273's five multi-construct targets; "
                     "every h1-test barcode in the source with a matching lane suffix",
               multi_from_full_test=multi_src, multi_in_h1test=sorted(multi_h1),
               guide_ids={t: v for t, v in multi_h1.items()},
               all_barcodes_found=ka4_all_found, suffix_equals_batch=ka4_suffix,
               ok=bool(multi_src == sorted(REF_GUIDE) and sorted(multi_h1) == sorted(REF_GUIDE)
                       and ka4_all_found and ka4_suffix))
    P(f"{ka4['label']}: {'PASS' if ka4['ok'] else 'FAIL'}  full-test multi {multi_src}; h1-test multi "
      f"{sorted(multi_h1)}; barcodes found {ka4_all_found}; suffix==batch {ka4_suffix}")
    if not ka4["ok"]:
        FAIL.append(ka4["label"])

    # --- group sums (273's layout)
    is_ctrl = lab == X.CTRL
    ntc_constructs = sorted(set(guide[is_ctrl].tolist()))
    ntc_par = {g: i % 2 for i, g in enumerate(ntc_constructs)}
    multi = multi_h1
    mpos = {t: i for i, t in enumerate(multi)}
    G_T, G_A, G_B, G_C, G_CA, G_CB = 0, 85, 170, 255, 256, 257
    G_GA, G_GB = 258, 258 + len(multi)
    G_GCA, G_GCB = 258 + 2 * len(multi), 258 + 2 * len(multi) + 1
    n_groups = G_GCB + 1

    def groups_full():
        return np.asarray([G_C if is_ctrl[i] else tpos[lab[i]] for i in range(n_cells)])

    def groups_half():
        return np.asarray([(G_CA if parity[i] == 0 else G_CB) if is_ctrl[i]
                           else (G_A if parity[i] == 0 else G_B) + tpos[lab[i]] for i in range(n_cells)])

    def groups_guide():
        out = np.full(n_cells, -1, dtype=np.int64)
        for i in range(n_cells):
            if is_ctrl[i]:
                out[i] = G_GCA if ntc_par[guide[i]] == 0 else G_GCB
            elif lab[i] in multi:
                out[i] = (G_GA if guide[i] == multi[lab[i]][0] else G_GB) + mpos[lab[i]]
        return out

    acc = np.zeros((n_groups, n_genes)); cnt = np.zeros(n_groups, dtype=np.int64)
    for gf in (groups_full, groups_half, groups_guide):
        a, c = X.stream_group_sums(X.H1TEST, gf(), n_groups)
        acc += a; cnt += c
    P(f"  sums: cells/target {cnt[G_T:G_A].min()}..{cnt[G_T:G_A].max()}; controls {cnt[G_C]} = {cnt[G_CA]} "
      f"+ {cnt[G_CB]}; {len(batches)} lanes; {time.time() - t0:.0f}s, peak {L.peak_gib():.2f} GiB")

    # --- KA-3: integer identity
    d_counts = int(np.abs(cnt[G_A:G_B] + cnt[G_B:G_C] - cnt[G_T:G_A]).max())
    d_sums = float(np.abs(acc[G_A:G_B] + acc[G_B:G_C] - acc[G_T:G_A]).max())
    d_ctrl = int(abs(cnt[G_CA] + cnt[G_CB] - cnt[G_C]))
    d_ctrl_s = float(np.abs(acc[G_CA] + acc[G_CB] - acc[G_C]).max())
    ka3 = dict(label="KA-3 lane half A + lane half B == the full row (counts and integer sums)",
               max_abs_count_diff=d_counts, max_abs_sum_diff=d_sums, ctrl_count_diff=d_ctrl,
               ctrl_sum_diff=d_ctrl_s, ok=bool(d_counts == 0 and d_sums == 0.0 and d_ctrl == 0 and d_ctrl_s == 0.0))
    P(f"{ka3['label']}: {'PASS' if ka3['ok'] else 'FAIL'}  counts |d| {d_counts}, sums |d| {d_sums:.1e}, "
      f"controls {d_ctrl} / {d_ctrl_s:.1e}")
    if not ka3["ok"]:
        FAIL.append(ka3["label"])

    # --- 273's residual space
    gpos = {g: i for i, g in enumerate(genes)}
    tcols = np.asarray([gpos[t] for t in targets], dtype=np.int64)
    scored = np.asarray(sorted(set(range(n_genes)) - set(tcols.tolist())), dtype=np.int64)
    assert scored.size == 17995

    def halves(sumA, sumB, cA, cB):
        eA = (L.prof(sumA) - L.prof(cA)[None, :])[:, scored]
        eB = (L.prof(sumB) - L.prof(cB)[None, :])[:, scored]
        return dict(full=L.row_cos(eA, eB),
                    resid=L.row_cos(eA - eA.mean(0, keepdims=True), eB - eB.mean(0, keepdims=True)))

    hb = halves(acc[G_A:G_B], acc[G_B:G_C], acc[G_CA], acc[G_CB])
    med_b = float(np.median(hb["resid"]))
    ka1 = dict(label="KA-1 273's lane-parity 200/200 residual median on h1-test", mine=med_b,
               reference=REF_BATCH_MEDIAN, abs_diff=abs(med_b - REF_BATCH_MEDIAN),
               mean=float(hb["resid"].mean()), n=int(hb["resid"].size),
               ok=bool(med_b == REF_BATCH_MEDIAN))
    P(f"{ka1['label']}: {'PASS' if ka1['ok'] else 'FAIL'}  {med_b:.17f} vs 273 {REF_BATCH_MEDIAN:.17f} "
      f"(|d| {ka1['abs_diff']:.1e}); mean {ka1['mean']:.6f} (273: 0.330950)")
    if not ka1["ok"]:
        FAIL.append(ka1["label"])

    # --- KA-2: 273's guide row (halves_sub: centre on the 85-target FULL-side mean of that side's control)
    gA = np.stack([acc[G_GA + mpos[t]] for t in multi])
    gB = np.stack([acc[G_GB + mpos[t]] for t in multi])
    cA, cB = acc[G_GCA], acc[G_GCB]
    eA = (L.prof(gA) - L.prof(cA)[None, :])[:, scored]
    eB = (L.prof(gB) - L.prof(cB)[None, :])[:, scored]
    mA = (L.prof(acc[G_T:G_A]) - L.prof(cA)[None, :])[:, scored].mean(0, keepdims=True)
    mB = (L.prof(acc[G_T:G_A]) - L.prof(cB)[None, :])[:, scored].mean(0, keepdims=True)
    g_resid = L.row_cos(eA - mA, eB - mB)
    g_full = L.row_cos(eA, eB)
    mt = list(multi)
    per = {t: dict(mine=float(g_resid[i]), ref=REF_GUIDE[t][0], abs_diff=abs(float(g_resid[i]) - REF_GUIDE[t][0]),
                   mine_full=float(g_full[i]), ref_full=REF_GUIDE_FULL[t],
                   cells=[int(cnt[G_GA + mpos[t]]), int(cnt[G_GB + mpos[t]])],
                   ref_cells=[REF_GUIDE[t][1], REF_GUIDE[t][2]]) for i, t in enumerate(mt)}
    med_g = float(np.median(g_resid))
    ok2 = (med_g == REF_GUIDE_MEDIAN
           and all(v["mine"] == v["ref"] and v["mine_full"] == v["ref_full"] and v["cells"] == v["ref_cells"]
                   for v in per.values()))
    ka2 = dict(label="KA-2 273's 5-target construct row (resid, full, cell counts) from h1-test",
               median=med_g, reference_median=REF_GUIDE_MEDIAN, per_target=per, ok=bool(ok2),
               max_abs_diff=float(max(v["abs_diff"] for v in per.values())))
    P(f"{ka2['label']}: {'PASS' if ka2['ok'] else 'FAIL'}  median {med_g:.17f} vs 273 {REF_GUIDE_MEDIAN:.17f}; "
      f"max per-target |d| {ka2['max_abs_diff']:.1e}")
    for t in mt:
        v = per[t]
        P(f"    {t:7s} resid {v['mine']:+.6f} (273 {v['ref']:+.6f})  full {v['mine_full']:+.6f} "
          f"(273 {v['ref_full']:+.6f})  cells {v['cells']} (273 {v['ref_cells']})")
    if not ka2["ok"]:
        FAIL.append(ka2["label"])

    out = dict(ka=[ka1, ka2, ka3, ka4], failures=FAIL, seconds=round(time.time() - t0, 1),
               peak_rss_gib=L.peak_gib(), h1_sha256=L.sha256_file(X.H1TEST),
               per_target_pdstat_source="lib285.per_target sha256 " + L.PER_TARGET_SHA)
    L.Logger.dump(P, out, X.WORK / "ka293.json")
    P(f"\ndone: {len(FAIL)} KA failures, {out['seconds']}s, peak {out['peak_rss_gib']:.2f} GiB")
    for f in FAIL:
        P(f"  *** {f}")
    P.save()
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
