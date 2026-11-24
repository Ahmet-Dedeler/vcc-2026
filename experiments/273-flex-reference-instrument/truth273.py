"""273 item 1-2 -- the h1-test TRUTH in the scorer's own space, the scorer's per-target `pds` statistic validated on a
cloud-scored prediction (KA-2), 159's reliability rows reproduced (KA-1), and the reference's split-half reliability
by experimental unit (batch / guide / random cells).

    PYTHONDONTWRITEBYTECODE=1 uv run python truth273.py

Reads (read-only): data/harness/h1-test.h5ad (0.68 GB, streamed in CSR row blocks), data/external/vcc2025/test/
adata_Test.h5ad (obs ONLY, for guide_id; the batch is the barcode suffix and is checked against the source column),
experiments/159-transfer-ceiling/work/shards-h1test.npz + reliability.json, experiments/095-clean-harness/.../pred.h5ad
(2.05 GB, streamed) + its run/results.csv, experiments/221-cloud-scorer/cloud-output/*/ka095/run/results.csv.
Writes work/truth273.npz, work/reliability273.{json,log}.

KNOWN-ANSWER ROWS (exit 1 on failure):
  KA-1  159's `reliability.curve` on 159's shards reproduces reliability.json's h1test rows (all 10 depths, 4 keys) at
        1e-12: full 0.3178 / resid 0.3300 at 200 cells per side.
  KA-2  the scorer's `pds_cosine`, recomputed here from the emitted counts of the 095 artifact (bulk log1p at TS
        50,000, the REAL control on both sides, the 85 target-gene columns removed, cosine, midrank ties, denominator
        n - 1), equals the cloud scorer's per-target values (221 ka095 results.csv == 095's own results.csv) and the
        mean 0.8655462184873949, to the digit.
"""
from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp
from anndata.io import read_elem

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib273 as X  # noqa: E402
import reliability as R159  # noqa: E402  (159's module; reads 159's shards read-only)

H1 = X.ROOT / "data/harness/h1-test.h5ad"
SRC = X.ROOT / "data/external/vcc2025/test/adata_Test.h5ad"
PRED095 = X.E095 / "pred.h5ad"
CHUNK = 4000
LOG: list[str] = []
FAIL: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def labels_of(h, path):
    o = h[path]
    if isinstance(o, h5py.Group) and "categories" in o:
        cats = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]])
        return cats[o["codes"][:]]
    v = o["values"][:] if isinstance(o, h5py.Group) else o[:]
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in v])


def stream_sums(path: Path, row_groups: np.ndarray, n_groups: int) -> tuple[np.ndarray, np.ndarray]:
    """Exact per-group count sums over the rows (row_groups[i] = group of row i, -1 = skip), streamed."""
    with h5py.File(path, "r") as h:
        n, g = (int(x) for x in h["X"].attrs["shape"])
        indptr = h["X/indptr"][:]
        acc = np.zeros((n_groups, g), dtype=np.float64)
        for r0 in range(0, n, CHUNK):
            r1 = min(r0 + CHUNK, n)
            s, e = int(indptr[r0]), int(indptr[r1])
            Xb = sp.csr_matrix((h["X/data"][s:e].astype(np.float64), h["X/indices"][s:e], indptr[r0:r1 + 1] - indptr[r0]),
                               shape=(r1 - r0, g))
            grp = row_groups[r0:r1]
            keep = grp >= 0
            if keep.any():
                oh = sp.csr_matrix((np.ones(int(keep.sum())), (grp[keep], np.flatnonzero(keep))), shape=(n_groups, r1 - r0))
                acc += np.asarray((oh @ Xb).todense())
    cnt = np.bincount(row_groups[row_groups >= 0], minlength=n_groups)
    return acc, cnt


def prof(P_):
    P_ = np.asarray(P_, dtype=np.float64)
    if P_.ndim == 1:
        return np.log1p(X.TS * P_ / P_.sum())
    tot = P_.sum(1, keepdims=True)
    return np.log1p(X.TS * P_ / np.where(tot > 0, tot, 1.0))


def rel_stats(r: np.ndarray, label: str) -> dict:
    r = np.asarray(r, dtype=np.float64)
    ok = np.isfinite(r)
    v = r[ok]
    sb = 2 * v / (1 + v)
    d = dict(label=label, n=int(ok.sum()), mean=float(v.mean()), median=float(np.median(v)), p10=float(np.percentile(v, 10)),
             p25=float(np.percentile(v, 25)), p75=float(np.percentile(v, 75)), p90=float(np.percentile(v, 90)),
             share_le0=float((v <= 0).mean()), share_ge005=float((v >= 0.05).mean()), share_ge020=float((v >= 0.2).mean()),
             sb_of_median=float(2 * np.median(v) / (1 + np.median(v))), sb_of_mean=float(2 * v.mean() / (1 + v.mean())),
             median_sb=float(np.median(sb)), mean_sb=float(sb.mean()))
    P(f"  {label:44s} n {d['n']:3d}  median {d['median']:+.4f}  mean {d['mean']:+.4f}  p10 {d['p10']:+.4f}  p90 {d['p90']:+.4f}  "
      f"share<=0 {d['share_le0']:.3f}  >=0.05 {d['share_ge005']:.3f}  SB(median) {d['sb_of_median']:.3f}")
    return d


def read_results_csv(p: Path, metric: str = "pds_cosine") -> dict[str, float]:
    out = {}
    with open(p) as fh:
        for row in csv.DictReader(fh):
            if row["metric"] == metric:
                out[row["perturbation"]] = float(row["value"])
    return out


def main() -> None:
    t0 = time.time()
    P(f"273 truth273.py -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}; h1-test sha256 {X.sha256_file(H1)[:12]}")

    # ------------------------------------------------------------------------------------------ KA-1: 159's rows
    ref = json.load(open(X.E159 / "work/reliability.json"))["h1test"]
    d159 = R159.load("h1test")
    rows159, meta159 = R159.curve(d159, verbose=False)
    ms = [r["m"] for r in rows159]
    worst = 0.0
    for key in R159.KEYS:
        mine = np.array([r["per_target"][key].mean() for r in rows159])
        theirs = np.asarray(ref[key]["r_mean"], dtype=np.float64)
        worst = max(worst, float(np.abs(mine - theirs).max()))
    ka1 = dict(label="KA-1 159 reliability rows (4 keys x 10 depths) reproduced from 159's shards", max_abs_diff=worst, ok=bool(worst < 1e-12),
               resid_disjoint_m200=float(rows159[-1]["per_target"]["resid_disjoint"].mean()),
               full_disjoint_m200=float(rows159[-1]["per_target"]["full_disjoint"].mean()), m=ms,
               resid_disjoint_r_mean=[float(r["per_target"]["resid_disjoint"].mean()) for r in rows159],
               full_disjoint_r_mean=[float(r["per_target"]["full_disjoint"].mean()) for r in rows159])
    P(f"{ka1['label']}: {'PASS' if ka1['ok'] else 'FAIL'} max|d| {worst:.2e}; resid/disjoint at m=200 {ka1['resid_disjoint_m200']:.6f} "
      f"(159: 0.3300), full/disjoint {ka1['full_disjoint_m200']:.6f} (159: 0.3178); SB(0.3300) = {2*0.33/(1.33):.4f}; "
      f"159's power-law rho_400 {ref['resid_disjoint']['power']['rho_400']:.4f}, linear {ref['resid_disjoint']['linear']['rho_400']:.4f}")
    if not ka1["ok"]:
        FAIL.append(ka1["label"])
    # per-target arrays from 159's shards: averaged over 24 pairings at m = 80 (267's Orion depth) and m = 200
    r159_80 = np.asarray(rows159[3]["per_target"]["resid_disjoint"]); assert rows159[3]["m"] == 80.0
    r159_200 = np.asarray(rows159[-1]["per_target"]["resid_disjoint"])
    r159_full_200 = np.asarray(rows159[-1]["per_target"]["full_disjoint"])
    targets159 = d159["targets"]

    # ------------------------------------------------------------------------------------------ h1-test: stream
    with h5py.File(H1, "r") as h:
        obs = read_elem(h["obs"])
        h1_genes = labels_of(h, "var/_index")
        n_cells, n_genes = (int(x) for x in h["X"].attrs["shape"])
    lab = obs["target_gene"].to_numpy().astype(str)
    barcodes = obs.index.to_numpy().astype(str)
    targets = np.asarray(sorted(set(lab) - {X.CTRL}))
    assert targets.size == 85 and np.array_equal(targets, targets159)
    tpos = {t: i for i, t in enumerate(targets)}
    # batch = barcode suffix (checked against the source's batch column); guide_id joined from the source obs
    with h5py.File(SRC, "r") as h:
        src_idx = labels_of(h, "obs/_index")
        src_batch = labels_of(h, "obs/batch")
        src_guide = labels_of(h, "obs/guide_id")
    pos = {b: i for i, b in enumerate(src_idx)}
    hit = np.asarray([pos[b] for b in barcodes], dtype=np.int64)
    suffix = np.asarray([b.split("-", 1)[1] for b in barcodes])
    assert np.array_equal(suffix, src_batch[hit]), "barcode suffix != source batch column"
    guide = src_guide[hit]
    batches = sorted(set(suffix))
    bidx = {b: i for i, b in enumerate(batches)}
    parity = np.asarray([bidx[b] % 2 for b in suffix])          # half A = even batch index, half B = odd
    P(f"h1-test: {n_cells} cells x {n_genes} genes, {targets.size} targets, {len(batches)} Flex batches; every barcode found in the 2025 test "
      f"source; batch split by parity of the sorted batch index ({(parity == 0).sum()} / {(parity == 1).sum()} cells)")
    # groups: 0..84 target full; 85..169 target A; 170..254 target B; 255 ctrl full; 256 ctrl A; 257 ctrl B
    # guide split (multi-construct targets only): constructs sorted, P1 -> A, P2 -> B; controls by NTC construct parity
    is_ctrl = lab == X.CTRL
    ntc_constructs = sorted(set(guide[is_ctrl]))
    ntc_par = {g: i % 2 for i, g in enumerate(ntc_constructs)}
    multi = {}
    for t in targets:
        gs = sorted(set(guide[lab == t]))
        if len(gs) > 1:
            multi[t] = gs
    G_T, G_A, G_B, G_C, G_CA, G_CB = 0, 85, 170, 255, 256, 257
    G_GA, G_GB, G_GCA, G_GCB = 258, 258 + len(multi), 258 + 2 * len(multi), 258 + 2 * len(multi) + 1
    n_groups = G_GCB + 1
    mpos = {t: i for i, t in enumerate(multi)}
    # one pass per group set (three cheap passes over a 0.68 GB file)
    def groups_full(): return np.asarray([G_C if is_ctrl[i] else tpos[lab[i]] for i in range(n_cells)])
    def groups_half():
        return np.asarray([(G_CA if parity[i] == 0 else G_CB) if is_ctrl[i] else (G_A if parity[i] == 0 else G_B) + tpos[lab[i]]
                           for i in range(n_cells)])
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
        a, c = stream_sums(H1, gf(), n_groups)
        acc += a; cnt += c
    P(f"  sums: cells per target {cnt[G_T:G_A].min()}..{cnt[G_T:G_A].max()}; batch halves A {cnt[G_A:G_B].min()}..{cnt[G_A:G_B].max()} / "
      f"B {cnt[G_B:G_C].min()}..{cnt[G_B:G_C].max()} (median {np.median(cnt[G_A:G_B]):.0f} / {np.median(cnt[G_B:G_C]):.0f}); controls "
      f"{cnt[G_C]} = {cnt[G_CA]} + {cnt[G_CB]}; guide split on {len(multi)} targets {dict((t, [int(cnt[G_GA + mpos[t]]), int(cnt[G_GB + mpos[t]])]) for t in multi)}, "
      f"NTC constructs {len(ntc_constructs)} -> {cnt[G_GCA]} / {cnt[G_GCB]} control cells; {time.time() - t0:.0f}s, peak {X.peak_gib():.2f} GiB")
    assert np.all(cnt[G_T:G_A] == 400)

    # ------------------------------------------------------------------------------------------ the scorer's truth
    b_ctrl = prof(acc[G_C])
    truth = prof(acc[G_T:G_A]) - b_ctrl[None, :]
    gpos = {g: i for i, g in enumerate(h1_genes)}
    tcols = np.asarray([gpos[t] for t in targets], dtype=np.int64)       # every target gene is on the h1 axis
    scored = np.asarray(sorted(set(range(n_genes)) - set(tcols.tolist())), dtype=np.int64)
    expr = np.expm1(b_ctrl)
    cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
    gate = (cpm > X.MIN_CPM) & (b_ctrl > 0)
    P(f"  truth: {truth.shape}, scored columns {scored.size} (= 18,080 - 85), render gate {int(gate.sum())} genes at cpm > {X.MIN_CPM}; "
      f"median truth norm on scored {np.median(np.linalg.norm(truth[:, scored], axis=1)):.4f}")

    # ------------------------------------------------------------------------------------------ KA-2: the statistic
    with h5py.File(PRED095, "r") as h:
        pobs = read_elem(h["obs"])
        pgenes = labels_of(h, "var/_index")
        pn = int(h["X"].attrs["shape"][0])
    assert np.array_equal(pgenes, h1_genes), "095 artifact gene axis != h1-test's"
    plab = pobs["target_gene"].to_numpy().astype(str)
    pgroups = np.asarray([tpos.get(l, -1) for l in plab])                 # the prediction's own controls are IGNORED (control_source real)
    psum, pcnt = stream_sums(PRED095, pgroups, 85)
    P(f"  095 artifact: {pn} rows, {int((pgroups >= 0).sum())} emitted cells over 85 targets ({pcnt.min()}..{pcnt.max()} each); {time.time() - t0:.0f}s, "
      f"peak {X.peak_gib():.2f} GiB")
    pred_eff = prof(psum) - b_ctrl[None, :]
    pds_loc, cos_loc = X.pds_own(pred_eff[:, scored], truth[:, scored])
    ka095 = read_results_csv(X.E221 / "ladder6-20260904/ka095/run/results.csv")
    r095 = read_results_csv(X.E095 / "run/results.csv")
    agg = 0.8655462184873949
    d_cloud = float(max(abs(pds_loc[i] - ka095[t]) for i, t in enumerate(targets)))
    d_local = float(max(abs(pds_loc[i] - r095[t]) for i, t in enumerate(targets)))
    d_mean = abs(float(pds_loc.mean()) - agg)
    ka2 = dict(label="KA-2 the scorer's per-target pds_cosine recomputed from the 095 artifact's emitted counts == cell-eval2's (cloud ka095 and the 2026-08-27 local run)",
               max_abs_diff_per_target_vs_cloud=d_cloud, max_abs_diff_per_target_vs_local_run=d_local, mean_local=float(pds_loc.mean()),
               mean_scorer=agg, abs_diff_mean=d_mean, ok=bool(d_cloud < 1e-12 and d_mean < 1e-12), n=85,
               scorer_recipe="bulk log1p(TS*P/sum P) TS=50000; effect vs the REAL control profile on both sides; all 85 target-gene columns removed "
                             "(exclusion_scope panel); cosine distance; midrank ties; PDS_p = 1 - rank/(n-1)")
    P(f"{ka2['label']}: {'PASS' if ka2['ok'] else 'FAIL'} max|d| per target {d_cloud:.2e} (vs 095's own run {d_local:.2e}); mean {pds_loc.mean():.16f} vs "
      f"scorer {agg:.16f} (|d| {d_mean:.2e})")
    if not ka2["ok"]:
        FAIL.append(ka2["label"])

    # ------------------------------------------------------------------------------------------ reliability by unit
    P("\nReliability of the reference by experimental unit (cos(e_A, e_B) per target on the 17,995 scored columns; resid = minus the mean over the 85 targets of THAT side)")
    rel: dict[str, dict] = {}
    arrays: dict[str, np.ndarray] = {}

    def halves(sumA, sumB, cA, cB, label):
        eA, eB = prof(sumA) - prof(cA)[None, :], prof(sumB) - prof(cB)[None, :]
        eA, eB = eA[:, scored], eB[:, scored]
        rA, rB = eA - eA.mean(0, keepdims=True), eB - eB.mean(0, keepdims=True)
        return dict(full=X.row_cos(eA, eB), resid=X.row_cos(rA, rB))

    hb = halves(acc[G_A:G_B], acc[G_B:G_C], acc[G_CA], acc[G_CB], "batch/disjoint")
    hbs = halves(acc[G_A:G_B], acc[G_B:G_C], acc[G_C], acc[G_C], "batch/shared")
    arrays["batch|resid"], arrays["batch|full"] = hb["resid"], hb["full"]
    arrays["batch_shared|resid"], arrays["batch_shared|full"] = hbs["resid"], hbs["full"]
    rel["batch|resid|disjoint"] = rel_stats(hb["resid"], "batch parity, resid, disjoint controls (200/200)")
    rel["batch|full|disjoint"] = rel_stats(hb["full"], "batch parity, full, disjoint controls")
    rel["batch|resid|shared"] = rel_stats(hbs["resid"], "batch parity, resid, shared control")
    rel["batch|full|shared"] = rel_stats(hbs["full"], "batch parity, full, shared control")
    # random cells: 159's shards (20 shards of 20 cells; disjoint control halves), one fixed pairing at 200/200 and the
    # 24-pairing averages at 80 and 200 per side
    z = d159
    mask159 = z["mask"]                    # 159 scored columns: its own genes minus the 85 targets (== `scored` up to order)
    assert np.array_equal(np.asarray([str(g) for g in z["genes"]]), h1_genes) and np.array_equal(mask159, scored)
    tsum, csum, kt, kc = z["tsum"], z["csum"], z["kt"], z["kc"]
    hr = halves(tsum[:, :kt // 2].sum(1), tsum[:, kt // 2:].sum(1), csum[:kc // 2].sum(0), csum[kc // 2:].sum(0), "random/fixed")
    arrays["rand_fixed|resid"], arrays["rand_fixed|full"] = hr["resid"], hr["full"]
    rel["rand_fixed|resid|disjoint"] = rel_stats(hr["resid"], "random cells (159 shards 0-9 vs 10-19), resid, disjoint")
    rel["rand_fixed|full|disjoint"] = rel_stats(hr["full"], "random cells (fixed pairing), full, disjoint")
    arrays["rand24_m200|resid"], arrays["rand24_m80|resid"], arrays["rand24_m200|full"] = r159_200, r159_80, r159_full_200
    rel["rand24_m200|resid"] = rel_stats(r159_200, "random cells, resid, mean of 24 pairings, 200/200 (159)")
    rel["rand24_m80|resid"] = rel_stats(r159_80, "random cells, resid, mean of 24 pairings, 80/80 (159; 267's Orion depth)")
    rel["rand24_m200|full"] = rel_stats(r159_full_200, "random cells, full, mean of 24 pairings, 200/200 (159)")
    # guide split: the multi-construct targets only (P1 half vs P2 half; controls by NTC construct parity), and the SAME
    # targets' batch-split values for the comparison
    gA = np.stack([acc[G_GA + mpos[t]] for t in multi]); gB = np.stack([acc[G_GB + mpos[t]] for t in multi])
    def halves_sub(sumA, sumB, cA, cB):
        # resid centring needs the panel mean: use ALL 85 targets' batch-half means of that side? No -- a guide half
        # exists for 5 targets only; centre on the 85-target mean of the corresponding FULL side (the depth bias is
        # the same order at 150-360 cells), reported as such.
        eA, eB = (prof(sumA) - prof(cA)[None, :])[:, scored], (prof(sumB) - prof(cB)[None, :])[:, scored]
        mA = (prof(acc[G_T:G_A]) - prof(cA)[None, :])[:, scored].mean(0, keepdims=True)
        mB = (prof(acc[G_T:G_A]) - prof(cB)[None, :])[:, scored].mean(0, keepdims=True)
        return dict(full=X.row_cos(eA, eB), resid=X.row_cos(eA - mA, eB - mB))
    hg = halves_sub(gA, gB, acc[G_GCA], acc[G_GCB])
    gi = np.asarray([tpos[t] for t in multi])
    rel["guide|resid"] = rel_stats(hg["resid"], f"GUIDE construct P1 vs P2 ({len(multi)} targets), resid")
    rel["guide|full"] = rel_stats(hg["full"], f"guide construct P1 vs P2 ({len(multi)} targets), full")
    rel["guide_targets_batch|resid"] = rel_stats(hb["resid"][gi], "  the same targets, batch split, resid")
    rel["guide_targets_rand|resid"] = rel_stats(r159_200[gi], "  the same targets, random 24-pairing 200/200, resid")
    rel["guide_cells"] = {t: dict(constructs=multi[t], cells=[int(cnt[G_GA + mpos[t]]), int(cnt[G_GB + mpos[t]])],
                                  resid=float(hg["resid"][i]), full=float(hg["full"][i]), batch_resid=float(hb["resid"][tpos[t]]))
                          for i, t in enumerate(multi)}
    # cell-count dependence within the batch split
    ncA, ncB = cnt[G_A:G_B], cnt[G_B:G_C]
    from scipy.stats import spearmanr
    rel["rho_resid_vs_min_cells_batch"] = float(spearmanr(hb["resid"], np.minimum(ncA, ncB)).statistic)
    rel["batch_vs_random_fixed_spearman"] = float(spearmanr(hb["resid"], hr["resid"]).statistic)
    rel["batch_vs_random24_spearman"] = float(spearmanr(hb["resid"], r159_200).statistic)
    P(f"  Spearman(batch resid, random-fixed resid) {rel['batch_vs_random_fixed_spearman']:+.3f}; (batch resid, 24-pairing random) {rel['batch_vs_random24_spearman']:+.3f}; "
      f"rho(batch resid, min cells per half) {rel['rho_resid_vs_min_cells_batch']:+.3f}")
    # 267's Orion / K562 batch-split numbers (R267 section 3; R271 D1) for the side-by-side table
    rel["reference_267"] = {
        "k562_gwps batch (81 cells/half)": dict(median=0.0471, mean=0.0654, p10=0.0079, p90=0.1619, share_le0=0.056, share_ge005=0.468, sb_median=0.090),
        "orion_hct116 batch (70 cells/half)": dict(median=0.0276, mean=0.0356, p10=-0.0494, p90=0.1191, share_le0=0.321, share_ge005=0.338, sb_median=0.054),
        "orion_hek293t batch (95 cells/half)": dict(median=0.0282, mean=0.0417, p10=-0.1218, p90=0.2140, share_le0=0.344, share_ge005=0.415, sb_median=0.055),
        "k562_gwps random (81)": dict(median=0.0450, share_le0=0.059), "orion_hct116 random (72)": dict(median=0.0293, share_le0=0.297),
        "orion_hek293t random (95)": dict(median=0.0706, share_le0=0.146),
        "159 orion-hct116-thinned random 75/side (resid mean)": 0.0286, "159 orion-hek293t-thinned random 75/side (resid mean)": 0.0637,
        "159 h1-test random 80/side (resid mean)": float(rows159[3]["per_target"]["resid_disjoint"].mean()),
        "159 h1-test random 200/side (resid mean)": float(rows159[-1]["per_target"]["resid_disjoint"].mean())}

    np.savez_compressed(X.WORK / "truth273.npz", targets=targets, h1_genes=h1_genes, b_ctrl=b_ctrl, truth=truth.astype(np.float64), scored=scored,
                        gate=gate, tcols=tcols, h1_sha256=X.sha256_file(H1), cells_A=ncA, cells_B=ncB, ctrl_cells=np.asarray([cnt[G_C], cnt[G_CA], cnt[G_CB]]),
                        b_ctrl_A=prof(acc[G_CA]), b_ctrl_B=prof(acc[G_CB]), pds095_local=pds_loc, cos095_local=cos_loc,
                        pds095_cloud=np.asarray([ka095[t] for t in targets]), **{f"rel|{k}": v for k, v in arrays.items()})
    out = dict(ka=[ka1, ka2], failures=FAIL, n_cells=n_cells, n_genes=n_genes, n_targets=int(targets.size), n_batches=len(batches),
               scored_columns=int(scored.size), render_gate_genes=int(gate.sum()), controls=dict(full=int(cnt[G_C]), A=int(cnt[G_CA]), B=int(cnt[G_CB])),
               cells_per_half=dict(A=dict(min=int(ncA.min()), median=float(np.median(ncA)), max=int(ncA.max())),
                                   B=dict(min=int(ncB.min()), median=float(np.median(ncB)), max=int(ncB.max()))),
               multi_construct_targets={t: v for t, v in multi.items()}, reliability=rel, inputs=dict(
                   h1_test=str(H1.relative_to(X.ROOT)), h1_sha256=X.sha256_file(H1), source_obs=str(SRC.relative_to(X.ROOT)),
                   shards159=str((X.E159 / "work/shards-h1test.npz").relative_to(X.ROOT)), pred095=str(PRED095.relative_to(X.ROOT)),
                   pred095_sha256=X.sha256_file(PRED095)), seconds=round(time.time() - t0, 1), peak_rss_gib=X.peak_gib())
    (X.WORK / "reliability273.json").write_text(json.dumps(out, indent=1, default=float))
    (X.WORK / "reliability273.log").write_text("\n".join(LOG) + "\n")
    P(f"\ndone: {len(FAIL)} KA failures, {out['seconds']}s, peak RSS {out['peak_rss_gib']:.2f} GiB")
    for f in FAIL:
        P(f"  *** {f}")
    if FAIL:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
