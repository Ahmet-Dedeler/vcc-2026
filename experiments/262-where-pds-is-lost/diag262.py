"""262 -- WHERE does `cd4svd5` lose `pds`, target by target, on one PLOLO fold.  MEASUREMENT ONLY.

    uv run python diag262.py --fold <k562_gwps | orion_hct116 | orion_hek293t | kolf2>

Nothing mathematical is new.  The construction is 257/261's, imported: `lib257` -> `lib246` (162's PLOLO,
218's residuals and `union_rows`, the shipped stack, `deflate`, `compose`, `score`).  The incumbent
`cd4svd5` is 261's `cd4svd5_parts`: CD4 deflated by `lib246.deflate(cd4, "svd", 5)`, RE-SOLVED by
`plolo.solve_fraction` at f = 0.34, merged with the H1-pair union at its shipped carrying-row fraction
(g = 0.28), one add, one preserve.  `h1union` is `lib246.stack_parts` (the un-deflated stack).

KNOWN-ANSWER ROWS (the run exits 1 on any failure):
  KA-1  `cd4svd5` endpoints (cos_agg, a23_cos, pds) == 257 `sweep257-<fold>.json` row `k5_s1` at 1e-9,
        and its per-target `pds` / cosine arrays == 257 `pertarget257-<fold>.npz` `k5_s1|pds` / `k5_s1|cos`
        AND 261 `pertarget261-<fold>.npz` `cd4svd5|pds` at max|d| == 0.0 (bit-identical).
  KA-2  `h1union` (STACK) likewise vs 257 row `STACK` and `STACK|pds` / `STACK|cos`.
  KA-3  this file's own per-target rank statistic, recomputed from the (prediction x truth) cosine matrix
        it saves, == the instrument's `_pds_per_target_a23` at max|d| < 1e-12 on every arm.

WHAT IS SAVED (`work/diag262-<fold>.npz`, `work/diag262-<fold>.json`), per target:
  * per-arm per-target `pds` and own cosine, delivered (a23, the instrument's) and at effect level
    (un-rendered), for arms cd4svd5 / h1union / BASE (no residual) / CD4ONLY (deflated CD4 at its resolved
    fraction, no union) / UNIONONLY (the union at its shipped fraction, no CD4);
  * the geometry of the loss: competitor cosine max / median / mean, beat matrix B[p, q] (how many of the
    three challenge frames rank truth q above truth p for prediction p), prediction-prediction and
    truth-truth nearest-neighbour and mean cosines, delivered prediction norm, truth norm;
  * structural covariates: which whitened sources carry the target (and PerturbMatch), CD4 condition count
    and row norm (complete rows), the per-row energy the k = 5 deflation removed, union carriers and row
    norm, the `k562_essential` flag, whitened-base row norm and its cosine to the truth, the merged
    residual's norm, its cosine to the truth and to the base, and the residual's share of the row.
Read-only everywhere except `experiments/262-where-pds-is-lost/work/`; the fold caches under `work/plolo/`
are byte copies of 257's (sha-listed in `work/plolo/SHA256SUMS.txt`).
"""
from __future__ import annotations

import argparse
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
E257 = ROOT / "experiments/257-soft-deflation-price"
E261 = ROOT / "experiments/261-union-deflation"
sys.path.insert(0, str(E257))
import lib257 as X  # noqa: E402  (imports lib246 as X.L -> plolo / resid218 / compose218)
import plolo  # noqa: E402

L = X.L
# lib246 pins plolo.WORK to 246's directory and lib257 to 257's; pin it to OURS so every cache read comes
# from this experiment's byte copies and nothing outside it is ever touched.
plolo.WORK = WORK / "plolo"
L.WORK = WORK
X.WORK = WORK
E062, E064 = plolo.E062, plolo.E064

TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
ESSENTIAL = {ln.strip() for ln in (ROOT / "data/donors/k562_essential/perts.txt").read_text().splitlines()
             if ln.strip()}


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def rank_stat(S):
    """The instrument's per-target rank statistic from a (prediction x truth) cosine matrix
    (`resid218._per_target`, verbatim arithmetic): the fraction of the other truths that are FARTHER from
    prediction p than p's own truth, ties at half."""
    n = S.shape[0]
    own = np.diag(S).copy()
    return ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1), own


def geometry(P, T):
    """Everything the loss is made of, for one frame: prediction rows P and truth rows T on the scored
    columns.  Returns per-target arrays and the beat matrix."""
    A, B = unit(P), unit(T)
    S = A @ B.T                                   # S[p, q] = cos(pred_p, truth_q)
    n = S.shape[0]
    own = np.diag(S).copy()
    off = ~np.eye(n, dtype=bool)
    Sm = np.where(off, S, np.nan)
    beats = (S > own[:, None]) & off               # truth q ranked ABOVE p's own truth for prediction p
    ties = (S == own[:, None]) & off
    Gp = A @ A.T
    Gt = B @ B.T
    Gpm = np.where(off, Gp, np.nan)
    Gtm = np.where(off, Gt, np.nan)
    pds, _ = rank_stat(S)
    return dict(
        pds=pds, own=own,
        comp_max=np.nanmax(Sm, axis=1), comp_med=np.nanmedian(Sm, axis=1), comp_mean=np.nanmean(Sm, axis=1),
        n_beat=beats.sum(1).astype(np.float64), n_tie=ties.sum(1).astype(np.float64),
        pred_nn=np.nanmax(Gpm, axis=1), pred_mean=np.nanmean(Gpm, axis=1),
        truth_nn=np.nanmax(Gtm, axis=1), truth_mean=np.nanmean(Gtm, axis=1),
        pred_norm=np.linalg.norm(P, axis=1), truth_norm=np.linalg.norm(T, axis=1),
        # the prediction's cosine to the truth's nearest-neighbour truth: does it point at the
        # RIGHT NEIGHBOURHOOD but the wrong member of it?
        own_vs_nn=own - np.nanmax(Sm, axis=1),
    ), beats.astype(np.int8), ties.astype(np.int8), S


def frame_geometry(final, fold, truth, mass="minnorm", render=True):
    """The instrument's a23 read on the three challenge frames (mean over A/B/C, as `score_effects`
    averages), or the effect-level read (render=False: the same rows un-rendered, one frame)."""
    sc = fold["scored"]
    acc, B_sum, T_sum = [], None, None
    S_by_ctx = {}
    ctxs = ("A", "B", "C") if render else ("effect",)
    for ctx in ctxs:
        if render:
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            P = plolo.deliver(final, b, gate, mass)[:, sc]
        else:
            P = np.asarray(final, dtype=np.float64)[:, sc]
        g, beats, ties, S = geometry(P, truth)
        acc.append(g)
        B_sum = beats.astype(np.int16) if B_sum is None else B_sum + beats
        T_sum = ties.astype(np.int16) if T_sum is None else T_sum + ties
        S_by_ctx[ctx] = S
    out = {k: np.mean([a[k] for a in acc], axis=0) for k in acc[0]}
    return out, B_sum, T_sum, S_by_ctx


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    failures: list[str] = []

    ref257 = {r["arm"]: r for r in json.loads((E257 / f"work/sweep257-{held}.json").read_text())["rows"]}
    pt257 = np.load(E257 / f"work/pertarget257-{held}.npz")
    pt261 = np.load(E261 / f"work/pertarget261-{held}.npz")

    st = L.stage(held)
    fold, targets = st["fold"], list(st["targets"])
    truth = np.asarray(st["truth"], dtype=np.float64)          # (n, scored)
    n = len(targets)
    print(f"[{held}] staged {n} targets; stack frac_cd4 {st['frac_cd4']:.6f} frac_x {st['frac_x']:.6f} "
          f"({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)

    # ------------------------------------------------------------------ the two objects + components
    d, aud_defl = L.deflate(st["cd4"], "svd", 5)
    fr, pu = L.solve_fraction(st["whitened"], d, L.CD4_F)
    arms = {
        "cd4svd5": [(d, float(fr)), (st["union"], st["frac_x"])],
        "h1union": L.stack_parts(st),
        "BASE": [],
        "CD4ONLY": [(d, float(fr))],
        "UNIONONLY": [(st["union"], st["frac_x"])],
    }
    ka_ref = {"cd4svd5": ("k5_s1", "k5_s1"), "h1union": ("STACK", "STACK")}

    arrays: dict[str, np.ndarray] = {}
    summary: dict[str, dict] = {}
    finals: dict[str, np.ndarray] = {}
    for name, parts in arms.items():
        t1 = time.time()
        final = L.compose(st, parts)
        finals[name] = final
        e = L.score(st, final)                                  # the instrument (a23 + a15 tags)
        pds_inst = np.asarray(e["_pds_per_target_a23"], dtype=np.float64)
        cos_inst = np.asarray(e["_cos_per_target_a23"], dtype=np.float64)
        g, B, Tt, S_ctx = frame_geometry(final, fold, truth)     # delivered, three frames
        ge, Be, _, _ = frame_geometry(final, fold, truth, render=False)   # effect level
        # KA-3: my rank statistic == the instrument's, bit for bit
        d3 = float(np.abs(g["pds"] - pds_inst).max())
        d3c = float(np.abs(g["own"] - cos_inst).max())
        if not (d3 < 1e-12 and d3c < 1e-12):
            failures.append(f"KA-3 {name}: per-target pds max|d| {d3:.3e}, cos {d3c:.3e}")
        row = dict(arm=name, **{k: float(v) for k, v in e.items() if not k.startswith("_")},
                   ka3_max_abs_diff=max(d3, d3c), seconds=round(time.time() - t1, 1))
        if name in ka_ref:
            r257, key = ka_ref[name]
            dd = max(abs(float(ref257[r257][k]) - float(e[k])) for k in KEYS)
            dp = float(np.abs(np.asarray(pt257[f"{key}|pds"]) - pds_inst).max())
            dc = float(np.abs(np.asarray(pt257[f"{key}|cos"]) - cos_inst).max())
            row["ka"] = dict(ref=f"257 {r257}", endpoints_max_abs_diff=dd, pertarget_pds_max_abs_diff=dp,
                             pertarget_cos_max_abs_diff=dc)
            if name == "cd4svd5":
                d61 = float(np.abs(np.asarray(pt261["cd4svd5|pds"]) - pds_inst).max())
                row["ka"]["pertarget_pds_vs_261_max_abs_diff"] = d61
                if d61 != 0.0:
                    failures.append(f"KA-1 cd4svd5 vs 261 pertarget: {d61:.3e}")
            ok = dd < TOL and dp == 0.0 and dc == 0.0
            row["ka"]["ok"] = bool(ok)
            if not ok:
                failures.append(f"KA {name} vs 257 {r257}: endpoints {dd:.3e}, pertarget pds {dp:.3e}, cos {dc:.3e}")
            print(f"  KA {name} vs 257 {r257}: endpoints {dd:.1e}, per-target pds {dp:.1e} cos {dc:.1e} "
                  f"-> {'PASS' if ok else 'FAIL'}", flush=True)
        summary[name] = row
        for k, v in g.items():
            arrays[f"{name}|{k}"] = np.asarray(v, dtype=np.float64)
        for k in ("pds", "own", "comp_max", "comp_med", "comp_mean", "n_beat", "pred_nn", "pred_mean"):
            arrays[f"{name}|effect_{k}"] = np.asarray(ge[k], dtype=np.float64)
        arrays[f"{name}|beats"] = B.astype(np.int8)
        arrays[f"{name}|ties"] = Tt.astype(np.int8)
        arrays[f"{name}|effect_beats"] = Be.astype(np.int8)
        arrays[f"{name}|S_A"] = S_ctx["A"].astype(np.float32)   # frame A cosine matrix, for the record
        print(f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
              f"effect-level pds {float(ge['pds'].mean()):.6f}  ({row['seconds']}s)", flush=True)

    # ------------------------------------------------------------------ structural covariates
    spec = plolo.FOLDS[held]
    banks = {nm: E062.Bank.load(nm, plolo.DONOR_PATHS[nm]) for nm in spec["sources"]}
    src_carried = {nm: np.asarray([t in b.pert_pos for t in targets]) for nm, b in banks.items()}
    pm = E064.load_perturbmatch(plolo.PM_PATH)
    pm_carried = np.asarray([("k562_gwps" in banks) and (t in pm.pert_pos) and (t in banks["k562_gwps"].pert_pos)
                             for t in targets])
    n_whitened = sum(v.astype(int) for v in src_carried.values())
    cd4 = st["cd4"]
    cd4_cnt = np.asarray(cd4.target_condition_count, dtype=np.int64)
    cd4_added = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0)
    cd4_norm = np.linalg.norm(cd4_added, axis=1)
    d_added = np.where(np.asarray(d.gene_votes) > 0, np.asarray(d.mean_effect), 0.0)
    d_norm = np.linalg.norm(d_added, axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cd4_energy_removed = np.where(cd4_norm > 0, 1.0 - d_norm**2 / np.maximum(cd4_norm**2, 1e-300), np.nan)
    u = st["union"]
    u_cnt = np.asarray(u.target_condition_count, dtype=np.int64)
    u_added = np.where(np.asarray(u.gene_votes) > 0, np.asarray(u.mean_effect), 0.0)
    u_norm = np.linalg.norm(u_added, axis=1)
    src_u = {nm: (np.asarray(L.source_residual(nm, tuple(targets), "gene_norm").target_condition_count) > 0)
             for nm in ("gse295214_h1", "h1_2025")}
    essential = np.asarray([t in ESSENTIAL for t in targets])

    sc = fold["scored"]
    W = np.asarray(st["whitened"], dtype=np.float64)
    Ws = W[:, sc]
    R = np.asarray(L.merged_many(arms["cd4svd5"]).mean_effect, dtype=np.float64)   # fr*d + frac_x*u
    Rs = R[:, sc]
    R0 = np.asarray(L.merged_many(arms["h1union"]).mean_effect, dtype=np.float64)   # frac_cd4*cd4 + frac_x*u
    R0s = R0[:, sc]
    Tu = unit(truth)

    def cos_rows(A, Bm):
        na = np.linalg.norm(A, axis=1)
        nb = np.linalg.norm(Bm, axis=1)
        return np.where((na > 0) & (nb > 0), (A * Bm).sum(1) / np.maximum(na * nb, 1e-300), np.nan)

    cov = dict(
        n_whitened=n_whitened.astype(np.float64),
        pm_carried=pm_carried.astype(np.float64),
        cd4_count=cd4_cnt.astype(np.float64),
        cd4_complete=(cd4_cnt == 3).astype(np.float64),
        cd4_norm=cd4_norm, cd4_deflated_norm=d_norm, cd4_energy_removed=cd4_energy_removed,
        union_count=u_cnt.astype(np.float64), union_norm=u_norm,
        union_gse295214=src_u["gse295214_h1"].astype(np.float64),
        union_h1_2025=src_u["h1_2025"].astype(np.float64),
        essential=essential.astype(np.float64),
        base_norm=np.linalg.norm(Ws, axis=1),
        base_cos_truth=cos_rows(Ws, truth),
        resid_norm=np.linalg.norm(Rs, axis=1),
        resid_cos_truth=cos_rows(Rs, truth),
        resid_cos_base=cos_rows(Rs, Ws),
        resid_share=np.linalg.norm(Rs, axis=1) / np.maximum(np.linalg.norm(Ws + Rs, axis=1), 1e-300),
        resid0_norm=np.linalg.norm(R0s, axis=1),
        resid0_cos_truth=cos_rows(R0s, truth),
        cd4_cos_truth=cos_rows(cd4_added[:, sc], truth),
        cd4_deflated_cos_truth=cos_rows(d_added[:, sc], truth),
        union_cos_truth=cos_rows(u_added[:, sc], truth),
        cd4_cos_base=cos_rows(cd4_added[:, sc], Ws),
        cd4_deflated_cos_base=cos_rows(d_added[:, sc], Ws),
        # the truth's nearest-neighbour cosine to the OTHER truths and the prediction's are in the arm
        # geometry above (identical truth side for every arm)
    )
    for nm, v in src_carried.items():
        cov[f"src_{nm}"] = v.astype(np.float64)
    for k, v in cov.items():
        arrays[f"cov|{k}"] = np.asarray(v, dtype=np.float64)
    arrays["targets"] = np.asarray(targets)

    meta = dict(fold=held, n_targets=n, sources=list(spec["sources"]), op=str(fold["op_path"]),
                frac_cd4_stack=st["frac_cd4"], frac_cd4_resolved=float(fr), frac_x=st["frac_x"],
                deflation=aud_defl, n_essential_on_panel=int(essential.sum()),
                n_whitened_hist={str(k): int((n_whitened == k).sum()) for k in range(4)},
                cd4_count_hist={str(k): int((cd4_cnt == k).sum()) for k in range(4)},
                union_count_hist={str(k): int((u_cnt == k).sum()) for k in range(3)},
                pm_carried=int(pm_carried.sum()),
                fold_cache=fold["audit"].get("cache"), failures=failures,
                seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    np.savez_compressed(WORK / f"diag262-{held}.npz", **arrays)
    (WORK / f"diag262-{held}.json").write_text(json.dumps(dict(meta=meta, arms=summary), indent=1, default=float))
    print(f"[{held}] done: {len(failures)} KA failures, {meta['seconds']}s, peak RSS {meta['peak_rss_gib']:.2f} GiB "
          f"-> {WORK / f'diag262-{held}.npz'}", flush=True)
    for f in failures:
        print(f"  *** {f}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
