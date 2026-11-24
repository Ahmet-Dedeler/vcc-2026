"""277 -- the class-conditional RENDERINGS on one PLOLO fold.  MEASUREMENT ONLY.

    PYTHONDONTWRITEBYTECODE=1 uv run python render277.py --fold <k562_gwps|orion_hct116|orion_hek293t|kolf2>

Constructions frozen in `work/DEFINITION-277.md` (sha `05a89391...`) before this file scored anything.
One fold per process (peak RSS ~5 GiB, the operator eigenbasis).

KNOWN-ANSWER ROWS (exit 1 on any failure):
  KA-1  V0 (`cd4svd5`) endpoints == 257 `sweep257-<fold>.json` row `k5_s1` at 1e-9 AND its per-target
        `pds` / cosine == 262 `diag262-<fold>.npz` `cd4svd5|pds` / `|cos` at max|d| == 0.0.
  KA-1b 262's `h1union` (the un-deflated stack) likewise, and 262's `cov|base_cos_truth` reproduced at 0.0.
  KA-3  V4 with CLASS = empty reproduces V0's composed matrix at max|d| == 0.0.

Writes `work/render277-<fold>.npz` (per-target arrays) and `work/render277-<fold>.json` (endpoints, the
class, the audits).  Nothing outside `work/` is written.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib277 as Z  # noqa: E402

L, X, plolo, O, E062 = Z.L, Z.X, Z.plolo, Z.O, Z.E062
TOL = 1e-9


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    fail: list[str] = []
    log: list[str] = []

    def P(s: str = "") -> None:
        print(s, flush=True)
        log.append(s)

    P(f"277 render277.py --fold {held} -- {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}")

    # ------------------------------------------------------------------------------------ stage
    st = L.stage(held)
    fold, targets = st["fold"], [str(t) for t in st["targets"]]
    truth = np.asarray(st["truth"], dtype=np.float64)
    W = np.asarray(st["whitened"], dtype=np.float64)
    sc = fold["scored"]
    n = len(targets)
    P(f"  staged {n} targets, {fold['wcols'].size} operator columns; frac_cd4 {st['frac_cd4']:.6f} "
      f"frac_x {st['frac_x']:.6f}  ({time.time()-t0:.0f}s, peak {Z.peak_gib():.2f} GiB)")

    # ------------------------------------------------------------------------------------ the class
    cl = Z.c_lolo(held, targets)
    n_elig = int(np.isfinite(cl).sum())
    masks: dict[str, np.ndarray] = {}
    masks["CLASS30"] = Z.bottom_class(cl, targets, 0.30)
    masks["CLASS20"] = Z.bottom_class(cl, targets, 0.20)
    masks["CLASS40"] = Z.bottom_class(cl, targets, 0.40)
    masks["TOPCONS30"] = Z.bottom_class(cl, targets, 0.30, top=True)
    n_class = int(masks["CLASS30"].sum())
    # 262's donor-donor agreement on this fold (276's `<fold>|dd`, aligned by target symbol)
    zp = np.load(Z.E276 / "work/predict276.npz", allow_pickle=True)
    dd_map = {str(t): float(v) for t, v in zip(zp[f"{held}|targets"], zp[f"{held}|dd"])}
    dd = np.asarray([dd_map.get(t, np.nan) for t in targets], dtype=np.float64)
    masks["AGREE30"] = Z.bottom_class(dd, targets, 0.30)
    # POST-HOC: the same control SIZE-MATCHED to CLASS30 (fewer targets have a donor-donor agreement, so
    # the frozen "bottom 30 % of eligible" rule gives a smaller set; both are reported)
    masks["AGREEN"] = Z.bottom_class(dd, targets, n_class / max(int(np.isfinite(dd).sum()), 1))
    for s, seed in enumerate((20260906, 20260907, 20260908)):
        masks[f"RAND30_s{s+1}"] = Z.random_class(cl, n_class, seed)
    P(f"  CLASS: {n_elig}/{n} eligible (finite C_lolo_{held}); CLASS30 n = {n_class} "
      f"(C_lolo cut {np.nanmax(cl[masks['CLASS30']]):+.5f}); CLASS20 {int(masks['CLASS20'].sum())}, "
      f"CLASS40 {int(masks['CLASS40'].sum())}, TOPCONS30 {int(masks['TOPCONS30'].sum())}, "
      f"AGREE30 {int(masks['AGREE30'].sum())} (overlap with CLASS30 "
      f"{int((masks['AGREE30'] & masks['CLASS30']).sum())})")

    # ------------------------------------------------------------------------------------ the objects
    cd4 = st["cd4"]
    d, aud5 = L.deflate(cd4, "svd", 5)
    carried = np.asarray(cd4.target_condition_count) > 0
    fr0, pu0 = L.solve_fraction(W, d, Z.CD4_F)
    P(f"  CD4: {int(carried.sum())}/{n} complete rows; deflation k = 5 removes "
      f"{aud5['energy_removed_frac']:.6f} of the carried energy (top-5 share of centred "
      f"{aud5['top_k_energy_frac_of_centred']:.6f}); V0 fraction {fr0:.6f} (per unit {pu0:.6f})")
    P(f"  CLASS30 n cd4-carried {int((masks['CLASS30'] & carried).sum())}/{n_class}; "
      f"AGREE30 {int((masks['AGREE30'] & carried).sum())}; TOPCONS30 "
      f"{int((masks['TOPCONS30'] & carried).sum())}")

    idn = X.Identity(st)
    base_cos = Z.row_cos(W[:, sc], truth)
    grp = Z.groups_262(base_cos)

    arms: dict[str, dict] = {}
    audits: dict[str, dict] = {}
    arrays: dict[str, np.ndarray] = {}
    finals: dict[str, np.ndarray] = {}

    def add(name: str, parts, *, keep_final: bool = False, **meta) -> dict:
        t1 = time.time()
        final = L.compose(st, parts)
        e = L.score(st, final)
        idm = idn.measure(final)
        pds_t = np.asarray(e["_pds_per_target_a23"], dtype=np.float64)
        cos_t = np.asarray(e["_cos_per_target_a23"], dtype=np.float64)
        arrays[f"{name}|pds"] = pds_t
        arrays[f"{name}|cos"] = cos_t
        for k in ("U", "M", "D"):
            arrays[f"{name}|{k}"] = np.asarray(idm[k], dtype=np.float64)
        row = dict(arm=name, **{k: float(v) for k, v in e.items() if not k.startswith("_")},
                   **{k: float(np.mean(v)) for k, v in idm.items()},
                   seconds=round(time.time() - t1, 1))
        arms[name] = row
        audits[name] = dict(meta)
        if keep_final:
            finals[name] = final
        P(f"  [{name:16s}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} | "
          f"nmae {idm['nmae_wald89']:.6f} delivered {idm['delivered_norm']:.5f}  ({row['seconds']}s)")
        return row

    # -------------------------------------------------------------- V0 + the KA rows
    add("V0_cd4svd5", [(d, fr0), (st["union"], st["frac_x"])], keep_final=True,
        fraction=fr0, per_unit=pu0, **{k: (v if not isinstance(v, np.ndarray) else None)
                                       for k, v in aud5.items()})
    add("h1union", L.stack_parts(st), fraction=st["frac_cd4"])
    add("BASE", [])

    ref257 = {r["arm"]: r for r in json.loads((Z.E257 / f"work/sweep257-{held}.json").read_text())["rows"]}
    pt262 = np.load(Z.E262 / f"work/diag262-{held}.npz", allow_pickle=True)
    ka: list[dict] = []
    for name, r257, k262 in (("V0_cd4svd5", "k5_s1", "cd4svd5"), ("h1union", "STACK", "h1union")):
        dd_end = max(abs(float(ref257[r257][k]) - float(arms[name][k])) for k in Z.KEYS)
        dp = float(np.abs(np.asarray(pt262[f"{k262}|pds"]) - arrays[f"{name}|pds"]).max())
        dc = float(np.abs(np.asarray(pt262[f"{k262}|own"]) - arrays[f"{name}|cos"]).max())
        ok = dd_end < TOL and dp == 0.0 and dc == 0.0
        ka.append(dict(row=f"KA-1 {name}", endpoints_vs_257=dd_end, pertarget_pds_vs_262=dp,
                       pertarget_cos_vs_262=dc, ok=bool(ok)))
        P(f"  KA-1 {name}: endpoints vs 257 {r257} {dd_end:.1e}; per-target pds vs 262 {dp:.1e}, "
          f"cos {dc:.1e} -> {'PASS' if ok else 'FAIL'}")
        if not ok:
            fail.append(f"KA-1 {name}")
    dbc = float(np.abs(np.asarray(pt262["cov|base_cos_truth"]) - base_cos).max())
    ka.append(dict(row="KA-1b base_cos_truth vs 262", value=dbc, ok=bool(dbc == 0.0)))
    P(f"  KA-1b base.truth cosine vs 262 max|d| {dbc:.1e} -> {'PASS' if dbc == 0.0 else 'FAIL'}")
    if dbc != 0.0:
        fail.append("KA-1b base_cos_truth")

    # -------------------------------------------------------------- V1 / V2 / V3 (+ the fixed-fraction diagnostics)
    def cd4_arm(name: str, m, *, resolve: bool = True):
        if resolve:
            fr, pu = L.solve_fraction(W, m, Z.CD4_F)
        else:
            fr, pu = fr0, pu0
        nrm = np.linalg.norm(np.where(np.asarray(m.gene_votes) > 0,
                                      np.asarray(m.mean_effect), 0.0), axis=1)
        return add(name, [(m, fr), (st["union"], st["frac_x"])], fraction=fr, per_unit=pu,
                   resolved=resolve, median_row_norm_all=float(np.median(nrm)),
                   carrying=int((nrm > 0).sum()))

    C30 = masks["CLASS30"]
    m1 = Z.mix_rows(cd4, d, C30, "V1")          # un-deflated ON the class, deflated elsewhere
    m2 = Z.mix_rows(d, cd4, C30, "V2")          # the mirror
    m3 = Z.zero_rows(d, C30, "V3")              # no CD4 residual on the class
    cd4_arm("V1_undefl", m1)
    cd4_arm("V1f_undefl", m1, resolve=False)
    cd4_arm("V2_mirror", m2)
    cd4_arm("V2f_mirror", m2, resolve=False)
    cd4_arm("V3_noresid", m3)
    cd4_arm("V3f_noresid", m3, resolve=False)

    # -------------------------------------------------------------- V5 boundary + the controls
    for tag, key in (("V5_20", "CLASS20"), ("V5_40", "CLASS40"),
                     ("C_top30", "TOPCONS30"), ("C_agree30", "AGREE30"), ("C_agreeN", "AGREEN"),
                     ("C_rand1", "RAND30_s1"), ("C_rand2", "RAND30_s2"), ("C_rand3", "RAND30_s3")):
        cd4_arm(tag, Z.mix_rows(cd4, d, masks[key], tag))

    # -------------------------------------------------------------- V4: the soft noise floor on the class only
    # kolf2 runs on the PRODUCTION operator, whose noise spectrum 272 saved under its own name
    spz_path = Z.E272 / (f"work/spectrum272-{held}.npz" if plolo.FOLDS[held]["op_kind"] != "production"
                         else "work/spectrum272-production.npz")
    if spz_path.exists():
        op = plolo.Operator(L.op_path_for(held))
        spz = np.load(spz_path)
        assert np.array_equal(spz["evals"], op.evals), "272's spectrum is not this fold's operator"
        nvec = np.asarray(spz["n"], dtype=np.float64)
        w_inc = (op.evals + Z.LAM) ** (-Z.GAMMA / 2.0)
        w_flr = (np.maximum(op.evals, 1.0 * nvec) + Z.LAM) ** (-Z.GAMMA / 2.0)
        pre = np.asarray(st["pre"], dtype=np.float64) * L.ALPHA
        # KA-3: an EMPTY class is `plolo.Operator.apply` bit for bit
        W_empty = Z.mixed_whiten(pre, fold["wcols"], op.evecs, w_flr, w_inc,
                                 np.zeros(n, dtype=bool))
        d3 = float(np.abs(W_empty - W).max())
        ka.append(dict(row="KA-3 mixed_whiten(class=empty) == Operator.apply", value=d3, ok=bool(d3 == 0.0)))
        P(f"  KA-3 mixed_whiten with CLASS = empty vs the staged whitened base: max|d| {d3:.1e} "
          f"-> {'PASS' if d3 == 0.0 else 'FAIL'}")
        if d3 != 0.0:
            fail.append("KA-3 mixed_whiten")
        W4 = Z.mixed_whiten(pre, fold["wcols"], op.evecs, w_flr, w_inc, C30)
        frac_damped = float(np.mean((w_flr / w_inc) ** 2))
        del op, W_empty
        st4 = dict(st)
        st4["whitened"] = W4
        fr4, pu4 = L.solve_fraction(W4, d, Z.CD4_F)
        frx4, _ = L.solve_fraction_carrying(W4, st["union"], L.H1_G)
        t1 = time.time()
        final4 = L.compose(st4, [(d, fr4), (st["union"], frx4)])
        e4 = L.score(st4, final4)
        # the identity is measured on the same fold objects (the truth side is unchanged)
        id4 = idn.measure(final4)
        arrays["V4_floor|pds"] = np.asarray(e4["_pds_per_target_a23"], dtype=np.float64)
        arrays["V4_floor|cos"] = np.asarray(e4["_cos_per_target_a23"], dtype=np.float64)
        for k in ("U", "M", "D"):
            arrays[f"V4_floor|{k}"] = np.asarray(id4[k], dtype=np.float64)
        arms["V4_floor"] = dict(arm="V4_floor", **{k: float(v) for k, v in e4.items() if not k.startswith("_")},
                                **{k: float(np.mean(v)) for k, v in id4.items()},
                                seconds=round(time.time() - t1, 1))
        bc4 = Z.row_cos(W4[:, sc], W[:, sc])
        audits["V4_floor"] = dict(fraction=fr4, per_unit=pu4, frac_x=frx4,
                                  mean_weight_ratio_sq=frac_damped,
                                  base_cos_to_incumbent_class=float(np.median(bc4[C30])),
                                  base_cos_to_incumbent_rest=float(np.median(bc4[~C30])))
        arrays["V4_floor|base_cos_to_incumbent"] = bc4
        P(f"  [{'V4_floor':16s}] cos_agg {e4['cos_agg']:.6f} macro {e4['a23_cos']:.6f} pds {e4['pds']:.6f} | "
          f"nmae {id4['nmae_wald89']:.6f} delivered {id4['delivered_norm']:.5f}; class base cosine to the "
          f"incumbent base {np.median(bc4[C30]):.5f} (rest {np.median(bc4[~C30]):.5f})")
        del W4, final4
    else:
        P(f"  V4 SKIPPED: no 272 noise spectrum for {held} ({spz_path} absent)")
        audits["V4_floor"] = dict(skipped=f"no spectrum272-{held}.npz")

    # ------------------------------------------------------------------------------------ save
    for k, m in masks.items():
        arrays[f"mask|{k}"] = m
    arrays["targets"] = np.asarray(targets)
    arrays["c_lolo"] = cl
    arrays["dd_agree"] = dd
    arrays["base_cos_truth"] = base_cos
    arrays["cd4_carried"] = carried
    for gname, gm in grp.items():
        arrays[f"grp|{gname}"] = gm

    meta = dict(fold=held, n_targets=n, n_eligible=n_elig, n_class=n_class,
                class_cut=float(np.nanmax(cl[C30])), frac_cd4_stack=st["frac_cd4"],
                frac_x=st["frac_x"], v0_fraction=fr0, deflation=aud5,
                op=str(L.op_path_for(held)), fold_cache=fold["audit"].get("cache"),
                seconds=round(time.time() - t0, 1), peak_rss_gib=Z.peak_gib(),
                ka=ka, failures=fail)
    np.savez_compressed(Z.WORK / f"render277-{held}.npz", **arrays)
    (Z.WORK / f"render277-{held}.json").write_text(
        json.dumps(dict(meta=meta, arms=arms, audits=audits, log=log), indent=1, default=float))
    P(f"[{held}] done: {len(fail)} KA failures, {meta['seconds']}s, peak RSS {meta['peak_rss_gib']:.2f} GiB")
    for f in fail:
        P(f"  *** {f}")
    if fail:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
