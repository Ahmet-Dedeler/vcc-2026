#!/usr/bin/env python3
"""400 harness: cap15 known answer + the shipped 399 construction on H1 (per-target six-bank neighbourhood tilt, K=25, shrink 400,
q80 protection) and its extension to UNCALLED genes with confident neighbourhood signs (397c), six members end to end.

Call sequence per target is 385's exact chain (profile / spread / rank-axis /
block + CallSet(400, 2200, 1.9) / block_stats / score); only the block function
and its kwargs vary by arm. Reads 323's reference and cap15 effects read-only;
writes only into experiments/391-consensus-block-signs/work/.

Arms: cap15 (production rc_hybrid_block; S0 gate), k1 (rc_hybrid_block_k K=1;
must match cap15), k3plain (S1), k3s125/k3s148/k3s200 (S2 dose), k3shuf
(shuffled-class control at 1.48x).

Usage: uv run python run389.py cap15|k1|k3plain|k3s125|k3s148|k3s200|k3shuf
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault(
    "NUMBA_CACHE_DIR",
    str(Path(__file__).resolve().parent / "work" / "numba-cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("POLARS_MAX_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np
import polars as pl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
BASE = ROOT / "experiments/323-current-family-validation"
for _p in (str(HERE), str(BASE), str(ROOT / "src")):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import run323 as R  # noqa: E402 -- read-only import, never edited
import evaluate323 as E  # noqa: E402 -- read-only import, never edited
from vcc2026 import rc_emit  # noqa: E402
from vcc2026.donor_transfer import emit_counts as EMIT_COUNTS  # noqa: E402
import k_emit400 as K  # noqa: E402

CLASSES3 = ((1 / 3, 16078.75), (1 / 3, 20190.5), (1 / 3, 25330.5))
CLASSES3W = None   # filled at runtime from the H1 control libraries (q10 / q50 / q90) for the class-spread arm
CALLSET = (400, 2200, 1.9)
SHUF_SEED = 20260913
PRIOR_BANKS = ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")   # CLEAN4 (371); no H1 data
KNN_K, KNN_M = 25, 400.0
BANK_DIRS = {"k562_gwps": ROOT / "data/donors/k562_gwps", "k562_essential": ROOT / "data/donors/k562_essential",
             "nadig_jurkat": ROOT / "data/donors/nadig_jurkat", "nadig_hepg2": ROOT / "experiments/166-donor-acquisition/donors/nadig_hepg2",
             "orion_hct116": ROOT / "data/donors/orion_hct116", "orion_hek293t": ROOT / "data/donors/orion_hek293t"}


def prior_signs_full(genes_full):
    """CLEAN4 consensus sign on the full reference gene axis: -1 down if mean P_down > 0.5, +1 up; uncovered -> -1."""
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str)); PB = cp["P_banks"].astype(float)
    P = PB[[names.index(b) for b in PRIOR_BANKS]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(str(g), -1) for g in genes_full])
    out = np.full(len(genes_full), -1.0)
    covered = col >= 0
    out[covered] = np.where(P[col[covered]] > 0.5, -1.0, 1.0)
    return out, float(covered.mean())


PRIOR_FULL = None
PRIOR_BY_TARGET, CONF_BY_TARGET = {}, {}


def knn_priors_full(genes_full, targets, fallback):
    """399's construction on the given gene axis: per-target sign (+ confidence) from the 25 nearest GWPS perturbations' votes in
    the six-bank family, shrunk toward the family marginal (m=400); default = family marginal sign, `fallback` where no vote."""
    off = [str(g) for g in genes_full]; nG = len(off)
    B = {}
    for b, d in BANK_DIRS.items():
        D = np.asarray(np.load(d / "delta.npy", mmap_mode="r"), dtype=np.float32)
        P = np.array([l.strip() for l in open(d / "perts.txt")]); G = np.array([l.strip() for l in open(d / "genes.txt")])
        gi = {g: i for i, g in enumerate(G)}; B[b] = (D, {p: i for i, p in enumerate(P)}, np.array([gi.get(g, -1) for g in off]), P)
    GW, gwp, _, GWP = B["k562_gwps"]
    GWc = GW - GW.mean(1, keepdims=True); GWn = GWc / np.maximum(np.linalg.norm(GWc, axis=1, keepdims=True), 1e-12)

    def votes(rows):
        dn = np.zeros(nG); nn = np.zeros(nG); names = [GWP[j] for j in rows]
        for b, (D, pi, gm, _) in B.items():
            ridx = [pi[p] for p in names if p in pi]
            if not ridx:
                continue
            cv = gm >= 0; sub = D[ridx][:, np.where(cv, gm, 0)]; nz = (sub != 0) & cv[None, :]
            dn += ((sub < 0) & nz).sum(0); nn += nz.sum(0)
        return dn, nn
    dn_m, nn_m = votes(np.arange(len(GWP))); P_marg = np.where(nn_m > 0, dn_m / np.maximum(nn_m, 1), np.nan)
    default = np.where(np.isfinite(P_marg), np.where(P_marg > 0.5, -1.0, 1.0), fallback)
    by_t, conf_t, n_rows = {}, {}, 0
    for t in map(str, targets):
        if t not in gwp:
            continue
        sim = GWn @ GWn[gwp[t]]; sim[gwp[t]] = -np.inf; nb = np.argsort(-sim)[:KNN_K]
        dn, nn = votes(nb); P = np.where(nn > 0, dn / np.maximum(nn, 1), np.nan)
        Ps = np.where(np.isfinite(P) & np.isfinite(P_marg), (nn * np.nan_to_num(P) + KNN_M * np.nan_to_num(P_marg)) / (nn + KNN_M),
                      np.where(np.isfinite(P), P, P_marg))
        by_t[t] = np.where(np.isfinite(Ps), np.where(Ps > 0.5, -1.0, 1.0), default)
        conf_t[t] = np.where(np.isfinite(Ps), np.abs(Ps - 0.5), 0.0); n_rows += 1
    return default, by_t, conf_t, n_rows


def builder389():
    with R.model_env("uniform"):
        b = R.load_module("builder389", ROOT / "scripts/make_submission_pm_cd4.py")
    assert b.H1_UNION_G is None and b.MASS_FIX is None
    return b


def native_cap15(builder_module):
    with np.load(BASE / "work/effects.npz", allow_pickle=False) as z:
        comp = json.loads((BASE / "work/compose.json").read_text())
        if builder_module.array_sha256(z["cap15"]) != comp["arms"]["cap15"]["final_sha256"]:
            raise ValueError("cap15 final effect hash drift vs compose.json")
        return z["cap15"].copy(), z["genes"].copy(), z["targets"].copy()


def kernel_check(emit_c):
    """One-time hist_z_k == wilcoxon_z proof; preserved, not re-run per arm."""
    out = WORK / "kernel389.json"
    if out.exists():
        return json.loads(out.read_text())
    w = K.check_kernel_equivalence(emit_c.cs, emit_c.tie)
    rec = {"max_abs_dz": w, "pass": bool(w < 1e-9)}
    R.write_json(out, rec)
    print(f"kernel equivalence: max|dz|={w:.3e}", flush=True)
    if not rec["pass"]:
        raise AssertionError("kernel equivalence failed")
    return rec


def block_for_arm(arm, profile, libs, emit_rng, spread, key, i, emit_c):
    cs = rc_emit.CallSet(*CALLSET)
    if arm == "cap15":
        return rc_emit.rc_hybrid_block(
            profile, libs, .1, emit_rng, spread=spread,
            overdispersion=4., control=emit_c, key=key,
            jitter_seed=1000 + i, call_set_first=True,
            call_set=cs, emit_counts=EMIT_COUNTS)
    kw = dict(spread=spread, overdispersion=4., control=emit_c, key=key,
              jitter_seed=1000 + i, call_set_first=True, call_set=cs,
              emit_counts=EMIT_COUNTS)
    if arm == "k1":
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng, classes=None,
                                   lfc_target=None, **kw)
    if arm == "k3plain":
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng,
                                   classes=list(CLASSES3), lfc_target=None, **kw)
    if arm == "k3c5bq80":
        gated = np.asarray(emit_c.keep, dtype=bool)
        spec = {"prior": np.ascontiguousarray(PRIOR_FULL[gated]), "shrink_floor": 0.002, "protect_q": 0.8}
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng,
                                   classes=list(CLASSES3), lfc_target=spec, **kw)
    if arm in ("k3knn", "k3knn_u20", "k3knn_u10", "k3knn_u05", "k3knn_q50", "k3knn_q65", "k3knn_q90", "k3knn_q100", "k3knn_wide"):
        gated = np.asarray(emit_c.keep, dtype=bool); t = str(TARGET_NAMES[i])
        psign = PRIOR_BY_TARGET.get(t, KNN_DEFAULT)
        pq = {"k3knn_q50": 0.5, "k3knn_q65": 0.65, "k3knn_q90": 0.9}.get(arm, 0.8)
        spec = {"prior": np.ascontiguousarray(psign[gated]), "shrink_floor": 0.002}
        if arm != "k3knn_q100":
            spec["protect_q"] = pq
        if arm in ("k3knn_u20", "k3knn_u10", "k3knn_u05"):
            tau = {"k3knn_u20": 0.2, "k3knn_u10": 0.1, "k3knn_u05": 0.05}[arm]
            conf = CONF_BY_TARGET.get(t, np.zeros(len(psign)))
            spec["extra_mask"] = np.ascontiguousarray((conf >= tau)[gated])
        classes = list(CLASSES3W) if arm == "k3knn_wide" else list(CLASSES3)
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng,
                                   classes=classes, lfc_target=spec, **kw)
    if arm in ("k3c5b", "k3c5bshuf", "k3down"):
        gated = np.asarray(emit_c.keep, dtype=bool)
        psign = PRIOR_FULL[gated] if arm != "k3down" else -np.ones(int(gated.sum()))
        spec = {"prior": np.ascontiguousarray(psign), "shrink_floor": 0.002}
        if arm == "k3c5bshuf":
            spec["permute"] = SHUF_SEED
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng,
                                   classes=list(CLASSES3), lfc_target=spec, **kw)
    if arm in ("k3s125", "k3s148", "k3s200", "k3shuf"):
        s = {"k3s125": 1.25, "k3s148": 1.48, "k3s200": 2.0, "k3shuf": 1.48}[arm]
        mask = K.prescribed_mask(key, emit_c, cs)
        spec = {"scale": s, "mask": np.ascontiguousarray(mask)}
        if arm == "k3shuf":
            spec["permute"] = SHUF_SEED
        return K.rc_hybrid_block_k(profile, libs, .1, emit_rng,
                                   classes=list(CLASSES3), lfc_target=spec, **kw)
    raise ValueError(f"unknown arm {arm}")


def emit_arm(name, effects_native, genes) -> dict:
    if (WORK / f"metrics-{name}.json").exists():
        raise FileExistsError(f"refusing to overwrite metrics-{name}.json")
    b = builder389()
    r, c, emit_c = E.load_reference()
    global PRIOR_FULL, PRIOR_BY_TARGET, CONF_BY_TARGET, KNN_DEFAULT, TARGET_NAMES
    PRIOR_FULL, cov = prior_signs_full(r["genes"])
    print(f"CLEAN4 prior coverage of the full axis {cov:.4f}", flush=True)
    TARGET_NAMES = [str(t) for t in r["targets"]]
    global CLASSES3W
    _q = np.quantile(np.asarray(r["libraries"], dtype=float), [0.1, 0.5, 0.9])
    CLASSES3W = ((1 / 3, float(_q[0])), (1 / 3, float(_q[1])), (1 / 3, float(_q[2])))
    if name.startswith("k3knn"):
        KNN_DEFAULT, PRIOR_BY_TARGET, CONF_BY_TARGET, n_rows = knn_priors_full(r["genes"], r["targets"], PRIOR_FULL)
        print(f"kNN priors: {n_rows}/{len(TARGET_NAMES)} targets with a row; default agrees with CLEAN4 on "
              f"{float((KNN_DEFAULT == PRIOR_FULL).mean()):.3f} of the axis", flush=True)
    if name != "cap15":
        kernel_check(emit_c)
    effects = R.map_axis(effects_native, genes, r["genes"])
    pick_rng, emit_rng = b.make_rngs(0)
    sums, jks, ps, qs, lfcs, digests, diags = [], [], [], [], [], [], []
    counts_sha = hashlib.sha256()
    t0 = time.monotonic()
    for i, target in enumerate(r["targets"]):
        profile = b.profile_from_effect(effects[i], r["b_ctrl"], min_cpm=10.)
        picks = pick_rng.choice(len(r["libraries"]), 400, replace=False)
        spread = b.emission_spread_vector(
            effects[i], r["b_ctrl"], min_cpm=10., tail_spread=1., keep_top=8000,
            keep_top_frac=None, keep_spread=.25, rank_score=None, spread_map="block",
            spread_lo=0., spread_hi=.60)
        gate = b.emission_rank_axis(effects[i], r["b_ctrl"], min_cpm=10.)[0]
        try:
            block, dg = block_for_arm(
                name, profile, r["libraries"][picks], emit_rng, spread,
                np.where(gate, effects[i], 0.), i, emit_c)
        except BaseException as exc:
            R.write_json(WORK / f"fail-{name}.json",
                         {"arm": name, "target": str(target), "index": i,
                          "error": f"{type(exc).__name__}: {exc}"})
            raise
        if K.LAST_PRIOR_REC:
            dg = {**dg, **K.LAST_PRIOR_REC.pop()}
        dg = {k: (float(v) if isinstance(v, (np.floating, np.integer)) else v)
              for k, v in dg.items()}
        dg["target"] = str(target)
        diags.append(dg)
        s, j, p, q, l = E.block_stats(block, c)
        sums.append(s)
        jks.append(j)
        ps.append(p)
        qs.append(q)
        lfcs.append(l)
        digests.append(hashlib.sha256(np.ascontiguousarray(s).tobytes()).hexdigest())
        counts_sha.update(np.asarray(block, dtype="<u4", order="C").tobytes())
        del block
        if (i + 1) % 20 == 0:
            print(f"{name} {i + 1}/{len(r['targets'])} ({time.monotonic() - t0:.1f}s)",
                  flush=True)
    panel_sha = hashlib.sha256("".join(digests).encode()).hexdigest()
    raw, rows, legs = E.score(r, c, sums, jks, ps, qs, lfcs)
    np.savez(WORK / f"stats-{name}.npz", targets=r["targets"], sums=np.asarray(sums),
             jk=np.asarray(jks), p=np.asarray(ps), q=np.asarray(qs), lfc=np.asarray(lfcs))
    pl.DataFrame(rows).write_csv(WORK / f"per-target-{name}.csv")
    R.write_json(WORK / f"diag-{name}.json", diags)
    out = {"arm": name, "raw": raw, "panel_total_digest_sha256": panel_sha,
           "emitted_count_digest_sha256": counts_sha.hexdigest(),
           "mse_legs": legs, "seconds": time.monotonic() - t0,
           "statistics_sha256": R.sha(WORK / f"stats-{name}.npz")}
    R.write_json(WORK / f"metrics-{name}.json", out)
    print(name, json.dumps({"raw": raw}), flush=True)
    return out


def known_answer(receipt: dict) -> dict:
    expected = json.loads((BASE / "work/metrics-cap15.json").read_text())
    diffs = {m: abs(receipt["raw"][m] - expected["raw"][m]) for m in E.MEMBERS}
    checks = {
        "emitted_count_digest_equal": receipt["emitted_count_digest_sha256"]
        == expected["emitted_count_digest_sha256"],
        "panel_total_digest_equal": receipt["panel_total_digest_sha256"]
        == expected["panel_total_digest_sha256"],
        "raw_metrics_within_tolerance": {m: diffs[m] <= E.TOLERANCE[m] for m in E.MEMBERS},
        "mse_legs_exact": receipt["mse_legs"] == expected["mse_legs"],
    }
    passed = bool(checks["emitted_count_digest_equal"] and checks["panel_total_digest_equal"]
                  and all(checks["raw_metrics_within_tolerance"].values())
                  and checks["mse_legs_exact"])
    result = {"pass": passed, "checks": checks, "absolute_differences": diffs,
              "tolerance": E.TOLERANCE,
              "expected_source": "experiments/323-current-family-validation/work/metrics-cap15.json"}
    R.write_json(WORK / f"known-answer-{receipt['arm']}.json", result)
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arm", choices=["cap15", "k1", "k3plain", "k3s125", "k3s148",
                                   "k3s200", "k3shuf", "k3c5b", "k3c5bshuf", "k3down", "k3c5bq80",
                                   "k3knn", "k3knn_u20", "k3knn_u10", "k3knn_u05",
                                   "k3knn_q50", "k3knn_q65", "k3knn_q90", "k3knn_q100", "k3knn_wide"])
    args = ap.parse_args()
    WORK.mkdir(parents=True, exist_ok=True)
    b = builder389()
    eff, genes, _ = native_cap15(b)
    receipt = emit_arm(args.arm, eff, genes)
    if args.arm in ("cap15", "k1"):
        ka = known_answer(receipt)
        print("known-answer:", json.dumps(ka, indent=2))
        if not ka["pass"]:
            raise AssertionError(f"{args.arm} did not reproduce 323; "
                                 f"see work/known-answer-{args.arm}.json")
    gc.collect()


if __name__ == "__main__":
    main()
