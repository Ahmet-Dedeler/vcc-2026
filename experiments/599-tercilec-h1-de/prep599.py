#!/usr/bin/env python3
"""599 prep (no truth is read): the H1 key, each mode's bulk / carried / q, and the EXACT ballast K every run will use.

The ballast count is fixed before any emission: run496.ballast_block takes K = min(16, max(1, ceil(S_c (1/ratio - 1) / LB)))
from ratio = min over live genes of profile_bulk / profile_de and the target's first 384 library picks, and the picks come
from the builder's pick stream alone (b.make_rngs(seed)[0], one 400-draw per target in order). So K per (mode, seed,
target) is computed here exactly, without the emitter and without h1-test's perturbed cells. It is the basis of the
registered prediction.

Known answers (each must hold or nothing is written):
  KA-bulk: mode ka's bulk (cached key, key410.weights 0.3 / 3) == run496.ps_bulk(eff, genes, targets) bit for bit;
  KA-K   : mode ka's K per target == the ballast_K that 435's base runs recorded (diag-base_s{0,1,2}_bps_prodem_mde0.1cap0.1
           min0.3coarse.json), every target, every seed.
Writes work/key599.npz (cache), work/prep599.json, work/prep599.npz (q and K arrays). ~6-7 GB (the reference + operator):
run through gate602.
    uv run python experiments/599-tercilec-h1-de/prep599.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/400-uncalled-knn"))
import common599 as C  # noqa: E402
import run400 as H  # noqa: E402  (read-only)

R, E = H.R, H.E
LB, KMAX, NCELL = 999_000, 16, 400
DIAG435 = ROOT / "experiments/435-zs-on-coarse/work/h1/diag-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.json"
SEEDS = (0, 1, 2)


def log(m):
    print(f"{time.strftime('%H:%M:%S')} {m}", flush=True)


def k_of(prof_de, prof_bk, libs):
    """run496.ballast_block's K, verbatim."""
    live = prof_de > 1e-9
    ratio = float(np.min(prof_bk[live] / prof_de[live]))
    s_c_est = float(np.sum(libs[: NCELL - KMAX]))
    need = s_c_est * (1.0 / max(ratio, 1e-6) - 1.0)
    return int(min(KMAX, max(1, np.ceil(need / LB)))), ratio


def main() -> int:
    b = H.builder389()
    eff, genes, targets = H.native_cap15(b)
    kz = C.fit_key(targets)
    log(f"key: {int(np.isfinite(kz['key']).sum())} finite, carried {int(kz['carried'].sum())}/{len(targets)}, "
        f"n_donors hist {np.bincount(np.asarray(kz['n_donors'], int)).tolist()}")
    sys.path.insert(0, str(ROOT / "experiments/496-uncalled-sign-hook"))
    import run496 as M  # noqa: PLC0415
    ka_bulk = C.mode_bulk("ka", eff, genes, targets, kz)[0]
    ref_bulk = M.ps_bulk(eff, genes, targets)
    ka_bulk_ok = bool(np.array_equal(ka_bulk, ref_bulk))
    log(f"KA-bulk (cached key == run496.ps_bulk): {ka_bulk_ok}")
    if not ka_bulk_ok:
        return 1
    del ref_bulk
    r, cc, emit_c = E.load_reference()
    assert [str(t) for t in r["targets"]] == [str(t) for t in targets], "reference target order != cap15 order"
    eff_de = R.map_axis(eff, genes, r["genes"])
    panel = set(map(str, targets))
    panel_mask = np.array([str(g) in panel for g in r["genes"]])
    prof_de = [b.profile_from_effect(eff_de[i], r["b_ctrl"], min_cpm=10.) for i in range(len(targets))]
    picks = {}
    for s in SEEDS:
        pr = b.make_rngs(s)[0]
        picks[s] = [pr.choice(len(r["libraries"]), NCELL, replace=False) for _ in range(len(targets))]
    out = {"targets": [str(t) for t in targets], "key": {"carried": int(kz["carried"].sum()),
                                                         "n_donors_hist": np.bincount(np.asarray(kz["n_donors"], int)).tolist()},
           "modes": {}, "known_answers": {"bulk": ka_bulk_ok}}
    arrays = {}
    for mode in C.MODES:
        eb, w, carried, c = C.mode_bulk(mode, eff, genes, targets, kz)
        qt, audit = C.mode_q(mode, eb, eff, genes, targets, carried)
        eff_bk = R.map_axis(eb, genes, r["genes"])
        K, ratio = {}, {}
        for s in SEEDS:
            ks, rs = [], []
            for i in range(len(targets)):
                if qt is None:
                    pb = b.profile_from_effect(eff_bk[i], r["b_ctrl"], min_cpm=10.)
                else:
                    pb = b.profile_from_effect(eff_bk[i], r["b_ctrl"], min_cpm=10., sharpen_q=float(qt[i]), protect=panel_mask)
                k, rt = k_of(prof_de[i], pb, r["libraries"][picks[s][i]])
                ks.append(k)
                rs.append(rt)
            K[s], ratio[s] = ks, rs
        rec = {"n_carried": int(carried.sum()), "c": c, "w_quantiles": np.quantile(w, [0, .1, .25, .5, .75, .9, 1]).round(4).tolist(),
               "w_at_lo": int((w[carried] <= C.MODES[mode]["lo"] + 1e-12).sum()), "w_at_hi": int((w[carried] >= C.MODES[mode]["hi"] - 1e-12).sum()),
               "audit": {k: v for k, v in audit.items() if k != "moved_vs_tercile"},
               "moved_vs_tercile": audit.get("moved_vs_tercile"),
               "K_mean": {s: float(np.mean(K[s])) for s in SEEDS}, "K_hist": {s: np.bincount(K[s], minlength=17).tolist() for s in SEEDS}}
        if qt is not None:
            rec["tier_counts"] = np.bincount(C.tiers_of(qt), minlength=3).tolist()
            rec["tier_counts_noncarried"] = np.bincount(C.tiers_of(qt)[~carried], minlength=3).tolist()
            arrays[f"q_{mode}"] = qt
        arrays[f"carried_{mode}"] = carried
        arrays[f"w_{mode}"] = w
        for s in SEEDS:
            arrays[f"K_{mode}_s{s}"] = np.asarray(K[s])
            arrays[f"ratio_{mode}_s{s}"] = np.asarray(ratio[s])
        out["modes"][mode] = rec
        log(f"{mode}: carried {rec['n_carried']}, tiers {rec.get('tier_counts')}, K mean "
            + ", ".join(f"s{s} {rec['K_mean'][s]:.3f}" for s in SEEDS))
    ka_k = {}
    for s in SEEDS:
        d = json.loads(Path(str(DIAG435).format(s=s)).read_text())
        rows = d if isinstance(d, list) else d.get("diags", d.get("per_target", []))
        want = [int(x["ballast_K"]) for x in rows]
        ka_k[s] = want == [int(k) for k in arrays[f"K_ka_s{s}"]]
    out["known_answers"]["K_vs_435"] = ka_k
    log(f"KA-K (mode ka K == 435's recorded ballast_K): {ka_k}")
    if not all(ka_k.values()):
        return 1
    # the comparison the prediction is built on
    cmp = {}
    for a_, b_ in (("tc", "t"), ("t", "s175"), ("tc", "s175")):
        d = {s: (arrays[f"K_{a_}_s{s}"] - arrays[f"K_{b_}_s{s}"]) for s in SEEDS}
        cmp[f"{a_}-{b_}"] = {"dK_mean": {s: float(d[s].mean()) for s in SEEDS},
                             "dK_mean_avg": float(np.mean([d[s].mean() for s in SEEDS])),
                             "n_targets_K_changed": {s: int((d[s] != 0).sum()) for s in SEEDS},
                             "q_changed": int((arrays[f"q_{a_}"] != arrays[f"q_{b_}"]).sum())}
    out["compare"] = cmp
    C.dump_json(C.WORK / "prep599.json", out)
    np.savez(C.WORK / "prep599.npz", **arrays)
    log(f"compare: {json.dumps(cmp)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
