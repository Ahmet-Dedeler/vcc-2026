"""257 sweep -- every cell of the softer-deflation grid on one PLOLO fold, resolved sizing.

    uv run python sweep257.py --fold <k562_gwps | orion_hct116 | orion_hek293t | kolf2>

Cells (all on top of the shipped stack; the deflated CD4 residual RE-SIZED by `solve_fraction`,
the package's construction, 249 section 2):
  STACK                     the shipped stack (KA: == 246 sweep246 `STACK`, 1e-9)
  a15                       CD4 alone at f = 0.10 (six246's identity reference row; KA in price257)
  k{3,5,8}_s{0.25,0.5,0.75,1}   R_c - s U_k U_k^T R_c            (k5_s1 IS `cd4svd5`)
  k10_s1                    the k = 10 projector (255's arm)
  k5_dir{1..5}              ONE of the top-5 directions removed alone (s = 1)
  randE_<cell>              energy-matched random subspace at that cell's removed energy, 246's
                            `rand` (one seed), resolved sizing; the k5_s1 one also at FIXED sizing
KNOWN-ANSWER rows, asserted at 1e-9 on (cos_agg, a23_cos, pds), the run exits 1 on any failure:
  STACK == sweep246 STACK;  k3_s1 / k5_s1 / k10_s1 == sweep246b B2_cd4_svd{3,5,10}_resolved (and the
  deflated arrays bit-identical to lib246.deflate('svd', k));  randE_k5_s1_fixed == sweep246b
  B2_cd4_randE5 (same dimension);  randE_k5_s1 (resolved) == 249 repro249-<fold> control.
Per row: PLOLO endpoints, 188's exact `nmae` identity (wald89), top-89 precision, delivered norm
(frame A, no mass fix; six246's convention); per-target pds/cos and the identity arrays to npz.
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
sys.path.insert(0, str(HERE))
import lib257 as X  # noqa: E402

L = X.L
TOL = 1e-9
S_GRID = (0.25, 0.5, 0.75, 1.0)
K_GRID = (3, 5, 8)
KEYS = ("cos_agg", "a23_cos", "pds")


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    ref_a = {r["arm"]: r for r in json.loads((X.E246 / f"work/sweep246-{held}.json").read_text())["rows"]}
    ref_b = {r["arm"]: r for r in json.loads((X.E246 / f"work/sweep246b-{held}.json").read_text())["rows"]}
    ref249 = json.loads((X.E249 / f"work/repro249-{held}.json").read_text())

    st = L.stage(held)
    ident = X.Identity(st)
    print(f"[{held}] staged {len(st['targets'])} targets; stack frac_cd4 {st['frac_cd4']:.6f} "
          f"frac_x {st['frac_x']:.6f} ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    failures: list[str] = []

    def run(name, parts, extra=None, ka=None, ka_ref=None):
        t1 = time.time()
        final = L.compose(st, parts)
        e = L.score(st, final)
        m = ident.measure(final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        row = dict(fold=held, arm=name, n_targets=len(st["targets"]),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"],
                   delivered_norm=m["delivered_norm"], audit=extra or {},
                   seconds=round(time.time() - t1, 1))
        msg = (f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
               f"q_w {e['a15_sgn_w']:.6f} nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} "
               f"({row['seconds']}s)")
        if ka_ref is not None:
            d = max(abs(float(ka_ref[k]) - float(e[k])) for k in KEYS)
            ok = d < TOL
            row["ka"] = dict(ref=ka, max_abs_diff=float(d), ok=bool(ok))
            msg += f"   KA {'PASS' if ok else 'FAIL'} vs {ka} ({d:.1e})"
            if not ok:
                failures.append(f"{name} vs {ka}: {d:.3e}")
        rows.append(row)
        print(msg, flush=True)
        return row

    # ------------------------------------------------------------------ reference rows
    run("STACK", L.stack_parts(st), dict(frac_cd4=st["frac_cd4"], frac_x=st["frac_x"], sizing="stack"),
        ka="sweep246 STACK", ka_ref=ref_a["STACK"])
    f010, _ = L.solve_fraction(st["whitened"], st["cd4"], 0.10)
    run("a15", [(st["cd4"], f010)], dict(fraction=float(f010), sizing="cd4 alone f=0.10"))

    # ------------------------------------------------------------------ the grid
    energy: dict[str, float] = {}
    for k in K_GRID + (10,):
        for s in S_GRID:
            if k == 10 and s != 1.0:
                continue
            d, aud = X.deflate_shrink(st["cd4"], k, s)
            if s == 1.0:
                d246, a246 = L.deflate(st["cd4"], "svd", k)
                same = (np.array_equal(np.asarray(d.mean_effect), np.asarray(d246.mean_effect))
                        and aud["energy_removed_frac"] == a246["energy_removed_frac"]
                        and aud["top_k_energy_frac_of_centred"] == a246["top_k_energy_frac_of_centred"])
                aud["bit_identical_to_lib246_svd"] = bool(same)
                print(f"  shrink k={k} s=1 vs lib246.deflate('svd', {k}): "
                      f"{'bit-identical' if same else 'DIFFERS'} (energy removed {aud['energy_removed_frac']:.9f})",
                      flush=True)
                if not same:
                    failures.append(f"shrink k={k} s=1 is not lib246.deflate svd{k}")
            parts, fr, pu = X.resolved_parts(st, d)
            name = f"k{k}_s{s:g}"
            ka = f"B2_cd4_svd{k}_resolved" if (s == 1.0 and k in (3, 5, 10)) else None
            run(name, parts, dict(**aud, fraction=fr, per_unit=pu, sizing="resolved"),
                ka=ka, ka_ref=(ref_b[ka] if ka else None))
            energy[name] = aud["energy_removed_frac"]
            if ka:
                fr_ref = float(ref_b[ka]["audit"]["fraction"])
                if abs(fr - fr_ref) >= 1e-12:
                    failures.append(f"{name}: solved fraction {fr:.12f} != 246 {fr_ref:.12f}")

    # ------------------------------------------------------------------ per-direction cells
    for j in range(5):
        d, aud = X.deflate_shrink(st["cd4"], 5, 1.0, dirs=[j])
        parts, fr, pu = X.resolved_parts(st, d)
        name = f"k5_dir{j+1}"
        run(name, parts, dict(**aud, fraction=fr, per_unit=pu, sizing="resolved"))
        energy[name] = aud["energy_removed_frac"]

    # ------------------------------------------------------------------ controls
    seed = L.RAND_SEED + 5                    # sweep246b: seed_i = 0, k = 5 (246's B2_cd4_randE5 stream)
    for name, efrac in energy.items():
        d, aud, kk = X.rand_energy_matched(st["cd4"], efrac, seed)
        ctl = dict(control="random_energy_matched", matched_to=name, matched_energy=float(efrac),
                   rand_dim=int(kk), seed=int(seed), **aud)
        if name == "k5_s1":
            run("randE_k5_s1_fixed", [(d, st["frac_cd4"]), (st["union"], st["frac_x"])],
                dict(**ctl, sizing="fixed"), ka="B2_cd4_randE5", ka_ref=ref_b["B2_cd4_randE5"])
            if kk != ref_b["B2_cd4_randE5"]["audit"]["k"]:
                failures.append(f"randE5 dimension {kk} != 246's {ref_b['B2_cd4_randE5']['audit']['k']}")
        parts, fr, pu = X.resolved_parts(st, d)
        run(f"randE_{name}", parts, dict(**ctl, sizing="resolved", fraction=fr, per_unit=pu),
            ka=("repro249 control (resolved)" if name == "k5_s1" else None),
            ka_ref=(ref249["control"]["resolved"] if name == "k5_s1" else None))

    meta = dict(fold=held, n_targets=len(st["targets"]), frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"],
                n_rows=len(rows), failures=failures, seconds=round(time.time() - t0, 1),
                peak_rss_gib=peak_gib(), fold_cache=st["fold"]["audit"].get("cache"))
    (X.WORK / f"sweep257-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=float))
    np.savez_compressed(X.WORK / f"pertarget257-{held}.npz", **per_target)
    np.savez_compressed(X.WORK / f"identity257-{held}.npz", **identity)
    print(f"[{held}] {len(rows)} rows, {len(failures)} KA failures, {meta['seconds']}s, "
          f"peak RSS {meta['peak_rss_gib']:.2f} GiB -> {X.WORK/f'sweep257-{held}.json'}", flush=True)
    for f in failures:
        print(f"  *** {f}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
