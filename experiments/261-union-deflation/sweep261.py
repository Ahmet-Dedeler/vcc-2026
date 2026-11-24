"""261 sweep -- every cell of the union-deflation grid on ONE PLOLO fold, on the `cd4svd5` base.

    uv run python sweep261.py --fold <k562_gwps | orion_hct116 | orion_hek293t | kolf2>

Rows (every union cell rides the cd4svd5 base: CD4 deflated k = 5, RESOLVED at f = 0.34):
  STACK, a15, k5_s0.75          257's known-answer rows (the shipped stack, CD4 alone at 0.10, 257's best cell)
  cd4svd5                       the base == 257 k5_s1 == 246 B2_cd4_svd5_resolved == 258/259 'cd4svd5 <- INCUMBENT'
  u_k0_g0.28                    the un-deflated union re-solved by the carrying-row rule (KA: == cd4svd5 at 0.0)
  u_c0_g0.28                    per-gene centring over the union's carrying rows only (rank-1 mean removed)
  u_k{1,2,3,5,8}_g0.28          THE ARM: top-k target-shared directions of the centred union removed, re-sized at 0.28
  u_k{0,3,5}_g{0.20,0.35}       the dose axis with and without the hub
  u_k{3,5}_g0.28_fixed          the projector at the SHIPPED union fraction (sizing not re-solved)
  randE_k{k}                    energy-matched random subspace (246 `rand`, one seed), resolved sizing
  shufP_k{k}                    shuffled projector (246b: fitted on deranged rows, applied to the real rows)
  split_k{k}[b]                 split-half projector (246c: fitted on one half of the union's rows, applied to
                                the other, both ways; two seeds)
KNOWN-ANSWER rows, asserted at 1e-9 on (cos_agg, a23_cos, pds) and exactly on the fractions; the run exits 1 on
any failure: STACK / a15 / k5_s0.75 / cd4svd5 vs 257's sweep257-<fold>.json (also nmae_wald89, top-89 precision,
delivered norm at 0.0); cd4svd5 vs 258's and 259's fold-<fold>.json (cos_agg, macro, pds, CD4 fraction, union
fraction, carrying rows); the deflated CD4 bit-identical to lib246.deflate('svd', 5); u_k0_g0.28 == cd4svd5 at 0.0.
Per row: PLOLO endpoints, 188's exact `nmae` identity (wald89, frame A), top-89 precision, delivered norm; the
per-target pds/cos arrays and the identity arrays to npz.
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
import lib261 as Y  # noqa: E402

L, X = Y.L, Y.X
TOL = 1e-9
KEYS = ("cos_agg", "a23_cos", "pds")
K_GRID = (1, 2, 3, 5, 8)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    held = ap.parse_args().fold
    t0 = time.time()
    ref257 = {r["arm"]: r for r in json.loads((Y.E257 / f"work/sweep257-{held}.json").read_text())["rows"]}
    ref258 = json.loads((Y.E258 / f"work/fold-{held}.json").read_text())["rows"]["cd4svd5 <- INCUMBENT"]
    ref259 = json.loads((Y.E259 / f"work/fold-{held}.json").read_text())["rows"]["cd4svd5 <- INCUMBENT"]

    st = L.stage(held)
    ident = X.Identity(st)
    targets = list(st["targets"])
    union = st["union"]
    n_carry = int(np.sum(np.asarray(union.target_condition_count) > 0))
    print(f"[{held}] staged {len(targets)} targets; stack frac_cd4 {st['frac_cd4']:.6f} frac_x {st['frac_x']:.6f}; "
          f"union carries {n_carry} rows ({time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB)", flush=True)
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    failures: list[str] = []
    ka_log: list[dict] = []

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka_log.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, ok=bool(ok)))
        if not ok:
            failures.append(f"{label}: got {got!r} want {want!r} |d| {d:.3e}")
        return ok

    def run(name, parts, extra=None, ka=None, ka_ref=None):
        t1 = time.time()
        final = L.compose(st, parts)
        e = L.score(st, final)
        m = ident.measure(final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        row = dict(fold=held, arm=name, n_targets=len(targets),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"],
                   delivered_norm=m["delivered_norm"], audit=extra or {},
                   seconds=round(time.time() - t1, 1))
        msg = (f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
               f"q_w {e['a15_sgn_w']:.6f} nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} ({row['seconds']}s)")
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

    # ------------------------------------------------------------------ 257's known-answer rows
    r_stack = run("STACK", L.stack_parts(st), dict(frac_cd4=st["frac_cd4"], frac_x=st["frac_x"], sizing="stack"),
                  ka="257 STACK", ka_ref=ref257["STACK"])
    f010, _ = L.solve_fraction(st["whitened"], st["cd4"], 0.10)
    r_a15 = run("a15", [(st["cd4"], f010)], dict(fraction=float(f010), sizing="cd4 alone f=0.10"),
                ka="257 a15", ka_ref=ref257["a15"])
    d5, fr5, pu5, aud5 = Y.cd4svd5_parts(st)
    d5b, aud5b = L.deflate(st["cd4"], "svd", 5)
    same = np.array_equal(np.asarray(d5.mean_effect), np.asarray(d5b.mean_effect))
    check("deflated CD4 bit-identical to lib246.deflate('svd', 5)", 1.0 if same else 0.0, 1.0, 0.0)
    r_base = run("cd4svd5", [(d5, fr5), (union, st["frac_x"])],
                 dict(**aud5, fraction=fr5, per_unit=pu5, frac_x=st["frac_x"], union_carrying=n_carry, sizing="resolved"),
                 ka="257 k5_s1", ka_ref=ref257["k5_s1"])
    d75, aud75 = X.deflate_shrink(st["cd4"], 5, 0.75)
    parts75, fr75, pu75 = X.resolved_parts(st, d75)
    r_75 = run("k5_s0.75", parts75, dict(**aud75, fraction=fr75, per_unit=pu75, sizing="resolved"),
               ka="257 k5_s0.75", ka_ref=ref257["k5_s0.75"])
    # the identity / precision / delivered-norm keys and the fractions against 257, exactly
    for mine, theirs in ((r_stack, "STACK"), (r_a15, "a15"), (r_base, "k5_s1"), (r_75, "k5_s0.75")):
        for key in ("nmae_wald89", "top89_precision", "delivered_norm", "a15_sgn_w"):
            check(f"257 {theirs}.{key} == {mine['arm']}", mine[key], ref257[theirs][key], 0.0)
    check("257 k5_s1 CD4 fraction == cd4svd5", fr5, ref257["k5_s1"]["audit"]["fraction"], 0.0)
    check("257 k5_s0.75 CD4 fraction", fr75, ref257["k5_s0.75"]["audit"]["fraction"], 0.0)
    # 258's and 259's cd4svd5 row (the builder's own path, apply_frozen_cd4_residual) at 0.0
    for tag, ref in (("258", ref258), ("259", ref259)):
        for key in ("cos_agg", "a23_cos", "pds", "a15_sgn_w", "top89_precision"):
            check(f"{tag} cd4svd5.{key}", r_base[key], ref[key], 0.0)
        check(f"{tag} cd4svd5.cd4_fraction", fr5, ref["cd4_fraction"], 0.0)
        check(f"{tag} cd4svd5.union_fraction == frac_x", st["frac_x"], ref["union_fraction"], 0.0)
        check(f"{tag} cd4svd5.union_carrying", n_carry, ref["union_carrying"], 0.0)
        check(f"{tag} cd4svd5.delivered_norm_A", r_base["delivered_norm"], ref["delivered_norm_A"], 0.0)
        check(f"{tag} cd4svd5.energy_removed", aud5["energy_removed_frac"], ref["energy_removed"], 0.0)

    # ------------------------------------------------------------------ the union's spectrum
    spec = Y.spectrum(union)
    print(f"[{held}] union: {spec['n_carried']} carrying rows, centred energy {spec['centred_energy_frac_of_raw']:.4f} of raw, "
          f"top shares {[round(x, 4) for x in spec['top_shares_of_centred'][:8]]}, cum {[round(x, 3) for x in spec['cum_shares_of_centred'][:8]]}",
          flush=True)

    # ------------------------------------------------------------------ the grid
    def union_row(name, u, g, aud, sizing="resolved", ka=None, ka_ref=None):
        parts, fu, pu = Y.union_parts(st, d5, fr5, u, g, sizing=sizing)
        return run(name, parts, dict(**aud, g=float(g), fraction=fu, per_unit=pu, frac_x_shipped=st["frac_x"],
                                     fraction_ratio=fu / st["frac_x"], cd4_fraction=fr5, sizing=sizing),
                   ka=ka, ka_ref=ka_ref)

    # k = 0 through the resolved path must BE the base (R1)
    u0, aud0 = L.deflate(union, "none")
    r_u0 = union_row("u_k0_g0.28", u0, Y.H1_G, aud0, ka="cd4svd5 (k=0 identity)", ka_ref=r_base)
    check("u_k0_g0.28 union fraction == shipped frac_x", r_u0["audit"]["fraction"], st["frac_x"], 0.0)
    for key in ("nmae_wald89", "top89_precision", "delivered_norm"):
        check(f"u_k0_g0.28.{key} == cd4svd5", r_u0[key], r_base[key], 0.0)
    same0 = np.array_equal(per_target["u_k0_g0.28|pds"], per_target["cd4svd5|pds"])
    check("u_k0_g0.28 per-target pds bit-identical to cd4svd5", 1.0 if same0 else 0.0, 1.0, 0.0)

    uc, audc = L.deflate(union, "centre")
    union_row("u_c0_g0.28", uc, Y.H1_G, audc)
    deflated: dict[int, tuple] = {}
    for k in K_GRID:
        uk, audk = L.deflate(union, "svd", k)
        deflated[k] = (uk, audk)
        union_row(f"u_k{k}_g0.28", uk, Y.H1_G, audk)
    for g in (0.20, 0.35):
        union_row(f"u_k0_g{g:g}", u0, g, aud0)
        for k in (3, 5):
            union_row(f"u_k{k}_g{g:g}", deflated[k][0], g, deflated[k][1])
    for k in (3, 5):
        union_row(f"u_k{k}_g0.28_fixed", deflated[k][0], Y.H1_G, deflated[k][1], sizing="fixed")

    # ------------------------------------------------------------------ controls (every k)
    seed_rand = Y.RAND_SEED + 5                  # 257 / 246's B2_cd4_randE5 stream, one seed
    for k in K_GRID:
        efrac = deflated[k][1]["energy_removed_frac"]
        dr, audr, kk = X.rand_energy_matched(union, efrac, seed_rand)
        union_row(f"randE_k{k}", dr, Y.H1_G, dict(control="random_energy_matched", matched_to=f"u_k{k}_g0.28",
                                                  matched_energy=float(efrac), rand_dim=int(kk), seed=int(seed_rand), **audr))
        ds, auds = Y.shuffled_projector(union, k, Y.SHUF_SEED + k)
        union_row(f"shufP_k{k}", ds, Y.H1_G, dict(control="shuffled_projector", **auds))
        for seed_i in (0, 1):
            dh, audh = Y.split_half_deflate(union, k, Y.RAND_SEED + 31 * seed_i + k)
            union_row(f"split_k{k}{'' if seed_i == 0 else 'b'}", dh, Y.H1_G,
                      dict(control="split_half_projector", seed_index=seed_i, **audh))

    n_ok = sum(1 for r in ka_log if r["ok"])
    meta = dict(fold=held, n_targets=len(targets), targets=targets, frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"],
                cd4svd5_fraction=fr5, union_carrying=n_carry, union_spectrum=spec, n_rows=len(rows),
                ka_rows=ka_log, ka_pass=n_ok, ka_total=len(ka_log), failures=failures,
                seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib(), fold_cache=st["fold"]["audit"].get("cache"))
    (Y.WORK / f"sweep261-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=float))
    np.savez_compressed(Y.WORK / f"pertarget261-{held}.npz", **per_target)
    np.savez_compressed(Y.WORK / f"identity261-{held}.npz", **identity)
    print(f"[{held}] {len(rows)} rows, known-answer {n_ok}/{len(ka_log)} (+{sum(1 for r in rows if 'ka' in r and r['ka']['ok'])}"
          f"/{sum(1 for r in rows if 'ka' in r)} row KAs), {len(failures)} failures, {meta['seconds']}s, "
          f"peak RSS {meta['peak_rss_gib']:.2f} GiB -> {Y.WORK/f'sweep261-{held}.json'}", flush=True)
    for f in failures:
        print(f"  *** {f}", flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
