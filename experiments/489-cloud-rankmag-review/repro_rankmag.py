#!/usr/bin/env python3
"""Independent checks of src/vcc2026/rank_mag.py (dn_k250_zs) on synthetic blocks.

No challenge data. Prints one line per check and writes work/repro.json.
Wilcoxon: rc_emit.wilcoxon_z (scanpy-identical tie-corrected z, the project's
exact_de323 kernel) AND cell_eval2.compute_de(backend='scanpy'), which is the
installed scorer's CPU path (pdex and gpudge are not installed here).
"""
from __future__ import annotations

import json
import sys
import traceback
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/423-rank-preserving-magnitude"))

from vcc2026 import rank_mag as RM  # noqa: E402
from vcc2026 import rc_emit  # noqa: E402
from vcc2026 import zero_repair as ZR  # noqa: E402
import rp423  # noqa: E402

OUT = Path(__file__).resolve().parent / "work" / "repro.json"
FINDINGS: list[dict] = []


def rec(name, ok, **info):
    row = {"check": name, "ok": bool(ok), **info}
    FINDINGS.append(row)
    flag = "PASS" if ok else "FAIL"
    extra = " ".join(f"{k}={v}" for k, v in info.items() if k != "detail")
    print(f"[{flag}] {name} {extra}", flush=True)
    if not ok and "detail" in info:
        print(f"       {info['detail']}", flush=True)
    return bool(ok)


def r1_and_tie(logged, cs, tie):
    """Full scanpy rank sum R1 and tie sum. logged is (G, n), cs is (G, m) sorted.

    Same kernel as rp423.r1_tie / rc_emit.wilcoxon_z (midrank, tie-corrected)."""
    if rp423._R1T is None:
        rp423._R1T = rp423._r1_tie_kernel()
    return rp423._R1T(np.ascontiguousarray(cs), np.ascontiguousarray(tie), np.ascontiguousarray(logged))


def logged_f64(block, keep):
    x = np.asarray(block, dtype=np.float64)
    scale = x.sum(1) / 1e6
    scale[scale == 0] = 1.0
    return np.log1p((x / scale[:, None])[:, keep])


def full_stats(block, cc):
    """R1, tie, z, p, q, lfc on the DE universe. z via rc_emit.wilcoxon_z."""
    keep = np.asarray(cc.keep, dtype=bool)
    logged = np.ascontiguousarray(logged_f64(block, keep).T)  # (G, n)
    if not hasattr(cc, "tie"):
        tie = rc_emit.ctrl_tiesums(np.ascontiguousarray(cc.cs))
    else:
        tie = np.ascontiguousarray(cc.tie)
    R1, TT = r1_and_tie(logged, cc.cs, tie)
    z = rc_emit.wilcoxon_z(np.ascontiguousarray(cc.cs), tie, logged)
    from scipy.stats import false_discovery_control, norm
    p = 2 * norm.sf(np.abs(z))
    q = false_discovery_control(p, method="bh")
    lfc = RM.block_lfc(block, cc)
    return {"R1": R1, "tie": TT, "z": z, "p": p, "q": q, "lfc": lfc}


class CC:
    def __init__(self, keep, cs, ref_lin):
        self.keep = np.asarray(keep, dtype=bool)
        self.cs = np.ascontiguousarray(cs, dtype=np.float64)
        self.ref_lin = np.asarray(ref_lin, dtype=np.float64)
        self.tie = rc_emit.ctrl_tiesums(self.cs)
        self.n_controls = int(cs.shape[1])
        self.min_cpm = 5.0


def controls_from_counts(C):
    """C (m, G) integer counts -> ScorerControl-like object using the module."""
    sc = RM.scorer_control_block(np.asarray(C, dtype=np.int64), min_cpm=5.0, chunk=1000)
    return CC(sc.keep, sc.cs, sc.ref_lin)


# ---------------------------------------------------------------- cell_bounds stress
def section_bounds():
    rng = np.random.default_rng(489)
    n_fail = 0
    n_try = 0
    n_tied = 0
    worst_gap = 0
    examples = []
    for trial in range(400):
        m = int(rng.integers(1, 80))
        n = int(rng.integers(1, 40))
        lib_c = rng.integers(500, 800_000, m)
        # sparse-ish counts so many exact zeros exist
        rate = rng.uniform(0, 0.02, m)
        cnt_c = rng.poisson(rate * lib_c)
        # a few heavy cells
        cnt_c[rng.choice(m, size=max(1, m // 10), replace=False)] = rng.integers(1, 500, size=max(1, m // 10))
        s_c = lib_c / 1e6
        c = np.sort(np.log1p(cnt_c / s_c))
        lib = rng.integers(500, 800_000, n).astype(np.float64)
        x = rng.integers(0, 400, n).astype(np.int64)
        # sprinkle exact structural copies of a control (count, lib) so ties happen
        if trial % 3 == 0:
            j = int(rng.integers(0, m))
            x[0] = int(cnt_c[j])
            lib[0] = float(lib_c[j])
        s = lib / 1e6
        n_try += n
        try:
            xl, xh, lt, rt, tied = RM.cell_bounds(x, s, c)
        except AssertionError as e:
            n_fail += 1
            if len(examples) < 8:
                examples.append(f"raise {e} trial {trial}")
            continue
        n_tied += int(tied.sum())
        # every integer in [xl, xh] (capped sample) must keep lt/rt; neighbours must leave
        for i in range(n):
            lo_i, hi_i = int(xl[i]), int(xh[i])
            if hi_i - lo_i > 12:
                grid = np.unique(np.concatenate([
                    np.arange(lo_i, min(lo_i + 4, hi_i) + 1),
                    np.arange(max(lo_i, hi_i - 4), hi_i + 1),
                    [int(x[i])],
                ]))
            else:
                grid = np.arange(lo_i, hi_i + 1)
            v0_lt, v0_rt = int(lt[i]), int(rt[i])
            for cand in grid:
                v = float(np.log1p(cand / s[i]))
                lt2 = int(np.searchsorted(c, v, "left"))
                rt2 = int(np.searchsorted(c, v, "right"))
                if lt2 != v0_lt or rt2 != v0_rt:
                    n_fail += 1
                    if len(examples) < 8:
                        examples.append(f"inside moved lt/rt cell {i} cand {cand} x {int(x[i])} "
                                        f"{v0_lt, v0_rt} -> {lt2, rt2}")
                    break
            # one below xl, if xl>0 and not a tied freeze that equals x
            if lo_i >= 1 and not (tied[i] and lo_i == int(x[i])):
                cand = lo_i - 1
                v = float(np.log1p(cand / s[i]))
                lt2 = int(np.searchsorted(c, v, "left"))
                rt2 = int(np.searchsorted(c, v, "right"))
                if lt2 == v0_lt and rt2 == v0_rt:
                    # allowed only if cand is still inside; that would mean xl is not tight
                    worst_gap += 1
            if hi_i < int(RM.BIG) and not tied[i]:
                cand = hi_i + 1
                v = float(np.log1p(np.float64(cand) / s[i]))
                lt2 = int(np.searchsorted(c, v, "left"))
                rt2 = int(np.searchsorted(c, v, "right"))
                if lt2 == v0_lt and rt2 == v0_rt:
                    worst_gap += 1
                    if len(examples) < 8:
                        examples.append(f"xh+1 still inside xh {hi_i} lib {lib[i]}")
    rec("cell_bounds random stress", n_fail == 0, trials=400, cells=n_try, tied=n_tied,
        interior_fails=n_fail, loose_edges=worst_gap, detail="; ".join(examples[:6]))


def section_f64_matches_celleval():
    """Per-cell log1p(count/(lib/1e6)) vs cell_eval2 _to_linear + scanpy log1p."""
    from cell_eval2.de_compute import _to_linear, _log1p_view
    rng = np.random.default_rng(7)
    n, g = 30, 12
    lib = rng.integers(3000, 40000, n)
    X = rng.poisson(0.01, (n, g))
    # plant a wide range of counts
    X[:, 0] = rng.integers(0, 50, n)
    X[0, :] = 0
    # force libraries
    # scale rows to exact libraries by dumping remainder on the last gene
    X = X.astype(np.int64)
    X[:, -1] = 0
    row = X.sum(1)
    X[:, -1] = np.maximum(lib - row, 0)
    # if row already exceeds lib, shrink
    over = X.sum(1) > lib
    X[over, -1] = 0
    # recompute: set last gene so sum == lib where possible
    X[:, -1] = 0
    slack = lib - X.sum(1)
    ok = slack >= 0
    X[ok, -1] = slack[ok]
    var = [f"g{i}" for i in range(g)]
    obs = pd.DataFrame({"target_gene": ["non-targeting"] * n})
    a = ad.AnnData(X=X.copy(), obs=obs, var=pd.DataFrame(index=var))
    linear = _to_linear(a, "counts", target_sum=1e6)
    logged = _log1p_view(linear)
    mine = logged_f64(X, np.ones(g, dtype=bool))
    # cell_eval2 linear is count/(lib/1e6); log1p of that
    diff_lin = np.max(np.abs(np.asarray(linear.X, dtype=np.float64) - (X / (X.sum(1) / 1e6)[:, None])))
    diff_log = np.max(np.abs(np.asarray(logged.X, dtype=np.float64) - mine))
    rec("cell_eval2 linear == count/(lib/1e6)", diff_lin == 0.0, max_abs=float(diff_lin))
    rec("cell_eval2 log1p == np.log1p of that", diff_log == 0.0, max_abs=float(diff_log))


def section_f32_boundary():
    """Cells parked on the f64 interval edge: does float32 log1p cross the control?"""
    rng = np.random.default_rng(11)
    cross = 0
    checked = 0
    examples = []
    for trial in range(200):
        m = 40
        lib_c = rng.integers(4000, 60000, m)
        cnt_c = rng.integers(0, 30, m)
        s_c = lib_c / 1e6
        c = np.sort(np.log1p(cnt_c / s_c))
        n = 8
        lib = rng.integers(4000, 60000, n).astype(np.float64)
        x = rng.integers(0, 40, n).astype(np.int64)
        s = lib / 1e6
        try:
            xl, xh, lt, rt, tied = RM.cell_bounds(x, s, c)
        except AssertionError:
            continue
        edges = []
        for i in range(n):
            if tied[i]:
                continue
            edges.append((i, int(xl[i])))
            if int(xh[i]) < int(RM.BIG):
                edges.append((i, int(xh[i])))
        for i, cand in edges:
            checked += 1
            v64 = float(np.log1p(cand / s[i]))
            # float32 path used by f32check429 / rowstore scaled_log1p
            v32 = float(np.log1p(np.float32(cand) * (np.float32(1e6) / np.float32(lib[i]))))
            c32 = c.astype(np.float32).astype(np.float64)  # quantise controls the same way? 
            # fairer: recompute control values in float32 from the same counts
            c32 = np.sort(np.log1p(cnt_c.astype(np.float32) * (np.float32(1e6) / lib_c.astype(np.float32))).astype(np.float64))
            lt64 = int(np.searchsorted(c, v64, "left"))
            rt64 = int(np.searchsorted(c, v64, "right"))
            lt32 = int(np.searchsorted(c32, v32, "left"))
            rt32 = int(np.searchsorted(c32, v32, "right"))
            if lt32 != lt64 or rt32 != rt64:
                cross += 1
                if len(examples) < 6:
                    examples.append(f"trial {trial} lib {lib[i]:.0f} count {cand} f64 ({lt64},{rt64}) f32 ({lt32},{rt32})")
    rec("f32 edge vs f64 interval", True, crossings=cross, edges=checked,
        detail="; ".join(examples) if examples else "none")
    # store the count on the finding via the last record
    FINDINGS[-1]["crossings"] = cross
    FINDINGS[-1]["edges"] = checked
    FINDINGS[-1]["examples"] = examples


# ---------------------------------------------------------------- synthetic production-like block
def make_block(seed=1, n_car=60, Kb=4, n_keep=4, m=400, holder_genes=3):
    """Mass-holder geometry the pass is written for.

    Controls: `n_keep` modest genes (counts 0..11) plus one library-holder gene, at libraries
    18k–26k, and one all-zero sub-gate gene. Max control CPM is 11 * 1e6 / 18000 ≈ 611, and
    the dropout puts exact zeros in the controls, so empty carriers are tied and frozen.

    Carriers: three library classes. Each down gene parks counts on two exact-9000 cells and
    two exact-22000 cells, both above that max, and parks the rest of the library on the
    sub-gate gene. The first two class-0 cells are lib+1, matching `class_cells`' remainder,
    and they carry none of the tested genes.
    """
    rng = np.random.default_rng(seed)
    G = n_keep + 2
    lib_c = rng.integers(18000, 26000, m)
    C = np.zeros((m, G), dtype=np.int64)
    for g in range(n_keep):
        C[:, g] = rng.integers(0, 12, m)
    C[:, -1] = 0
    C[:, -2] = lib_c - C.sum(1)
    if np.any(C[:, -2] <= 0):
        raise RuntimeError("control library remainder negative")
    cc = controls_from_counts(C)
    if bool(cc.keep[-1]) or not bool(np.all(cc.keep[:-1])):
        raise RuntimeError(f"gate keep={cc.keep.tolist()}")
    base = np.array([9000, 14000, 22000])
    lib = np.repeat(base, n_car // 3 + 1)[:n_car].astype(np.int64)
    lib[0] += 1
    lib[1] += 1
    X = np.zeros((n_car, G), dtype=np.int64)
    low = np.flatnonzero(lib == 9000)
    high = np.flatnonzero(lib == 22000)
    # 25 counts at lib 9000 is CPM 2778; 20 at lib 22000 is CPM 909. Both sit above ~611.
    for g in range(holder_genes):
        X[low[:2], g] = 25
        X[high[:2], g] = 20
    if n_keep > holder_genes:
        X[::5, holder_genes] = 3
    X[:, -1] = lib - X.sum(1)
    if np.any(X < 0):
        raise RuntimeError("carrier library oversubscribed")
    B = np.zeros((Kb, G), dtype=np.int64)
    B[:, -1] = 900_000
    for g in range(holder_genes):
        B[:, g] = 5
        B[:, -1] -= 5
    return np.vstack([X, B]), n_car, cc, C


def section_invariants():
    block, n_car, cc, _C = make_block()
    before = full_stats(block, cc)
    out, diag = RM.apply(block, n_car, cc, "dn_k250_zs")
    after = full_stats(out, cc)
    rec("totals bit-identical", np.array_equal(out.sum(0), block.sum(0)))
    rec("libraries bit-identical", np.array_equal(out.sum(1), block.sum(1)))
    rec("ballast rows bit-identical", np.array_equal(out[n_car:], block[n_car:]))
    rec("non-negative", int(out.min()) >= 0, mn=int(out.min()))
    dR = after["R1"] - before["R1"]
    rec("full R1 identical (scanpy midrank sum)", np.all(dR == 0),
        n_changed=int(np.sum(dR != 0)), max_abs=float(np.max(np.abs(dR)) if dR.size else 0))
    # control-relative sum, what the module asserts
    keep = cc.keep
    s = block.sum(1) / 1e6
    rel0 = []
    rel1 = []
    gidx = np.flatnonzero(keep)
    for j, g in enumerate(gidx):
        rel0.append(RM.ranks_of(block[:, g], s, cc.cs[j]).sum())
        rel1.append(RM.ranks_of(out[:, g], s, cc.cs[j]).sum())
    rel0 = np.array(rel0)
    rel1 = np.array(rel1)
    rec("control-relative ranks_of sum identical", np.array_equal(rel0, rel1))
    # per-cell, carriers only, on moved genes
    per_cell_bad = 0
    for j, g in enumerate(gidx):
        a = RM.ranks_of(block[:n_car, g], s[:n_car], cc.cs[j])
        b = RM.ranks_of(out[:n_car, g], s[:n_car], cc.cs[j])
        per_cell_bad += int(np.any(a != b))
    rec("per-cell carrier ranks_of identical", per_cell_bad == 0, genes_changed=per_cell_bad)
    dt = after["tie"] - before["tie"]
    dp = np.log(np.maximum(after["p"], 1e-300)) - np.log(np.maximum(before["p"], 1e-300))
    calls0 = before["q"] < 0.05
    calls1 = after["q"] < 0.05
    rec("tie term may move", True, n_tie_changed=int(np.sum(dt != 0)),
        max_abs_dlogp=float(np.max(np.abs(dp)) if dp.size else 0),
        call_flips=int(np.sum(calls0 != calls1)))
    fin = np.isfinite(before["lfc"]) & np.isfinite(after["lfc"])
    rec("no lfc sign flip", not np.any(np.sign(before["lfc"][fin]) != np.sign(after["lfc"][fin])))
    # down genes that moved did not rise
    moved = np.flatnonzero(np.any(out[:, gidx] != block[:, gidx], axis=0))
    rose = []
    for j in moved:
        if before["lfc"][j] < 0 and after["lfc"][j] > before["lfc"][j] + 1e-9:
            rose.append(int(j))
    rec("moved down genes did not rise", len(rose) == 0, rose=rose,
        n_genes_changed=int(moved.size), n_moved_diag=diag.get("zs_n_moved"))
    # k = 1 is a no-op
    out1, _ = RM.apply(block, n_car, cc, "dn_k100_zs")
    rec("k=1 is bit-identical", np.array_equal(out1, block))
    # determinism
    outb, diagb = RM.apply(block, n_car, cc, "dn_k250_zs")
    rec("deterministic", np.array_equal(out, outb) and diag == diagb)
    # vendor matches rp423 zs
    l = RM.block_lfc(block, cc)
    xr, dr = rp423.postpass(block, n_car, cc, l, rp423.parse_variant("dn_k250_zs"))
    rec("matches rp423.postpass zs", np.array_equal(xr, out))
    # k independence on this mass-holder block
    outs = {}
    fracs = {}
    for name in ("dn_k200_zs", "dn_k250_zs", "dn_k1000_zs"):
        o, d = RM.postpass(block, n_car, cc, l, RM.parse_spec(name))
        outs[name] = o
        fracs[name] = (d["frac_massw"], d["n_moved"], d["counts_moved"])
    same = np.array_equal(outs["dn_k200_zs"], outs["dn_k250_zs"]) and np.array_equal(
        outs["dn_k250_zs"], outs["dn_k1000_zs"])
    rec("k=2, 2.5, 10 identical on this block (capacity-bound)", same, fracs=str(fracs))
    # feasible fraction
    t, fr, ll, dA = RM.postpass(block, n_car, cc, l, RM.parse_spec("dn_k250_zs"))[1]["_fr_by_gene"]
    if t.size:
        w = np.abs(dA) / max(np.abs(dA).sum(), 1e-300)
        rec("feasible fraction", True, n_touch=int(t.size),
            mass_weighted=float((fr * w).sum()), max_fr=float(fr.max()), min_fr=float(fr.min()),
            n_reached_95=int(np.sum(fr >= 0.95)), n_negative_fr=int(np.sum(fr < -0.01)))
    # infeasible: all carriers one library
    b2 = block.copy()
    # collapse carrier libraries by moving subgate... easier to rebuild a 1-class block
    return diag


def section_slack_reaches_target():
    """A mild down gene inside one wide control gap, with room to reach 2× and 2.5× but not 10×.

    Controls are lib 20_000, gene 0 alternating 0 and 40 (CPM 0 and 2_000, mean 1_000). Gene 1 is
    the sub-gate (all zeros). Gene 2 holds the control library and is absent from the block, so
    its arithmetic mean is 0 and `select` skips it (`A > 0` is required).
    """
    m = 300
    C = np.zeros((m, 3), dtype=np.int64)
    C[:, 0] = np.where(np.arange(m) % 2 == 0, 0, 40)
    C[:, 2] = 20_000 - C[:, 0]
    cc = controls_from_counts(C)
    n_low, n_high = 8, 8
    n_car = n_low + n_high
    libs = np.array([10_000] * n_low + [25_000] * n_high, dtype=np.int64)
    X = np.zeros((n_car + 2, 3), dtype=np.int64)
    # CPM 1_200 on the low class, 800 on the high class, both strictly inside (0, 2_000).
    X[:n_low, 0] = 12
    X[n_low:n_car, 0] = 20
    X[:n_car, 1] = libs - X[:n_car, 0]
    X[n_car:, 1] = 800_000
    l = RM.block_lfc(X, cc)
    info = {"lfc0": None if l.size == 0 else float(l[0]), "keep": str(cc.keep.tolist()),
            "ref0": None if not cc.keep[0] else float(cc.ref_lin[0])}
    if not (cc.keep[0] and not cc.keep[1] and np.isfinite(l[0]) and l[0] < 0):
        rec("slack gene is down and gated", False, **info)
        return
    outs, diags, lfcs = {}, {}, {}
    for name, spec in (("k2", "dn_k200_zs"), ("k25", "dn_k250_zs"), ("k10", "dn_k1000_zs")):
        o, d = RM.postpass(X, n_car, cc, l, RM.parse_spec(spec))
        outs[name] = o
        diags[name] = d
        lfcs[name] = float(RM.block_lfc(o, cc)[0])
    rec("slack gene: k=2, 2.5 and 10 are not the same block",
        (not np.array_equal(outs["k2"], outs["k25"])) and (not np.array_equal(outs["k25"], outs["k10"])),
        l_base=float(l[0]), **{f"l_{k}": v for k, v in lfcs.items()},
        target_k2=float(2 * l[0]), target_k25=float(2.5 * l[0]), target_k10=float(10 * l[0]),
        frac_k2=diags["k2"]["frac_massw"], frac_k25=diags["k25"]["frac_massw"],
        frac_k10=diags["k10"]["frac_massw"],
        moved_k2=diags["k2"]["counts_moved"], moved_k25=diags["k25"]["counts_moved"],
        moved_k10=diags["k10"]["counts_moved"])
    st0, st1 = full_stats(X, cc), full_stats(outs["k25"], cc)
    rec("slack gene: full R1 held, tie term recorded",
        np.all(st0["R1"] == st1["R1"]),
        dR1=float(np.max(np.abs(st1["R1"] - st0["R1"]))),
        n_tie_changed=int(np.sum(st0["tie"] != st1["tie"])),
        max_abs_dlogp=float(np.max(np.abs(np.log(np.maximum(st1["p"], 1e-300))
                                          - np.log(np.maximum(st0["p"], 1e-300))))),
        k2_overshoot=float(lfcs["k2"] - 2 * l[0]))
    # overshoot of the k*l target, in lfc units, on the gene that can pay for it
    rec("slack gene k=2.5 lands on the 2.5x target",
        abs(lfcs["k25"] - 2.5 * float(l[0])) < abs(lfcs["k10"] - 2.5 * float(l[0])),
        l_k25=lfcs["k25"], target=float(2.5 * l[0]),
        abs_overshoot=abs(lfcs["k25"] - 2.5 * float(l[0])),
        l_k10=lfcs["k10"], target_k10=float(10 * l[0]))


def _two_gene(ctrl_counts, ctrl_lib=20_000, n_low=10, n_high=10, place_low=4):
    """Gene 0 tested, gene 1 sub-gate (zero in every control), gene 2 holds the control library.

    Carrier libraries sit on gene 1, so gene 2 has block mean 0 and is not selected (`A > 0`)."""
    ctrl_counts = np.asarray(ctrl_counts, dtype=np.int64)
    m = len(ctrl_counts)
    C = np.zeros((m, 3), dtype=np.int64)
    C[:, 0] = ctrl_counts
    C[:, 2] = int(ctrl_lib) - ctrl_counts
    if np.any(C[:, 2] < 0):
        raise RuntimeError("control count exceeds ctrl_lib")
    cc = controls_from_counts(C)
    if bool(cc.keep[1]) or not bool(cc.keep[0]):
        raise RuntimeError(f"two-gene gate {cc.keep.tolist()} ref {cc.ref_lin.tolist()}")
    n_car = n_low + n_high
    libs = np.array([8000] * n_low + [24_000] * n_high, dtype=np.int64)
    X = np.zeros((n_car + 2, 3), dtype=np.int64)
    X[:n_car, 1] = libs
    X[:place_low, 0] = 1
    X[:place_low, 1] -= 1
    X[n_car:, 1] = 900_000
    return X, n_car, cc


def section_zero_repair():
    """Repair places 1 count per chosen carrier. A 1-count can leave its cell only when
    that cell sits strictly below every control value (no control zero, and no control
    with a smaller CPM)."""
    rng = np.random.default_rng(5)
    # Case A: 30% of controls are exact zeros. 1-count cells are frozen.
    ctrl = np.where(rng.random(400) < 0.30, 0, rng.integers(1, 8, 400))
    X, n_car, cc = _two_gene(ctrl, place_low=0)
    X0 = X.copy()
    rep, zd = ZR.zero_repair(X0, n_car, cc.keep, cc.ref_lin, -1,
                             np.random.default_rng(1), floor_log2=-5.0)
    out, diag = RM.apply(rep, n_car, cc, "dn_k250_zs")
    l_rep = RM.block_lfc(rep, cc)
    l_zs = RM.block_lfc(out, cc)
    changed_cols = np.flatnonzero(np.any(out != rep, axis=0)).tolist()
    rec("repair with control zeros: repaired gene counts frozen",
        np.array_equal(out[:, 0], rep[:, 0]),
        keep0=bool(cc.keep[0]), zero_diag=str(zd),
        l_repair=float(l_rep[0]), l_zs=float(l_zs[0]),
        gene0_moved=int(np.abs(out[:, 0] - rep[:, 0]).sum()),
        n_moved_other=diag.get("zs_n_moved"), changed_cols=str(changed_cols))

    # Case B: every control count is 4 at lib 15_000 -> CPM 266.7.
    # A low-lib carrier (lib 8_000) with 1 count has CPM 125, below every control.
    X, n_car, cc = _two_gene(np.full(300, 4), place_low=4)
    l_rep = float(RM.block_lfc(X, cc)[0])
    out, diag = RM.apply(X, n_car, cc, "dn_k250_zs")
    l_zs = float(RM.block_lfc(out, cc)[0])
    st, st2 = full_stats(X, cc), full_stats(out, cc)
    j0 = 0  # gene 0 is the first kept gene if it is gated
    moved_b = int(np.abs(out[:, 0] - X[:, 0]).sum())
    dR1 = float(st2["R1"][j0] - st["R1"][j0]) if cc.keep[0] else None
    rec("no control below the 1-count: zs moves it and R1 holds",
        moved_b > 0 and dR1 == 0 and l_zs < l_rep,
        keep0=bool(cc.keep[0]), ref=float(cc.ref_lin[0]) if cc.keep[0] else None,
        l_before=l_rep, l_after=l_zs, dl=(l_zs - l_rep),
        gene0_l1_moved=moved_b,
        low_before=X[:10, 0].tolist(), low_after=out[:10, 0].tolist(),
        high_before=X[10:20, 0].tolist(), high_after=out[10:20, 0].tolist(),
        dR1=dR1,
        dlogp=float(np.log(max(st2["p"][j0], 1e-300)) - np.log(max(st["p"][j0], 1e-300))) if cc.keep[0] else None)

    # Same regime through the real repair. Control CPM = 2000, so the floor -5 asks for
    # more single counts than there are carriers: every carrier gets a 1-count, and the
    # low-library ones sit below every control.
    X0, n_car, cc = _two_gene(np.full(300, 40), place_low=0)
    rep, zd = ZR.zero_repair(X0, n_car, cc.keep, cc.ref_lin, -1,
                             np.random.default_rng(2), floor_log2=-5.0)
    out2, _ = RM.apply(rep, n_car, cc, "dn_k250_zs")
    low_b = int(rep[:10, 0].sum())
    low_a = int(out2[:10, 0].sum())
    rec("no-zero controls: zero_repair then zs moves low-lib 1-counts",
        low_a < low_b and int(out2[:, 0].sum()) == int(rep[:, 0].sum()),
        zero_diag=str(zd),
        l_repair=float(RM.block_lfc(rep, cc)[0]),
        l_zs=float(RM.block_lfc(out2, cc)[0]),
        gene0_moved=int(np.abs(out2[:, 0] - rep[:, 0]).sum()),
        low_before=rep[:10, 0].tolist(), low_after=out2[:10, 0].tolist(),
        high_before=rep[10:20, 0].tolist(), high_after=out2[10:20, 0].tolist())


def section_lib_plus_one():
    """`_levels` rounds cw/cw.max() to 4 decimals. lib and lib+1 differ by ~1e-4 and split."""
    libs = np.array([9000.0, 9001.0, 14000.0, 22000.0])
    cw = 1.0 / (400 * libs / 1e6)
    groups, gw = RM._levels(cw)
    keys = np.round(cw / cw.max(), 4)
    rec("lib and lib+1 are different _levels groups", len(groups) == 4,
        n_groups=len(groups), keys=str(keys.tolist()), gw=str(np.round(gw, 8).tolist()))
    # both cells above the control max, so the split is the only reason a pair exists
    m = 80
    C = np.zeros((m, 3), dtype=np.int64)
    C[:, 0] = 4
    C[:, 2] = 20_000 - 4
    cc = controls_from_counts(C)
    n_car = 8
    libs_i = np.array([9000] * 4 + [9001] * 4, dtype=np.int64)
    # CPM ~889, above the control value 200. Forty ballast rows of zeros pull the mean
    # below the control, so the gene is down and the only legal pair is lib vs lib+1.
    X = np.zeros((n_car + 40, 3), dtype=np.int64)
    X[:n_car, 0] = 8
    X[:n_car, 1] = libs_i - 8
    X[n_car:, 1] = 500_000
    l = RM.block_lfc(X, cc)
    out, d = RM.postpass(X, n_car, cc, l, RM.parse_spec("dn_k250_zs"))
    rec("lib vs lib+1 above the max: counts move, lfc barely changes",
        int(np.abs(out[:, 0] - X[:, 0]).sum()) > 0 and l[0] < 0,
        l0=float(l[0]), l1=float(RM.block_lfc(out, cc)[0]),
        counts_moved=int(np.abs(out[:, 0] - X[:, 0]).sum()),
        low_after=out[:4, 0].tolist(), high_after=out[4:8, 0].tolist(),
        n_moved=d["n_moved"],
        dlfc=float(RM.block_lfc(out, cc)[0] - l[0]))
    st0, st1 = full_stats(X, cc), full_stats(out, cc)
    j = int(np.flatnonzero(cc.keep)[0])
    rec("lib vs lib+1: R1 held, tie term may move",
        st0["R1"][j] == st1["R1"][j],
        dR1=float(st1["R1"][j] - st0["R1"][j]),
        dtie=float(st1["tie"][j] - st0["tie"][j]),
        dlogp=float(np.log(max(st1["p"][j], 1e-300)) - np.log(max(st0["p"][j], 1e-300))))


def section_sign_guard_and_infeasible():
    """One library class: a down gene has a target and sub-gate budget, but no legal pair."""
    m = 80
    ctrl = np.where(np.arange(m) % 5 == 0, 0, 8)
    C = np.column_stack([ctrl, np.zeros(m, dtype=np.int64), 20_000 - ctrl])
    cc = controls_from_counts(C)
    n_car = 12
    libs = np.full(n_car, 10000, dtype=np.int64)
    X = np.zeros((n_car + 2, 3), dtype=np.int64)
    # four cells just heavy enough to be down in the mean (diluted by eight zeros) and,
    # if the interval allows, movable — but every donor shares one library, so no pair exists
    X[:4, 0] = 2
    X[:n_car, 1] = libs - X[:n_car, 0]
    X[n_car:, 1] = 500_000
    l = RM.block_lfc(X, cc)
    out, d = RM.postpass(X, n_car, cc, l, RM.parse_spec("dn_k250_zs"))
    rec("one library class: infeasible down target leaves the block unchanged",
        np.array_equal(out, X) and l[0] < 0, n_moved=d["n_moved"], n_touch=d["n_touch"],
        frac_zero=d["frac_zero"], l0=float(l[0]))
    # tied zeros: carriers at 0, controls include 0
    X2 = np.zeros_like(X)
    X2[:n_car, 1] = libs
    X2[n_car:, 1] = 500_000
    l2 = RM.block_lfc(X2, cc)
    out2, d2 = RM.postpass(X2, n_car, cc, l2, RM.parse_spec("dn_k250_zs"))
    rec("all-zero carriers stay unchanged when controls contain zero",
        np.array_equal(out2, X2), l0=float(l2[0]), n_moved=d2["n_moved"])


def section_asserts_fire():
    """The final asserts are reachable: a hand-broken subgate pool must raise, not return."""
    block, n_car, cc, _ = make_block(seed=2)
    # Monkeypatch by calling _subgate_zs with an inconsistent delta on a copy — the public
    # apply path should not return a block that violates libraries. We instead check that
    # postpass's own asserts reject a corrupted internal state by driving _subgate_zs.
    X = block.copy()
    sidx = np.flatnonzero(~cc.keep)
    bad = np.zeros(n_car, dtype=np.int64)
    bad[0] = 10**9  # asks for more subgate than exists
    raised = False
    try:
        RM._subgate_zs(X, sidx, bad)
    except AssertionError as e:
        raised = True
        msg = str(e)
    rec("subgate shortfall raises", raised, msg=msg if raised else "")
    # public apply on a normal block does return, and the returned block satisfies the claims
    out, _ = RM.apply(block, n_car, cc, "dn_k250_zs")
    rec("apply returns a block satisfying the four invariants",
        np.array_equal(out.sum(0), block.sum(0)) and np.array_equal(out.sum(1), block.sum(1))
        and np.array_equal(out[n_car:], block[n_car:]) and out.min() >= 0)


def section_celleval_de():
    """Run the installed scorer on one small block, before and after the pass."""
    from cell_eval2.de_compute import compute_de
    block, n_car, cc, C = make_block(seed=4, n_car=40, Kb=2, n_keep=4, m=120, holder_genes=3)
    out, diag = RM.apply(block, n_car, cc, "dn_k250_zs")
    genes = [f"g{i}" for i in range(block.shape[1])]

    def frame(pred, ctrl):
        X = np.vstack([ctrl, pred])
        labels = ["non-targeting"] * len(ctrl) + ["T"] * len(pred)
        a = ad.AnnData(X=np.asarray(X, dtype=np.float64),
                       obs=pd.DataFrame({"target_gene": labels}),
                       var=pd.DataFrame(index=genes))
        df = compute_de(
            a, backend="scanpy", groupby="target_gene", reference="non-targeting",
            mean_calc="arithmetic", epsilon=1e-9, input_type="counts", target_sum=1e6,
            fdr_scope="per_pert", filter_gene_min_cpm_cell=5.0, threads=1,
        )
        df = df.filter(df["target"] == "T").sort("feature")
        return df

    d0 = frame(block, C)
    d1 = frame(out, C)
    # join on feature
    a = d0.to_pandas().set_index("feature")
    b = d1.to_pandas().set_index("feature")
    common = a.index.intersection(b.index)
    dp = np.abs(np.log(np.maximum(a.loc[common, "p_value"].to_numpy(), 1e-300))
                - np.log(np.maximum(b.loc[common, "p_value"].to_numpy(), 1e-300)))
    calls0 = a.loc[common, "p_adj"].to_numpy() < 0.05
    calls1 = b.loc[common, "p_adj"].to_numpy() < 0.05
    dl = np.abs(a.loc[common, "log2_fold_change"].to_numpy() - b.loc[common, "log2_fold_change"].to_numpy())
    # sign of lfc
    s0 = np.sign(a.loc[common, "log2_fold_change"].to_numpy())
    s1 = np.sign(b.loc[common, "log2_fold_change"].to_numpy())
    rec("cell_eval2 scanpy: no call flip", int(np.sum(calls0 != calls1)) == 0,
        n_genes=int(len(common)), call_flips=int(np.sum(calls0 != calls1)),
        max_abs_dlogp=float(dp.max() if len(dp) else 0),
        max_abs_dlfc=float(dl.max() if len(dl) else 0),
        sign_flips=int(np.sum(s0 != s1)),
        n_moved=diag.get("zs_n_moved"))
    # compare our z-p to cell_eval2 p on the BEFORE block (calibration)
    st = full_stats(block, cc)
    # features that are gated
    gidx = np.flatnonzero(cc.keep)
    ours = []
    theirs = []
    for j, g in enumerate(gidx):
        name = genes[g]
        if name in a.index:
            ours.append(st["p"][j])
            theirs.append(float(a.loc[name, "p_value"]))
    if ours:
        ratio = np.array(theirs) / np.maximum(np.array(ours), 1e-300)
        rec("our normal-approx p vs cell_eval2 p (before)", True,
            n=len(ours), median_ratio=float(np.median(ratio)),
            max_abs_log_ratio=float(np.max(np.abs(np.log(np.maximum(ratio, 1e-300))))))


def section_selector_and_sha():
    import hashlib
    path = ROOT / "src/vcc2026/rank_mag.py"
    sha = hashlib.sha256(path.read_bytes()).hexdigest()
    pin = "ed6728327a15f039d89b05bfc8f65281f138242ea914544116d54ec99ecde638"
    rec("module sha equals builder pin", sha == pin, sha=sha[:16])
    # selector: dn keeps only l<0; up keeps only l>0; k parsing
    l = np.array([-2.0, -0.1, 0.0, 0.4, 3.0, np.nan])
    ok = np.array([True, True, True, True, True, False])
    dn = RM.select("dn", l, ok)
    up = RM.select("up", l, ok)
    rec("dn selects only finite negative", bool(dn.tolist() == [True, True, False, False, False, False]))
    rec("up selects only finite positive", bool(up.tolist() == [False, False, False, True, True, False]))
    spec = RM.parse_spec("dn_k250_zs")
    rec("dn_k250_zs parses to k=2.5", spec == {"sel": "dn", "k": 2.5, "k_up": None, "mode": "zs"}, spec=str(spec))
    # builder accepts other selectors; module does too
    ok_variants = []
    for v in ("up_k250_zs", "all_k500_zs", "q80dn_k200_zs"):
        ok_variants.append(RM.parse_spec(v)["sel"])
    rec("unpriced selectors parse", ok_variants == ["up", "all", "q80dn"], got=str(ok_variants))
    # target formula: A_t corresponds to k * lfc
    ref = np.array([10.0])
    A = np.array([10.0 * 2 ** (-0.4)])
    lfc = np.log2((A + 1e-9) / (ref + 1e-9))
    k = 2.5
    A_t = (ref + 1e-9) * np.power(2.0, k * lfc) - 1e-9
    # shrink story: if delivered is 0.6 * design lfc
    design = -0.8
    delivered = 0.6 * design
    A_del = (ref + 1e-9) * 2 ** delivered - 1e-9
    l_del = np.log2((A_del + 1e-9) / (ref + 1e-9))
    A_tgt = (ref + 1e-9) * np.power(2.0, 2.5 * l_del) - 1e-9
    l_tgt = float(np.log2((A_tgt + 1e-9) / (ref + 1e-9))[0])
    rec("target is 2.5x the handed lfc (shrink example)",
        abs(l_tgt - 2.5 * delivered) < 1e-6,
        delivered=float(delivered), target=l_tgt, design=design,
        target_over_design=float(l_tgt / design))


def section_jackknife():
    """Does the pass always raise the block jackknife? Use the builder's formula."""
    TS = 50_000.0

    def jk(block):
        x = np.asarray(block, dtype=np.float64)
        p = x.sum(0)
        s = p.sum()
        lib = x.sum(1)
        live = p > 0
        v = np.log1p(TS * (p[live][None, :] - x[:, live]) / (s - lib)[:, None])
        n = x.shape[0]
        return float((n - 1) / n * ((v - v.mean(0)) ** 2).sum())

    deltas = []
    for seed in range(8):
        block, n_car, cc, _ = make_block(seed=seed)
        out, _ = RM.apply(block, n_car, cc, "dn_k250_zs")
        deltas.append(jk(out) / jk(block) - 1.0)
    deltas = np.array(deltas)
    rec("jackknife relative change on 8 synthetic blocks", True,
        min=float(deltas.min()), median=float(np.median(deltas)), max=float(deltas.max()),
        n_decreased=int(np.sum(deltas < -1e-9)))


def section_gate_ulp():
    """Can count/(lib/1e6) and count*(1e6/lib) put a gene on opposite sides of CPM 5?"""
    rng = np.random.default_rng(13)
    found = None
    for trial in range(2000):
        n = 50
        lib = rng.integers(8000, 40000, n)
        # choose a count vector whose mean CPM is near 5
        # CPM_i = c_i * 1e6 / lib_i
        c = rng.integers(0, 5, n)
        a = (c / (lib / 1e6)).mean()
        b = (c * (1e6 / lib)).mean()
        if (a > 5) != (b > 5):
            found = (float(a), float(b), trial)
            break
        # also try to nudge one count
    rec("emitter vs scorer gate can disagree near CPM 5", True,
        found=found is not None, example=str(found))


def main():
    sections = [
        section_bounds,
        section_f64_matches_celleval,
        section_f32_boundary,
        section_invariants,
        section_slack_reaches_target,
        section_lib_plus_one,
        section_zero_repair,
        section_sign_guard_and_infeasible,
        section_asserts_fire,
        section_celleval_de,
        section_selector_and_sha,
        section_jackknife,
        section_gate_ulp,
    ]
    for fn in sections:
        print(f"\n== {fn.__name__} ==", flush=True)
        try:
            fn()
        except Exception:
            traceback.print_exc()
            rec(fn.__name__ + " crashed", False, detail=traceback.format_exc().splitlines()[-1])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(FINDINGS, indent=2, default=str))
    n_fail = sum(1 for r in FINDINGS if not r["ok"])
    print(f"\n{len(FINDINGS)} checks, {n_fail} failed", flush=True)
    return n_fail


if __name__ == "__main__":
    sys.exit(main())
