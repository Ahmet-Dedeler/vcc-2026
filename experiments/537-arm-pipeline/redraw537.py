#!/usr/bin/env python3
"""537 — RD: does a bulk-only step re-draw pass 2? (RESULTS §2.6)

§2's B1 pairs are bulk-only board pairs at the same builder seed. They measure the DE re-roll of today's law only if a bulk
step re-draws the carriers; if the two builds shared their draws, B1 would under-read the re-roll.

Mechanism (read from the code): the builder draws every carrier from ONE emission stream for the whole file
(`make_submission_pm_cd4.py`: `pick_rng, emission_rng = make_rngs(SEED)`, then per target `ballast_count(profile,
profile_bulk, …)` sets K and `carrier_libs = libraries[picks][: cells_per_pert - ballast_k]` are emitted from
`emission_rng`); the H1 harness does the same (`run496.ballast_block`, one `emit_rng` for all 85 targets). A bulk step
changes K for some targets; each K change changes how many carrier cells are drawn, so every LATER target's carriers come
from a shifted stream.

Known answer (exact): for a bulk-only step, targets BEFORE the first K change draw their carriers from the same stream
state with the same DE profile and libraries, so the carriers' own jackknife statistic `jk_target` (a function of the
carrier counts) must be identical there; after it, pass 1 (`n_pred_pass1`) and the carriers should differ on most targets.
(`carrier_total` is not a test: every emitted cell carries its library total, so it moves only with K.)
For the DE knobs the same diagnostics show whether pass 1 is shared (`n_pred_pass1` equal on every target).

Measures (f 1 truth pool, the harness's own; per-target member rows from price528.rows_at):
  - K changes per seed (count, first index); carriers (jk_target) and pass 1 (n_pred_pass1) equal before / after it;
  - rho_m: the pooled seed-correlation of the per-target member values between the arm and its base at the SAME seed
    (1 = shared draws, 0 = independent re-draws), for the bulk step and, for comparison, for DE knobs;
  - per-target lfc correlation over genes, same seed vs across seeds;
  - the overall member deltas per seed (their mean is the step's own DE leakage, their sd the noise), 85 targets and
    scaled to the board's 900 target-rows.
usage: uv run python experiments/537-arm-pipeline/redraw537.py   (writes work/redraw537.json)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import price537 as P7  # noqa: E402  (load, stats_path, TAG, P5 = price528, L = lib521, regime_pool, SGN)

PAIRS = [
    ("bulk-only step (ps bulk scale ^1.2)", "cpm-bs1.2", "base"),
    ("DE knob: CLEAN4 hook", "cpm+c4h", "base"),
    ("DE knob: call-propensity block (kept)", "lcap-kept", "base"),
    ("DE knob: K 2,000", "cpm-n1600", "base"),
    ("DE knob: hook on the block", "lcap-kept+c4h", "lcap-kept"),
]
MEM = ("fid", "jac", "reach", "nmae")
SC = np.sqrt(85 / 900)


def diag(mode, seed):
    t = P7.TAG.format(seed=seed)
    if mode == "base":
        return json.loads((P7.ST435 / f"diag-{t}.json").read_text())
    for d in P7.DIRS:
        p = d / f"h1_{mode}" / f"diag-{t}.json"
        if p.exists():
            return json.loads(p.read_text())
    return None


def main():
    seeds = (0, 1, 2)
    states = {s: P7.L.load_state(s) for s in seeds}
    pools = {s: P7.regime_pool(states[s], "1.0")[0] for s in seeds}
    rows, lfcs = {}, {}

    def get(mode, s):
        if (mode, s) not in rows:
            A = P7.load(mode, s)
            if A is None:
                return None
            rows[(mode, s)] = P7.P5.rows_at(states[s], A, pools[s])
            lfcs[(mode, s)] = A["lfc"]
        return rows[(mode, s)]

    out = {}
    for label, a, b in PAIRS:
        ss = [s for s in seeds if get(a, s) is not None and get(b, s) is not None and diag(a, s) and diag(b, s)]
        if len(ss) < 2:
            print(f"{label}: {len(ss)} seed(s) -- skipped", flush=True)
            continue
        rec = {"arm": a, "base": b, "seeds": ss, "per_seed": {}}
        first = {}
        for s in ss:
            da, db = diag(a, s), diag(b, s)
            ka = np.array([d["ballast_K"] for d in da])
            kb = np.array([d["ballast_K"] for d in db])
            ja = np.array([d["jk_target"] for d in da])
            jb = np.array([d["jk_target"] for d in db])
            na = np.array([d["n_pred_pass1"] for d in da])
            nb = np.array([d["n_pred_pass1"] for d in db])
            ch = np.flatnonzero(ka != kb)
            f = int(ch[0]) if ch.size else len(ka)
            first[s] = f
            Ra, Rb = rows[(a, s)], rows[(b, s)]
            dm = np.nanmean(Ra, 1) - np.nanmean(Rb, 1)
            rec["per_seed"][s] = dict(
                n_K_changed=int(ch.size), first_K_change=f, K_mean_arm=float(ka.mean()), K_mean_base=float(kb.mean()),
                carriers_equal_before=None if f == 0 else bool(np.array_equal(ja[:f], jb[:f])),
                carriers_equal_after_frac=float((ja[f:] == jb[f:]).mean()) if f < len(ja) else None,
                pass1_equal_before=None if f == 0 else bool(np.array_equal(na[:f], nb[:f])),
                pass1_equal_after_frac=float((na[f:] == nb[f:]).mean()) if f < len(na) else None,
                # the stream-shift signature: after the first K change, targets whose own K did NOT change still re-draw
                sameK_after=int((ka[f:] == kb[f:]).sum()),
                sameK_after_pass1_equal=int(((ka[f:] == kb[f:]) & (na[f:] == nb[f:])).sum()),
                d_members=dict(zip(MEM, dm.tolist())), d_overall=float(P7.SGN @ dm))
        # the seed-correlation of per-target members, same seed (targets after the latest first-K-change: re-drawn rows)
        lo = max(first.values())
        early = min(first.values())
        rho = {}
        for part, idx in (("late", np.arange(lo, 85)), ("early", np.arange(0, early)), ("all", np.arange(85))):
            if idx.size < 3:
                continue
            rho[part] = {"n_targets": int(idx.size)}
            for mi, m in enumerate(MEM):
                A = np.stack([rows[(a, s)][mi, idx] for s in ss])      # seeds x targets
                B = np.stack([rows[(b, s)][mi, idx] for s in ss])
                ok = np.all(np.isfinite(A) & np.isfinite(B), 0)
                A, B = A[:, ok], B[:, ok]
                Ac, Bc = A - A.mean(0), B - B.mean(0)
                cov, va, vb = (Ac * Bc).sum(), (Ac * Ac).sum(), (Bc * Bc).sum()
                rho[part][m] = float(cov / np.sqrt(va * vb)) if va > 0 and vb > 0 else None
        rec["rho_same_seed"] = rho
        # per-target lfc correlation over genes: same seed (arm vs base) vs across seeds (base s vs base s')
        same, cross = [], []
        tt = range(lo, 85) if lo < 85 else range(85)          # DE knobs change no K: all targets
        for s in ss:
            for t in tt:
                x, y = lfcs[(a, s)][t], lfcs[(b, s)][t]
                ok = np.isfinite(x) & np.isfinite(y)
                same.append(np.corrcoef(x[ok], y[ok])[0, 1])
        for s1, s2 in ((ss[0], ss[1]),) + (((ss[1], ss[2]),) if len(ss) > 2 else ()):
            for t in tt:
                x, y = lfcs[(b, s1)][t], lfcs[(b, s2)][t]
                ok = np.isfinite(x) & np.isfinite(y)
                cross.append(np.corrcoef(x[ok], y[ok])[0, 1])
        rec["lfc_corr"] = dict(same_seed_median=float(np.nanmedian(same)), cross_seed_median=float(np.nanmedian(cross)))
        # overall deltas: mean over seeds = the step's own leakage, sd over seeds = noise (85 targets; x sqrt(85/900))
        dov = np.array([rec["per_seed"][s]["d_overall"] for s in ss])
        lev_a = np.array([float(P7.SGN @ np.nanmean(rows[(a, s)], 1)) for s in ss])
        lev_b = np.array([float(P7.SGN @ np.nanmean(rows[(b, s)], 1)) for s in ss])
        rec["overall"] = dict(d_mean=float(dov.mean()), d_sd=float(dov.std(ddof=1)), d_sd_board=float(dov.std(ddof=1) * SC),
                              level_sd_arm=float(lev_a.std(ddof=1)), level_sd_base=float(lev_b.std(ddof=1)),
                              independent_pred_sd=float(np.sqrt(lev_a.var(ddof=1) + lev_b.var(ddof=1))))
        out[label] = rec
        ps = rec["per_seed"]
        print(f"{label}: seeds {ss}", flush=True)
        for s in ss:
            r = ps[s]
            print(f"   s{s}: K changed on {r['n_K_changed']}/85 targets, first at {r['first_K_change']} "
                  f"(mean K {r['K_mean_arm']:.2f} vs {r['K_mean_base']:.2f}); carriers (jk) equal before: {r['carriers_equal_before']}, "
                  f"after: {r['carriers_equal_after_frac']}; pass 1 equal before: {r['pass1_equal_before']}, after: "
                  f"{r['pass1_equal_after_frac']} (same-K targets after it: {r['sameK_after']}, of which pass 1 equal: "
                  f"{r['sameK_after_pass1_equal']}); d overall {r['d_overall']:+.5f} "
                  + " ".join(f"{m} {r['d_members'][m]:+.5f}" for m in MEM), flush=True)
        for part, v in rho.items():
            print(f"   same-seed rho ({part}, {v['n_targets']} targets): " + " ".join(
                f"{m} {v[m]:+.2f}" if v[m] is not None else f"{m} n/a" for m in MEM), flush=True)
        o = rec["overall"]
        print(f"   lfc corr (late targets, or all when no K changed): same seed {rec['lfc_corr']['same_seed_median']:.3f} vs across seeds "
              f"{rec['lfc_corr']['cross_seed_median']:.3f}; overall d {o['d_mean']:+.5f} sd {o['d_sd']:.5f} "
              f"(board {o['d_sd_board']:.5f}); independent re-draw predicts sd {o['independent_pred_sd']:.5f}", flush=True)
    (HERE / "work/redraw537.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
