#!/usr/bin/env python3
"""528 — price the real-emitter H1 runs (drive528.sh) on all six members, at 521's board regimes, with target bootstraps,
and measure deliverability (does the emitter make the prescribed block genes significant, by count stratum).

Per arm and seed: the run's own stats (p, q, lfc; pre-zs tag and the zs tag) against 435's production stats of the same
seed (the base; drive528's KA proves the harness reproduces them bit for bit). Members use 521's exact functions (lib521:
fid / jac / nmae / reach rows, the REALIZED order: sig first, q, p, |lfc|, name) on the truth's pools at f 1.0 / 0.5 / 0.35
(PowerSurface) and cmp (table521.cmp_counts). pds / mse come from the harness's own scorer (f 1 only) and must be pinned.
Units: overall = sum_m INV_m / 6 * d_m (192's board-fitted 1/span; nmae sign flipped), `policy521.value`; board columns =
global521.board_cols: fid and jac at 0.72x, reach x 1 | x 0 | x -0.76 (docs/H1-BOARD-CALIBRATION.md's call-set rule).
CIs: 4,000 target resamples of the seed-averaged per-target rows (paired against the seed-averaged base rows).

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/528-callprop-block/price528.py [--ka-only]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
sys.path.insert(0, str(HERE))
import lib521 as L  # noqa: E402
from global521 import board_cols  # noqa: E402
from policy521 import SGN, means, per_target_contrib, top3_share  # noqa: E402
from table521 import regime_pool  # noqa: E402

TAG = "base_s{seed}_bps_prodem_mde0.1cap0.1min0.3coarse"
ST435 = ROOT / "experiments/435-zs-on-coarse/work/h1"
REG = ("1.0", "0.5", "0.35", "cmp")
MEMBERS = ("fid", "jac", "reach", "nmae")
NB = 4000
LIB_MED = 20190.5          # h1-test's median control library (feat528_H1 lib_median)


def stats_path(mode, seed, zs):
    t = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
    return (ST435 if mode == "base" else HERE / f"work/h1_{mode}") / f"stats-{t}.npz"


def metrics_path(mode, seed, zs):
    t = TAG.format(seed=seed) + ("_dnk250zs" if zs else "")
    return (ST435 if mode == "base" else HERE / f"work/h1_{mode}") / f"metrics-{t}.json"


def load_arm(mode, seed, zs):
    p = stats_path(mode, seed, zs)
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=True)
    return dict(p=np.asarray(z["p"], float), q=np.asarray(z["q"], float), lfc=np.asarray(z["lfc"], float),
                sums=np.asarray(z["sums"]), raw=json.load(open(metrics_path(mode, seed, zs)))["raw"])


def rows_at(S, A, pool):
    """(4, T) per-target rows [fid, jac, reach, nmae] of run A on a truth pool, realized order."""
    St = dict(S, oq=A["q"], op=A["p"], ol=A["lfc"])
    call = (A["q"] < L.ALPHA) & ~S["own"]
    adj = L.adjud(S)
    rl = S["rl"]
    fr = L.fid_rows(call, A["lfc"], pool, adj, rl)
    jr = L.jac_rows(call, pool)
    nr = L.nmae_rows(A["lfc"], pool, rl)
    rr = np.array([L.reach_row(L.realized_order_pos(St, i), pool[i], adj[i], A["lfc"][i], rl[i]) for i in range(call.shape[0])])
    return np.vstack([fr, jr, rr, nr])


def val(R, R0, idx=None):
    return float(SGN @ (means(R, idx) - means(R0, idx)))


def boot(R, R0, seed=528):
    rng = np.random.default_rng(seed)
    n = R.shape[1]
    v = np.array([val(R, R0, rng.integers(0, n, n)) for _ in range(NB)])
    return float(np.quantile(v, .025)), float(np.quantile(v, .975)), float(np.mean(v <= 0))


def boot_member(R, R0, m, seed=528):
    rng = np.random.default_rng(seed)
    n = R.shape[1]
    v = []
    for _ in range(NB):
        idx = rng.integers(0, n, n)
        v.append(np.nanmean(R[m, idx]) - np.nanmean(R0[m, idx]))
    return float(np.quantile(v, .025)), float(np.quantile(v, .975))


def parse_mode(mode):
    """(key, keep, decoupled, n_cpm, hook) of a run528 mode; `base` is 435's production run."""
    import re
    if mode == "base":
        return "cpm", 0, False, 2200, False
    m = re.fullmatch(r"([a-z0-9]+)(?:-k([0-9]+))?(-dec|-kept)?(?:-n([0-9]+))?(\+mzg)?", mode)
    return m.group(1), int(m.group(2) or 0), bool(m.group(3)), int(m.group(4) or 2200), bool(m.group(5))


def base_of(mode):
    """the arm's comparison run: 435's production run, or 528's own cpm run with the same size / hook."""
    key, keep, dec, n, hook = parse_mode(mode)
    if n == 2200 and not hook:
        return "base"
    return "cpm" + (f"-n{n}" if n != 2200 else "") + ("+mzg" if hook else "")


def block_of(mode, S, keys):
    """the prescribed (CALLED) block of an arm on H1 per target: the CPM block of its size for cpm/base, else
    block_from_key (the -dec arms call the key's block too; only their steered block stays on CPM)."""
    from k_emit528 import block_from_key  # noqa: PLC0415
    key, keep, dec, n, hook = parse_mode(mode)
    T = S["head"].shape[0]
    if key == "cpm":
        if n == 2200:
            return S["block"]
        return np.array([block_from_key(S["head"][i], S["ref_lin"], n, S["ref_lin"], n) for i in range(T)])
    return np.array([block_from_key(S["head"][i], S["ref_lin"], n, keys[key][i], keep) for i in range(T)])


def deliverability(S, A, block, base_block):
    """realized call rate of the prescribed genes by class and control-count stratum; unsuppressible rate outside."""
    own = S["own"]
    D = (A["q"] < L.ALPHA) & ~own
    head = S["head"] & ~own
    blk = block & ~own
    C = head | blk
    cnt = S["ref_lin"] * LIB_MED / 1e6
    edges = [0, 1.5, 2.0, 3.0, 5.0, 10.0, np.inf]
    cls = {"head": head, "block kept": blk & base_block, "block swapped-in": blk & ~base_block,
           "outside (complement)": ~C & ~own, "swapped-out (now complement)": base_block & ~blk & ~S["head"] & ~own}
    out = {}
    for name, M in cls.items():
        n = int(M.sum())
        r = dict(n_per_target=n / M.shape[0], call_rate=float((D & M).sum() / max(n, 1)))
        strata = {}
        for lo, hi in zip(edges[:-1], edges[1:]):
            g = (cnt >= lo) & (cnt < hi)
            Mg = M & g[None, :]
            if Mg.sum():
                strata[f"[{lo},{hi})"] = dict(n=int(Mg.sum()), call_rate=float((D & Mg).sum() / Mg.sum()))
        r["by_count_per_cell"] = strata
        out[name] = r
    out["missed_per_target"] = float((C & ~D).sum() / C.shape[0])
    out["unsuppressible_per_target"] = float((D & ~C).sum() / C.shape[0])
    return out, C


def sign_acc(S, A, M):
    """share of genes in M whose delivered lfc sign matches the truth's (adjudicable genes); and on the truth's gate."""
    adj = L.adjud(S)
    ok = (np.sign(A["lfc"]) == np.sign(S["rl"])) & (A["lfc"] != 0)
    gate = (S["rq"] < L.ALPHA) & np.isfinite(S["rl"]) & (S["rl"] != 0) & ~S["own"]
    a = M & adj
    g = M & gate
    return dict(n_per_target=float(a.sum() / M.shape[0]), acc_all=float(ok[a].sum() / max(a.sum(), 1)),
                n_gate_per_target=float(g.sum() / M.shape[0]), acc_gate=float(ok[g].sum() / max(g.sum(), 1)))


def main():
    ka_only = "--ka-only" in sys.argv
    states = {s: L.load_state(s) for s in (0, 1, 2)}
    pools = {(s, r): regime_pool(states[s], r)[0] for s in (0, 1, 2) for r in REG}
    keys = dict(np.load(HERE / "work/keys528.npz", allow_pickle=False))
    out = {"KA": {}}
    # KA: 521's realized-order members on 435's own stats == its metrics json, and the cpm run == 435 exactly
    for s in (0, 1, 2):
        A = load_arm("base", s, False)
        R = rows_at(states[s], A, pools[(s, "1.0")])
        got = dict(zip(MEMBERS, (float(np.nanmean(R[m])) for m in range(4))))
        ok = all(abs(got[m] - A["raw"][m]) < 1e-12 for m in MEMBERS)
        out["KA"][f"base_s{s}_members_eq_435_metrics"] = ok
        print(f"KA base s{s}: members from stats == 435 metrics json: {ok}  {got}", flush=True)
        assert ok
    for zs in (False, True):
        A, B = load_arm("cpm", 0, zs), load_arm("base", 0, zs)
        if A is not None:
            eq = all(np.array_equal(A[k], B[k], equal_nan=True) for k in ("p", "q", "lfc", "sums"))
            out["KA"][f"cpm_s0{'_zs' if zs else ''}_eq_435"] = eq
            print(f"KA cpm s0 {'zs' if zs else 'pre-zs'} == 435: {eq}", flush=True)
    if ka_only:
        return
    arms = sorted(p.name[3:] for p in (HERE / "work").glob("h1_*")
                  if p.is_dir() and not p.name[3:].startswith("cpm") and not p.name.endswith("-v3check"))
    res = {}
    for arm in arms:
        for zs in (False, True):
            b0 = base_of(arm)
            seeds = [s for s in (0, 1, 2) if load_arm(arm, s, zs) is not None and load_arm(b0, s, zs) is not None]
            if not seeds:
                continue
            lab = arm + (" +zs" if zs else "")
            r_arm = {"seeds": seeds, "base": b0}
            # pds / mse pinned?
            pm = {}
            for s in seeds:
                A, B = load_arm(arm, s, zs), load_arm(b0, s, zs)
                pm[s] = dict(pds=A["raw"]["pds"] - B["raw"]["pds"], mse=A["raw"]["mse"] - B["raw"]["mse"],
                             sums_equal=bool(np.array_equal(A["sums"], B["sums"])))
            r_arm["pds_mse"] = pm
            for reg in REG:
                per_seed = []
                RA, RB = [], []
                for s in seeds:
                    A, B = load_arm(arm, s, zs), load_arm(b0, s, zs)
                    Ra, Rb = rows_at(states[s], A, pools[(s, reg)]), rows_at(states[s], B, pools[(s, reg)])
                    dm = means(Ra) - means(Rb)
                    per_seed.append(dict(value=val(Ra, Rb), d=dict(zip(MEMBERS, dm.tolist())), board=list(board_cols(dm))))
                    RA.append(Ra)
                    RB.append(Rb)
                Ra, Rb = np.nanmean(np.stack(RA), 0), np.nanmean(np.stack(RB), 0)
                dm = np.array([np.mean([p["d"][m] for p in per_seed]) for m in MEMBERS])
                lo, hi, p0 = boot(Ra, Rb)
                ci_m = {m: boot_member(Ra, Rb, i) for i, m in enumerate(MEMBERS)}
                fj = float(SGN[0] * dm[0] + SGN[1] * dm[1])
                contrib = per_target_contrib(Ra, Rb)
                r_arm[reg] = dict(per_seed=per_seed, value_mean=float(np.mean([p["value"] for p in per_seed])),
                                  top3_share=top3_share(contrib), n_targets_gaining=int((contrib > 0).sum()),
                                  value_sd=float(np.std([p["value"] for p in per_seed])), ci=[lo, hi], p_le0=p0,
                                  d_mean=dict(zip(MEMBERS, dm.tolist())), d_ci=ci_m, fid_jac_overall=fj,
                                  board=list(board_cols(dm)))
                print(f"{lab:16s} {reg:5s} overall {r_arm[reg]['value_mean']:+.5f} (seeds {', '.join(f'{p['value']:+.5f}' for p in per_seed)}) "
                      f"CI [{lo:+.5f}, {hi:+.5f}] P<=0 {p0:.3f} | fid {dm[0]:+.5f} jac {dm[1]:+.5f} reach {dm[2]:+.5f} "
                      f"nmae {dm[3]:+.5f} | fid+jac {fj:+.5f} | board [{r_arm[reg]['board'][0]:+.5f} | {r_arm[reg]['board'][1]:+.5f} | "
                      f"{r_arm[reg]['board'][2]:+.5f}]", flush=True)
            # deliverability (f-independent), seed 0, pre-zs (zs never changes calls)
            if not zs:
                s = seeds[0]
                A, B = load_arm(arm, s, False), load_arm(b0, s, False)
                blk = block_of(arm, states[s], keys)
                bblk = block_of(b0, states[s], keys) if b0 != "base" else states[s]["block"]
                dl, C = deliverability(states[s], A, blk, bblk)
                dl0, C0 = deliverability(states[s], B, bblk, bblk)
                own = states[s]["own"]
                D = (A["q"] < L.ALPHA) & ~own
                sw_in = blk & ~bblk & ~own
                sw_out = bblk & ~blk & ~states[s]["head"] & ~own
                DB = (B["q"] < L.ALPHA) & ~own
                sa = {"swapped-in (arm, called)": sign_acc(states[s], A, sw_in & D),
                      "swapped-out (base, called)": sign_acc(states[s], B, sw_out & DB),
                      "kept block (arm, called)": sign_acc(states[s], A, blk & bblk & D),
                      "kept block (base, called)": sign_acc(states[s], B, blk & bblk & DB),
                      "head (arm, called)": sign_acc(states[s], A, states[s]["head"] & ~own & D),
                      "head (base, called)": sign_acc(states[s], B, states[s]["head"] & ~own & DB)}
                # the cost of calls that do not materialise: jac / fid rows with the prescribed set delivered exactly
                pool = pools[(s, "1.0")]
                adj = L.adjud(states[s])
                jr_real, jr_pres = L.jac_rows(D, pool), L.jac_rows(C & ~own, pool)
                fr_real = L.fid_rows(D, A["lfc"], pool, adj, states[s]["rl"])
                fr_pres = L.fid_rows(C & ~own, A["lfc"], pool, adj, states[s]["rl"])
                r_arm["deliverability"] = dict(arm=dl, base=dl0, sign_accuracy=sa,
                                               realized_vs_prescribed=dict(jac=[float(np.mean(jr_real)), float(np.mean(jr_pres))],
                                                                           fid=[float(np.nanmean(fr_real)), float(np.nanmean(fr_pres))]))
                print(f"   deliverability {arm} s{s}: missed/target {dl['missed_per_target']:.1f} (base {dl0['missed_per_target']:.1f}), "
                      f"unsuppressible/target {dl['unsuppressible_per_target']:.1f} (base {dl0['unsuppressible_per_target']:.1f}); "
                      f"call rate head {dl['head']['call_rate']:.4f} kept {dl['block kept']['call_rate']:.4f} "
                      f"swapped-in {dl['block swapped-in']['call_rate']:.4f} (n/target {dl['block swapped-in']['n_per_target']:.0f}), "
                      f"swapped-out now {dl['swapped-out (now complement)']['call_rate']:.4f}", flush=True)
                for k, v in dl["block swapped-in"]["by_count_per_cell"].items():
                    print(f"      swapped-in count/cell {k:12s} n {v['n']:6d} call rate {v['call_rate']:.4f}", flush=True)
                for k, v in sa.items():
                    print(f"      sign acc {k:28s} n/target {v['n_per_target']:7.1f} acc {v['acc_all']:.4f} | gate n/target "
                          f"{v['n_gate_per_target']:6.1f} acc {v['acc_gate']:.4f}", flush=True)
            res[lab] = r_arm
    out["arms"] = res
    (HERE / "work/price528.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
