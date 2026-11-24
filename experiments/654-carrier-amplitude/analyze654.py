#!/usr/bin/env python3
"""654 — read the carrier-amplitude arms against PREREGISTRATION.md (b0bbfa20).

Arms: x080 / base / x125 (run654.py cpm <design> <seed> ... --zs dn_k250_zs), seeds 0-2; zs stats.
base s2 is read from 630's stored stats (same construction; ka654 checks s0 vs 537 and s1 vs 630 bit for bit).
1. Bulk pin: per-gene block totals (`sums`) equal to base's -> pds/mse pinned (else the arm is void); harness pds/mse raw.
2. Members at 521's four regimes (score625d's code path: realised order, table521.regime_pool), overall-eq by
   lib521.overall_delta (192's INV), paired per seed vs base.
3. 631's regions (decomp631: lib521.load_state's head / block / uncalled x the BASE's delivered direction, lib521.pools at
   f 1 and 0.35): the nmae numerator share per region (sums to the total nmae) and the mean delivered |lfc| per region.
4. Calls per target (q < 0.05, own target excluded).
5. The preregistered verdicts: P1, P2 (bins at f 0.35 AND 0.5, all seeds the same sign), P3.
usage: uv run python experiments/654-carrier-amplitude/analyze654.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import table521 as TB  # noqa: E402

TAG = "base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
MINE = HERE / "work/h1_cpm"
W630 = ROOT / "experiments/630-rank-neutral-up/work/h1"
ARMS = ("x080", "base", "x125")
SEEDS = (0, 1, 2)
ZR = "--zr" in sys.argv          # Amendment A1: the production-faithful construction (--zero-repair -5)
OUT = "analyze654zr" if ZR else "analyze654"


def _asreg(design, s, kind, suf, ext):
    name = f"{kind}-{TAG.format(s=s).replace('base', design, 1)}{suf}.{ext}"
    p = MINE / name
    if not p.exists() and design == "base":
        p = W630 / name
    return p


def path(design, s, kind="stats", suf="_dnk250zs", ext="npz"):
    if not ZR:
        return _asreg(design, s, kind, suf, ext)
    hits = sorted(MINE.glob(f"{kind}-{design}_s{s}_bps*zr-5*{suf}.{ext}"))   # run496 names them bps_zr-5_prodem_...
    if hits:
        assert len(hits) == 1, hits
        return hits[0]
    if design == "x080":           # A1: zero repair is an exact no-op without zero-total gated genes (checked in load)
        return _asreg(design, s, kind, suf, ext)
    return MINE / f"MISSING-{design}-zr-s{s}"


def load(design, s):
    z = np.load(path(design, s), allow_pickle=True)
    out = {k: np.asarray(z[k], float) for k in ("q", "p", "lfc")}
    if ZR and "zr-5" not in path(design, s).name:
        assert not (out["lfc"] < -30).any(), f"{design} s{s} has zero-total genes: run its zr arm"
    out["sums"] = np.asarray(z["sums"])
    m = path(design, s, "metrics", ext="json")
    out["raw"] = json.loads(m.read_text())["raw"] if m.exists() else None
    return out


S0 = L.load_state(0)
OWN, RL = S0["own"], S0["rl"]
ADJ = L.adjud(S0)
POOLS = {reg: TB.regime_pool(S0, reg)[0] for reg in TB.REGIMES}


def members(run, reg):
    pool = POOLS[reg]
    call = (run["q"] < 0.05) & ~OWN
    St = dict(S0)
    St["oq"], St["op"], St["ol"] = run["q"], run["p"], run["lfc"]
    rr = np.array([L.reach_row(L.realized_order_pos(St, i), pool[i], ADJ[i], run["lfc"][i], RL[i])
                   for i in range(len(S0["targets"]))])
    return dict(fid=float(np.nanmean(L.fid_rows(call, run["lfc"], pool, ADJ, RL))), jac=float(np.mean(L.jac_rows(call, pool))),
                reach=float(np.nanmean(rr)), nmae=float(np.nanmean(L.nmae_rows(run["lfc"], pool, RL))))


def regions_table(run, base, f):
    """decomp631's regions with the BASE's delivered direction, on lib521.nmae_rows' exact conventions (targets with >= 10
    finite-truth pool genes, lfc nan/inf -> 0), so the region shares of a target sum to its nmae exactly and their mean over
    the scored targets is members()'s nmae (asserted). decomp631's own estimator keeps small-pool targets and does not."""
    pool = L.pools(S0, f)
    P = np.asarray(run["lfc"], float)
    Pb = np.nan_to_num(np.asarray(base["lfc"], float), nan=0.0, posinf=0.0, neginf=0.0)
    regs = {}
    for rname, R in (("head", S0["head"]), ("block", S0["block"]), ("uncalled", ~S0["C"])):
        for dname, D in (("UP", Pb > 0), ("DOWN", Pb < 0), ("zero", Pb == 0)):
            regs[f"{rname}|{dname}"] = R & D
    share = {k: [] for k in regs}
    mag = {k: [] for k in regs}
    tot = []
    for i in range(P.shape[0]):
        sel = pool[i] & np.isfinite(RL[i])
        if sel.sum() < 10:
            continue
        den = np.abs(RL[i, sel]).sum()
        if not np.isfinite(den) or den == 0:
            continue
        h = np.nan_to_num(P[i], nan=0.0, posinf=0.0, neginf=0.0)
        err = np.where(sel, np.abs(h - np.nan_to_num(RL[i])), 0.0)
        tot.append(err.sum() / den)
        for k, M in regs.items():
            m = sel & M[i]
            share[k].append(err[m].sum() / den)
            mag[k].extend(np.abs(h[m]).tolist())
    ref = float(np.nanmean(L.nmae_rows(P, pool, RL)))
    assert abs(float(np.mean(tot)) - ref) < 1e-12, (np.mean(tot), ref)
    assert abs(sum(float(np.mean(v)) for v in share.values()) - ref) < 1e-12
    return {k: dict(share=float(np.mean(share[k])), mag=float(np.mean(mag[k])) if mag[k] else float("nan"),
                    cells=len(mag[k]) / max(len(tot), 1)) for k in regs}


def main() -> int:
    lines, res = [], {"arms": {}, "deltas": {}, "regions": {}}

    def say(s=""):
        print(s, flush=True)
        lines.append(s)

    runs = {(a, s): load(a, s) for a in ARMS for s in SEEDS if path(a, s).exists()}
    missing = [(a, s) for a in ARMS for s in SEEDS if (a, s) not in runs]
    if missing:
        say(f"MISSING runs: {missing}")
    say("== 1. bulk pin (sums == base) and harness raw pds / mse, f 1")
    for (a, s), r in sorted(runs.items()):
        b = runs.get(("base", s))
        pin = None if b is None else bool(np.array_equal(r["sums"], b["sums"]))
        raw = r["raw"] or {}
        say(f"  {a:5s} s{s}: sums==base {pin}  pds {raw.get('pds', float('nan')):.6f}  mse {raw.get('mse', float('nan')):.6f}  "
            f"calls/target {float(((r['q'] < 0.05) & ~OWN).sum(1).mean()):.1f}")
        res["arms"][f"{a}|{s}"] = dict(sums_equal_base=pin, raw=raw, calls=float(((r['q'] < 0.05) & ~OWN).sum(1).mean()))
    say("\n== 2. members vs base (paired per seed); raw deltas and overall-eq (192 INV); negative nmae = better")
    for a in ("x080", "x125"):
        for reg in TB.REGIMES:
            rows = []
            for s in SEEDS:
                if (a, s) not in runs or ("base", s) not in runs:
                    continue
                m1, m0 = members(runs[(a, s)], reg), members(runs[("base", s)], reg)
                d = {k: m1[k] - m0[k] for k in m1}
                d["overall"] = L.overall_delta(m1, m0)
                rows.append(d)
                res["deltas"][f"{a}|{reg}|{s}"] = d
            if rows:
                mean = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
                sd = float(np.std([r["overall"] for r in rows], ddof=1)) if len(rows) > 1 else float("nan")
                say(f"  {a} f {reg:4s}: overall {mean['overall']:+.5f} (sd {sd:.5f}, seeds {[round(r['overall'], 5) for r in rows]})  "
                    f"fid {mean['fid']:+.5f} jac {mean['jac']:+.5f} reach {mean['reach']:+.5f} nmae {mean['nmae']:+.5f}")
                res["deltas"][f"{a}|{reg}|mean"] = mean | {"sd_overall": sd, "n": len(rows)}
    say("\n== 3. 631 regions (BASE's delivered direction): nmae share per arm and mean delivered |lfc|, mean over seeds")
    for f in (1.0, 0.35):
        say(f"  f {f:g}: region|dir    cells   share x080 / base / x125      mean|lfc| x080 / base / x125   (|lfc| rel x080, x125)")
        tabs = {a: [regions_table(runs[(a, s)], runs[("base", s)], f) for s in SEEDS if (a, s) in runs and ("base", s) in runs]
                for a in ARMS}
        if not tabs["base"]:
            say("         (no base runs yet)")
            continue
        for key in tabs["base"][0]:
            sh = {a: float(np.mean([t[key]["share"] for t in tabs[a]])) for a in ARMS if tabs[a]}
            mg = {a: float(np.nanmean([t[key]["mag"] for t in tabs[a]])) for a in ARMS if tabs[a]}
            cells = float(np.mean([t[key]["cells"] for t in tabs["base"]]))
            rel = {a: (mg[a] / mg["base"] - 1.0) if mg.get("base") else float("nan") for a in ("x080", "x125") if a in mg}
            say(f"         {key:15s} {cells:6.0f}   {sh.get('x080', np.nan):.4f} / {sh['base']:.4f} / {sh.get('x125', np.nan):.4f}"
                f"      {mg.get('x080', np.nan):.4f} / {mg['base']:.4f} / {mg.get('x125', np.nan):.4f}"
                f"     ({rel.get('x080', np.nan):+.1%}, {rel.get('x125', np.nan):+.1%})")
            res["regions"][f"{f:g}|{key}"] = dict(share=sh, mag=mg, cells=cells, mag_rel=rel)
    say("\n== 4. preregistered verdicts")
    g = lambda a, reg, k="overall": res["deltas"].get(f"{a}|{reg}|mean", {}).get(k, float("nan"))  # noqa: E731
    p1 = g("x125", "0.35") < 0 and g("x125", "0.35", "nmae") > 0 and g("x125", "0.35", "reach") < 0
    say(f"  P1 (x1.25 worse at f 0.35, nmae up, reach down): {'HOLDS' if p1 else 'FAILS'} "
        f"(overall {g('x125', '0.35'):+.5f}, nmae {g('x125', '0.35', 'nmae'):+.5f}, reach {g('x125', '0.35', 'reach'):+.5f})")
    seeds_ok = all(np.sign(res['deltas'].get(f'x080|{reg}|{s}', {}).get('overall', np.nan)) == np.sign(g('x080', reg))
                   for reg in ('0.35', '0.5') for s in SEEDS)
    lo = min(g("x080", "0.35"), g("x080", "0.5"))
    hi = max(g("x080", "0.35"), g("x080", "0.5"))
    if seeds_ok and lo >= 0.0010:
        b = "CANDIDATE (>= +0.0010 at f 0.35 and 0.5, every seed the same sign)"
    elif seeds_ok and hi <= 0:
        b = "CLOSED (<= 0: amplitude is at or past its optimum below x1 too)"
    else:
        b = "INCONCLUSIVE (between, or the seeds disagree)"
    say(f"  P2 sign bet (x0.8 better, +0.0003 ... +0.0015 at f 0.35): x0.8 - x1 overall f 0.35 {g('x080', '0.35'):+.5f}, "
        f"f 0.5 {g('x080', '0.5'):+.5f} -> bet {'WON' if g('x080', '0.35') > 0 else 'LOST'}; BIN: {b}")
    hd = res["regions"].get("1|head|DOWN", {}).get("mag", {})
    bu = res["regions"].get("1|block|UP", {}).get("mag", {})
    uu = res["regions"].get("1|uncalled|UP", {}).get("mag", {})
    rel = lambda m: abs(m.get("x125", np.nan) / m.get("x080", np.nan) - 1.0) if m else float("nan")  # noqa: E731
    p3 = rel(hd) < 0.05 and (rel(bu) > 0.10 or rel(uu) > 0.10)
    say(f"  P3 (head DOWN |lfc| moves < 5 % x0.8 -> x1.25; block or uncalled UP > 10 %): {'HOLDS' if p3 else 'FAILS'} "
        f"(head DOWN {rel(hd):.1%}, block UP {rel(bu):.1%}, uncalled UP {rel(uu):.1%}; f 1)")
    res["verdicts"] = dict(P1=p1, P2_bin=b, P2_bet_won=bool(g("x080", "0.35") > 0), P3=p3)
    res["construction"] = "production-faithful (--zero-repair -5)" if ZR else "as registered (no zero repair)"
    (HERE / f"work/{OUT}.json").write_text(json.dumps(res, indent=1, default=float))
    (HERE / f"work/{OUT}.out").write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
