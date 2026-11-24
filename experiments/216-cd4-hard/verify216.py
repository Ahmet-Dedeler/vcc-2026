"""216 Phase 0c — verify 215's Phase-2 windows independently, before building on them.

The mandate is explicit: "Verify its Phase-2 windows yourself before building on them."  This
recomputes them from the raw snapshots with 197's identity and NOTHING from 215's code:

    cos_floor = (1.023611 - raw_mse) / 2.074495          (vacuous above raw_mse 1.023611)

and re-derives 215's amplitude discriminant on the two headline windows, so that the x2.06 /
x2.07 this whole experiment is measured against is a number this experiment owns.

Read-only.  The snapshots in `data/leaderboard/` and 215's `work/lb-live.json` are inputs;
nothing is written outside this experiment and no network call is made.

    uv run python verify216.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
LB = ROOT / "data/leaderboard"
LIVE = ROOT / "experiments/215-leader-decomposition/work/lb-live.json"

VACUOUS_ABOVE = 1.023611
SLOPE = 2.074495
LEAK_P90 = 0.0236          # 215 §1.4 V2, the conservative leak bound |dpds|/|dmse|
RATE_PDS = 0.6479          # 215's one fitted constant

EXPECT = {
    ("cqawesome", "from"): 0.06427, ("cqawesome", "to"): 0.13211,
    ("jw_13288", "from"): 0.06333, ("jw_13288", "to"): 0.13117,
    ("Aginglab.com", "to"): 0.09887,
}


def floor_of(raw_mse):
    return (VACUOUS_ABOVE - raw_mse) / SLOPE


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def rows_of(payload):
    if isinstance(payload, dict):
        for k in ("entries", "leaderboard", "results", "data", "rows", "items"):
            if k in payload and isinstance(payload[k], list):
                return payload[k]
        for v in payload.values():
            if isinstance(v, list) and v and isinstance(v[0], dict):
                return v
        return []
    return payload if isinstance(payload, list) else []


def get(row, *names):
    for n in names:
        if n in row and row[n] is not None:
            return row[n]
        for k in row:
            if k.lower().replace("_", "") == n.lower().replace("_", ""):
                return row[k]
    return None


def main() -> None:
    obs = {}
    files = sorted(LB.glob("*.json")) + ([LIVE] if LIVE.exists() else [])
    for p in files:
        try:
            payload = json.loads(p.read_text())
        except Exception:                                        # noqa: BLE001
            continue
        for r in rows_of(payload):
            if not isinstance(r, dict):
                continue
            team = get(r, "team", "team_name", "teamName", "name", "user", "username")
            model = get(r, "model", "model_name", "modelName", "message", "description")
            mse = num(get(r, "mse_raw", "rawMse", "expr_mse_unbiased_capped_norm", "mse"))
            pds = num(get(r, "pds_raw", "rawPds", "pds_cosine", "pds"))
            if team is None or mse is None or pds is None:
                continue
            # some payloads carry the SCALED mse; the raw one is the large number
            key = (str(team), str(model))
            prev = obs.get(key)
            cand = dict(team=str(team), model=str(model), raw_mse=mse, raw_pds=pds,
                        src=p.name)
            if prev is None or p.name > prev["src"]:
                obs[key] = cand

    teams = {}
    for (team, model), r in obs.items():
        teams.setdefault(team, []).append(r)

    rep = {"n_snapshots": len(files), "n_observations": len(obs), "n_teams": len(teams)}
    checks = []
    for team in ("cqawesome", "jw_13288", "Aginglab.com", "Ahmet Dedeler"):
        rows = teams.get(team, [])
        informative = [r for r in rows if r["raw_mse"] < VACUOUS_ABOVE]
        if not informative:
            checks.append(dict(team=team, error="no informative row found",
                               n_rows=len(rows)))
            continue
        floors = sorted((floor_of(r["raw_mse"]), r["model"], r["raw_pds"], r["raw_mse"])
                        for r in informative)
        lo, hi = floors[0], floors[-1]
        d_pds = hi[2] - lo[2]
        d_mse = hi[3] - lo[3]
        leak = LEAK_P90 * abs(d_mse)
        checks.append(dict(
            team=team, n_rows=len(rows), n_informative=len(informative),
            floor_min=lo[0], floor_min_model=lo[1],
            floor_max=hi[0], floor_max_model=hi[1],
            ratio=hi[0] / lo[0] if lo[0] > 0 else None,
            d_raw_pds=d_pds, d_raw_mse=d_mse,
            conservative_leak_bound=leak,
            times_bound=(abs(d_pds) / leak) if leak > 0 else None,
            rho=((d_pds / RATE_PDS) / (-d_mse / SLOPE)) if d_mse else None))
        print(f"[{team}] informative rows {len(informative)}/{len(rows)}  "
              f"floor {lo[0]:.5f} ({lo[1]}) -> {hi[0]:.5f} ({hi[1]})  "
              f"x{hi[0]/lo[0]:.3f}   d_pds {d_pds:+.5f}  "
              f"x bound {abs(d_pds)/leak if leak else float('nan'):.1f}  "
              f"rho {checks[-1]['rho']:.2f}" if checks[-1].get("rho") is not None else "")

    rep["teams"] = checks
    gates = []
    for c in checks:
        for side, key in (("from", "floor_min"), ("to", "floor_max")):
            e = EXPECT.get((c["team"], side))
            if e is None:
                continue
            ok = abs(c[key] - e) <= 5e-5
            gates.append(dict(team=c["team"], side=side, expected=e,
                              measured=round(c[key], 5), passed=bool(ok)))
            print(f"[{'PASS' if ok else 'FAIL'}] 215 {c['team']} {side}: "
                  f"expected {e} measured {c[key]:.5f}")
    rep["gates"] = gates
    rep["all_passed"] = bool(all(g["passed"] for g in gates)) and bool(gates)

    cq = [c for c in checks if c["team"] == "cqawesome"]
    jw = [c for c in checks if c["team"] == "jw_13288"]
    if cq and jw:
        rep["headline"] = dict(
            cqawesome_ratio=cq[0]["ratio"], jw_ratio=jw[0]["ratio"],
            mean_ratio=float(np.mean([cq[0]["ratio"], jw[0]["ratio"]])),
            leader_move_pct=100.0 * (float(np.mean([cq[0]["ratio"], jw[0]["ratio"]])) - 1.0))
        print(f"\nverified leader move: cqawesome x{cq[0]['ratio']:.4f}, "
              f"jw_13288 x{jw[0]['ratio']:.4f}  (215 published x2.06 / x2.07)")
    (WORK / "verify216.json").write_text(json.dumps(rep, indent=2))


if __name__ == "__main__":
    main()
