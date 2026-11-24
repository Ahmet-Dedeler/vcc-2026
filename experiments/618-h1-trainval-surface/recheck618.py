#!/usr/bin/env python3
"""618: the registered 613d re-check on the fresh draw (TRUSTGATE.md, 50f6708c), run only after the gate's delivery-only
verdict is PASS and the hook runs (hk_c4h / hk_id / hk_kg, 3 seeds) exist.

1. Writes the tv read map: h1-test's reference / key / compose / 435 base stats (seeds 0-2) / 480 gate table ->
   the h1-trainval equivalents (ref_tv, tv/effects.npz + compose.json, g2mde's stats = 435's exact flags, gate-h1_tv).
2. Runs, inside 616's sandbox (map mode; the instruments' code is unchanged):
     613/block613d.py   -> 613d's c4h − base and KGID(v3) − c4h  (h1 row = trainval now)
     613/block613d5.py  -> 613d's ungated v3 (idall3's form) − c4h (h1 row)
     618/score618.py    -> the real emitter's contrasts at 521's regimes (3 seeds)
   Map completeness, measured (616's SANDBOX616_TRACE, 19:33-19:58Z): block613d and block613d5 read exactly these
   H1-target-specific files (323 reference / effects, 435's base stats s0, 480's gate-h1_full_ka). Everything else they
   read is target-independent: the donor tables, 371's gene-level CLEAN4 prior, and the Orion / K562 gates for the other rows.
3. Compares per contrast and regime: PASS iff |613d − real| <= 0.0006 raw nmae wherever >= 20 targets are scored.
   -> work/recheck618.json, and a table on stdout.
    uv run python recheck618.py
"""
import json
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
W = HERE / "work"
SB = ROOT / "experiments/616-clean-orion-folds"
TOL = 0.0006
T435 = "experiments/435-zs-on-coarse/work/h1/stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.npz"
MAP = {
    "experiments/323-current-family-validation/work/reference.npz": str(W / "ref_tv/reference.npz"),
    "experiments/323-current-family-validation/work/effects.npz": str(W / "tv/effects.npz"),
    "experiments/323-current-family-validation/work/compose.json": str(W / "tv/compose.json"),
    "experiments/480-nmae-instrument/work/gate-h1_full_ka.npz": str(W / "gate-h1_tv.npz"),
    **{T435.format(s=s): str(W / f"g2mde/stats-base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse.npz") for s in (0, 1, 2)},
}


def sandbox(tag, script):
    env = dict(os.environ, SANDBOX616_MAP=str(W / "recheck618_map.json"))
    r = subprocess.run([sys.executable, str(SB / "sandbox616.py"), tag, "map", str(script)], env=env, cwd=ROOT,
                       capture_output=True, text=True)
    (W / "logs").mkdir(exist_ok=True)
    (W / f"logs/recheck-{tag}.log").write_text(r.stdout + "\n" + r.stderr)
    if r.returncode != 0:
        raise SystemExit(f"{tag} failed (rc {r.returncode}); see work/logs/recheck-{tag}.log")
    return SB / "work/sandbox" / tag


def main():
    missing = [v for v in MAP.values() if not Path(v).exists()]
    if missing:
        raise SystemExit(f"tv inputs missing: {missing}")
    (W / "recheck618_map.json").write_text(json.dumps(MAP, indent=1))
    d = sandbox("recheck618-613d", ROOT / "experiments/613-block-prior-validation/block613d.py")
    j613d = json.loads((d / "experiments/613-block-prior-validation/work/block613d.json").read_text())["h1"]
    d5 = sandbox("recheck618-613d5", ROOT / "experiments/613-block-prior-validation/block613d5.py")
    t5 = (d5 / "experiments/613-block-prior-validation/work/block613d5.out").read_text()
    m = re.search(r"== h1: f 1: .*?ungated v3−c4h ([-+0-9.e]+).*?f 0\.5: .*?ungated v3−c4h ([-+0-9.e]+).*?"
                  r"f 0\.35: .*?ungated v3−c4h ([-+0-9.e]+)", t5)
    if not m:
        raise SystemExit("could not read 613d5's h1 ungated row")
    pred = {
        "c4h-base": {f: j613d[f]["c4h_minus_base"] for f in ("1.0", "0.5", "0.35")},
        "idkg3-c4h": {f: j613d[f]["kgid_minus_c4h"] for f in ("1.0", "0.5", "0.35")},
        "idall3-c4h": dict(zip(("1.0", "0.5", "0.35"), map(float, m.groups()))),
        "n_scored": {f: j613d[f]["n_scored"] for f in ("1.0", "0.5", "0.35")},
    }
    ds = sandbox("recheck618-score", HERE / "score618.py")
    real = json.loads((ds / "experiments/618-h1-trainval-surface/work/score618.json").read_text())
    rows, ok = [], True
    for c in ("c4h-base", "idkg3-c4h", "idall3-c4h"):
        for f in ("1.0", "0.5", "0.35"):
            rv = real.get(c, {}).get(f, {}).get("nmae")
            pv = pred[c][f]
            n = pred["n_scored"][f]
            if rv is None:
                rows.append(dict(contrast=c, f=f, pred=pv, real=None, n_scored=n, verdict="MISSING"))
                ok = False
                continue
            judged = n >= 20
            good = abs(pv - rv) <= TOL
            ok &= good or not judged
            rows.append(dict(contrast=c, f=f, pred=pv, real=rv, err=pv - rv, n_scored=n,
                             verdict=("PASS" if good else "FAIL") if judged else "not judged (< 20 targets)"))
    for r in rows:
        print(f"{r['contrast']:11s} f {r['f']:4s} 613d {r['pred']:+.5f}  real {r['real'] if r['real'] is None else format(r['real'], '+.5f')}  "
              f"n {r['n_scored']}  {r['verdict']}")
    print(f"613d out-of-panel re-check (|613d − real| <= {TOL}): {'PASS' if ok else 'FAIL'}")
    (W / "recheck618.json").write_text(json.dumps({"pass": ok, "tol": TOL, "rows": rows, "pred": pred}, indent=1))


if __name__ == "__main__":
    main()
