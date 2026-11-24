#!/usr/bin/env python3
"""517 — turn the 09-30 board read (458d vs 446b, 458g vs 458d) into (a) the implied loss of Arc's programme information on
A/B/C, r, on H1's attenuation curves (grid517: Arc with a random fraction r of genes' signs flipped), and (b) the expected
price of the programme-sourced hook arms at that r (hookgrid517 / grid517), before any hook entry exists.

usage: uv run python experiments/517-flex-programme/predict517.py --fid-458d 0.0021 [--fid-458g-vs-458d -0.0003]
          [--nmae-458d -0.0012] [--t-fid 0.67] [--t-hook 0.35,1.76] [--inflation]
  --t-fid      board / H1 fid ratio for a sign prior of unchanged quality (392: 0.67)
  --t-hook     board / H1 nmae band for the hook (until 496a reads: 446bja's 0.35 … the MDE family's 1.76; then
               496a's board Δnmae / −0.00146)
  --inflation  shrink H1's same-line fid curve by the share of Arc's fid-like edge that is shared-control noise (halves517:
               0.0071 of 0.0301)
Reads work/grid517.json and work/hookgrid517.json (3-seed means). Prints r from each input and the hook arms' prices.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def curve(rows, key, **filt):
    rs = sorted({r["r"] for r in rows if r.get("r") is not None})
    out = []
    for r in rs:
        v = [x[key] for x in rows if x.get("r") == r and all(x.get(k) == w for k, w in filt.items())]
        out.append((r, float(np.mean(v))))
    return np.array(out)


def invert(cur, y):
    """r at which the (monotone decreasing in r) curve equals y, linear interpolation; clamped to the grid's range."""
    r, v = cur[:, 0], cur[:, 1]
    if y >= v[0]:
        return 0.0
    if y <= v[-1]:
        return float(r[-1])
    return float(np.interp(-y, -v, r))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fid-458d", type=float, required=True, help="board raw Δfid, 458d − 446b")
    ap.add_argument("--fid-458g-vs-458d", type=float, default=None, help="board raw Δfid, 458g − 458d")
    ap.add_argument("--nmae-458d", type=float, default=None, help="board raw Δnmae, 458d − 446b (a check only)")
    ap.add_argument("--t-fid", type=float, default=0.67)
    ap.add_argument("--t-hook", default="0.35,1.76")
    ap.add_argument("--inflation", action="store_true")
    a = ap.parse_args()
    g = json.loads((HERE / "work/grid517.json").read_text())["rows"]
    h = json.loads((HERE / "work/hookgrid517.json").read_text())
    base = dict(w=1.0, q=0.8, tau=None)
    fid_all = curve(g, "fid", gated=False, **base)
    fid_g33 = curve(g, "fid", gated=True, **base)
    nmae_all = curve(g, "nmae", gated=False, **base)
    shrink = (1 - 0.0071 / 0.0301) if a.inflation else 1.0
    fid_all[:, 1] *= shrink; fid_g33[:, 1] *= shrink
    r_d = invert(fid_all, a.fid_458d / a.t_fid)
    print(f"458d Δfid {a.fid_458d:+.4f} → H1-equivalent {a.fid_458d / a.t_fid:+.4f} → r ≈ {r_d:.3f} (ungated fid curve"
          f"{', control-noise shrunk' if a.inflation else ''}: {np.round(fid_all[:, 1], 4).tolist()} at r {fid_all[:, 0].tolist()})")
    rs = [r_d]
    if a.fid_458g_vs_458d is not None:
        fid_g = a.fid_458d + a.fid_458g_vs_458d
        r_g = invert(fid_g33, fid_g / a.t_fid)
        print(f"458g Δfid vs 446b {fid_g:+.4f} → r ≈ {r_g:.3f} (gated curve)")
        rs.append(r_g)
    if a.nmae_458d is not None:
        # nmae decreases (improves) as quality rises: invert on -nmae
        cur = nmae_all.copy(); cur[:, 1] = -cur[:, 1]
        r_n = invert(cur, -a.nmae_458d / 1.0)
        print(f"458d Δnmae {a.nmae_458d:+.4f} → r ≈ {r_n:.3f} on the block nmae curve at transfer 1 (a check: the block's nmae "
              f"carries a kf magnitude term and 392's nmae read inverted, so fid decides)")
    r = float(np.mean(rs))
    lo, hi = (float(x) for x in a.t_hook.split(","))
    print(f"\nusing r = {r:.3f}")

    def hook_at(tau, gated):
        cur = []
        for rr in (0.0, 0.1, 0.2, 0.3, 0.4):
            v = [x["nmae"] for x in h if x["src"] == "arc" and x["r"] == rr and x["tau"] == tau and x["gated"] == gated]
            cur.append((rr, float(np.mean(v))))
        cur = np.array(cur)
        return float(np.interp(r, cur[:, 0], cur[:, 1]))
    c4 = float(np.mean([x["nmae"] for x in h if x["src"] == "clean4"]))
    print(f"H1 nmae of 496's CLEAN4 hook: {c4:+.5f} (board band {c4 * lo:+.4f} … {c4 * hi:+.4f})")
    for tau in (0.05, 0.1):
        for gated in (True, False):
            v = hook_at(tau, gated)
            print(f"Arc hook tau {tau} {'gated (le1)' if gated else 'ungated':11s} at r {r:.2f}: H1-equivalent nmae {v:+.5f} "
                  f"→ board {v * lo:+.4f} … {v * hi:+.4f} ({'beats' if v < c4 else 'does NOT beat'} the CLEAN4 hook)")
    mzp = HERE / "work/mzgrid517.json"
    if mzp.exists():
        mz = json.loads(mzp.read_text())
        cur = np.array([(rr, float(np.mean([x["nmae"] for x in mz if x["arm"] == "hook_g" and x["r"] == rr])))
                        for rr in (0.0, 0.1, 0.2, 0.3, 0.4)])
        v = float(np.interp(r, cur[:, 0], cur[:, 1]))
        print(f"meanz_ne hook (the recommended source, $MZ) gated at r {r:.2f}: H1-equivalent nmae {v:+.5f} → board "
              f"{v * lo:+.4f} … {v * hi:+.4f} ({'beats' if v < c4 else 'does NOT beat'} the CLEAN4 hook)")
    print("\nBlock prior at this r (H1 overall-equiv, q0.8): "
          f"gated {float(np.interp(r, *curve(g, 'ov', gated=True, **base).T)):+.4f}, "
          f"ungated {float(np.interp(r, *curve(g, 'ov', gated=False, **base).T)):+.4f}")


if __name__ == "__main__":
    main()
