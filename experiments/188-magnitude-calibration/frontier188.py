"""188 — the consolidated screen: every magnitude SHAPE, on the corrected `nmae` pool.

Two corrections force this table to be recomputed rather than read off `budget188.py`:

  1. `pds` must be the mean over the THREE challenge frames A/B/C, which is 162's
     pre-committed predictor and the quantity 169's `+0.019452` is.  `budget188.py` used frame
     A alone.
  2. `nmae*` must be scored on a pool that resembles the scorer's own DE call.  `pool188.py`
     measures on h1-test, where the reference's BH<0.05 set is on disk, that the call is
     reproduced 0.7219 by `|lfc|*sqrt(cpm)` at that row's own n, 0.6406 by `|h|`, and only
     **0.3796** by `|lfc|` — and that its median size is **89**, matching audit 181's
     board-implied median `n_real` of 90-110.  The `top-2000 |lfc|` pool `budget188.py` used is
     the wrong construction AND the wrong size, and it inverts the sign of the L1 budget term.

KNOWN-ANSWER ROWS
  * `identity` reproduces 169's cached PLOLO a23 rows on all four folds and all three frames;
  * `q1.15` must read **+0.019452** scaled and `q1.15L` **+0.017631** (169 §2, §6.4).
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
sys.path.insert(0, HERE)
sys.path.insert(0, E169)
from budget188 import FOLDS, PDS_SPAN, deliver, load, pds_of        # noqa: E402
from mapfit import Map, apply_map, apply_map_lfc, to_lfc            # noqa: E402
from normfit188 import apply_norm                                    # noqa: E402
from pool188 import topn, umd                                        # noqa: E402
from tilt188 import sharpen_l1, tilt                                 # noqa: E402
from sharp import sharpen, sharpen_lfc                                # noqa: E402

NMAE_SPAN = 0.60          # AGENTS.md §4 range 0.57-0.63; 169 §7 uses the conservative end
POOL_N = 89               # h1-test's median n_real; audit 181's board-implied 90-110


def main() -> None:
    z = np.load(os.path.join(WORK, "chalframes.npz"))
    cg = z["gA"] & z["gB"] & z["gC"]
    frames = {c: (z[f"b{c}"], z[f"g{c}"]) for c in "ABC"}
    bA, gA = frames["A"]
    expr = np.expm1(bA)
    cpmA = 1e6 * expr / float(expr.sum())
    RAW = {k: Map.from_dict(v) for k, v in
           json.loads(open(os.path.join(WORK,
                                        "budget188.json")).read())["maps"].items()}
    NRM = {k: Map.from_dict(v) for k, v in
           json.loads(open(os.path.join(WORK, "maps_norm.json")).read()).items()}
    p169 = {}
    for r in json.loads(open(os.path.join(E169, "work/plolo169.json")).read()):
        p169[(r["fold"], r["arm"], round(r["q"], 2), r["space"])] = r["pds"]

    def arms(P0, held, b, g):
        """Every arm, built with THIS render frame's own gate — which is what 169's
        `plolo169.py` does (`tf(final, q, gate)` with `gate` the frame's).  Applying the A
        gate in the B and C frames reads 0.8% low on `q1.15` and would fail the gate below."""
        cpm = 1e6 * np.expm1(b) / float(np.expm1(b).sum())
        out = [("identity", lambda: np.where(g[None, :], P0, 0.0))]
        for q in (1.05, 1.10, 1.15, 1.20, 1.25):
            out.append((f"q{q:.2f}", lambda q=q: sharpen(P0, q, g, mode="row")))
            out.append((f"q{q:.2f}L", lambda q=q: sharpen_lfc(P0, q, b, g, mode="row")))
        for t in (-0.20, -0.10, 0.10, 0.20):
            out.append((f"tilt{t:+.2f}", lambda t=t: tilt(P0, g, cpm, t, 1.0)))
            out.append((f"q1.15_tilt{t:+.2f}", lambda t=t: tilt(P0, g, cpm, t, 1.15)))
        out.append(("q1.15_L1norm", lambda: sharpen_l1(P0, 1.15, g)))
        for k in ("median", "abs", "de", "budget", "debudget"):
            out.append((f"map_{k}", lambda k=k: apply_map(P0, RAW[f"h|{k}|{held}"], g)))
            out.append((f"map_{k}L",
                        lambda k=k: apply_map_lfc(P0, RAW[f"lfc|{k}|{held}"], b, g)))
        for k in ("median", "budget"):
            out.append((f"nmap_{k}", lambda k=k: apply_norm(P0, NRM[f"{k}|{held}"], g)))
        return out

    acc, kagate = {}, []
    for held in FOLDS:
        F = load(held)
        sc = F["scored"]
        Tr = F["truth"][:, sc]
        mk = np.zeros(cg.size, bool)
        mk[sc] = True
        c = np.flatnonzero(mk & cg)
        bh = F["b"]
        cpmf = 1e6 * np.expm1(bh[c]) / float(np.expm1(bh).sum())
        LR = to_lfc(F["truth"], bh, gA)[:, c]
        pool = topn(np.abs(LR) * np.sqrt(cpmf)[None, :], np.full(LR.shape[0], POOL_N))
        P0 = F["pred"]
        names = [n for n, _ in arms(P0, held, bA, gA)]
        got = {n: dict(pds=[], U=None, M=None, D=None, rA=0.0, ov=0.0) for n in names}
        for ctx in "ABC":
            b, g = frames[ctx]
            base_top = None
            for name, fn in arms(P0, held, b, g):
                Q = fn()
                Dl = deliver(np.where(g[None, :], Q, 0.0), b, g)
                got[name]["pds"].append(pds_of(Dl[:, sc], Tr))
                kt = min(8000, int(g.sum()))
                tops = [set(np.argsort(-np.abs(Q[i]))[:kt].tolist()) for i in range(3)]
                if name == "identity":
                    base_top = tops
                got[name]["ov"] += float(np.mean([len(x & y) / kt
                                                  for x, y in zip(base_top, tops)])) / 3.0
                if ctx == "A":
                    got[name]["rA"] = float(
                        np.median(np.linalg.norm(Q, axis=1)) /
                        np.median(np.linalg.norm(np.where(g[None, :], P0, 0.0), axis=1)))
                    U, M, D = umd(to_lfc(Dl, b, g)[:, c], LR, pool)
                    got[name].update(U=U, M=M, D=D)
                    if name == "identity":
                        kagate.append((held, got[name]["pds"][0],
                                       p169.get((held, "a23", 1.00, "chalA/delivered"))))
        for name in names:
            r = got[name]
            acc.setdefault(name, []).append(dict(pds=float(np.mean(r["pds"])), U=r["U"],
                                                 M=r["M"], D=r["D"], rA=r["rA"], ov=r["ov"],
                                                 pdsA=r["pds"][0]))
        print(f"  [{held}] {len(acc)} arms", flush=True)

    print(f"\n--- KNOWN-ANSWER GATE: identity vs 169's plolo169.json chalA rows ---")
    worst = max(abs(a - b) for _, a, b in kagate)
    for f, a, b in kagate:
        print(f"  {f:<14} here {a:.6f}  169 {b:.6f}  |d| {abs(a-b):.2e}")
    print(f"  worst {worst:.2e}  {'PASS' if worst < 1e-9 else 'FAIL'}")
    if worst >= 1e-9:
        sys.exit("ABORT")

    base = acc["identity"]
    U0 = np.concatenate([x["U"] for x in base]); M0 = np.concatenate([x["M"] for x in base])
    D0 = np.concatenate([x["D"] for x in base])
    nm0 = float(np.mean(1 + (U0 - 2 * M0) / D0))
    print(f"\n{'=' * 118}\nTHE FRONTIER — PLOLO, 4 folds x 3 challenge frames for `pds`; "
          f"`nmae*` on the wald-89 pool (frame A)")
    print(f"  {'arm':<18}{'d pds raw':>11}{'d pds scaled':>13}{'d nmae*':>10}{'dU/D':>9}"
          f"{'e':>8}{'pds ovr':>10}{'nmae ovr':>10}{'SUM(2)':>10}{'rA':>8}{'ktop ov':>9}")
    rows = {}
    for name, v in acc.items():
        U = np.concatenate([x["U"] for x in v]); M = np.concatenate([x["M"] for x in v])
        nm = float(np.mean(1 + (U - 2 * M) / D0))
        du = float(np.mean((U - U0) / D0)); dm = float(np.mean((M - M0) / D0))
        dp = float(np.mean([a["pds"] - b["pds"] for a, b in zip(v, base)]))
        po = dp / PDS_SPAN * 1.011 / 6
        no = -(nm - nm0) / NMAE_SPAN / 6
        rows[name] = dict(dpds=dp, dpds_scaled=dp / PDS_SPAN, dnmae=nm - nm0, dU=du, dM=dm,
                          pds_ovr=po, nmae_ovr=no, rA=float(np.mean([x["rA"] for x in v])),
                          ov=float(np.mean([x["ov"] for x in v])),
                          per_fold=[a["pds"] - b["pds"] for a, b in zip(v, base)])
        print(f"  {name:<18}{dp:>+11.6f}{dp/PDS_SPAN:>+13.6f}{nm-nm0:>+10.5f}{du:>+9.5f}"
              f"{(dm/du if abs(du) > 1e-12 else float('nan')):>8.2f}{po:>+10.6f}"
              f"{no:>+10.6f}{po+no:>+10.6f}{rows[name]['rA']:>8.5f}{rows[name]['ov']:>9.5f}")

    print(f"\n  169's own values for the two known arms: q1.15 pds +0.019452 scaled, "
          f"q1.15L +0.017631 scaled")
    with open(os.path.join(WORK, "frontier188.json"), "w") as fh:
        json.dump(rows, fh, indent=1, default=float)
    print(f"\nwrote {WORK}/frontier188.json")


if __name__ == "__main__":
    sys.exit(main())
