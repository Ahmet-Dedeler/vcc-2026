"""264 -- where does the three-donor mean's delivered amplitude go?  Per donor on ONE fold: raw and whitened row norms,
the whitening amplification (un-preserved median norm after / before), the energy fraction on frame A's CPM gate,
and the delivered norm of each donor's rows alone; the same for the two-donor and three-donor means.
    PYTHONDONTWRITEBYTECODE=1 uv run python diag264.py --fold orion_hct116 -> work/diag264-<fold>.{json,log}"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib264 as Z  # noqa: E402

Y, X, L, plolo = Z.Y, Z.X, Z.L, Z.plolo
ap = argparse.ArgumentParser()
ap.add_argument("--fold", default="orion_hct116")
args = ap.parse_args()
held = args.fold
t0 = time.time()
st = Z.stage3(held)
wc, sc = st["wcols"], st["sc"]
op = st["op"]
bA = plolo.challenge_bctrl("A")
eA = np.expm1(bA)
gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
gate_wc = gA[wc]
held_gate = np.asarray(st["fold"]["gate"], bool)[wc]
lines, out = [], {}


def P(s=""):
    print(s, flush=True)
    lines.append(s)


w = (op.evals + Z.LAM) ** (-Z.GAMMA / 2.0)
order = np.argsort(op.evals)[::-1]


def describe(name, rows, avail=None):
    R = rows[:, wc]
    live = np.linalg.norm(R, axis=1) > 0 if avail is None else avail
    R = R[live]
    nr = np.linalg.norm(R, axis=1)
    coef = R @ op.evecs                                     # coordinates in the eigenbasis
    T = coef * w                                             # whitened (un-preserved), same norm as T @ evecs.T
    nt = np.linalg.norm(T, axis=1)
    e_raw_gate = float((R[:, gate_wc] ** 2).sum() / (R * R).sum())
    e_raw_held = float((R[:, held_gate] ** 2).sum() / (R * R).sum())
    Tw = T @ op.evecs.T
    e_wht_gate = float((Tw[:, gate_wc] ** 2).sum() / (Tw * Tw).sum())
    # energy by eigenvalue decile of the operator (top = the donors' dominant directions)
    e_coef = (coef ** 2).sum(0)
    e_T = (T ** 2).sum(0)
    dec_raw = [float(e_coef[order[k * len(order) // 10:(k + 1) * len(order) // 10]].sum() / e_coef.sum()) for k in range(10)]
    dec_wht = [float(e_T[order[k * len(order) // 10:(k + 1) * len(order) // 10]].sum() / e_T.sum()) for k in range(10)]
    # delivered norm of these rows alone after the preserve (as Operator.apply would scale them), frame A, no mass fix
    full = np.zeros((R.shape[0], rows.shape[1]))
    full[:, wc] = Tw * (np.median(nr) / np.median(nt))
    H = plolo.deliver(full, bA, gA, None)
    nd = np.linalg.norm(H[:, wc], axis=1)
    d = dict(n=int(live.sum()), raw_median=float(np.median(nr)), raw_p90=float(np.percentile(nr, 90)),
             wht_unpreserved_median=float(np.median(nt)), amplification=float(np.median(nt) / np.median(nr)),
             wht_p90_over_median=float(np.percentile(nt, 90) / np.median(nt)), raw_p90_over_median=float(np.percentile(nr, 90) / np.median(nr)),
             energy_on_gateA_raw=e_raw_gate, energy_on_held_gate_raw=e_raw_held, energy_on_gateA_wht=e_wht_gate,
             delivered_median=float(np.median(nd)), delivered_over_effect=float(np.median(nd) / np.median(nr)),
             energy_by_operator_eigen_decile_raw=dec_raw, energy_by_operator_eigen_decile_wht=dec_wht)
    P(f"{name:28s} n {d['n']:3d}  raw med {d['raw_median']:.3f} (p90/med {d['raw_p90_over_median']:.2f})  wht unpreserved med {d['wht_unpreserved_median']:.3f} "
      f"amp x{d['amplification']:.2f} (p90/med {d['wht_p90_over_median']:.2f})  energy on gateA raw {100*e_raw_gate:.1f} % / wht {100*e_wht_gate:.1f} % / held gate raw {100*e_raw_held:.1f} %  "
      f"delivered/effect {d['delivered_over_effect']:.3f}")
    P(f"{'':28s} energy by operator eigen-decile (top->bottom) raw {' '.join(f'{100*x:.0f}' for x in dec_raw)} | wht {' '.join(f'{100*x:.0f}' for x in dec_wht)}")
    return d


P(f"[{held}] donors {st['donors3']}; gate A {int(gA.sum())} genes ({int(gate_wc.sum())} on the operator axis); held gate {int(held_gate.sum())} on axis; {time.time()-t0:.0f}s")
for j, nm in enumerate(st["donors3"]):
    out[nm] = describe(nm, st["rows3"][j], st["avail3"][j])
three = st["avail3"].all(0)
out["mean2 (cached pre)"] = describe("mean2 (cached pre)", st["pre"], three)
out["mean3 (a)"] = describe("mean3 (a)", st["pre_a"], three)
# the three whitened components of the three-donor mean on the three-donor targets, pre-preserve
P("\ncomponent norms inside the whitened three-donor mean (three-donor targets; W is linear, pre-preserve):")
comp = {}
for j, nm in enumerate(st["donors3"]):
    R = st["rows3"][j][three][:, wc] / 3.0
    T = (R @ op.evecs) * w
    comp[nm] = float(np.median(np.linalg.norm(T, axis=1)))
    P(f"  {nm:20s} median |W r_j / 3| = {comp[nm]:.3f}")
out["components"] = comp
# cosine between each donor's whitened row and the whitened three-donor mean row
P("\nmedian cosine of each donor's whitened row with the whitened mean row (three-donor targets):")
Wa = st["whitened_a"][three][:, wc]
for j, nm in enumerate(st["donors3"]):
    c = Z.row_cos(st["wht_rows3"][j][three][:, wc], Wa)
    P(f"  {nm:20s} {np.nanmedian(c):+.3f}")
out["seconds"] = round(time.time() - t0, 1)
out["peak_rss_gib"] = Z.peak_gib()
(Z.WORK / f"diag264-{held}.json").write_text(json.dumps(out, indent=1, default=float))
(Z.WORK / f"diag264-{held}.log").write_text("\n".join(lines) + "\n")
P(f"done {out['seconds']}s peak {out['peak_rss_gib']:.2f} GiB")
