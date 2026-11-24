"""521 — validate the reconstruction on the shipped H1 state (435 production stats, seeds 0-2):
  (1) the reconstructed prescribed set C = head(400 by |key|) | block(2,200 by control CPM) against the realized calls;
  (2) the four DE members with the REALIZED order equal 435's metrics json (known answer, every seed);
  (3) the MODEL order (calls in realized order, then |key| + jitter) reproduces the realized reach: the tail model the
      policy pricer uses for genes that leave the call set.
    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/check521.py
"""
from __future__ import annotations

import json
import sys

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

import lib521 as L  # noqa: E402


def main():
    out = {}
    for seed in (0, 1, 2):
        S = L.load_state(seed)
        sig = S["oq"] < L.ALPHA
        own = S["own"]
        C = S["C"]
        n_sig = (sig & ~own).sum(1)
        uns = (sig & ~C & ~own).sum(1)
        missed = (C & ~sig & ~own).sum(1)
        want = json.load(open(str(L.MET435).format(seed=seed)))["raw"]
        m_real = L.members(S, sig & ~own, S["ol"], 1.0, order="realized")
        m_model = L.members(S, sig & ~own, S["ol"], 1.0, order="model", rows=True)
        m_real_rows = L.members(S, sig & ~own, S["ol"], 1.0, order="realized", rows=True)["rows"]
        ka = {m: abs(m_real[m] - want[m]) < 1e-12 for m in ("fid", "jac", "reach", "nmae")}
        drr = m_model["rows"]["reach"] - m_real_rows["reach"]
        print(f"seed {seed}: n_calls median {np.median(n_sig):.0f} (min {n_sig.min()}, max {n_sig.max()}); prescribed "
              f"{C.sum(1).min()}..{C.sum(1).max()}; unsuppressible median {np.median(uns):.0f} mean {uns.mean():.1f} "
              f"max {uns.max()}; missed median {np.median(missed):.0f} max {missed.max()}")
        print(f"   KA realized order: " + " ".join(f"{m} {m_real[m]:.10f} ({'PASS' if ka[m] else 'FAIL'})" for m in ka))
        print(f"   model order reach {m_model['reach']:.6f} vs realized {m_real['reach']:.6f} (d {m_model['reach'] - m_real['reach']:+.6f}); "
              f"targets differing {int(np.sum(np.abs(np.nan_to_num(drr)) > 1e-12))}, max |d| {np.nanmax(np.abs(drr)):.4f}")
        # head / block sign accuracy on today's calls vs the truth's lfc sign (adjudicable)
        adj = L.adjud(S)
        ok = (np.sign(S["ol"]) == np.sign(S["rl"])) & adj
        for nm, m in (("head", S["head"] & sig), ("block", S["block"] & sig & ~S["head"]), ("uns", sig & ~C)):
            mm = m & adj & ~own
            print(f"   calls in {nm:5s}: {mm.sum() / 85:.0f}/target, sign acc vs truth lfc {ok[mm].sum() / max(mm.sum(), 1):.4f}")
        out[seed] = dict(ka=ka, reach_model=m_model["reach"], reach_real=m_real["reach"], uns_mean=float(uns.mean()),
                         missed_max=int(missed.max()))
        if not all(ka.values()):
            raise SystemExit("known answer failed")
    (L.WORK / "check521.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
