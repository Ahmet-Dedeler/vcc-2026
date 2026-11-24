#!/usr/bin/env python3
"""537 — how much of the call-propensity block's swapped-in set the CLEAN4 hook (tau 0.05) reaches, on A/B/C and on H1.

The hook signs complement genes (not head, not the steered block) where |P_down - 0.5| >= 0.05; under the block's `kept`
steering the swapped-in genes are called but unsteered, so the hook reaches them (fid). This compares the coverage on A/B/C
(528's abc artifact; each context's gated 5-CPM axis; top 2,200 by the key vs by control CPM, head ignored) with H1's
(per target, head excluded, 528's keys528 `lcap`, the production CallSet). Prior: the production CLEAN4 file
(396's prior_ctx_challenge.npz P_down) on A/B/C; 371's consensus mean (== it at tau 0.05, 517's KA) on H1.
Writes work/coverage537.json. Light (< 1 GB).
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
from vcc2026.rc_emit_k import block_from_key  # noqa: E402


def cover(P, sel):
    ok = np.isfinite(P) & (np.abs(P - 0.5) >= 0.05)
    return float(ok[sel].mean()), (float((P[sel & ok] > 0.5).mean()) if (sel & ok).any() else float("nan"))


def main():
    out = {}
    pr = np.load(ROOT / "experiments/396-jurkat-tilt-A/work/prior_ctx_challenge.npz", allow_pickle=True)
    P, gpr = pr["P_down"], pr["genes"].astype(str)
    cb = np.load(ROOT / "experiments/528-callprop-block/work/callblock528_abc.npz", allow_pickle=True)
    assert (cb["genes"].astype(str) == gpr).all()
    for c in "ABC":
        f = np.load(ROOT / f"experiments/528-callprop-block/work/feat528_{c}.npz", allow_pickle=True)
        assert (f["genes"].astype(str) == gpr).all()
        keep, cpm, key = f["keep"].astype(bool), f["mean_cpm"], cb[f"key_{c}"]
        idx = np.flatnonzero(keep)
        tc = set(idx[np.argsort(-cpm[idx], kind="stable")[:2200]].tolist())
        tk = set(idx[np.argsort(-key[idx], kind="stable")[:2200]].tolist())
        m = lambda s_: np.isin(np.arange(len(gpr)), list(s_))  # noqa: E731
        sin, sout, kept = m(tk - tc), m(tc - tk), m(tk & tc)
        comp = keep & ~sin & ~kept & ~sout
        out[c] = dict(n_swapped_in=int(sin.sum()), swapped_in=cover(P, sin), swapped_out=cover(P, sout)[0],
                      kept_block=cover(P, kept)[0], complement=cover(P, comp)[0])
        print(f"{c}: swapped-in {sin.sum()} covered {out[c]['swapped_in'][0]:.3f} (down {out[c]['swapped_in'][1]:.3f}); "
              f"swapped-out {out[c]['swapped_out']:.3f}; kept block {out[c]['kept_block']:.3f}; complement {out[c]['complement']:.3f}")
    S = L.load_state(0)
    key = np.load(ROOT / "experiments/528-callprop-block/work/keys528.npz")["lcap"]
    ref = np.load(L.REF, allow_pickle=False)
    GF, KEEPH = ref["genes"].astype(str), ref["keep"].astype(bool)
    cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
    names = list(cp["bank_names"].astype(str))
    Pb = cp["P_banks"].astype(float)[[names.index(b) for b in ("k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2")]].mean(0)
    pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
    col = np.array([pos.get(g, -1) for g in GF])
    Pf = np.full(len(GF), np.nan)
    Pf[col >= 0] = Pb[col[col >= 0]]
    Pk = Pf[KEEPH]
    ci, co = [], []
    for i in range(85):
        nb = block_from_key(S["head"][i], S["ref_lin"], 2200, key[i], 0)
        ci.append(cover(Pk, nb & ~S["block"][i] & ~S["head"][i])[0])
        co.append(cover(Pk, S["block"][i] & ~nb)[0])
    out["H1"] = dict(swapped_in_median=float(np.median(ci)), swapped_out_median=float(np.median(co)))
    print(f"H1: swapped-in coverage median {np.median(ci):.3f}; swapped-out {np.median(co):.3f}")
    (HERE / "work/coverage537.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
