#!/usr/bin/env python3
"""630 — WHERE the UP pass's capacity goes: the binding condition for every gene-round it does not deliver.

Runs src_rank_up.up_pass2 (the builder module, unmodified semantics) on the first N seed-0 H1 blocks after the zs pass,
with each `continue` branch tagged by a distinct counter (the function source is copied and only the counter names are
rewritten; a KA asserts the instrumented copy returns the identical block). Reports per target:
  no_carrier      no carrier has sub-gate budget left
  no_pos_donor    no ballast row raises the mean (never on production blocks)
  promo_budget    the chosen recipient's sub-gate budget < the promotion cost
  promo_donor     donors' in-budget capacity < the promotion cost
  k_zero          nothing left to move
  revert_u        the gene's net U change left its bounds
  revert_cmax     the recipient ended below the control max
and the sub-gate budget: total in carriers, used, and the share of carriers drained to 0.
  uv run python experiments/630-rank-neutral-up/diagcap630.py <variant> [<n_targets>]
Analysis only, ~4 GB (the H1 reference).
"""
from __future__ import annotations

import glob
import inspect
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rnup630 as U  # noqa: E402
import src_rank_up as M  # noqa: E402


def instrumented():
    src = inspect.getsource(M.up_pass2)
    tags = iter(["no_carrier", "no_pos_donor", "promo", "k_zero"])
    out, n = [], 0
    for line in src.splitlines():
        if 'st["n_no_candidate"] += 1' in line:
            tag = next(tags)
            n += 1
            line = line.replace('st["n_no_candidate"] += 1', f'st["n_no_candidate"] += 1; st.setdefault("x_{tag}", 0); '
                                f'st["x_{tag}"] += 1')
        out.append(line)
    src = "\n".join(out)
    assert n == 4, n
    # split the promotion branch into budget vs donor
    src = src.replace('st["x_promo"] += 1',
                      'st.setdefault("x_promo_budget", 0); st.setdefault("x_promo_donor", 0); '
                      'st["x_promo_budget" if k_promo > int(sub_left[a]) else "x_promo_donor"] += 1')
    a = '            st["n_reverted"] += 1\n'
    assert src.count(a) == 1
    src = src.replace(a, a + '            st.setdefault("x_revert_u", 0); st.setdefault("x_revert_cmax", 0)\n'
                      '            st["x_revert_u" if (du < lo_j or du > hi_j) else "x_revert_cmax"] += 1\n')
    a = "    X[:, gidx] = Xk\n"
    assert src.count(a) == 1
    src = src.replace(a, a + "    st['x_sub_init'] = int(X[np.ix_(car, sidx)].sum()); st['x_sub_left'] = int(sub_left.sum()); "
                      "st['x_drained'] = float((sub_left == 0).mean())\n")
    ns = dict(M.__dict__)
    exec(compile(src, "<up_pass2_instrumented>", "exec"), ns)
    return ns["up_pass2"]


def main() -> int:
    variant = sys.argv[1] if len(sys.argv) > 1 else "b_k100_q20r3p10"
    n_t = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    spec = M.parse_spec(variant)
    f2 = instrumented()
    r, c, _ = U.E.load_reference()
    keep = np.asarray(c.keep, dtype=bool)
    kpos = np.cumsum(keep) - 1
    gpos = {str(g): i for i, g in enumerate(r["genes"])}
    lib_med = float(np.median(r["libraries"]))
    tie = M.scorer_tie(c)
    mde = 0.24 / np.sqrt(np.maximum(np.asarray(c.ref_lin, float) * lib_med / 1e6, 1e-12))
    files = sorted(glob.glob(str(HERE / f"work/blocks/{U.TAG.format(seed=0)}_*.npz")))[:n_t]
    zspec = U.ZS.parse_spec("dn_k250_zs")
    keys = ["x_no_carrier", "x_no_pos_donor", "x_promo_budget", "x_promo_donor", "x_k_zero", "x_revert_u", "x_revert_cmax"]
    rows = []
    for f in files:
        X0, n_car, target = U.load_block(f)
        _, _, _, _, l0 = U.E.block_stats(X0, c)
        Xz, _ = U.ZS.postpass(X0, n_car, c, l0, zspec)
        p0, q0, lz = M.de_stats(Xz, c, tie)
        sel = (q0 >= 0.05) & np.isfinite(lz) & (lz > 0)
        full = gpos.get(target, -1)
        if full >= 0 and keep[full]:
            sel[kpos[full]] = False
        t = np.maximum(lz, spec["K"] * mde)
        v3 = {k: spec[k] for k in ("quantum", "rounds", "rec_pool") if k in spec}
        Xa, sa = M.up_pass2(Xz, n_car, c, lz, t, sel, cross=spec["cross"], **v3)
        Xb, sb = f2(Xz, n_car, c, lz, t, sel, cross=spec["cross"], **v3)
        assert np.array_equal(Xa, Xb), "instrumented copy differs from the module"
        lib = Xz.sum(1)
        rows.append({**{k: sb.get(k, 0) for k in keys}, "sel": sb["n_sel"], "moved": sb["n_moved"], "frac": sb["frac_got"],
                     "sub_init": sb["x_sub_init"], "sub_used": sb["x_sub_init"] - sb["x_sub_left"],
                     "drained": sb["x_drained"], "car_lib_med": float(np.median(lib[:n_car])),
                     "sub_share": sb["x_sub_init"] / float(lib[:n_car].sum())})
        rw = rows[-1]
        print(f"{target:10s} sel {rw['sel']:4d} moved {rw['moved']:4d} frac {rw['frac']:.3f} | "
              + " ".join(f"{k[2:]} {rw[k]}" for k in keys)
              + f" | sub-gate in carriers {rw['sub_init']} ({rw['sub_share']:.3f} of carrier counts), used {rw['sub_used']}, "
              f"drained {rw['drained']:.2f}", flush=True)
    m = {k: float(np.mean([r_[k] for r_ in rows])) for k in rows[0]}
    print("\nMEAN " + " ".join(f"{k} {v:.3g}" for k, v in m.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
