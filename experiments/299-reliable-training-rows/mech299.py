"""299 -- the MECHANISM row: does selecting reliable training rows move the operator's spectrum?

272's statistic, computed identically for the incumbent (cached fold operator) and for the reliability-selected
refit, in one process on one fold: with an INDEPENDENT noise estimate N built from the donors' batch-parity half
differences (`op272.noise_spectrum`), the share of eigen-directions on which n_i exceeds lambda_i (s_i < 0.5, i.e.
the direction is noise-dominated) and the share of the WHITENED PREDICTION ENERGY of the fold's `pre` that those
directions carry.  The incumbent reads 77-79 % of directions / 87-92 % of the energy (R272 section 0.10).

    PYTHONDONTWRITEBYTECODE=1 uv run python mech299.py orion_hct116 sel:q0.5
"""
from __future__ import annotations

import gc
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/272-noise-corrected-operator"))
import run299 as RU  # noqa: E402
import op272 as O  # noqa: E402

plolo, L = RU.plolo, RU.L
ALPHA, GAMMA, LAM = RU.ALPHA, RU.GAMMA, RU.LAM
P = RU.P


def report(tag: str, evals, evecs, n_i, s, e_raw, e_wht) -> dict:
    rep = O.spectrum_report(np.asarray(evals), s, n_i, e_raw, e_wht, GAMMA)
    P(f"  {tag:16s} p {rep['p']}  s<0.5 on {100*rep['share_noise_dominated']:.1f} % of directions "
      f"(s=0 {100*rep['share_s_zero']:.1f} %) carrying {100*rep['whitened_energy_share_noise_dominated']:.1f} % of the "
      f"whitened prediction energy (raw share {100*rep['raw_energy_share_noise_dominated']:.1f} %); "
      f"sum n / sum lambda {rep['total_noise_over_total_lambda']:.3f}; energy-weighted mean s {rep['energy_weighted_mean_s']:.4f}")
    return rep


def main(held: str, arm: str) -> None:
    t0 = time.time()
    st = RU.stage299(held)
    donors = list(st["sources"])
    wc = st["wcols"]
    targets = st["targets"]
    P(f"[mech {held} / {arm}] donors {donors}, wcols {wc.size}, {time.time()-t0:.0f}s")

    H = {line: O.load_half_rows(line, targets) for line in donors}
    unit, renorm = O.basis_unit_factor(st["raw_full"])
    D_list = [O.half_differences(H[line])[:, wc] * unit for line in donors]
    pre_wc = np.asarray(st["pre"], dtype=np.float64)[:, wc] * ALPHA
    out: dict = dict(fold=held, arm=arm, donors=donors, unit_factor=float(unit), renorm=renorm,
                     n_half_rows={line: int(H[line]["ok"].sum()) for line in donors})

    # ---------------- the incumbent (cached fold operator)
    op = plolo.Operator(L.op_path_for(held))
    meta = json.loads(Path(str(op.path)).with_suffix(".json").read_text()) if hasattr(op, "path") else None
    reg = plolo.OPERATOR_REGISTRY[f"held_{held}"]
    jm = (ROOT / reg["path"]).with_suffix(".json")
    brs = json.loads(jm.read_text())["basis_rows_by_source"] if jm.exists() else {d: 1.0 for d in donors}
    w_inc = [float(brs[d]) for d in donors]
    n_i, s, _ = O.noise_spectrum(op.evals, op.evecs, D_list, w_inc)
    e_raw, e_wht = O.whitened_energy(pre_wc, op.evecs, op.evals, GAMMA)
    out["incumbent"] = report("incumbent", op.evals, op.evecs, n_i, s, e_raw, e_wht)
    out["incumbent"]["basis_weights"] = w_inc
    del op, n_i, s, e_raw, e_wht
    gc.collect()

    # ---------------- the selected refit
    rel = {nm: RU.load_rel(nm) for nm in donors}
    banks = [st["banks"][nm] for nm in donors]
    kind, q = arm.split(":")[0], float(arm.split(":")[1][1:])
    seed = 299000 + (["k562_gwps", "orion_hct116", "orion_hek293t"].index(held))
    rows_by, auds = {}, {}
    for b in banks:
        idx, a = RU.select_rows(b, rel[b.name], kind, q, seed)
        rows_by[b.name] = idx
        auds[b.name] = a
    evals, evecs, fmeta = RU.fit299(banks, rows_by, st["raw_full"], wc)
    w_sel = [float(fmeta["basis_rows_by_source"][d]) for d in donors]
    n_i, s, _ = O.noise_spectrum(evals, evecs, D_list, w_sel)
    e_raw, e_wht = O.whitened_energy(pre_wc, evecs, evals, GAMMA)
    out[arm] = report(arm, evals, evecs, n_i, s, e_raw, e_wht)
    out[arm]["basis_weights"] = w_sel
    out["fit"] = fmeta
    out["selection"] = auds
    # the eigenvalue spectra side by side (top decile / tail), and the operator's own conditioning
    out["eigen"] = dict(sel_max=float(evals.max()), sel_min=float(evals.min()),
                        sel_trace=float(evals.sum()), sel_top10_share=float(np.sort(evals)[-int(0.1 * evals.size):].sum() / evals.sum()))
    P(f"  eigen: sel max {out['eigen']['sel_max']:.6f} trace {out['eigen']['sel_trace']:.4f} "
      f"top-decile share {100*out['eigen']['sel_top10_share']:.2f} %")
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = RU.rss()
    (WORK / f"mech299-{held}.json").write_text(json.dumps(out, indent=1, default=float))
    (WORK / f"mech299-{held}.log").write_text("\n".join(RU.LOG) + "\n")
    P(f"[mech {held}] done {(time.time()-t0)/60:.1f} min, peak {RU.rss():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "sel:q0.5")
