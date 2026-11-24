"""300 item 5 -- the MECHANISM row: 272's statistic (the share of eigen-directions whose INDEPENDENT noise estimate
exceeds the eigenvalue, and the share of the WHITENED PREDICTION ENERGY those directions carry) for the incumbent
and for each named arm, on one fold.  The incumbent must reproduce R272 section 0.10 (76.8 % of directions at
s = 0 carrying 86.8 % of the whitened energy on `orion_hct116`); 299's `sel:q0.5` read 99.9 %.

    PYTHONDONTWRITEBYTECODE=1 uv run python mech300.py orion_hct116 tilt:d2.5 spread:d1.0
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
import run300 as RU  # noqa: E402
import op272 as O  # noqa: E402

plolo, L = RU.plolo, RU.L
ALPHA, GAMMA, LAM = RU.ALPHA, RU.GAMMA, RU.LAM
P = RU.P


def report(tag: str, evals, evecs, n_i, s, e_raw, e_wht) -> dict:
    rep = O.spectrum_report(np.asarray(evals), s, n_i, e_raw, e_wht, GAMMA)
    P(f"  {tag:18s} p {rep['p']}  s<0.5 on {100*rep['share_noise_dominated']:.1f} % of directions "
      f"(s=0 {100*rep['share_s_zero']:.1f} %) carrying {100*rep['whitened_energy_share_noise_dominated']:.1f} % of the "
      f"whitened prediction energy (raw share {100*rep['raw_energy_share_noise_dominated']:.1f} %); "
      f"sum n / sum lambda {rep['total_noise_over_total_lambda']:.3f}; energy-weighted mean s {rep['energy_weighted_mean_s']:.4f}")
    return rep


def main(held: str, want: list[str]) -> None:
    t0 = time.time()
    st = RU.stage300(held)
    donors = list(st["sources"])
    banks = [st["banks"][nm] for nm in donors]
    st["fit_banks"] = banks
    wc = st["wcols"]
    targets = st["targets"]
    P(f"[mech {held}] arms {want}, donors {donors}, wcols {wc.size}, {time.time()-t0:.0f}s")

    H = {line: O.load_half_rows(line, targets) for line in donors}
    unit, renorm = O.basis_unit_factor(st["raw_full"])
    D_list = [O.half_differences(H[line])[:, wc] * unit for line in donors]
    pre_wc = np.asarray(st["pre"], dtype=np.float64)[:, wc] * ALPHA
    out: dict = dict(fold=held, arms=want, donors=donors, unit_factor=float(unit), renorm=renorm,
                     n_half_rows={line: int(H[line]["ok"].sum()) for line in donors})

    op = plolo.Operator(L.op_path_for(held))
    reg = plolo.OPERATOR_REGISTRY[f"held_{held}"]
    jm = (ROOT / reg["path"]).with_suffix(".json")
    brs = json.loads(jm.read_text())["basis_rows_by_source"] if jm.exists() else {d: 1.0 for d in donors}
    w_inc = [float(brs[d]) for d in donors]
    n_i, s, _ = O.noise_spectrum(op.evals, op.evecs, D_list, w_inc)
    e_raw, e_wht = O.whitened_energy(pre_wc, op.evecs, op.evals, GAMMA)
    out["incumbent"] = report("incumbent", op.evals, op.evecs, n_i, s, e_raw, e_wht)
    V50 = None
    if any(a.startswith("spread") for a in want):
        top = np.argsort(op.evals)[-RU.N_PC:][::-1]
        V50 = np.ascontiguousarray(np.asarray(op.evecs, dtype=np.float64)[:, top])
    del op, n_i, s, e_raw, e_wht
    gc.collect()

    rel = {nm: RU.load_rel(nm) for nm in donors}
    full_rows = {b.name: np.arange(len(b.perts), dtype=np.int64) for b in banks}
    U = {b.name: RU.tilt_u(b, rel[b.name], full_rows[b.name]) for b in banks}
    V = {}
    if V50 is not None:
        for b in banks:
            V[b.name], _ = RU.spread_v(st, b, full_rows[b.name], wc, V50)
        del V50
        gc.collect()

    def flush() -> None:
        out["seconds"] = round(time.time() - t0, 1)
        out["peak_rss_gib"] = RU.rss()
        (WORK / f"mech300-{held}.json").write_text(json.dumps(out, indent=1, default=float))
        (WORK / f"mech300-{held}.log").write_text("\n".join(RU.LOG) + "\n")

    flush()
    for arm in want:
        kind, dose = arm.split(":")[0], float(arm.split(":")[1][1:])
        rows_by = RU.prune_rows(banks) if kind == "prune" else full_rows
        q = U if kind in ("tilt", "rand") else V
        if kind == "prune":
            q = {b.name: RU.tilt_u(b, rel[b.name], rows_by[b.name]) for b in banks}
        wts = None if dose == 0 else {b.name: RU.dose_weights(q[b.name], dose) for b in banks}
        evals, evecs, fmeta = RU.fit300(banks, rows_by, st["raw_full"], wc, weights_by=wts)
        w_arm = [float(fmeta["basis_rows_by_source"][d]) for d in donors]
        n_i, s, _ = O.noise_spectrum(evals, evecs, D_list, w_arm)
        e_raw, e_wht = O.whitened_energy(pre_wc, evecs, evals, GAMMA)
        out[arm] = report(arm, evals, evecs, n_i, s, e_raw, e_wht)
        out[arm]["fit"] = {k: v for k, v in fmeta.items() if k != "raw_medians"}
        del evals, evecs, n_i, s, e_raw, e_wht
        gc.collect()
        flush()
    P(f"[mech {held}] done {(time.time()-t0)/60:.1f} min, peak {RU.rss():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
