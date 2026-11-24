"""216 — analyse the extended dose sweep: the curve, its turnover, and the ceiling.

The question the mandate poses is arithmetic once the curve is measured.  The board's two
leaders moved their cosine by **x2.06 and x2.07** (215 §3.2), i.e. **+106%**.  This asks
what the CD4 dose can deliver at its own maximum, on the same relative scale, against the
shipped f = 0.10.

Also produces:
  * the turnover, located by a local quadratic fit through the three grid points around
    the argmax rather than read off the grid;
  * paired-TARGET bootstrap CIs on every dose's Delta-cosine against the shipped point
    (204's instrument: +-0.20% relative resolution);
  * `pds` in units of 162's LOO sigma (0.001494 raw) at every dose -- falsifier 2;
  * the target-shuffled control at every dose it was run -- falsifier 1;
  * per-fold curves including leak-free `kolf2` -- falsifier 3.

    uv run python analyse216.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
WORK = HERE / "work"
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
SPAN_PDS = 0.452
SIGMA_PDS = 0.001494          # 162's leave-one-out sigma on raw PLOLO pds
GATE_1_SIGMA = 1.645          # 162 Gate 1
SHIPPED_F = 0.10
BOOT = 6000
SEED = 20260904


def load():
    rows = json.loads((WORK / "sweep216.json").read_text())
    by = {}
    for r in rows:
        by[(r["fold"], r["arm"])] = r
    return rows, by


def macro(by, arm, key):
    return float(np.mean([by[(h, arm)][key] for h in FOLDS]))


def turnover(fs, ys):
    """Locate the maximum by a quadratic through the argmax and its two neighbours."""
    i = int(np.argmax(ys))
    if i == 0 or i == len(ys) - 1:
        return float(fs[i]), float(ys[i]), "grid endpoint"
    x = np.log(np.asarray(fs[i - 1:i + 2], dtype=float))
    y = np.asarray(ys[i - 1:i + 2], dtype=float)
    a, b, c = np.polyfit(x, y, 2)
    if a >= 0:
        return float(fs[i]), float(ys[i]), "grid (no interior fit)"
    xs = -b / (2 * a)
    return float(np.exp(xs)), float(a * xs * xs + b * xs + c), "quadratic in log f"


def main() -> None:
    rows, by = load()
    doses = sorted({r["f"] for r in rows if r["kind"] == "cd4"})
    shuf_doses = sorted({r["f"] for r in rows if r["kind"] == "shuffled"})
    rep: dict[str, object] = {}

    # ---------------------------------------------------------------- the curve
    rng = np.random.default_rng(SEED)
    base_arm = "f%g" % SHIPPED_F
    off_arm = "f0"
    n_t = {h: len(by[(h, base_arm)]["_cos_pt"]) for h in FOLDS}
    boot_idx = {h: rng.integers(0, n_t[h], size=(BOOT, n_t[h])) for h in FOLDS}

    def boot_delta(arm_a, arm_b, key="_cos_pt"):
        """Paired-target bootstrap of the macro (across folds) mean difference."""
        acc = np.zeros(BOOT)
        for h in FOLDS:
            a = np.asarray(by[(h, arm_a)][key])
            b = np.asarray(by[(h, arm_b)][key])
            d = a - b
            acc += d[boot_idx[h]].mean(axis=1)
        acc /= len(FOLDS)
        return float(np.percentile(acc, 2.5)), float(np.percentile(acc, 97.5)), \
            float((acc > 0).mean())

    curve = []
    for f in doses:
        arm = "f%g" % f
        row = dict(f=f,
                   fraction=float(np.mean([by[(h, arm)]["fraction"] for h in FOLDS])),
                   cos_agg=macro(by, arm, "cos_agg"),
                   cos_a15=macro(by, arm, "a15_cos"),
                   cos_a23=macro(by, arm, "a23_cos"),
                   pds=macro(by, arm, "pds"),
                   pds_a15=macro(by, arm, "a15_pds"),
                   q_w=macro(by, arm, "a15_sgn_w"),
                   rr=macro(by, arm, "a23_rr"))
        row["d_pds_vs_shipped"] = row["pds"] - macro(by, base_arm, "pds")
        row["pds_sigma_vs_shipped"] = row["d_pds_vs_shipped"] / SIGMA_PDS
        row["d_pds_vs_off"] = row["pds"] - macro(by, off_arm, "pds")
        for h in FOLDS:
            row[f"cos_agg_{h}"] = by[(h, arm)]["cos_agg"]
            row[f"pds_{h}"] = by[(h, arm)]["pds"]
        if f > 0:
            lo, hi, p = boot_delta(arm, base_arm)
            row["dcos_vs_shipped_ci"] = [lo, hi]
            row["dcos_vs_shipped_P>0"] = p
        curve.append(row)
    rep["curve"] = curve

    c0 = [c for c in curve if c["f"] == 0][0]
    cs = [c for c in curve if c["f"] == SHIPPED_F][0]
    fs = [c["f"] for c in curve if c["f"] > 0]

    rep["turnover"] = {}
    for key in ("cos_agg", "cos_a15", "pds", "q_w"):
        ys = [c[key] for c in curve if c["f"] > 0]
        f_star, y_star, how = turnover(fs, ys)
        grid_i = int(np.argmax(ys))
        rep["turnover"][key] = dict(
            f_star=f_star, value=y_star, method=how,
            grid_argmax_f=fs[grid_i], grid_argmax_value=ys[grid_i],
            rel_vs_shipped_pct=100.0 * (ys[grid_i] / cs[key] - 1.0),
            rel_vs_off_pct=100.0 * (ys[grid_i] / c0[key] - 1.0))

    # ------------------------------------------- the arithmetic the mandate asks for
    best_cos = max(c["cos_agg"] for c in curve)
    best_cos15 = max(c["cos_a15"] for c in curve)
    rep["ceiling_vs_leader_doubling"] = dict(
        leader_move_relative=2.06,
        leader_move_pct=106.0,
        cd4_dose_max_cos_agg_rel_to_shipped=best_cos / cs["cos_agg"],
        cd4_dose_max_cos_agg_pct=100.0 * (best_cos / cs["cos_agg"] - 1.0),
        cd4_dose_max_cos_a15_rel_to_shipped=best_cos15 / cs["cos_a15"],
        cd4_dose_max_cos_a15_pct=100.0 * (best_cos15 / cs["cos_a15"] - 1.0),
        cd4_full_on_off_rel=best_cos / c0["cos_agg"],
        fraction_of_leader_move_explained=(
            (best_cos / cs["cos_agg"] - 1.0) / 1.06),
        shipped_cosine_board=0.079249,
        leader_cosine_board=0.13211,
        board_cosine_needed_rel=0.13211 / 0.079249)

    # ---------------------------------------------------------------- falsifier 1
    shuf = []
    for f in shuf_doses:
        real, sh = "f%g" % f, "shuf%g" % f
        shuf.append(dict(
            f=f,
            real_dcos=macro(by, real, "cos_agg") - c0["cos_agg"],
            shuf_dcos=macro(by, sh, "cos_agg") - c0["cos_agg"],
            real_dpds=macro(by, real, "pds") - c0["pds"],
            shuf_dpds=macro(by, sh, "pds") - c0["pds"],
            shuf_dcos_ci=list(boot_delta(sh, off_arm)[:2]),
            gap=macro(by, real, "cos_agg") - macro(by, sh, "cos_agg")))
    rep["falsifier1_shuffled"] = shuf
    rep["falsifier1_verdict"] = dict(
        n_doses=len(shuf),
        shuffled_negative_on_cos_at_all_doses=bool(all(s["shuf_dcos"] < 0 for s in shuf)),
        shuffled_negative_on_pds_at_all_doses=bool(all(s["shuf_dpds"] < 0 for s in shuf)),
        gap_monotone_increasing=bool(all(
            shuf[i]["gap"] <= shuf[i + 1]["gap"] for i in range(len(shuf) - 1))),
        max_gap=max(s["gap"] for s in shuf))

    # ---------------------------------------------------------------- falsifier 2
    guard = []
    for c in curve:
        if c["f"] == 0:
            continue
        guard.append(dict(f=c["f"], d_pds_raw=c["d_pds_vs_shipped"],
                          sigma=c["pds_sigma_vs_shipped"],
                          within_gate1=bool(abs(c["pds_sigma_vs_shipped"]) < GATE_1_SIGMA)))
    rep["falsifier2_pds_guard"] = dict(
        sigma_raw=SIGMA_PDS, gate1=GATE_1_SIGMA, rows=guard,
        largest_f_within_gate1=max((g["f"] for g in guard if g["within_gate1"]),
                                   default=None))

    # ---------------------------------------------------------------- falsifier 3
    per_fold = {}
    for h in FOLDS:
        ys = [by[(h, "f%g" % f)]["cos_agg"] for f in fs]
        ps = [by[(h, "f%g" % f)]["pds"] for f in fs]
        i = int(np.argmax(ys))
        j = int(np.argmax(ps))
        per_fold[h] = dict(
            cos_argmax_f=fs[i], cos_argmax=ys[i],
            cos_at_shipped=by[(h, base_arm)]["cos_agg"],
            cos_rel_pct=100.0 * (ys[i] / by[(h, base_arm)]["cos_agg"] - 1.0),
            pds_argmax_f=fs[j], pds_argmax=ps[j],
            sign_of_gain_at_shipped=float(np.sign(
                by[(h, base_arm)]["cos_agg"] - by[(h, off_arm)]["cos_agg"])))
    rep["falsifier3_per_fold"] = per_fold
    rep["falsifier3_verdict"] = dict(
        all_folds_positive_at_shipped=bool(
            all(v["sign_of_gain_at_shipped"] > 0 for v in per_fold.values())),
        kolf2_leak_free_cos_argmax_f=per_fold["kolf2"]["cos_argmax_f"],
        argmax_f_range=[min(v["cos_argmax_f"] for v in per_fold.values()),
                        max(v["cos_argmax_f"] for v in per_fold.values())])

    (WORK / "analyse216.json").write_text(json.dumps(rep, indent=2))

    # ------------------------------------------------------------------- printout
    print(f"{'f':>7} {'fraction':>9} {'cos_agg':>9} {'rel%':>7} {'cos_a15':>9} "
          f"{'rel%':>7} {'pds':>9} {'d_pds/sig':>10} {'q_w':>8}")
    for c in curve:
        print(f"{c['f']:>7.2f} {c['fraction']:>9.4f} {c['cos_agg']:>9.6f} "
              f"{100*(c['cos_agg']/cs['cos_agg']-1):>7.1f} {c['cos_a15']:>9.6f} "
              f"{100*(c['cos_a15']/cs['cos_a15']-1):>7.1f} {c['pds']:>9.6f} "
              f"{c['pds_sigma_vs_shipped']:>10.2f} {c['q_w']:>8.6f}")
    print()
    print(json.dumps({k: rep[k] for k in ("turnover", "ceiling_vs_leader_doubling",
                                          "falsifier1_verdict", "falsifier3_verdict")},
                     indent=2))
    print("largest f inside Gate-1:",
          rep["falsifier2_pds_guard"]["largest_f_within_gate1"])


if __name__ == "__main__":
    main()
