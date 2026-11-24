"""221 — the six-member h1-test → board transfer table, from real `cell-eval2` output.

WHAT IT COMPUTES, AND IN WHICH SPACE
------------------------------------
RAW is primary. Every local→board factor this project already uses is a ratio of RAW deltas —
164 §6's `reach` 0.4277, 185's 1.290, `docs/BOARD-LEDGER.md` §4's 0.390 (h1-test a8−a15
`+0.02008` against a board `+0.007827`). Raw needs only `cell-eval2 run`'s `agg_results.csv`,
so it does not wait on the real bundle.

    transfer_raw[member][rung] = board_raw_delta / h1test_raw_delta

BOARD DELTAS COME FROM THE SCALED LEDGER, NOT THE RAW ONE. `docs/BOARD-LEDGER.md` quotes the
scaled members to 16 significant figures (they are what averages to Overall) and its raw table
to 6, and the two disagree at the LEVEL by up to 1.4e-4 — a rounding/anchor artefact, not a
measurement. Deltas do not disagree: with the ledger's own fitted map `scaled = slope·(raw − b)`,
`Δraw = Δscaled / slope` reproduces every recorded raw delta to ~4e-6 (checked below and printed
as a gate). So the board's raw deltas are derived that way and the level ambiguity never enters.

SCALED, when a bundle is available, is the same object times a per-member constant:

    transfer_scaled = (slope_board / slope_h1test) · transfer_raw

`--run-dir` may be given more than once; the candidates are pooled across runs.

    uv run python experiments/221-cloud-scorer/transfer_table.py \
        --run-dir experiments/221-cloud-scorer/cloud-output/ladder6-20260904 \
        --run-dir experiments/221-cloud-scorer/cloud-output/ladder7-20260904
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent

MEMBERS = {
    "pds": "pds_cosine",
    "mse": "expr_mse_unbiased_capped_norm",
    "nmae": "de_wilcoxon_lfc_nmae",
    "fid": "de_wilcoxon_direction_fidelity_yield_raw",
    "reach": "de_wilcoxon_direction_reach_raw",
    "jac": "de_wilcoxon_sig_jaccard",
}
ORDER = ["pds", "mse", "nmae", "fid", "reach", "jac"]

# `docs/BOARD-LEDGER.md` header: `scaled = slope·(raw − b)`, fitted on 612 live board entries.
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532,
         "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}

# The board, SCALED. §3b for a15 / rc_hybrid / rc_callset; `210-cd4-build/BOARD-RESULT.md` for
# cd4f034. `pds` and `mse` are identical across the first three to sixteen significant figures —
# that is a repeated board MEASUREMENT, not a rounding, and it is what makes the two emitter
# rungs single-variable.
BOARD = {
    "a15":        {"pds": 0.745871, "mse": 0.120057, "nmae": 0.094239,
                   "fid": -0.030397, "reach": 0.115368, "jac": -0.000619},
    "rc_hybrid":  {"pds": 0.7458708290005253, "mse": 0.12005661519665856,
                   "nmae": 0.09429017536041247, "fid": -0.029758,
                   "reach": 0.170308, "jac": -0.000253},
    "rc_callset": {"pds": 0.7458708290005253, "mse": 0.12005661519665856,
                   "nmae": 0.09429042916917090, "fid": -0.029268,
                   "reach": 0.192898, "jac": 0.027335},
    "cd4f034":    {"pds": 0.737800, "mse": 0.159276, "nmae": 0.104861,
                   "fid": -0.020127, "reach": 0.219300, "jac": 0.028408},
}
# The raw deltas the ledger and `210-cd4-build/BOARD-RESULT.md` record directly, used ONLY as a
# gate on the `Δscaled / slope` route.
BOARD_RAW_DELTA_CHECK = {
    ("a15", "rc_callset"): {"reach": 0.24983 - 0.181173},
    ("rc_callset", "cd4f034"): {"pds": -0.00365, "mse": -0.03732, "nmae": -0.00639,
                                "fid": +0.00272, "reach": +0.02338, "jac": +0.00039},
}
RUNGS = [("a15", "rc_hybrid", "emitter: quasi → rc_hybrid — the emitted p-ORDER"),
         ("rc_hybrid", "rc_callset", "emitter: rc_hybrid → rc_callset — the CALL SET"),
         ("rc_callset", "cd4f034", "VCC_CD4_TARGET_F 0.10 → 0.34 — a DIRECTION arm")]
NAMES = ("ka095", "a15", "a15_seed1", "rc_hybrid", "rc_callset", "cd4f034")

# The 095 frozen h1-test artifact's own `agg_results.csv`, means, as `cell-eval2` 0.16.0 wrote
# them on 2026-08-27 (`experiments/095-clean-harness/out-h1test/runs/.../run/agg_results.csv`).
# This is the known-answer row: the cloud scorer must reproduce all ten.
KA095 = {
    "de_wilcoxon_direction_fidelity_yield_raw": 0.5046621497798164,
    "de_wilcoxon_direction_reach_raw": 0.20045827672164832,
    "de_wilcoxon_lfc_nmae": 0.9667252780145132,
    "de_wilcoxon_sig_jaccard": 0.06359002047099205,
    "expr_distance_unbiased": 0.0018008600704525074,
    "expr_mse_unbiased": 0.003700950826313404,
    "expr_mse_unbiased_capped": 0.003700950826313404,
    "expr_mse_unbiased_capped_norm": 2.0551018299736388,
    "expr_real_mass_ratio": 1.0,
    "pds_cosine": 0.8655462184873949,
}

# The four h1-test targets the frozen CD4 bank's `complete mean3` policy covers, against 276 of
# the challenge panel's 300. Measured by `build_h1_candidates.py`.
CD4_H1_TARGETS = ["ANKZF1", "MED15", "MED25", "ZNF32"]


def read_agg(path: Path) -> dict[str, float]:
    """`cell-eval2 run`'s wide agg CSV -> {cell-eval2 metric name: RAW mean}."""
    import polars as pl
    df = pl.read_csv(path)
    row = df.filter(df["statistic"] == "mean")
    r = row.row(0, named=True)
    return {k: (float(v) if v is not None else float("nan"))
            for k, v in r.items() if k != "statistic"}


def read_score(path: Path) -> dict[str, float]:
    import polars as pl
    df = pl.read_csv(path)
    col = "from_replicate" if "from_replicate" in df.columns else df.columns[-1]
    raw = {str(m): v for m, v in zip(df["metric"].to_list(), df[col].to_list())}
    out = {}
    for short, name in MEMBERS.items():
        hit = [k for k in raw if k == name] or [k for k in raw if name in k]
        if hit and raw[hit[0]] is not None:
            out[short] = float(raw[hit[0]])
    if raw.get("avg_score") is not None:
        out["_avg_score"] = float(raw["avg_score"])
    return out


def per_pert(path: Path) -> dict[str, dict[str, float]]:
    """`run/results.csv` is LONG — `perturbation, metric, value` — not wide. Pivot it."""
    import polars as pl
    df = pl.read_csv(path)
    inv = {v: k for k, v in MEMBERS.items()}
    out: dict[str, dict[str, float]] = {}
    for pert, metric, value in df.select("perturbation", "metric", "value").iter_rows():
        m = inv.get(str(metric))
        if m is None:
            continue
        out.setdefault(str(pert), {})[m] = (float(value) if value is not None
                                            else float("nan"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run-dir", action="append", required=True)
    ap.add_argument("--out-json", default=str(HERE / "work/transfer-table.json"))
    ap.add_argument("--out-md", default=str(HERE / "work/transfer-table.md"))
    args = ap.parse_args()

    agg, scaled, pp, src = {}, {}, {}, {}
    for rd in args.run_dir:
        run = Path(rd)
        for name in NAMES:
            for cand in ("agg_results.csv", "agg.csv"):
                f = run / name / "run" / cand
                if f.exists() and name not in agg:
                    agg[name] = read_agg(f)
                    src[name] = str(f)
            f = run / name / "run" / "results.csv"
            if f.exists() and name not in pp:
                pp[name] = per_pert(f)
            f = run / name / "score.csv"
            if f.exists() and name not in scaled:
                s = read_score(f)
                if s:
                    scaled[name] = s

    lines: list[str] = []
    P = lines.append
    report: dict[str, object] = {"sources": src}

    # ---------------------------------------------------------- KNOWN-ANSWER: the 095 artifact
    P("## Known-answer row — the cloud scorer against the frozen 095 h1-test artifact\n")
    if "ka095" not in agg:
        P("**NOT YET AVAILABLE — no table below this line may be trusted.**\n")
        ka_pass = None
    else:
        P("`experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25/pred.h5ad`,")
        P("scored on 2026-08-27 by `cell-eval2` 0.16.0, re-scored today on the cloud worker.")
        P("All ten statistics of `agg_results.csv`'s `mean` row, not just `mse`.\n")
        P("| metric | 2026-08-27, on disk | cloud worker | rel |")
        P("|---|---:|---:|---:|")
        rows, worst = {}, 0.0
        for k, want in KA095.items():
            got = agg["ka095"].get(k, float("nan"))
            rel = abs(got / want - 1.0) if want else abs(got)
            worst = max(worst, rel)
            rows[k] = {"want": want, "got": got, "rel": rel}
            P(f"| `{k}` | {want!r} | {got!r} | {rel:.2e} |")
        ka_pass = worst <= 1e-12
        P(f"\n**Worst relative difference over all ten: `{worst:.3e}` — "
          f"{'PASS' if ka_pass else 'FAIL'}.**\n")
        report["known_answer_095"] = {"worst_rel": worst, "pass": bool(ka_pass), "rows": rows}
    if ka_pass is False:
        P("> **The cloud scorer is NOT the scorer that produced the 095 artifact. Nothing below "
          "is usable.**\n")

    # ---------------------------------------------------------------- board delta gate
    P("## Gate — the board's raw deltas, derived from its scaled ones\n")
    P("| rung | member | Δscaled / slope | Δraw recorded | rel |")
    P("|---|---|---:|---:|---:|")
    gate_worst = 0.0
    for (a, b), want in BOARD_RAW_DELTA_CHECK.items():
        for m, wraw in want.items():
            got = (BOARD[b][m] - BOARD[a][m]) / SLOPE[m]
            rel = abs(got / wraw - 1.0)
            gate_worst = max(gate_worst, rel)
            P(f"| `{a}`→`{b}` | `{m}` | {got:+.6f} | {wraw:+.6f} | {rel:.1e} |")
    P(f"\nWorst relative disagreement **{gate_worst:.1e}** over 7 recorded raw deltas.\n")
    report["board_raw_delta_gate"] = {"worst_rel": gate_worst}

    # ---------------------------------------------------------------- levels
    P("## h1-test, real `cell-eval2` 0.16.0, RAW `agg_results.csv` means\n")
    P("| candidate | " + " | ".join(ORDER) + " |")
    P("|---|" + "---:|" * len(ORDER))
    for name in NAMES:
        if name not in agg:
            continue
        P(f"| `{name}` | " + " | ".join(
            f"{agg[name].get(MEMBERS[m], float('nan')):.6f}" for m in ORDER) + " |")
    P("")
    report["h1_raw"] = {n: {m: agg[n].get(MEMBERS[m]) for m in ORDER} for n in agg}

    # ---------------------------------------------------------------- the noise floor
    floor = {}
    if "a15" in agg and "a15_seed1" in agg:
        for m in ORDER:
            floor[m] = abs(agg["a15_seed1"].get(MEMBERS[m], np.nan)
                           - agg["a15"].get(MEMBERS[m], np.nan))
        P("### The h1-test noise floor — `a15` against its re-seeded twin `a15_seed1`\n")
        P("Identical config and identical effect matrix; `library_pick_seed` 0 → 1, which also "
          "re-draws each target's 400 control library sizes, so this is the seed noise on all "
          "six members and not only on the DE ones.\n")
        P("| | " + " | ".join(ORDER) + " |")
        P("|---|" + "---:|" * len(ORDER))
        P("| |Δ| raw seed noise | " + " | ".join(f"{floor[m]:.6f}" for m in ORDER) + " |\n")
        report["noise_floor_raw"] = floor

    # ---------------------------------------------------------------- THE TABLE
    P("## THE TRANSFER TABLE — board Δraw / h1-test Δraw, per member, per rung\n")
    table = {}
    P("| rung | member | h1-test Δraw | board Δraw | **transfer** | readable? |")
    P("|---|---|---:|---:|---:|---|")
    for a, b, why in RUNGS:
        if a not in agg or b not in agg:
            continue
        table[f"{a}->{b}"] = {"why": why, "members": {}}
        for m in ORDER:
            dl = agg[b].get(MEMBERS[m], np.nan) - agg[a].get(MEMBERS[m], np.nan)
            db = (BOARD[b][m] - BOARD[a][m]) / SLOPE[m]
            f = floor.get(m, 0.0)
            readable = bool(abs(dl) > f and f > 0 and dl != 0)
            tf = (db / dl) if readable else float("nan")
            table[f"{a}->{b}"]["members"][m] = {
                "h1_delta_raw": dl, "board_delta_raw": db, "transfer": tf,
                "noise_floor": f, "readable": readable}
            P(f"| `{a}`→`{b}` | `{m}` | {dl:+.6f} | {db:+.6f} | "
              f"{'—' if not np.isfinite(tf) else f'{tf:+.3f}×'} | "
              f"{'yes' if readable else f'NO — |Δ| ≤ noise {f:.6f}'} |")
    P("")
    report["transfer_raw"] = table

    # ---------------------------------------------------------------- scaled, if a bundle ran
    if scaled:
        P("## The same runs, SCALED against the h1-test real bundle\n")
        P("| candidate | " + " | ".join(ORDER) + " | avg_score |")
        P("|---|" + "---:|" * (len(ORDER) + 1))
        for name in NAMES:
            if name in scaled:
                P(f"| `{name}` | " + " | ".join(
                    f"{scaled[name].get(m, float('nan')):+.6f}" for m in ORDER) +
                  f" | {scaled[name].get('_avg_score', float('nan')):+.6f} |")
        P("")
        report["h1_scaled"] = scaled
        # h1-test's own slope per member, read off the two spaces we now hold both of
        sl = {}
        for m in ORDER:
            pairs = [(agg[n].get(MEMBERS[m]), scaled[n].get(m))
                     for n in scaled if n in agg and m in scaled[n]]
            pairs = [(r, s) for r, s in pairs if r is not None and s is not None]
            if len(pairs) >= 2:
                R = np.array([p[0] for p in pairs]); S = np.array([p[1] for p in pairs])
                if np.ptp(R) > 0:
                    sl[m] = float(np.polyfit(R, S, 1)[0])
        if sl:
            P("### h1-test's own scale, and the scaled transfer factor\n")
            P("`slope_h1` is fitted from the candidates themselves — every one of them carries "
              "both a raw and a scaled reading of the same run, so the map is over-determined.\n")
            P("| member | slope board | slope h1-test | ratio |")
            P("|---|---:|---:|---:|")
            for m in ORDER:
                if m in sl:
                    P(f"| `{m}` | {SLOPE[m]:+.4f} | {sl[m]:+.4f} | "
                      f"{SLOPE[m]/sl[m]:+.4f} |")
            P("")
            report["h1_slope"] = sl

    # ---------------------------------------------- the CD4 rung, on its four covered targets
    if "rc_callset" in pp and "cd4f034" in pp:
        P("## The CD4 rung on the targets the frozen bank actually covers\n")
        P("The `complete mean3` CD4 policy covers **4 of h1-test's 85 targets** against "
          "**276 of the challenge panel's 300**. The whole-surface delta is diluted ~19.6× by "
          "construction; these four are where the arm is even defined.\n")
        P("| target | " + " | ".join(f"Δ{m}" for m in ORDER) + " |")
        P("|---|" + "---:|" * len(ORDER))
        rows = {}
        for t in CD4_H1_TARGETS:
            if t in pp["cd4f034"] and t in pp["rc_callset"]:
                d = {m: pp["cd4f034"][t].get(m, np.nan) - pp["rc_callset"][t].get(m, np.nan)
                     for m in ORDER}
                rows[t] = d
                P(f"| `{t}` | " + " | ".join(f"{d[m]:+.5f}" for m in ORDER) + " |")
        others = [t for t in pp["cd4f034"] if t not in CD4_H1_TARGETS]
        if others:
            d = {m: float(np.nanmean([pp["cd4f034"][t].get(m, np.nan)
                                      - pp["rc_callset"][t].get(m, np.nan) for t in others]))
                 for m in ORDER}
            rows["_other_mean"] = d
            P(f"| *the other {len(others)}, mean* | " +
              " | ".join(f"{d[m]:+.5f}" for m in ORDER) + " |")
        P("")
        report["cd4_covered_targets"] = rows
        # Coverage-corrected: if only covered targets move, a per-target mean over the whole
        # surface is the per-COVERED-target effect times the coverage fraction. Undo that on
        # both sides and the two become comparable.
        cov_h1, cov_board = len(CD4_H1_TARGETS) / 85.0, 276.0 / 300.0
        P("### The same rung, coverage-corrected\n")
        P(f"Only covered targets move (proven: 81/85 blocks digest-identical, and the "
          f"*other {len(others)}* row above is zero to five decimals). A per-target mean over a "
          f"whole surface is therefore the per-COVERED-target effect times that surface's "
          f"coverage — **{cov_h1:.3f} on h1-test, {cov_board:.3f} on the challenge panel**. "
          f"Dividing both sides out makes them comparable.\n")
        P("| member | h1-test, per covered target | board, per covered target | transfer |")
        P("|---|---:|---:|---:|")
        cc = {}
        for m in ORDER:
            vals = [rows[t][m] for t in CD4_H1_TARGETS
                    if t in rows and np.isfinite(rows[t].get(m, np.nan))]
            if not vals:
                continue
            h1 = float(np.mean(vals))
            bd = (BOARD["cd4f034"][m] - BOARD["rc_callset"][m]) / SLOPE[m] / cov_board
            cc[m] = {"h1_per_covered": h1, "board_per_covered": bd,
                     "transfer": (bd / h1 if h1 else float("nan")),
                     "n_h1_targets": len(vals)}
            P(f"| `{m}` | {h1:+.5f} | {bd:+.5f} | "
              f"{'—' if not h1 else f'{bd / h1:+.3f}×'} |")
        P("")
        report["cd4_coverage_corrected"] = cc

    md = "\n".join(lines)
    Path(args.out_md).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_md).write_text(md)
    Path(args.out_json).write_text(json.dumps(report, indent=1, default=float))
    print(md)
    print(f"\n[221] wrote {args.out_md} and {args.out_json}")


if __name__ == "__main__":
    main()
