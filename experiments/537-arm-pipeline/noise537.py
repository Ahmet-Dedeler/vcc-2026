#!/usr/bin/env python3
"""537 — how noisy is the DE channel on the board for TODAY's construction?

The rule in use (AGENTS 11.1 rule 13, FINDINGS §1): "the ordering channel carries sd ~0.0045 overall between digest-identical
arms". Its source (230 §6) is the spread of board Δreach among arms with IDENTICAL per-gene totals but DIFFERENT emission laws
(a21/a22 resample, rc_hybrid, rc_callset: library CV 0.019 -> 0.382, n_pred 4,420 -> 1,141). That is a spread of designs,
not the re-roll noise of one design. This script measures the re-roll noise of the current design two ways:

1. Board pairs of ours that keep the DE design and change only the bulk (the ballast count K changes, so the carriers are
   re-drawn on other cells: a re-roll plus a small ballast effect on nmae) — read directly from data/board-entries.
2. H1: 435's production construction (446bj2z's DE design) on 3 emission seeds; per-target member rows with 521's exact
   functions (528's price machinery, realized order) at 521's regimes; the board-size re-roll sd of a paired difference
   = sqrt(2) x per-target seed sd / sqrt(900 scored target-rows) (A/B/C: 3 contexts x ~300 targets).
Units: raw members and overall (192's INV / 6: fid 3.32482, jac 2.72081, reach 1.12927, nmae 1.64787).
Writes work/noise537.json. Light (< 1 GB).
"""
from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/528-callprop-block"))
sys.path.insert(0, str(ROOT / "experiments/521-strength-nreal"))
import lib521 as L  # noqa: E402
import price528 as P5  # noqa: E402  (read-only: rows_at, load_arm)
from table521 import regime_pool  # noqa: E402

INV = dict(fid=3.32482, jac=2.72081, reach=1.12927, nmae=1.64787, pds=2.2138, mse=1.051)
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}


def board():
    L_ = json.load(open(ROOT / "submissions/SUBMITTED_LEDGER.json"))
    name = {v.get("entry_id"): Path(v.get("package", "")).name.split("-")[0] for v in L_.values()}
    E = {}
    for f in glob.glob(str(ROOT / "data/board-entries/*.json")):
        d = json.load(open(f))
        e = d.get("entry", d)
        E[name.get(e["entry_id"], e["entry_id"])] = {m: e[k] for m, k in RAW.items()} | {"overall": e["score_avg"]}
    # pairs that keep the DE design (same prior, MDE, repair, zs) and change only the bulk; the ballast re-draws carriers
    pairs = [("410b", "417a", "q 1.25 (bulk)"), ("417a", "417e", "stack (bulk)"), ("433h", "433k", "q 1.5 + re-level (bulk)"),
             ("433k", "446b", "q 1.75 (bulk)"), ("446b", "446bj2", "LODO key (bulk)"), ("446b", "446bd2", "DLD-1 union (bulk)"),
             ("453dr", "457dr", "CD4 x1.5 on the October path (bulk)")]
    out = []
    for a, b, what in pairs:
        if a in E and b in E:
            d = {m: E[b][m] - E[a][m] for m in RAW}
            out.append({"pair": f"{b} - {a}", "what": what, **{m: d[m] for m in RAW},
                        "reach_overall": d["reach"] * INV["reach"] / 6, "fid_overall": d["fid"] * INV["fid"] / 6,
                        "jac_overall": d["jac"] * INV["jac"] / 6})
    return out


def h1():
    res = {}
    states = {s: L.load_state(s) for s in (0, 1, 2)}
    for reg in ("1.0", "0.5", "0.35", "cmp"):
        R = []
        for s in (0, 1, 2):
            A = P5.load_arm("base", s, True)            # 435 with zs = 446bj2z's DE construction on H1
            pool = regime_pool(states[s], reg)[0]
            R.append(P5.rows_at(states[s], A, pool))   # (4, 85): fid, jac, reach, nmae
        R = np.stack(R)                                 # (3 seeds, 4, T)
        sd_t = np.nanstd(R, axis=0, ddof=1)             # per-target seed sd (4, T)
        n_ok = np.isfinite(R).all(0).sum(1)
        mean_sd = np.sqrt(np.nanmean(sd_t ** 2, axis=1))   # rms per-target sd per member
        panel_sd_h1 = np.nanstd(np.nanmean(R, axis=2), axis=0, ddof=1)   # sd of the 85-target panel mean across seeds
        board_diff = np.sqrt(2) * mean_sd / np.sqrt(900.0)
        res[reg] = {m: dict(per_target_sd=float(mean_sd[i]), h1_panel_sd=float(panel_sd_h1[i]), n_targets=int(n_ok[i]),
                            board_diff_sd_raw=float(board_diff[i]), board_diff_sd_overall=float(board_diff[i] * INV[m] / 6))
                    for i, m in enumerate(("fid", "jac", "reach", "nmae"))}
        tot = float(np.sqrt(sum(res[reg][m]["board_diff_sd_overall"] ** 2 for m in ("fid", "jac", "reach", "nmae"))))
        res[reg]["all_four_overall"] = tot
    return res


def main():
    b = board()
    print("Board pairs that keep the DE design (bulk-only; carriers re-drawn through the ballast):")
    print(f"{'pair':18s} {'what':38s} {'dfid':>9s} {'djac':>9s} {'dreach':>9s} {'dnmae':>9s} | reach overall")
    for r in b:
        print(f"{r['pair']:18s} {r['what']:38s} {r['fid']:+.6f} {r['jac']:+.6f} {r['reach']:+.6f} {r['nmae']:+.6f} | {r['reach_overall']:+.6f}")
    rr = np.array([r["reach"] for r in b])
    print(f"reach |d| raw: median {np.median(np.abs(rr)):.5f}, rms {np.sqrt(np.mean(rr ** 2)):.5f}  -> overall rms "
          f"{np.sqrt(np.mean(rr ** 2)) * INV['reach'] / 6:.5f}")
    for m in ("fid", "jac", "nmae"):
        v = np.array([r[m] for r in b])
        print(f"{m} |d| raw: median {np.median(np.abs(v)):.6f}, rms {np.sqrt(np.mean(v ** 2)):.6f} -> overall rms "
              f"{np.sqrt(np.mean(v ** 2)) * INV[m] / 6:.6f}")
    h = h1()
    print("\nH1 (435 + zs, 3 seeds): per-target seed sd -> board paired-difference sd (900 target-rows)")
    for reg, v in h.items():
        print(f"  f {reg:4s}: " + "  ".join(f"{m} raw {v[m]['board_diff_sd_raw']:.6f} (ov {v[m]['board_diff_sd_overall']:.6f})"
                                           for m in ("fid", "jac", "reach", "nmae")) + f" | all four overall {v['all_four_overall']:.6f}")
    (HERE / "work/noise537.json").write_text(json.dumps({"board_pairs": b, "h1": h}, indent=1))


if __name__ == "__main__":
    main()
