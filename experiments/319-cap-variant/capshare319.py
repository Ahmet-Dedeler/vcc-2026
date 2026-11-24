"""319 -- how much of the REAL 300-target panel the cap actually touches.

  PYTHONDONTWRITEBYTECODE=1 uv run python experiments/319-cap-variant/capshare319.py

`_normrel_pre_whiten_effects` weights the THREE production donors (`k562_gwps`, `orion_hct116`,
`orion_hek293t`) per panel target; the CD4 residual, the H1 union and PerturbMatch enter elsewhere
and are not weighted by this rule.  Production builds that model with `gene_norm=False, eb=False`
(`MB:4972-4977`), so `model.s[name]` is exactly `Donor.load(path).decompose()[1]` and the ratio

    r(t, j) = ||s_j[row(t)][dcols_j]|| / median over the panel rows j carries

can be read off the donor banks alone, without the 20 GiB build.  Each donor is loaded, measured
and released in turn (`Donor.load` + `decompose` are the real functions, so the numbers are the
build's own, not a re-implementation).

KNOWN-ANSWER ROW: the same `donor_panel_ratios` is first run on 313's fixture and must reproduce
`work/test319.json`'s per-donor block, which was produced from the real `MultiDonorModel`, exactly.

Read-only outside `experiments/319-cap-variant/`.  No build, no bank write, no submission.
"""
from __future__ import annotations

import json
import resource
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E313 = ROOT / "experiments/313-normrel-parity"

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(E313))
sys.path.insert(0, str(HERE))

import lib319 as L  # noqa: E402
import parity313 as P313  # noqa: E402

EXP = -1.5
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def main() -> None:
    from vcc2026.donor_transfer import Donor
    from build_submission import ChallengeSpec

    KA: list[dict] = []

    def ka(label, got, tol, want=0.0):
        d = abs(float(got) - float(want))
        ok = bool(d <= tol)
        KA.append(dict(label=label, value=float(got), abs_diff=d, tol=tol, ok=ok))
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.6e} (tol {tol:.0e})")
        return ok

    # ------------------------------------------------------------------ KA: the fixture, same code
    ref = json.loads((WORK / "test319.json").read_text())["cap_reach_fixture"]["per_donor"]
    fx = P313.load_fixture()
    fmodel = P313.build_model(fx)
    ftargets = [str(t) for t in fx["targets"]]
    worst = 0.0
    for name in fmodel.names:
        r, med = L.donor_panel_ratios(
            fmodel.s[name], fmodel.dcols[name], fmodel.pert_index[name], ftargets)
        worst = max(worst, abs(med - ref[name]["panel_median_norm"]),
                    abs(float(np.median(list(r.values()))) - ref[name]["median"]),
                    abs(float(sum(1 for v in r.values() if v > 1)) - ref[name]["n_gt1"]),
                    abs(float(max(r.values())) - ref[name]["max"]))
    ka("KA-C  donor_panel_ratios reproduces test319.json's fixture block from the real MultiDonorModel",
       worst, 0.0)
    del fmodel, fx

    # ------------------------------------------------------------------ the real panel
    spec = ChallengeSpec.load(ROOT / "data/challenge")
    panel = [str(t) for t in spec.perts]
    official = [str(g) for g in spec.genes]
    tindex = {g: i for i, g in enumerate(official)}
    P(f"[319] real bundle data/challenge: {len(panel)} panel targets, {len(official)} official genes")

    donor_paths = {
        "k562_gwps": ROOT / "data/donors/k562_gwps",
        "orion_hct116": ROOT / "data/donors/orion_hct116",
        "orion_hek293t": ROOT / "data/donors/orion_hek293t",
    }
    ratios: dict[str, dict[str, float]] = {}
    medians: dict[str, float] = {}
    per_donor: dict[str, dict] = {}
    for name, path in donor_paths.items():
        d = Donor.load(path)
        dcols = np.asarray([j for j, g in enumerate(map(str, d.genes)) if str(g) in tindex],
                           dtype=np.int64)
        pert_index = {str(p): i for i, p in enumerate(d.perts)}
        # production builds the model with gene_norm=False, eb=False, so `model.s[name]` is
        # exactly `decompose()[1]`.  The row-wise form is used to keep the two 2.6 GB Orion banks
        # from needing a 2.6 GB copy beside them; on the small k562 bank both forms are run and
        # must agree bitwise (KA-C2).
        m = d.delta.mean(axis=0)
        r, med = L.donor_panel_ratios_rows(d.delta, m, dcols, pert_index, panel)
        if name == "k562_gwps":
            m_full, s_full = d.decompose()
            r_full, med_full = L.donor_panel_ratios(s_full, dcols, pert_index, panel)
            worst_rows = max([abs(r[t] - r_full[t]) for t in r] + [abs(med - med_full)])
            ka("KA-C2 the row-wise ratio path == the full `decompose()` path on a real donor bank "
               "(k562_gwps, 269 panel rows)", worst_rows, 0.0)
            del m_full, s_full, r_full
        ratios[name] = r
        medians[name] = med
        vals = np.array(sorted(r.values()))
        per_donor[name] = dict(
            n_panel_targets=int(vals.size), panel_median_norm=med,
            median_r=float(np.median(vals)), n_gt1=int((vals > 1).sum()),
            share_gt1=float((vals > 1).mean()) if vals.size else float("nan"),
            min=float(vals.min()), max=float(vals.max()),
            q05=float(np.quantile(vals, 0.05)), q25=float(np.quantile(vals, 0.25)),
            q75=float(np.quantile(vals, 0.75)), q95=float(np.quantile(vals, 0.95)),
            n_genes=int(d.delta.shape[1]), n_shared_genes=int(dcols.size),
            n_perts=int(d.delta.shape[0]))
        P(f"[319]   {name:14s} {vals.size:3d}/{len(panel)} panel targets carried; panel median row "
          f"norm {med:.6f}; r median {np.median(vals):.6f} range [{vals.min():.4f}, {vals.max():.4f}]; "
          f"r>1 on {int((vals > 1).sum())}/{vals.size} = {100 * (vals > 1).mean():.1f} %; "
          f"q05 {np.quantile(vals, 0.05):.4f} q25 {np.quantile(vals, 0.25):.4f} "
          f"q75 {np.quantile(vals, 0.75):.4f} q95 {np.quantile(vals, 0.95):.4f}; "
          f"{dcols.size}/{d.delta.shape[1]} genes on the official axis")
        del d, m
        P(f"[319]     peak RSS so far {peak_gib():.2f} GiB")

    names = list(donor_paths)
    carried = {nm: set(ratios[nm]) for nm in names}
    uni = {nm: 1.0 for nm in names}          # production's default donor weights

    cells = [(t, nm, ratios[nm][t]) for nm in names for t in ratios[nm]]
    touched = [c for c in cells if c[2] > 1.0]
    share_cells = len(touched) / max(len(cells), 1)

    W_cap = L.rule_weights(ratios, uni, panel, EXP, "cap", carried)
    W_plain = L.rule_weights(ratios, uni, panel, EXP, "plain", carried)
    W_uni = L.rule_weights(ratios, uni, panel, 0.0, "plain", carried)

    ok = np.isfinite(W_cap).all(1) & np.isfinite(W_plain).all(1)
    d_cap_vs_uni = np.abs(W_cap - W_uni).max(1)
    d_plain_vs_uni = np.abs(W_plain - W_uni).max(1)
    d_cap_vs_plain = np.abs(W_cap - W_plain).max(1)

    per_target = []
    for i, t in enumerate(panel):
        have = [nm for nm in names if t in carried[nm]]
        c = sum(1 for nm in have if ratios[nm][t] > 1.0)
        per_target.append(dict(target=t, n_donors=len(have), n_capped=c,
                               max_share_move_vs_uniform=float(d_cap_vs_uni[i]),
                               max_share_move_vs_plain=float(d_cap_vs_plain[i])))
    no_cap = [p["target"] for p in per_target if p["n_donors"] >= 1 and p["n_capped"] == 0]
    all_cap = [p["target"] for p in per_target if p["n_donors"] >= 1
               and p["n_capped"] == p["n_donors"]]
    # a target whose carried donors are ALL uncapped keeps uniform weights -> the row is the
    # incumbent's, bit for bit (the norm match is then the identity).
    P("")
    P(f"[319] REAL PANEL cap reach: {len(touched)}/{len(cells)} (target, donor) cells have r > 1 "
      f"= {100 * share_cells:.2f} %.  Each donor's own panel median r is 1.0 BY CONSTRUCTION, so "
      f"about half of every donor's cells sit above it; the cap is not a rare event.")
    P(f"[319]   targets where NO carried donor is capped -> the row reverts to the INCUMBENT "
      f"exactly: {len(no_cap)}/{len(panel)}")
    P(f"[319]   targets where EVERY carried donor is capped -> the cap is a pure rescale of the "
      f"plain rule's ratios, but not of the row: {len(all_cap)}/{len(panel)}")
    hist = {}
    for p in per_target:
        hist[f"{p['n_capped']}/{p['n_donors']}"] = hist.get(f"{p['n_capped']}/{p['n_donors']}", 0) + 1
    P(f"[319]   capped-donors-per-target histogram (capped/carried): "
      f"{dict(sorted(hist.items()))}")
    P(f"[319]   share move vs the INCUMBENT's uniform weights, per target (max over donors): "
      f"cap median {np.median(d_cap_vs_uni[ok]):.5f} mean {np.mean(d_cap_vs_uni[ok]):.5f} "
      f"max {d_cap_vs_uni[ok].max():.5f}  |  plain median {np.median(d_plain_vs_uni[ok]):.5f} "
      f"mean {np.mean(d_plain_vs_uni[ok]):.5f} max {d_plain_vs_uni[ok].max():.5f}")
    P(f"[319]   cap vs plain, per target (max over donors): median "
      f"{np.median(d_cap_vs_plain[ok]):.5f} mean {np.mean(d_cap_vs_plain[ok]):.5f} "
      f"max {d_cap_vs_plain[ok].max():.5f}")
    ratio_move = float(np.median(d_cap_vs_uni[ok]) / max(np.median(d_plain_vs_uni[ok]), 1e-300))
    P(f"[319]   => the cap keeps {100 * ratio_move:.1f} % of the shipped arm's median weight move "
      f"away from uniform.")

    out = dict(
        bundle=str((ROOT / "data/challenge").relative_to(ROOT)),
        exponent=EXP,
        n_panel_targets=len(panel), n_official_genes=len(official),
        known_answers=KA,
        per_donor=per_donor,
        cells=len(cells), touched=len(touched), share_cells_touched=share_cells,
        targets_no_capped_donor=no_cap, targets_all_donors_capped=len(all_cap),
        capped_per_target_histogram=dict(sorted(hist.items())),
        share_move=dict(
            cap_vs_uniform=dict(median=float(np.median(d_cap_vs_uni[ok])),
                                mean=float(np.mean(d_cap_vs_uni[ok])),
                                max=float(d_cap_vs_uni[ok].max())),
            plain_vs_uniform=dict(median=float(np.median(d_plain_vs_uni[ok])),
                                  mean=float(np.mean(d_plain_vs_uni[ok])),
                                  max=float(d_plain_vs_uni[ok].max())),
            cap_vs_plain=dict(median=float(np.median(d_cap_vs_plain[ok])),
                              mean=float(np.mean(d_cap_vs_plain[ok])),
                              max=float(d_cap_vs_plain[ok].max())),
            cap_over_plain_median_ratio=ratio_move),
        per_target=per_target,
        peak_gib=peak_gib(),
    )
    (WORK / "capshare319.json").write_text(json.dumps(out, indent=2))
    (WORK / "capshare319.log").write_text("\n".join(LINES) + "\n")
    P("")
    P(f"[319] KA rows {sum(1 for k in KA if k['ok'])}/{len(KA)} pass; peak {peak_gib():.2f} GiB")
    if not all(k["ok"] for k in KA):
        raise SystemExit("KNOWN-ANSWER FAILURE")


if __name__ == "__main__":
    main()
