"""254 -- the PRODUCTION-panel partial residual the builder staged == 246's `lib246.cd4_partial`
on the same 300 targets, bit for bit; and the staged complete-row object == `plolo.cd4_for`.

Reads `work/stage/partial-composed.npz` (and `partial-partial.npz`, which must be identical), and
`work/stage/cd4-composed.npz`; recomputes both objects with 246's / 162's code on the official
panel order (`data/challenge/pert_counts.csv`, the order the builder uses).  Also re-derives the
partial's solved fraction from the saved residual and the frame's recorded base norm.

    uv run python experiments/254-partial-rows/partialcheck254.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
S = WORK / "stage"
sys.path.insert(0, str(ROOT / "experiments/246-direction-hypotheses"))
import lib246 as L  # noqa: E402
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"


def main() -> None:
    zc = np.load(S / "partial-composed.npz")
    zp = np.load(S / "partial-partial.npz")
    zcd4 = np.load(S / "cd4-composed.npz")
    j = json.loads((S / "stage-composed.json").read_text())
    targets = [str(t) for t in zc["targets"]]
    panel = [line.split(",")[0] for line in (ROOT / "data/challenge/pert_counts.csv").read_text().splitlines()[1:] if line.strip()]
    checks = {}
    checks["staged targets == pert_counts.csv order"] = targets == panel and list(plolo.PANEL) == panel
    checks["partial-composed == partial-partial (the object does not depend on the deflation flag)"] = all(
        np.array_equal(zc[k], zp[k]) for k in ("mean_effect", "gene_votes", "target_condition_count"))
    part246 = L.cd4_partial(tuple(panel), equalise=False)
    checks["staged partial residual == lib246.cd4_partial(production targets) mean_effect, bit for bit"] = np.array_equal(
        zc["mean_effect"], np.asarray(part246.mean_effect))
    checks["... gene_votes and target_condition_count bit for bit"] = (
        np.array_equal(zc["gene_votes"], np.asarray(part246.gene_votes))
        and np.array_equal(zc["target_condition_count"], np.asarray(part246.target_condition_count)))
    info = part246.condition_details[0]
    block = j["audit"]["cd4_partial"]
    checks["row census: 15 partial rows (4 one-condition, 11 two-condition), 15 carrying"] = (
        info["n_partial"] == 15 == block["rows"] == block["carrying_rows"]
        and info["n_one"] == 4 == block["rows_one_condition"] and info["n_two"] == 11 == block["rows_two_conditions"])
    checks["median complete-row norm in the meta == lib246's"] = block["median_complete_row_norm"] == info["median_complete_norm"]
    cd4_162 = plolo.cd4_for(tuple(panel))
    checks["staged complete-row CD4 object == plolo.cd4_for(production targets), bit for bit"] = all(
        np.array_equal(zcd4[k], np.asarray(getattr(cd4_162, k))) for k in ("mean_effect", "gene_votes", "target_condition_count"))
    checks["the partial's support is disjoint from the complete-row object's support"] = not np.any(
        (zc["gene_votes"] > 0) & (zcd4["gene_votes"] > 0))
    # the fraction, re-derived from the saved residual and the recorded base norm
    frame = block["residual_frame"]
    added = np.where(zc["gene_votes"] > 0, zc["mean_effect"], 0.0)
    norms = np.linalg.norm(added, axis=1)
    live = norms[norms > 0]
    per_unit = float(np.median(live)) / frame["whitened_base_median_row_norm"]
    solved = 0.34 / per_unit
    checks["the partial's solved fraction re-derives from the saved residual and the recorded base norm (== to 1e-15)"] = (
        abs(solved - frame["solved_fraction"]) < 1e-15 and float(np.median(norms)) == 0.0
        and frame["shipped_median_rule_would_solve"] == 0.0)
    checks["the whitened base's median row norm is the shipped one (4.229295604454299)"] = (
        frame["whitened_base_median_row_norm"] == 4.229295604454299)
    partial_targets = [panel[i] for i in np.flatnonzero(zc["gene_votes"].sum(1) > 0)]
    out = dict(checks=checks, partial_targets=partial_targets, n_partial=int(info["n_partial"]), n_one=int(info["n_one"]),
               n_two=int(info["n_two"]), median_complete_norm=float(info["median_complete_norm"]),
               carrying_median_norm=float(np.median(live)), solved_fraction=frame["solved_fraction"], solved_rederived=solved,
               per_row_norms={t: float(norms[panel.index(t)]) for t in partial_targets},
               per_row_conditions={t: int(zc["target_condition_count"][panel.index(t)]) for t in partial_targets})
    (WORK / "partialcheck254.json").write_text(json.dumps(out, indent=1))
    for k, v in checks.items():
        print(f"[{'PASS' if v else 'FAIL'}] {k}")
    print(f"partial targets ({len(partial_targets)}): {partial_targets}")
    print(f"carrying-row median norm {np.median(live):.6f}, complete-row median norm {info['median_complete_norm']:.6f}, "
          f"solved fraction {frame['solved_fraction']!r} (re-derived {solved!r})")
    ok = all(checks.values())
    print(f"PARTIALCHECK {'PASS' if ok else 'FAIL'} {sum(checks.values())}/{len(checks)}")
    raise SystemExit(0 if ok else 1)


if __name__ == "__main__":
    main()
