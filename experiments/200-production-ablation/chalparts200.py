"""200 — the CHALLENGE panel's effect matrix, decomposed per donor, for the `mse` channel.

`mse` cannot be priced on PLOLO: 162 §3 measured its leave-one-arm-out RMS at 0.112 against a
predict-zero null of 0.044 -- 2.5x WORSE than predicting nothing -- because a cross-line truth
wants an amplitude 17-46x smaller than production ships.  The only instrument that has ever
retrodicted a board `mse` move is 096's decomposition on `data/challenge/context_{A,B,C}`,

    mse - 1 = ( s^2 A - 2 s X + Kc - min(s^2 B, claim) ) / sumD                          (*)

and `A`, `B` and `claim` are computable from our own side on the real panel.  For an EMISSION
arm the effect matrix does not move and 096's own `dh` (read off the shipped v11 group sums)
is the right object -- that is what `mse200.py`'s known-answer gates use.  For a GEOMETRY arm
it does move, so the panel's effect matrix has to be rebuilt with the component removed.

This module builds it through the same E062 path `folds200.py` uses, with all three
production donors and no held-out library, and caches the per-donor pieces.

VALIDATION, and it is not optional: production's own audit numbers are recorded in
`submissions/a15-cd4f010.vcc.meta.json` and this reconstruction has to land on them --

    effects.whitening.before_operator_median_norm      4.122892101713494
    effects.whitening.after_all_median_norm            4.2292956044543
    effects.cd4.residual_frame.f_per_unit_fraction     1.055664536904655
    effects.cd4.fraction (solved at target_f = 0.10)   0.09472706196346517
    main_source_support.per_donor_target_counts        269 / 290 / 294
    effects.cd4.target_condition_count_histogram["3"]  276

`report()` prints all six against the artifact.  E062 zeroes the 300 panel coordinates and
production does not (`before_all` 4.2381 > `before_operator` 4.1229 is exactly that mass), so
the two differ off the operator axis by construction; the check is on the operator axis, where
they must agree, and the residual is reported rather than assumed.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
P162 = os.path.join(ROOT, "experiments/162-retrodicting-validator")
sys.path.insert(0, P162)
sys.path.insert(0, HERE)

import plolo                                                          # noqa: E402
from plolo import (CD4_NAMES, CD4_ROOT, DONOR_PATHS, E062, E064, OFFICIAL,  # noqa: E402
                   PANEL, PANEL_SET, PM_PATH)

PROD_SOURCES = ("k562_gwps", "orion_hct116", "orion_hek293t")
os.makedirs(WORK, exist_ok=True)


def parts(limit_gib: float = 24.0) -> dict:
    """Per-donor effects on the 300-target challenge panel, plus PerturbMatch and CD4."""
    cache = os.path.join(WORK, "parts-challenge.npz")
    pos = {g: i for i, g in enumerate(OFFICIAL)}
    panel_positions = np.asarray([pos[g] for g in PANEL], dtype=np.int64)
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        out = dict(
            targets=tuple(str(t) for t in z["targets"]),
            src={n: z[f"src_{n}"].astype(np.float64) for n in PROD_SOURCES},
            avail={n: z[f"avail_{n}"] for n in PROD_SOURCES},
            mean_on_axis={n: z[f"mean_{n}"] for n in PROD_SOURCES},
            pm_effects=z["pm_effects"].astype(np.float64), pm_avail=z["pm_avail"],
            shared=z["shared"], wcols=z["wcols"], panel_positions=panel_positions,
            audit=json.loads(str(z["audit"])))
        print("[challenge] parts cache hit", flush=True)
        return out

    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in PROD_SOURCES}
    prepared = {n: E062.prepare_main_source(b, PANEL, OFFICIAL, panel_positions,
                                            chunk_rows=256, max_rss_gib=limit_gib)
                for n, b in banks.items()}
    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm_prep = E062.prepare_main_source(pm_bank, PANEL, OFFICIAL, panel_positions,
                                       chunk_rows=256, max_rss_gib=limit_gib)
    shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET)
    wcols = E064.operator_columns([banks[n] for n in PROD_SOURCES], OFFICIAL, PANEL_SET)

    audit = dict(targets=len(PANEL), wcols=int(len(wcols)),
                 shared_columns=int(np.asarray(shared).size),
                 per_donor_target_counts={n: int(prepared[n].target_available.sum())
                                          for n in PROD_SOURCES},
                 pm_available=int(pm_prep.target_available.sum()),
                 raw_median_norm={n: float(prepared[n].raw_median_norm)
                                  for n in PROD_SOURCES})
    np.savez_compressed(
        cache, targets=np.asarray(PANEL), audit=json.dumps(audit), wcols=wcols,
        pm_effects=pm_prep.effects.astype(np.float32), pm_avail=pm_prep.target_available,
        shared=np.asarray(shared),
        **{f"src_{n}": prepared[n].effects.astype(np.float32) for n in PROD_SOURCES},
        **{f"avail_{n}": prepared[n].target_available for n in PROD_SOURCES},
        **{f"mean_{n}": prepared[n].mean_on_axis for n in PROD_SOURCES})
    print(f"[challenge] parts built: {json.dumps(audit)}", flush=True)
    return dict(targets=tuple(PANEL), src={n: prepared[n].effects for n in PROD_SOURCES},
                avail={n: prepared[n].target_available for n in PROD_SOURCES},
                mean_on_axis={n: prepared[n].mean_on_axis for n in PROD_SOURCES},
                pm_effects=pm_prep.effects, pm_avail=pm_prep.target_available,
                shared=np.asarray(shared), wcols=wcols,
                panel_positions=panel_positions, audit=audit)


def cd4() -> object:
    """CD4 residual prepared against the challenge panel (production's `complete mean3`)."""
    cache = os.path.join(WORK, "cd4-challenge.npz")
    if os.path.exists(cache):
        z = np.load(cache)
        return E062.CD4Prepared(mode="gene_norm", mean_effect=z["mean_effect"],
                                gene_votes=z["gene_votes"],
                                target_condition_count=z["tcc"], condition_details=[])
    banks = [E062.Bank.load(n, CD4_ROOT / n, require_cd4_scale=True) for n in CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, PANEL, OFFICIAL, PANEL_SET, "gene_norm",
                                 chunk_rows=256, max_rss_gib=24.0)
    complete = np.asarray(prep.target_condition_count) == len(CD4_NAMES)
    me = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    gv = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    me[~complete] = 0.0
    gv[~complete] = 0
    tcc = np.asarray(prep.target_condition_count, dtype=np.uint8).copy()
    np.savez_compressed(cache, mean_effect=me, gene_votes=gv, tcc=tcc)
    return E062.CD4Prepared(mode="gene_norm", mean_effect=me, gene_votes=gv,
                            target_condition_count=tcc, condition_details=[])


ARTIFACT = dict(before_operator_median_norm=4.122892101713494,
                after_all_median_norm=4.2292956044543,
                f_per_unit_fraction=1.055664536904655,
                cd4_solved_fraction=0.09472706196346517,
                per_donor_target_counts={"k562_gwps": 269, "orion_hct116": 290,
                                         "orion_hek293t": 294},
                cd4_complete_targets=276)


def report() -> dict:
    """The six known-answer rows against `submissions/a15-cd4f010.vcc.meta.json`."""
    from folds200 import consensus
    from plolo import PROD_OPS, Operator, solve_fraction
    p = parts()
    c = cd4()
    pre = consensus(p, pm_blend=0.50, weights=None) * 0.982
    op = Operator(PROD_OPS / "whitener_union_pruned.npz")
    before_op = float(np.median(np.linalg.norm(pre[:, p["wcols"]], axis=1)))
    whitened = op.apply(pre, p["wcols"], 0.65, 0.0)
    after_all = float(np.median(np.linalg.norm(whitened, axis=1)))
    frac, per_unit = solve_fraction(whitened, c, 0.10)
    del op
    got = dict(before_operator_median_norm=before_op, after_all_median_norm=after_all,
               f_per_unit_fraction=per_unit, cd4_solved_fraction=frac,
               per_donor_target_counts={n: int(np.asarray(p["avail"][n]).sum())
                                        for n in PROD_SOURCES},
               cd4_complete_targets=int((np.asarray(c.target_condition_count) == 3).sum()))
    print(f"\n{'=' * 96}\nKNOWN-ANSWER: reconstruction vs the SHIPPED artifact's own meta")
    print(f"{'quantity':<34}{'artifact':>22}{'here':>22}{'rel':>10}")
    for k in ("before_operator_median_norm", "after_all_median_norm",
              "f_per_unit_fraction", "cd4_solved_fraction"):
        a, b = ARTIFACT[k], got[k]
        print(f"{k:<34}{a:>22.12f}{b:>22.12f}{(b / a - 1) * 100:>9.4f}%")
    print(f"{'per_donor_target_counts':<34}{str(ARTIFACT['per_donor_target_counts']):>44}")
    print(f"{'':<34}{str(got['per_donor_target_counts']):>44}")
    print(f"{'cd4_complete_targets':<34}{ARTIFACT['cd4_complete_targets']:>22}"
          f"{got['cd4_complete_targets']:>22}")
    print(f"{'operator columns':<34}{17942:>22}{int(len(p['wcols'])):>22}")
    with open(os.path.join(WORK, "chalparts200.json"), "w") as fh:
        json.dump(dict(artifact=ARTIFACT, got=got, audit=p["audit"]), fh, indent=1,
                  default=float)
    return got


if __name__ == "__main__":
    report()
