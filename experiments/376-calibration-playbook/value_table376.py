#!/usr/bin/env python3
"""376 step 3 — per-line value table for October.

For every recipient bank with >= 500 perturbations, donors = the three production
banks minus the recipient if it is one of them (October-faithful: on Oct 22 the
donors are the production set). Runs 372's two-split within-bank held-out sign
protocol via build_calibration_artifact's own functions (no code duplication), plus
the retrieval change at beta=0.35 on 341's saved frame (read-only, no refit; n/a
for lines outside 341's six).

Writes only into this experiment's work/. Serial, one recipient at a time.
"""
from __future__ import annotations
import gc, json, os, resource, sys, time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(ROOT / "experiments/341-known-line-calibration"))
sys.path.insert(0, str(ROOT / "experiments/334-text-conditioned-transfer"))
import build_calibration_artifact as B
import core341 as C341
import core334 as C334

PROD_DONORS = ["data/donors/k562_gwps", "data/donors/orion_hct116",
               "data/donors/orion_hek293t"]
# Fixed order (deterministic seeds); name -> held-out seed.
RECIPIENTS = [
    ("data/donors/k562_gwps", 3760910),
    ("data/donors/orion_hct116", 3760911),
    ("data/donors/orion_hek293t", 3760912),
    ("data/donors/nadig_jurkat", 3760913),
    ("data/donors/rpe1", 3760914),
    ("data/donors/k562_essential", 3760915),
    ("experiments/166-donor-acquisition/donors/nadig_hepg2", 3760916),
]
EXCLUDED = {"experiments/166-donor-acquisition/donors/frangieh_melanoma": 234,
            "data/donors/h1_2025": 254}
PANEL = ROOT / "data/challenge/pert_counts.csv"
AXIS = ROOT / "data/challenge/gene_names.csv"
BETA = 0.35
WORK341 = ROOT / "experiments/341-known-line-calibration/work"
LINES341 = ["k562_gwps", "orion_hct116", "orion_hek293t", "rpe1",
            "nadig_jurkat", "nadig_hepg2"]


def blend(d, p, beta):
    """374's verbatim fixed-norm blend (346's rule)."""
    nd = np.linalg.norm(d, axis=1, keepdims=True)
    np_ = np.linalg.norm(p, axis=1, keepdims=True)
    assert (nd > 0).all() and (np_ > 0).all()
    m = (1 - beta) * d / nd + beta * p / np_
    nm = np.linalg.norm(m, axis=1, keepdims=True)
    return m * (nd / nm)


def mem_ok() -> bool:
    return True  # caller gates on memory_pressure before launching this script


def main() -> None:
    started, clock = datetime.now(timezone.utc).isoformat(), time.monotonic()
    work = HERE / "work"
    work.mkdir(parents=True, exist_ok=True)

    # ---- K3 gate: 341's Jurkat weighted-sign gain vs diagonal from saved preds --
    short = {r.split("/")[-1]: r for r, _ in RECIPIENTS}
    frame = {}
    for ln in LINES341:
        z = np.load(WORK341 / ln / "predictions.npz")
        frame[ln] = {k: np.asarray(z[k], dtype=np.float64)
                     for k in ("truth", "crossline_diag", "same_line_full")}
        assert list(z["targets"]) == json.loads(
            (ROOT / "experiments/341-known-line-calibration/selection.json")
            .read_text())["targets"], f"341 target order drift for {ln}"
    j = "nadig_jurkat"
    sj = C341.sign_endpoints(frame[j]["same_line_full"], frame[j]["truth"])
    dj = C341.sign_endpoints(frame[j]["crossline_diag"], frame[j]["truth"])
    k3 = float(sj["weighted_sign"].mean() - dj["weighted_sign"].mean())
    print(f"K3 jurkat wsign same-crossdiag: {k3:.8f} (record 0.03380218)")
    # Record is quoted to 8 decimals (341 RESULTS / 374 print), so the honest gate
    # is half-ulp of the 8th decimal, not 1e-9.
    k3_ok = abs(k3 - 0.03380218) < 5e-9
    print("K3 GATE:", "PASS" if k3_ok else "FAIL", flush=True)
    if not k3_ok:
        (work / "value_table.json").write_text(json.dumps(
            {"k3": k3, "k3_ok": False}, indent=2))
        sys.exit(1)

    # ---- retrieval change at beta=0.35 on 341's frame (read-only) ---------------
    retr = {}
    for ln in LINES341:
        d0, p0, tt = frame[ln]["crossline_diag"], frame[ln]["same_line_full"], \
            frame[ln]["truth"]
        base = C334.endpoints(d0, tt)["retrieval"]
        cur = C334.endpoints(blend(d0, p0, BETA), tt)["retrieval"]
        iv = C341.interval((cur - base)[None, :])
        retr[ln] = dict(n_targets=int(len(tt)),
                        retr_diag=float(base.mean()),
                        retr_b035=float(cur.mean()),
                        delta=float(iv["mean"]),
                        target_ci95=iv["target_ci95"])

    # ---- per-recipient held-out ------------------------------------------------
    rows = []
    for rdir, seed in RECIPIENTS:
        t0 = time.monotonic()
        recip = (ROOT / rdir).resolve()
        donors = [(ROOT / d).resolve() for d in PROD_DONORS]
        dropped = None
        if recip in donors:
            donors = [d for d in donors if d != recip]
            dropped = rdir
        axes = B.build_axes(recip, donors, PANEL, AXIS)
        n_train, n_cov = len(axes["training"]), len(axes["covered"])
        print(f"== {rdir}: donors={[str(d.relative_to(ROOT)) for d in donors]} "
              f"training={n_train} covered={n_cov} query={len(axes['test'])} "
              f"recip∩panel={len(axes['recip_panel_overlap'])}", flush=True)
        direction, basal, rms = B.assemble(axes, recip, donors, poison=False)
        ho = B.heldout_two_split(axes, direction, basal, rms, len(donors), seed)
        del direction, basal, rms
        gc.collect()
        key = rdir.split("/")[-1]
        rows.append(dict(recipient=rdir, donors=[str(d) for d in donors],
                         dropped_recipient_from_donors=dropped,
                         n_training=n_train, n_covered=n_cov,
                         n_query=len(axes["test"]),
                         n_missing_consensus=len(axes["missing"]),
                         n_recipient_panel_overlap=len(axes["recip_panel_overlap"]),
                         heldout=ho,
                         retrieval_b035_on_341frame=retr.get(key, "n/a"),
                         seconds=time.monotonic() - t0,
                         max_rss_bytes_macos=resource.getrusage(
                             resource.RUSAGE_SELF).ru_maxrss))
        for split, v in ho.items():
            print(f"   {split}: gain {v['gain_same_minus_diag']:+.5f} "
                  f"CI [{v['ci95'][0]:+.5f},{v['ci95'][1]:+.5f}] "
                  f"shuf {v['gain_shuf_minus_diag']:+.5f} "
                  f"frac {v['frac_improved']:.3f} top64 {v['top64_gain']:+.5f}",
                  flush=True)

    out = dict(started_utc=started, finished_utc=datetime.now(timezone.utc).isoformat(),
               seconds=time.monotonic() - clock, beta=BETA, k3=k3, k3_ok=True,
               donor_policy="production three minus recipient-if-member",
               excluded_banks={k: f"{v} perts (<500)" for k, v in EXCLUDED.items()},
               recipients=rows)
    (work / "value_table.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"wrote work/value_table.json in {time.monotonic() - clock:.0f}s", flush=True)


if __name__ == "__main__":
    main()
