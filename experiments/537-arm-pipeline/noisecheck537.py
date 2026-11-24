#!/usr/bin/env python3
"""537 — the checker preregistered in PREREGISTRATION-noise.md (written 2026-09-30 13:45Z, before any 10-01 entry).

usage: uv run python experiments/537-arm-pipeline/noisecheck537.py <label>=<arm id>:<ref id> [...]
  e.g. 460dr=<460dr entry>:q3VuTImh4FRGKXH2waCp  479p2g=<id>:YvVru3DAF89Jiyr34olR ...
Reads only data/board-entries/<id>.json. Bands (raw): |d fid| <= 0.0001, |d jac| <= 0.00005, |d reach| <= 0.0040; nmae is
reported. Verdict: CONFIRMED if >= 5 scored pairs (or all but one when fewer than 6) are inside all three bands AND the rms
of the new reach deltas <= 0.0026; REFUTED if that rms >= 0.0050 or >= 3 pairs break the fid or jac band; else UNRESOLVED.
Also prints the pooled reach rms with B1's seven pairs (RESULTS §2.2) and its chi-square 95 % interval.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[2]
RAW = {"pds": "pds_cosine", "mse": "expr_mse_unbiased_capped_norm", "nmae": "de_wilcoxon_lfc_nmae",
       "fid": "de_wilcoxon_direction_fidelity_yield_raw", "reach": "de_wilcoxon_direction_reach_raw",
       "jac": "de_wilcoxon_sig_jaccard"}
BAND = {"fid": 0.0001, "jac": 0.00005, "reach": 0.0040}
B1_REACH = [0.001127, -0.000504, 0.000230, 0.000700, -0.001963, -0.002314, 0.000711]   # RESULTS §2.2 (noise537b.json)


def entry(i):
    d = json.loads((ROOT / f"data/board-entries/{i}.json").read_text())
    return d.get("entry", d)


def main():
    rows = []
    for arg in sys.argv[1:]:
        lab, ids = arg.split("=", 1)
        a, b = ids.split(":")
        ea, eb = entry(a), entry(b)
        d = {m: ea[k] - eb[k] for m, k in RAW.items()}
        inside = all(abs(d[m]) <= BAND[m] for m in BAND)
        rows.append((lab, d, inside))
        print(f"{lab:10s} {a} - {b}: fid {d['fid']:+.6f} jac {d['jac']:+.6f} reach {d['reach']:+.6f} nmae {d['nmae']:+.6f} "
              f"(pds {d['pds']:+.6f} mse {d['mse']:+.6f}) -> {'inside' if inside else 'OUTSIDE'} the bands")
    if not rows:
        raise SystemExit("no pairs given")
    n = len(rows)
    n_in = sum(r[2] for r in rows)
    rr = np.array([r[1]["reach"] for r in rows])
    rms = float(np.sqrt(np.mean(rr ** 2)))
    n_fj = sum(1 for r in rows if abs(r[1]["fid"]) > BAND["fid"] or abs(r[1]["jac"]) > BAND["jac"])
    need = 5 if n >= 6 else n - 1
    if rms >= 0.0050 or n_fj >= 3:
        verdict = "REFUTED"
    elif n_in >= need and rms <= 0.0026:
        verdict = "CONFIRMED"
    else:
        verdict = "UNRESOLVED"
    pooled = np.concatenate([np.array(B1_REACH), rr])
    k = pooled.size
    prms = float(np.sqrt(np.mean(pooled ** 2)))
    lo, hi = prms * np.sqrt(k / chi2.ppf(0.975, k)), prms * np.sqrt(k / chi2.ppf(0.025, k))
    print(f"pairs {n}, inside {n_in} (need {need}); new reach rms {rms:.5f}; fid/jac breaks {n_fj} -> {verdict}")
    print(f"pooled with B1 (n {k}): reach rms {prms:.5f} raw [{lo:.5f}, {hi:.5f}] = {prms * 1.12927 / 6:.5f} overall")


if __name__ == "__main__":
    main()
