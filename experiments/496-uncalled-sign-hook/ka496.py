"""496: known answer and reading.
KA: the hook emitter without --mde-add-sign == 435's production-emitter coarse run (seed 0), stats bit-equal.
Reading: each --mde-add-sign run vs 435's coarse base at the same seed. Raw members from the metrics files (pds, mse must
be bit-identical: the hook moves the A-channel only); per-target DE members from the stats via 490's score_de (the
harness's own tables, asserted equal to the metrics' raw), overall-equiv with 435's board slopes, target bootstrap
(2,000 resamples) per seed and pooled over seeds. Diagnostics: unsuppressible calls, clip fraction, fd, flips.
usage: uv run python experiments/496-uncalled-sign-hook/ka496.py
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
W = ROOT / "experiments/496-uncalled-sign-hook/work/h1"; W435 = ROOT / "experiments/435-zs-on-coarse/work/h1"
sys.path.insert(0, str(ROOT / "experiments/490-de-emitter-review"))
import score_de  # noqa: E402

MEM = ["pds", "mse", "nmae", "fid", "reach", "jac"]
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
DE = ["nmae", "fid", "reach", "jac"]
BASE = "base_s{s}_bps_prodem_mde0.1cap0.1min0.3coarse"
ARM = "base_s{s}_bps_hookem_mde0.1cap0.1min0.3coarseas{tau}"


def eq(a, b):
    A, B = np.load(a), np.load(b)
    bad = [k for k in sorted(set(A.files) & set(B.files)) if not np.array_equal(A[k], B[k])]
    return (not bad) and set(A.files) == set(B.files), bad


ka = W / "stats-base_s0_bps_hookem_mde0.1cap0.1min0.3coarse.npz"
if ka.exists():
    e, bad = eq(ka, W435 / f"stats-{BASE.format(s=0)}.npz")
    print("KA hook emitter, no sign source == 435 production coarse (seed 0):", "PASS (bit-equal)" if e else f"FAIL {bad}")
    if not e:
        raise SystemExit(1)
else:
    print("KA PENDING")
zr = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
keep = np.asarray(zr["keep"], dtype=bool); T = zr["targets"].astype(str); G = np.asarray(zr["genes"])[keep].astype(str)


def per_target(tag, wdir):
    st = np.load(wdir / f"stats-{tag}.npz", allow_pickle=False)
    raw, per = score_de.score(T, G, st["p"], st["q"], st["lfc"])
    want = json.loads((wdir / f"metrics-{tag}.json").read_text())["raw"]
    assert max(abs(raw[m] - want[m]) for m in raw) <= 1e-9, f"{tag}: score_de != the harness metrics"
    return raw, per, want


rng = np.random.default_rng(496)
out = {}
for tau in ("0.05", "0"):
    pooled = []
    for s in (0, 1, 2):
        a_tag, b_tag = ARM.format(s=s, tau=tau), BASE.format(s=s)
        if not (W / f"metrics-{a_tag}.json").exists():
            continue
        ra, pa, wa = per_target(a_tag, W)
        rb, pb, wb = per_target(b_tag, W435)
        d = {m: wa[m] - wb[m] for m in MEM}
        falsifier = "OK" if abs(d["pds"]) < 1e-12 and abs(d["mse"]) < 1e-12 else "FIRED (pds/mse moved)"
        ov = sum(SLOPE[m] * d[m] for m in MEM) / 6
        D = np.array([[pa[m].get(t, np.nan) - pb[m].get(t, np.nan) for t in T] for m in DE])      # 4 x targets
        pooled.append(D)
        def ov_of(ii, D=D):
            return float(sum(SLOPE[m] * np.nanmean(D[k][ii]) for k, m in enumerate(DE)) / 6)
        bo = np.array([ov_of(rng.integers(0, len(T), len(T))) for _ in range(2000)])
        ma, mb = json.loads((W / f"metrics-{a_tag}.json").read_text()), json.loads((W435 / f"metrics-{b_tag}.json").read_text())
        dg = json.loads((W / f"diag-{a_tag}.json").read_text())
        dg = dg if isinstance(dg, list) else dg.get("per_target", dg.get("diags", []))
        flips = [x.get("mde_n_sign_flip") for x in dg if isinstance(x, dict) and x.get("mde_n_sign_flip") is not None]
        row = {"d": {m: round(d[m], 6) for m in MEM}, "overall_equiv": round(ov, 6),
               "ci95": [round(float(np.quantile(bo, .025)), 6), round(float(np.quantile(bo, .975)), 6)],
               "falsifier_pds_mse": falsifier,
               "unsuppressible": [mb.get("unsuppressible_mean"), ma.get("unsuppressible_mean")],
               "clip_frac": [mb.get("clip_frac_mean"), ma.get("clip_frac_mean")],
               "fd": [mb.get("mde_fd_median"), ma.get("mde_fd_median")], "fu": [mb.get("mde_fu_median"), ma.get("mde_fu_median")],
               "flips_median": float(np.median(flips)) if flips else None}
        out[f"tau{tau}_s{s}"] = row
        print(f"tau {tau} seed {s}: " + " ".join(f"{m} {d[m]:+.5f}" for m in MEM) + f" | overall {ov:+.5f} {row['ci95']}"
              f" | pds/mse {falsifier} | unsupp {row['unsuppressible']} clip {row['clip_frac']} fd {row['fd']} flips {row['flips_median']}")
    if len(pooled) > 1:
        P = np.nanmean(np.stack(pooled), 0)
        def ovp(ii):
            return float(sum(SLOPE[m] * np.nanmean(P[k][ii]) for k, m in enumerate(DE)) / 6)
        bo = np.array([ovp(rng.integers(0, len(T), len(T))) for _ in range(4000)])
        pt = ovp(np.arange(len(T)))
        out[f"tau{tau}_pooled"] = {"overall_equiv": round(pt, 6), "ci95": [round(float(np.quantile(bo, .025)), 6), round(float(np.quantile(bo, .975)), 6)],
                                   "p_le_0": round(float((bo <= 0).mean()), 4), "seeds": len(pooled),
                                   **{m: round(float(np.nanmean(P[k])), 6) for k, m in enumerate(DE)}}
        print(f"tau {tau} pooled over {len(pooled)} seeds (target bootstrap): overall {pt:+.5f} {out[f'tau{tau}_pooled']['ci95']} "
              f"P(<=0) {out[f'tau{tau}_pooled']['p_le_0']} | " + " ".join(f"{m} {out[f'tau{tau}_pooled'][m]:+.5f}" for m in DE))
(W.parent / "ka496.json").write_text(json.dumps(out, indent=1))
