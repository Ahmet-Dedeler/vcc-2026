"""436 known answers and the six-member reading.
KA1: em436 (production emitter copy + mde_sign + pure diagnostics), option OFF == 435's prodem coarse stats, bit-equal.
KA2: em436, option on with an EMPTY target set (all-zero vectors on every block) == the same, bit-equal.
Then per seed: donor-signed (`_dsall`, or the variant given as argv[1]) vs the same-seed option-off run, all six raw
members on all targets, on the essential-covered targets and on the others (per-target CSV; nmae on the 58 gated
targets: 20 covered / 38 not), plus mde_lambda / prior_lambda per block and the delivered-LFC change.
Writes work/ka436.json (KA) and work/tab436-<variant>.json."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import polars as pl

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
W = HERE / "work/h1"; W435 = ROOT / "experiments/435-zs-on-coarse/work/h1"
BASE = "mde0.1cap0.1min0.3coarse"
VAR = sys.argv[1] if len(sys.argv) > 1 else "dsall"
MEM = ["pds", "mse", "nmae", "fid", "reach", "jac"]
SLOPE = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532, "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}   # 435's board slopes
PT = {"pds": "pds_cosine", "nmae": "de_wilcoxon_lfc_nmae", "fid": "de_wilcoxon_direction_fidelity_yield_raw",
      "reach": "de_wilcoxon_direction_reach_raw", "jac": "de_wilcoxon_sig_jaccard"}


def eq(a, b):
    A, B = np.load(a), np.load(b)
    bad = [k for k in sorted(set(A.files) & set(B.files)) if not np.array_equal(A[k], B[k])]
    return (not bad) and set(A.files) == set(B.files), bad


ka = {}
for s in (0, 1, 2):
    ref = W435 / f"stats-base_s{s}_bps_prodem_{BASE}.npz"
    for nm, f in ((f"KA1 s{s} off", W / f"stats-base_s{s}_bps_em436_{BASE}.npz"),
                  (f"KA2 s{s} empty", W / f"stats-base_s{s}_bps_em436_{BASE}_dsempty.npz")):
        if not f.exists() or not ref.exists():
            continue
        e, bad = eq(f, ref)
        ka[nm] = "PASS (bit-equal)" if e else f"FAIL {bad}"
        print(nm, "vs 435", ref.name, "->", ka[nm])
# KA3: the emitter after the `capped` edit (v2), plain mode, re-run of seed 0's donor-signed arm == the v1 output
f3, r3 = W / f"stats-base_s0_bps_em436_{BASE}_dsall_ka3.npz", W / f"stats-base_s0_bps_em436_{BASE}_dsall.npz"
if f3.exists() and r3.exists():
    e, bad = eq(f3, r3)
    ka["KA3 s0 dsall v2-emitter == v1"] = "PASS (bit-equal)" if e else f"FAIL {bad}"
    print("KA3 s0 dsall (emitter v2, plain) vs v1 ->", ka["KA3 s0 dsall v2-emitter == v1"])
if ka:
    (HERE / "work/ka436.json").write_text(json.dumps(ka, indent=1))
if any(v.startswith("FAIL") for v in ka.values()):
    raise SystemExit("known answer FAILED: read nothing below")

sys.path.insert(0, str(HERE))
import signs436 as S6  # noqa: E402

rows, tab = [], {"variant": VAR, "seeds": {}}
for s in (0, 1, 2):
    off_m = W / f"metrics-base_s{s}_bps_em436_{BASE}.json"
    on_m = W / f"metrics-base_s{s}_bps_em436_{BASE}_{VAR}.json"
    if not (off_m.exists() and on_m.exists()):
        continue
    mo, mn = json.loads(off_m.read_text()), json.loads(on_m.read_text())
    cov = set(mn["sign_info"]["covered_targets"])
    po = pl.read_csv(W / f"per-target-base_s{s}_bps_em436_{BASE}.csv")
    pn = pl.read_csv(W / f"per-target-base_s{s}_bps_em436_{BASE}_{VAR}.csv")
    j = po.join(pn, on=["perturbation", "metric"], suffix="_on")
    rec = {"all": {k: mn["raw"][k] - mo["raw"][k] for k in MEM}}
    for nm_, keep in (("covered", True), ("uncovered", False)):
        d = {}
        for k, pk in PT.items():
            sub = j.filter((pl.col("metric") == pk) & (pl.col("perturbation").is_in(list(cov)) == keep))
            v = (sub["value_on"] - sub["value"]).to_numpy()
            d[k] = float(np.nanmean(v)) if v.size else float("nan")
            d[f"n_{k}"] = int(np.isfinite(sub["value"].to_numpy()).sum())
            d[f"moved_{k}"] = int((np.abs(v) > 0).sum())
        rec[nm_] = d
    # the overall-mean nmae change decomposed: covered share x covered mean
    nm_all = j.filter(pl.col("metric") == PT["nmae"])
    rec["nmae_decomp_check"] = float(np.nanmean((nm_all["value_on"] - nm_all["value"]).to_numpy()))
    rec["overall_all"] = sum(SLOPE[k] * rec["all"][k] for k in MEM) / 6
    # blocks: mde_lambda / prior_lambda, capacity diagnostics
    do = json.loads((W / f"diag-base_s{s}_bps_em436_{BASE}.json").read_text())
    dn = json.loads((W / f"diag-base_s{s}_bps_em436_{BASE}_{VAR}.json").read_text())
    blk = []
    for a_, b_ in zip(do, dn):
        assert a_["target"] == b_["target"]
        blk.append({"target": a_["target"], "covered": a_["target"] in cov,
                    "mde_lambda_off": a_.get("mde_lambda"), "mde_lambda_on": b_.get("mde_lambda"),
                    "prior_lambda_off": a_.get("prior_lambda"), "prior_lambda_on": b_.get("prior_lambda"),
                    "dclass_off": a_.get("mde_dclass_max"), "dclass_on": b_.get("mde_dclass_max"),
                    "dclass_full_on": b_.get("mde_dclass_full"), "slack_on": b_.get("mde_slack_min"),
                    "overfill_on": b_.get("mde_overfill_frac"), "fu_on": b_.get("mde_fu"), "fd_on": b_.get("mde_fd"),
                    "fu_off": a_.get("mde_fu"), "fd_off": a_.get("mde_fd"),
                    "n_up_off": a_.get("mde_n_up"), "n_dn_off": a_.get("mde_n_dn"),
                    "n_up_on": b_.get("mde_n_up"), "n_dn_on": b_.get("mde_n_dn"),
                    "sign_n": b_.get("mde_sign_n"), "sign_rev": b_.get("mde_sign_n_rev"), "sign_rev_blk": b_.get("mde_sign_n_rev_blk"),
                    "unsupp_off": a_.get("call_set_unsuppressible"), "unsupp_on": b_.get("call_set_unsuppressible"),
                    "clip_off": a_.get("clip_frac"), "clip_on": b_.get("clip_frac"),
                    "calls_off": a_.get("call_set_delivered"), "calls_on": b_.get("call_set_delivered"),
                    "cs_jac_off": a_.get("call_set_jaccard"), "cs_jac_on": b_.get("call_set_jaccard"),
                    "fu0_on": b_.get("mde_fu0"), "fd0_on": b_.get("mde_fd0")})
    rec["blocks"] = blk
    # the call-leak read (435's dose-ladder criterion: shipped dose = ~5 unsuppressible calls/target)
    cb = [b for b in blk if b["covered"]]
    rec["leak"] = {"covered_blocks": len(cb),
                   "unsupp_mean_off": float(np.mean([b["unsupp_off"] for b in cb])), "unsupp_mean_on": float(np.mean([b["unsupp_on"] for b in cb])),
                   "blocks_unsupp_gt50_on": int(sum(b["unsupp_on"] > 50 for b in cb)), "blocks_unsupp_gt50_off": int(sum(b["unsupp_off"] > 50 for b in cb)),
                   "extra_calls_sum": int(sum(b["calls_on"] - b["calls_off"] for b in cb)),
                   "extra_calls_max": int(max(b["calls_on"] - b["calls_off"] for b in cb)),
                   "mde_lambda_lt1_on": int(sum((b["mde_lambda_on"] if b["mde_lambda_on"] is not None else 1) < 1 for b in cb)),
                   "mde_lambda_lt1_off": int(sum((b["mde_lambda_off"] if b["mde_lambda_off"] is not None else 1) < 1 for b in cb)),
                   "prior_lambda_lt1_on": int(sum((b["prior_lambda_on"] if b["prior_lambda_on"] is not None else 1) < 1 for b in cb))}
    # stats-level identity of the uncovered targets (does the shared rng stream shift after a covered block?)
    zo, zn = np.load(W / f"stats-base_s{s}_bps_em436_{BASE}.npz"), np.load(W / f"stats-base_s{s}_bps_em436_{BASE}_{VAR}.npz")
    tg = zo["targets"].astype(str)
    ident = [bool(np.array_equal(zo["sums"][k], zn["sums"][k]) and np.array_equal(zo["lfc"][k], zn["lfc"][k], equal_nan=True))
             for k in range(len(tg))]
    rec["uncovered_bit_identical"] = int(sum(i for i, t in zip(ident, tg) if t not in cov))
    rec["uncovered_n"] = int(sum(t not in cov for t in tg))
    rec["covered_sums_identical"] = int(sum(bool(np.array_equal(zo["sums"][k], zn["sums"][k])) for k, t in enumerate(tg) if t in cov))
    # delivered LFC change on covered targets vs the intended donor-signed add (on the genes the donor sign reverses)
    gk = np.load(ROOT / "experiments/434-posterior-median-magnitude/work/allcells434.npz", allow_pickle=True)["genes"].astype(str)
    Ssig, _ = S6.donor_signs(tg, gk, mode=("allx" if "allx" in VAR else "all" if "all" in VAR else "ess"),
                             min_agree=int(VAR.split("ag")[1]) if "ag" in VAR else 1)
    lo, ln = np.nan_to_num(zo["lfc"]), np.nan_to_num(zn["lfc"])
    dl = []
    for k, t in enumerate(tg):
        if t in cov:
            rev = (Ssig[k] != 0) & (Ssig[k] != np.sign(lo[k])) & (lo[k] != 0)
            agr = (Ssig[k] != 0) & (Ssig[k] == np.sign(lo[k]))
            dl.append({"target": t, "rev_genes": int(rev.sum()),
                       "rev_mean_signed_move": float(np.mean(Ssig[k][rev] * (ln[k][rev] - lo[k][rev]))) if rev.any() else 0.0,
                       "agree_mean_signed_move": float(np.mean(Ssig[k][agr] * (ln[k][agr] - lo[k][agr]))) if agr.any() else 0.0,
                       "delivered_sign_changes": int((np.sign(ln[k]) != np.sign(lo[k])).sum())})
    rec["delivered_lfc"] = dl
    tab["seeds"][s] = rec
    rows.append(rec)
    c, u = rec["covered"], rec["uncovered"]
    lam_c = [b["mde_lambda_on"] for b in blk if b["covered"]]
    print(f"seed {s} {VAR}: all " + " ".join(f"{k} {rec['all'][k]:+.5f}" for k in MEM) + f"  -> overall {rec['overall_all']:+.5f}")
    print(f"   covered ({c['n_nmae']} gated / {c['n_pds']} all): " + " ".join(f"{k} {c[k]:+.5f}" for k in PT)
          + f" | uncovered ({u['n_nmae']} / {u['n_pds']}): " + " ".join(f"{k} {u[k]:+.5f}" for k in PT))
    print(f"   uncovered bit-identical {rec['uncovered_bit_identical']}/{rec['uncovered_n']}; covered mde_lambda on "
          f"median {np.nanmedian([x if x is not None else np.nan for x in lam_c]):.3f} min {np.nanmin([x if x is not None else np.nan for x in lam_c]):.3f}; "
          f"prior_lambda<1 on/off {sum((b['prior_lambda_on'] or 0) < 1 for b in blk)}/{sum((b['prior_lambda_off'] or 0) < 1 for b in blk)}")
    print(f"   leak: {rec['leak']}")
if rows:
    for part in ("all", "covered", "uncovered"):
        ks = MEM if part == "all" else list(PT)
        arr = {k: np.array([r[part][k] for r in rows]) for k in ks}
        tab[f"mean_{part}"] = {k: float(v.mean()) for k, v in arr.items()}
        tab[f"sd_{part}"] = {k: float(v.std(ddof=1)) if len(v) > 1 else 0.0 for k, v in arr.items()}
        print(f"{part:9s} mean±sd: " + " ".join(f"{k} {arr[k].mean():+.5f}±{(arr[k].std(ddof=1) if len(rows) > 1 else 0):.5f}" for k in ks))
    ov = np.array([r["overall_all"] for r in rows])
    tab["overall_all_mean"] = float(ov.mean()); tab["overall_all_sd"] = float(ov.std(ddof=1)) if len(ov) > 1 else 0.0
    print(f"overall (H1, 435 slopes): {ov.mean():+.5f} ± {tab['overall_all_sd']:.5f}")
(HERE / f"work/tab436-{VAR}.json").write_text(json.dumps(tab, indent=1, default=float))
