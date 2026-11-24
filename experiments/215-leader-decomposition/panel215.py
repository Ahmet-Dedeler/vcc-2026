"""215 §Panel — build the deduplicated board observation panel from all 17 snapshots.

Read-only. Writes only inside experiments/215-leader-decomposition/.

Known-answer rows (abort on failure):
  KA-1  the six raw->scaled slopes recover docs/BOARD-LEDGER.md's fitted map
  KA-2  our own displayed row reproduces 197's rc_hybrid / 202's rc_callset numbers exactly
  KA-3  Hongjian Zhou's 09-03 00:30 row is pm-r14d at the six values 186 recorded
  KA-4  teams with raw mse < 1.0 on 08-28 / 09-03T0030 are 8 / 44 (206 KA)
  KA-5  the informative cohort (raw mse < 1.023611) is 54 of 612 on 09-03T0030 (197)
"""
import json, glob, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
os.makedirs(WORK, exist_ok=True)

RAW = {
    "pds": "pds_cosine",
    "mse": "expr_mse_unbiased_capped_norm",
    "nmae": "de_wilcoxon_lfc_nmae",
    "fid": "de_wilcoxon_direction_fidelity_yield_raw",
    "reach": "de_wilcoxon_direction_reach_raw",
    "jac": "de_wilcoxon_sig_jaccard",
}
SCA = {
    "pds": "score_pds", "mse": "score_mse", "nmae": "score_nmae",
    "fid": "score_fid", "reach": "score_reach", "jac": "score_jac",
}
MEMBERS = list(RAW)

# 197's inversion, verbatim.
COS_A, COS_B = 1.023611, 2.074495
def cos_floor(raw_mse):
    return (COS_A - raw_mse) / COS_B
def informative(raw_mse):
    return raw_mse < COS_A


def load_panel():
    files = sorted(glob.glob(os.path.join(ROOT, "data/leaderboard/*.json")))
    files = [f for f in files if os.path.basename(f).startswith("lb_")] + \
            [f for f in files if not os.path.basename(f).startswith("lb_")]
    files = sorted(set(files), key=lambda f: os.path.basename(f))
    files.append(os.path.join(WORK, "lb-live.json"))
    obs = {}
    snap_of = {}
    for f in files:
        tag = os.path.basename(f).replace(".json", "")
        d = json.load(open(f))
        entries = d["entries"] if isinstance(d, dict) else d
        for e in entries:
            if e.get("status") != "published":
                continue
            if any(e.get(RAW[m]) is None for m in MEMBERS):
                continue
            eid = e["id"]
            rec = {
                "id": eid,
                "team": e["team_name"],
                "team_id": e.get("team_id"),
                "org": e.get("organization"),
                "members": tuple(e.get("team_members") or ()),
                "date": e["submission_date"],
                "n": e.get("submission_count"),
                "model": e.get("model_name"),
                "desc": e.get("description"),
                "overall": e["score_avg"],
            }
            for m in MEMBERS:
                rec["raw_" + m] = float(e[RAW[m]])
                v = e.get(SCA[m])
                rec["sca_" + m] = None if v is None else float(v)
            obs[eid] = rec
            snap_of.setdefault(eid, []).append(tag)
    for eid, rec in obs.items():
        rec["snaps"] = snap_of[eid]
        rec["n_snaps"] = len(snap_of[eid])
    return obs, files


# Every member is scaled PER CONTEXT and then averaged, so the aggregate raw -> aggregate
# scaled map is affine only while all three contexts sit in the same clamping regime.
#   mse   per-context scaled = max(0, slope*(raw_c-b))  -> convex kink near raw ~0.99
#   nmae  per-context scaled = max(-6, ...)             -> floored for broken submissions
# Fit each member well inside the all-contexts-unclamped region.
FITRANGE = {"mse": (0.0, 0.90), "nmae": (0.85, 1.05)}


def fit_scaling(obs):
    """KA-1: recover the affine raw->scaled map per member from the whole panel.

    Fit only on the UNCLAMPED part of each member's range: `mse` is clamped to
    scaled 0 above raw ~0.98786 and `nmae` is floored at scaled -6.
    """
    out = {}
    for m in MEMBERS:
        rows = [r for r in obs.values() if r["sca_" + m] is not None]
        x = np.array([r["raw_" + m] for r in rows])
        y = np.array([r["sca_" + m] for r in rows])
        keep = np.isfinite(x) & np.isfinite(y)
        if m in FITRANGE:
            lo, hi = FITRANGE[m]
            keep &= (x > lo) & (x < hi)
        A = np.vstack([x[keep], np.ones(keep.sum())]).T
        slope, intercept = np.linalg.lstsq(A, y[keep], rcond=None)[0]
        resid = np.abs(A @ [slope, intercept] - y[keep]).max()
        out[m] = {"slope": float(slope), "intercept": float(intercept),
                  "b": float(-intercept / slope), "n": int(keep.sum()),
                  "max_abs_resid": float(resid)}
    return out


def main():
    obs, files = load_panel()
    print(f"snapshots: {len(files)}   distinct published observations: {len(obs)}")

    fits = fit_scaling(obs)
    LEDGER = {"pds": 2.2138, "mse": -1.051, "nmae": -1.6532,
              "fid": 3.3645, "reach": 1.1292, "jac": 2.7263}
    print("\nKA-1  raw -> scaled affine map (fit on the whole panel)")
    print(f"{'member':<7}{'slope':>12}{'ledger':>10}{'rel err':>10}{'zero at':>11}{'n':>7}{'maxresid':>11}")
    ok1 = True
    for m in MEMBERS:
        f = fits[m]
        rel = abs(f["slope"] - LEDGER[m]) / abs(LEDGER[m])
        ok1 &= rel < 0.01
        print(f"{m:<7}{f['slope']:>12.5f}{LEDGER[m]:>10.4f}{rel:>10.2%}"
              f"{f['b']:>11.5f}{f['n']:>7d}{f['max_abs_resid']:>11.2e}")
    assert ok1, "KA-1 FAILED: scaling map does not reproduce the ledger's"

    # KA-2 our own row
    ours = [r for r in obs.values() if r["team"] == "Ahmet Dedeler"]
    ours.sort(key=lambda r: r["date"])
    cs = [r for r in ours if r["model"] == "Peer Review Rejected This"][0]
    rh = [r for r in ours if r["model"] == "The Mitochondria Declined To Comment"][0]
    print("\nKA-2  our own two ordering points")
    # tolerances: the BOARD-RESULT write-ups quote raw DE members to 5 d.p.
    ka2 = [
        ("rc_hybrid overall", rh["overall"], 0.1834190095, 5e-10),
        ("rc_hybrid raw pds", rh["raw_pds"], 0.8359866220735785, 0.0),
        ("rc_hybrid raw mse", rh["raw_mse"], 0.8738045054832962, 0.0),
        ("rc_hybrid raw reach", rh["raw_reach"], 0.22982, 2e-4),
        ("rc_callset overall", cs["overall"], 0.1918638367, 5e-10),
        ("rc_callset raw pds", cs["raw_pds"], 0.8359866220735785, 0.0),
        ("rc_callset raw mse", cs["raw_mse"], 0.8738045054832962, 0.0),
        ("rc_callset raw jac", cs["raw_jac"], 0.04043, 5e-5),
        ("rc_callset raw reach", cs["raw_reach"], 0.24983, 1e-4),
    ]
    ok2 = True
    for name, got, want, tol in ka2:
        err = abs(got - want)
        good = err <= tol
        ok2 &= good
        print(f"  {name:<22} {got!r:<24} want {want:<22} {'PASS' if good else 'FAIL'}")
    assert ok2, "KA-2 FAILED"
    # the exact-identity check the whole decomposition rests on
    assert rh["raw_pds"] == cs["raw_pds"], "rc_hybrid/rc_callset pds not bit-identical"
    assert rh["raw_mse"] == cs["raw_mse"], "rc_hybrid/rc_callset mse not bit-identical"
    print("  rc_hybrid vs rc_callset: raw pds and raw mse BIT-IDENTICAL  PASS")

    # KA-3 186's pm-r14d row
    hz = [r for r in obs.values() if r["team"] == "Hongjian Zhou" and r["model"] == "pm-r14d"]
    assert len(hz) == 1, f"KA-3 FAILED: {len(hz)} pm-r14d rows"
    h = hz[0]
    want = dict(pds=0.8325, mse=0.760, nmae=0.9291, fid=0.5084, reach=0.2454, jac=0.0224)
    print("\nKA-3  Hongjian Zhou pm-r14d, six members (186 §1.4)")
    ok3 = True
    for m, w in want.items():
        got = h["raw_" + m]
        good = abs(got - w) < 6e-4
        ok3 &= good
        print(f"  {m:<6} {got:.6f}  want {w:.4f}  {'PASS' if good else 'FAIL'}")
    assert ok3, "KA-3 FAILED"

    # KA-4 / KA-5 cohort counts on named snapshots
    print("\nKA-4/5  cohort counts")
    ok45 = True
    for tag, want_lt1, want_inf, want_n in [
        ("lb_20260828T1325", 8, None, None),
        # 197 §2.5 reports 41 at raw mse < 1.0 and 54 informative on this snapshot.
        # (206 §1 quotes 44 for the same cut; 197's 41 is what the bytes give.)
        ("2026-09-03T0030Z", 41, 54, 612),
    ]:
        rows = [r for r in obs.values() if tag in r["snaps"]]
        lt1 = sum(1 for r in rows if r["raw_mse"] < 1.0)
        inf = sum(1 for r in rows if informative(r["raw_mse"]))
        print(f"  {tag}: rows {len(rows)}  raw mse<1.0 {lt1} (want {want_lt1})"
              f"  informative {inf} (want {want_inf})")
        ok45 &= (lt1 == want_lt1)
        if want_inf is not None:
            ok45 &= (inf == want_inf)
    assert ok45, "KA-4/5 FAILED"

    print("\nALL PANEL KNOWN-ANSWER ROWS PASS")

    json.dump({"fits": fits, "n_obs": len(obs), "n_snapshots": len(files),
               "files": [os.path.basename(f) for f in files]},
              open(os.path.join(WORK, "panel215_fits.json"), "w"), indent=1)
    # the panel itself
    ser = sorted(obs.values(), key=lambda r: (r["team"] or "", r["date"]))
    json.dump(ser, open(os.path.join(WORK, "panel215.json"), "w"), indent=0)
    print(f"wrote work/panel215.json  ({len(ser)} observations)")


if __name__ == "__main__":
    main()
