"""234 Phase 0 — the mandatory known-answer rows.  Nothing downstream runs if one fails.

KA-1  232's SIGN-ACCURACY EXCHANGE RATE and the oracle `nmae` response slope, re-derived
      from 232's stored `arm232.json` by an INDEPENDENT implementation of the six-member
      pricing (the constants are the published spans; the code is written fresh here).
      Targets: composite +0.003990 / +0.003785 per +0.01 `q_w` (k_cos 1.384 / 1.000),
      raw d`nmae` -0.002573 per +0.01 with r = -0.995.

KA-2  228's `e = dM/dU = 0.072` for 224's per-gene magnitude gain, and its d`nmae` +0.0373,
      recomputed from 228's stored `adel228.json` through 188's identity `1 + U/D - 2M/D`.
      This is the number the SIGN route has to beat, so it is a load-bearing row.

KA-3  220 §7's measured key correlations with `1[our sign is right]` inside the reference DE
      pool on h1-test: `|H|` +0.0743 (AUC 0.5350), 189's `theta` +0.002, control expression
      +0.0080 (AUC 0.5049).  These are the SAME keys 234 is about to price as flip keys, so
      the row is not a formality -- if the keys do not reproduce, the screen is meaningless.
      Gated by 220's own `check_old` (rule 10) on the a15 instrument.

    uv run python phase0_234.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common234 import (KA220_AUC, KA220_CORR, KA228_E, KA232_EXCHANGE,  # noqa: E402
                       QW_NEEDED, ROOT, WORK)

CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")

# ---- published spans (AGENTS.md §4 / 228 / 232 `price232.py` header).  Constants, not code.
SPAN_PDS, CALIB = 0.452, 1.011
SPAN_MSE, DMSE_DCOS = 0.950162566, 2.074495
SPAN_NMAE, SPAN_FID, SPAN_JAC = 0.5963, 0.30, 0.365
OURS_COS, OURS_JAC_RAW = 0.09946902740184964, 0.04095656108311888


def six_member(base_rows, arm_rows, k_cos: float) -> dict:
    """An independent re-implementation of 228/232's pricing routes.

    `pds`  d pds / span / 6 * calib          `mse`  197's inversion through rel cosine
    `nmae` -d nmae / span / 6                `fid`  d call-set sign precision / span / 6
    `jac`  RELATIVE, rescaled to the board's own raw `jac`
    """
    def m(fn):
        return float(np.mean([fn(a, b) for a, b in zip(arm_rows, base_rows)]))

    rel = m(lambda a, b: (a["cos_agg"] - b["cos_agg"]) / abs(b["cos_agg"]))
    o_pds = m(lambda a, b: a["pds"] - b["pds"]) / SPAN_PDS / 6.0 * CALIB
    o_mse = (DMSE_DCOS * OURS_COS * rel * k_cos) / SPAN_MSE / 6.0
    o_nmae = -m(lambda a, b: a["nmae_pools"]["wald89"] - b["nmae_pools"]["wald89"]) \
        / SPAN_NMAE / 6.0
    o_fid = m(lambda a, b: a["fid_precision"] - b["fid_precision"]) / SPAN_FID / 6.0
    o_jac = m(lambda a, b: a["jac"]["89"] / b["jac"]["89"] - 1.0) * OURS_JAC_RAW \
        / SPAN_JAC / 6.0
    return dict(rel_cos=rel, o_pds=o_pds, o_mse=o_mse, o_nmae=o_nmae, o_fid=o_fid,
                o_jac=o_jac, composite=o_pds + o_mse + o_nmae + o_fid + o_jac,
                d_qw=m(lambda a, b: a["q_w_pool"]["wald89"] - b["q_w_pool"]["wald89"]),
                d_nmae=m(lambda a, b: a["nmae_pools"]["wald89"] - b["nmae_pools"]["wald89"]))


def ka1() -> dict:
    rows = json.loads((ROOT / "experiments/232-sign-vote/work/arm232.json").read_text())
    R = {(r["fold"], r["arm"]): r for r in rows}
    orc = sorted({r["arm"] for r in rows
                  if ("oracle@" in r["arm"] or "antioracle@" in r["arm"])
                  and all((f, r["arm"]) in R for f in CANCER)})
    base = [R[(f, "BASE")] for f in CANCER]
    x, y = [], {"composite_k1384": [], "composite_k1000": [], "d_nmae": [],
                "o_pds": [], "o_nmae": []}
    for arm in orc:
        a = [R[(f, arm)] for f in CANCER]
        p384, p100 = six_member(base, a, 1.384), six_member(base, a, 1.000)
        x.append(p384["d_qw"])
        y["composite_k1384"].append(p384["composite"])
        y["composite_k1000"].append(p100["composite"])
        y["d_nmae"].append(p384["d_nmae"])
        y["o_pds"].append(p384["o_pds"])
        y["o_nmae"].append(p384["o_nmae"])
    x = np.asarray(x)
    out = {"n_oracle_rows": len(orc), "arms": orc}
    for k, v in y.items():
        v = np.asarray(v)
        slope = float(np.linalg.lstsq(x[:, None], v, rcond=None)[0][0])
        out[k] = dict(per_0p01_qw=slope * 0.01, r=float(np.corrcoef(x, v)[0, 1]))
    ok = (abs(out["composite_k1384"]["per_0p01_qw"] - KA232_EXCHANGE["composite_k1384"]) < 5e-6
          and abs(out["composite_k1000"]["per_0p01_qw"]
                  - KA232_EXCHANGE["composite_k1000"]) < 5e-6
          and abs(out["d_nmae"]["per_0p01_qw"] - KA232_EXCHANGE["d_nmae"]) < 5e-6
          and abs(out["d_nmae"]["r"] - KA232_EXCHANGE["r_nmae"]) < 5e-3)
    print(f"[KA-1] 232 exchange rate over {len(orc)} oracle/antioracle doses, "
          "independent pricing code")
    for k, exp in (("composite_k1384", KA232_EXCHANGE["composite_k1384"]),
                   ("composite_k1000", KA232_EXCHANGE["composite_k1000"]),
                   ("d_nmae", KA232_EXCHANGE["d_nmae"]),
                   ("o_pds", KA232_EXCHANGE["o_pds"]),
                   ("o_nmae", KA232_EXCHANGE["o_nmae"])):
        print(f"        {k:<18} {out[k]['per_0p01_qw']:>+11.6f} per +0.01 q_w "
              f"(232: {exp:+.6f})   r {out[k]['r']:+.3f}")
    print(f"        -> gate +0.0019 needs d q_w = +{QW_NEEDED['k1.384']:.5f} (k1.384) .. "
          f"+{QW_NEEDED['k1.000']:.5f} (k1.000)   {'PASS' if ok else 'FAIL'}", flush=True)
    out["PASS"] = bool(ok)
    assert ok, "KA-1 FAILED: 232's exchange rate does not reproduce"
    return out


def ka2() -> dict:
    d = json.loads((ROOT / "experiments/228-theta-gain-plolo/work/adel228.json").read_text())
    got = {}
    for arm in ("BASE", "gain_reg@0.25"):
        r = [x for x in d if x["arm"] == arm and x["fold"] in CANCER]
        assert len(r) == 3, f"{arm}: {len(r)} cancer rows"
        U = float(np.mean([x["nmae_U"]["wald89"] for x in r]))
        M = float(np.mean([x["nmae_M"]["wald89"] for x in r]))
        nm = float(np.mean([x["nmae_pools"]["wald89"] for x in r]))
        assert abs((1.0 + U - 2.0 * M) - nm) < 5e-6, f"188 identity broken for {arm}"
        got[arm] = dict(U=U, M=M, nmae=nm)
    dU = got["gain_reg@0.25"]["U"] - got["BASE"]["U"]
    dM = got["gain_reg@0.25"]["M"] - got["BASE"]["M"]
    dn = got["gain_reg@0.25"]["nmae"] - got["BASE"]["nmae"]
    e = dM / dU
    ok = abs(e - KA228_E) < 0.002 and abs(dn - 0.0373) < 0.0005
    print(f"[KA-2] 228's magnitude arm: dU {dU:+.5f}  dM {dM:+.5f}  e = {e:.4f} "
          f"(228: {KA228_E})  d_nmae {dn:+.4f} (228: +0.0373)  {'PASS' if ok else 'FAIL'}",
          flush=True)
    out = dict(dU=dU, dM=dM, e=e, d_nmae=dn, PASS=bool(ok), **got)
    assert ok, "KA-2 FAILED"
    return out


def ka3() -> dict:
    sys.path.insert(0, str(ROOT / "experiments/220-ordering-keys-retry"))
    import common220 as C220                                        # noqa: E402
    T = C220.Instrument(C220.DE_RC, "rc_callset")
    old = C220.Instrument(C220.DE_A15, "a15")
    C220.check_old(old)                                             # rule 10 gate
    pool = T.pool
    y = T.correct.astype(np.float64)

    def auc(sig):
        s = np.asarray(sig)[pool]
        t = y[pool] > 0.5
        o = np.argsort(s, kind="stable")
        sv, r = s[o], np.empty(s.size)
        i = 0
        while i < s.size:
            j = i
            while j + 1 < s.size and sv[j + 1] == sv[i]:
                j += 1
            r[o[i:j + 1]] = 0.5 * (i + 1 + j + 1)
            i = j + 1
        n1, n0 = int(t.sum()), int((~t).sum())
        return float((r[t].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

    sigs = dict(absH=np.abs(T.H),
                theta=np.broadcast_to(T.theta[None, :], T.H.shape),
                expr=np.broadcast_to(T.b_ctrl[None, :], T.H.shape),
                key189=C220.key189(T.H, T.theta, 0.5))
    out = {}
    for k, s in sigs.items():
        out[k] = dict(corr=float(np.corrcoef(np.asarray(s)[pool], y[pool])[0, 1]),
                      auc=auc(s))
    ok = (abs(out["absH"]["corr"] - KA220_CORR["absH"]) < 5e-4
          and abs(out["absH"]["auc"] - KA220_AUC["absH"]) < 5e-4
          and abs(out["theta"]["corr"] - KA220_CORR["theta"]) < 5e-4
          and abs(out["expr"]["corr"] - KA220_CORR["expr"]) < 5e-4
          and abs(out["expr"]["auc"] - KA220_AUC["expr"]) < 5e-4)
    print("[KA-3] 220 §7 key correlations with `1[our sign is right]`, h1-test rc_callset "
          f"pool ({int(pool.sum()):,} pairs)")
    for k in ("absH", "theta", "expr", "key189"):
        exp = KA220_CORR.get(k)
        print(f"        {k:<8} corr {out[k]['corr']:+.4f}   AUC {out[k]['auc']:.4f}"
              + (f"   (220: {exp:+.4f})" if exp is not None else "   (not in 220)"))
    print(f"        required by 220's spec: corr ~ +0.30 / AUC ~ 0.639   "
          f"{'PASS' if ok else 'FAIL'}", flush=True)
    out["PASS"] = bool(ok)
    out["n_pool_pairs"] = int(pool.sum())
    out["base_q_unweighted"] = float(y[pool].mean())
    assert ok, "KA-3 FAILED: 220's key correlations do not reproduce"
    return out


def main() -> None:
    t0 = time.time()
    audit = dict(KA_1=ka1(), KA_2=ka2(), KA_3=ka3(),
                 qw_needed=QW_NEEDED, seconds=None)
    audit["seconds"] = time.time() - t0
    (WORK / "phase0_234.json").write_text(json.dumps(audit, indent=1))
    print(f"\n-> {WORK / 'phase0_234.json'}  ({audit['seconds']:.0f}s)  "
          "ALL KNOWN-ANSWER ROWS PASS")


if __name__ == "__main__":
    main()
