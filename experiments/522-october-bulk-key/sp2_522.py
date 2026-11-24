"""522 — does the src-strength rider survive 479's per-target sharpening? The October recipe will carry S_P2 (stercile
2.5 / 2.0 / 1.75 by the unkeyed stack row's non-panel norm; 546's default partner) or P2g (key terciles). Arms per
surface, all exact (fold522 with a per-target q vector; rows sharpened per q group, as comp546o / arms540c):
  pct + S_P2, pct x S0.25 + S_P2, pct x S0.5 + S_P2, cap 6 + S_P2 (and the uniform-q pct, pct x S0.25 references)
KNOWN ANSWERS (abort unless |d| < 1e-12): October-like folds: `p96.6 + S_P2` (P 96.647) == 546o's row (L.cmp);
A/B/C folds: `cap 6 + S_P2` == 540c's `V1_0.5_6 + S_P2` row (cos_agg, pds).

    uv run python experiments/522-october-bulk-key/sp2_522.py        (~1-2 GB)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import fold522 as FD  # noqa: E402
import lib438 as L  # noqa: E402
from lib438 import SLOTS  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

WORK = HERE / "work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
P2 = (2.5, 2.0, 1.75)
TOL = 1e-12


def terc(v, spec):
    t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
    return np.asarray([spec[k] for k in t], float)


def score_q(S, w, qv):
    """S.score with a per-target q vector (rows sharpened per q group; the row mode is per-row, so grouping is exact)."""
    if S.kind == "oct":
        bulk = L.key_scale(S.fin, np.asarray(w, float))[0]
        rea, K = {}, []
        for c, (b, gate) in S.frames.items():
            F = np.empty_like(bulk)
            for q in np.unique(qv):
                sel = qv == q
                F[sel] = L.sharpen(bulk[sel], float(q), gate, mode="row")
            R = np.empty_like(F)
            for i in range(S.n):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                p_re, _, k, _, _ = L.realise_k(S.P_DE[c][i], p_bk, 16)
                R[i] = np.log1p(TS * p_re) - b
                K.append(k)
            rea[c] = R
        s, arr = L.row_score(rea, S.scored, S.truth, S.betas, S.nonpanel)
    else:
        from lib420 import TS as TS4, key_scale, profile_from_effect as pfe, realise_k, sharpen
        from comp420 import score_rows
        F0 = S.F
        bulk = key_scale(F0.fST, np.asarray(w, float))[0]
        rea, K = {}, []
        for c, (b, g) in F0.fr.items():
            X = np.empty_like(bulk)
            for q in np.unique(qv):
                sel = qv == q
                X[sel] = sharpen(bulk[sel], float(q), g, mode="row")
            R = np.empty_like(X)
            for i in range(F0.n):
                p_bk = pfe(X[i], b, g, "minnorm")
                p_re, _, k, _, _ = realise_k(F0.P_DE[c][i], p_bk, 16)
                R[i] = np.log1p(TS4 * p_re) - b
                K.append(k)
            rea[c] = R
        s, arr = score_rows(rea, F0.scored, F0.truth, F0.betas, F0.nonpanel)
    s["K_mean"] = float(np.mean(K))
    return s, arr


def run(panel, held):
    t0 = time.time()
    S = FD.load(panel, held)
    z = FD.zsrc(S)
    kv, car = S.kv, S.carried
    if S.kind == "oct":
        strength = np.linalg.norm(S.fin[:, S.nonpanel], axis=1)
    else:
        strength = np.linalg.norm(S.F.fST[:, S.F.nonpanel], axis=1)
    qs = terc(strength, P2)
    q1 = np.full(S.n, 1.75)
    w3 = FD.wclip(kv, car, 0.5, 3.0)
    w6 = FD.wclip(kv, car, 0.5, 6.0)
    wp = FD.wclip(kv, car, 0.5, FD.pct_hi(kv, car, FD.P_PROD))
    ka = {}
    if S.kind == "oct":
        w546 = FD.wclip(kv, car, 0.5, FD.pct_hi(kv, car, FD.P_546))
        s, _ = score_q(S, w546, qs)
        ref = json.loads((FD.W546 / f"comp546o-{panel}-{held}.json").read_text())["arms"]["p96.6 + S_P2"]["realised"]["16"]
        d = L.cmp(s, ref)
        assert d < TOL, f"KA p96.6 + S_P2 vs 546o {panel}/{held}: {d}"
        ka["p96.6 + S_P2 == 546o"] = d
    else:
        s, _ = score_q(S, w6, qs)
        ref = json.loads((W540 / f"arms540c-{held}.json").read_text())["arms"]["V1_0.5_6 + S_P2"]["realised"]["16"]
        d = max(abs(s["cos_agg"] - ref["cos_agg"]), abs(s["pds"] - ref["pds"]))
        assert d < TOL, f"KA cap6 + S_P2 vs 540c {held}: {d}"
        ka["cap6 + S_P2 == 540c"] = d
    arms = {"cap3": (w3, q1), "pct": (wp, q1), "pctxS0.25": (wp * np.exp(0.25 * z), q1),
            "cap6 + S_P2": (w6, qs), "pct + S_P2": (wp, qs), "pctxS0.25 + S_P2": (wp * np.exp(0.25 * z), qs),
            "pctxS0.5 + S_P2": (wp * np.exp(0.5 * z), qs)}
    res = dict(panel=panel, fold=held, ka=ka, arms={})
    base = None
    for nm, (w, qv) in arms.items():
        s, _ = score_q(S, w, qv)
        if base is None:
            base = s
        res["arms"][nm] = dict(cos_agg=s["cos_agg"], pds=s["pds"], r=s["cos_agg"] / base["cos_agg"] - 1,
                               dpds=s["pds"] - base["pds"], K_mean=s["K_mean"])
    a = res["arms"]
    res["rider_under_S_P2"] = a["pctxS0.25 + S_P2"]["cos_agg"] / a["pct + S_P2"]["cos_agg"] - 1
    res["rider_plain"] = a["pctxS0.25"]["cos_agg"] / a["pct"]["cos_agg"] - 1
    res["S_P2_on_pct"] = a["pct + S_P2"]["cos_agg"] / a["pct"]["cos_agg"] - 1
    print(f"[{panel}/{held}] KA {ka}; S_P2 on pct {100 * res['S_P2_on_pct']:+.2f} %; rider plain {100 * res['rider_plain']:+.2f} %, "
          f"under S_P2 {100 * res['rider_under_S_P2']:+.2f} %; dpds rider|S_P2 {a['pctxS0.25 + S_P2']['pds'] - a['pct + S_P2']['pds']:+.5f} "
          f"({time.time() - t0:.0f}s)", flush=True)
    return res


def main():
    out = {}
    for p, f in FD.surfaces(("abc", "efp", "efr", "real")):
        if p == "abc" and f == "kolf2":
            continue
        out[f"{p}/{f}"] = run(p, f)
        (WORK / "sp2_522.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
