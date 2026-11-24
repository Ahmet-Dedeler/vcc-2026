"""521 known answers (abort on failure).

KA1  192's per-perturbation sizing table at f = 1.0 with 192's OWN code (frontier192.build_keys / PowerSurface /
     topk_call / order_key, imported read-only; nothing is written into 192's directory):
         ORACLE    n_pred_p = m * n_real_p, FIT key, m 2.0  ->  d(fid+jac) overall +0.03608 (192 §3.3: +0.0361)
         PREDICTED n_pred_p = m * n_hat_p,  FIT key, m 2.0  ->  d(fid+jac) overall -0.09717 (192 §3.3: -0.0972)
     plus the full m sweep of both rows against work/frontier192.log's printed values (5 decimals).
KA2  513's four H1 DE members recomputed from 435's production-emitter stats with walk513's own functions
     (fid 0.5089291277653674, jac 0.07469853420125332, nmae 0.9427040517722616, reach 0.26925697759272765).

    PYTHONDONTWRITEBYTECODE=1 uv run python experiments/521-strength-nreal/ka521.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E192 = ROOT / "experiments/192-predicted-callset"
sys.path.insert(0, str(E192))
sys.path.insert(0, str(ROOT / "experiments/513-reach-walk-anatomy"))


def ka1():
    import frontier192 as F  # noqa: PLC0415  (read-only import)
    s, keys, absH, sig, coef, _ = F.build_keys()
    nreal = sig.sum(1)
    T = absH.shape[0]
    # 192's n_real_p predictor, verbatim (frontier192.main)
    hn = np.linalg.norm(absH, axis=1)
    feats = np.stack([np.log(hn), np.log(np.abs(absH).sum(1)), np.log((absH > np.quantile(absH, 0.99)).sum(1) + 1)], 1)
    y = np.log(nreal + 1.0)
    pred = np.zeros(T)
    for p in range(T):
        tr = np.setdiff1d(np.arange(T), [p])
        A = np.c_[feats[tr], np.ones(T - 1)]
        beta, *_ = np.linalg.lstsq(A, y[tr], rcond=None)
        pred[p] = np.r_[feats[p], 1.0] @ beta
    ok = nreal > 0
    r = float(np.corrcoef(pred[ok], y[ok])[0, 1])
    npred_hat = np.exp(pred) - 1.0
    f = 1.0
    sf = F.PowerSurface(F.A15, f)
    nr = np.array([int(sf.real_sig[p].sum()) for p in range(T)])
    f0, j0 = sf.fid(), sf.jac()
    # the printed values in 192's own log, f = 1.0 block, FIT key
    log = (E192 / "work/frontier192.log").read_text().split("\n")
    i0 = next(i for i, l in enumerate(log) if l.strip().startswith("f = 1.0"))
    want = {}
    for l in log[i0:i0 + 40]:
        m_ = re.match(r"\s+(ORACLE n_real_p|PREDICTED n_hat_p)\s+m=(\S+)\s+FIT \(cpm,nz,cv2\) LOTO\s+fid (\S+) jac (\S+) "
                      r"reach (\S+)\s+d\(fid\+jac\) overall (\S+)", l)
        if m_:
            want[(m_.group(1), float(m_.group(2)))] = tuple(float(x) for x in m_.group(3, 4, 5, 6))
    rows, ok_all = [], True
    for lab, base in (("ORACLE n_real_p", nr.astype(float)), ("PREDICTED n_hat_p", npred_hat * f)):
        for m in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0):
            call = F.topk_call(keys["FIT (cpm,nz,cv2) LOTO"], np.round(m * base).astype(int))
            fv, jv = sf.fid(call=call), sf.jac(call=call)
            rv = sf.reach(key=F.order_key(call, absH))
            d = ((fv - f0) * F.INV["fid"] + (jv - j0) * F.INV["jac"]) / 6
            w = want[(lab, m)]
            good = all(abs(a - b) <= 6e-6 for a, b in zip((fv, jv, rv, d), w))
            ok_all &= good
            rows.append(dict(row=lab, m=m, fid=fv, jac=jv, reach=rv, d_overall=d, want=w, ok=good))
            print(f"KA1 {lab:18s} m={m:<4} fid {fv:.5f} jac {jv:.5f} reach {rv:.5f} d(fid+jac) {d:+.5f}  "
                  f"192 log {w[3]:+.5f}  {'PASS' if good else 'FAIL'}")
    print(f"KA1 192 predictor LOTO Pearson {r:+.4f} (192: +0.2401)")
    ok_all &= abs(r - 0.2401) < 1e-4
    return ok_all, dict(rows=rows, pearson_192=r, fid0=f0, jac0=j0)


def ka2():
    import walk513 as W  # noqa: PLC0415  (read-only import; its main() is NOT called: it writes into 513's dir)
    d = W.load()
    T, G = d["targets"], d["genes"]
    fid, jac, reach, nm = [], [], [], []
    for k, t in enumerate(T):
        own = G == t
        rs = (d["rq"][k] < W.ALPHA) & ~own
        adj = np.isfinite(d["rl"][k]) & (d["rl"][k] != 0) & ~own
        os_ = (d["oq"][k] < W.ALPHA) & ~own
        n_real = int(rs.sum()); calls = os_ & adj; n_pred = int(calls.sum())
        kk = int((np.sign(d["ol"][k][calls]) == np.sign(d["rl"][k][calls])).sum())
        if n_real > 0 or n_pred > 0:
            fid.append(kk / max(n_pred, n_real))
        u = (rs | os_).sum(); jac.append((rs & os_).sum() / u if u else 1.0)
        gate = rs & np.isfinite(d["rl"][k])
        if gate.sum() >= 10:
            lh = np.nan_to_num(d["ol"][k][gate], nan=0.0, posinf=0.0, neginf=0.0)
            nm.append(np.abs(lh - d["rl"][k][gate]).sum() / np.abs(d["rl"][k][gate]).sum())
        rr = W.reach_row(np.nonzero(rs)[0], W.order_keys(d["oq"][k], d["op"][k], d["ol"][k], G), d["rl"][k], d["ol"][k])
        if rr is not None:
            reach.append(rr[0] / rr[1])
    got = dict(fid=float(np.mean(fid)), jac=float(np.mean(jac)), nmae=float(np.mean(nm)), reach=float(np.mean(reach)))
    want = json.load(open(W.MET))["raw"]
    res = {m: dict(got=got[m], want=want[m], ok=abs(got[m] - want[m]) < 1e-12) for m in got}
    for m, v in res.items():
        print(f"KA2 {m:6s} {v['got']:.16f} vs 435/513 {v['want']:.16f}  {'PASS' if v['ok'] else 'FAIL'}")
    return all(v["ok"] for v in res.values()), res


def main():
    ok1, r1 = ka1()
    ok2, r2 = ka2()
    out = dict(ka1_pass=ok1, ka2_pass=ok2, ka1=r1, ka2=r2)
    (HERE / "work/ka521.json").write_text(json.dumps(out, indent=1, default=float))
    print("KNOWN ANSWERS", "PASS" if ok1 and ok2 else "FAIL")
    if not (ok1 and ok2):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
