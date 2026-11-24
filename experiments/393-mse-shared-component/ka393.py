#!/usr/bin/env python3
"""393 KA reproduction (own code, writes only into this experiment's work/).

KA-i: 187's C_p-free LS solve on a8+a10+a11 -> X = 0.019240, sumD = 0.229141.
KA-ii: 211's per-context claim saturation q_c(a15) = 1.000000/0.999904/1.000000,
  recomputed from the stored per-row C_hat values in chat211.json with own code.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", str(Path(__file__).resolve().parent / "work" / "numba-cache"))
os.environ.setdefault("NUMBA_NUM_THREADS", "4")
os.environ.setdefault("OMP_NUM_THREADS", "4")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "experiments/096-mse-floor"))

CTX = ("A", "B", "C")
ALPHA_PROD = 1.0796
BOARD_RAW = {"a8": 1.362666960949, "a10": 1.3524148, "a11": 0.9099893}
ALPHA = {"a8": 1.0796, "a10": 1.0796, "a11": 0.838}


def emin(chat, mu):
    return np.minimum(np.asarray(chat, float), mu)


def claim_of(chat_rows, w_row, cp):
    return float((emin(chat_rows, cp) * w_row).sum())


def main() -> None:
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))

    # ---- KA-i ----
    from mech_a12_audit import cp_proxies, measure
    W096 = ROOT / "experiments/096-mse-floor/work"
    cps = cp_proxies()["ntc"]
    cp16 = {c: cps[c][0] for c in CTX}
    arms = {}
    for arm in ("a8", "a10", "a11"):
        rows = measure(str(W096 / f"{arm}-stats.npz"), cp16)
        per = {}
        for r in rows:
            w_row = r["K"] / (r["Cp"] * r["n"])
            per[r["ctx"]] = dict(A=r["A"], B=r["B"], Kc=r["Kc"],
                                 Chat=np.asarray(r["Chat"]), w=w_row, Cp=r["Cp"])
        arms[arm] = per
    terms, ss, ys = [], [], []
    for a in ("a8", "a10", "a11"):
        t = {}
        for c in CTX:
            p = arms[a][c]
            w = np.full(p["Chat"].size, p["w"])
            t[c] = dict(A=p["A"], B=p["B"], Kc=p["Kc"],
                        claim=claim_of(p["Chat"], w, p["Cp"]))
        terms.append(t)
        ss.append(ALPHA[a] / ALPHA_PROD)
        ys.append(BOARD_RAW[a])
    M = np.array([[-2.0 * s, -(y - 1.0)] for s, y in zip(ss, ys)])
    b = np.array([-np.mean([t[c]["A"] + t[c]["Kc"] - min(t[c]["claim"], t[c]["B"]) for c in CTX])
                  for t in terms])
    sol, *_ = np.linalg.lstsq(M, b, rcond=None)
    X0, sumD0 = float(sol[0]), float(sol[1])
    ref187 = json.loads((ROOT / "experiments/211-amplitude-reopened/work/mse187.json").read_text())
    P(f"  full precision: X0={X0!r} (187 rerun {ref187['X0']!r}, "
      f"|d|={abs(X0 - ref187['X0']):.2e}); sumD0={sumD0!r} (187 {ref187['sumD0']!r}, "
      f"|d|={abs(sumD0 - ref187['sumD0']):.2e})")
    ok_x = abs(X0 - 0.019240) < 2e-5
    ok_s = abs(sumD0 - 0.229141) < 2e-4
    P(f"[KA-i] X = {X0:.6f} (187: 0.019240) {'PASS' if ok_x else 'FAIL'}; "
      f"sumD = {sumD0:.6f} (187: 0.229141) {'PASS' if ok_s else 'FAIL'}")
    # 211 prints residuals at 2 s.f. (+6.0e-5/-6.1e-5/+0.2e-5); full-precision check
    # below uses 6.5e-5 and the equality to 211's printed digits is the real gate.
    for t, s, y in zip(terms, ss, ys):
        num = np.mean([t[c]["A"] + t[c]["Kc"] - min(t[c]["claim"], t[c]["B"]) for c in CTX])
        r = 1.0 + (num - 2.0 * s * X0) / sumD0
        P(f"       residual {r - y:+.2e} (tol 6.5e-5) {'PASS' if abs(r - y) <= 6.5e-5 else 'FAIL'}")

    # ---- KA-ii: q_c from stored per-row C_hat ----
    ch = json.loads((ROOT / "experiments/211-amplitude-reopened/work/chat211.json").read_text())
    cp = {c: ch["C_p"][c] for c in CTX}
    P("\n[KA-ii] q_c = mean_p min(C_hat_p, C_p)/C_p recomputed from stored rows:")
    ok_all = True
    q = {}
    for c in CTX:
        block = ch["chat"][f"{c}|quasi"]
        # 211's published q is the MAXT24 sample: sel = linspace(0, maxkey, 24)
        # (chat211.py; the stored dict is the 8+24 union, the record quotes MAXT24)
        maxkey = max(int(k) for k in block.keys())
        sel24 = np.linspace(0, maxkey, 24).astype(int)
        P(f"  ctx {c}: stored union rows={len(block)} maxkey={maxkey} "
          f"sel24={sel24[:4].tolist()}...{sel24[-2:].tolist()}")
        vals = np.array([block[str(i)] for i in sel24])
        q[c] = float(np.minimum(vals, cp[c]).mean() / cp[c])
        P(f"  ctx {c}: n={len(vals)} C_p={cp[c]:.4f} q={q[c]:.6f}")
    exp = {"A": 1.000000, "B": 0.999904, "C": 1.000000}
    for c in CTX:
        ok = abs(q[c] - exp[c]) < 2e-6
        ok_all &= ok
        P(f"  ctx {c}: recomputed {q[c]:.6f} vs 211 {exp[c]:.6f} {'PASS' if ok else 'FAIL'}")

    res = dict(X=X0, sumD=sumD0, ka_i_pass=bool(ok_x and ok_s), q=q, ka_ii_pass=bool(ok_all))
    (WORK / "ka393.json").write_text(json.dumps(res, indent=1))
    (WORK / "ka393.log").write_text("\n".join(lines))
    P("\nOVERALL " + ("PASS" if (ok_x and ok_s and ok_all) else "FAIL"))
    if not (ok_x and ok_s and ok_all):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
