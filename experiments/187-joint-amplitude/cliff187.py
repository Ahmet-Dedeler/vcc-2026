"""187 — where are the per-context rho cliffs, really?  Joint constraint from the two board
amplitude deltas that straddle them.

Two board deltas are amplitude-only moves whose sign and size depend on the cliff positions:
    a12 -> a14 : alpha 1.0108 -> 0.982 on the a12 profile   (board d raw -0.018574)
    a23 -> knee: A .98234 / B 1.01291 / C 1.04786 on a23    (board d scaled -0.002259)

mse187.py showed that no single (level bias, CV) hypothesis common to all three contexts fits
both.  This script lets the C_p error be PER CONTEXT — a multiplicative factor f_c on the NTC
proxy (cv = 0), or a CV — and maps the region of (f_B, f_C) that reproduces both deltas.  f_A
is pinned by a14 (which sits on A's proxy cliff to 0.02%) and is scanned separately.

Reads work/mse187.json (per-context terms for the a12 family and for a23) — no stats reload.

    uv run python cliff187.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from mse187 import (ALPHA, ALPHA_PROD, BOARD_RAW, CTX, H1_CV, arm_terms, claim_of, raw_of,
                    score, solve_X, terms_180)

HERE = Path(__file__).resolve().parent
J = json.loads((HERE / "work/mse187.json").read_text())
X0, sumD0 = J["X0"], J["sumD0"]
BOARD_D14 = BOARD_RAW["a14"] - BOARD_RAW["a12"]
BOARD_DKNEE = -0.002259


def per12():
    per = {}
    for c in CTX:
        d = J["arms"]["a12"][c]
        per[c] = dict(A=d["A"], B=d["B"], Kc=d["Kc"], Cp=d["Cp"], w=d["w"], n=d["n"],
                      Chat=np.asarray(d["Chat"]))
    return per


def per23():
    per = {}
    for c in CTX:
        d = J["a23"][c]
        per[c] = {k: (np.asarray(v) if isinstance(v, list) else v) for k, v in d.items()}
    return per


P12, P23 = per12(), per23()
S12, S14, F14 = ALPHA["a12"] / ALPHA_PROD, ALPHA["a14"] / ALPHA_PROD, ALPHA["a14"] / ALPHA["a12"]
KAPPA = {c: P23[c]["kappa"] for c in CTX}


def terms_pc(per, f, cv, s=1.0, key="Cp"):
    """per-context bias f[c] and cv[c] on the proxy C_p (a12 family: `Cp`; a23: `Cp46`)."""
    out = {}
    for c in CTX:
        p = per[c]
        chat = p["Chat"]
        w = np.full(chat.size, (p["w"] if "w" in p else p["weight_sum"] / chat.size))
        out[c] = dict(A=s * s * p["A"], B=s * s * p["B"], Kc=p["Kc"],
                      claim=claim_of(chat, w, p[key], f[c], cv[c]))
    return out


def terms_23(f, cv, which="before"):
    out = {}
    for c in CTX:
        p = P23[c]
        chat = p["Chat"] if which == "before" else p["Chat_after"]
        w = np.full(chat.size, p["weight_sum"] / chat.size)
        out[c] = dict(A=p["A"] if which == "before" else p["A_after"],
                      B=p["B"] if which == "before" else p["B_after"], Kc=p["Kc"],
                      claim=claim_of(chat, w, p["Cp46"], f[c], cv[c]))
    return out


def deltas(f, cv):
    """(pred a12->a14 raw delta, pred a23->knee scaled delta, a23 cliffs in alpha, X23)."""
    r12, _ = raw_of(terms_pc(P12, f, cv), X0, sumD0, S12)
    r14, _ = raw_of(terms_pc(P12, f, cv, s=F14), X0, sumD0, S14)
    t23 = terms_23(f, cv)
    tk = terms_23(f, cv, which="after")
    X23 = solve_X(t23, BOARD_RAW["a23"], sumD0)
    _, p0 = raw_of(t23, X23, sumD0, 1.0)
    _, p1 = raw_of(tk, X23, sumD0, 1.0, kappa=KAPPA)
    dk = score(list(p1.values())) - score(list(p0.values()))
    cliffs = {c: 0.982 * float(np.sqrt(t23[c]["claim"] / t23[c]["B"])) for c in CTX}
    cliffs12 = {c: ALPHA["a12"] * float(np.sqrt(terms_pc(P12, f, cv)[c]["claim"] / P12[c]["B"])) for c in CTX}
    return r14 - r12, dk, cliffs, cliffs12, X23, {c: p1[c] - p0[c] for c in CTX}


def main():
    lines = []
    P = lambda *a: (print(*a, flush=True), lines.append(" ".join(str(x) for x in a)))
    one = {c: 1.0 for c in CTX}
    zero = {c: 0.0 for c in CTX}

    P("=== single-context factor scans (cv = 0 elsewhere): error = pred - board ===")
    P(f"  {'ctx':>4}{'f_c':>8}{'e(a12->a14)':>13}{'e(knee)':>10}{'a23 cliff A/B/C':>28}{'a12 cliff A/B/C':>28}")
    scans = {}
    for c in CTX:
        scans[c] = []
        for f in np.arange(0.950, 1.0301, 0.005):
            fc = dict(one); fc[c] = float(f)
            d14, dk, cl, cl12, X23, per = deltas(fc, zero)
            scans[c].append((float(f), d14 - BOARD_D14, dk - BOARD_DKNEE))
            P(f"  {c:>4}{f:>8.3f}{d14-BOARD_D14:>+13.5f}{dk-BOARD_DKNEE:>+10.5f}   "
              + "/".join(f"{cl[x]:.4f}" for x in CTX) + "   " + "/".join(f"{cl12[x]:.4f}" for x in CTX))

    P("\n=== per-context CV scans (bias 1 everywhere) ===")
    P(f"  {'ctx':>4}{'cv_c':>8}{'e(a12->a14)':>13}{'e(knee)':>10}{'a23 cliff A/B/C':>28}")
    for c in CTX:
        for cv in (0.0, 0.02, 0.04, H1_CV, 0.07, 0.10):
            cvc = dict(zero); cvc[c] = float(cv)
            d14, dk, cl, cl12, X23, per = deltas(one, cvc)
            P(f"  {c:>4}{cv:>8.3f}{d14-BOARD_D14:>+13.5f}{dk-BOARD_DKNEE:>+10.5f}   "
              + "/".join(f"{cl[x]:.4f}" for x in CTX))

    # 2-D map over (f_B, f_C), f_A = 1, cv = 0: where do BOTH deltas fit to +-0.001?
    P("\n=== (f_B, f_C) map, f_A = 1, cv = 0: cells with |e14| < 0.001 raw AND |e_knee| < 0.001 scaled ===")
    fits = []
    grid = np.arange(0.960, 1.0201, 0.0025)
    for fB in grid:
        for fC in grid:
            f = {"A": 1.0, "B": float(fB), "C": float(fC)}
            d14, dk, cl, cl12, X23, per = deltas(f, zero)
            e14, ek = d14 - BOARD_D14, dk - BOARD_DKNEE
            if abs(e14) < 0.001 and abs(ek) < 0.001:
                fits.append(dict(fB=float(fB), fC=float(fC), e14=e14, ek=ek, cliffs=cl, cliffs12=cl12,
                                 per_ctx=per))
    for r in fits:
        P(f"  f_B {r['fB']:.4f} f_C {r['fC']:.4f}  e14 {r['e14']:+.5f} e_knee {r['ek']:+.5f}  "
          f"a23 cliffs " + "/".join(f"{r['cliffs'][x]:.4f}" for x in CTX) +
          "  knee per-ctx d raw " + " ".join(f"{x} {r['per_ctx'][x]:+.5f}" for x in CTX))
    if fits:
        fB = np.array([r["fB"] for r in fits]); fC = np.array([r["fC"] for r in fits])
        P(f"  region: f_B in [{fB.min():.4f}, {fB.max():.4f}], f_C in [{fC.min():.4f}, {fC.max():.4f}], "
          f"{len(fits)} cells of {grid.size**2}")

    # f_A from a14 alone (a14 sits on A's proxy cliff): scan f_A with f_B = f_C = 1
    P("\n=== f_A scan (a14 sits on A's proxy cliff; only a14's level and the delta respond) ===")
    for fA in np.arange(0.97, 1.0301, 0.005):
        f = {"A": float(fA), "B": 1.0, "C": 1.0}
        r12, _ = raw_of(terms_pc(P12, f, zero), X0, sumD0, S12)
        r14, _ = raw_of(terms_pc(P12, f, zero, s=F14), X0, sumD0, S14)
        P(f"  f_A {fA:.3f}  a12 {r12-BOARD_RAW['a12']:+.5f}  a14 {r14-BOARD_RAW['a14']:+.5f}  delta err {(r14-r12)-BOARD_D14:+.5f}")

    (HERE / "work/cliff187.json").write_text(json.dumps(dict(scans=scans, fits=fits), indent=1, default=float))
    (HERE / "work/cliff187.log").write_text("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
