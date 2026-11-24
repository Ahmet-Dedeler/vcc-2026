"""534 Part C -- efr's re-level factors from the FOLD-ROUTE energies (438 section 2.5: efr has no 3-donor operator). Light.

price438.fold_energies + foldlevel438, verbatim arithmetic: each design's realised render energy averaged over efr's
three cancer folds (comp534c's first, unlevelled pass), anchored at 417a's measured energy through A/B/C's own fold-route
417a energy (438's stored abc rows, read only); ntclevel's rule s = sqrt(rho* claim / B), rho* at sd 0.01.
Claims: ntclevel417's f = 1 claims from the efr bundle's NTC cells (comp534c cached them).

KNOWN ANSWERS: the claims == 438's efr_s1 claims; the factors of h1c1|K3 and h2c1|V1_5_3+sh1.5 == 438's efr
level438.json (relative |d| <= 1e-12).

    uv run python foldlevel534c.py efr        -> work/partC/efr/level534c.json
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401
import json  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

import designs534c as D  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
W438 = ROOT / "experiments/438-october-stack-repricing/work"
W420 = ROOT / "experiments/420-bulk-stack/work"
SLOTS = ("A", "B", "C")
REF = "h1c1|K3+sh1.25"


def fold_energies(panel, sfx=""):
    """{design: AB (anchored, unlevelled) and, where levelled rows exist, AB_level} over the three cancer folds."""
    a420 = json.loads((W420 / "abc420.json").read_text())
    M0 = {c: a420["measured"]["417a"][c] for c in SLOTS}
    R0 = {f: json.loads((W438 / "abc" / f"comp438-{f}-held.json").read_text())["arms"] for f in D.CANCER}
    R = {f: json.loads((HERE / panel / f"comp534c{sfx}-{f}.json").read_text())["arms"] for f in D.CANCER}
    mE = lambda RR, n: {c: {x: float(np.mean([RR[f][n]["realised"]["16"]["AB"][c][x] for f in D.CANCER])) for x in ("A", "B")}
                        for c in SLOTS}
    E0 = mE(R0, REF)
    out = {}
    for n in R[D.CANCER[0]]:
        if n.endswith(":L") or not R[D.CANCER[0]][n].get("realised"):
            continue
        e = mE(R, n)
        d = dict(AB={c: {x: M0[c][x] * e[c][x] / E0[c][x] for x in ("A", "B")} for c in SLOTS})
        if all(f"{n}:L" in R[f] for f in D.CANCER):
            el = mE(R, f"{n}:L")
            d["AB_level"] = {c: {x: M0[c][x] * el[c][x] / E0[c][x] for x in ("A", "B")} for c in SLOTS}
        out[n] = d
    return out


def main(panel, sfx=""):
    assert panel in D.FOLD_ROUTE
    from scipy.stats import norm
    claims = json.loads((HERE / panel / f"comp534c{sfx}-k562_gwps.json").read_text())["claims"]
    c438 = json.loads((W438 / "efr_s1" / "full438.json").read_text())["claims"]
    dcl = max(abs(claims[c] - c438[c]) for c in SLOTS)
    assert dcl == 0.0, f"KA-CLAIMS: {dcl}"
    sd = 0.01
    rho_t = float(np.exp(2.0 * sd * norm.ppf(0.31 / (0.31 + 1.78))))
    E = fold_energies(panel, sfx)
    fac = {n: {c: float(np.sqrt(rho_t * claims[c] / e["AB"][c]["B"])) for c in SLOTS} for n, e in E.items()}
    l438 = json.loads((W438 / panel / "level438.json").read_text())
    worst = {n: max(abs(fac[n][c] / l438["factors"][n][c] - 1) for c in SLOTS) for n in ("h1c1|K3", D.NEAREST_438)}
    assert "h1c1|none" not in fac and "h1c1|none" not in l438["factors"]      # design-only: no fold-route energy (as 438)
    for n, v in worst.items():
        print(f"KA factor {n:22s} vs 438's efr level438.json: rel |d| {v:.1e}")
        assert v <= 1e-12, f"KA FAILED {n}: {v}"
    out = dict(panel=panel, route="fold", sd=sd, rho_target=rho_t, claims=claims, claims_vs_438_efr_s1=dcl, factors=fac,
               ka=worst, energies=E)
    if sfx:                                                   # KA-PRIMARY: the primary factors, exactly
        pl = json.loads((HERE / panel / "level534c.json").read_text())["factors"]
        assert all(fac[n] == pl[n] for n in pl), "KA-PRIMARY factors"
        out["ka"]["primary"] = f"{len(pl)} designs' factors == the primary run's"
    (HERE / panel / f"level534c{sfx}.json").write_text(json.dumps(out, indent=1))
    print(f"{panel}: {len(fac)} designs; B0 factors " + "/".join(f"{fac['B0'][c]:.4f}" for c in SLOTS))


if __name__ == "__main__":
    main(sys.argv[1], "pk" if "--pk" in sys.argv else "")
