"""Is HepG2 a good transfer partner, or are common-essential targets easy?

`independence.py` reads `nadig_hepg2 x orion_hct116` at resid same +0.1091 against
`orion_hct116 x orion_hek293t` at +0.0291 -- a 3.7x difference that would say the new line is a
far better partner than either shipped Orion line is for the other.

**That comparison is confounded and must not be quoted as it stands.** The HepG2 pairs can only
be evaluated on HepG2's own targets, which are the Nadig DepMap common-essential set; the
Orion-Orion pair is evaluated on 3,000 targets drawn from ~17,700 mostly non-essential ones.
159 section 5.3 already records that common-essential genes transfer best -- their responses are
large, stereotyped growth-arrest programmes. So the difference could be entirely about WHICH
targets are scored.

The control: re-run `orion_hct116 x orion_hek293t` restricted to **exactly the target set the
HepG2 pair uses**, on the same gene axis rule. If the Orion pair jumps to HepG2-like values, the
target set explains it and HepG2 is unremarkable. If it stays near +0.029, the line is.

Writes `work/confound.json`.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

import independence as I

HERE = Path(__file__).resolve().parent


def restricted(dA, dB, targets, rng):
    """159's pair statistic on a FIXED target list (no subsampling)."""
    pa = {p: i for i, p in enumerate(dA["perts"])}
    pb = {p: i for i, p in enumerate(dB["perts"])}
    sh = sorted(set(targets) & set(pa) & set(pb))
    ga = {g: i for i, g in enumerate(dA["genes"])}
    gb = {g: i for i, g in enumerate(dB["genes"])}
    shg = [g for g in sorted(set(ga) & set(gb)) if g not in set(sh)]
    A = dA["delta"][[pa[p] for p in sh]][:, [ga[g] for g in shg]].astype(np.float64)
    B = dB["delta"][[pb[p] for p in sh]][:, [gb[g] for g in shg]].astype(np.float64)
    out = {"n_perts": len(sh), "n_genes": len(shg)}
    for space in ("full", "resid"):
        X = A - A.mean(0, keepdims=True) if space == "resid" else A
        Y = B - B.mean(0, keepdims=True) if space == "resid" else B
        a, b = I.unit(X), I.unit(Y)
        S = a @ b.T
        same, diff = float(np.diag(S).mean()), I.offdiag_mean(S)
        out[space] = {"same": same, "diff": diff,
                      "ratio": same / diff if abs(diff) > 1e-6 else float("nan")}
        if space == "resid":
            out[space]["same_permuted"] = float(
                (a * b[rng.permutation(len(sh))]).sum(1).mean())
    return out


def main() -> None:
    rng = np.random.default_rng(I.SEED)
    L = {k: I.load_donor(rel) for k, (rel, _) in I.DONORS.items()}
    L = {k: v for k, v in L.items() if v is not None}

    hep = set(L["nadig_hepg2"]["perts"])
    hct, hek = L["orion_hct116"], L["orion_hek293t"]

    # the exact target set the HepG2 x Orion pairs are scored on
    tgt = sorted(hep & set(hct["perts"]) & set(hek["perts"]))
    print(f"HepG2 targets present in BOTH Orion lines: {len(tgt)}\n")

    rows = {}
    rows["orion_hct116|orion_hek293t ON HEPG2 TARGETS"] = restricted(hct, hek, tgt, rng)
    rows["nadig_hepg2|orion_hct116"] = restricted(L["nadig_hepg2"], hct, tgt, rng)
    rows["nadig_hepg2|orion_hek293t"] = restricted(L["nadig_hepg2"], hek, tgt, rng)
    # and the same Orion pair on a size-matched random draw of its OWN target universe,
    # so target COUNT is held fixed and only target IDENTITY changes
    own = sorted(set(hct["perts"]) & set(hek["perts"]) - hep)
    draw = [own[i] for i in sorted(rng.choice(len(own), len(tgt), replace=False))]
    rows[f"orion_hct116|orion_hek293t ON {len(tgt)} NON-HEPG2 TARGETS"] = restricted(
        hct, hek, draw, rng)

    print(f"{'pair / target set':<52}{'n_pert':>7}{'same':>9}{'diff':>9}"
          f"{'RATIO':>8}{'resid same':>12}{'resid null':>12}")
    for k, r in rows.items():
        print(f"{k:<52}{r['n_perts']:7d}{r['full']['same']:+9.4f}{r['full']['diff']:+9.4f}"
              f"{r['full']['ratio']:8.2f}{r['resid']['same']:+12.4f}"
              f"{r['resid']['same_permuted']:+12.4f}")

    a = rows["orion_hct116|orion_hek293t ON HEPG2 TARGETS"]["resid"]["same"]
    b = rows[f"orion_hct116|orion_hek293t ON {len(tgt)} NON-HEPG2 TARGETS"]["resid"]["same"]
    c = rows["nadig_hepg2|orion_hct116"]["resid"]["same"]
    print(f"\n  target-set effect  (Orion pair, essential vs not): {a:+.4f} vs {b:+.4f} "
          f"= x{a/b:.2f}")
    print(f"  line effect        (HepG2 vs Orion, SAME targets):  {c:+.4f} vs {a:+.4f} "
          f"= x{c/a:.2f}")

    (HERE / "work/confound.json").write_text(json.dumps(rows, indent=2))
    print("\nwrote work/confound.json")


if __name__ == "__main__":
    main()
