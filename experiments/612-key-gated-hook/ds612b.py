"""612 ds-b: ds612's forms on the four folds at 521's REGIMES.

At regime f, a target's pool is the truth's top round(f * n_real) gate genes by the truth's own p (stable sort; own target
already off the gate), and a row is scored only with >= 10 pool genes (cell-eval2's nmae rule). f 1 is ds612's gate.
These are the regimes in which H1's real emitter showed the plain hook's STRONG-half shrink helping (612/diag612b: +0.0011
raw nmae at f 0.35, seed 0). ds612's thinned H1 surfaces thin the CELLS instead, which asks a different question.

IDw_C4s - kgID is exactly the value of the strong-half CLEAN4 shrink: the two forms differ only on strong-half cells.
Analysis only; < 3 GB; seconds. Output: work/ds612b.out, work/ds612b.json.
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import ds612 as D  # noqa: E402  (puts 480/537/622 on the path)
from common480 import boot_delta  # noqa: E402

SURFS = ["h1", "orion_hct116", "orion_hek293t", "k562_gwps"]
REGIMES = (1.0, 0.5, 0.35)
FORMS = ["plain", "kgC4", "kgID", "IDw_C4s", "plainID", "kgID_nf", "kgID+mshC4s_0.5"]


def pvals(name, T):
    if name == "h1":
        p = np.load(D.ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=True)["p"]
    else:
        p = np.load(D.ROOT / f"experiments/480-nmae-instrument/work/gate-{name}.npz", allow_pickle=True)["p"]
    p = np.asarray(p, float)
    assert p.shape == T.shape, (name, p.shape, T.shape)
    return p


def pool_T(T, p, f, min_gate=10):
    """521's PowerSurface pool on the fold truth: top round(f n) gate genes by p; rows with < min_gate dropped."""
    if f == 1.0:
        return T
    out = np.full(T.shape, np.nan)
    for i in range(T.shape[0]):
        g = np.isfinite(T[i])
        k = int(round(f * int(g.sum())))
        if k < min_gate:
            continue
        idx = np.argsort(np.where(g & np.isfinite(p[i]), p[i], np.inf), kind="stable")[:k]
        out[i, idx] = T[i, idx]
    return out


def main() -> None:
    BETA, prof = D.setup()
    lines, out = [], {}

    def say(s=""):
        print(s)
        lines.append(s)

    for name in SURFS:
        X = D.surface_forms(name, BETA, prof)
        T, P = X["T"], X["P"]
        p = pvals(name, T)
        out[name] = {}
        say(f"\n== {name}: Δnmae (+ = better) [95% target bootstrap]; columns = 521 regimes")
        hdr = f"  {'form':16s}" + "".join(f" | f {f:<4g} vs nohook            vs kgID              " for f in REGIMES)
        say(hdr)
        res = {}
        for f in REGIMES:
            Tf = pool_T(T, p, f)
            rows = {k: D.nmae_cells(v, Tf)[0] for k, v in P.items()}
            n = int(np.isfinite(rows["nohook"]).sum())
            res[f] = {"n": n, "vs_nohook": {k: boot_delta(rows[k], rows["nohook"]) for k in FORMS},
                      "vs_kgID": {k: boot_delta(rows[k], rows["kgID"]) for k in FORMS if k != "kgID"}}
        say("  scored targets: " + ", ".join(f"f {f:g}: {res[f]['n']}" for f in REGIMES))
        for k in FORMS:
            cells = []
            for f in REGIMES:
                a = res[f]["vs_nohook"][k]
                b = res[f]["vs_kgID"].get(k)
                cells.append(f"{a[0]:+.5f} [{a[1]:+.4f},{a[2]:+.4f}] "
                             + (f"{b[0]:+.5f} [{b[1]:+.4f},{b[2]:+.4f}]" if b else " " * 27))
            say(f"  {k:16s} " + " | ".join(cells))
        s = [res[f]["vs_kgID"]["IDw_C4s"] for f in REGIMES]
        say("  strong-half CLEAN4 shrink value (IDw_C4s - kgID): "
            + ", ".join(f"f {f:g} {v[0]:+.5f} [{v[1]:+.5f},{v[2]:+.5f}]" for f, v in zip(REGIMES, s)))
        out[name] = {str(f): res[f] for f in REGIMES}
    (HERE / "work/ds612b.out").write_text("\n".join(lines) + "\n")
    (HERE / "work/ds612b.json").write_text(json.dumps({"BETA": BETA, "surfaces": out}, indent=1))
    say("\nwrote work/ds612b.out, work/ds612b.json")


if __name__ == "__main__":
    main()
