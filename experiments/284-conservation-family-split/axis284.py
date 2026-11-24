"""284 — POST-HOC control the mandate does not ask for but the result needs: is the reproducible
remainder the perturbed gene's OWN knockdown?

279's statistic is the cosine of two lines' residual response rows on a common 2,739-gene axis
(`lib276.load_source_path`); it does not drop the perturbed gene's own column. A gene that sits ON
that axis therefore contributes its own on-target knockdown to every pair's cosine — a reproducible,
gene-specific, biologically uninteresting component. Genes NOT on the axis cannot: for them the
cosine is downstream response only.

  PYTHONDONTWRITEBYTECODE=1 uv run python axis284.py

Outputs: work/axis284.json, work/axis284.log
"""
from __future__ import annotations

import hashlib
import json
import sys
import time

import numpy as np

import lib284 as L

T0 = time.time()
LOG: list[str] = []


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def main() -> None:
    h = hashlib.sha256((L.HERE / "PREREGISTRATION.md").read_bytes()).hexdigest()
    if h != (L.WORK / "PREREGISTRATION.sha256").read_text().split()[0]:
        sys.exit("preregistration digest mismatch — aborting")
    held = L.load_held()
    sp = np.load(L.WORK / "splits284.npz", allow_pickle=True)
    comp = sp["comp"]
    rem = np.load(L.WORK / "rem284.npz", allow_pickle=True)
    C1, C2, r1, r2, ok = rem["C1"], rem["C2"], rem["r1"], rem["r2"], rem["ok"]
    z = np.load(L.HELD279, allow_pickle=True)
    axis = set(np.asarray(z["axis"]).astype(str).tolist())
    gene = held["gene"].to_numpy().astype(str)
    on = np.array([g in axis for g in gene])
    P(f"axis {len(axis)} genes; held genes ON the axis {int(on.sum())} "
      f"({on.mean():.1%}); of the {int(ok.sum())} remainder-test genes, "
      f"{int((ok & on).sum())} are on the axis and {int((ok & ~on).sum())} are not")

    out = {"n_axis": len(axis), "n_held_on_axis": int(on.sum()), "strata": {}}
    for lab, m in (("ON the axis (self-knockdown possible)", ok & on),
                   ("OFF the axis (downstream response only)", ok & ~on)):
        b = L.boot_rho(r1[m], r2[m], comp[m], np.random.default_rng(L.SEED + 801))
        c = L.boot_rho(C1[m], C2[m], comp[m], np.random.default_rng(L.SEED + 802))
        out["strata"][lab] = dict(n=int(m.sum()), remainder=b, ceiling=c,
                                  median_C1=float(np.nanmedian(C1[m])),
                                  median_C2=float(np.nanmedian(C2[m])))
        P(f"  {lab:42s} n {int(m.sum()):5d}  remainder {b['point']:+.4f} "
          f"[{b['lo']:+.4f}, {b['hi']:+.4f}]  ceiling {c['point']:+.4f} "
          f"[{c['lo']:+.4f}, {c['hi']:+.4f}]  median C {np.nanmedian(C1[m]):+.4f} / "
          f"{np.nanmedian(C2[m]):+.4f}")
    out["log"] = LOG
    (L.WORK / "axis284.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "axis284.log").write_text("\n".join(LOG))
    P("wrote work/axis284.json")


if __name__ == "__main__":
    main()
