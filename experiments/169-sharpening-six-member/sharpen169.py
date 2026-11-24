"""169 — the refined q grid on the two held Orion lines, both renormalisations, mass fix both ways.

Completes three things 168 left open:

  1. the knee is inside `[1.0, 1.3]` and 168 only sampled 1.15 and 1.25;
  2. `nmae` had ONE surface (HCT116).  h1-test's second surface is priced in `de169.py`; this
     script adds the two renormalisations and the mass-fix interaction, which are the two
     things that could make 168's number not transfer;
  3. whether sharpening composes with or double-counts a23's mass fix
     (`VCC_MASS_FIX=minnorm`, shipped, board `pds +0.001737 / mse +0.001525`).

KNOWN-ANSWER ROW (rule 10's corollary): `median` renormalisation with the mass fix ON must
reproduce 168's `work/sharpen.json` -- worst-surface d pds +0.01145 at q=1.15, +0.01799 at
q=1.25, exact d nmae +0.0036 / +0.0327 on HCT116 -- because that is literally the same
computation.  The script asserts it and aborts if it does not.

Deterministic.  No emitter, no cell-eval2, no build, no submission.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, HERE)
from debench import Surface                                     # noqa: E402
from massfix import mass_fix, production_effect                 # noqa: E402
from screen import endpoints, render                            # noqa: E402
from sharp import sharpen, sharpen_lfc                          # noqa: E402

QS = (1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30, 1.50)
EPS = 1e-9
DE_TABLE = {"hct116": os.path.join(REPO, "experiments/095-clean-harness/de_thinned.npz")}
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")

# 168 `work/sharpen.json`, median renormalisation + mass fix, HCT116 / HEK293T.
KNOWN = {1.15: dict(dpds_hct=+0.01145, dpds_hek=+0.01359, dnmae_hct=+0.0036),
         1.25: dict(dpds_hct=+0.01799, dpds_hek=+0.02198, dnmae_hct=+0.0327),
         1.50: dict(dpds_hct=+0.02944, dpds_hek=+0.03950, dnmae_hct=+0.2094)}


def de_endpoints(surface_name, eff, b, ptar):
    """Exact `nmae` + sign accuracy via the cached reference DE table (168's function)."""
    path = DE_TABLE.get(surface_name)
    if path is None or not os.path.exists(path):
        return {}
    s = Surface(path)
    z = np.load(FROZEN % surface_name, allow_pickle=True)
    ctrl_counts = np.load(os.path.join(REPO, "experiments/156-delivered-geometry/work",
                                       f"truth-{surface_name}.npz"),
                          allow_pickle=True)["ctrl_counts"].astype(np.float64)
    ctrl_cpm = 1e6 * ctrl_counts / ctrl_counts.sum()
    lin = np.clip(np.expm1(b[None, :] + eff), 0.0, None)
    p = lin / lin.sum(1, keepdims=True)
    L = np.log2((1e6 * p + EPS) / (ctrl_cpm[None, :] + EPS))
    gpos = {g: i for i, g in enumerate(s.genes)}
    cols = np.asarray([gpos.get(str(g), -1) for g in z["genes"]])
    prow = {pp: i for i, pp in enumerate(s.perts)}
    rows = np.asarray([prow.get(pp, -1) for pp in ptar])
    ok, keep = cols >= 0, rows >= 0
    M = np.full((len(s.perts), len(s.genes)), np.nan)
    M[np.ix_(rows[keep], cols[ok])] = L[np.ix_(np.flatnonzero(keep), np.flatnonzero(ok))]
    M = np.where(np.isfinite(M), M, 0.0)
    s.pred_lfc = M
    return {"nmae_exact": s.nmae(), "sign_acc": s.sign_accuracy()}


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    for surface in ("hct116", "hek293t"):
        D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
        print(f"\n{'=' * 122}\n{surface}: {D.shape[0]} targets, {int(gate.sum())} gated genes, "
              f"{len(mask)} scored genes")
        for mode in ("median", "row", "lfc", "lfcL"):
            for mf in (True, False):
                tag = f"{mode}/{'massfix' if mf else 'raw'}"
                print(f"\n  --- {tag} ---")
                print(f"  {'q':>5} {'pds':>9} {'d pds':>9} {'cos':>9} {'rownorm sd':>11} "
                      f"{'sse*':>10} {'gstar':>8} {'nmae':>9} {'d nmae':>9} {'sign':>7}")
                base, rows = None, []
                for q in QS:
                    if mode in ('lfc', 'lfcL'):
                        h = sharpen_lfc(D, q, b, gate, mode='row',
                                        restore='h' if mode == 'lfc' else 'lfc')
                    else:
                        h = sharpen(D, q, gate, mode=mode)
                    eff = mass_fix(h, b, gate) if mf else h
                    scored = render(eff, b)
                    e = endpoints(scored[:, mask], T[:, mask])
                    de = de_endpoints(surface, eff, b, ptar)
                    nr = np.linalg.norm(scored[:, mask], axis=1)
                    rec = dict(q=q, rn_sd=float(nr.std() / nr.mean()), **e, **de)
                    if base is None:
                        base = rec
                    rows.append(rec)
                    nm = de.get("nmae_exact", float("nan"))
                    print(f"  {q:>5.2f} {e['pds']:9.5f} {e['pds'] - base['pds']:+9.5f} "
                          f"{e['cos']:9.5f} {rec['rn_sd']:11.5f} {e['sse']:10.6f} "
                          f"{e['gstar']:8.4f} {nm:9.4f} "
                          f"{nm - base.get('nmae_exact', nm):+9.4f} "
                          f"{de.get('sign_acc', float('nan')):7.4f}")
                out[f"{surface}|{tag}"] = rows

    with open(os.path.join(WORK, "sharpen169.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)

    # ---------------- known-answer gate against 168 ----------------
    print(f"\n{'=' * 122}\nKNOWN-ANSWER GATE vs experiment 168 (median renorm + mass fix)")
    bad = []
    for q, want in KNOWN.items():
        h = out["hct116|median/massfix"]
        k = out["hek293t|median/massfix"]
        i = QS.index(q)
        got = dict(dpds_hct=h[i]["pds"] - h[0]["pds"], dpds_hek=k[i]["pds"] - k[0]["pds"],
                   dnmae_hct=h[i]["nmae_exact"] - h[0]["nmae_exact"])
        for key, w in want.items():
            g = got[key]
            ok = abs(g - w) <= max(2e-4, 0.02 * abs(w))
            print(f"  q={q:.2f} {key:<10} got {g:+.5f}  expect {w:+.5f}  "
                  f"{'PASS' if ok else 'FAIL'}")
            if not ok:
                bad.append((q, key))
    if bad:
        sys.exit(f"ABORT: 168 reproduction failed for {bad}")

    # ---------------- worst-surface tables ----------------
    for mode in ("median", "row", "lfc", "lfcL"):
        for mf in (True, False):
            tag = f"{mode}/{'massfix' if mf else 'raw'}"
            h, k = out[f"hct116|{tag}"], out[f"hek293t|{tag}"]
            print(f"\n{'=' * 122}\nWORST-SURFACE, {tag}")
            print(f"  {'q':>5} {'worst d pds':>12} {'hct116':>10} {'hek293t':>10} "
                  f"{'sse* hct':>10} {'sse* hek':>10} {'d nmae hct':>11} {'rn sd hct':>10}")
            for i, q in enumerate(QS):
                dh, dk = h[i]["pds"] - h[0]["pds"], k[i]["pds"] - k[0]["pds"]
                print(f"  {q:>5.2f} {min(dh, dk):+12.5f} {dh:+10.5f} {dk:+10.5f} "
                      f"{h[i]['sse']:10.6f} {k[i]['sse']:10.6f} "
                      f"{h[i]['nmae_exact'] - h[0]['nmae_exact']:+11.4f} "
                      f"{h[i]['rn_sd']:10.5f}")

    # ---------------- does sharpening double-count the mass fix? ----------------
    print(f"\n{'=' * 122}\nCOMPOSITION WITH a23's MASS FIX (board: pds +0.001737)")
    print("  If sharpening merely re-did what the mass fix does, the mass fix's own gain would")
    print("  SHRINK as q rises.  d(massfix) at each q, per surface:")
    for mode in ("median", "row", "lfc", "lfcL"):
        print(f"  {'q':>5}  " + "  ".join(
            f"{s + '/' + mode:>22}" for s in ("hct116", "hek293t")))
        for i, q in enumerate(QS):
            cells = []
            for s in ("hct116", "hek293t"):
                a = out[f"{s}|{mode}/massfix"][i]["pds"]
                r = out[f"{s}|{mode}/raw"][i]["pds"]
                cells.append(f"{a - r:+22.5f}")
            print(f"  {q:>5.2f}  " + "  ".join(cells))


if __name__ == "__main__":
    sys.exit(main())
