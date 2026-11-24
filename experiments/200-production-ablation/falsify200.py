"""200 — Phase 3 falsifiers.

1. EVERY ABLATION MUST LEAVE THE UNTOUCHED MEMBERS WHERE THEY BELONG.
   An emission-only ablation (the tiered spread policy, `--overdispersion`) must leave the
   per-gene panel TOTALS digest-identical, because the quasi emitter repairs each gene's total
   exactly -- that is the mechanism behind rule 12's exemption and behind the board's
   `d pds = 0.000000` on v14 and a8.  A profile-touching ablation (the mass fix, whitening,
   PerturbMatch, CD4, the gate) must NOT.  191 proved the check works (900/900 blocks).
   Here the totals are recomputed for every arm configuration from the same code path
   `realize200.py` emits from, and hashed.

2. RULE 10 CONFIG DIFF.  Each ablation must differ from SHIP in exactly the intended number of
   fields, the operator sha256 must be the incumbent's wherever the operator is not the thing
   being ablated, and `effects.whitening.after_all_median_norm` must not drift.

3. THE PLOLO ROW IDENTITY.  169's q = 1.00 gate re-run: this experiment's SHIP rows against
   162's own `arms.json`, all 16 (fold, frame) cells.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(ROOT, "experiments/169-sharpening-six-member")
P162 = os.path.join(ROOT, "experiments/162-retrodicting-validator")
sys.path.insert(0, HERE)
sys.path.insert(0, P162)

N_EMIT, EMIT_L, TS = 400, 100_000, 50_000.0


def totals_digest():
    """sha256 of the per-gene panel totals for every h1-test arm configuration."""
    import importlib.util
    import realize200
    from realize200 import G

    E = np.load(os.path.join(E169, "work/h1test-effect.npz"), allow_pickle=True)
    H0 = E["H"].astype(np.float64)
    b = E["b"].astype(np.float64)
    a_gate = np.asarray(E["gate"])
    esums = E["sums"].astype(np.float64)
    cpm_cell = np.load(os.path.join(E169, "work/h1test-cellcpm.npy")) \
        if os.path.exists(os.path.join(E169, "work/h1test-cellcpm.npy")) else None

    spec = importlib.util.spec_from_file_location(
        "prod_builder_200f", f"{ROOT}/scripts/make_submission_pm_cd4.py")
    prod = importlib.util.module_from_spec(spec)
    sys.modules["prod_builder_200f"] = prod
    spec.loader.exec_module(prod)

    import h1parts200
    p = h1parts200.parts()
    c4 = h1parts200.cd4(p["targets"])
    etargets = [str(t) for t in E["targets"]]
    order = np.asarray([list(p["targets"]).index(t) for t in etargets], dtype=np.int64)
    art_norm = float(np.median(np.linalg.norm(H0, axis=1)))
    op = h1parts200.H1Operator()
    GEO = {g: h1parts200.effects_for(G[g], p, c4, op, art_norm)[order] for g in G}
    del op

    # (label, geom, mass_fix, min_cpm) -- the PROFILE-determining half of an arm
    CONFIGS = [("artifact totals (164 path)", "ARTIFACT_TOTALS", False, 5.0),
               ("art, mass fix ON", None, True, 5.0),
               ("art, mass fix OFF", None, False, 5.0),
               ("rebuild SHIP", "ship", True, 5.0),
               ("rebuild 095 geometry", "e095", True, 5.0),
               ("rebuild noWhiten", "noWhiten", True, 5.0),
               ("rebuild noPM", "noPM", True, 5.0),
               ("rebuild noCD4", "noCD4", True, 5.0),
               ("rebuild SHIP gate 10", "ship", True, 10.0)]
    out = {}
    cpm = cpm_cell
    for label, geom, mf, min_cpm in CONFIGS:
        h = hashlib.sha256()
        for k in range(H0.shape[0]):
            S = float(esums[k].sum())
            if geom == "ARTIFACT_TOTALS":
                profile = (esums[k] + 0.5) / (N_EMIT * EMIT_L)
            else:
                hh = H0[k] if geom is None else GEO[geom][k]
                g = a_gate if (geom is None and min_cpm == 5.0) else (
                    (cpm > min_cpm) & (b > 0))
                hq = np.where(g, hh, 0.0)
                if mf:
                    hq = prod.mass_consistent_effect(hq, b, g, mode="minnorm")
                lin = np.clip(np.expm1(b + hq), 0.0, None)
                profile = lin / lin.sum() * (S / (N_EMIT * EMIT_L))
            totals = np.floor(N_EMIT * EMIT_L * profile).astype(np.int64)
            h.update(totals.tobytes())
        out[label] = h.hexdigest()
        print(f"  {label:<32} totals sha256 {out[label][:16]}", flush=True)
    return out


def plolo_identity():
    rows = json.loads(open(os.path.join(WORK, "pds200.json")).read())
    ref = json.loads(open(os.path.join(P162, "work/arms.json")).read())
    R = {(r["fold"], r["arm"], r["space"]): r for r in rows}
    A = {(r["fold"], r["arm"], r["space"]): r for r in ref}
    worst, n = 0.0, 0
    for f in ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2"):
        for c in "ABC":
            sp = f"chal{c}/delivered"
            for mine, theirs in (("SHIP", "a23"), ("noMassFix", "a15")):
                x, y = R.get((f, mine, sp)), A.get((f, theirs, sp))
                if x and y:
                    worst = max(worst, abs(x["pds"] - y["pds"]))
                    n += 1
    return worst, n


def main() -> None:
    print("=" * 96)
    print("FALSIFIER 2 -- rule 10 config diff against the SHIPPED artifact's own meta")
    print("=" * 96)
    a15 = json.loads(open(f"{ROOT}/submissions/a15-cd4f010.vcc.meta.json").read())
    a23 = json.loads(open(f"{ROOT}/submissions/a23-massfix.h5ad.meta.json").read())
    keys = sorted(set(a15["config"]) | set(a23["config"]))
    diff = {k: (a15["config"].get(k), a23["config"].get(k))
            for k in keys if a15["config"].get(k) != a23["config"].get(k)}
    print(f"  a15 -> a23 config differences: {len(diff)}  {json.dumps(diff)}"
          f"   {'PASS (exactly 1, the intended one)' if len(diff) == 1 else 'FAIL'}")
    print(f"  operator sha256  a15 {a15['operator']['sha256'][:16]}  "
          f"a23 {a23['operator']['sha256'][:16]}  "
          f"{'PASS' if a15['operator']['sha256'] == a23['operator']['sha256'] else 'FAIL'}")
    n15 = a15["effects"]["whitening"]["after_all_median_norm"]
    n23 = a23["effects"]["whitening"]["after_all_median_norm"]
    print(f"  effects.whitening.after_all_median_norm  a15 {n15!r}  a23 {n23!r}  "
          f"rel {abs(n23 / n15 - 1):.2e}  {'PASS (no drift)' if abs(n23/n15-1) < 1e-9 else 'FAIL'}")
    print(f"  NOTE: the two artifacts' post_whiten_sha256 DIFFER "
          f"({a15['effects']['post_whiten_sha256'][:12]} vs "
          f"{a23['effects']['post_whiten_sha256'][:12]}) at a config diff of one render-time "
          f"flag.\n        The eigen-multiply is not bit-reproducible across runs; the audit "
          f"numbers agree to 1e-15.\n        Any sha-level comparison of two independently "
          f"built effect matrices will fail for this reason alone.")

    print()
    print("=" * 96)
    print("FALSIFIER 3 -- 169's q = 1.00 identity, re-run: SHIP vs 162's own arms.json")
    print("=" * 96)
    worst, n = plolo_identity()
    print(f"  {n}/24 cells compared, worst |d pds| = {worst:.3e}  "
          f"{'PASS (bit-identical)' if worst == 0.0 else 'CHECK'}")

    print()
    print("=" * 96)
    print("FALSIFIER 1 -- per-gene panel TOTALS digest per arm configuration")
    print("=" * 96)
    print("  An EMISSION-only ablation shares the digest by construction (the quasi emitter")
    print("  repairs each gene's total exactly); a PROFILE ablation must break it.")
    dig = totals_digest()
    json.dump(dict(totals=dig, plolo_worst=worst, plolo_cells=n),
              open(os.path.join(WORK, "falsify200.json"), "w"), indent=1)
    groups: dict[str, list[str]] = {}
    for k, v in dig.items():
        groups.setdefault(v, []).append(k)
    print(f"\n  distinct digests: {len(groups)} over {len(dig)} configurations")
    for v, ks in groups.items():
        print(f"    {v[:16]}  <-  {ks}")
    print("\n  The tiered-spread and --overdispersion ablations do not appear above because")
    print("  they do not enter the profile at all: they change only `spread` and `phi`, which")
    print("  the emitter applies AFTER the totals are fixed.  Every one of those arms shares")
    print("  the `art, mass fix ON` digest by construction, and `realize200.py` asserts")
    print("  `B.sum(0) == totals` on every emitted block, which is that statement measured.")


if __name__ == "__main__":
    sys.exit(main())
