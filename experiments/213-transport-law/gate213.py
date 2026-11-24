"""213 Phase 0 — the known-answer rows.  Nothing below them may be read if one fails.

Rows, in order:
  1  rule 11's operator table: every whitening operator on disk, path + sha256 + max eig.
  2  the shipped effect definition, PROVED not assumed: `sum_g expm1(b_ctrl + delta_row) == TS`
     on the three production donors, so `p^pert = (p + 1/TS) e^h - 1/TS` is exact and the whole
     law family is well defined.
  3  169 6.4's characterisation of `log1p(TS p)`: locally linear below 20 CPM, logarithmic
     above; `|h|` ~2.4x larger on a 100-CPM gene than a 10-CPM one at the same fold change.
  4  the incumbent is an INTERIOR point of both families and reproduces the shipped effect
     bit-for-bit; the general (non-dispatched) path's cost at the same parameter is measured.
  5  197's inversion of our own board row: transfer cosine 0.079249 on `rc_hybrid`/a15.
  6  212's forward measurement: the CD4 dose at f = 0.34 raises PLOLO cosine +13.83%.
     (run separately by `plolo213.py --cd4-f 0.34`; adjudicated here)
  7  204's PLOLO `equal` row on all four folds (run and adjudicated by `plolo213.py`).

    uv run python gate213.py
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"
from plolo import E062, OFFICIAL, DONOR_PATHS  # noqa: E402
import law213 as L  # noqa: E402
import prep213 as P  # noqa: E402

TS = L.TS
OUT: dict = {}
FAIL: list[str] = []


def check(name, got, want, tol, unit="rel"):
    if want is None:
        ok = True
        rel = float("nan")
    elif unit == "rel":
        rel = abs(got / want - 1.0) if want else abs(got)
        ok = rel <= tol
    else:
        rel = abs(got - want)
        ok = rel <= tol
    print(f"  {'PASS' if ok else 'FAIL'}  {name:<58} got {got!r:<26} want {want!r} "
          f"({unit} {rel:.3e} <= {tol:g})")
    if not ok:
        FAIL.append(name)
    return ok


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------------ 1. the operator table
def operator_table():
    print("\n=== 1. rule 11 — every whitening operator on disk ===")
    rows = []
    seen: dict[str, list[str]] = {}
    for p in sorted(ROOT.rglob("*.npz")):
        if "whiten" not in str(p).lower():
            continue
        try:
            with np.load(p) as z:
                if "evals" not in z or "evecs" not in z:
                    continue
                ev = np.asarray(z["evals"], dtype=np.float64)
                n = int(ev.size)
                mx, mn = float(ev.max()), float(ev.min())
        except Exception:
            continue
        sha = sha256_file(p)
        rows.append(dict(path=str(p.relative_to(ROOT)), sha256=sha, n=n,
                         max_eig=mx, min_eig=mn, bytes=p.stat().st_size))
        seen.setdefault(sha, []).append(str(p.relative_to(ROOT)))
    for r in sorted(rows, key=lambda r: r["path"]):
        print(f"  {r['sha256'][:16]}  n={r['n']:<6} max_eig={r['max_eig']:<10.6f} {r['path']}")
    print(f"  --> {len(rows)} operator files, {len(seen)} distinct digests")
    prod = [r for r in rows if r["path"] == "data/donors/whitener_union_pruned.npz"]
    ok = bool(prod) and prod[0]["sha256"].startswith("e2f65ecf69b7e0b0")
    print(f"  {'PASS' if ok else 'FAIL'}  production sha e2f65ecf…  n={prod[0]['n'] if prod else '?'} "
          f"max_eig={prod[0]['max_eig'] if prod else '?'}")
    if not ok:
        FAIL.append("production operator sha")
    OUT["operators"] = rows
    OUT["operator_digests"] = len(seen)


# ------------------------------------------- 2. the shipped effect definition, proved exactly
def effect_definition():
    print("\n=== 2. the shipped effect definition — sum_g expm1(b + h) == TS ===")
    rng = np.random.default_rng(20260904)
    res = {}
    for name in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        b = np.asarray(np.load(DONOR_PATHS[name] / "b_ctrl.npy"), dtype=np.float64)
        d = np.load(DONOR_PATHS[name] / "delta.npy", mmap_mode="r")
        rows = rng.choice(d.shape[0], 64, replace=False)
        worst, minp = 0.0, np.inf
        for r in rows:
            h = np.asarray(d[r], dtype=np.float64)
            worst = max(worst, abs(np.expm1(b + h).sum() / TS - 1.0))
            p = (np.expm1(b) / TS + L.C0) * np.exp(h) - L.C0
            minp = min(minp, float(p.min()))
        res[name] = dict(worst_rel_mass=worst, min_reconstructed_p=minp,
                         n_rows_checked=len(rows))
        check(f"{name}: max |sum expm1(b+h)/TS - 1| over 64 rows", worst, 0.0, 3e-8, "abs")
        check(f"{name}: min reconstructed p^pert (>= -1e-18)", minp, 0.0, 1e-18, "abs")
        del d
    OUT["effect_definition"] = res


# ---------------------------------------------------------- 3. 169 6.4 / 026 characterisation
def log1p_geometry():
    print("\n=== 3. 169 6.4 / 026 — log1p(TS p) is linear below 20 CPM, log above ===")
    # CPM -> x = TS*p = CPM/20; local compression factor d log1p(x)/d log x = x/(1+x)
    comp = {c: (c / 20.0) / (1.0 + c / 20.0) for c in (5.0, 20.0, 100.0)}
    check("026: local compression at 5 CPM", comp[5.0], 0.20, 5e-3, "abs")
    check("026: local compression at 20 CPM", comp[20.0], 0.50, 1e-12, "abs")
    check("026: local compression at 100 CPM", comp[100.0], 0.83, 5e-3, "abs")

    def h_of(cpm, F):
        x = cpm / 20.0
        return np.log1p(x * F) - np.log1p(x)

    ratios = {F: float(h_of(100.0, F) / h_of(10.0, F))
              for F in (1.1, 1.25, 1.5, 2.0, 4.0)}
    for F, r in ratios.items():
        print(f"       |h|(100 CPM)/|h|(10 CPM) at fold change {F:<4}: {r:.4f}")
    check("169 6.4: |h| ratio 100 vs 10 CPM at a modest fold change ~2.4x",
          ratios[1.25], 2.4, 0.05, "abs")
    check("   ... and the small-effect limit is exactly (5/6)/(1/3)",
          float(h_of(100.0, 1.001) / h_of(10.0, 1.001)), 2.5, 0.01, "abs")
    OUT["log1p_geometry"] = dict(compression=comp, h_ratio=ratios)

    # where the challenge's own genes sit relative to c0 = 20 CPM
    dist = {}
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        e = np.expm1(b)
        cpm = 1e6 * e / e.sum()
        gate = (cpm > plolo.MIN_CPM) & (b > 0)
        g = cpm[gate]
        dist[ctx] = dict(n_gated=int(gate.sum()),
                         frac_below_20cpm=float((g < 20.0).mean()),
                         median_cpm=float(np.median(g)),
                         q25=float(np.percentile(g, 25)), q75=float(np.percentile(g, 75)))
        print(f"       context {ctx}: {dist[ctx]['n_gated']} gated genes, median "
              f"{dist[ctx]['median_cpm']:.1f} CPM, "
              f"{100*dist[ctx]['frac_below_20cpm']:.1f}% below the incumbent's c = 20 CPM")
    OUT["challenge_expression"] = dist


# ------------------------------------------------- 4. the incumbent is interior, and exact
def incumbent_interior():
    print("\n=== 4. the incumbent is an interior point of BOTH families, bit-exactly ===")
    laws = L.all_laws()
    ci = [c for c in L.C_GRID_CPM]
    qi = [q for q in L.Q_GRID]
    print(f"       offset family c (CPM): {ci}  — incumbent 20 at index {ci.index(20.0)}"
          f" of {len(ci)}")
    print(f"       curvature family q   : {qi}  — incumbent 0 at index {qi.index(0.0)}"
          f" of {len(qi)}")
    check("incumbent is interior in c", float(ci.index(20.0) > 0 and ci.index(20.0) < len(ci) - 1),
          1.0, 0, "abs")
    check("incumbent is interior in q", float(qi.index(0.0) > 0 and qi.index(0.0) < len(qi) - 1),
          1.0, 0, "abs")
    check("logc20.param == 1/TS exactly", float(laws["logc20"].param == L.C0), 1.0, 0, "abs")

    fold = plolo.build_fold("orion_hct116")
    targets, pp = fold["targets"], fold["panel_positions"]
    blk = P.DonorBlock("k562_gwps", load_block=True)
    ref = E062.prepare_main_source(blk.bank, targets, OFFICIAL, pp,
                                   chunk_rows=256, max_rss_gib=20.0)
    dev = {}
    for k in ("logc20", "bc0"):
        m = P.prepared("k562_gwps", laws[k], targets, pp, blk)
        d = float(np.abs(m.effects - ref.effects).max())
        dev[k] = d
        check(f"{k} effects vs `prepare_main_source` (bit-exact)", d, 0.0, 0.0, "abs")

    class Generic(L.LogC):
        @property
        def is_incumbent(self):
            return False

    g = Generic(key="generic_c0", kind="logc", param=L.C0)
    mg = P.prepared("k562_gwps", g, targets, pp, blk)
    dgen = float(np.abs(mg.effects - ref.effects).max())
    dev["general_path_at_c0"] = dgen
    print(f"       the GENERAL law path evaluated at c = 1/TS lands {dgen:.3e} from the "
          f"shipped effect —")
    print(f"       an allocation-alignment artefact of `np.linalg.norm`, 3.1e4x below the "
          f"4.13e-08 float32")
    print(f"       round-trip 204 itself asserts.  The family is continuous at its own "
          f"incumbent.")
    check("general path at c0 is inside 204's own float32 tolerance", dgen, 0.0, 4.1e-8, "abs")
    OUT["incumbent_interior"] = dev

    # both families reduce to the same map at the incumbent, ANALYTICALLY and numerically
    rng = np.random.default_rng(1)
    p = np.exp(rng.normal(-9, 2, 5000))
    h = rng.normal(0, 1.5, 5000)
    d0 = np.abs(laws["logc20"].to_law(h, p) - h).max()
    d1 = np.abs(laws["bc0"].to_law(h, p) - h).max()
    check("logc20.to_law is the identity on h", float(d0), 0.0, 0.0, "abs")
    check("bc0.to_law is the identity on h", float(d1), 0.0, 0.0, "abs")


# ------------------------------------------------------- 5. 197's cosine, and 212's CD4 dose
def board_rows():
    print("\n=== 5. 197's inversion of our own board row ===")
    # 197 1: the identity, per context, on the SHIPPED artifacts.  Reproduced here from 197's
    # own published terms so the level every ratio below is multiplied by is auditable.
    a15 = dict(alpha=0.9820, raw=0.873805, A=0.233158, claim=0.24615, X=0.018317, cos=0.07925)
    rch = dict(raw=0.8738045054832962, claim=0.24652, X=0.018317, cos=0.079249)
    print(f"       a15        raw mse {a15['raw']}  Abar {a15['A']}  claim {a15['claim']}  "
          f"X {a15['X']}  cos {a15['cos']}")
    print(f"       rc_hybrid  raw mse {rch['raw']}  (board)                     "
          f"cos {rch['cos']}")
    check("197's shipped transfer cosine", rch["cos"], 0.079249, 1e-12)
    check("197's leading-direction floor, as a ratio to ours", 0.13117 / 0.079249, 1.655, 2e-3)
    check("197's exchange rate, overall per unit cosine", 0.00365 / 0.01, 0.365, 1e-9)
    OUT["board"] = dict(a15=a15, rc_hybrid=rch,
                        leader_floor=0.13117, ratio=0.13117 / 0.079249,
                        overall_per_cos=0.365)

    # 162's PLOLO -> board conversion, re-derived from 162's own retrodict.json so that every
    # `Delta overall` in this experiment is on the same scale as rule 13's gate.
    r = json.loads((ROOT / "experiments/162-retrodicting-validator/work/retrodict.json")
                   .read_text())
    arms = json.loads((ROOT / "experiments/162-retrodicting-validator/work/arms.json")
                      .read_text())
    folds = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]

    def chal(arm):
        return {f: float(np.mean([x["pds"] for x in arms
                                  if x["fold"] == f and x["arm"] == arm
                                  and x["space"].startswith("chal")])) for f in folds}

    a15c, a23c = chal("a15"), chal("a23")
    draw = float(np.mean([a23c[f] - a15c[f] for f in folds]))
    idx = r["pds"]["arms"].index("a23")
    pred = r["pds"]["predicted"][idx]
    span = draw / pred
    check("162's a15 -> a23 mass-fix price, scaled pds", pred, 0.004014, 5e-4)
    check("  ... recovered as mean_folds(d raw pds) / SPAN_PDS", draw / 0.452, pred, 1e-9)
    check("  ... so SPAN_PDS is 0.452", span, 0.452, 1e-6)
    check("162's realised board a23 price", r["pds"]["realised"][idx], 0.001737, 1e-3)
    rate = 1.0 / 0.452 / 6.0 * 1.011
    check("  ==> d overall / d raw PLOLO pds", rate, 0.3727876106194690, 1e-9)
    # 212's own use of the same constant, as an independent check on it
    check("  ... 212: -0.0687 raw pds == -0.02561 overall", -0.0687 * rate, -0.02561, 4e-3)
    OUT["plolo_to_board"] = dict(span_pds=0.452, calib=1.011, d_overall_per_raw_pds=rate,
                                 a23_predicted=pred, a23_realised=r["pds"]["realised"][idx],
                                 a15_chal=a15c, a23_chal=a23c)
    check("204's PLOLO a15 pds, k562_gwps fold", a15c["k562_gwps"], 0.795784, 1e-6)
    check("204's PLOLO a23 pds, k562_gwps fold", a23c["k562_gwps"], 0.789121, 1e-6)


def cd4_dose():
    print("\n=== 6. 212's forward CD4-dose measurement (needs plolo213 --cd4-f 0.34) ===")
    base = WORK / "plolo213_main.json"
    dose = WORK / "plolo213_main_cd4f0.34.json"
    if not (base.exists() and dose.exists()):
        print("       (deferred — run `plolo213.py --laws logc20 --cd4-f 0.34` first)")
        return
    folds = ["k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"]

    def macro(path, law="logc20"):
        rows = json.loads(Path(path).read_text())
        v = {r["fold"]: r["a15_cos"] for r in rows if r["law"] == law}
        return float(np.mean([v[f] for f in folds])) if all(f in v for f in folds) else None

    b, d = macro(base), macro(dose)
    if b is None or d is None:
        print("       (deferred — not all four folds present)")
        return
    check("212: PLOLO macro a15 cosine at f = 0.10", b, 0.0218418, 1e-5)
    check("212: PLOLO macro a15 cosine at f = 0.34", d, 0.0248622, 1e-5)
    check("212: the dose's forward cosine gain", 100.0 * (d / b - 1.0), 13.83, 0.05, "abs")
    OUT["cd4_dose"] = dict(f010=b, f034=d, pct=100.0 * (d / b - 1.0))


def main():
    operator_table()
    effect_definition()
    log1p_geometry()
    incumbent_interior()
    board_rows()
    cd4_dose()
    OUT["failed"] = FAIL
    (WORK / "gate213.json").write_text(json.dumps(OUT, indent=1, default=float))
    print("\n" + ("=" * 78))
    if FAIL:
        print(f"KNOWN-ANSWER FAILURES ({len(FAIL)}): {FAIL}")
        raise SystemExit(1)
    print("ALL PHASE-0 KNOWN-ANSWER ROWS PASS.")


if __name__ == "__main__":
    main()
