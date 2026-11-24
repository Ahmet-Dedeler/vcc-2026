"""193 — screen alternative eigenvalue maps on PLOLO (162), a15 config, challenge frames.

The arm is a *pure operator swap*: everything else in the a15 pipeline is held fixed
(alpha 0.982, CD4 target f 0.10, no mass fix, the same frozen operator EIGENBASIS), and only
the eigenvalue map `f(ev)` changes.  `plolo.Operator.apply` renormalises the median row norm
afterwards, so this is a shape-only arm by construction.

KNOWN-ANSWER GATE (AGENTS.md 11.1 rule 10)
------------------------------------------
Map `prod(lam=0)` IS a15's operator.  Its four folds x three challenge frames must reproduce
162's stored `work/arms2.json` rows for arm `a15`.  `pds` is a rank statistic and 162 measured
it bit-identical across processes on all 165 of its rows, so the gate is EQUALITY at 1e-9
absolute; anything larger is a wiring fault, not float noise.  The run aborts if it fails.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import importlib.util
import json
import os
import resource
import sys
import time
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")
os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "8")

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E162 = ROOT / "experiments/162-retrodicting-validator"
sys.path.insert(0, str(E162))
sys.path.insert(0, str(ROOT / "src"))

import plolo                                        # noqa: E402  (162's harness, verbatim)
from plolo import FOLDS, build_fold, cd4_for, deliver, endpoints, solve_fraction  # noqa: E402
from run_arms import CONFIGS                        # noqa: E402

sys.path.insert(0, str(HERE))
import maps193                                      # noqa: E402

A15 = CONFIGS["a15"]
FOLD_ORDER = ("orion_hct116", "orion_hek293t", "k562_gwps", "kolf2")
SPAN_PDS = 0.452           # 162's board span for pds


def rss_gib() -> float:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 2**30 if sys.platform == "darwin" else peak / 2**20


def sha256_file(path: Path, buf: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while chunk := fh.read(buf):
            h.update(chunk)
    return h.hexdigest()


def cd4_cached(targets):
    key = hashlib.sha256("\n".join(targets).encode()).hexdigest()[:16]
    path = WORK / f"cd4-{key}.npz"
    if path.exists():
        z = np.load(path)
        return plolo.E062.CD4Prepared(
            mode=str(z["mode"]), mean_effect=z["mean_effect"], gene_votes=z["gene_votes"],
            target_condition_count=z["tcc"], condition_details=())
    c = cd4_for(targets)
    np.savez(path, mode=c.mode, mean_effect=np.asarray(c.mean_effect),
             gene_votes=np.asarray(c.gene_votes),
             tcc=np.asarray(c.target_condition_count))
    return c


def frames_for(fold):
    out = []
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        out.append((f"chal{ctx}", b, (cpm > plolo.MIN_CPM) & (b > 0)))
    return out


def run_fold(held: str, op_override: Path | None, maps_filter, n_basis_override=None,
             registry: str = "a"):
    """Score every map on one PLOLO fold.  Returns list of row dicts."""
    fold = build_fold(held)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    wcols = np.asarray(fold["wcols"])
    cd4 = cd4_cached(fold["targets"])
    frames = frames_for(fold)

    op_path = Path(op_override) if op_override else Path(fold["op_path"])
    op_sha = sha256_file(op_path)
    with np.load(op_path) as z:
        evals = np.asarray(z["evals"], dtype=np.float64)
        evecs = np.asarray(z["evecs"], dtype=np.float64)
    assert evecs.shape[0] == len(wcols) == len(evals), (
        f"operator axis {evecs.shape} vs fold wcols {len(wcols)}")
    meta_path = op_path.with_suffix(".json")
    n_basis = n_basis_override
    if n_basis is None and meta_path.exists():
        m = json.loads(meta_path.read_text())
        n_basis = int(m.get("basis", [0, 0])[0]) or int(sum(
            m.get("basis_rows_by_source", {}).values()))
    assert n_basis, f"basis row count unknown for {op_path}"
    print(f"[{held}] operator {op_path.relative_to(ROOT)} sha {op_sha[:8]} "
          f"p={len(evals)} n={n_basis} maxev={evals.max():.4f}", flush=True)

    # --- the part of the pipeline that does not depend on f(ev) -------------------------
    pre = np.asarray(fold["pre"], dtype=np.float64) * A15["alpha"]
    sub = pre[:, wcols]
    before = float(np.median(np.linalg.norm(sub, axis=1))) or 1.0
    Y = sub @ evecs                     # (n_targets, p), reused by every map
    del sub
    gc.collect()

    builder = maps193.build_maps if registry == "a" else maps193.build_maps_b
    fmaps = builder(evals, n_basis, gamma=A15["gamma"])
    rows = []
    for name, w in fmaps.items():
        if maps_filter and name not in maps_filter:
            continue
        t0 = time.time()
        t = (Y * w) @ evecs.T
        after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
        whitened = pre.copy()
        whitened[:, wcols] = t * (before / after)
        del t
        frac, per_unit = solve_fraction(whitened, cd4, A15["cd4_f"])
        final, _ = plolo.E062.add_post_whiten_residual(whitened, cd4, frac, wcols)
        final[:, fold["panel_positions"]] = 0.0
        del whitened
        rec = dict(fold=held, map=name, op=str(op_path.relative_to(ROOT)), op_sha256=op_sha,
                   n_basis=int(n_basis), p=int(len(evals)), cd4_fraction=float(frac))
        e = endpoints(final[:, fold["scored"]], truth)
        rec["intended_pds"] = e["pds"]
        rec["intended_cos"] = e["cos"]
        pds = []
        for fname, b, gate in frames:
            mat = deliver(final, b, gate, A15["mass"])
            e = endpoints(mat[:, fold["scored"]], truth)
            rec[f"{fname}_pds"] = e["pds"]
            rec[f"{fname}_cos"] = e["cos"]
            rec[f"{fname}_sgn_w"] = e["sgn_w"]
            rec[f"{fname}_sse"] = e["sse"]
            pds.append(e["pds"])
            del mat
        rec["chal_pds"] = float(np.mean(pds))
        rows.append(rec)
        del final
        gc.collect()
        print(f"  {held:<14} {name:<18} chal_pds {rec['chal_pds']:.6f} "
              f"intended {rec['intended_pds']:.6f}  ({time.time()-t0:.1f}s, "
              f"peakRSS {rss_gib():.1f}G)", flush=True)
    del evecs, Y
    gc.collect()
    return rows


def gate(rows):
    """Known-answer gate: `prod(lam=0)` must reproduce 162's stored a15 rows."""
    ref = {(r["fold"], r["arm"], r["space"]): r for r in
           json.loads((E162 / "work/arms2.json").read_text())}
    worst = 0.0
    checked = 0
    for r in rows:
        if r["map"] != "prod(lam=0)" or r.get("op_sha256_is_default") is False:
            pass
        if r["map"] != "prod(lam=0)":
            continue
        for c in "ABC":
            k = (r["fold"], "a15", f"chal{c}/delivered")
            if k not in ref:
                continue
            d = abs(r[f"chal{c}_pds"] - ref[k]["pds"])
            worst = max(worst, d)
            checked += 1
    return worst, checked


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", nargs="*", default=list(FOLD_ORDER))
    ap.add_argument("--maps", nargs="*", default=None)
    ap.add_argument("--op-override", default=None,
                    help="path to an alternative operator (same eigen-axis) for ALL folds")
    ap.add_argument("--n-basis", type=int, default=None)
    ap.add_argument("--tag", default="phase1")
    ap.add_argument("--registry", default="a", choices=["a", "b"])
    ap.add_argument("--no-gate", action="store_true")
    args = ap.parse_args()

    out_path = WORK / f"{args.tag}.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["map"], r["op_sha256"]) for r in rows}
    for held in args.folds:
        new = run_fold(held, args.op_override, set(args.maps) if args.maps else None,
                       args.n_basis, args.registry)
        new = [r for r in new if (r["fold"], r["map"], r["op_sha256"]) not in done]
        rows.extend(new)
        out_path.write_text(json.dumps(rows, indent=1))

    if not args.no_gate and args.op_override is None:
        worst, checked = gate(rows)
        print(f"\nKNOWN-ANSWER GATE: prod(lam=0) vs 162 a15, {checked} cells, "
              f"worst |dpds| = {worst:.3e}")
        if checked and worst > 1e-9:
            raise SystemExit("GATE FAILED — the harness does not reproduce a15")
        print("GATE PASSED" if checked else "GATE NOT EVALUABLE (no folds compared)")
    print(f"wrote {out_path}")


if __name__ == "__main__":
    main()
