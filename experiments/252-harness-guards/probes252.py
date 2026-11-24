"""252 — the guard probes.  ONE module, run twice:

    uv run python probes252.py --mode blind    # against the PRE-252 harness (run before the edit)
    uv run python probes252.py --mode fire     # against the guarded harness (run after the edit)

Each probe constructs the broken input 242 described and records what the harness does with
it: in `blind` mode the silent value, in `fire` mode the raise (verbatim message) -- and every
probe carries a PASS twin, the same call on a correct input, whose numbers must be identical
across the two modes.  The two json outputs (`work/probes-blind.json`, `work/probes-fire.json`)
are the BLIND_pre / FIRE rows AGENTS.md §11.1 rule 10 asks for; `compare252.py` diffs them.

Items:  #5/#4  operator identity (plolo.Operator + E064.load_frozen_operator)
        trap   the Oct-22 axis mismatch at apply time (233's pess operator on the abc axis)
        #13    the on-axis preserve's `or 1.0` at 40% coverage
        #14    endpoints() on an all-zero truth block
        #6     build_fold's cache-validity test (sources / op / PM_BLEND / MIN_CPM)
        zero   the silent zero: telohaec / gse295214_h1 on the k562_gwps fold at f = 0.34

Everything reads its fold caches from THIS experiment's `work/plolo/` (byte copies of 162's)
and writes nothing outside `experiments/252-harness-guards/`.
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import resource
import sys
import time
import traceback
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E162 = ROOT / "experiments/162-retrodicting-validator"
E218 = ROOT / "experiments/218-independent-donors"
_pre = argparse.ArgumentParser(add_help=False)
_pre.add_argument("--harness", default=str(E162),
                  help="directory holding the plolo.py to probe (default: the live 162 dir; "
                       "work/pre252/faux/experiments/162-retrodicting-validator for the byte copy)")
HARNESS = Path(_pre.parse_known_args()[0].harness).resolve()
sys.path.insert(0, str(HARNESS))
sys.path.insert(0, str(E218))

import plolo  # noqa: E402

assert Path(plolo.__file__).resolve().parent == HARNESS, plolo.__file__
import resid218  # noqa: E402  (re-points plolo.WORK at 218's dir on import)
import compose218  # noqa: E402

plolo.WORK = WORK / "plolo"            # 252's byte copies of the fold caches: read-only hits
from plolo import (E062, E064, FOLDS, HELD_OPS, OFFICIAL, PANEL_SET, PROD_OPS,  # noqa: E402
                   Operator, build_fold, cd4_for, endpoints, solve_fraction)

ALPHA, GAMMA, LAM = 0.982, 0.65, 0.0
TAMPERED = WORK / "probe/held_orion_hct116_TAMPERED.npz"
OP233_PESS = ROOT / "experiments/233-october-retention/work/ops/pess__held_orion_hct116.npz"
OP095 = ROOT / "experiments/095-clean-harness/out/whitener_union.npz"


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def digest(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def attempt(fn):
    """Run `fn`; record its value or the exception (class + verbatim message)."""
    t0 = time.time()
    try:
        v = fn()
        out = dict(raised=None, value=v)
    except Exception as e:  # noqa: BLE001 -- the whole point is to record what was raised
        out = dict(raised=type(e).__name__, message=str(e),
                   is_value_error=isinstance(e, ValueError))
        if type(e).__name__ in ("TypeError", "AttributeError"):
            out["trace_tail"] = traceback.format_exc().splitlines()[-1]
    out["seconds"] = round(time.time() - t0, 2)
    return out


def op_summary(op) -> dict:
    return dict(name=getattr(op, "name", "<no `name` attribute in this harness>"),
                sha256=getattr(op, "sha256", None),
                axis_sha256=getattr(op, "axis_sha256", None),
                n_cols=int(op.evals.shape[0]), max_eigenvalue=float(op.evals.max()))


def load_op(path, **kw):
    op = Operator(path, **kw)
    s = op_summary(op)
    del op
    gc.collect()
    return s


# ------------------------------------------------------------------------ probes
def probe_identity(mode: str) -> dict:
    out = {}
    for label, path in (("095_union_maxeig_8.036", OP095),
                        ("donors_union_maxeig_6.257", PROD_OPS / "whitener_union.npz"),
                        ("donors_union_pruned_maxeig_5.895", PROD_OPS / "whitener_union_pruned.npz"),
                        ("TAMPERED_held_orion_hct116", TAMPERED)):
        out[f"Operator({label})"] = attempt(lambda p=path: load_op(p))
        print(f"  Operator({label}): {out[f'Operator({label})']}", flush=True)
    # a swapped file under an explicit expectation, and the correct file under the same one
    out["Operator(union.npz, name='union_pruned')  [SWAP]"] = attempt(
        lambda: load_op(PROD_OPS / "whitener_union.npz", name="union_pruned"))
    out["Operator(union_pruned.npz, name='union_pruned')  [PASS]"] = attempt(
        lambda: load_op(PROD_OPS / "whitener_union_pruned.npz", name="union_pruned"))
    out["Operator(TAMPERED, name='held_orion_hct116')  [SWAP]"] = attempt(
        lambda: load_op(TAMPERED, name="held_orion_hct116"))
    for k in list(out)[-3:]:
        print(f"  {k}: {out[k]}", flush=True)
    # every fold's own operator, called the way the callers call it
    for held, spec in FOLDS.items():
        p = Path(spec["op"]) if spec["op_kind"] != "production" else PROD_OPS / "whitener_union_pruned.npz"
        out[f"Operator(FOLDS[{held}].op)  [PASS]"] = attempt(lambda p=p: load_op(p))
        print(f"  Operator(FOLDS[{held}].op): {out[f'Operator(FOLDS[{held}].op)  [PASS]']}", flush=True)
    return out


def probe_064(mode: str) -> dict:
    """E064.load_frozen_operator called EXACTLY as plolo.build_fold calls it (expected = own hash)."""
    out = {}
    spec = FOLDS["orion_hct116"]
    banks = [E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in spec["sources"]]

    def load(path, expected):
        op, wcols = E064.load_frozen_operator(Path(path), banks, OFFICIAL, PANEL_SET,
                                              expected_sha256=expected)
        r = dict(n_wcols=int(len(wcols)), evecs_shape=list(op.evecs.shape),
                 metadata_donors=op.metadata.get("donors"))
        del op
        gc.collect()
        return r

    genuine = Path(spec["op"])
    out["load_frozen_operator(TAMPERED, expected=sha256_file(TAMPERED))  [TAUTOLOGY]"] = attempt(
        lambda: load(TAMPERED, E062.sha256_file(TAMPERED)))
    out["load_frozen_operator(genuine, expected=sha256_file(genuine))  [PASS]"] = attempt(
        lambda: load(genuine, E062.sha256_file(genuine)))
    out["load_frozen_operator(genuine, expected=<pinned 7cecbf7f...>)  [PASS]"] = attempt(
        lambda: load(genuine, E064.EXPECTED_WHITENER_SHA256["orion_hct116"]))
    out["load_frozen_operator(TAMPERED, expected=<pinned 7cecbf7f...>)  [FIRE, old check]"] = attempt(
        lambda: load(TAMPERED, E064.EXPECTED_WHITENER_SHA256["orion_hct116"]))
    out["load_frozen_operator(genuine, expected='')  [empty expectation]"] = attempt(
        lambda: load(genuine, ""))
    for k, v in out.items():
        print(f"  {k}: {v}", flush=True)
    return out


def probe_axis_trap(mode: str) -> dict:
    """233's pess-panel operator (17,942 cols, self-described) on the abc axis: the Oct-22 trap."""
    out = {}
    fold = build_fold("orion_hct116")
    wcols_abc = np.asarray(fold["wcols"])
    # 233's fold cache stores no `wcols`; rebuild the pess axis the way fold233/064 do and
    # check it against the digest 233 recorded in the operator's own json at fit time.
    pess_panel = tuple(ln.split(",")[0].strip() for ln in
                       (ROOT / "experiments/199-def-runbook/work/def-bundle/pert_counts.csv")
                       .read_text().splitlines()[1:] if ln.strip())
    assert len(pess_panel) == 300
    banks = [E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in FOLDS["orion_hct116"]["sources"]]
    wcols_pess = E064.operator_columns(banks, OFFICIAL, set(pess_panel))
    meta233 = json.loads(OP233_PESS.with_suffix(".json").read_text())
    out["pess axis reconstructed == 233 json axis_sha256 (index digest)"] = bool(
        hashlib.sha256(wcols_pess.tobytes()).hexdigest() == meta233["axis_sha256"])
    assert out["pess axis reconstructed == 233 json axis_sha256 (index digest)"]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    out["n_cols"] = dict(abc=int(len(wcols_abc)), pess=int(len(wcols_pess)),
                         same_length=bool(len(wcols_abc) == len(wcols_pess)),
                         columns_identical=bool(np.array_equal(wcols_abc, wcols_pess)))
    op = Operator(OP233_PESS)
    out["Operator(233 pess op)"] = op_summary(op)
    out["apply(pess op, ABC wcols)  [TRAP]"] = attempt(
        lambda: dict(out_sha256=digest(op.apply(pre * ALPHA, wcols_abc, GAMMA, LAM))))
    out["apply(pess op, PESS wcols)  [PASS]"] = attempt(
        lambda: dict(out_sha256=digest(op.apply(pre * ALPHA, wcols_pess, GAMMA, LAM))))
    if mode == "fire":
        out["apply(pess op, ABC wcols, check_axis=False)  [named bypass, = blind]"] = attempt(
            lambda: dict(out_sha256=digest(op.apply(pre * ALPHA, wcols_abc, GAMMA, LAM,
                                                    check_axis=False))))
    del op
    gc.collect()
    for k, v in out.items():
        print(f"  {k}: {v}", flush=True)
    return out


def probe_preserve_and_zero(mode: str) -> dict:
    """#13 (40% coverage through Operator.apply) and the silent zero, on the k562_gwps fold."""
    out = {}
    held = "k562_gwps"
    fold = build_fold(held)
    wcols = fold["wcols"]
    pre = np.asarray(fold["pre"], dtype=np.float64)
    n = pre.shape[0]
    op = Operator(Path(FOLDS[held]["op"]))
    out["operator"] = op_summary(op)

    rng = np.random.default_rng(252)
    keep40 = np.sort(rng.choice(n, int(round(0.40 * n)), replace=False))
    keep60 = np.sort(rng.choice(n, int(round(0.60 * n)), replace=False))
    for label, keep in (("40pct", keep40), ("60pct", keep60)):
        M = np.zeros_like(pre)
        M[keep] = pre[keep]
        before_raw = float(np.median(np.linalg.norm(M[:, wcols], axis=1)))

        def run(M=M, before_raw=before_raw, keep=keep):
            W = op.apply(M * ALPHA, wcols, GAMMA, LAM)
            return dict(out_sha256=digest(W), carrying=int(len(keep)), n_rows=int(n),
                        median_row_norm_before_raw=before_raw,
                        median_row_norm_out=float(np.median(np.linalg.norm(W[:, wcols], axis=1))),
                        preserve_ratio_implied=float(
                            np.median(np.linalg.norm(W[keep][:, wcols], axis=1))
                            / np.median(np.linalg.norm(M[keep][:, wcols], axis=1))))
        out[f"apply({label} coverage)  [{'FIRE' if label == '40pct' else 'PASS'}]"] = attempt(run)
        print(f"  apply({label}): {out[f'apply({label} coverage)  [{'FIRE' if label == '40pct' else 'PASS'}]']}", flush=True)

    whitened = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    out["apply(full pre)  [PASS]"] = dict(out_sha256=digest(whitened))
    del op
    gc.collect()

    # ---- #14 endpoints on an all-zero truth block (and on the real one)
    sc = fold["scored"]
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
    P = whitened[:, sc]
    e_real = endpoints(P, truth)
    out["endpoints(real truth)  [PASS]"] = {k: float(v) for k, v in e_real.items()}
    e_zero = attempt(lambda: {k: float(v) for k, v in endpoints(P, np.zeros_like(truth)).items()})
    out["endpoints(all-zero truth)  [FIRE]"] = e_zero
    if mode == "fire":
        out["mean_endpoints([real, zero])"] = attempt(
            lambda: plolo.mean_endpoints([e_real, endpoints(P, np.zeros_like(truth))]))
    print(f"  endpoints(real): sgn {e_real['sgn']:.6f} sgn_w {e_real['sgn_w']:.6f}", flush=True)
    print(f"  endpoints(zero): {e_zero}", flush=True)

    # ---- the silent zero
    targets = fold["targets"]
    cd4 = cd4_for(targets)
    resid218.PANEL_SET = set(plolo.PANEL)
    resid218._PREP_CACHE.clear()

    def size(res, **kw):
        frac, per_unit = solve_fraction(whitened, res, 0.34, **kw)
        return dict(fraction=float(frac), per_unit=float(per_unit))

    def coverage(res):
        norms = np.linalg.norm(np.where(np.asarray(res.gene_votes) > 0,
                                        np.asarray(res.mean_effect), 0.0), axis=1)
        return dict(carrying=int((norms > 0).sum()), n_rows=int(norms.size),
                    median_all=float(np.median(norms)),
                    median_carrying=float(np.median(norms[norms > 0])) if (norms > 0).any() else 0.0)

    out["cd4 coverage"] = coverage(cd4)
    out["solve_fraction(cd4, 0.34)  [PASS]"] = attempt(lambda: size(cd4))
    print(f"  cd4: {out['cd4 coverage']} -> {out['solve_fraction(cd4, 0.34)  [PASS]']}", flush=True)
    for src in ("telohaec", "gse295214_h1"):
        res = resid218.source_residual(src, targets, "gene_norm")
        out[f"{src} coverage"] = coverage(res)
        out[f"solve_fraction({src}, 0.34)  [FIRE]"] = attempt(lambda r=res: size(r))
        ref = compose218.solve_fraction_carrying(whitened, res, 0.34)
        out[f"compose218.solve_fraction_carrying({src}, 0.34)"] = dict(fraction=float(ref[0]),
                                                                       per_unit=float(ref[1]))
        if mode == "fire":
            out[f"solve_fraction({src}, 0.34, sizing='carrying')"] = attempt(
                lambda r=res: size(r, sizing="carrying"))
            out[f"solve_fraction({src}, 0.34, sizing='legacy_all_rows')"] = attempt(
                lambda r=res: size(r, sizing="legacy_all_rows"))
        for k in [k for k in out if src in k]:
            print(f"  {k}: {out[k]}", flush=True)
    if mode == "fire":
        out["solve_fraction(cd4, 0.34, sizing='carrying')"] = attempt(lambda: size(cd4, sizing="carrying"))
        out["solve_fraction(cd4, 0.34, sizing='legacy_all_rows')  [= default where per_unit > 0]"] = attempt(
            lambda: size(cd4, sizing="legacy_all_rows"))
    return out


def probe_cache(mode: str) -> dict:
    """#6: the same cache file under changed inputs.  A rebuild is short-circuited by a sentinel
    so the probe costs nothing: reaching `prepare_main_source` == the cache was (rightly) refused."""
    out = {}
    held = "k562_gwps"
    spec = FOLDS[held]
    saved = dict(spec)
    saved_pm, saved_cpm = plolo.PM_BLEND, plolo.MIN_CPM

    class RebuildAttempted(RuntimeError):
        pass

    def sentinel(*a, **k):
        raise RebuildAttempted("REBUILD ATTEMPTED: the cache was not used")

    original = E062.prepare_main_source
    E062.prepare_main_source = sentinel

    def call(label):
        r = attempt(lambda: (lambda f: dict(audit_sources=f["audit"]["sources"],
                                            audit_op=f["audit"]["op"],
                                            returned_sources=list(f["sources"]),
                                            returned_op=str(Path(f["op_path"]).relative_to(ROOT)),
                                            audit_cache=f["audit"].get("cache"),
                                            n_targets=len(f["targets"])))(build_fold(held)))
        out[label] = r
        print(f"  {label}: {r}", flush=True)

    try:
        call("build_fold(unchanged)  [PASS: hit]")
        FOLDS[held]["sources"] = ("orion_hek293t",)
        call("build_fold(sources changed to ('orion_hek293t',))  [FIRE: must miss]")
        FOLDS[held].update(saved)
        FOLDS[held]["op"] = HELD_OPS / "held_orion_hek293t.npz"
        call("build_fold(op changed to held_orion_hek293t)  [FIRE: must miss]")
        FOLDS[held].update(saved)
        plolo.PM_BLEND = 0.25
        call("build_fold(PM_BLEND 0.5 -> 0.25)  [FIRE: must miss]")
        plolo.PM_BLEND = saved_pm
        # a MIN_CPM that CHANGES this fold's gate (every gene the k562 bank carries is above
        # 5 CPM, so 10 -> 5 leaves the gate identical and a legacy cache is rightly still valid)
        zc = np.load(plolo.WORK / f"fold-{held}.npz", allow_pickle=True)
        b_ctrl = np.asarray(zc["b_ctrl"], dtype=np.float64)
        expr = np.expm1(b_ctrl)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate10 = np.asarray(zc["gate"], dtype=bool)
        alt = next(v for v in (5.0, 20.0, 50.0, 100.0, 200.0, 500.0)
                   if not np.array_equal((cpm > v) & (b_ctrl > 0), gate10))
        out["MIN_CPM alternative that changes the gate"] = dict(
            min_cpm=alt, gated_at_10=int(gate10.sum()), gated_at_alt=int(((cpm > alt) & (b_ctrl > 0)).sum()),
            gate_identical_at_5=bool(np.array_equal((cpm > 5.0) & (b_ctrl > 0), gate10)))
        plolo.MIN_CPM = alt
        call(f"build_fold(MIN_CPM 10 -> {alt:g})  [FIRE: must miss]")
        plolo.MIN_CPM = saved_cpm
        call("build_fold(restored)  [PASS: hit]")
    finally:
        E062.prepare_main_source = original
        FOLDS[held].update(saved)
        plolo.PM_BLEND, plolo.MIN_CPM = saved_pm, saved_cpm
    if mode == "fire":
        out["fold_cache_status(unchanged)"] = attempt(lambda: plolo.fold_cache_status(held))
        FOLDS[held]["sources"] = ("orion_hek293t",)
        out["fold_cache_status(sources changed)"] = attempt(lambda: plolo.fold_cache_status(held))
        FOLDS[held].update(saved)
        for k in list(out)[-2:]:
            print(f"  {k}: {out[k]}", flush=True)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(parents=[_pre])
    ap.add_argument("--mode", choices=("blind", "fire"), required=True)
    args = ap.parse_args()
    t0 = time.time()
    out = dict(mode=args.mode, harness=str(Path(plolo.__file__).resolve()),
               plolo_sha256=E062.sha256_file(Path(plolo.__file__).resolve()),
               e064_file=str(Path(E064.__file__).resolve()),
               e064_sha256=E062.sha256_file(Path(E064.__file__).resolve()),
               tampered=json.loads((WORK / "probe/tampered.json").read_text()))
    print(f"MODE {args.mode}: plolo {out['plolo_sha256'][:12]} e064 {out['e064_sha256'][:12]}", flush=True)
    for name, fn in (("identity", probe_identity), ("e064", probe_064),
                     ("axis_trap", probe_axis_trap), ("preserve_endpoints_zero", probe_preserve_and_zero),
                     ("cache", probe_cache)):
        print(f"\n=== {name} ===", flush=True)
        out[name] = fn(args.mode)
        out[name]["_peak_gib"] = round(peak_gib(), 2)
        gc.collect()
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gib"] = round(peak_gib(), 2)
    path = WORK / f"probes-{args.mode}.json"
    path.write_text(json.dumps(out, indent=1, default=str))
    print(f"\nwrote {path} ({out['seconds']}s, peak {out['peak_rss_gib']} GiB)")


if __name__ == "__main__":
    main()
