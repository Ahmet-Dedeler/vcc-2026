"""213 Phase 2 — sweep the TRANSPORT LAW on PLOLO, with the COSINE as the primary endpoint.

WHAT VARIES, AND WHAT DOES NOT
------------------------------
The law varies in BOTH of its halves at once, because that is what a law is:

  * donor side   -- the effect is re-expressed as `g = f_c(p^pert_don) - f_c(p^ctrl_don)`
                    BEFORE 062's centring and `EFFECT_SCALE/raw_median` normalisation
                    (`prep213.prepared`), so the consensus, the operator and the CD4 solve all
                    act in the law's own space;
  * render side  -- the delivered abundance is `f_c^-1(f_c(p^ctrl_rec) + g)` instead of
                    `expm1(b_rec + g)/TS` (`law213.deliver_law`).

Everything else is 162's/204's byte-unchanged code, imported not copied: the four folds, the
frozen per-fold operators, `alpha` 0.982, `gamma` 0.65, `lam` 0.0, the PerturbMatch 50:50
blend, the CD4 residual re-solved at every arm exactly as production does, the three challenge
frames A/B/C, the read-back `log1p(TS p) - b_ctrl`, and the truth.

THE READ-BACK IS NOT A FREE PARAMETER.  `cell-eval2` scores log-pseudobulk, so every arm is
scored in `log1p(TS p)` space whatever space it transported in.  Only the map from (donor
effect, recipient basal) to `p` is ours.

KNOWN-ANSWER ROW (rule 10's corollary).  `logc20` and `bc0` are the shipped law reached from
the two different families, and both must reproduce 204's published `equal` row EXACTLY:

    k562_gwps      a23_pds 0.7891213079583496   a15_pds 0.7957840906989219
                   a23_cos 0.03146796694667673  q 0.5016431146906033  q_w 0.5048799395717629
    orion_hct116   a23_pds 0.7918227737342401   a23_cos 0.026040412606551395
    orion_hek293t  a23_pds 0.7773753414904073   a23_cos 0.02247089282968584
    kolf2          a23_pds 0.6101280255692021   a23_cos 0.007559021017607683

    uv run python plolo213.py --laws all
"""
from __future__ import annotations

import argparse
import gc
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(ROOT / "experiments/204-repression-weighted-mixing"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
assert (plolo.WORK / "fold-kolf2.npz").exists(), plolo.WORK

import score204  # noqa: E402
score204.WORK = plolo.WORK.parent

from plolo import (E062, E064, FOLDS, OFFICIAL, PANEL_SET, PM_PATH, PROD_OPS,  # noqa: E402
                   Operator, build_fold, cd4_for, deliver, endpoints, solve_fraction)
from score204 import PROD, _pds_per_target  # noqa: E402

import law213 as L  # noqa: E402
import prep213 as P  # noqa: E402

# rule 11: a path and a sha, never a description.  The production operator is pinned in full
# (212 1); the three held operators are pinned on 162's/212's recorded 16-hex prefix, because
# no full digest for them is published anywhere in the repo.  Nine near-identically-named
# operators exist on disk, so the assert runs BEFORE the file is loaded.
OP_SHA_FULL = {
    "whitener_union_pruned.npz":
        "e2f65ecf69b7e0b0b817d02d850c7be07a78b442644a239c4b0e520cda1db123",
}
OP_SHA_PREFIX = {
    "held_k562_gwps.npz": "13d988fcc3b129d2",
    "held_orion_hct116.npz": "7cecbf7f36e40d22",
    "held_orion_hek293t.npz": "97e486e0d1a805e1",
    "whitener_union_pruned.npz": "e2f65ecf69b7e0b0",
}

KA = {
    "k562_gwps": dict(a23_pds=0.7891213079583496, a15_pds=0.7957840906989219,
                      a23_cos=0.03146796694667673, a15_sgn=0.5016431146906033,
                      a15_sgn_w=0.5048799395717629),
    "orion_hct116": dict(a23_pds=0.7918227737342401, a23_cos=0.026040412606551395),
    "orion_hek293t": dict(a23_pds=0.7773753414904073, a23_cos=0.02247089282968584),
    "kolf2": dict(a23_pds=0.6101280255692021, a23_cos=0.007559021017607683),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for blk in iter(lambda: fh.read(1 << 22), b""):
            h.update(blk)
    return h.hexdigest()


_META: dict[str, P.DonorBlock] = {}


def meta(name: str) -> P.DonorBlock:
    if name not in _META:
        _META[name] = P.DonorBlock(name, load_block=False)
    return _META[name]


def law_sources(held: str, fold: dict, law, use_pm: bool = True):
    """The fold's sources in LAW space, with the PerturbMatch blend exactly as 162.

    PerturbMatch is a limma-voom log-fold-change, not a `log1p(TS p)` difference, so it has no
    law-space form and is left alone.  `--no-pm` prices that impurity.
    """
    spec = FOLDS[held]
    out = []
    for n in spec["sources"]:
        sp = P.prepared(n, law, fold["targets"], fold["panel_positions"], meta(n))
        if n == "k562_gwps" and use_pm:
            pm_bank = E064.load_perturbmatch(PM_PATH)
            pm = E062.prepare_main_source(pm_bank, fold["targets"], OFFICIAL,
                                          fold["panel_positions"], chunk_rows=256,
                                          max_rss_gib=20.0)
            shared = E064.shared_k562_pm_columns(meta(n).bank, pm_bank, OFFICIAL, PANEL_SET)
            eff = np.asarray(sp.effects, dtype=np.float64).copy()
            both = sp.target_available & pm.target_available
            rows = np.flatnonzero(both)
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    (1 - plolo.PM_BLEND) * sp.effects[np.ix_(rows, shared)]
                    + plolo.PM_BLEND * pm.effects[np.ix_(rows, shared)])
            sp = E062.PreparedSource(
                name="k562_gwps_pm50", effects=eff,
                target_available=np.asarray(sp.target_available, bool).copy(),
                mean_on_axis=np.asarray(sp.mean_on_axis, float).copy(),
                raw_median_norm=float(sp.raw_median_norm),
                mapped_gene_count=int(sp.mapped_gene_count),
                duplicate_symbols=tuple(sp.duplicate_symbols))
        out.append((n, sp))
    return out


def score_law(fold, pre, op, cd4, truth_cols, truth, law, gamma=None, alpha=None,
              cd4_f=None) -> dict:
    """204's `score()` with the transport law as an argument."""
    eff = pre * (PROD["alpha"] if alpha is None else alpha)
    whitened = op.apply(eff, fold["wcols"], PROD["gamma"] if gamma is None else gamma,
                        PROD["lam"])
    frac, per_unit = solve_fraction(whitened, cd4, PROD["cd4_f"] if cd4_f is None else cd4_f)
    final, _ = E062.add_post_whiten_residual(whitened, cd4, frac, fold["wcols"])
    final[:, fold["panel_positions"]] = 0.0
    out: dict = {}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        if mass is not None and law.kind == "rank":
            continue                       # no analytic density; a15 geometry only
        acc, per, perc = [], [], []
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            dl = L.deliver_law(final, b, gate, law, mass)
            acc.append(endpoints(dl[:, truth_cols], truth))
            p, c = _pds_per_target(dl[:, truth_cols], truth)
            per.append(p)
            perc.append(c)
        for k in acc[0]:
            out[f"{tag}_{k}"] = float(np.mean([a[k] for a in acc]))
        out[f"_pds_per_target_{tag}"] = np.mean(per, axis=0).tolist()
        out[f"_cos_per_target_{tag}"] = np.mean(perc, axis=0).tolist()
    for k in ("pds", "cos", "sgn", "sgn_w", "rr", "sse"):
        out[k] = out.get(f"a23_{k}", out.get(f"a15_{k}"))
    out["median_row_norm_final"] = float(np.median(np.linalg.norm(final, axis=1)))
    out["median_row_norm_pre"] = float(np.median(np.linalg.norm(pre, axis=1)))
    out["cd4_solved_fraction"] = float(frac)
    out["cd4_f_per_unit"] = float(per_unit)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--laws", nargs="*", default=["all"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--no-pm", action="store_true")
    ap.add_argument("--cd4-f", type=float, default=None)
    ap.add_argument("--alpha", type=float, default=None)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    laws = L.all_laws()
    keys = list(laws) if args.laws == ["all"] else args.laws
    tag = "nopm" if args.no_pm else "main"
    if args.cd4_f is not None:
        tag += f"_cd4f{args.cd4_f:g}"
    if args.alpha is not None:
        tag += f"_a{args.alpha:g}"
    out_path = Path(args.out) if args.out else WORK / f"plolo213_{tag}.json"
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["law"]) for r in rows}
    ops_seen = {}

    for held in args.folds:
        fold = build_fold(held)
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        cd4 = cd4_for(fold["targets"])
        op_path = Path(fold["op_path"] if FOLDS[held]["op_kind"] != "production"
                       else PROD_OPS / "whitener_union_pruned.npz")
        sha = sha256_file(op_path)
        full = OP_SHA_FULL.get(op_path.name)
        pref = OP_SHA_PREFIX.get(op_path.name)
        if (full is not None and sha != full) or (pref is not None
                                                  and not sha.startswith(pref)):
            raise SystemExit(f"operator sha drift on {op_path}: {sha}")
        op = Operator(op_path)
        ops_seen[held] = dict(path=str(op_path.relative_to(ROOT)), sha256=sha,
                              max_eig=float(op.evals.max()), min_eig=float(op.evals.min()),
                              n=int(op.evals.size))
        print(f"[{held}] operator {op_path.name} sha {sha[:16]} "
              f"eig [{op.evals.min():.3e}, {op.evals.max():.4f}] n={op.evals.size} "
              f"targets={len(fold['targets'])}", flush=True)
        for key in keys:
            if (held, key) in done:
                continue
            t0 = time.time()
            srcs = law_sources(held, fold, laws[key], use_pm=not args.no_pm)
            pre, aud = score204.weighted_baseline(srcs, fold["targets"], "equal",
                                                 {n: 1.0 for n, _ in srcs})
            e = score_law(fold, pre, op, cd4, fold["scored"], truth, laws[key],
                          cd4_f=args.cd4_f, alpha=args.alpha)
            rows.append(dict(fold=held, law=key, arm="equal", use_pm=not args.no_pm,
                             cd4_f=args.cd4_f, alpha_override=args.alpha,
                             n_targets=len(fold["targets"]), **laws[key].asdict(),
                             **e, operator=ops_seen[held], baseline_audit=aud,
                             seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=1))
            print(f"[{held}/{key}] a15cos {e['a15_cos']:.7f} "
                  f"a15pds {e['a15_pds']:.6f} "
                  f"a23cos {e.get('a23_cos', float('nan')):.7f} "
                  f"a23pds {e.get('a23_pds', float('nan')):.6f} "
                  f"q_w {e['a15_sgn_w']:.6f} f {e['cd4_solved_fraction']:.6f} "
                  f"({time.time()-t0:.0f}s)", flush=True)
            del srcs, pre
            gc.collect()
        del op, cd4
        gc.collect()

    if (not args.no_pm and args.cd4_f is None and args.alpha is None
            and set(args.folds) == set(FOLDS)):
        print("\n===== KNOWN-ANSWER ROW: logc20 and bc0 == 204's published `equal` row =====")
        ok = True
        for ident in ("logc20", "bc0"):
            if ident not in keys:
                continue
            for held, want in KA.items():
                got = [r for r in rows if r["fold"] == held and r["law"] == ident]
                if not got:
                    continue
                got = got[0]
                for k, v in want.items():
                    rel = abs(got[k] / v - 1.0)
                    good = rel <= 1e-12
                    ok &= good
                    print(f"  {ident:<8} {held:<15} {k:<10} {'PASS' if good else 'FAIL'}  "
                          f"{got[k]!r}  rel {rel:.2e}")
        (WORK / f"plolo213_ka_{tag}.json").write_text(json.dumps(
            {"passed": bool(ok), "operators": ops_seen, "laws": keys}, indent=1))
        if not ok:
            raise SystemExit("PLOLO KNOWN-ANSWER ROW FAILED — not 204's instrument. STOP.")
        print("  ==> PASSED.\n", flush=True)


if __name__ == "__main__":
    main()
