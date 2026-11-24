"""250 — the A/B/C bridge between the two instruments.

(a) 246's instrument exactly (062's held operators, 246's own fold caches, the shipped CD4
    banks, today's panel) driven through 250's functions: STACK, `B2_cd4_svd5` (fixed),
    `B2_cd4_svd5_resolved`, `B2_cd4_randE5` (fixed), `B2_cd4_svd5_SHUFDEF` (fixed),
    `B3_cd4_svd5_splithalf` (fixed), `C2_partial@0.34`, `E_cd4_svd5+partial@0.34` must
    reproduce 246's json rows to 1e-9 (expected 0.0e+00) -- so every deflation / control /
    composition function in `lib250` is asserted bit-faithful to 246's before it is pointed at
    a disjoint panel.  Abort on failure.
(b) the SAME objects through 233's re-fitted A/B/C operator (essential-pruned, 12,000-row cap:
    the recipe every disjoint panel's operator is fitted by), so the difference between the
    two instruments on today's panel is isolated to the operator recipe.

    /usr/bin/time -l uv run python bridge250.py --fold k562_gwps      (also orion_*, kolf2)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gate250 import wait_for_memory  # noqa: E402
import lib250 as L  # noqa: E402
import plolo  # noqa: E402
from plolo import Operator, cd4_for, solve_fraction  # noqa: E402

TOL = 1e-9
HELD_OPS = L.ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners"
PROD_OP = L.ROOT / "data/donors/whitener_union_pruned.npz"
CACHE246 = L.ROOT / "experiments/246-direction-hypotheses/work/plolo"
SHIPPED_CD4 = L.ROOT / "data/external/cd4t/cloud-compact/cd4_compact"
W246 = L.ROOT / "experiments/246-direction-hypotheses/work"


def ref246(held: str) -> dict:
    rows = {}
    for tag in ("sweep246", "sweep246b", "sweep246c"):
        d = json.loads((W246 / f"{tag}-{held}.json").read_text())
        for r in d["rows"]:
            rows.setdefault(r["arm"], r)
    return rows


def stage246(held: str):
    """246's staging: today's panel, plolo's fold cache (246's byte copy), the shipped banks."""
    abc = L.B.read_panel("abc")
    plolo.PANEL = abc
    plolo.PANEL_SET = set(abc)
    plolo.CD4_ROOT = SHIPPED_CD4
    plolo.WORK = CACHE246                                    # read-only cache hit
    if held == "kolf2":
        plolo.FOLDS["kolf2"] = dict(sources=L.B.FOLD_SOURCES["kolf2"], op=PROD_OP,
                                    op_kind="production")
    else:
        plolo.FOLDS[held]["op"] = HELD_OPS / f"held_{held}.npz"
        plolo.FOLDS[held]["op_kind"] = "held"
        plolo.FOLDS[held]["sources"] = L.B.FOLD_SOURCES[held]
    L.resid218.PANEL_SET = set(abc)
    L.resid218._PREP_CACHE.clear()
    before = sorted(p.name for p in CACHE246.iterdir())
    fold = plolo.build_fold(held)
    assert sorted(p.name for p in CACHE246.iterdir()) == before, "wrote into 246's cache"
    plolo.WORK = L.S.SHARED                                  # bctrl caches: 250's copy
    targets = fold["targets"]
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
    pre = np.asarray(fold["pre"], dtype=np.float64)          # stored float32 by plolo
    cd4 = cd4_for(targets)
    h1 = L.union_rows([L.resid218.source_residual(n, targets, "gene_norm") for n in L.UNION])
    return dict(fold=fold, targets=targets, truth=truth, pre=pre, cd4=cd4, h1=h1,
                wcols=np.asarray(fold["wcols"]), panel_positions=fold["panel_positions"])


def whiten_with(base: dict, op_path: Path) -> dict:
    op = Operator(op_path)
    whitened = op.apply(base["pre"] * L.ALPHA, base["wcols"], L.GAMMA, L.LAM)
    del op
    frac_u, pu_u = L.solve_fraction_carrying(whitened, base["h1"], L.UNION_G)
    st = dict(base)
    st.update(whitened=whitened, frac_u=float(frac_u), pu_u=float(pu_u), op=str(op_path))
    return st


def arms(st: dict, ref: dict | None, label: str, held: str, part) -> list[dict]:
    out = []
    cd4, whitened = st["cd4"], st["whitened"]
    frac_cd4, _ = solve_fraction(whitened, cd4, L.CD4_F)

    def run(name, cd4_res, frac, extra=(), audit=None, ka=None):
        t1 = time.time()
        final, preserve = L.compose(st, cd4_res, frac, extra)
        e = L.score(st, final)
        row = dict(fold=held, operator=label, arm=name, n_targets=len(st["targets"]),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   audit=dict(frac_cd4=float(frac), preserve=preserve, **(audit or {})),
                   seconds=round(time.time() - t1, 1))
        msg = (f"[{held}/{label}/{name}] cos_agg {e['cos_agg']:.9f} macro {e['a23_cos']:.9f} "
               f"pds {e['pds']:.9f} q_w {e['a15_sgn_w']:.6f}")
        if ka is not None and ref is not None:
            r = ref[ka]
            d = max(abs(r["cos_agg"] - e["cos_agg"]), abs(r["a23_cos"] - e["a23_cos"]),
                    abs(r["pds"] - e["pds"]), abs(r["a15_sgn_w"] - e["a15_sgn_w"]))
            row["ka"] = dict(ref=ka, max_abs_diff=float(d), ok=bool(d < TOL))
            msg += f"   KA {'PASS' if d < TOL else 'FAIL'} vs 246 {ka} ({d:.1e})"
            assert d < TOL, f"KA FAILED {name} vs 246 {ka}: {d:.3e}"
        out.append(row)
        print(msg, flush=True)

    run("STACK", cd4, frac_cd4, ka="STACK")
    d5, a5 = L.deflate(cd4, "svd", 5)
    run("B2_cd4_svd5", d5, frac_cd4, audit=dict(sizing="fixed", **a5), ka="B2_cd4_svd5")
    fr5, pu5 = solve_fraction(whitened, d5, L.CD4_F)
    run("B2_cd4_svd5_resolved", d5, fr5, audit=dict(sizing="resolved", per_unit=float(pu5), **a5),
        ka="B2_cd4_svd5_resolved")
    for k in (3, 10):
        d, a = L.deflate(cd4, "svd", k)
        fr, pu = solve_fraction(whitened, d, L.CD4_F)
        run(f"B2_cd4_svd{k}_resolved", d, fr, audit=dict(sizing="resolved", **a),
            ka=f"B2_cd4_svd{k}_resolved")
    d, a, kk = L.rand_energy_matched(cd4, a5["energy_removed_frac"], L.RAND_SEED + 5)
    run("B2_cd4_randE5", d, frac_cd4, audit=dict(sizing="fixed", rand_dim=int(kk), **a),
        ka="B2_cd4_randE5")
    d, a = L.shuffled_projector(cd4, 5, L.SHUF_SEED + 5)
    run("B2_cd4_svd5_SHUFDEF", d, frac_cd4, audit=dict(sizing="fixed", **a),
        ka="B2_cd4_svd5_SHUFDEF")
    d, a = L.split_half_deflate(cd4, 5, L.RAND_SEED + 5)
    run("B3_cd4_svd5_splithalf", d, frac_cd4, audit=dict(sizing="fixed", **a),
        ka="B3_cd4_svd5_splithalf")
    if part is not None and part.condition_details[0]["n_partial"] > 0:
        frp, _ = L.solve_fraction_carrying(whitened, part, 0.34)
        run("C2_partial@0.34", cd4, frac_cd4, extra=[(part, frp)],
            audit=dict(fraction_partial=float(frp)), ka="C2_partial@0.34")
        run("E_cd4_svd5+partial@0.34", d5, frac_cd4, extra=[(part, frp)],
            audit=dict(sizing="fixed", fraction_partial=float(frp)),
            ka="E_cd4_svd5+partial@0.34")
    # the pinned production basis on 246's own folds (in-sample reference for Phase 3)
    prod, _ = L.production_cd4_abc()
    V5, _U, b5 = L.gene_basis(prod, 5)
    d, a = L.deflate_gene_basis(cd4, V5, "ABCprod")
    fr, pu = solve_fraction(whitened, d, L.CD4_F)
    run("pinABC5_resolved", d, fr, audit=dict(sizing="resolved", **a))
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True, choices=("k562_gwps", "orion_hct116",
                                                     "orion_hek293t", "kolf2"))
    args = ap.parse_args()
    held = args.fold
    t0 = time.time()
    print(f"[{held}] gate clear: {wait_for_memory(tag=f' bridge/{held}')}", flush=True)
    ref = ref246(held)
    base = stage246(held)
    part = L.cd4_partial("abc", base["targets"], cd4_root=SHIPPED_CD4)
    rows = []
    # (a) 246's operator
    op246 = PROD_OP if held == "kolf2" else HELD_OPS / f"held_{held}.npz"
    st = whiten_with(base, op246)
    rows += arms(st, ref, "op062" if held != "kolf2" else "production", held, part)
    del st
    # (b) 233's re-fitted A/B/C operator (238's re-fit of the production one for kolf2)
    op233 = L.B.op_path("abc", held)
    st = whiten_with(base, op233)
    rows += arms(st, None, "refit233", held, part)
    del st
    (L.WORK / f"bridge250-{held}.json").write_text(json.dumps(
        dict(meta=dict(fold=held, n_targets=len(base["targets"]), op246=str(op246),
                       op233=str(op233), seconds=round(time.time() - t0, 1),
                       peak_rss_gib=L.peak_gib()), rows=rows), indent=1))
    print(f"[{held}] bridge done: {len(rows)} rows in {time.time() - t0:.0f}s, "
          f"peak {L.peak_gib():.2f} GiB", flush=True)


if __name__ == "__main__":
    main()
