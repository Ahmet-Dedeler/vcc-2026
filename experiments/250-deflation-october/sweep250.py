"""250 — one (panel, fold) per process: known-answer rows, then the k ladder, the controls,
the pinned A/B/C basis and the `+partial` composition.

    /usr/bin/time -l uv run python sweep250.py --panel real --fold k562_gwps

KNOWN-ANSWER ROWS (abort on failure, tolerance 1e-12 on cos_agg / cos / pds / sgn_w):
  BASE233   233's `build_final` (verbatim) on the panel's re-fitted operator == 233's stored
            BASE row (`score233.json`; 238's `six238.json` for efp/efr), AND 250's `compose`
            reproduces that matrix bit-for-bit (`np.array_equal`).
  null233   the deranged BASE == 233's stored `null` row (233's seed; informational on efp/efr,
            238 used the same seed and construction, six238.py:223-231).

SIZING: `allrows` (production A/B/C, f = 0.34) on `abc`; `carrying` (247's rule at g_cd4,
what `def_pipeline.DEF_ENV` ships on D/E/F) on every disjoint panel.  Every deflation arm is
RESOLVED (the fraction re-solved on the deflated object, 249's construction) unless its name
says `_fixed`.
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

TOL = 1e-12
KS = (2, 3, 5, 7, 10, 15)
KS_SPLIT = (3, 5, 10)
KS_PIN = (3, 5, 10)
REF = {
    "abc": L.ROOT / "experiments/233-october-retention/work/score233.json",
    "real": L.ROOT / "experiments/233-october-retention/work/score233.json",
    "pess": L.ROOT / "experiments/233-october-retention/work/score233.json",
    "efp": L.ROOT / "experiments/238-essential-free-panel/work/six238.json",
    "efr": L.ROOT / "experiments/238-essential-free-panel/work/six238.json",
}


def reference_rows(panel: str, held: str) -> tuple[dict, str]:
    rows = json.loads(REF[panel].read_text())
    ref = {r["arm"]: r for r in rows if r["panel"] == panel and r["fold"] == held
           and r["arm"] in ("BASE", "null")}
    return ref, str(REF[panel].relative_to(L.ROOT))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel", required=True, choices=list(L.B.BUNDLES))
    ap.add_argument("--fold", required=True, choices=list(L.B.FOLD_SOURCES))
    ap.add_argument("--sizing", default="auto", choices=("auto", "allrows", "carrying"))
    args = ap.parse_args()
    panel, held = args.panel, args.fold
    sizing = args.sizing if args.sizing != "auto" else ("allrows" if panel == "abc" else "carrying")
    t0 = time.time()
    gate_msg = wait_for_memory(tag=f" {panel}/{held}")
    print(f"[{panel}/{held}] gate clear: {gate_msg}", flush=True)

    st = L.stage(panel, held)
    targets = st["targets"]
    n = len(targets)
    cd4, whitened = st["cd4"], st["whitened"]
    ref, ref_file = reference_rows(panel, held)
    print(f"[{panel}/{held}] staged {n} targets, sizing={sizing}, op {st['opinfo']['path']} "
          f"axis {st['opinfo']['axis_sha256'][:12]}, union frac {st['frac_u']:.6f} "
          f"({time.time() - t0:.0f}s, peak {L.peak_gib():.2f} GiB)", flush=True)

    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    finals: dict[str, tuple[np.ndarray, float]] = {}
    fails = 0

    def run(name, cd4_res, frac_cd4, extra=(), audit=None, keep=False, ka=None,
            final_override=None):
        nonlocal fails
        t1 = time.time()
        if final_override is None:
            final, preserve = L.compose(st, cd4_res, frac_cd4, extra)
        else:
            final, preserve = final_override
        e = L.score(st, final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        row = dict(panel=panel, fold=held, arm=name, n_targets=n, sizing=sizing,
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   audit=dict(frac_cd4=float(frac_cd4), preserve=float(preserve),
                              **(audit or {})),
                   seconds=round(time.time() - t1, 1))
        msg = (f"[{panel}/{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} "
               f"pds {e['pds']:.6f} q_w {e['a15_sgn_w']:.6f} frac {frac_cd4:.5f} "
               f"({row['seconds']}s)")
        if ka is not None:
            r = ref.get(ka)
            if r is None:
                msg += f"   KA: no reference row {ka} in {ref_file}"
                row["ka"] = dict(ref=ka, source=ref_file, missing=True)
            else:
                d = max(abs(r["cos_agg"] - e["cos_agg"]), abs(r["cos"] - e["cos"]),
                        abs(r["pds"] - e["pds"]), abs(r["sgn_w"] - e["sgn_w"]))
                strict = True   # 238 derived its null with 233 seed 233233 (six238.py:223-231)
                ok = d < TOL
                row["ka"] = dict(ref=ka, source=ref_file, max_abs_diff=float(d), ok=bool(ok),
                                 strict=strict)
                msg += f"   KA {'PASS' if ok else 'FAIL'} vs {ka} ({d:.1e}){'' if strict else ' [informational]'}"
                if not ok and strict:
                    fails += 1
        rows.append(row)
        if keep:
            finals[name] = (final, preserve)
        print(msg, flush=True)
        return row

    # ------------------------------------------------------------- KNOWN-ANSWER ROWS
    frac_all, pu_all = L.solve_fraction(whitened, cd4, L.CD4_F)
    base233, pres233 = L.compose(st, cd4, frac_all)
    assert np.array_equal(base233, st["final233"]), "250 compose != 233 build_final (bitwise)"
    assert frac_all == st["aud233"]["frac_cd4"] and st["frac_u"] == st["aud233"]["frac_union"]
    print(f"[{panel}/{held}] KA compose==build_final bitwise PASS; frac_cd4 {frac_all:.9f} "
          f"frac_union {st['frac_u']:.9f}", flush=True)
    run("BASE233", cd4, frac_all, audit=dict(rule="allrows@0.34", per_unit=float(pu_all),
                                             sizing233=dict(st["aud233"])),
        final_override=(base233, pres233), ka="BASE", keep=True)
    run("null233", cd4, frac_all, audit=dict(rule="allrows@0.34", derange_seed=L.NULL_SEED),
        final_override=(L.derange(base233), pres233), ka="null")
    if fails:
        raise SystemExit(f"[{panel}/{held}] KNOWN-ANSWER FAILURE ({fails}) -- aborting")

    # ------------------------------------------------------------- the panel's stack
    frac0, pu0 = L.solve_cd4(whitened, cd4, sizing)
    norms0 = L.residual_norms(cd4)
    run("BASE", cd4, frac0, audit=dict(family="base", rule=sizing, per_unit=float(pu0),
                                       **norms0), keep=True)
    if sizing == "allrows":
        assert np.array_equal(finals["BASE"][0], base233), "BASE != BASE233 on allrows"
        fr_c, pu_c = L.solve_fraction_carrying(whitened, cd4, L.G_CD4)
        run("BASE_carrying", cd4, fr_c, audit=dict(family="base", rule="carrying@g_cd4",
                                                    per_unit=float(pu_c), **norms0))

    # ------------------------------------------------------------- Phase 1: the k ladder
    energy: dict[int, float] = {}
    d5 = aud5 = fr5 = None
    for k in KS:
        d, aud = L.deflate(cd4, "svd", k)
        aud.update(L.residual_norms(d))
        energy[k] = aud["energy_removed_frac"]
        fr, pu = L.solve_cd4(whitened, d, sizing)
        run(f"svd{k}", d, fr, audit=dict(family="ladder", sizing_mode="resolved",
                                          per_unit=float(pu), **aud), keep=(k == 5))
        if k == 5:
            d5, aud5, fr5 = d, aud, fr
    run("svd5_fixed", d5, frac0, audit=dict(family="ladder", sizing_mode="fixed", **aud5))
    other = "allrows" if sizing == "carrying" else "carrying"
    fro, puo = L.solve_cd4(whitened, d5, other)
    run(f"svd5_{other}", d5, fro, audit=dict(family="ladder", sizing_mode="resolved", rule=other,
                                              per_unit=float(puo), **aud5))

    # ------------------------------------------------------------- Phase 2: the controls
    for i in (0, 1):
        tag = "" if i == 0 else "b"
        d, aud, kk = L.rand_energy_matched(cd4, energy[5], L.RAND_SEED + 7 * i + 5)
        fr, pu = L.solve_cd4(whitened, d, sizing)
        run(f"randE5{tag}", d, fr, audit=dict(family="ctrl", control="random_energy_matched",
                                              matched_to_k=5, rand_dim=int(kk), seed=i,
                                              sizing_mode="resolved", per_unit=float(pu), **aud))
        run(f"randE5{tag}_fixed", d, frac0, audit=dict(family="ctrl",
                                                      control="random_energy_matched",
                                                      matched_to_k=5, rand_dim=int(kk), seed=i,
                                                      sizing_mode="fixed", **aud))
    d, aud = L.shuffled_projector(cd4, 5, L.SHUF_SEED + 5)
    fr, pu = L.solve_cd4(whitened, d, sizing)
    run("svd5_SHUFDEF", d, fr, audit=dict(family="ctrl", control="shuffled_projector",
                                          sizing_mode="resolved", per_unit=float(pu), **aud))
    for k in KS_SPLIT:
        for i in (0, 1):
            d, aud = L.split_half_deflate(cd4, k, L.RAND_SEED + 31 * i + k)
            fr, pu = L.solve_cd4(whitened, d, sizing)
            run(f"svd{k}_splithalf{'' if i == 0 else 'b'}", d, fr,
                audit=dict(family="ctrl", control="split_half_projector", seed=i,
                           sizing_mode="resolved", per_unit=float(pu), **aud))

    # ------------------------------------------------------------- Phase 3: the pinned A/B/C basis
    prod, prod_targets = L.production_cd4_abc()
    for k in KS_PIN:
        V, _U, baud = L.gene_basis(prod, k)
        d, aud = L.deflate_gene_basis(cd4, V, "ABCprod")
        fr, pu = L.solve_cd4(whitened, d, sizing)
        run(f"pinABC{k}", d, fr, audit=dict(family="pin", basis="abc_production_291rows",
                                            basis_uk_sha256=baud["basis_u_k_sha256"],
                                            basis_top_k_energy=baud["top_k_energy_frac_of_centred"],
                                            sizing_mode="resolved", per_unit=float(pu), **aud),
            keep=(k == 5))
        if k == 5:
            run("pinABC5_fixed", d, frac0, audit=dict(family="pin", basis="abc_production_291rows",
                                                      sizing_mode="fixed", **aud))
    V5own, _U5, baud5 = L.gene_basis(cd4, 5)
    V10own, _U10, baud10 = L.gene_basis(cd4, 10)

    # ------------------------------------------------------------- Phase 4: the partial rows
    part = L.cd4_partial(panel, targets)
    info = dict(part.condition_details[0])
    partial_rows = np.flatnonzero(L._carried(part)).tolist()
    if info["n_partial"] > 0:
        frp, pup = L.solve_fraction_carrying(whitened, part, 0.34)
        assert frp > 0, "VACUOUS partial-CD4 residual"
        pn = L.residual_norms(part)
        run("partial@0.34", cd4, frac0, extra=[(part, frp)],
            audit=dict(family="partial", dose=0.34, fraction_partial=float(frp),
                       per_unit_partial=float(pup), **info, **pn))
        run("svd5+partial@0.34", d5, fr5, extra=[(part, frp)],
            audit=dict(family="partial", k=5, dose=0.34, fraction_partial=float(frp),
                       sizing_mode="resolved", **info))
        ps = L.resid218.shuffled(part, seed=L.SHUF_SEED)
        frs, _ = L.solve_fraction_carrying(whitened, ps, 0.34)
        run("partial_SHUF@0.34", cd4, frac0, extra=[(ps, frs)],
            audit=dict(family="partial_ctrl", control="shuffled", dose=0.34,
                       fraction_partial=float(frs), **info))
    else:
        print(f"[{panel}/{held}] partial rows: VACUOUS (0 partial rows on this fold)", flush=True)

    # ------------------------------------------------------------- non-carrying rows, exactly
    nc = ~L._carried(cd4)
    noncarry = {}
    f0, s0 = finals["BASE"]
    for name in ("svd5", "pinABC5"):
        f1, s1 = finals[name]
        ratio = s1 / s0
        if nc.any():
            diff = float(np.abs(f1[nc] - ratio * f0[nc]).max())
            scale = float(np.abs(f0[nc]).max())
            noncarry[name] = dict(n_noncarry=int(nc.sum()), preserve_ratio=float(ratio),
                                  max_abs_diff=diff, max_rel=diff / scale if scale else 0.0)
        else:
            noncarry[name] = dict(n_noncarry=0)
    print(f"[{panel}/{held}] non-carrying rows: {json.dumps(noncarry)}", flush=True)

    meta = dict(panel=panel, fold=held, sizing=sizing, n_targets=n, targets=list(targets),
                carried_cd4=L._carried(cd4).tolist(), cd4_norm_pos=(np.linalg.norm(
                    np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0),
                    axis=1) > 0).tolist(),
                carried_union=L._carried(st["h1"]).tolist(), partial_rows=partial_rows,
                frac_u=st["frac_u"], pu_u=st["pu_u"], g_cd4=L.G_CD4, energy_at_k=energy,
                own_basis5=baud5, own_basis10=baud10, noncarrying=noncarry,
                opinfo=st["opinfo"], bind=st["bind"], fold_audit=st["fold"]["audit"],
                fold_cache=st["fold"]["cache"], gate=gate_msg, ka_failures=fails,
                seconds=round(time.time() - t0, 1), peak_rss_gib=L.peak_gib())
    (L.WORK / f"sweep250-{panel}-{held}.json").write_text(
        json.dumps(dict(meta=meta, rows=rows), indent=1))
    np.savez_compressed(L.WORK / f"pertarget250-{panel}-{held}.npz", **per_target)
    np.savez_compressed(L.WORK / f"basis250-{panel}-{held}.npz", V5=V5own, V10=V10own,
                        audit5=json.dumps(baud5), audit10=json.dumps(baud10))
    print(f"[{panel}/{held}] {len(rows)} arms in {time.time() - t0:.0f}s, peak RSS "
          f"{L.peak_gib():.2f} GiB, KA failures {fails}", flush=True)


if __name__ == "__main__":
    main()
