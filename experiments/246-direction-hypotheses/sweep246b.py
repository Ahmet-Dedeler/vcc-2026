"""246 Phase 2b — the extension the first fold asked for.

Fold `k562_gwps` (sweep246.py) said: the deflated SELF-residual is negative on cosine on top
of the shipped stack at every (k, dose) -- the stack already carries that relaxation -- but
deflating the shipped CD4 residual's own shared subspace buys `pds` (+0.0094 at k = 5 against
+0.0002 / -0.0031 for two random 5-dim subspaces) at a cosine cost (-2.4 %).  This script
measures that properly:

  B2_cd4_svd<k>            k in 1..20, shipped fraction, one add, one preserve
  B2_cd4_svd<k>_resolved   same, CD4's fraction re-solved on the deflated object
  B2_cd4_randE<k>          ENERGY-MATCHED random control: a random subspace of the centred
                           row space whose dimension is chosen so the energy it removes
                           equals what svd<k> removes (two seeds).  209's control, matched.
  B2_cd4_svd<k>_SHUFDEF    the deflation fitted on a target-SHUFFLED copy of the residual and
                           its (row-permuted) projector applied to the real one -- removes a
                           k-dim subspace of the same rank fitted to the wrong rows.
  D_<mode><k>@0.1          MECHANISM ROW: the deflated self-residual on the BARE whitened base
                           (no CD4, no H1), against 218's `self0.1` (+5.35 % cos, -0.0097 pds
                           on the cancer folds).  Does deflation keep the cosine and drop the
                           `pds` cost when the relaxation is NOT already in the stack?
  C2_partial@g             family C again with the row-matched control built correctly: the
                           consensus `pre` on the partial rows' OWN gene support (sweep246's
                           `restrict_rows(self_res, ...)` was vacuous because the self-residual
                           is masked to CD4's complete-row support).
  E_cd4_svd<k>+partial     the two positives composed (k = 3, 5).

    uv run python sweep246b.py --fold <held>
"""
from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib246 as L  # noqa: E402

TOL = 1e-9
KS = (1, 2, 3, 4, 5, 7, 10, 15, 20)
KS_RESOLVED = (2, 3, 5, 7, 10)
KS_RANDE = (3, 5, 10)
KS_SHUFDEF = (3, 5, 10)
D_MODES = [("none", 0), ("centre", 0), ("svd", 1), ("svd", 2), ("svd", 5), ("svd", 10),
           ("svd", 20)]
DOSES_C = (0.17, 0.34)


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def rand_energy_matched(res, target_frac: float, seed: int, kmax: int = 260):
    """Random subspace of the centred carried-row space, dimension chosen by bisection so
    the removed energy fraction is >= target_frac (smallest such k)."""
    lo, hi = 1, kmax
    best = None
    while lo <= hi:
        mid = (lo + hi) // 2
        d, aud = L.deflate(res, "rand", mid, seed=seed)
        if aud["energy_removed_frac"] >= target_frac:
            best = (d, aud, mid)
            hi = mid - 1
        else:
            lo = mid + 1
    if best is None:
        d, aud = L.deflate(res, "rand", kmax, seed=seed)
        best = (d, aud, kmax)
    return best


def shuffled_projector(res, k: int, seed: int):
    """Fit the top-k projector on a row-DERANGED copy, apply it to the real residual."""
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(L._carried(res))
    R = eff[rows]
    mean = R.mean(axis=0, keepdims=True)
    Rc = R - mean
    Rc[~sup[rows]] = 0.0
    rng = np.random.default_rng(seed)
    perm = np.arange(rows.size)
    while True:
        rng.shuffle(perm)
        if not np.any(perm == np.arange(rows.size)):
            break
    # the projector fitted on the shuffled matrix is the SAME subspace (row order does not
    # change the row space), so the honest control is the projector's ROW-side action
    # U_k U_k^T permuted: apply the left projector of the shuffled fit to the un-shuffled rows.
    G = Rc[perm] @ Rc[perm].T
    w, U = np.linalg.eigh(G)
    Uk = U[:, np.argsort(w)[::-1][:k]]          # fitted on shuffled row order
    Rk = Rc - Uk @ (Uk.T @ Rc)                  # applied to the REAL row order
    Rk[~sup[rows]] = 0.0
    e0 = float((R * R).sum())
    eff[rows] = Rk
    audit = dict(mode="svd_shufdef", k=int(k), n_carried=int(rows.size),
                 energy_removed_frac=float(1.0 - (Rk * Rk).sum() / e0))
    return L._as_res(res, eff, "svd_shufdef", k), audit


def partial_row_matched_self(pre, part):
    """`pre` (the un-whitened consensus) on exactly the partial rows and their gene support."""
    votes = np.asarray(part.gene_votes)
    eff = np.where(votes > 0, np.asarray(pre, dtype=np.float64), 0.0)
    return L.E062.CD4Prepared(mode="self_on_partial_rows", mean_effect=eff,
                              gene_votes=votes.copy(),
                              target_condition_count=np.asarray(
                                  part.target_condition_count).copy(),
                              condition_details=part.condition_details)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    args = ap.parse_args()
    held = args.fold
    t0 = time.time()
    prev = json.loads((L.WORK / f"sweep246-{held}.json").read_text())
    prev_rows = {r["arm"]: r for r in prev["rows"]}

    st = L.stage(held)
    fold, targets = st["fold"], st["targets"]
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}

    def run(name, parts, extra=None, ka=None):
        t1 = time.time()
        final = L.compose(st, parts)
        e = L.score(st, final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        row = dict(fold=held, arm=name, n_targets=len(targets),
                   **{k: v for k, v in e.items() if not k.startswith("_")},
                   audit=extra or {}, seconds=round(time.time() - t1, 1))
        rows.append(row)
        msg = (f"[{held}/{name}] cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} "
               f"pds {e['pds']:.6f} q_w {e['a15_sgn_w']:.6f} ({row['seconds']}s)")
        if ka is not None:
            r = prev_rows[ka]
            d = max(abs(r["cos_agg"] - e["cos_agg"]), abs(r["a23_cos"] - e["a23_cos"]),
                    abs(r["pds"] - e["pds"]))
            assert d < TOL, f"KA FAILED {name} vs sweep246 {ka}: {d:.3e}"
            msg += f"   KA PASS vs {ka} ({d:.1e})"
        print(msg, flush=True)
        return row

    # known answers: this process must reproduce sweep246's STACK and B_cd4_svd5 rows
    run("STACK", L.stack_parts(st), dict(frac_cd4=st["frac_cd4"], frac_x=st["frac_x"]),
        ka="STACK")
    d5, a5 = L.deflate(st["cd4"], "svd", 5)
    run("B2_cd4_svd5", [(d5, st["frac_cd4"]), (st["union"], st["frac_x"])],
        dict(family="B2", sizing="fixed", **a5, **L.residual_norms(d5)), ka="B_cd4_svd5")

    # ---------------------------------------------------------------- B2: the k sweep
    energy_at_k = {5: a5["energy_removed_frac"]}
    for k in KS:
        if k == 5:
            continue
        d, aud = L.deflate(st["cd4"], "svd", k)
        energy_at_k[k] = aud["energy_removed_frac"]
        run(f"B2_cd4_svd{k}", [(d, st["frac_cd4"]), (st["union"], st["frac_x"])],
            dict(family="B2", sizing="fixed", **aud, **L.residual_norms(d)))
    for k in KS_RESOLVED:
        d, aud = L.deflate(st["cd4"], "svd", k)
        fr, per_unit = L.solve_fraction(st["whitened"], d, L.CD4_F)
        run(f"B2_cd4_svd{k}_resolved", [(d, fr), (st["union"], st["frac_x"])],
            dict(family="B2", sizing="resolved", fraction=float(fr),
                 per_unit=float(per_unit), **aud))
    for k in KS_RANDE:
        for seed_i in (0, 1):
            d, aud, kk = rand_energy_matched(st["cd4"], energy_at_k[k],
                                             L.RAND_SEED + 7 * seed_i + k)
            run(f"B2_cd4_randE{k}{'' if seed_i == 0 else 'b'}",
                [(d, st["frac_cd4"]), (st["union"], st["frac_x"])],
                dict(family="B2_ctrl", control="random_energy_matched", matched_to_k=k,
                     rand_dim=int(kk), seed=seed_i, sizing="fixed", **aud))
    for k in KS_SHUFDEF:
        d, aud = shuffled_projector(st["cd4"], k, L.SHUF_SEED + k)
        run(f"B2_cd4_svd{k}_SHUFDEF", [(d, st["frac_cd4"]), (st["union"], st["frac_x"])],
            dict(family="B2_ctrl", control="shuffled_projector", sizing="fixed", **aud))

    # ---------------------------------------------------------------- D: mechanism, bare base
    self_res = L.self_residual(st["cd4"], st["pre"])
    for mode, k in D_MODES:
        d, aud = L.deflate(self_res, mode, k)
        fr, _ = L.solve_fraction(st["whitened"], d, 0.10)
        run(f"D_{mode}{k if k else ''}@0.1", [(d, fr)],
            dict(family="D", dose=0.10, fraction=float(fr), **aud, **L.residual_norms(d)),
            ka=("self0.1_bare" if mode == "none" else None))
    for seed_i in (0, 1):
        d, aud, kk = rand_energy_matched(self_res, L.deflate(self_res, "svd", 5)[1][
            "energy_removed_frac"], L.RAND_SEED + 11 * seed_i)
        fr, _ = L.solve_fraction(st["whitened"], d, 0.10)
        run(f"D_randE5{'' if seed_i == 0 else 'b'}@0.1", [(d, fr)],
            dict(family="D_ctrl", control="random_energy_matched", rand_dim=int(kk),
                 dose=0.10, fraction=float(fr), **aud))
    ds = L.shuffled(L.deflate(self_res, "svd", 5)[0], seed=L.SHUF_SEED)
    fr, _ = L.solve_fraction(st["whitened"], ds, 0.10)
    run("D_svd5_SHUF@0.1", [(ds, fr)], dict(family="D_ctrl", control="shuffled", dose=0.10,
                                            fraction=float(fr)))

    # ---------------------------------------------------------------- C2 + E
    part = L.cd4_partial(targets, equalise=False)
    info = part.condition_details[0]
    if info["n_partial"] > 0:
        ctrl = partial_row_matched_self(st["pre"], part)
        for g in DOSES_C:
            fr, _ = L.solve_fraction_carrying(st["whitened"], part, g)
            run(f"C2_partial@{g:g}", L.stack_parts(st) + [(part, fr)],
                dict(family="C", dose=g, fraction=float(fr), **info), ka=f"C_partial@{g:g}")
            frc, _ = L.solve_fraction_carrying(st["whitened"], ctrl, g)
            assert frc > 0, "row-matched control is vacuous"
            run(f"C2_partial_ROWSELF@{g:g}", L.stack_parts(st) + [(ctrl, frc)],
                dict(family="C_ctrl", control="row_matched_self_on_support", dose=g,
                     fraction=float(frc), **info, **L.residual_norms(ctrl)))
        for k in (3, 5):
            d, aud = L.deflate(st["cd4"], "svd", k)
            fr, _ = L.solve_fraction_carrying(st["whitened"], part, 0.34)
            run(f"E_cd4_svd{k}+partial@0.34",
                [(d, st["frac_cd4"]), (st["union"], st["frac_x"]), (part, fr)],
                dict(family="E", k=k, dose=0.34, fraction=float(fr), **info))

    meta = dict(fold=held, n_targets=len(targets), energy_at_k=energy_at_k,
                seconds=round(time.time() - t0, 1), peak_rss_gib=peak_gib())
    (L.WORK / f"sweep246b-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows),
                                                              indent=1))
    np.savez_compressed(L.WORK / f"pertarget246b-{held}.npz", **per_target)
    print(f"[{held}] {len(rows)} arms in {time.time()-t0:.0f}s, peak RSS {peak_gib():.2f} GiB",
          flush=True)


if __name__ == "__main__":
    main()
