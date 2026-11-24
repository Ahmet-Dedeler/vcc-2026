"""246 Phase 2c — is CD4's shared subspace a STABLE property of the source, or fitted to the
rows it is applied to?  Split-half projector transfer: fit the top-k subspace on a random
half of CD4's carried rows and apply it to the OTHER half (and vice versa), so no row is
deflated by a basis fitted on itself.  If the `pds` gain survives, the subspace is a property
of CD4 (and transfers to an October panel refit); if it collapses, it is in-sample structure.

Also the CD4 dose-reduction comparison the mandate's `pds`-must-not-fall rule invites: the
same `pds` bought on the DOSE axis (216) costs 20x the cosine -- measured here on the same
folds at the same time: `cd4 f` in {0.28, 0.20} with the union at 0.28.

    uv run python sweep246c.py --fold <held>
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
KS = (2, 3, 5, 10)


def split_half_deflate(res, k: int, seed: int):
    eff = np.asarray(res.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(res.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(L._carried(res))
    R = eff[rows]
    rng = np.random.default_rng(seed)
    perm = rng.permutation(rows.size)
    halves = (perm[: rows.size // 2], perm[rows.size // 2:])
    e0 = float((R * R).sum())
    out = R.copy()
    for fit, apply in ((halves[0], halves[1]), (halves[1], halves[0])):
        Rf = R[fit]
        mean = Rf.mean(axis=0, keepdims=True)
        Rfc = Rf - mean
        Rfc[~sup[rows][fit]] = 0.0
        # right-singular basis of the FIT half:  V_k = Rfc^T U_k / s_k
        G = Rfc @ Rfc.T
        w, U = np.linalg.eigh(G)
        idx = np.argsort(w)[::-1][:k]
        Uk, sk = U[:, idx], np.sqrt(np.maximum(w[idx], 1e-300))
        Vk = (Rfc.T @ Uk) / sk[None, :]                     # (G, k), orthonormal columns
        Ra = R[apply] - mean
        Ra[~sup[rows][apply]] = 0.0
        Ra = Ra - (Ra @ Vk) @ Vk.T
        Ra[~sup[rows][apply]] = 0.0
        out[apply] = Ra
    eff[rows] = out
    audit = dict(mode="svd_splithalf", k=int(k), n_carried=int(rows.size),
                 energy_removed_frac=float(1.0 - (out * out).sum() / e0))
    return L._as_res(res, eff, "svd_splithalf", k), audit


def peak_gib():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fold", required=True)
    args = ap.parse_args()
    held = args.fold
    t0 = time.time()
    prev = {r["arm"]: r for r in json.loads(
        (L.WORK / f"sweep246-{held}.json").read_text())["rows"]}
    st = L.stage(held)
    targets = st["targets"]
    rows, per_target = [], {}

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
            r = prev[ka]
            d = max(abs(r["cos_agg"] - e["cos_agg"]), abs(r["a23_cos"] - e["a23_cos"]),
                    abs(r["pds"] - e["pds"]))
            assert d < TOL, f"KA FAILED {name} vs {ka}: {d:.3e}"
            msg += f"   KA PASS vs {ka} ({d:.1e})"
        print(msg, flush=True)

    run("STACK", L.stack_parts(st), ka="STACK")
    for k in KS:
        for seed_i in (0, 1):
            d, aud = split_half_deflate(st["cd4"], k, L.RAND_SEED + 31 * seed_i + k)
            run(f"B3_cd4_svd{k}_splithalf{'' if seed_i == 0 else 'b'}",
                [(d, st["frac_cd4"]), (st["union"], st["frac_x"])],
                dict(family="B3", control="split_half_projector", seed=seed_i,
                     sizing="fixed", **aud))
    # the dose axis, same folds, same base, for the exchange-rate comparison
    for f in (0.28, 0.20):
        fr, _ = L.solve_fraction(st["whitened"], st["cd4"], f)
        run(f"DOSE_cd4@{f:g}", [(st["cd4"], fr), (st["union"], st["frac_x"])],
            dict(family="DOSE", cd4_f=f, fraction=float(fr)))

    (L.WORK / f"sweep246c-{held}.json").write_text(json.dumps(
        dict(meta=dict(fold=held, seconds=round(time.time() - t0, 1),
                       peak_rss_gib=peak_gib()), rows=rows), indent=1))
    np.savez_compressed(L.WORK / f"pertarget246c-{held}.npz", **per_target)
    print(f"[{held}] {len(rows)} arms in {time.time()-t0:.0f}s, peak {peak_gib():.2f} GiB",
          flush=True)


if __name__ == "__main__":
    main()
