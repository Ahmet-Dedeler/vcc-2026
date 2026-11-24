"""290 step 3 — CONTROL C3 (PREREG §7): the 9 hallmark indicators replaced by random gene sets of
matched size drawn from the 2,739-gene axis, 20 draws, the rest of the feature block unchanged.

Reloads the 14 half-1 banks (unit rows only) because the projections must be recomputed per draw.

    PYTHONDONTWRITEBYTECODE=1 uv run python ctrl290.py

Writes work/ctrl290.json, work/ctrl290.log.
"""
from __future__ import annotations

import json
import resource
import time

import numpy as np
import pandas as pd

import lib290 as L

T0 = time.time()
LOG: list[str] = []
N_DRAW = 20


def P(msg: str = "") -> None:
    line = f"[{time.time() - T0:7.1f}s] {msg}" if msg else ""
    print(line, flush=True)
    LOG.append(line)


def main() -> None:
    h = L.check_prereg()
    P(f"prereg {h} MATCH")

    per_src, genes_held, keys, axis = L.scan_sources()
    lines = [per_src[k]["line"] for k in keys]
    nG, nA = len(genes_held), len(axis)
    h1 = [i for i, ln in enumerate(lines) if ln in L.HALF1_LINES]
    P(f"half-1 sources: {[keys[i] for i in h1]}")

    _, _, _, psizes = L.programme_directions(axis)
    sizes = [psizes[p]["n_on_axis"] for p in L.PROGRAMMES]
    prolif_size = len({g for p in L.PROLIF for g in L.hallmark_sets()[p]} & set(axis))
    P(f"matched sizes on the axis: {dict(zip([p.replace('HALLMARK_', '') for p in L.PROGRAMMES], sizes))}"
      f"; proliferation axis {prolif_size}")

    idxs, Us = {}, {}
    for i in h1:
        k = keys[i]
        idx, U, _, _, _ = L.load_unit(per_src[k]["rel"], axis, genes_held)
        idxs[k], Us[k] = idx, U
        P(f"  loaded {k:22s} {idx.size:6d} rows")

    # ---- the unit, target and folds (rebuilt exactly as fit290) --------------------------------
    held = L.load_held()
    X21 = held[L.FEATS].to_numpy(float)
    comp = np.load(L.SPLITS284, allow_pickle=True)["comp"]
    fold = L.component_folds(comp, 5)
    ref = np.load(L.REM284, allow_pickle=True)
    ok = np.asarray(ref["ok"])
    y = np.asarray(ref["r2"])
    u = np.flatnonzero(ok)
    yu, foldu, compu = y[u], fold[u], comp[u]
    feat = pd.read_csv(L.WORK / "feat290.csv")
    assert list(feat["gene"].astype(str)) == list(held["gene"].astype(str))
    P(f"unit {len(u)} genes; frozen remainder rho {L.spearman(np.asarray(ref['r1'])[u], yu):+.4f}")

    keep = ["logmag", "breadth", "gen_cos", "gen_abscos", "gen_pc5", "gen_pc20", "ontarget"]
    Xkeep = feat.loc[u, [f"{c}_h1" for c in keep]].to_numpy(float)
    on_axis = feat.loc[u, "on_axis"].to_numpy(float)[:, None]

    rng = np.random.default_rng(L.SEED + 91)
    rows = []
    for d in range(N_DRAW):
        Wr = np.zeros((len(L.PROGRAMMES), nA))
        for j, s in enumerate(sizes):
            Wr[j, rng.choice(nA, size=s, replace=False)] = 1.0
            Wr[j] /= np.linalg.norm(Wr[j])
        wp = np.zeros(nA)
        wp[rng.choice(nA, size=prolif_size, replace=False)] = 1.0
        wp /= np.linalg.norm(wp)

        PE = np.full((len(L.PROGRAMMES), len(h1), nG), np.nan, np.float32)
        PP = np.full((len(h1), nG), np.nan, np.float32)
        for q, i in enumerate(h1):
            k = keys[i]
            proj = Us[k] @ Wr.T.astype(np.float32)
            for j in range(len(L.PROGRAMMES)):
                PE[j, q, idxs[k]] = proj[:, j]
            PP[q, idxs[k]] = Us[k] @ wp.astype(np.float32)
        with np.errstate(invalid="ignore"):
            pe = np.stack([np.nanmean(PE[j], axis=0) for j in range(len(L.PROGRAMMES))], axis=1)
            pab = np.nanmean(np.abs(PP), axis=0)
        Xd = np.column_stack([Xkeep, pe[u], pab[u][:, None], on_axis])
        pr, _ = L.oof_ridge(Xd, yu, foldu)
        rho = L.spearman(yu, pr)
        rho_pab = L.spearman(pab[u], yu)
        rows.append(dict(draw=d, rho=rho, rho_prolif_abs=rho_pab))
        P(f"  draw {d:2d}: model rho {rho:+.4f}   random-prolif |proj| rho {rho_pab:+.4f} "
          f"({time.time() - T0:.0f}s)")

    r = np.array([x["rho"] for x in rows])
    rp = np.array([x["rho_prolif_abs"] for x in rows])
    fitj = json.loads((L.WORK / "fit290.json").read_text())
    obs = fitj["primary"]["component"]["response_ridge"]["rho"]["point"]
    obs_p = fitj["per_feature"]["pe_prolif_abs"]["point"]
    out = dict(prereg_sha256=h, n_draws=N_DRAW, matched_sizes=dict(zip(L.PROGRAMMES, sizes)),
               prolif_size=prolif_size, draws=rows,
               model=dict(observed=obs, mean=float(r.mean()), sd=float(r.std()),
                          lo=float(r.min()), hi=float(r.max()),
                          n_reaching=int((r >= obs).sum())),
               prolif_abs=dict(observed=obs_p, mean=float(rp.mean()), sd=float(rp.std()),
                               lo=float(rp.min()), hi=float(rp.max()),
                               n_at_least_as_negative=int((rp <= obs_p).sum())),
               seconds=time.time() - T0,
               peak_rss_gib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 30)
    P(f"\nCONTROL C3: model with random matched programmes {r.mean():+.4f} +- {r.std():.4f} "
      f"[{r.min():+.4f}, {r.max():+.4f}]; observed {obs:+.4f}; {int((r >= obs).sum())}/{N_DRAW} reach it")
    P(f"           random |proj| onto a matched random axis {rp.mean():+.4f} +- {rp.std():.4f} "
      f"[{rp.min():+.4f}, {rp.max():+.4f}]; observed pe_prolif_abs {obs_p:+.4f}; "
      f"{int((rp <= obs_p).sum())}/{N_DRAW} at least as negative")
    (L.WORK / "ctrl290.json").write_text(json.dumps(out, indent=1, default=float))
    (L.WORK / "ctrl290.log").write_text("\n".join(LOG))
    P(f"peak RSS {out['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    main()
