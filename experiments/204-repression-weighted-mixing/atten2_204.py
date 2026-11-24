"""204 Phase 1b — the attenuation estimate with GENE COVERAGE controlled, plus a
known-answer row for the estimator itself.

`atten204.py` measured `a_i` two ways and they disagreed in opposite directions.  One
confound is structural and was not controlled there: K562-gwps carries 7,679 of the
18,533 official genes and each Orion line carries 18,106.  Every inner product
`<d_i, d_X>` lives only on the genes BOTH carry, so a ratio `G_iX / G_jX` compares two
different gene sets, and 039 measured that the coverage difference alone is worth 0.019
macro PDS through a different mechanism.  Coverage is a separate lever from attenuation
and must not be smuggled into it.

So every quantity here is computed on ONE gene set: the columns carried by all three
production donors (and, for an instrument route, by the instrument too).

KNOWN-ANSWER ROW (rule 10's corollary).  The Gram/triple estimator is handed a synthetic
donor built from a real one with a KNOWN attenuation `c` and fresh independent noise; it
must return `c`.  A screen whose estimator has never been shown a right answer cannot
tell a real 1.7x from a bug.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(HERE))

from prep204 import DONORS, prepared, shared_targets  # noqa: E402
from plolo import DONOR_PATHS, PANEL_SET  # noqa: E402

EFFECT_SCALE = 6.817
SEED = 20260903


def triple(Gm: np.ndarray) -> np.ndarray:
    """a_i (relative, index 0 = 1) from the 3x3 Gram matrix of donor effects."""
    i, j, k = 0, 1, 2
    u = np.array([
        np.sqrt(Gm[i, j] * Gm[i, k] / Gm[j, k]),
        np.sqrt(Gm[j, i] * Gm[j, k] / Gm[i, k]),
        np.sqrt(Gm[k, i] * Gm[k, j] / Gm[i, j])])
    return u


def gram_matrix(E: list[np.ndarray], rows: np.ndarray, cols: np.ndarray) -> np.ndarray:
    n = len(E)
    Gm = np.zeros((n, n))
    sub = [e[np.ix_(rows, cols)] for e in E]
    for a in range(n):
        for b in range(a, n):
            Gm[a, b] = Gm[b, a] = float((sub[a] * sub[b]).sum())
    return Gm


def main() -> None:
    tg = shared_targets()
    panel_mask = np.array([t in PANEL_SET for t in tg])
    P = {d: prepared(d, tg, "shared") for d in DONORS}
    E = [P[d]["effects"] for d in DONORS]
    n_t = len(tg)

    # ---- the one gene set: columns every production donor actually measures -----------
    support = [np.abs(e).sum(0) > 0 for e in E]
    common = np.flatnonzero(support[0] & support[1] & support[2])
    each = {d: int(support[i].sum()) for i, d in enumerate(DONORS)}
    allcols = np.arange(E[0].shape[1])
    print(f"gene support: {each}, common {len(common)}", flush=True)

    out: dict = dict(n_targets=n_t, n_panel=int(panel_mask.sum()),
                     gene_support=each, n_common_genes=int(len(common)))

    # ---- production-frame row norms, the thing the median-norm rescale actually leaves -
    out["median_row_norm_production_frame"] = {}
    for i, d in enumerate(DONORS):
        out["median_row_norm_production_frame"][d] = dict(
            all_official=float(np.median(np.linalg.norm(E[i], axis=1))),
            common_genes=float(np.median(np.linalg.norm(E[i][:, common], axis=1))),
            raw_median_norm_own_axis=P[d]["raw_median_norm"])

    # ---- known-answer row: recover a planted attenuation -------------------------------
    rng = np.random.default_rng(SEED)
    ka = {}
    for c in (0.5, 0.8):
        # synthetic donor 3 = c * (shared component proxy) + fresh independent noise of
        # the same energy the real donor 3 carries.  The shared proxy is donor 0's row,
        # so the planted a_2/a_0 ratio is exactly `c` after matching row scale.
        S = E[0][:, common]
        noise = rng.normal(size=S.shape) * float(np.sqrt(
            ((E[2][:, common] ** 2).sum() - (E[0][:, common] * E[2][:, common]).sum() ** 2
             / (S ** 2).sum()) / S.size))
        synth = c * S + noise
        Gm = np.zeros((3, 3))
        mats = [E[0][:, common], E[1][:, common], synth]
        for a in range(3):
            for b in range(a, 3):
                Gm[a, b] = Gm[b, a] = float((mats[a] * mats[b]).sum())
        u = triple(Gm)
        ka[f"c={c}"] = dict(recovered=float(u[2] / u[0]), planted=c,
                            rel_err=float(u[2] / u[0] / c - 1.0))
        print(f"[known-answer] planted a_2/a_0 = {c} -> recovered {u[2]/u[0]:.4f} "
              f"({100*(u[2]/u[0]/c-1):+.2f}%)", flush=True)
    out["known_answer_planted_attenuation"] = ka

    # ---- route B on the common gene set, with a target bootstrap ----------------------
    routeB = {}
    for cols_name, cols in (("common_genes", common), ("all_official", allcols)):
        rec = {}
        for scope, m in (("shared", np.ones(n_t, bool)), ("panel", panel_mask)):
            rows = np.flatnonzero(m)
            Gm = gram_matrix(E, rows, cols)
            u = triple(Gm)
            a = u / u[0]
            tot = np.diag(Gm)
            noise = tot - u ** 2
            # per-target bootstrap
            boot = []
            rb = np.random.default_rng(SEED)
            sub = [e[np.ix_(rows, cols)] for e in E]
            pair = np.zeros((len(rows), 3, 3))
            for x in range(3):
                for y in range(3):
                    pair[:, x, y] = (sub[x] * sub[y]).sum(1)
            for _ in range(2000):
                idx = rb.integers(0, len(rows), len(rows))
                Gb = pair[idx].sum(0)
                ub = triple(Gb)
                boot.append(ub / ub[0])
            boot = np.array(boot)
            rec[scope] = dict(
                a_rel={d: float(a[i]) for i, d in enumerate(DONORS)},
                a_rel_ci={d: [float(np.percentile(boot[:, i], 2.5)),
                              float(np.percentile(boot[:, i], 97.5))]
                          for i, d in enumerate(DONORS)},
                snr={d: float(u[i] ** 2 / noise[i]) for i, d in enumerate(DONORS)},
                signal_energy={d: float(u[i] ** 2) for i, d in enumerate(DONORS)},
                noise_energy={d: float(noise[i]) for i, d in enumerate(DONORS)},
                w_minvar_rel={d: float((u[i] / noise[i]) / (u[0] / noise[0]))
                              for i, d in enumerate(DONORS)},
                gram=Gm.tolist())
            print(f"[triple/{cols_name}/{scope}] a_rel=" +
                  " ".join(f"{a[i]:.4f}" for i in range(3)) + "  snr=" +
                  " ".join(f"{u[i]**2/noise[i]:.4f}" for i in range(3)) + "  w_mv=" +
                  " ".join(f"{(u[i]/noise[i])/(u[0]/noise[0]):.4f}" for i in range(3)),
                  flush=True)
        routeB[cols_name] = rec
    out["routeB_triple_gene_controlled"] = routeB

    # ---- instrument route on the instrument's own common set --------------------------
    inst = {}
    for iname, ipath in (("kolf2", ROOT / "data/external/kolf2"),
                         ("h1_2025", ROOT / "data/donors/h1_2025"),
                         ("rpe1", ROOT / "data/donors/rpe1"),
                         ("nadig_jurkat", ROOT / "data/donors/nadig_jurkat"),
                         ("k562_essential", ROOT / "data/donors/k562_essential")):
        names = set(p.strip() for p in (ipath / "perts.txt").read_text().splitlines()
                    if p.strip())
        keep = [t for t in tg if t in names]
        if len(keep) < 40:
            continue
        DONOR_PATHS[iname] = ipath
        pi = prepared(iname, tuple(keep), f"inst-{iname}")
        X, ok = pi["effects"], pi["target_available"]
        isup = np.abs(X).sum(0) > 0
        cols = np.flatnonzero(support[0] & support[1] & support[2] & isup)
        rows = np.array([tg.index(t) for t in keep])[ok]
        Xs = X[ok][:, cols]
        g = np.array([float((E[i][rows][:, cols] * Xs).sum()) for i in range(3)])
        # bootstrap over targets
        per = np.array([[float((E[i][r][cols] * Xs[j]).sum()) for i in range(3)]
                        for j, r in enumerate(rows)])
        rb = np.random.default_rng(SEED)
        bt = np.array([(lambda s: s / s[0])(per[rb.integers(0, len(rows), len(rows))].sum(0))
                       for _ in range(2000)])
        inst[iname] = dict(
            n_targets=int(len(rows)), n_common_genes=int(len(cols)),
            a_rel={d: float(g[i] / g[0]) for i, d in enumerate(DONORS)},
            a_rel_ci={d: [float(np.percentile(bt[:, i], 2.5)),
                          float(np.percentile(bt[:, i], 97.5))]
                      for i, d in enumerate(DONORS)},
            same_line_as_k562=(iname == "k562_essential"))
        print(f"[inst/{iname}] n={len(rows)} genes={len(cols)} a_rel=" +
              " ".join(f"{g[i]/g[0]:.4f}" for i in range(3)) + "  ci_hct=" +
              f"[{np.percentile(bt[:,1],2.5):.3f},{np.percentile(bt[:,1],97.5):.3f}]"
              f" ci_hek=[{np.percentile(bt[:,2],2.5):.3f},{np.percentile(bt[:,2],97.5):.3f}]",
              flush=True)
    out["routeB_instrument_gene_controlled"] = inst

    (WORK / "atten2_204.json").write_text(json.dumps(out, indent=2))
    print("wrote work/atten2_204.json")


if __name__ == "__main__":
    main()
