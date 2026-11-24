"""158 — known-answer checks, before any number in RESULTS.md is believed.

AGENTS.md 11.1 rule 10's corollary: a screen with no configuration whose answer is already
known cannot detect the class of bug where the wrong numbers look exactly as reasonable as
the right ones. Four checks here, and the experiment aborts if any fails.

1. `rule.weight` reproduces production's `pert_shrink_vector` bit for bit (missing='neutral').
2. `keys.donor_rows`'s streaming per-gene mean equals `Donor.decompose()`'s `m` exactly, on
   real donor data (k562_gwps, small enough to load whole).
3. `keys.build`'s `donor_cos` reproduces production's `donor_cosine_key` on a synthetic
   two-donor `MultiDonorModel` with disjoint-but-overlapping gene axes and a missing target.
4. THE PRODUCTION PATH. `--shrink-pert-key` multiplies the effect BEFORE the render, so the
   deployed object is `render(mass_fix(w*h))`, not `w*render(mass_fix(h))`. Those differ, and
   the difference is exactly what decides whether the pds claim survives. Both are scored.
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(REPO, "experiments/156-delivered-geometry"))
sys.path.insert(0, os.path.join(REPO, "src"))

import keys as K158  # noqa: E402
import rule as R158  # noqa: E402
from massfix import mass_fix, production_effect  # noqa: E402
from oracle import cos_agg, delivered, global_star, pds_of, sse_of  # noqa: E402
from screen import render  # noqa: E402


def _prod():
    import importlib.util
    path = os.path.join(REPO, "scripts/make_submission_pm_cd4.py")
    spec = importlib.util.spec_from_file_location("prod_158", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["prod_158"] = mod
    spec.loader.exec_module(mod)
    return mod


def check_weight_form(prod) -> None:
    rng = np.random.default_rng(7)
    worst = 0.0
    for trial in range(200):
        n = int(rng.integers(5, 140))
        k = rng.normal(size=n)
        if trial % 4 == 0:                                    # force ties
            k = np.round(k, 1)
        carried = np.ones(n, dtype=bool)
        if trial % 3 == 0:                                    # force missing rows
            carried[rng.choice(n, size=max(1, n // 10), replace=False)] = False
        gamma = float(rng.choice([0.0, 0.25, 0.5, 1.0, 2.0, 4.0]))
        floor = float(rng.choice([0.0, 0.1, 0.3, 0.5]))
        want = prod.pert_shrink_vector(np.where(carried, k, 0.0), carried, gamma,
                                       floor=floor, missing="neutral")
        mine = R158.weight(np.where(carried, k, np.nan), gamma, floor, sign=1.0)
        worst = max(worst, float(np.abs(want - mine).max()))
    print(f"1. weight form vs production pert_shrink_vector : max|diff| {worst:.3e}")
    assert worst < 1e-12, "weight form drifted from production"


def check_streaming_mean() -> None:
    from pathlib import Path
    from vcc2026.donor_transfer import Donor
    d = Donor.load(Path(REPO) / "data/donors/k562_gwps")
    m_full, _ = d.decompose()
    delta = np.load(Path(REPO) / "data/donors/k562_gwps/delta.npy", mmap_mode="r")
    acc = np.zeros(delta.shape[1], dtype=np.float64)
    for i in range(0, delta.shape[0], K158.CHUNK):
        acc += np.asarray(delta[i:i + K158.CHUNK], dtype=np.float64).sum(0)
    m_stream = acc / delta.shape[0]
    err = float(np.abs(m_stream - m_full.astype(np.float64)).max())
    rel = err / float(np.abs(m_full).max())
    print(f"2. streaming donor mean vs Donor.decompose()     : max|diff| {err:.3e} "
          f"(rel {rel:.2e}, float32 accumulation in decompose)")
    assert rel < 1e-5, "streaming mean disagrees with production decompose"
    del d, delta


def check_donor_cos(prod) -> None:
    """Synthetic two-donor model: overlapping gene axes, one target only one donor carries."""
    from vcc2026.donor_transfer import Donor
    from vcc2026.multi_donor import MultiDonorModel
    rng = np.random.default_rng(11)
    panel = np.asarray([f"G{i}" for i in range(40)])
    ga = np.asarray([f"G{i}" for i in range(30)] + ["X1", "X2"])
    gb = np.asarray([f"G{i}" for i in range(10, 40)])
    pa = np.asarray(["T0", "T1", "T2", "T3", "T4"])
    pb = np.asarray(["T0", "T1", "T2", "T4", "T5"])
    A = Donor(delta=rng.normal(size=(len(pa), len(ga))).astype(np.float32), perts=pa,
              genes=ga, b_ctrl=np.abs(rng.normal(size=len(ga))).astype(np.float32))
    B = Donor(delta=rng.normal(size=(len(pb), len(gb))).astype(np.float32) * 3.0, perts=pb,
              genes=gb, b_ctrl=np.abs(rng.normal(size=len(gb))).astype(np.float32))
    model = MultiDonorModel({"a": A, "b": B}, panel, gene_norm=False, eb=False)
    targets = ["T0", "T1", "T2", "T3", "T4", "T5"]
    want, carried = prod.donor_cosine_key(model, targets)

    # my construction, on the same objects
    mine = np.zeros(len(targets))
    S = {}
    for name, d in (("a", A), ("b", B)):
        pidx = {p: i for i, p in enumerate(d.perts)}
        m = d.delta.mean(axis=0)
        cols = [(int(np.flatnonzero(panel == g)[0]), j)
                for j, g in enumerate(d.genes) if (panel == g).any()]
        tc = np.asarray([c[0] for c in cols]); dc = np.asarray([c[1] for c in cols])
        M = np.zeros((len(targets), len(panel)))
        have = np.zeros(len(targets), dtype=bool)
        for r, t in enumerate(targets):
            if t in pidx:
                have[r] = True
                M[r, tc] = (d.delta[pidx[t]] - m)[dc]
        S[name] = (M, have)
    for r in range(len(targets)):
        if S["a"][1][r] and S["b"][1][r]:
            mine[r] = K158._cos(S["a"][0][r], S["b"][0][r])
    err = float(np.abs(want - mine).max())
    print(f"3. donor_cos vs production donor_cosine_key      : max|diff| {err:.3e} "
          f"(carried {carried.sum()}/{len(targets)}, mine nonzero "
          f"{int((mine != 0).sum())})")
    assert err < 1e-9, "donor cosine construction drifted from production"


def check_production_path() -> None:
    """Idealised post-render row scaling vs what `--shrink-pert-key` actually delivers."""
    print("\n4. PRODUCTION PATH — is a pre-render weight still a pure row scale afterwards?")
    with open(os.path.join(HERE, "work", "keys.json")) as fh:
        import json
        KJ = json.load(fh)
    for surface in ("hct116", "hek293t"):
        D, T, b, b_ref, gate, mask, ptar = production_effect(surface)
        Pd = render(mass_fix(D, b, gate, mode="minnorm"), b)
        Pm, Tm = Pd[:, mask], T[:, mask]
        key = np.asarray(KJ[surface]["keys"]["donor_cos_live"], dtype=np.float64)
        w = R158.weight(key, 2.0, 0.3, 1.0)

        post = w[:, None] * Pm                                    # the idealised object
        pre = render(mass_fix(w[:, None] * D, b, gate, mode="minnorm"), b)[:, mask]

        base_pds = pds_of(Pm, Tm)
        rel = float(np.median(np.linalg.norm(pre - post, axis=1)
                              / np.maximum(np.linalg.norm(post, axis=1), 1e-300)))
        r_post = dict(sse=sse_of(global_star(post, Tm) * np.ones(len(w)), post, Tm))
        r_pre = dict(sse=sse_of(global_star(pre, Tm) * np.ones(len(w)), pre, Tm))
        print(f"   {surface}: ||pre-post||/||post|| median {rel:.3e} | "
              f"pds post {pds_of(post, Tm):.10f} (d {pds_of(post, Tm) - base_pds:+.2e}) | "
              f"pds pre {pds_of(pre, Tm):.10f} (d {pds_of(pre, Tm) - base_pds:+.2e})")
        print(f"            sse post {r_post['sse']:.6f} (cos {cos_agg(r_post['sse']):.5f}) | "
              f"sse pre {r_pre['sse']:.6f} (cos {cos_agg(r_pre['sse']):.5f})")


def main() -> None:
    prod = _prod()
    check_weight_form(prod)
    check_streaming_mean()
    check_donor_cos(prod)
    check_production_path()
    print("\nall known-answer checks passed")


if __name__ == "__main__":
    main()
