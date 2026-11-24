"""158 — build-time-computable per-target keys.

Every key here is a function of DONORS, the recipient context's CONTROL profile, our own
predicted effect, and frozen external tables. None of them touches the held surface's
perturbation truth. The two Orion surfaces use disjoint donor sets
(`experiments/104-context-heldout-credit-adapter/work/<surface>/manifest.json`,
`baseline.sources`), so a key computed for HCT116 never sees HCT116 responses:

    hct116  <- k562_gwps + orion_hek293t
    hek293t <- k562_gwps + orion_hct116

DONOR PROCESSING MUST MATCH THE BUILD. `scripts/make_submission_pm_cd4.py` constructs
`MultiDonorModel(donors, spec.genes, gene_norm=False, eb=False)`, so the donor vector for
target `p` is `delta[p] - mean_over_all_perts(delta)` mapped onto the 18,533-gene panel axis
with zeros off the donor's own axis -- exactly what `donor_cosine_key` and `donor_agreement`
consume. `MultiDonorModel._renorm` cancels out of every pairwise cosine and out of every
ratio used below, so it is not applied; the geometric-mean keys ARE renorm-sensitive and are
built with it explicitly.

Memory: Orion's `delta.npy` is 2.7 GB. It is memory-mapped and the per-gene mean is
accumulated in row chunks, so peak RSS stays near 0.5 GB rather than the 5.4 GB a
`Donor.load()` + `decompose()` would cost.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
DONORS = os.path.join(REPO, "data/donors")
PM = os.path.join(REPO, "data/external/perturbmatch")
FROZEN = os.path.join(REPO, "experiments/104-context-heldout-credit-adapter/work/%s/frozen.npz")

# from each surface's manifest, key `baseline.sources`
SOURCES = {
    "hct116": ("k562_gwps", "orion_hek293t"),
    "hek293t": ("k562_gwps", "orion_hct116"),
}
CHUNK = 2000


def _lines(path: Path) -> np.ndarray:
    """`Donor.load`'s reader, verbatim, so indices line up with production's."""
    return np.asarray(path.read_text().split("\n"))


def donor_rows(name: str, targets: list[str], panel_genes: np.ndarray):
    """`s[p]` for each requested target, on the panel gene axis, zeros off the donor axis.

    Returns `(S, carried, meta)` where `S` is (n_targets, n_panel_genes), `carried` marks
    targets the donor screened, and `meta` carries the per-target own-gene knockdown depth
    plus the donor-axis effect norm (both computed on the donor's OWN axis, before the
    panel projection loses genes).
    """
    d = Path(DONORS) / name
    perts = _lines(d / "perts.txt")
    genes = _lines(d / "genes.txt")
    delta = np.load(d / "delta.npy", mmap_mode="r")
    assert delta.shape == (len(perts), len(genes)), (delta.shape, len(perts), len(genes))

    # streaming per-gene mean over ALL perturbations -- `Donor.decompose`'s `m`
    acc = np.zeros(delta.shape[1], dtype=np.float64)
    for i in range(0, delta.shape[0], CHUNK):
        acc += np.asarray(delta[i:i + CHUNK], dtype=np.float64).sum(0)
    m = acc / delta.shape[0]

    pindex = {str(p): i for i, p in enumerate(perts)}
    gpos = {str(g): i for i, g in enumerate(genes)}
    tcols, dcols = [], []
    for j, g in enumerate(genes):
        t = np.flatnonzero(panel_genes == g)
        if t.size:
            tcols.append(int(t[0])); dcols.append(j)
    tcols = np.asarray(tcols, dtype=np.int64)
    dcols = np.asarray(dcols, dtype=np.int64)

    S = np.zeros((len(targets), len(panel_genes)), dtype=np.float64)
    carried = np.zeros(len(targets), dtype=bool)
    kd = np.full(len(targets), np.nan)
    dnorm = np.full(len(targets), np.nan)
    for r, t in enumerate(targets):
        i = pindex.get(str(t))
        if i is None:
            continue
        carried[r] = True
        s_full = np.asarray(delta[i], dtype=np.float64) - m
        S[r, tcols] = s_full[dcols]
        dnorm[r] = float(np.linalg.norm(s_full))
        own = gpos.get(str(t))
        if own is not None:
            kd[r] = float(s_full[own])

    # MultiDonorModel._renorm's numerator input: median on-axis effect norm
    scale = float(np.median(np.linalg.norm(np.asarray(delta[::37])[:, dcols]
                                          - m[dcols], axis=1)))
    return S, carried, dict(kd=kd, dnorm=dnorm, scale=scale,
                            n_perts=int(delta.shape[0]), n_axis_shared=int(len(dcols)))


def _cos(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(a @ b) / den if den > 1e-12 else 0.0


def build(surface: str, targets: list[str], panel_genes: np.ndarray,
          gate: np.ndarray, mask: np.ndarray, ctrl_cpm: np.ndarray,
          bctrl: np.ndarray, delivered: np.ndarray) -> dict[str, np.ndarray]:
    """Every candidate key for one surface. `delivered` is our a23 effect on `mask`."""
    names = SOURCES[surface]
    S, C, M = {}, {}, {}
    for n in names:
        S[n], C[n], M[n] = donor_rows(n, targets, panel_genes)
    renorm_ref = float(np.median([M[n]["scale"] for n in names]))
    renorm = {n: renorm_ref / (M[n]["scale"] or 1.0) for n in names}

    n_t = len(targets)
    a, b = names
    scored = np.zeros(len(panel_genes), dtype=bool)
    scored[mask] = True
    live = gate & scored                    # the coordinates production lets move AND scores

    out: dict[str, np.ndarray] = {}

    # ---- whole-vector cross-donor DIRECTION agreement (production's `donor_cosine_key`)
    out["donor_cos"] = np.asarray([_cos(S[a][i], S[b][i]) for i in range(n_t)])
    out["donor_cos_live"] = np.asarray(
        [_cos(S[a][i][live], S[b][i][live]) for i in range(n_t)])
    both = (np.abs(S[a]).sum(0) > 0) & (np.abs(S[b]).sum(0) > 0)
    out["donor_cos_shared"] = np.asarray(
        [_cos(S[a][i][both], S[b][i][both]) for i in range(n_t)])

    # ---- per-gene `donor_agreement` reduced to one number per target, three ways
    A = np.zeros((n_t, len(panel_genes)))
    for i in range(n_t):
        st = np.vstack([renorm[a] * S[a][i], renorm[b] * S[b][i]])
        unan = (np.abs(np.sign(st).sum(0)) == 2) & np.all(st != 0, axis=0)
        geo = np.exp(np.log(np.abs(st) + 1e-300).mean(0))
        A[i] = np.where(unan, geo, 0.0)
    Alive = A[:, live]
    out["agree_frac"] = (Alive > 0).mean(1)                     # how much of the row agrees
    out["agree_norm"] = np.linalg.norm(Alive, axis=1)           # how strongly it agrees
    mean_mag = np.abs(np.vstack([[renorm[a] * S[a][i], renorm[b] * S[b][i]]
                                 for i in range(n_t)]).reshape(n_t, 2, -1)).mean(1)[:, live]
    out["agree_share"] = Alive.sum(1) / np.maximum(mean_mag.sum(1), 1e-300)

    # ---- donor magnitude structure
    na = np.linalg.norm(S[a][:, live], axis=1)
    nb = np.linalg.norm(S[b][:, live], axis=1)
    out["donor_norm_geo"] = np.sqrt(np.maximum(renorm[a] * na * renorm[b] * nb, 0.0))
    out["donor_norm_ratio"] = (np.minimum(renorm[a] * na, renorm[b] * nb)
                               / np.maximum(np.maximum(renorm[a] * na, renorm[b] * nb), 1e-300))

    # ---- the target's own measured knockdown depth in each donor (more negative = stronger)
    kd = np.vstack([M[a]["kd"], M[b]["kd"]])
    out["kd_worst"] = -np.nanmin(-kd, axis=0)     # the shallower of the two knockdowns
    out["kd_mean"] = np.nanmean(kd, axis=0)
    for j, n in enumerate(names):
        out[f"kd_{n}"] = kd[j]

    # ---- recipient-context basal expression of the target gene itself
    own = np.asarray([int(np.flatnonzero(panel_genes == t)[0])
                      if (panel_genes == t).any() else -1 for t in targets])
    cpm = np.where(own >= 0, ctrl_cpm[np.maximum(own, 0)], np.nan)
    out["basal_log_cpm"] = np.log10(np.maximum(cpm, 1e-3))

    # ---- our own predicted object (build-time by construction)
    out["row_norm"] = np.linalg.norm(delivered, axis=1)
    out["row_nnz"] = (delivered != 0).sum(1).astype(float)

    # ---- external support flags
    pm_perts = {str(x) for x in _lines(Path(PM) / "perts.txt")}
    out["pm_support"] = np.asarray([1.0 if str(t) in pm_perts else 0.0 for t in targets])

    # ---- experiment 104's context-held-out fitted credit score, reduced per target
    z = np.load(FROZEN % surface, allow_pickle=True)
    for tag in ("score_agreement", "score_full"):
        sc = z[tag].astype(np.float64)[:, live]
        out[f"{tag}_mean"] = sc.mean(1)
        w = np.abs(delivered)
        wl = w / np.maximum(w.sum(1, keepdims=True), 1e-300)
        # weight by our own amplitude: the credit that actually gets delivered
        full = np.zeros((n_t, len(panel_genes)))
        full[:, live] = sc
        out[f"{tag}_wmean"] = (wl * full[:, mask]).sum(1)

    meta = {n: {k: v for k, v in M[n].items() if k not in ("kd", "dnorm")} for n in names}
    meta["renorm"] = renorm
    return out, meta


def main() -> None:
    import sys
    sys.path.insert(0, HERE)
    from oracle import SURFACES, delivered as load_delivered

    work = os.path.join(HERE, "work")
    os.makedirs(work, exist_ok=True)
    store = {}
    for surface in SURFACES:
        P, T, ptar, aux = load_delivered(surface)
        z = np.load(FROZEN % surface, allow_pickle=True)
        genes = np.asarray([str(x) for x in z["genes"]])
        k, meta = build(surface, [str(x) for x in ptar], genes, aux["gate"],
                        aux["mask"], z["ctrl_cpm"].astype(np.float64),
                        z["bctrl"].astype(np.float64), P)
        store[surface] = {"targets": [str(x) for x in ptar],
                         "keys": {name: v.tolist() for name, v in k.items()},
                         "meta": meta}
        print(f"{surface}: {len(k)} keys over {len(ptar)} targets; donors {SOURCES[surface]}")
        for name, v in k.items():
            finite = np.isfinite(v)
            print(f"   {name:24s} finite {finite.sum():3d}/{len(v)}  "
                  f"[{np.nanmin(v):+.4g}, {np.nanmax(v):+.4g}]  "
                  f"uniq {len(np.unique(v[finite])):3d}")
    with open(os.path.join(work, "keys.json"), "w") as fh:
        json.dump(store, fh, indent=1)
    print(f"\nwrote {os.path.join(work, 'keys.json')}")


if __name__ == "__main__":
    main()
