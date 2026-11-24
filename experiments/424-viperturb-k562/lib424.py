"""424 — shared machinery: VIPerturb-seq K562 (data/donors/viperturb_k562) as a donor in 420's bulk instruments.

Everything is 420's (lib420 -> lib418 -> sens415 / plolo / 414): the PP composition, 410's key clip (0.3, 3) held fixed,
417a's sharpening q 1.25, 414's ballast realisation at KMAX 16, the a23 read-back with the board-emulated pds. What is
new here is only how the donor BASE (`pre`, before whitening) is formed:

  base  (417a)  normrel-cap(-1.5) mean over the fold's sources, PM50 correction at the uniform K562 share
  (i)   extra   VIPerturb as one more source at equal weight in that mean (normrel median over its own panel rows)
  (ii)  avg     one K562 row: K562 slot := 0.5 (K562-GWPS + PM50 corr) + 0.5 VIPerturb where both carry the target;
                VIPerturb alone where only it carries the target (the fill); K562 slot's normrel median recomputed
  (ii-) avgonly (ii) without the fill (the decomposition row: (ii) = (ii-) + fill, to first order)
  (ii') ivw     (ii) with cell-count weights n_vip / (n_vip + n_gwps) instead of 0.5 (inverse-variance)
  (iii) fill    K562 slot := VIPerturb only on targets K562-GWPS lacks (slot's normrel median recomputed, as a
                builder would: this moves the cap weight of other K562 rows slightly)
  (iii0) fillpure  (iii) with the K562 slot's normrel median held at its K562-GWPS value: only the fill rows change
  CTRL  shuf    (i) with the VIPerturb rows permuted across targets (must lose)

VIPerturb rows are prepared exactly as every production donor is (E062.prepare_main_source / sens415.full_rows): the
donor's own genome-wide column mean (all 17,489 perts) subtracted, scaled by EFFECT_SCALE / the median centred row norm,
mapped to the official axis, panel columns kept (normrel norms include them). The operator, CD4 / H1 / partial residuals
and the key are the fold's (held fixed). The K562-GWPS fold holds K562 out as truth, so VIPerturb (also K562) is
excluded there by the leakage rule: every arm equals the base on that fold by construction.

Nothing here writes outside experiments/424-viperturb-k562 (bytecode off; fold caches are APFS clones in ./work/plolo).
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import json  # noqa: E402
import os  # noqa: E402
import shutil  # noqa: E402
import subprocess  # noqa: E402
from pathlib import Path  # noqa: E402
from types import SimpleNamespace  # noqa: E402

import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
W420 = ROOT / "experiments/420-bulk-stack"
sys.path.insert(0, str(W420))
import lib420  # noqa: E402
from lib420 import (ALPHA, E062, FOLDS, GAMMA, LAM, PROD_OPS, UNION, W414, W418, Donors, Operator, betas_for,  # noqa: E402,F401
                    build_fold, cd4_frozen_and_partial, compose_arm, compose_pp, de_profiles, deflate, design_and_realised,
                    energy, endpoints, frames, key_scale, pre_normrel, shared_energy_key, source_residual, union_rows,
                    weights)
import plolo  # noqa: E402

# ---- fold caches: clones of 420's (read-only use), so nothing is ever written into another experiment's dir
MYPLOLO = WORK / "plolo"
if not MYPLOLO.exists():
    MYPLOLO.mkdir(parents=True)
    for f in (W420 / "work/plolo").iterdir():
        subprocess.run(["cp", "-c", str(f), str(MYPLOLO / f.name)], check=True)
plolo.WORK = MYPLOLO
VIP = "viperturb_k562"
VIP_TABLE = os.environ.get("VCC424_VIP", "")     # "" = the primary full-axis table; else an axis variant in ./work/donors
VIP_PATH = (WORK / "donors" / f"{VIP}_{VIP_TABLE}") if VIP_TABLE else (ROOT / "data/donors" / VIP)
TAG = f"_{VIP_TABLE}" if VIP_TABLE else ""
plolo.DONOR_PATHS[VIP] = VIP_PATH       # for size410.donor_rows (the key-refit secondary only); in-process only
PANEL, OFFICIAL = plolo.PANEL, plolo.OFFICIAL
MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t")
Q = 1.25
ARMS = ("i:extra", "ii:avg", "ii-:avgonly", "ii':ivw", "iii:fill", "iii0:fillpure", "CTRL:shuf")


def vip_full_rows():
    """sens415.full_rows for VIPerturb (panel x official, panel columns kept) + availability; cached in ./work."""
    cache = WORK / f"full-{VIP}{TAG}.npz"
    if cache.exists():
        z = np.load(cache)
        return z["effects"], z["available"], json.loads(str(z["audit"]))
    bank = E062.Bank.load(VIP, VIP_PATH)
    mp = E062.map_bank_to_axis(bank, OFFICIAL)
    mean_all = E062.stream_column_mean(bank, mp.all_source_columns, chunk_rows=256, max_rss_gib=8.0)
    raw_median, _, _ = E062.centered_statistics(bank, mp.all_source_columns, mean_all, need_gene_scale=False,
                                                chunk_rows=256, max_rss_gib=8.0)
    at = {int(s): i for i, s in enumerate(mp.all_source_columns.tolist())}
    umean = np.asarray([mean_all[at[int(s)]] for s in mp.unique_source_columns])
    delta = np.load(bank.delta_path, mmap_mode="r")
    eff = np.zeros((len(PANEL), len(OFFICIAL)))
    avail = np.zeros(len(PANEL), bool)
    for i, t in enumerate(PANEL):
        r = bank.pert_pos.get(t)
        if r is None:
            continue
        eff[i, mp.unique_target_positions] = (np.asarray(delta[r, mp.unique_source_columns], float) - umean) \
            * (E062.EFFECT_SCALE / raw_median)
        avail[i] = True
    del delta
    audit = dict(raw_median=float(raw_median), mapped_genes=int(len(mp.unique_target_positions)),
                 duplicates=list(mp.duplicate_symbols), n_perts=int(bank.shape[0]), panel_available=int(avail.sum()))
    np.savez(cache, effects=eff, available=avail, audit=json.dumps(audit))
    return eff, avail, audit


def vip_cells():
    prov = json.loads((ROOT / "data/donors" / VIP / "provenance.json").read_text())   # cells: identical in every variant
    return prov["cells_per_target"]


def gwps_cells():
    """K562-GWPS cells per perturbation (its meta.json keeps only quantiles; the counts come from the h5ad obs)."""
    cache = WORK / "gwps_cells.json"
    if cache.exists():
        return json.loads(cache.read_text())
    import h5py
    path = ROOT / "data/external/scperturb/ReplogleWeissman2022_K562_gwps.h5ad"
    with h5py.File(path, "r") as h:
        node = h["obs/perturbation"]
        if isinstance(node, h5py.Group):
            cats = [c.decode() if isinstance(c, bytes) else str(c) for c in node["categories"][:]]
            codes = node["codes"][:]
            labels = np.asarray(cats, dtype=object)[codes]
        else:
            labels = np.asarray([c.decode() if isinstance(c, bytes) else str(c) for c in node[:]], dtype=object)
    u, c = np.unique(labels.astype(str), return_counts=True)
    out = {str(k): int(v) for k, v in zip(u, c)}
    cache.write_text(json.dumps(out))
    return out


def arm_view(dn, targets, sources, arm, rng=None):
    """(view, subset) for pre_normrel: a Donors-like namespace with the arm's F / A / med / pmcorr, and the source list.
    The base is (dn, sources) itself; VIPerturb is never used on a fold whose held line is K562 (the caller checks)."""
    assert "k562_gwps" not in getattr(dn, "held", ()), "leakage"
    effv, availv, _ = vip_full_rows()
    idx = {t: i for i, t in enumerate(PANEL)}
    ii = np.asarray([idx[t] for t in targets])
    Fv, Av = effv[ii].copy(), availv[ii].copy()
    F, A, med = dict(dn.F), dict(dn.A), dict(dn.med)
    pm = {k: v.copy() for k, v in dn.pmcorr.items()}
    medof = lambda M, a: float(np.median(np.linalg.norm(M[a], axis=1)))
    if arm in ("i:extra", "CTRL:shuf"):
        if arm == "CTRL:shuf":
            perm = rng.permutation(len(targets))
            Fv, Av = Fv[perm], Av[perm]
        F[VIP], A[VIP], med[VIP] = Fv, Av, medof(Fv, Av)
        return SimpleNamespace(F=F, A=A, med=med, pmcorr=pm), list(sources) + [VIP]
    assert "k562_gwps" in sources
    Fk, Ak = dn.F["k562_gwps"].copy(), dn.A["k562_gwps"].copy()
    both, only = Ak & Av, (~Ak) & Av
    if arm in ("ii:avg", "ii-:avgonly", "ii':ivw"):
        if arm == "ii':ivw":
            nv, ng = vip_cells(), gwps_cells()
            wv = np.asarray([nv.get(t, 0) / max(nv.get(t, 0) + ng.get(t, 0), 1) for t in targets])
        else:
            wv = np.full(len(targets), 0.5)
        wv = np.where(both, wv, 0.0)
        Fk[both] = (1 - wv[both, None]) * Fk[both] + wv[both, None] * Fv[both]
        for k in pm:
            pm[k][both] = (1 - wv[both, None]) * pm[k][both]
    if arm in ("ii:avg", "ii':ivw", "iii:fill", "iii0:fillpure"):
        Fk[only] = Fv[only]
        Ak = Ak | only
        for k in pm:
            assert not np.any(pm[k][only]), "PM correction on a target K562-GWPS lacks"
    elif arm != "ii-:avgonly":
        raise ValueError(arm)
    # iii0:fillpure keeps the K562 slot's normrel median at its K562-GWPS value, so ONLY the fill rows change
    F["k562_gwps"], A["k562_gwps"] = Fk, Ak
    if arm != "iii0:fillpure":
        med["k562_gwps"] = medof(Fk, Ak)
    return SimpleNamespace(F=F, A=A, med=med, pmcorr=pm), list(sources)


def arm_stats(dn, targets, sources, arm):
    """How many fold targets the arm touches (for the table)."""
    _, availv, _ = vip_full_rows()
    idx = {t: i for i, t in enumerate(PANEL)}
    Av = availv[np.asarray([idx[t] for t in targets])]
    Ak = dn.A["k562_gwps"] if "k562_gwps" in sources else np.zeros(len(targets), bool)
    return dict(n_targets=len(targets), vip=int(Av.sum()), both=int((Ak & Av).sum()), fill=int((~Ak & Av).sum()),
                gain_source=int(Av.sum()) if arm in ("i:extra", "CTRL:shuf") else int((~Ak & Av).sum()))


def score_rows(dls, scored, truth, betas, nonpanel, lams=(0.5, 1.0)):
    """comp420.score_rows (identical arithmetic)."""
    out, X, Aa, pdsp, pdsl, AB = {}, [], [], [], [], {}
    for ctx in "ABC":
        dl = dls[ctx]
        P = dl[:, scored]
        e = endpoints(P, truth)
        e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
        for lam in lams:
            e[f"pds_lam{lam:g}"] = endpoints((dl + lam * betas[ctx][None, :])[:, scored], truth)["pds"]
        out[ctx] = e
        X.append((P * truth).sum(1))
        Aa.append((P * P).sum(1))
        Tu = truth / np.maximum(np.linalg.norm(truth, axis=1, keepdims=True), 1e-300)
        for Pm, acc in ((P, pdsp), ((dl + 0.5 * betas[ctx][None, :])[:, scored], pdsl)):
            Pu = Pm / np.maximum(np.linalg.norm(Pm, axis=1, keepdims=True), 1e-300)
            S = Pu @ Tu.T
            own = np.diag(S).copy()
            n = S.shape[0]
            acc.append(((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1))
        H = dl[:, nonpanel]
        a_ = float((H * H).sum(1).mean())
        Hc = H - H.mean(0, keepdims=True)
        b_ = float((Hc * Hc).sum(1).mean())
        AB[ctx] = dict(A=a_, B=b_, BoverA=b_ / a_)
    keys = ("pds", "cos", "sse", "cos_agg", "sgn_w", "pds_lam0.5", "pds_lam1")
    summ = {k: float(np.mean([out[c][k] for c in "ABC"])) for k in keys}
    summ["per_frame"] = {c: {k: out[c][k] for k in ("pds", "cos_agg", "pds_lam0.5")} for c in "ABC"}
    summ["AB"] = AB
    arrays = dict(x=np.asarray(X), a=np.asarray(Aa), pds=np.asarray(pdsp), pdsl=np.asarray(pdsl))
    return summ, arrays


def cmp(s, r, fields=("cos_agg", "pds", "pds_lam0.5", "pds_lam1")):
    d = max(abs(s[k] - r[k]) for k in fields)
    return max(d, *[abs(s["AB"][c][k] / r["AB"][c][k] - 1) for c in "ABC" for k in ("A", "B")])
