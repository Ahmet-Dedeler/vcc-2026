"""263 -- per-target donor selection on the whitened base, with BUILD-TIME-OBSERVABLE arbiters.  MEASUREMENT ONLY.

Construction (nothing mathematical is new downstream of the base):
  * the fold's per-donor prepared rows are rebuilt exactly as `plolo.build_fold` builds them (`E062.prepare_main_source`
    per donor, the PerturbMatch 0.5 blend on `k562_gwps`'s shared columns, `E062.build_baseline`); the rebuilt mean cast
    to float32 must equal the cached `pre` BIT FOR BIT (KA-0), so the per-donor rows are the ones the incumbent's base
    was averaged from;
  * a per-target selection / weighting replaces the mean row by the same combination of the donors' un-whitened rows
    (targets the arbiter does not touch keep the cached mean row verbatim), the selected matrix is whitened by the
    fold's operator through `plolo.Operator.apply` itself (the production path: linear whitening on the operator
    columns, ONE global preserve), so the trivial selection is `cd4svd5` at 0.0e+00 (R1);
  * everything downstream is 261's `cd4svd5_parts`: CD4 deflated k = 5 (`lib246.deflate(cd4, "svd", 5)`), RE-SOLVED at
    f = 0.34 by `plolo.solve_fraction` against the arm's own base; the H1-pair union re-solved by 218's carrying-row rule
    at g = 0.28 against the same base; one add, one preserve (`lib246.compose`); scored by `lib246.score` (a23 + a15
    tags, the scorer's per-target rank statistic); 188's exact `nmae` identity through `lib257.Identity`.

Read-only everywhere except `experiments/263-donor-selection-arbiters/work/`; the fold caches under `work/plolo/` are
byte copies of 262's (= 257's = 261's), sha-listed in `work/plolo/SHA256SUMS.txt`.
"""
from __future__ import annotations

import resource
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E257 = ROOT / "experiments/257-soft-deflation-price"
E262 = ROOT / "experiments/262-where-pds-is-lost"
sys.path.insert(0, str(E257))
import lib257 as X  # noqa: E402  (imports lib246 as X.L -> plolo / resid218 / compose218)
import plolo  # noqa: E402
from compose218 import solve_fraction_carrying  # noqa: E402

L = X.L
# lib246 pins plolo.WORK to 246's directory and lib257 to 257's; pin it to OURS so every cache read comes from this
# experiment's byte copies and nothing outside it is ever touched.
plolo.WORK = WORK / "plolo"
L.WORK = WORK
X.WORK = WORK
E062, E064 = plolo.E062, plolo.E064
OFFICIAL, PANEL_SET = plolo.OFFICIAL, plolo.PANEL_SET
ALPHA, GAMMA, LAM = L.ALPHA, L.GAMMA, L.LAM
CD4_F, H1_G = L.CD4_F, L.H1_G
CANCER = L.CANCER
FOLDS = ("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2")
SHUF_SEED = L.SHUF_SEED      # 20726
RAND_SEED = L.RAND_SEED      # 20260904
KNN = 5


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


# ------------------------------------------------------------------------------------ per-donor rows
def prepared_sources(held: str, fold: dict):
    """`plolo.build_fold`'s per-source PreparedSource list, verbatim (including the PerturbMatch blend on the
    k562_gwps rows), so `E062.build_baseline` on it is the cached `pre`."""
    spec = plolo.FOLDS[held]
    banks = {n: E062.Bank.load(n, plolo.DONOR_PATHS[n]) for n in spec["sources"]}
    targets, panel_positions = fold["targets"], fold["panel_positions"]
    prepared = {n: E062.prepare_main_source(b, targets, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=40.0)
                for n, b in banks.items()}
    pm_n = 0
    if "k562_gwps" in banks:
        pm_bank = E064.load_perturbmatch(plolo.PM_PATH)
        pm_prepared = E062.prepare_main_source(pm_bank, targets, OFFICIAL, panel_positions, chunk_rows=256, max_rss_gib=40.0)
        shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET)
    source_prepared = []
    for n in spec["sources"]:
        if n == "k562_gwps":
            eff = np.asarray(prepared[n].effects, dtype=np.float64).copy()
            both = prepared[n].target_available & pm_prepared.target_available
            rows = np.flatnonzero(both)
            pm_n = int(both.sum())
            if len(rows):
                eff[np.ix_(rows, shared)] = (
                    (1 - plolo.PM_BLEND) * prepared[n].effects[np.ix_(rows, shared)]
                    + plolo.PM_BLEND * pm_prepared.effects[np.ix_(rows, shared)])
            source_prepared.append(E062.PreparedSource(
                name="k562_gwps_pm50", effects=eff,
                target_available=np.asarray(prepared[n].target_available, bool).copy(),
                mean_on_axis=np.asarray(prepared[n].mean_on_axis, float).copy(),
                raw_median_norm=float(prepared[n].raw_median_norm),
                mapped_gene_count=int(prepared[n].mapped_gene_count),
                duplicate_symbols=tuple(prepared[n].duplicate_symbols)))
        else:
            source_prepared.append(prepared[n])
    return source_prepared, banks, pm_n


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def row_cos(A, B):
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(1) / (na * nb)
    return np.where((na > 0) & (nb > 0), c, np.nan)


# --------------------------------------------------------------------------------------------- staging
def stage(held: str) -> dict:
    """Everything an arm needs.  The operator STAYS RESIDENT (2.6 GB) so every arm's base is whitened through
    `Operator.apply` itself; nothing is approximated by linearity."""
    fold = plolo.build_fold(held)
    spec = plolo.FOLDS[held]
    targets = list(fold["targets"])
    n = len(targets)
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    sc = np.asarray(fold["scored"], dtype=np.int64)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
    pre_cache = np.asarray(fold["pre"], dtype=np.float64)

    sources, banks, pm_n = prepared_sources(held, fold)
    donors = [s.name for s in sources]
    rows = [np.asarray(s.effects, dtype=np.float64) for s in sources]
    avail = np.stack([np.asarray(s.target_available, bool) for s in sources])        # (J, n)
    pre64, votes = E062.build_baseline(sources)
    ka0_pre = float(np.abs(pre64.astype(np.float32).astype(np.float64) - pre_cache).max())
    ka0_votes = int(np.abs(votes - avail.sum(0)).max())

    cd4 = plolo.cd4_for(targets)
    union = L.union_rows([L.source_residual(nm, targets, "gene_norm") for nm in ("gse295214_h1", "h1_2025")])

    # kolf2's rows as an EXTERNAL arbiter on the cancer folds (on disk, not a production donor; never the truth there)
    kolf_rows = kolf_avail = None
    if held != "kolf2":
        kb = E062.Bank.load("kolf2", plolo.DONOR_PATHS["kolf2"])
        kp = E062.prepare_main_source(kb, targets, OFFICIAL, fold["panel_positions"], chunk_rows=256, max_rss_gib=40.0)
        kolf_rows = np.asarray(kp.effects, dtype=np.float64)
        kolf_avail = np.asarray(kp.target_available, bool)

    # basal profiles (log-pseudobulk on the official axis, renormalised as plolo does)
    b_held = np.asarray(fold["b_ctrl"], dtype=np.float64)
    b_donor = [plolo._b_ctrl_on_axis(banks[nm])[0] for nm in spec["sources"]]

    op = plolo.Operator(L.op_path_for(held))
    whitened = op.apply(pre_cache * ALPHA, wcols, GAMMA, LAM)              # == lib246.stage's base, bit for bit
    wht_rows = [op.apply(r * ALPHA, wcols, GAMMA, LAM) for r in rows]      # each donor's rows, whitened alone
    frac_cd4, _ = plolo.solve_fraction(whitened, cd4, CD4_F)
    frac_x, _ = solve_fraction_carrying(whitened, union, H1_G)
    d5, aud5 = L.deflate(cd4, "svd", 5)
    st = dict(held=held, fold=fold, targets=targets, n=n, wcols=wcols, sc=sc, truth=truth, pre=pre_cache,
              whitened=whitened, cd4=cd4, union=union, frac_cd4=float(frac_cd4), frac_x=float(frac_x),
              d5=d5, aud5=aud5, donors=donors, sources=list(spec["sources"]), rows=rows, avail=avail,
              wht_rows=wht_rows, votes=votes, kolf_rows=kolf_rows, kolf_avail=kolf_avail,
              b_held=b_held, b_donor=b_donor, op=op, pm_n=pm_n,
              ka0=dict(pre_float32_max_abs_diff=ka0_pre, votes_max_abs_diff=ka0_votes))
    return st


# ---------------------------------------------------------------------------------------- the arms
def combine(st: dict, W, *, explicit_mean: bool = False, norm_match: bool = False) -> np.ndarray:
    """The selected un-whitened base: the cached mean row where W's row is NaN, else sum_j W[t, j] * row_j[t]
    (float32-cast, as the cache was).  `norm_match=True` rescales every combined row to the cached mean row's norm
    (the DIRECTION-only version of the selection: the global preserve and the residual sizing then see the same row
    norms as the incumbent, so the arm is shape-only).  `explicit_mean=True` ignores W and recomputes EVERY covered
    row as `build_baseline` does (sum over the carrying donors / count) -- the R1 known-answer path."""
    pre = st["pre"].copy()
    rows, avail = st["rows"], st["avail"]
    J, n = avail.shape
    if explicit_mean:
        acc = np.zeros_like(pre)
        cnt = np.zeros(n, dtype=np.int64)
        for j in range(J):
            m = avail[j]
            acc[m] += rows[j][m]
            cnt += m.astype(np.int64)
        cov = cnt > 0
        acc[cov] /= cnt[cov, None]
        pre[cov] = acc[cov].astype(np.float32).astype(np.float64)
        return pre
    use = np.isfinite(W).all(1) & (avail.sum(0) >= 1)
    idx = np.flatnonzero(use)
    if idx.size:
        Wu = W[idx].copy()
        Wu[~avail[:, idx].T] = 0.0
        s = Wu.sum(1)
        assert np.all(s > 0), "a selected target has no available donor with positive weight"
        Wu /= s[:, None]
        acc = np.zeros((idx.size, pre.shape[1]))
        for j in range(J):
            acc += Wu[:, j][:, None] * rows[j][idx]
        if norm_match:
            na = np.linalg.norm(acc, axis=1)
            nm = np.linalg.norm(pre[idx], axis=1)
            ok = (na > 0) & (nm > 0)
            acc[ok] *= (nm[ok] / na[ok])[:, None]
        pre[idx] = acc.astype(np.float32).astype(np.float64)
    return pre


def base_from(st: dict, pre_sel: np.ndarray) -> np.ndarray:
    return st["op"].apply(pre_sel * ALPHA, st["wcols"], GAMMA, LAM)


def arm_parts(st: dict, whitened: np.ndarray):
    """261's `cd4svd5_parts` against THIS base: the deflated CD4 re-solved at f = 0.34 (all-rows rule), the union
    re-solved by the carrying-row rule at g = 0.28."""
    fr, pu = plolo.solve_fraction(whitened, st["d5"], CD4_F)
    fx, px = solve_fraction_carrying(whitened, st["union"], H1_G)
    return [(st["d5"], float(fr)), (st["union"], float(fx))], dict(cd4_fraction=float(fr), cd4_per_unit=float(pu),
                                                                  union_fraction=float(fx), union_per_unit=float(px))


def evaluate(st: dict, whitened: np.ndarray, ident) -> tuple[dict, dict]:
    parts, aud = arm_parts(st, whitened)
    st_arm = dict(st, whitened=whitened)
    final = L.compose(st_arm, parts)
    e = L.score(st_arm, final)
    m = ident.measure(final)
    aud["base_median_norm_wcols"] = float(np.median(np.linalg.norm(whitened[:, st["wcols"]], axis=1)))
    return e, dict(**m, audit=aud)


# ----------------------------------------------------------------------------------- selections
def one_hot(st: dict, pick: np.ndarray) -> np.ndarray:
    """pick[t] = donor index, or -1 = keep the mean row."""
    n, J = st["n"], len(st["rows"])
    W = np.full((n, J), np.nan)
    for t in np.flatnonzero(pick >= 0):
        W[t] = 0.0
        W[t, pick[t]] = 1.0
    return W


def softmax_weights(st: dict, S: np.ndarray, beta: float, gate: np.ndarray) -> np.ndarray:
    """w_j ∝ exp(beta * z_j) over the available donors, z = the arbiter score centred per target and divided by the
    SD of the centred scores over the targets with every donor's score defined (one scale per fold)."""
    n, J = S.shape
    ok = np.isfinite(S).all(1) & st["avail"].all(0)
    Sc = S - np.nanmean(S, axis=1, keepdims=True)
    scale = float(np.std(Sc[ok])) if ok.sum() > 2 else 1.0
    scale = scale if scale > 0 else 1.0
    W = np.full((n, J), np.nan)
    for t in np.flatnonzero(gate & ok):
        z = Sc[t] / scale
        w = np.exp(beta * (z - z.max()))
        W[t] = w / w.sum()
    return W


def pick_from_scores(S: np.ndarray, avail: np.ndarray, gate: np.ndarray, *, reverse: bool = False) -> np.ndarray:
    """argmax (or argmin) of the score over the available donors, on the gated targets; -1 elsewhere or where fewer
    than two donors have a defined score."""
    n, J = S.shape
    pick = np.full(n, -1, dtype=np.int64)
    for t in np.flatnonzero(gate):
        s = np.where(avail[:, t] & np.isfinite(S[t]), S[t], np.nan)
        if np.isfinite(s).sum() < 2:
            continue
        pick[t] = int(np.nanargmin(s) if reverse else np.nanargmax(s))
    return pick


def random_pick(avail: np.ndarray, gate: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    J, n = avail.shape
    pick = np.full(n, -1, dtype=np.int64)
    for t in np.flatnonzero(gate):
        cand = np.flatnonzero(avail[:, t])
        if cand.size >= 2:
            pick[t] = int(rng.choice(cand))
    return pick


def knn_self_consistency(rows_j, avail_j, other, other_avail, k: int = KNN) -> np.ndarray:
    """For each target t: the k nearest targets to t in the OTHER space (cosine, among targets the other space
    carries), then the mean cosine of donor j's row for t to donor j's rows for those neighbours (those j carries)."""
    n = rows_j.shape[0]
    Uo = unit(other)
    Go = Uo @ Uo.T
    Uj = unit(rows_j)
    Gj = Uj @ Uj.T
    out = np.full(n, np.nan)
    for t in range(n):
        if not (avail_j[t] and other_avail[t]):
            continue
        g = Go[t].copy()
        g[t] = -np.inf
        g[~other_avail] = -np.inf
        nb = np.argsort(-g)[:k]
        nb = nb[np.isfinite(g[nb]) & avail_j[nb]]
        if nb.size:
            out[t] = float(Gj[t, nb].mean())
    return out
