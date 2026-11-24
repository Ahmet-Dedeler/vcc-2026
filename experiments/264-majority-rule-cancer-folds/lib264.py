"""264 -- the three-donor MAJORITY rule on the CANCER folds, with `kolf2` as the THIRD donor.  MEASUREMENT ONLY.

Everything is 263's construction, imported and not re-derived (`lib263`: 162's PLOLO fold caches, the per-donor
prepared rows rebuilt exactly as `plolo.build_fold` builds them, `combine` / `base_from` / `arm_parts` / `evaluate` =
261's `cd4svd5_parts` re-solved against the arm's own base, 188's exact `nmae` identity through `lib257.Identity`).
What is new here, and only here:

  * on each cancer fold `kolf2` is prepared AS A DONOR through the same call `build_fold` uses for every donor
    (`E062.prepare_main_source(bank, targets, OFFICIAL, panel_positions)`: per-gene centring over the bank's own
    perturbations, one global scale to a median centred row norm of EFFECT_SCALE = 6.817, panel columns zeroed) --
    the same call 263's `stage` makes for its `kolf2` arbiter rows -- and RESTRICTED to the
    fold's operator columns (the operator was fitted on the two production donors and is NOT re-fitted; a column the
    operator does not carry cannot be whitened, so the third donor contributes on the operator axis only; the energy
    outside it is measured and a sensitivity row keeps it);
  * the THREE-donor base (a) is `E062.build_baseline` over [production donor 1, production donor 2, kolf2] (float32-cast
    as the cache is); on the `kolf2` fold the three donors are the fold's own three and (a) is the cached `pre`
    bit for bit (263 KA-0), so that fold runs through the identical path as the known-answer fold;
  * `combine3` is 263's `combine` with three rows and the three-row mean as the untouched / norm-reference row;
  * the pickers: majority (drop the donor whose un-whitened row on the operator columns has the lowest mean cosine
    to the other two), reverse (drop the highest), random drop (3 seeds), always-drop-<donor>, oracle drop (the donor
    whose row has the lowest cosine to the held truth on the scored columns -- a ceiling reference, never an arm).

Read-only everywhere except `experiments/264-majority-rule-cancer-folds/work/`; the fold caches under `work/plolo/`
are byte copies of 263's (= 262's = 257's = 261's), sha-listed in `work/plolo/SHA256SUMS.txt`.
"""
from __future__ import annotations

import hashlib
import resource
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E263 = ROOT / "experiments/263-donor-selection-arbiters"
sys.path.insert(0, str(E263))
import lib263 as Y  # noqa: E402  (imports lib257 as Y.X, lib246 as Y.L, plolo as Y.plolo)

X, L, plolo = Y.X, Y.L, Y.plolo
# lib263 pins plolo.WORK / L.WORK / X.WORK to 263's directory at import; pin every one of them to OURS so each cache
# read comes from this experiment's byte copies and nothing outside it is ever touched.
plolo.WORK = WORK / "plolo"
Y.WORK = WORK
L.WORK = WORK
X.WORK = WORK
E062 = plolo.E062
OFFICIAL = plolo.OFFICIAL
ALPHA, GAMMA, LAM = Y.ALPHA, Y.GAMMA, Y.LAM
CANCER = Y.CANCER
FOLDS = Y.FOLDS
SHUF_SEED = Y.SHUF_SEED
PAIRS = ((0, 1), (0, 2), (1, 2))


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_array(a: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()


def row_cos(A, B):
    return Y.row_cos(A, B)


# --------------------------------------------------------------------------------------------- staging
def stage_base(held: str) -> tuple[dict, list, object]:
    """263's `stage`, verbatim, but keeping the per-donor PreparedSource list and kolf2's PreparedSource (one preparation
    pass instead of two).  The operator STAYS RESIDENT so every arm's base is whitened through `Operator.apply`."""
    fold = plolo.build_fold(held)
    spec = plolo.FOLDS[held]
    targets = list(fold["targets"])
    n = len(targets)
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    sc = np.asarray(fold["scored"], dtype=np.int64)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
    pre_cache = np.asarray(fold["pre"], dtype=np.float64)
    sources, banks, pm_n = Y.prepared_sources(held, fold)
    donors = [s.name for s in sources]
    rows = [np.asarray(s.effects, dtype=np.float64) for s in sources]
    avail = np.stack([np.asarray(s.target_available, bool) for s in sources])
    pre64, votes = E062.build_baseline(sources)
    ka0_pre = float(np.abs(pre64.astype(np.float32).astype(np.float64) - pre_cache).max())
    ka0_votes = int(np.abs(votes - avail.sum(0)).max())
    cd4 = plolo.cd4_for(targets)
    union = L.union_rows([L.source_residual(nm, targets, "gene_norm") for nm in ("gse295214_h1", "h1_2025")])
    kp = None
    kolf_rows = kolf_avail = None
    if held != "kolf2":
        kb = E062.Bank.load("kolf2", plolo.DONOR_PATHS["kolf2"])
        kp = E062.prepare_main_source(kb, targets, OFFICIAL, fold["panel_positions"], chunk_rows=256, max_rss_gib=40.0)
        kolf_rows = np.asarray(kp.effects, dtype=np.float64)
        kolf_avail = np.asarray(kp.target_available, bool)
    b_held = np.asarray(fold["b_ctrl"], dtype=np.float64)
    b_donor = [plolo._b_ctrl_on_axis(banks[nm])[0] for nm in spec["sources"]]
    op = plolo.Operator(L.op_path_for(held))
    whitened = op.apply(pre_cache * ALPHA, wcols, GAMMA, LAM)
    wht_rows = [op.apply(r * ALPHA, wcols, GAMMA, LAM) for r in rows]
    frac_cd4, _ = plolo.solve_fraction(whitened, cd4, Y.CD4_F)
    frac_x, _ = Y.solve_fraction_carrying(whitened, union, Y.H1_G)
    d5, aud5 = L.deflate(cd4, "svd", 5)
    st = dict(held=held, fold=fold, targets=targets, n=n, wcols=wcols, sc=sc, truth=truth, pre=pre_cache,
              whitened=whitened, cd4=cd4, union=union, frac_cd4=float(frac_cd4), frac_x=float(frac_x),
              d5=d5, aud5=aud5, donors=donors, sources=list(spec["sources"]), rows=rows, avail=avail,
              wht_rows=wht_rows, votes=votes, kolf_rows=kolf_rows, kolf_avail=kolf_avail,
              b_held=b_held, b_donor=b_donor, op=op, pm_n=pm_n,
              ka0=dict(pre_float32_max_abs_diff=ka0_pre, votes_max_abs_diff=ka0_votes))
    return st, sources, kp


def whitening_audit(st: dict, rows: np.ndarray, avail: np.ndarray) -> dict:
    """How the operator treats these rows: un-preserved whitened median norm over the raw median (the amplification),
    the energy on frame A's CPM gate before / after whitening, the energy in the operator's bottom eigen-decile, and
    the delivered / effect norm ratio of the rows alone (frame A, no mass fix)."""
    op, wc = st["op"], st["wcols"]
    bA = plolo.challenge_bctrl("A")
    eA = np.expm1(bA)
    gA = (1e6 * eA / float(eA.sum()) > plolo.MIN_CPM) & (bA > 0)
    R = rows[avail][:, wc]
    nr = np.linalg.norm(R, axis=1)
    w = (op.evals + LAM) ** (-GAMMA / 2.0)
    order = np.argsort(op.evals)[::-1]
    T = (R @ op.evecs) * w
    nt = np.linalg.norm(T, axis=1)
    Tw = T @ op.evecs.T
    e_T = (T ** 2).sum(0)
    bottom = order[9 * len(order) // 10:]
    full = np.zeros((R.shape[0], rows.shape[1]))
    full[:, wc] = Tw * (np.median(nr) / np.median(nt))
    H = plolo.deliver(full, bA, gA, None)
    return dict(n=int(avail.sum()), raw_median=float(np.median(nr)), wht_unpreserved_median=float(np.median(nt)),
                amplification=float(np.median(nt) / np.median(nr)),
                energy_on_gateA_raw=float((R[:, gA[wc]] ** 2).sum() / (R * R).sum()),
                energy_on_gateA_wht=float((Tw[:, gA[wc]] ** 2).sum() / (Tw * Tw).sum()),
                energy_bottom_eigen_decile_wht=float(e_T[bottom].sum() / e_T.sum()),
                delivered_over_effect=float(np.median(np.linalg.norm(H[:, wc], axis=1)) / np.median(nr)))


def kolf2_as_donor(st: dict, kp, *, weight: float = 1.0) -> tuple[E062.PreparedSource, dict]:
    """`kolf2` prepared as a donor for THIS fold's targets (the `build_fold` call, already made in `stage_base`),
    restricted to the fold's operator columns and scaled by `weight` (1.0 = the mandate's equal mean; the
    `matched` mode passes the fold-constant scalar that equalises its whitened amplitude with the production donors').
    Returns the PreparedSource and an audit (R3 material)."""
    eff = np.asarray(kp.effects, dtype=np.float64)
    avail = np.asarray(kp.target_available, bool)
    # R3: the same rows are what 263's kolf2 fold consumed as its TRUTH (`E062.raw_truth` = effects * raw_median /
    # EFFECT_SCALE, float32 in the cache), target by target
    zk = np.load(plolo.WORK / "fold-kolf2.npz", allow_pickle=True)
    truth_k = np.asarray(zk["truth"])                                  # float32, as cached
    pos_k = {str(t): i for i, t in enumerate(zk["targets"])}
    common = [i for i, t in enumerate(st["targets"]) if avail[i] and t in pos_k]
    mine = (eff[common] * (kp.raw_median_norm / E062.EFFECT_SCALE)).astype(np.float32)
    theirs = truth_k[[pos_k[st["targets"][i]] for i in common]]
    r3 = dict(n_common=len(common), n_kolf_carried=int(avail.sum()),
              max_abs_diff_float32=float(np.abs(mine.astype(np.float64) - theirs.astype(np.float64)).max()) if common else float("nan"),
              rows_bit_identical=int(sum(np.array_equal(mine[k], theirs[k]) for k in range(len(common)))),
              sha256_mine=sha256_array(mine), sha256_kolf2_fold_truth=sha256_array(theirs),
              raw_median_norm=float(kp.raw_median_norm), effect_scale=float(E062.EFFECT_SCALE),
              sha256_effects_float64=sha256_array(eff),
              kolf2_targets_not_in_fold=[str(t) for t in zk["targets"] if str(t) not in set(st["targets"])],
              fold_targets_kolf2_lacks=[st["targets"][i] for i in range(len(st["targets"])) if not avail[i]])
    wcols = np.asarray(st["wcols"], dtype=np.int64)
    inside = np.zeros(eff.shape[1], dtype=bool)
    inside[wcols] = True
    carried_cols = np.flatnonzero(np.abs(eff).sum(0) > 0)
    outside_cols = carried_cols[~inside[carried_cols]]
    e_all = float((eff * eff).sum())
    e_out = float((eff[:, outside_cols] ** 2).sum()) if outside_cols.size else 0.0
    eff_r = eff.copy()
    eff_r[:, ~inside] = 0.0
    eff_r *= weight
    moa = np.asarray(kp.mean_on_axis, dtype=np.float64).copy()
    moa[~inside] = 0.0
    moa *= weight
    restrict = dict(n_carried_cols=int(carried_cols.size), n_cols_outside_operator=int(outside_cols.size),
                    energy_frac_outside_operator=(e_out / e_all) if e_all > 0 else 0.0, n_operator_cols=int(wcols.size),
                    kolf2_weight=float(weight), sha256_effects_restricted_weighted=sha256_array(eff_r))
    src = E062.PreparedSource(name="kolf2", effects=eff_r, target_available=avail.copy(), mean_on_axis=moa,
                              raw_median_norm=float(kp.raw_median_norm), mapped_gene_count=int(kp.mapped_gene_count),
                              duplicate_symbols=tuple(kp.duplicate_symbols))
    src_unrestricted = E062.PreparedSource(name="kolf2_unrestricted", effects=eff * weight, target_available=avail.copy(),
                                           mean_on_axis=np.asarray(kp.mean_on_axis, dtype=np.float64) * weight,
                                           raw_median_norm=float(kp.raw_median_norm), mapped_gene_count=int(kp.mapped_gene_count),
                                           duplicate_symbols=tuple(kp.duplicate_symbols))
    return src, dict(r3=r3, restrict=restrict, unrestricted=src_unrestricted)


def stage3(held: str, kolf_weight: str = "equal") -> dict:
    """263's staging (two production donors on a cancer fold; three on `kolf2`), plus the third donor and the
    three-donor base (a).  `st['pre']` / `st['whitened']` stay the INCUMBENT's (c = cd4svd5's base).
    kolf_weight: 'equal' = the mandate's construction (weight 1.0); 'matched' = a FOLD-CONSTANT scalar on kolf2's rows
    that equalises its un-preserved whitened median row norm with the mean of the two production donors' (the
    sensitivity construction; the raw-space selection is scale-invariant, so the picks are identical)."""
    st, sources2, kp = stage_base(held)
    wa = {nm: whitening_audit(st, st["rows"][j], st["avail"][j]) for j, nm in enumerate(st["donors"])}
    if held == "kolf2":
        sources3 = list(sources2)
        aud3 = dict(third_donor=None, kolf2_weight=1.0, note="the fold's own three production donors; kolf2 is the truth")
        unres = None
    else:
        wk = whitening_audit(st, st["kolf_rows"], st["kolf_avail"])
        wa["kolf2 (weight 1)"] = wk
        if kolf_weight == "matched":
            weight = float(np.mean([wa[nm]["wht_unpreserved_median"] for nm in st["donors"]]) / wk["wht_unpreserved_median"])
        elif kolf_weight == "equal":
            weight = 1.0
        else:
            raise ValueError(kolf_weight)
        ksrc, aud = kolf2_as_donor(st, kp, weight=weight)
        sources3 = list(sources2) + [ksrc]
        aud3 = dict(third_donor="kolf2", kolf_weight_mode=kolf_weight, **aud["r3"], **aud["restrict"])
        unres = aud["unrestricted"]
    pre3_64, votes3 = E062.build_baseline(sources3)
    pre_a = pre3_64.astype(np.float32).astype(np.float64)
    rows3 = [np.asarray(s.effects, dtype=np.float64) for s in sources3]
    avail3 = np.stack([np.asarray(s.target_available, bool) for s in sources3])
    donors3 = [s.name for s in sources3]
    op = st["op"]
    wht_rows3 = list(st["wht_rows"]) + ([op.apply(rows3[-1] * ALPHA, st["wcols"], GAMMA, LAM)] if len(rows3) > len(st["wht_rows"]) else [])
    whitened_a = op.apply(pre_a * ALPHA, st["wcols"], GAMMA, LAM)
    three = avail3.all(0)
    wa["mean2 (cached pre), three-donor targets"] = whitening_audit(st, st["pre"], three)
    wa["mean3 (a), three-donor targets"] = whitening_audit(st, pre_a, three)
    if held != "kolf2":
        wa["kolf2 (as used)"] = whitening_audit(st, rows3[-1], avail3[-1])
    changed_vs_cache = np.abs(pre_a - st["pre"]).max(1) > 0
    filler = votes3 == 0
    st3 = dict(st, sources3=sources3, rows3=rows3, avail3=avail3, donors3=donors3, votes3=votes3, pre_a=pre_a,
               whitened_a=whitened_a, wht_rows3=wht_rows3, kolf_unrestricted=unres, whitening_audit=wa,
               aud3=dict(**aud3, n_three_donor_targets=int(three.sum()), n_changed_vs_cached_pre=int(changed_vs_cache.sum()),
                         n_filler_rows=int(filler.sum()), votes3_hist={int(k): int(v) for k, v in zip(*np.unique(votes3, return_counts=True))},
                         pre_a_equals_cached_pre=bool(np.array_equal(pre_a, st["pre"]))))
    return st3


# ---------------------------------------------------------------------------------------- combining
def combine3(st3: dict, W, *, norm_match: bool = False, norm_ref=None) -> np.ndarray:
    """263's `combine` on the THREE donor rows: the three-donor mean row (a) where W's row is NaN, else
    sum_j W[t, j] * row_j[t] over the donors carrying t (renormalised), float32-cast as the cache was; `norm_match`
    rescales the combined row to the reference row's norm (default: the three-row mean, i.e. (a); `norm_ref` = the
    cached incumbent `pre` for the `:nmc` cells)."""
    pre = st3["pre_a"].copy()
    rows, avail = st3["rows3"], st3["avail3"]
    ref = st3["pre_a"] if norm_ref is None else norm_ref
    use = np.isfinite(W).all(1) & (avail.sum(0) >= 1)
    idx = np.flatnonzero(use)
    if idx.size:
        Wu = W[idx].copy()
        Wu[~avail[:, idx].T] = 0.0
        s = Wu.sum(1)
        assert np.all(s > 0), "a selected target has no available donor with positive weight"
        Wu /= s[:, None]
        acc = np.zeros((idx.size, pre.shape[1]))
        for j in range(len(rows)):
            acc += Wu[:, j][:, None] * rows[j][idx]
        if norm_match:
            na = np.linalg.norm(acc, axis=1)
            nm = np.linalg.norm(ref[idx], axis=1)
            ok = (na > 0) & (nm > 0)
            acc[ok] *= (nm[ok] / na[ok])[:, None]
        pre[idx] = acc.astype(np.float32).astype(np.float64)
    return pre


def drop_weights(n: int, J: int, drop: np.ndarray) -> np.ndarray:
    """W with the dropped donor at 0 and the others at 1 on targets where drop >= 0; NaN (keep (a)) elsewhere."""
    W = np.full((n, J), np.nan)
    for t in np.flatnonzero(drop >= 0):
        W[t] = 1.0
        W[t, drop[t]] = 0.0
    return W


# ------------------------------------------------------------------------------------------ pickers
def pair_cosines(rows_on_cols, avail3):
    """(n, 3) cosine of the donor rows for pairs (0,1), (0,2), (1,2); NaN where either donor lacks the target."""
    n = avail3.shape[1]
    out = np.full((n, len(PAIRS)), np.nan)
    for k, (a, b) in enumerate(PAIRS):
        c = row_cos(rows_on_cols[a], rows_on_cols[b])
        c[~(avail3[a] & avail3[b])] = np.nan
        out[:, k] = c
    return out


def majority_drop(pairs_cos: np.ndarray, avail3: np.ndarray, gate: np.ndarray, *, reverse: bool = False) -> np.ndarray:
    """263's rule verbatim: per target carried by all three (and gated), each donor's score = the mean of its two
    pairwise cosines; drop the donor with the LOWEST score (reverse: the highest).  Returns the dropped index or -1."""
    n, J = avail3.shape[1], avail3.shape[0]
    drop = np.full(n, -1, dtype=np.int64)
    for t in np.flatnonzero(gate):
        if not avail3[:, t].all() or not np.isfinite(pairs_cos[t]).all():
            continue
        c = {PAIRS[k]: pairs_cos[t, k] for k in range(len(PAIRS))}
        a = np.array([np.mean([c[pr] for pr in PAIRS if j in pr]) for j in range(J)])
        drop[t] = int(np.argmax(a) if reverse else np.argmin(a))
    return drop


def random_drop(avail3: np.ndarray, gate: np.ndarray, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = avail3.shape[1]
    drop = np.full(n, -1, dtype=np.int64)
    for t in np.flatnonzero(gate):
        if avail3[:, t].all():
            drop[t] = int(rng.integers(0, avail3.shape[0]))
    return drop


def fixed_drop(avail3: np.ndarray, gate: np.ndarray, j: int) -> np.ndarray:
    n = avail3.shape[1]
    drop = np.full(n, -1, dtype=np.int64)
    m = gate & avail3.all(0)
    drop[m] = j
    return drop


def oracle_drop(cos_truth: np.ndarray, avail3: np.ndarray, gate: np.ndarray) -> np.ndarray:
    """Drop the donor whose row has the LOWEST cosine to the held truth (a ceiling reference for the same operation)."""
    n = avail3.shape[1]
    drop = np.full(n, -1, dtype=np.int64)
    for t in np.flatnonzero(gate):
        if avail3[:, t].all() and np.isfinite(cos_truth[t]).all():
            drop[t] = int(np.argmin(cos_truth[t]))
    return drop
