"""299 -- fit the fold's whitening operator on a CHOSEN SET OF TRAINING ROWS and read the arm through `cd4svd5`'s
construction.  ONE FOLD PER PROCESS.  Every arm's row is written as soon as it is measured, so a killed process
loses only the arm in flight.

    PYTHONDONTWRITEBYTECODE=1 uv run python run299.py orion_hct116 cd4svd5 STACK full sel:q0.5 sel:q0.25 rand:q0.5 mag:q0.5 rev:q0.5

Arms
  cd4svd5      the incumbent through the CACHED fold operator (KA-3 / KA-4 against 263)
  STACK        h1union (the un-deflated CD4 at f = 0.34 + the pair union at g = 0.28) -- the BR249 retrodiction leg
  full         the operator REFIT on every row of both donor banks (KA-2: must reproduce the cached operator)
  sel:q<Q>     refit on the TOP-Q fraction of each donor's non-panel bank rows by split-half reliability + every
               panel-target row the bank carries
  rev:q<Q>     the BOTTOM-Q by reliability          mag:q<Q>  the TOP-Q by magnitude
  rand:q<Q>[:s<seed>]   a uniform random subset of the same size

Everything downstream of the base is `cd4svd5` verbatim (lib263 / lib246 / lib257).  ONE VARIABLE: the row set.
"""
from __future__ import annotations

import gc
import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/299-reliable-training-rows"
WORK = HERE / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
E272 = ROOT / "experiments/272-noise-corrected-operator"
sys.path.insert(0, str(E263))
sys.path.insert(0, str(E272))
sys.path.insert(0, str(E271))
import lib263 as Y  # noqa: E402
import op272 as O  # noqa: E402

Y.plolo.WORK = WORK / "plolo"
Y.L.WORK = WORK
Y.X.WORK = WORK
Y.WORK = WORK
plolo, L, X = Y.plolo, Y.L, Y.X
E062 = plolo.E062
OFFICIAL, PANEL, PANEL_SET = plolo.OFFICIAL, list(plolo.PANEL), plolo.PANEL_SET
ALPHA, GAMMA, LAM = Y.ALPHA, Y.GAMMA, Y.LAM
KEYS = ("cos_agg", "a23_cos", "pds")
TOL = 1e-9
CHUNK = 2048
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_self() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------ staging (lib263.stage, minus
# the kolf2 arbiter rows and the per-donor whitened rows, which no arm here uses)
def stage299(held: str) -> dict:
    fold = plolo.build_fold(held)
    spec = plolo.FOLDS[held]
    targets = list(fold["targets"])
    n = len(targets)
    wcols = np.asarray(fold["wcols"], dtype=np.int64)
    sc = np.asarray(fold["scored"], dtype=np.int64)
    truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
    pre_cache = np.asarray(fold["pre"], dtype=np.float64)

    banks = {nm: E062.Bank.load(nm, plolo.DONOR_PATHS[nm]) for nm in spec["sources"]}
    prepared = {nm: E062.prepare_main_source(b, targets, OFFICIAL, fold["panel_positions"], chunk_rows=256, max_rss_gib=40.0)
                for nm, b in banks.items()}
    raw_full = {nm: float(p.raw_median_norm) for nm, p in prepared.items()}
    pm_n = 0
    if "k562_gwps" in banks:
        pm_bank = plolo.E064.load_perturbmatch(plolo.PM_PATH)
        pm_prepared = E062.prepare_main_source(pm_bank, targets, OFFICIAL, fold["panel_positions"], chunk_rows=256, max_rss_gib=40.0)
        shared = plolo.E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, OFFICIAL, PANEL_SET)
    sources = []
    for nm in spec["sources"]:
        if nm == "k562_gwps":
            eff = np.asarray(prepared[nm].effects, dtype=np.float64).copy()
            both = prepared[nm].target_available & pm_prepared.target_available
            rws = np.flatnonzero(both)
            pm_n = int(both.sum())
            if len(rws):
                eff[np.ix_(rws, shared)] = ((1 - plolo.PM_BLEND) * prepared[nm].effects[np.ix_(rws, shared)]
                                            + plolo.PM_BLEND * pm_prepared.effects[np.ix_(rws, shared)])
            sources.append(E062.PreparedSource(name="k562_gwps_pm50", effects=eff,
                                               target_available=np.asarray(prepared[nm].target_available, bool).copy(),
                                               mean_on_axis=np.asarray(prepared[nm].mean_on_axis, float).copy(),
                                               raw_median_norm=float(prepared[nm].raw_median_norm),
                                               mapped_gene_count=int(prepared[nm].mapped_gene_count),
                                               duplicate_symbols=tuple(prepared[nm].duplicate_symbols)))
        else:
            sources.append(prepared[nm])
    del prepared
    if "k562_gwps" in banks:
        del pm_prepared, pm_bank
    gc.collect()

    avail = np.stack([np.asarray(s.target_available, bool) for s in sources])
    pre64, votes = E062.build_baseline(sources)
    ka0_pre = float(np.abs(pre64.astype(np.float32).astype(np.float64) - pre_cache).max())
    del pre64
    cd4 = plolo.cd4_for(targets)
    union = L.union_rows([L.source_residual(nm, targets, "gene_norm") for nm in ("gse295214_h1", "h1_2025")])
    op = plolo.Operator(L.op_path_for(held))
    whitened = op.apply(pre_cache * ALPHA, wcols, GAMMA, LAM)
    d5, aud5 = L.deflate(cd4, "svd", 5)
    st = dict(held=held, fold=fold, targets=targets, n=n, wcols=wcols, sc=sc, truth=truth, pre=pre_cache,
              whitened=whitened, cd4=cd4, union=union, d5=d5, aud5=aud5,
              donors=[s.name for s in sources], sources=list(spec["sources"]), avail=avail, votes=votes,
              rows=[np.asarray(s.effects, dtype=np.float64) for s in sources], pm_n=pm_n,
              banks=banks, raw_full=raw_full, ka0=dict(pre_float32_max_abs_diff=ka0_pre))
    del op, sources
    gc.collect()
    return st


# ------------------------------------------------------------------------------------ the selection
def load_rel(line: str) -> dict:
    z = np.load(WORK / "rel" / f"{line}.npz", allow_pickle=True)
    return dict(names=[str(x) for x in z["names"]], r=np.asarray(z["r_nopanel"], dtype=np.float64),
                mag=np.asarray(z["mag"], dtype=np.float64), is_panel=np.asarray(z["is_panel"], bool))


def select_rows(bank, rel: dict, kind: str, q: float, seed: int) -> tuple[np.ndarray, dict]:
    """Row indices INTO THE REAL BANK for this arm, plus the audit.  `kind` in {full, sel, rev, mag, rand}."""
    pos = bank.pert_pos
    n_bank = bank.shape[0]
    names, r, mg, ip = rel["names"], rel["r"], rel["mag"], rel["is_panel"]
    if kind == "full":
        return np.arange(n_bank, dtype=np.int64), dict(
            kind=kind, q=1.0, seed=-1, n=int(n_bank), n_pool=0, n_selected_from_pool=0, n_panel=int(ip.sum()),
            n_bank=int(n_bank), key_median_selected=float("nan"), r_median_selected=float(np.nanmedian(r)),
            r_median_pool=float(np.nanmedian(r)), mag_median_selected=float(np.nanmedian(mg)),
            mag_median_pool=float(np.nanmedian(mg)))
    have = np.asarray([n in pos for n in names], bool)
    pool_m = have & (~ip) & np.isfinite(r)
    panel_m = have & ip
    pool = np.flatnonzero(pool_m)
    k = int(np.ceil(q * pool.size))
    if kind == "rand":
        rng = np.random.default_rng(seed)
        take = pool[rng.permutation(pool.size)[:k]]
    else:
        key = r[pool] if kind in ("sel", "rev") else mg[pool]
        nm = np.asarray([names[i] for i in pool])
        order = np.lexsort((nm, -key)) if kind in ("sel", "mag") else np.lexsort((nm, key))
        take = pool[order[:k]]
    chosen = np.concatenate([take, np.flatnonzero(panel_m)])
    idx = np.unique(np.asarray([pos[names[i]] for i in chosen], dtype=np.int64))
    aud = dict(kind=kind, q=q, seed=seed, n=int(idx.size), n_pool=int(pool.size), n_selected_from_pool=int(k),
               n_panel=int(panel_m.sum()), n_bank=int(n_bank),
               key_median_selected=float(np.median((r if kind != "mag" else mg)[take])) if take.size else float("nan"),
               r_median_selected=float(np.median(r[take])) if take.size else float("nan"),
               r_median_pool=float(np.median(r[pool])), mag_median_selected=float(np.median(mg[take])) if take.size else float("nan"),
               mag_median_pool=float(np.median(mg[pool])))
    return idx, aud


# ------------------------------------------------------------------------------------ the fit (062's recipe, rows subset)
def row_weights(bank, rel: dict, kind: str, seed: int) -> tuple[np.ndarray, np.ndarray, dict]:
    """POST-HOC (not preregistered) rank-preserving arms: EVERY row of the bank teaches the operator, but each row is
    weighted by its split-half reliability.  Row count, donor proportions, axis, gamma / lambda are the incumbent's,
    so the covariance keeps the incumbent's rank -- the confound the q < 1 arms cannot avoid.

      wgt:mild     w = 0.5 + u,  u = the row's rank quantile of R in [0, 1]  -- a BOUNDED monotone tilt toward the
                   reliable rows that leaves the effective sample size at ~96 % of n, so any movement is composition
                   and not a loss of rank (the q < 1 arms cannot separate the two)
      wgt:mildrev  w = 1.5 - u                     the reversed tilt, at the same weight distribution
      wgt:R        w ∝ clip(R, 0.005, 0.99)        the strong dose (n_eff 34-44 % of n)
      wgt:rev      w ∝ clip(1 - R, 0.01, 0.995)    the reversed strong dose
      wgt:shuf     the wgt:R weights permuted across the bank's rows (the composition control)
    A row with no reliability takes the bank's median R.  Weights are renormalised per bank to sum to the bank's row
    count, so the covariance is the incumbent's donor mixture exactly.
    """
    import est271 as EST
    pos, n_bank = bank.pert_pos, bank.shape[0]
    R = np.full(n_bank, np.nan)
    for nm, rv in zip(rel["names"], rel["r"]):
        j = pos.get(nm)
        if j is not None:
            R[j] = rv
    R = EST.spearman_brown(R)
    med = float(np.nanmedian(R))
    R = np.where(np.isfinite(R), R, med)
    Rc = np.clip(R, EST.R_MIN, EST.R_MAX)
    u = np.argsort(np.argsort(R, kind="stable")).astype(np.float64) / max(n_bank - 1, 1)
    if kind in ("wgtR", "wgtshuf"):
        w = Rc.copy()
    elif kind == "wgtrev":
        w = np.clip(1.0 - R, 1.0 - EST.R_MAX, 1.0 - EST.R_MIN)
    elif kind == "wgtmild":
        w = 0.5 + u
    elif kind == "wgtmildrev":
        w = 1.5 - u
    elif kind == "wgtsteep":                       # the same tilt at 4.75x the mild dose (w in [0.1, 1.9])
        w = 0.1 + 1.8 * u
    elif kind == "wgtsteeprev":
        w = 1.9 - 1.8 * u
    else:
        raise ValueError(kind)
    if kind == "wgtshuf":
        w = w[np.random.default_rng(seed).permutation(n_bank)]
    w = w * (n_bank / w.sum())
    neff = float(w.sum() ** 2 / (w * w).sum())
    aud = dict(kind=kind, seed=seed, n=int(n_bank), n_with_r=int(np.isfinite(R).sum()), median_R=med,
               n_eff=neff, n_eff_over_n=neff / n_bank, w_min=float(w.min()), w_max=float(w.max()),
               w_p10=float(np.percentile(w, 10)), w_p90=float(np.percentile(w, 90)))
    return np.arange(n_bank, dtype=np.int64), w, aud


def fit299(banks: list, rows_by: dict, raw_medians: dict, wcols: np.ndarray, *, weights_by: dict | None = None, log=P) -> tuple[np.ndarray, np.ndarray, dict]:
    """`op272.fit_operator` restricted to a row index set per bank, held in memory (no npz).  The chosen rows are
    centred on THEIR OWN column mean (062's recipe applied to the chosen set; identical to the bank mean when the
    set is every row).  scipy `eigh(overwrite_a=True, driver='evd')` for the peak-RSS cap; KA-2 validates it."""
    t0 = time.time()
    p = int(len(wcols))
    ref = float(np.median([raw_medians[b.name] for b in banks]))
    C = np.zeros((p, p), dtype=np.float64)
    total = 0
    rows_by_name, mean_shift = {}, {}
    for bank in banks:
        idx = np.asarray(rows_by[bank.name], dtype=np.int64)
        w = None if weights_by is None else np.asarray(weights_by[bank.name], dtype=np.float64)
        wsum = float(idx.size) if w is None else float(w.sum())
        src_cols, op_cols = E062.last_source_columns(bank, OFFICIAL, np.asarray(wcols, dtype=np.int64))
        delta = np.load(bank.delta_path, mmap_mode="r")
        acc = np.zeros(src_cols.size, dtype=np.float64)
        for lo in range(0, idx.size, CHUNK):
            blk = np.asarray(delta[idx[lo:lo + CHUNK]][:, src_cols], dtype=np.float64)
            acc += blk.sum(0) if w is None else (blk * w[lo:lo + CHUNK, None]).sum(0)
            del blk
        mean = acc / wsum
        mult = ref / raw_medians[bank.name]
        for lo in range(0, idx.size, CHUNK):
            hi = min(lo + CHUNK, idx.size)
            blk = np.zeros((hi - lo, p), dtype=np.float64)
            blk[:, op_cols] = (np.asarray(delta[idx[lo:hi]][:, src_cols], dtype=np.float64) - mean) * mult
            if w is not None:
                blk *= np.sqrt(w[lo:hi])[:, None]
            C += blk.T @ blk
            del blk
        del delta
        total += wsum
        rows_by_name[bank.name] = int(idx.size)
        mean_shift[bank.name] = float(np.linalg.norm(mean))
        log(f"    gram {bank.name}: {idx.size} rows x {src_cols.size} cols (mult {mult:.6f}, sum w {wsum:.1f}), "
            f"{time.time()-t0:.0f}s, rss {rss():.2f} GiB")
    C /= total
    log(f"    covariance {C.shape} over {total} rows; eigh (scipy evd, overwrite) ... rss {rss():.2f} GiB")
    evals, evecs = sla.eigh(C, overwrite_a=True, driver="evd")
    del C
    evals = np.maximum(np.asarray(evals, dtype=np.float64), float(np.max(evals)) * 1e-12)
    evecs = evecs.astype(np.float32)          # 062 stores float32; plolo.Operator loads it back to float64
    gc.collect()
    evecs = np.asarray(evecs, dtype=np.float64)
    gc.collect()
    floor = float(evals.max()) * 1e-12
    meta = dict(n_cols=p, basis_rows_by_source=rows_by_name, basis_total_rows=int(total), reference_norm=ref,
                raw_medians={b.name: raw_medians[b.name] for b in banks}, mean_norm_by_source=mean_shift,
                n_floored=int((evals <= floor * (1 + 1e-9)).sum()),
                trace=float(evals.sum()), top_decile_share=float(np.sort(evals)[-(p // 10):].sum() / evals.sum()),
                max_eigenvalue=float(evals.max()), min_eigenvalue=float(evals.min()),
                evals_sha256=hashlib.sha256(np.ascontiguousarray(evals).tobytes()).hexdigest()[:16],
                seconds=round(time.time() - t0, 1), peak_rss_gib=rss())
    log(f"    fitted: max eig {meta['max_eigenvalue']:.6f}, evals sha {meta['evals_sha256']}, {meta['seconds']}s, peak {rss():.2f} GiB")
    return evals, evecs, meta


# ------------------------------------------------------------------------------------ the measurement
def main(held: str, arms: list[str]) -> None:
    t0 = time.time()
    sweep = WORK / f"sweep299-{held}.json"
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    if sweep.exists():
        old = json.loads(sweep.read_text())
        rows = old["rows"]
        done = {r["arm"] for r in rows}
        z = np.load(WORK / f"pertarget299-{held}.npz")
        per_target = {k: z[k] for k in z.files if k not in ("targets", "avail")}
        z = np.load(WORK / f"identity299-{held}.npz")
        identity = {k: z[k] for k in z.files}
    else:
        done = set()
    arms = [a for a in arms if a not in done]
    P(f"[{held}] run299.py sha256 {sha256_self()[:16]}; arms {arms}; already done {sorted(done)}")

    st = stage299(held)
    ident = X.Identity(st)
    n, wc, sc = st["n"], st["wcols"], st["sc"]
    targets = st["targets"]
    pre = np.asarray(st["pre"], dtype=np.float64)
    donors = list(st["sources"])
    P(f"[{held}] staged {n} targets, donors {donors}, wcols {wc.size}, {time.time()-t0:.0f}s, rss {rss():.2f} GiB")
    ka_log: list[dict] = []
    failures: list[str] = []

    def check(label, got, want, tol):
        d = abs(float(got) - float(want))
        ok = d <= tol
        ka_log.append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
        if not ok:
            failures.append(label)
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")

    check("KA-0 build_baseline(per-donor rows) float32 == cached pre", st["ka0"]["pre_float32_max_abs_diff"], 0.0, 0.0)
    ref263 = {r["arm"]: r for r in json.loads((E263 / f"work/sweep263-{held}.json").read_text())["rows"]}
    pt263 = np.load(E263 / f"work/pertarget263-{held}.npz")
    rel = {nm: load_rel(nm) for nm in donors}
    banks = [st["banks"][nm] for nm in donors]

    def flush():
        (WORK / f"sweep299-{held}.json").write_text(json.dumps(dict(
            meta=dict(fold=held, n_targets=n, targets=targets, donors=donors, run299_sha256=sha256_self(),
                      op272_sha256=O.sha256_self(), n_wcols=int(wc.size), pm_n=st["pm_n"],
                      ka=ka_log, ka_pass=sum(1 for k in ka_log if k["ok"]), ka_total=len(ka_log), failures=failures,
                      seconds=round(time.time() - t0, 1), peak_rss_gib=rss()), rows=rows), indent=1, default=str))
        np.savez_compressed(WORK / f"pertarget299-{held}.npz", targets=np.asarray(targets),
                            avail=st["avail"].astype(np.int8), **per_target)
        np.savez_compressed(WORK / f"identity299-{held}.npz", **identity)

    def evaluate(name: str, base: np.ndarray, aud_extra: dict) -> dict:
        t1 = time.time()
        parts, aud = Y.arm_parts(st, base)
        st_arm = dict(st, whitened=base)
        final = L.compose(st_arm, parts)
        e = L.score(st_arm, final)
        m = ident.measure(final)
        per_target[f"{name}|pds"] = np.asarray(e["_pds_per_target_a23"])
        per_target[f"{name}|cos"] = np.asarray(e["_cos_per_target_a23"])
        for key in ("U", "M", "D"):
            identity[f"{name}|{key}"] = m[key]
        aud["base_median_norm_wcols"] = float(np.median(np.linalg.norm(base[:, wc], axis=1)))
        aud["base_cos_to_incumbent_base_median"] = float(np.nanmedian(Y.row_cos(base[:, wc], st["whitened"][:, wc])))
        row = dict(fold=held, arm=name, n_targets=n, **{k: v for k, v in e.items() if not k.startswith("_")},
                   nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"], delivered_norm=m["delivered_norm"],
                   audit=dict(**aud, **aud_extra), seconds=round(time.time() - t1, 1))
        rows.append(row)
        P(f"  {name:12s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
          f"nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} base~inc {aud['base_cos_to_incumbent_base_median']:.4f} "
          f"({row['seconds']}s, rss {rss():.2f} GiB)")
        return row

    for arm in arms:
        try:
            if arm == "cd4svd5":
                r = evaluate(arm, st["whitened"], dict(control="the incumbent, cached operator"))
                d = max(abs(float(ref263[arm][k]) - float(r[k])) for k in KEYS)
                check(f"KA-3 {arm} == 263 on {KEYS}", d, 0.0, TOL)
                check("KA-4 cd4svd5 per-target pds bit-identical to 263",
                      float(np.abs(per_target["cd4svd5|pds"] - pt263["cd4svd5|pds"]).max()), 0.0, 0.0)
            elif arm == "STACK":
                parts = [(st["cd4"], None), (st["union"], None)]
                fr, _ = plolo.solve_fraction(st["whitened"], st["cd4"], Y.CD4_F)
                fx, _ = Y.solve_fraction_carrying(st["whitened"], st["union"], Y.H1_G)
                st_arm = dict(st, whitened=st["whitened"])
                final = L.compose(st_arm, [(st["cd4"], float(fr)), (st["union"], float(fx))])
                e = L.score(st_arm, final)
                m = ident.measure(final)
                per_target[f"{arm}|pds"] = np.asarray(e["_pds_per_target_a23"])
                per_target[f"{arm}|cos"] = np.asarray(e["_cos_per_target_a23"])
                for key in ("U", "M", "D"):
                    identity[f"{arm}|{key}"] = m[key]
                rows.append(dict(fold=held, arm=arm, n_targets=n, **{k: v for k, v in e.items() if not k.startswith("_")},
                                 nmae_wald89=m["nmae_wald89"], top89_precision=m["top89_precision"],
                                 delivered_norm=m["delivered_norm"],
                                 audit=dict(cd4_fraction=float(fr), union_fraction=float(fx),
                                            arm="h1union: un-deflated CD4 at f=0.34 + pair union at g=0.28"), seconds=0.0))
                P(f"  {arm:12s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} nmae89 {m['nmae_wald89']:.5f}")
                d = max(abs(float(ref263[arm][k]) - float(e[k])) for k in KEYS) if arm in ref263 else float("nan")
                if arm in ref263:
                    check(f"KA-3b {arm} == 263 on {KEYS}", d, 0.0, TOL)
            else:
                spec = arm.split(":")
                kind = spec[0]
                seed = 299000 + (["k562_gwps", "orion_hct116", "orion_hek293t"].index(held) if held in
                                 ("k562_gwps", "orion_hct116", "orion_hek293t") else 9)
                for tok in spec[2:]:
                    if tok.startswith("s"):
                        seed = int(tok[1:])
                rows_by, auds, wts = {}, {}, None
                if kind == "wgt":
                    wts = {}
                    for b in banks:
                        idx, w, a = row_weights(b, rel[b.name], f"wgt{spec[1]}", seed)
                        rows_by[b.name], wts[b.name], auds[b.name] = idx, w, a
                        P(f"    {arm} {b.name}: {a['n']} rows, all of them; n_eff {a['n_eff']:.0f} "
                          f"({100*a['n_eff_over_n']:.1f} % of n), w in [{a['w_min']:.4f}, {a['w_max']:.4f}], "
                          f"p10 {a['w_p10']:.4f} p90 {a['w_p90']:.4f}, median R {a['median_R']:+.4f}")
                else:
                    q = 1.0 if kind == "full" else float(spec[1][1:])
                    for b in banks:
                        idx, a = select_rows(b, rel[b.name], kind, q, seed)
                        rows_by[b.name] = idx
                        auds[b.name] = a
                        P(f"    {arm} {b.name}: {a['n']} rows (pool {a['n_pool']}, took {a['n_selected_from_pool']}, "
                          f"panel {a['n_panel']}), r median selected {a['r_median_selected']:+.4f} vs pool {a['r_median_pool']:+.4f}, "
                          f"mag median {a['mag_median_selected']:.3f} vs {a['mag_median_pool']:.3f}")
                evals, evecs, fmeta = fit299(banks, rows_by, st["raw_full"], wc, weights_by=wts)
                dop = O.DampedOperator(evals, evecs, None, 0.0, arm)
                base = O.preserve(pre * ALPHA, dop.transform(pre[:, wc] * ALPHA, GAMMA, LAM), wc)
                if kind == "full":
                    regp = plolo.OPERATOR_REGISTRY[f"held_{held}"]
                    with np.load(ROOT / regp["path"]) as z:
                        ev_c = np.sort(np.asarray(z["evals"], dtype=np.float64))
                    ev_r = np.sort(evals)
                    relerr = np.abs(ev_r - ev_c) / np.maximum(np.abs(ev_c), 1e-300)
                    big = ev_c > ev_c.max() * 1e-8
                    check("KA-2a eigenvalues of the full-row refit == cached operator (max rel, eig > 1e-8 max)",
                          float(relerr[big].max()), 0.0, 1e-6)
                    rr = np.linalg.norm((base - st["whitened"])[:, wc], axis=1) / np.linalg.norm(st["whitened"][:, wc], axis=1)
                    check("KA-2e full-row refit applied to pre == the cached operator's base (max row rel)", float(rr.max()), 0.0, 1e-4)
                del dop, evecs
                gc.collect()
                r = evaluate(arm, base, dict(fit=fmeta, selection=auds))
                del base
                gc.collect()
                if kind == "full":
                    d = max(abs(float(ref263["cd4svd5"][k]) - float(r[k])) for k in KEYS)
                    check("KA-2f full == cd4svd5 on (cos_agg, macro, pds)", d, 0.0, 1e-5)
            flush()
        except Exception as exc:                      # keep every completed arm
            failures.append(f"{arm}: {type(exc).__name__}: {exc}")
            P(f"  ARM FAILED {arm}: {type(exc).__name__}: {exc}")
            flush()
            raise
    flush()
    (WORK / f"run299-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{held}] done: {len(rows)} rows, KA {sum(1 for k in ka_log if k['ok'])}/{len(ka_log)}, failures {failures}, "
      f"{(time.time()-t0)/60:.1f} min, peak {rss():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
