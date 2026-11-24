"""300 -- the SPREAD mechanism: fit the fold's whitening operator with rank-preserving row WEIGHTS (or on 062's
production prune row set) and read the arm through `cd4svd5`'s construction.  ONE FOLD PER PROCESS.  Every arm's
row is written as soon as it is measured, so a killed process loses only the arm in flight.

    PYTHONDONTWRITEBYTECODE=1 uv run python run300.py orion_hct116 cd4svd5 STACK tilt:d1.0 tilt:d-1.0 ...

Arms (PREREGISTRATION.md section 1; every arm keeps the incumbent's row count and donor proportions except `prune:*`)
  cd4svd5           the incumbent through the CACHED fold operator (KA-3 / KA-4 against 263)
  STACK             h1union -- the BR249 retrodiction leg
  tilt:d<D>         A: w = 1 + D*(0.5 - u), u = the rank quantile of the row's split-half reliability.
                    D>0 tilts AWAY from reliable rows; D<0 is the MIRROR control.  D=0 is the un-weighted refit
                    (KA-2a/2e/2f against the cached operator).  D=1.0/-1.0/1.8/-1.8 are 299's
                    wgt:mildrev / wgt:mild / wgt:steeprev / wgt:steep (KA-W1/W2/W3).
  spread:d<D>       C: w = 1 + D*(0.5 - v), v = the rank quantile of the row's energy share on the top 50
                    eigenvectors of the fold's own CACHED (un-weighted) operator.  D>0 up-weights LOW-leverage
                    rows.  Reliability is not used.
  prune:d<D>        B: the same tilt as `tilt` but on 062's PRODUCTION row set (12,000-row per-donor cap at
                    seed 0, then the common-essential rows dropped).  D=0 is the prune alone.
  rand:d<D>:s<SEED> D(ii): the tilt weights with u PERMUTED across the bank's rows.

Everything downstream of the base is `cd4svd5` verbatim (lib263 / lib246 / lib257).  ONE VARIABLE: what teaches
the operator.
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
HERE = ROOT / "experiments/300-spread-tilt"
WORK = HERE / "work"
E262 = ROOT / "experiments/262-pds-loss-diagnosis"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
E272 = ROOT / "experiments/272-noise-corrected-operator"
E299 = ROOT / "experiments/299-reliable-training-rows"
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
N_PC = 50                                    # arm C: the "leading directions" of the bank
PRUNE_CAP = 12000                            # 062 / scripts/build_union_whitener.py
PRUNE_SEED = 0
DONOR_ORDER = ("k562_gwps", "orion_hct116", "orion_hek293t")
ESSENTIAL = ROOT / "data/donors/k562_essential/perts.txt"
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_self() -> str:
    return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()


# ------------------------------------------------------------------------------------ staging (299's stage299)
def stage300(held: str) -> dict:
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


def load_rel(line: str) -> dict:
    z = np.load(WORK / "rel" / f"{line}.npz", allow_pickle=True)
    return dict(names=[str(x) for x in z["names"]], r=np.asarray(z["r_nopanel"], dtype=np.float64),
                mag=np.asarray(z["mag"], dtype=np.float64), is_panel=np.asarray(z["is_panel"], bool))


# ------------------------------------------------------------------------------------ 062's production row set
def prune_rows(fold_banks: list) -> dict[str, np.ndarray]:
    """062's basis rows, reproduced exactly: a 12,000-row cap per donor drawn by ONE `default_rng(0)` visiting
    the donors in DONOR_ORDER (donors under the cap consume no random numbers), then every row whose perturbation
    is in the common-essential list dropped.  The RNG stream is consumed for ALL THREE production donors so that
    a two-bank fold sees the SAME subsample the production three-bank operator saw."""
    ess = set(ESSENTIAL.read_text().split("\n")) - {""}
    rng = np.random.default_rng(PRUNE_SEED)
    capped: dict[str, np.ndarray] = {}
    for nm in DONOR_ORDER:
        b = E062.Bank.load(nm, plolo.DONOR_PATHS[nm])
        n = len(b.perts)
        capped[nm] = np.sort(rng.choice(n, PRUNE_CAP, False)) if n > PRUNE_CAP else np.arange(n, dtype=np.int64)
    out = {}
    for b in fold_banks:
        idx = capped[b.name]
        keep = np.asarray([b.perts[i] not in ess for i in idx], bool)
        out[b.name] = np.asarray(idx[keep], dtype=np.int64)
    return out


# ------------------------------------------------------------------------------------ weights
def _rank_quantile(x: np.ndarray) -> np.ndarray:
    n = x.size
    return np.argsort(np.argsort(x, kind="stable")).astype(np.float64) / max(n - 1, 1)


def tilt_u(bank, rel: dict, rows: np.ndarray) -> np.ndarray:
    """The rank quantile of the Spearman-Brown-corrected split-half reliability, within the bank's chosen rows.
    299's `row_weights` verbatim (median R for a row with no measurement), restricted to `rows`."""
    import est271 as EST
    pos, n_bank = bank.pert_pos, len(bank.perts)
    R = np.full(n_bank, np.nan)
    for nm, rv in zip(rel["names"], rel["r"]):
        j = pos.get(nm)
        if j is not None:
            R[j] = rv
    R = EST.spearman_brown(R)
    R = R[rows]
    med = float(np.nanmedian(R))
    R = np.where(np.isfinite(R), R, med)
    return _rank_quantile(R)


def spread_v(st: dict, bank, rows: np.ndarray, wcols: np.ndarray, V50: np.ndarray) -> tuple[np.ndarray, dict]:
    """Arm C's axis: the row's ENERGY SHARE on the top-`N_PC` eigenvectors of the fold's own cached (un-weighted)
    operator, measured on exactly the rows the un-weighted fit sees (062's centring and `mult`), then rank-
    quantiled within the bank.  Low share = low leverage on the bank's leading directions."""
    ref = float(np.median([st["raw_full"][b.name] for b in st["fit_banks"]]))
    mult = ref / st["raw_full"][bank.name]
    src_cols, op_cols = E062.last_source_columns(bank, OFFICIAL, np.asarray(wcols, dtype=np.int64))
    delta = np.load(bank.delta_path, mmap_mode="r")
    acc = np.zeros(src_cols.size, dtype=np.float64)
    for lo in range(0, rows.size, CHUNK):
        acc += np.asarray(delta[rows[lo:lo + CHUNK]][:, src_cols], dtype=np.float64).sum(0)
    mean = acc / rows.size
    Vsub = np.ascontiguousarray(V50[op_cols])                # (n_src, N_PC)
    share = np.empty(rows.size, dtype=np.float64)
    tot = np.empty(rows.size, dtype=np.float64)
    for lo in range(0, rows.size, CHUNK):
        hi = min(lo + CHUNK, rows.size)
        z = (np.asarray(delta[rows[lo:hi]][:, src_cols], dtype=np.float64) - mean) * mult
        proj = z @ Vsub
        tot[lo:hi] = np.einsum("ij,ij->i", z, z)
        share[lo:hi] = np.einsum("ij,ij->i", proj, proj)
        del z, proj
    del delta
    s = share / np.maximum(tot, 1e-300)
    aud = dict(n=int(rows.size), n_pc=N_PC, share_median=float(np.median(s)), share_p10=float(np.percentile(s, 10)),
               share_p90=float(np.percentile(s, 90)), share_min=float(s.min()), share_max=float(s.max()))
    return _rank_quantile(s), aud


def dose_weights(q: np.ndarray, d: float) -> np.ndarray:
    """w = 1 + d*(0.5 - q).  Sum(0.5 - q) == 0 exactly over a rank quantile, so sum(w) == n exactly: the row
    count, donor proportions and coverage are the un-weighted fit's on every dose."""
    return 1.0 + d * (0.5 - q)


# ------------------------------------------------------------------------------------ the fit (062's recipe)
def fit300(banks: list, rows_by: dict, raw_medians: dict, wcols: np.ndarray, *, weights_by: dict | None = None, log=P):
    """299's `fit299` verbatim: `op272.fit_operator` restricted to a row index set per bank and an optional row
    weight, held in memory.  The chosen rows are centred on THEIR OWN (weighted) column mean.  scipy
    `eigh(overwrite_a=True, driver='evd')` for the peak-RSS cap; KA-2 validates it."""
    t0 = time.time()
    p = int(len(wcols))
    ref = float(np.median([raw_medians[b.name] for b in banks]))
    C = np.zeros((p, p), dtype=np.float64)
    total = 0.0
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
    meta = dict(n_cols=p, basis_rows_by_source=rows_by_name, basis_total_rows=float(total), reference_norm=ref,
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
    sweep = WORK / f"sweep300-{held}.json"
    rows: list[dict] = []
    per_target: dict[str, np.ndarray] = {}
    identity: dict[str, np.ndarray] = {}
    if sweep.exists():
        old = json.loads(sweep.read_text())
        rows = old["rows"]
        done = {r["arm"] for r in rows}
        z = np.load(WORK / f"pertarget300-{held}.npz")
        per_target = {k: z[k] for k in z.files if k not in ("targets", "avail")}
        z = np.load(WORK / f"identity300-{held}.npz")
        identity = {k: z[k] for k in z.files}
    else:
        done = set()
    arms = [a for a in arms if a not in done]
    P(f"[{held}] run300.py sha256 {sha256_self()[:16]}; arms {arms}; already done {sorted(done)}")

    st = stage300(held)
    ident = X.Identity(st)
    n, wc, sc = st["n"], st["wcols"], st["sc"]
    targets = st["targets"]
    pre = np.asarray(st["pre"], dtype=np.float64)
    donors = list(st["sources"])
    banks = [st["banks"][nm] for nm in donors]
    st["fit_banks"] = banks
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
    ref299 = {r["arm"]: r for r in json.loads((E299 / f"work/sweep299-{held}.json").read_text())["rows"]}
    rel = {nm: load_rel(nm) for nm in donors}

    # caches that are computed at most once per process
    cache: dict = {}

    def full_rows():
        return {b.name: np.arange(len(b.perts), dtype=np.int64) for b in banks}

    def get_prune():
        if "prune" not in cache:
            cache["prune"] = prune_rows(banks)
            for b in banks:
                P(f"    062 prune {b.name}: {cache['prune'][b.name].size} rows kept of {len(b.perts)}")
        return cache["prune"]

    def get_u(tag: str, rows_by: dict):
        key = f"u:{tag}"
        if key not in cache:
            cache[key] = {b.name: tilt_u(b, rel[b.name], rows_by[b.name]) for b in banks}
        return cache[key]

    def get_v(rows_by: dict):
        if "v" not in cache:
            regp = plolo.OPERATOR_REGISTRY[f"held_{held}"]
            with np.load(ROOT / regp["path"]) as z:
                ev = np.asarray(z["evals"], dtype=np.float64)
                top = np.argsort(ev)[-N_PC:][::-1]
                V50 = np.asarray(z["evecs"], dtype=np.float64)[:, top]
            P(f"    arm C basis: top {N_PC} eigenvectors of the cached operator, eigenvalues "
              f"{ev[top][0]:.4f} .. {ev[top][-1]:.4f} ({100*ev[top].sum()/ev.sum():.2f} % of the trace)")
            out, aud = {}, {}
            for b in banks:
                out[b.name], aud[b.name] = spread_v(st, b, rows_by[b.name], wc, V50)
                P(f"    arm C {b.name}: energy share on the top {N_PC} PCs median {aud[b.name]['share_median']:.4f} "
                  f"(p10 {aud[b.name]['share_p10']:.4f}, p90 {aud[b.name]['share_p90']:.4f})")
            del V50
            gc.collect()
            cache["v"] = (out, aud)
        return cache["v"]

    def flush():
        (WORK / f"sweep300-{held}.json").write_text(json.dumps(dict(
            meta=dict(fold=held, n_targets=n, targets=targets, donors=donors, run300_sha256=sha256_self(),
                      op272_sha256=O.sha256_self(), n_wcols=int(wc.size), pm_n=st["pm_n"],
                      ka=ka_log, ka_pass=sum(1 for k in ka_log if k["ok"]), ka_total=len(ka_log), failures=failures,
                      seconds=round(time.time() - t0, 1), peak_rss_gib=rss()), rows=rows), indent=1, default=str))
        np.savez_compressed(WORK / f"pertarget300-{held}.npz", targets=np.asarray(targets),
                            avail=st["avail"].astype(np.int8), **per_target)
        np.savez_compressed(WORK / f"identity300-{held}.npz", **identity)

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
        P(f"  {name:18s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} "
          f"nmae89 {m['nmae_wald89']:.5f} dn {m['delivered_norm']:.4f} base~inc {aud['base_cos_to_incumbent_base_median']:.4f} "
          f"({row['seconds']}s, rss {rss():.2f} GiB)")
        return row

    def ka_vs_299(arm299: str, arm: str, r: dict) -> None:
        """the mandate's item 1: 299's wgt rows reproduced through THIS experiment's plumbing at <= 1e-9"""
        if arm299 not in ref299 or "cd4svd5" not in ref299:
            P(f"  KA SKIP {arm299}: not in 299's sweep for {held}")
            return
        d299 = float(ref299[arm299]["pds"]) - float(ref299["cd4svd5"]["pds"])
        d300 = float(r["pds"]) - float([x for x in rows if x["arm"] == "cd4svd5"][0]["pds"])
        check(f"KA-W {arm} d pds == 299's {arm299} ({d299:+.8f})", d300, d299, TOL)

    KA299 = {"tilt:d1.0": "wgt:mildrev", "tilt:d-1.0": "wgt:mild", "tilt:d1.8": "wgt:steeprev",
             "tilt:d-1.8": "wgt:steep"}

    for arm in arms:
        try:
            if arm == "cd4svd5":
                r = evaluate(arm, st["whitened"], dict(control="the incumbent, cached operator"))
                d = max(abs(float(ref263[arm][k]) - float(r[k])) for k in KEYS)
                check(f"KA-3 {arm} == 263 on {KEYS}", d, 0.0, TOL)
                check("KA-4 cd4svd5 per-target pds bit-identical to 263",
                      float(np.abs(per_target["cd4svd5|pds"] - pt263["cd4svd5|pds"]).max()), 0.0, 0.0)
            elif arm == "STACK":
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
                P(f"  {arm:18s} cos_agg {e['cos_agg']:.6f} macro {e['a23_cos']:.6f} pds {e['pds']:.6f} nmae89 {m['nmae_wald89']:.5f}")
                if arm in ref263:
                    check(f"KA-3b {arm} == 263 on {KEYS}",
                          max(abs(float(ref263[arm][k]) - float(e[k])) for k in KEYS), 0.0, TOL)
            else:
                spec = arm.split(":")
                kind = spec[0]
                d = float(spec[1][1:])
                seed = 300001
                for tok in spec[2:]:
                    if tok.startswith("s"):
                        seed = int(tok[1:])
                if kind == "prune":
                    rows_by = get_prune()
                    q = get_u("prune", rows_by)
                elif kind in ("tilt", "rand"):
                    rows_by = full_rows()
                    q = get_u("full", rows_by)
                    if kind == "rand":
                        rg = np.random.default_rng(seed)
                        q = {b.name: q[b.name][rg.permutation(q[b.name].size)] for b in banks}
                elif kind == "spread":
                    rows_by = full_rows()
                    q, vaud = get_v(rows_by)
                else:
                    raise ValueError(arm)
                wts, auds = None, {}
                if d != 0.0:
                    wts = {}
                    for b in banks:
                        w = dose_weights(q[b.name], d)
                        wts[b.name] = w
                        neff = float(w.sum() ** 2 / (w * w).sum())
                        auds[b.name] = dict(kind=kind, dose=d, seed=seed, n=int(w.size), n_eff=neff,
                                            n_eff_over_n=neff / w.size, w_min=float(w.min()), w_max=float(w.max()),
                                            w_sum=float(w.sum()))
                        P(f"    {arm} {b.name}: {w.size} rows, sum w {w.sum():.4f} (n {w.size}), n_eff {neff:.0f} "
                          f"({100*neff/w.size:.1f} %), w in [{w.min():.4f}, {w.max():.4f}]")
                else:
                    for b in banks:
                        auds[b.name] = dict(kind=kind, dose=0.0, n=int(rows_by[b.name].size))
                evals, evecs, fmeta = fit300(banks, rows_by, st["raw_full"], wc, weights_by=wts)
                dop = O.DampedOperator(evals, evecs, None, 0.0, arm)
                base = O.preserve(pre * ALPHA, dop.transform(pre[:, wc] * ALPHA, GAMMA, LAM), wc)
                if kind == "tilt" and d == 0.0:
                    regp = plolo.OPERATOR_REGISTRY[f"held_{held}"]
                    with np.load(ROOT / regp["path"]) as z:
                        ev_c = np.sort(np.asarray(z["evals"], dtype=np.float64))
                    ev_r = np.sort(evals)
                    relerr = np.abs(ev_r - ev_c) / np.maximum(np.abs(ev_c), 1e-300)
                    big = ev_c > ev_c.max() * 1e-8
                    check("KA-2a eigenvalues of the un-weighted refit == cached operator (max rel, eig > 1e-8 max)",
                          float(relerr[big].max()), 0.0, 1e-6)
                    rr = np.linalg.norm((base - st["whitened"])[:, wc], axis=1) / np.linalg.norm(st["whitened"][:, wc], axis=1)
                    check("KA-2e un-weighted refit applied to pre == the cached operator's base (max row rel)", float(rr.max()), 0.0, 1e-4)
                del dop, evecs
                gc.collect()
                extra = dict(fit=fmeta, selection=auds)
                if kind == "spread":
                    extra["spread_axis"] = vaud
                r = evaluate(arm, base, extra)
                del base
                gc.collect()
                if kind == "tilt" and d == 0.0:
                    check("KA-2f un-weighted refit == cd4svd5 on (cos_agg, macro, pds)",
                          max(abs(float(ref263["cd4svd5"][k]) - float(r[k])) for k in KEYS), 0.0, 1e-5)
                if arm in KA299:
                    ka_vs_299(KA299[arm], arm, r)
            flush()
        except Exception as exc:                      # keep every completed arm
            failures.append(f"{arm}: {type(exc).__name__}: {exc}")
            P(f"  ARM FAILED {arm}: {type(exc).__name__}: {exc}")
            flush()
            (WORK / f"run300-{held}.log").write_text("\n".join(LOG) + "\n")
            raise
        (WORK / f"run300-{held}.log").write_text("\n".join(LOG) + "\n")
    flush()
    (WORK / f"run300-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{held}] done: {len(rows)} rows, KA {sum(1 for k in ka_log if k['ok'])}/{len(ka_log)}, failures {failures}, "
      f"{(time.time()-t0)/60:.1f} min, peak {rss():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
