"""286 STAGE C -- the DIRECTION-RESOLVING reference: the same cross-product arms composed through 273's h1-test
instrument with the H1 union OFF (base + H1-CD4 deflated k = 5 re-solved at f = 0.34), per-target `pds` and own cosine
against the 10x Flex, 400-cell truth, with 265's paired target bootstrap.

    PYTHONDONTWRITEBYTECODE=1 uv run python stageC286.py

Everything is 273's construction verbatim (lib273): 170's staged `cand_h1`, `union_pruned`, the PerturbMatch 0.5 blend
variant 273 selected, the re-extracted H1 CD4 bank, `plolo.deliver` + the scorer's read-back, `resid218._per_target`.
The only new object is the donor row: line j's row filtered by the cross-product estimator fitted in Stage B.

Known answers: BASE and INCUMBENT reproduce 273's `instrument273.json` arms; the E0 rebuild reproduces INCUMBENT.
LEAKAGE NOTE: the Stage B fits exclude the 300 A/B/C panel targets but not h1-test's 85 targets; the overlap with each
line's fitted bank rows is measured and reported below (order 85 / 17,500 rows, i.e. ~0.5 % of a covariance fitted on
17.5 k rows).
"""
from __future__ import annotations

import gc
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/286-crossproduct-estimator"
WORK = HERE / "work"
E273 = ROOT / "experiments/273-flex-reference-instrument"
sys.path.insert(0, str(E273))
sys.path.insert(0, str(HERE))
import lib273 as X  # noqa: E402
import est286 as E  # noqa: E402
from boot265 import boot  # noqa: E402

plolo, L, E062, E064 = X.plolo, X.L, X.E062, X.E064
OFFICIAL, PANEL_SET = X.OFFICIAL, X.PANEL_SET
ALPHA, GAMMA, LAM, CD4_F = X.ALPHA, X.GAMMA, X.LAM, X.CD4_F
CD4_H1 = E273 / "work/cd4/cd4_compact_h1"
LOG: list[str] = []
OUT: dict = dict(est286_sha256=E.sha256_self(), ka=[], failures=[])
rng = np.random.default_rng(20260906)


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def rss() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def check(label, got, want, tol):
    d = abs(float(got) - float(want))
    ok = d <= tol
    OUT["ka"].append(dict(label=label, got=float(got), want=float(want), abs_diff=d, tol=tol, ok=bool(ok)))
    if not ok:
        OUT["failures"].append(f"{label}: |d| {d:.3e}")
    P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.2e}")
    return ok


PHASE = sys.argv[1] if len(sys.argv) > 1 else "all"      # "filters" | "score" | "all"
FILT_NPZ = WORK / "stageC-filtered.npz"

t0 = time.time()
st = X.staged170()
tr = X.load_truth()
targets = st["h1_targets"]
n = targets.size
cols, cand = st["cols"], st["cand_h1"]
tset = set(targets.tolist())
inst = json.loads((E273 / "work/instrument273.json").read_text())
pm_variant = inst["donors"]["pm_blend_variant"]
P(f"286 Stage C -- {n} h1-test targets, {cols.size} operator columns; 273 PM blend variant {pm_variant!r}; "
  f"rss {rss():.2f} GiB")

# ------------------------------------------------------------------------------ the donor rows (273's construction)
rows273 = {ln: np.load(E273 / f"work/half/rows273-{ln}.npz", allow_pickle=True) for ln in X.DONORS}
real = {ln: np.asarray(rows273[ln]["real"], dtype=np.float64) for ln in X.DONORS}
avail = np.stack([np.asarray(rows273[ln]["avail"], bool) for ln in X.DONORS])
banks = {ln: E062.Bank.load(ln, plolo.DONOR_PATHS[ln]) for ln in X.DONORS}
pm_bank = E064.load_perturbmatch(plolo.PM_PATH)
no_panel = np.zeros(0, dtype=np.int64)
pm_prep = E062.prepare_main_source(pm_bank, list(targets), OFFICIAL, no_panel, chunk_rows=256, max_rss_gib=9.0)
shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank,
                                     OFFICIAL, PANEL_SET if pm_variant == "panel excluded" else set())
both = avail[0] & np.asarray(pm_prep.target_available, bool)
k_rows = real["k562_gwps"].copy()
r_ = np.flatnonzero(both)
k_rows[np.ix_(r_, shared)] = ((1 - plolo.PM_BLEND) * real["k562_gwps"][np.ix_(r_, shared)]
                              + plolo.PM_BLEND * np.asarray(pm_prep.effects)[np.ix_(r_, shared)])
donor_rows = {"k562_gwps": k_rows, "orion_hct116": real["orion_hct116"], "orion_hek293t": real["orion_hek293t"]}
del pm_bank, pm_prep, banks
gc.collect()

# ------------------------------------------------------------------------------ the filters (2 eigenproblems / line)
# Two phases in separate processes: the eigendecompositions raise this process's PEAK RSS above E062's
# `enforce_memory_limit` guard, so the filtered rows are written to disk and the scoring runs fresh.
filtered: dict[str, dict[str, np.ndarray]] = {ln: {} for ln in X.DONORS}
lk = {}
if PHASE == "score":
    z = np.load(FILT_NPZ, allow_pickle=True)
    for k in z.files:
        if "|" in k:
            ln, arm = k.split("|", 1)
            filtered[ln][arm] = np.asarray(z[k], dtype=np.float64)
    lk = json.loads(str(z["leakage_json"]))
    P(f"  loaded {FILT_NPZ.name}: {sum(len(v) for v in filtered.values())} filtered row blocks")
for ln in (X.DONORS if PHASE != "score" else ()):
    t1 = time.time()
    meta = json.loads((WORK / "moments" / f"{ln}.json").read_text())
    assert meta["est286_sha256"] == E.sha256_self()
    cj = np.load(WORK / "moments" / f"{ln}_cols.npy")
    Sj = np.load(WORK / "moments" / f"{ln}_S.npy")
    Nj = np.load(WORK / "moments" / f"{ln}_N.npy")
    lam, V, a_eig = E.psd_eigh(Sj, log=P)
    del Sj
    gc.collect()
    lam_top = lam[::-1][:100].copy()
    V_top = np.ascontiguousarray(V[:, ::-1][:, :100], dtype=np.float32)
    Sp = E.s_plus(lam, V)
    del V
    gc.collect()
    Lf, Q, k, a_w = E.wiener_factors(Sp, Nj, 1.0, log=P)
    del Sp, Nj
    gc.collect()
    Xj = donor_rows[ln][:, cj].copy()
    for phi in E.PHIS:
        full = donor_rows[ln].copy()
        full[:, cj] = E.apply_e1(Xj, Lf, Q, k, float(phi))
        filtered[ln][f"e1:phi{phi:g}"] = full
    full = donor_rows[ln].copy()
    full[:, cj] = E.apply_e1(Xj, Lf, Q, k, 0.10, mirror=True)
    filtered[ln]["e1:phi0.1:mirror"] = full
    for r in E.RANKS:
        full = donor_rows[ln].copy()
        full[:, cj] = E.apply_e2(Xj, Lf, lam_top, V_top, int(r))
        filtered[ln][f"e2:r{r}"] = full
    # the h1-target overlap with this line's FITTED bank rows (the leakage note)
    perts = set((ROOT / f"data/donors/{ln}/perts.txt").read_text().split("\n"))
    lk[ln] = dict(h1_targets_in_bank=int(sum(1 for t in targets if str(t) in perts)), n_fit=int(meta["n_rows"]),
                  share=float(sum(1 for t in targets if str(t) in perts) / meta["n_rows"]),
                  eig=a_eig, wiener=a_w)
    P(f"  [{ln}] filters built in {time.time()-t1:.0f}s; h1 targets inside the fitted bank rows "
      f"{lk[ln]['h1_targets_in_bank']}/{lk[ln]['n_fit']} ({100*lk[ln]['share']:.2f} %); rss {rss():.2f} GiB")
    del Lf, Q, k, lam, V_top, Xj
    gc.collect()
OUT["leakage"] = lk
ARMS = sorted(filtered[X.DONORS[0]])
if PHASE == "filters":
    np.savez_compressed(FILT_NPZ, leakage_json=json.dumps(lk, default=float),
                        **{f"{ln}|{arm}": v.astype(np.float32) for ln in X.DONORS for arm, v in filtered[ln].items()})
    P(f"[Stage C filters] wrote {FILT_NPZ.name}, {len(ARMS)} arms x 3 lines, "
      f"{(time.time()-t0)/60:.1f} min, peak rss {rss():.2f} GiB")
    raise SystemExit(0)

# ------------------------------------------------------------------------------ the operator + the CD4 part
op = plolo.Operator(plolo.PROD_OPS / "whitener_union_pruned.npz", name="union_pruned")
base = op.apply(cand, cols, GAMMA, LAM)
banks_h1 = [E062.Bank.load(nm, CD4_H1 / nm, require_cd4_scale=True) for nm in X.CD4_NAMES]
prep = E062.prepare_cd4_mean(banks_h1, list(targets), OFFICIAL, tset, "gene_norm", chunk_rows=256, max_rss_gib=40.0)
cd4h1 = X.complete_mean3(prep)
d5, aud5 = L.deflate(cd4h1, "svd", 5)
fr5, pu5 = X.solve_fraction(base, d5, CD4_F)
P(f"  base + H1 CD4 svd5 at f = 0.34: fraction {fr5:.6f} (per unit {pu5:.6f}); rss {rss():.2f} GiB")

arms: dict[str, dict] = {}
aud: dict[str, dict] = {}


def add(name: str, final: np.ndarray, **meta) -> dict:
    e = X.score_h1(final, tr)
    arms[name] = e
    aud[name] = dict(meta, median_row_norm=float(np.median(np.linalg.norm(final, axis=1))))
    P(f"  [{name:22s}] pds a23 {e['pds']:.6f}  own cos a23 {e['cos']:.5f}  cos_agg {e['a23_cos_agg']:.5f}  "
      f"delivered {e['a23_delivered_norm']:.4f}  ({time.time()-t0:.0f}s)")
    return e


add("BASE", base)
add("INCUMBENT", X.compose(base, [(d5, fr5)], cols), fraction=fr5, per_unit=pu5)
check("KA-C1 BASE pds == 273's instrument273 BASE", arms["BASE"]["pds"], inst["arms"]["BASE"]["pds"], 1e-12)
check("KA-C2 INCUMBENT pds == 273's instrument273 INCUMBENT", arms["INCUMBENT"]["pds"], inst["arms"]["INCUMBENT"]["pds"], 1e-12)


def build_pre(src_rows: dict[str, np.ndarray]) -> np.ndarray:
    src = []
    for j, ln in enumerate(X.DONORS):
        src.append(E062.PreparedSource(name=ln, effects=src_rows[ln], target_available=avail[j].copy(),
                                       mean_on_axis=np.asarray(rows273[ln]["mean_on_axis"]),
                                       raw_median_norm=float(rows273[ln]["raw_median_real"]),
                                       mapped_gene_count=0, duplicate_symbols=()))
    pre64, _ = E062.build_baseline(src)
    return pre64


pre64 = build_pre(donor_rows)
check("KA-C0 the rebuilt equal-weight mean row == 170's cand_h1 / alpha (273's KA-0)",
      float(np.abs(pre64 * ALPHA - cand).max() / np.abs(cand).max()), 0.0, 1e-8)


def donor_arm(name: str, pre_sel: np.ndarray, **meta):
    b = op.apply(pre_sel * ALPHA, cols, GAMMA, LAM)
    fr, pu = X.solve_fraction(b, d5, CD4_F)
    return add(name, X.compose(b, [(d5, fr)], cols), fraction=fr, per_unit=pu,
               rows_changed=int((np.abs(pre_sel - pre64).max(1) > 0).sum()),
               base_cos_to_incumbent_base=float(np.median(X.row_cos(b[:, cols], base[:, cols]))), **meta)


donor_arm("E0:rebuilt", pre64, note="the explicit equal-weight mean of the rebuilt rows (== INCUMBENT if KA-C0 holds)")
check("KA-C3 E0:rebuilt pds == INCUMBENT", arms["E0:rebuilt"]["pds"], arms["INCUMBENT"]["pds"], 1e-9)
for arm in ARMS:
    for tag, nm in ((":nm", True), (":lit", False)):
        sel = {ln: filtered[ln][arm] for ln in X.DONORS}
        pre = build_pre(sel)
        if nm:
            cov = avail.any(0)
            pre[cov] = E.norm_match(pre[cov], pre64[cov])
        donor_arm(arm + tag, pre, estimator=arm, norm_matched=nm)
        del pre, sel
        gc.collect()

# ------------------------------------------------------------------------------ the bootstrap vs BASE and INCUMBENT
OUT["arms"] = {k: X.summary(v) for k, v in arms.items()}
OUT["audits"] = aud
OUT["boot"] = {}
for name in arms:
    if name in ("BASE", "INCUMBENT"):
        continue
    rec = {}
    for ref in ("INCUMBENT", "BASE"):
        rec[f"pds vs {ref}"] = boot(arms[name]["pds_a23"] - arms[ref]["pds_a23"], rng)
        rec[f"cos vs {ref}"] = boot(arms[name]["cos_a23"] - arms[ref]["cos_a23"], rng)
    OUT["boot"][name] = rec
    P(f"  {name:22s} vs INCUMBENT: d pds {rec['pds vs INCUMBENT']['point']:+.5f} "
      f"[{rec['pds vs INCUMBENT']['lo']:+.5f}, {rec['pds vs INCUMBENT']['hi']:+.5f}]  d own cos "
      f"{rec['cos vs INCUMBENT']['point']:+.5f} [{rec['cos vs INCUMBENT']['lo']:+.5f}, {rec['cos vs INCUMBENT']['hi']:+.5f}]")

arrays = {f"{k}|{kk}": v for k, e in arms.items() for kk, v in e.items() if isinstance(v, np.ndarray)}
arrays["targets"] = targets
np.savez_compressed(WORK / "stageC286.npz", **arrays)
OUT["seconds"] = round(time.time() - t0, 1)
OUT["peak_rss_gib"] = rss()
(WORK / "stageC286.json").write_text(json.dumps(OUT, indent=1, default=float))
(WORK / "stageC286.log").write_text("\n".join(LOG) + "\n")
P(f"\n[Stage C] done {(time.time()-t0)/60:.1f} min, peak rss {rss():.2f} GiB, KA "
  f"{sum(1 for k in OUT['ka'] if k['ok'])}/{len(OUT['ka'])}, failures {OUT['failures']}")
