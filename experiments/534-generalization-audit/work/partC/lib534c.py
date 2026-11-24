"""534 Part C -- Part B's lib534.py (sha256 eec51cebaa3d0bd0eebc85e821064069694216e886c68935907197065f05886a) with WORK pointed at work/partC/ and three
base TOGGLES (PREREGISTRATION-C section 2), each read at call time and inert at its production value:
  NORMREL_MODE, NORMREL_EXP   the per-target donor weights of the base (lib418.pre_normrel; production cap, -1.5)
                              -> used by pre_base() and lodo_bases_base() below (the base AND the LODO key follow)
  PARTIAL_G                   the CD4 partial rows' sizing (production 0.5); compose() now passes it to compose_arm
                              on the abc-rules path too (lib418's default argument binds 0.5 at definition time)
Everything else is lib438 verbatim. Part B's docstring follows.

534 Part B -- a COPY of 438's lib438.py (sha256 44aec6bd6be31e5486327653ff776233d9cf1b2cec872f9e6c9d9b4615c798df) with WORK and every cache path pointed at
experiments/534-generalization-audit/work/partB/ (preregistration 534-B section 5). The only changes: this paragraph,
the write guard imported first (guard534: any write under the repo outside 534 raises), and HERE/ROOT/WORK below.
Everything else is lib438 verbatim, so the arithmetic is 438's. 438's docstring follows.

438 -- 420's bulk-stack pricer with the PANEL as an argument (a bundle's `pert_counts.csv`), not `data/challenge`.

Everything that depends on the panel follows it:
  * the PLOLO fold's target set          targets = panel n held-library coverage (`build_fold438`)
  * the operator axis                     union(the fold's source donors' genes) - panel, asserted against the operator's
                                          own axis record (plolo.Operator._check_axis) -- the Oct-22 trap fires here
  * donor coverage / normrel medians      `Donors438` (sens415.Donors's arithmetic, rows for ANY target list)
  * the CD4 / H1-union / partial residuals  the panel's own CD4 banks; H1 union with the panel's columns masked
  * the key                               410's formula and 415's LODO-V1 on the fold (and the full panel) of THIS panel
  * the frames                            b_ctrl / gate from the bundle's own control cells
and the composition RULES follow the bundle the way the October driver does:
  rules "abc"  (ABC_ENV)  CD4 all-rows sizing at f 0.34 x m, SVD-5 deflation re-fitted on the rows   == 418/420 exactly
  rules "def"  (DEF_ENV)  CD4 carrying-row sizing at g 0.3443655278649641 x m (`VCC_CD4_TARGET_G`), deflation by the
                          pinned A/B/C gene basis V5_abc (`VCC_CD4_DEFLATE_BASIS`), everything else identical.
420's arithmetic is kept: compose -> key_scale -> sharpen -> ballast realisation (KMAX 16) -> a23 read-back -> cos_agg,
pds, board-emulated pds (414's beta at lam 0.5), render A/B.

Nothing here writes outside experiments/438-october-stack-repricing/work (bytecode writing is disabled BEFORE any foreign
import; every foreign module that repoints plolo.WORK at import is followed by a reset to our own directory).
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import guard534  # noqa: E402,F401  (534: refuse writes outside experiments/534-generalization-audit)
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from pathlib import Path  # noqa: E402

os.environ.setdefault("OPENBLAS_NUM_THREADS", "6")
os.environ.setdefault("OMP_NUM_THREADS", "6")
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent                  # 534 Part C: experiments/534-generalization-audit/work/partC
ROOT = HERE.parents[3]
WORK = HERE                                            # 534 Part C: every lib cache (panels, bankstats, frames) lives here
WORK.mkdir(parents=True, exist_ok=True)
assert WORK == ROOT / "experiments/534-generalization-audit/work/partC", WORK
W420 = ROOT / "experiments/420-bulk-stack/work"
W418 = ROOT / "experiments/418-bulk-composition/work"
W414 = ROOT / "experiments/414-bulk-reprice/work"
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))

import lib420  # noqa: E402  (imports lib418 -> 162/218/228/410/169/415/414 machinery; each repoints plolo.WORK)
import lib418  # noqa: E402
import plolo  # noqa: E402
import arms414  # noqa: E402
import sens415  # noqa: E402
import resid218  # noqa: E402
from lib418 import (ALPHA, GAMMA, LAM, UNION, UNION_G, E062, Operator, compose_arm, key_scale, lodo_bases,  # noqa: E402
                    lodo_key, pre_normrel, sharpen, shared_energy_key, weights, donor_rows, source_residual, union_rows,
                    deflate, cd4_frozen_and_partial, solve_fraction, endpoints, deliver)
from lib420 import realise_k, KMAXES  # noqa: E402
from compose218 import solve_fraction_carrying  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

plolo.WORK = WORK / "plolo_unbound"      # nothing may be cached through plolo until bind() points it at a panel dir
E064 = plolo.E064
OFFICIAL = plolo.OFFICIAL
OFF_POS = {g: i for i, g in enumerate(OFFICIAL)}
MAIN = ("k562_gwps", "orion_hct116", "orion_hek293t")
CD4_F = 0.34
CD4_G_DEF = 0.3443655278649641          # DEF_ENV VCC_CD4_TARGET_G (== 1 x; 432 T6: x1.5 is 0.5165482917974462)
PARTIAL_G = 0.5                         # 534 Part C toggle (production 0.5; 0 = the partial rows are not added)
NORMREL_MODE, NORMREL_EXP = "cap", -1.5  # 534 Part C toggles (production cap, -1.5; lib418.pre_normrel's modes)
PRODUCTION = dict(NORMREL_MODE="cap", NORMREL_EXP=-1.5, PARTIAL_G=0.5)
V5_PATH = ROOT / "data/donors/cd4_deflation_basis/V5_abc.npz"
V5_SHA = "ba7ade3fc4f59de5c2de5d8bf270f243ab8e8ddfd0728be9e2c8c71c9b18122b"
V5_V_SHA = "dabf88f7f95ffc4dac109de640ad29b943773ec4e45c5ddb2a4231b0003cce61"
HELD_OPS = ROOT / "experiments/062-cd4-post-whiten/cloud-output/whiteners"
PROD_OP = ROOT / "data/donors/whitener_union_pruned.npz"
FOLD_SOURCES = {
    "k562_gwps": ("orion_hct116", "orion_hek293t"),
    "orion_hct116": ("k562_gwps", "orion_hek293t"),
    "orion_hek293t": ("k562_gwps", "orion_hct116"),
    "kolf2": MAIN,
}
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
SLOTS = ("A", "B", "C")                 # frame slots; a bundle's own context labels map onto them in order

R233 = ROOT / "experiments/233-october-retention/work"
R238 = ROOT / "experiments/238-essential-free-panel/work"
R432 = ROOT / "experiments/432-def-rehearsal-new-knobs/work"
R244 = ROOT / "experiments/244-october-rehearsal/work/run1"

# ---------------------------------------------------------------------------------------------- panel registry
# Each entry: the bundle; the composition rules the October driver would apply; the panel's CD4 banks; the held-out
# operators PLOLO's cancer folds need; the 3-donor operator of the full panel (the builder's); the measured artifacts
# whose render energy is a known answer for the cross-panel energy route (432's ntclevel output).
PANELS = {
    "abc": dict(bundle=ROOT / "data/challenge", rules="abc",
                cd4_root=ROOT / "data/external/cd4t/cloud-compact/cd4_compact",
                fold_ops={"k562_gwps": HELD_OPS / "held_k562_gwps.npz", "orion_hct116": HELD_OPS / "held_orion_hct116.npz",
                          "orion_hek293t": HELD_OPS / "held_orion_hek293t.npz", "kolf2": PROD_OP},
                full_op=PROD_OP, folds=("k562_gwps", "orion_hct116", "orion_hek293t", "kolf2"), synthetic=False,
                meta=ROOT / "experiments/410-pertscale-revisit/work/build_b/a-spreadsheet-achieved-enlightenment.h5ad.meta.json"),
    # A/B/C's panel under the D/E/F composition rules: separates "the rules changed" from "the panel changed"
    "abcdef": dict(bundle=ROOT / "data/challenge", rules="def",
                   cd4_root=ROOT / "data/external/cd4t/cloud-compact/cd4_compact",
                   fold_ops={"k562_gwps": HELD_OPS / "held_k562_gwps.npz", "orion_hct116": HELD_OPS / "held_orion_hct116.npz",
                             "orion_hek293t": HELD_OPS / "held_orion_hek293t.npz"},
                   full_op=PROD_OP, folds=CANCER, synthetic=False, meta=None),
    # 199's realistic draw (K562 n Orion n CD4, 34 essential genes): 432 built the stack here -> measured energies
    "real": dict(bundle=ROOT / "experiments/199-def-runbook/work/def-bundle-r", rules="def",
                 cd4_root=R233 / "cd4/real/cd4_compact",
                 fold_ops={h: R233 / f"ops/real__held_{h}.npz" for h in CANCER},
                 full_op=R432 / "base/whitener_union_pruned_def.npz", folds=CANCER, synthetic=True,
                 meta=R432 / "b0/b0.h5ad.meta.json", key_ref=R432 / "key/pert_scale_def432_lo0.3_hi3.npz",
                 measured={"b0": "b0/b0.ntclevel.json", "b1": "b1/b1.ntclevel.json"}),
    # 238's essential-free pessimistic draw (K562 146/300): 432's ef_b1 is the stack built here
    "efp": dict(bundle=R238 / "ef-bundle", rules="def", cd4_root=R238 / "cd4/efp/cd4_compact",
                fold_ops={h: R238 / f"ops/efp__held_{h}.npz" for h in CANCER},
                full_op=R244 / "whitener_union_pruned_def.npz", folds=CANCER, synthetic=True,
                meta=R432 / "ef_b1/ef_b1.h5ad.meta.json", measured={"ef_b1": "ef_b1/ef_b1.ntclevel.json"}),
    # 238's essential-free realistic draw (K562 300/300). No 3-donor operator exists for it (a fit peaks ~13 GB, over
    # this mandate's 10 GB cap): the full-panel route uses the same-axis 2-donor surrogate, validated on abc/real/efp.
    # mechanism probe: efp with the CD4 deflation re-fitted on the panel's own rows (250's rejected day-of re-fit)
    "efpsvd": dict(bundle=R238 / "ef-bundle", rules="defsvd", cd4_root=R238 / "cd4/efp/cd4_compact",
                   fold_ops={h: R238 / f"ops/efp__held_{h}.npz" for h in CANCER},
                   full_op=R244 / "whitener_union_pruned_def.npz", folds=CANCER, synthetic=True, meta=None),
    "efr": dict(bundle=R238 / "ef-bundle-r", rules="def", cd4_root=R238 / "cd4/efr/cd4_compact",
                fold_ops={h: R238 / f"ops/efr__held_{h}.npz" for h in CANCER},
                full_op=None, surrogate_ops=(R238 / "ops/efr__held_orion_hct116.npz", R238 / "ops/efr__held_orion_hek293t.npz"),
                folds=CANCER, synthetic=True, meta=None),
}
SURROGATES = {  # the same-axis held-out-Orion operators of each panel that HAS a real 3-donor operator
    "abc": (HELD_OPS / "held_orion_hct116.npz", HELD_OPS / "held_orion_hek293t.npz"),
    "real": (R233 / "ops/real__held_orion_hct116.npz", R233 / "ops/real__held_orion_hek293t.npz"),
    "efp": (R238 / "ops/efp__held_orion_hct116.npz", R238 / "ops/efp__held_orion_hek293t.npz"),
}


def sha256(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_csv_col(path: Path, col: str) -> tuple[str, ...]:
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip()]
    hdr = [c.strip() for c in lines[0].split(",")]
    j = hdr.index(col)
    return tuple(ln.split(",")[j].strip() for ln in lines[1:])


# ---------------------------------------------------------------------------------------------- the bound panel
@dataclass
class Ctx:
    key: str
    bundle: Path
    rules: str
    panel: tuple
    panel_set: set
    pp_pos: np.ndarray
    nonpanel: np.ndarray
    contexts: tuple
    cd4_root: Path
    fold_ops: dict
    full_op: Path | None
    folds: tuple
    work: Path
    spec: dict = field(default_factory=dict)


def bundle_frames(bundle: Path, contexts, work: Path):
    """b_ctrl / gate per context from the bundle's own control cells (plolo.challenge_bctrl's arithmetic)."""
    import h5py
    work.mkdir(parents=True, exist_ok=True)
    out = {}
    for slot, ctx in zip(SLOTS, contexts):
        cache = work / f"bctrl-{ctx}.npy"
        if cache.exists():
            b = np.load(cache)
        else:
            with h5py.File(bundle / f"context_{ctx}.h5ad", "r") as h:
                idx = h["X/indices"][:]
                dat = h["X/data"][:].astype(np.float64)
                node = h["var/_index"]
                raw = node["values"][:] if isinstance(node, h5py.Group) else node[:]
                genes = [g.decode() if isinstance(g, bytes) else str(g) for g in raw]
            assert tuple(genes) == OFFICIAL, f"{bundle}/context_{ctx}: gene axis differs from gene_names.csv"
            totals = np.bincount(idx, weights=dat, minlength=len(OFFICIAL))
            b = np.log1p(TS * totals / totals.sum())
            np.save(cache, b)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        out[slot] = (b, (cpm > plolo.MIN_CPM) & (b > 0))
    return out


def bind(key: str, spec: dict | None = None, tag: str | None = None) -> Ctx:
    """Point every panel-dependent global of the imported machinery at `key`'s panel. One panel per process."""
    if spec is None:
        if key in PANELS:
            spec = PANELS[key]
        else:                                   # a day-of panel: dayof438.py writes work/<key>/spec.json
            spec = json.loads((WORK / key / "spec.json").read_text())
            spec["fold_ops"] = {h: Path(v) for h, v in spec.get("fold_ops", {}).items()}
    spec = dict(spec)
    bundle = Path(spec["bundle"]).resolve()
    official = read_csv_col(bundle / "gene_names.csv", "gene_name")
    assert official == OFFICIAL, f"{bundle}: gene_names.csv differs from the official axis"
    panel = read_csv_col(bundle / "pert_counts.csv", "target_gene")
    assert len(panel) == 300 and len(set(panel)) == 300, f"{bundle}: {len(panel)} targets"
    man = json.loads((bundle / "manifest.json").read_text()) if (bundle / "manifest.json").exists() else {}
    contexts = tuple(man.get("contexts", ("A", "B", "C")))
    assert len(contexts) == 3
    work = WORK / (tag or key)
    work.mkdir(parents=True, exist_ok=True)
    pset = set(panel)
    pp_pos = np.asarray([OFF_POS[g] for g in panel if g in OFF_POS], dtype=np.int64)
    nonpanel = np.ones(len(OFFICIAL), bool)
    nonpanel[pp_pos] = False
    # --- rebind the module globals the imported machinery reads at CALL time
    plolo.PANEL, plolo.PANEL_SET = panel, pset
    plolo.WORK = work
    plolo.CD4_ROOT = Path(spec["cd4_root"])
    sens415.PANEL, sens415.PANEL_SET, sens415.CD4_ROOT = panel, pset, Path(spec["cd4_root"])
    resid218.PANEL_SET = pset
    resid218._PREP_CACHE.clear()
    # --- frames from the bundle's controls; a synthetic bundle's controls are A/B/C relabelled (asserted)
    fr = bundle_frames(bundle, contexts, work)
    if spec.get("synthetic") or tuple(contexts) == SLOTS:
        ref = bundle_frames(ROOT / "data/challenge", SLOTS, WORK / "abc_frames")
        for s in SLOTS:
            assert np.array_equal(fr[s][0], ref[s][0]) and np.array_equal(fr[s][1], ref[s][1]), \
                f"{key}: frame {s} != A/B/C's (a synthetic bundle must carry A/B/C's control cells)"
    arms414.FRAMES.clear()
    arms414.FRAMES.update(fr)
    return Ctx(key=key, bundle=bundle, rules=spec["rules"], panel=panel, panel_set=pset, pp_pos=pp_pos,
               nonpanel=nonpanel, contexts=contexts, cd4_root=Path(spec["cd4_root"]),
               fold_ops={k: Path(v) for k, v in spec.get("fold_ops", {}).items()},
               full_op=Path(spec["full_op"]) if spec.get("full_op") else None,
               folds=tuple(spec.get("folds", CANCER)), work=work, spec=spec)


# ---------------------------------------------------------------------------------------------- donor rows, any panel
BANKS: dict = {}


def bank(name: str):
    if name not in BANKS:
        BANKS[name] = (E064.load_perturbmatch(plolo.PM_PATH) if name == "perturbmatch"
                       else E062.Bank.load(name, plolo.DONOR_PATHS[name]))
    return BANKS[name]


def bank_stats(name: str):
    """Panel-INDEPENDENT bank statistics: prepare_main_source's column mean over all mapped columns (chunk 256) and the
    median centred row norm. Cached once for every panel (they depend on the bank alone)."""
    cache = WORK / "bankstats" / f"{name}.npz"
    b = bank(name)
    mp = E062.map_bank_to_axis(b, OFFICIAL)
    if cache.exists():
        z = np.load(cache)
        assert np.array_equal(z["all_source_columns"], mp.all_source_columns)
        return mp, z["mean_all"], float(z["raw_median"])
    mean_all = E062.stream_column_mean(b, mp.all_source_columns, chunk_rows=256, max_rss_gib=9.0)
    raw_median, _, _ = E062.centered_statistics(b, mp.all_source_columns, mean_all, need_gene_scale=False,
                                                chunk_rows=256, max_rss_gib=9.0)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez(cache, all_source_columns=mp.all_source_columns, mean_all=mean_all, raw_median=raw_median)
    return mp, mean_all, float(raw_median)


def rows(name: str, targets):
    """(effects with panel columns KEPT, available, mean_on_axis (unique columns), raw_median) for any target list.
    == sens415.full_rows / E062.prepare_main_source before its panel zeroing (asserted on A/B/C by the callers)."""
    b = bank(name)
    mp, mean_all, raw_median = bank_stats(name)
    at = {int(s): i for i, s in enumerate(mp.all_source_columns.tolist())}
    umean = np.asarray([mean_all[at[int(s)]] for s in mp.unique_source_columns], dtype=np.float64)
    mean_on_axis = np.zeros(len(OFFICIAL))
    mean_on_axis[mp.unique_target_positions] = umean
    delta = np.load(b.delta_path, mmap_mode="r")
    eff = np.zeros((len(targets), len(OFFICIAL)))
    avail = np.zeros(len(targets), bool)
    for i, t in enumerate(targets):
        r = b.pert_pos.get(str(t))
        if r is None:
            continue
        eff[i, mp.unique_target_positions] = (np.asarray(delta[r, mp.unique_source_columns], np.float64) - umean) \
            * (E062.EFFECT_SCALE / raw_median)
        avail[i] = True
    return eff, avail, mean_on_axis, raw_median, mp


def prepared(name: str, targets, pp_pos) -> "E062.PreparedSource":
    """E062.prepare_main_source(bank, targets, OFFICIAL, panel_positions) from the cached bank statistics."""
    eff, avail, moa, rmn, mp = rows(name, targets)
    eff[:, pp_pos] = 0.0
    moa = moa.copy()
    moa[pp_pos] = 0.0
    return E062.PreparedSource(name=name, effects=eff, target_available=avail, mean_on_axis=moa,
                               raw_median_norm=rmn, mapped_gene_count=len(mp.unique_target_positions),
                               duplicate_symbols=mp.duplicate_symbols)


class Donors438:
    """sens415.Donors for ANY target list: full rows (panel columns kept), availability, per-donor normrel medians over
    the targets each donor carries in THIS list, and the production PM50 correction (061's broad genes). `pre` is
    sens415.Donors.pre verbatim (normrel cap -1.5 via lib418.pre_normrel)."""

    pre = sens415.Donors.pre

    def __init__(self, targets, sources):
        self.sources = list(sources)
        self.F, self.A = {}, {}
        for s in sources:
            f, a, *_ = rows(s, targets)
            self.F[s], self.A[s] = f, a
        fpm, apm, *_ = rows("perturbmatch", targets)
        fk, ak = (self.F["k562_gwps"], self.A["k562_gwps"]) if "k562_gwps" in self.F else rows("k562_gwps", targets)[:2]
        both = apm & ak
        obs = np.zeros(len(OFFICIAL), bool)
        obs[E062.map_bank_to_axis(bank("perturbmatch"), OFFICIAL).unique_target_positions] = True
        corr = plolo.PM_BLEND * (fpm - fk)
        corr[:, ~obs] = 0.0
        corr[~both] = 0.0
        self.pmcorr = {"prod": corr}
        self.med = {s: float(np.median(np.linalg.norm(self.F[s][self.A[s]], axis=1))) for s in sources}


# ---------------------------------------------------------------------------------------------- PLOLO fold, any panel
def fold_operator(ctx: Ctx, held: str, fold_ops: str):
    """(operator path, its axis wcols). 'held' = the fold's held-out operator (PLOLO proper); 'prod' = the panel's 3-donor
    operator on its 3-donor axis (leaky: the held library is in its basis) -- the day-of shortcut, priced as a variant."""
    if fold_ops == "held":
        path = ctx.fold_ops[held]
        srcs = FOLD_SOURCES[held]
    elif fold_ops == "prod":
        path = ctx.full_op
        srcs = MAIN
    else:
        raise ValueError(fold_ops)
    wcols = E064.operator_columns([bank(n) for n in srcs], OFFICIAL, ctx.panel_set)
    return Path(path), wcols


def build_fold438(ctx: Ctx, held: str, fold_ops: str = "held") -> dict:
    """plolo.build_fold's arithmetic for ANY panel (targets = panel n held coverage), cached in this panel's work dir.
    pre/truth are returned float32-rounded, exactly as plolo's cache hit returns them (so a fresh build and a cached
    build are the same object)."""
    sources = FOLD_SOURCES[held]
    op_path, wcols = fold_operator(ctx, held, fold_ops)
    cache = ctx.work / f"fold-{held}-{fold_ops}.npz"
    held_bank = bank(held)
    targets = tuple(g for g in ctx.panel if g in held_bank.pert_pos)
    if cache.exists():
        z = np.load(cache, allow_pickle=True)
        assert tuple(str(t) for t in z["targets"]) == targets and np.array_equal(z["wcols"], wcols)
        out = {k: z[k] for k in ("pre", "truth", "b_ctrl", "gate", "scored", "wcols")}
        audit = json.loads(str(z["audit"]))
    else:
        prep = {n: prepared(n, targets, ctx.pp_pos) for n in sources}
        pm = prepared("perturbmatch", targets, ctx.pp_pos)
        src = []
        pm_n = 0
        for n in sources:
            if n == "k562_gwps":
                shared = E064.shared_k562_pm_columns(bank("k562_gwps"), bank("perturbmatch"), OFFICIAL, ctx.panel_set)
                eff = np.asarray(prep[n].effects, np.float64).copy()
                both = prep[n].target_available & pm.target_available
                rws = np.flatnonzero(both)
                pm_n = int(both.sum())
                if len(rws):
                    eff[np.ix_(rws, shared)] = ((1 - plolo.PM_BLEND) * prep[n].effects[np.ix_(rws, shared)]
                                                + plolo.PM_BLEND * pm.effects[np.ix_(rws, shared)])
                src.append(E062.PreparedSource(name="k562_gwps_pm50", effects=eff,
                                               target_available=np.asarray(prep[n].target_available, bool).copy(),
                                               mean_on_axis=np.asarray(prep[n].mean_on_axis, float).copy(),
                                               raw_median_norm=float(prep[n].raw_median_norm),
                                               mapped_gene_count=int(prep[n].mapped_gene_count),
                                               duplicate_symbols=tuple(prep[n].duplicate_symbols)))
            else:
                src.append(prep[n])
        pre, votes = E062.build_baseline(src)
        hp = prepared(held, targets, ctx.pp_pos)
        assert hp.target_available.all(), "a held target has no truth row"
        truth = E062.raw_truth(hp)
        b_ctrl, held_cols = plolo._b_ctrl_on_axis(held_bank)
        expr = np.expm1(b_ctrl)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate = (cpm > plolo.MIN_CPM) & (b_ctrl > 0)
        sm = np.zeros(len(OFFICIAL), bool)
        sm[held_cols] = True
        sm[ctx.pp_pos] = False
        out = dict(pre=pre.astype(np.float32), truth=truth.astype(np.float32), b_ctrl=b_ctrl, gate=gate,
                   scored=np.flatnonzero(sm), wcols=wcols)
        audit = dict(held=held, sources=list(sources), n_targets=len(targets),
                     panel_coverage=f"{len(targets)}/300", op=str(op_path.relative_to(ROOT)), fold_ops=fold_ops,
                     op_cols=int(len(wcols)), pm_supported_targets=pm_n,
                     carriers={str(v): int((votes == v).sum()) for v in np.unique(votes)},
                     source_coverage={n: int(prep[n].target_available.sum()) for n in sources})
        np.savez_compressed(cache, targets=np.asarray(targets), audit=json.dumps(audit), **out)
    op = plolo.Operator(op_path)
    op._check_axis(wcols)                                  # refuses an operator fitted for another panel / axis
    in_axis = sorted(set(targets) & {OFFICIAL[i] for i in wcols})
    assert not in_axis, f"{len(in_axis)} panel genes on the operator axis"
    return dict(held=held, targets=targets, pre=np.asarray(out["pre"], np.float64),
                truth=np.asarray(out["truth"], np.float64), b_ctrl=out["b_ctrl"], gate=out["gate"],
                scored=np.asarray(out["scored"]), wcols=np.asarray(out["wcols"]), sources=sources,
                panel_positions=ctx.pp_pos, op=op, op_path=op_path, audit=audit)


# ---------------------------------------------------------------------------------------------- composition, by rules
def load_v5():
    """The builder's pinned `V5_abc` basis (path + npz sha256, rule 11)."""
    assert sha256(V5_PATH) == V5_SHA, "V5_abc npz sha"
    V = np.asarray(np.load(V5_PATH)["V"], np.float64)
    assert V.shape == (len(OFFICIAL), 5) and np.abs(V.T @ V - np.eye(5)).max() < 1e-10, "V5_abc not orthonormal"
    return V


def deflate_pinned(cd4, V):
    """The builder's `deflate_cd4_residual_pinned` (250's lib250.deflate_gene_basis), verbatim algebra."""
    eff = np.asarray(cd4.mean_effect, np.float64).copy()
    sup = np.asarray(cd4.gene_votes) > 0
    eff[~sup] = 0.0
    rws = np.flatnonzero(np.asarray(cd4.target_condition_count) > 0)
    R = eff[rws]
    e0 = float((R * R).sum())
    Rc = R - R.mean(axis=0, keepdims=True)
    Rc[~sup[rws]] = 0.0
    ec = float((Rc * Rc).sum())
    proj = Rc @ V
    captured = float((proj * proj).sum() / max(ec, 1e-300))
    Rk = Rc - proj @ V.T
    Rk[~sup[rws]] = 0.0
    e1 = float((Rk * Rk).sum())
    eff[rws] = Rk
    audit = dict(mode="pinned", n_carried=int(rws.size), energy_removed_frac=float(1.0 - e1 / e0),
                 captured_frac_of_centred=captured, centred_energy_frac_of_raw=float(ec / e0))
    return E062.CD4Prepared(mode=cd4.mode, mean_effect=eff, gene_votes=np.asarray(cd4.gene_votes).copy(),
                            target_condition_count=np.asarray(cd4.target_condition_count).copy(),
                            condition_details=cd4.condition_details), audit


def residuals(ctx: Ctx, targets):
    """(deflated CD4 main residual, CD4 partial rows, H1 union, deflation audit) for `targets` under the panel's rules."""
    frozen, partial = cd4_frozen_and_partial(list(targets))
    if ctx.rules in ("abc", "defsvd"):          # defsvd: D/E/F sizing, but SVD-5 re-fitted on THIS panel's rows (a mechanism probe)
        cd4, aud = deflate(frozen, 5)
    else:
        cd4, aud = deflate_pinned(frozen, load_v5())
    h1 = union_rows([source_residual(n, list(targets), "gene_norm") for n in UNION])
    return cd4, partial, h1, aud


def compose(ctx: Ctx, whitened, cd4, partial, h1, wcols, pp_pos, m_cd4=1.0, m_h1=1.0):
    """lib418.compose_arm (rules abc: bit-identical call) or its D/E/F twin (CD4 sized on its CARRYING rows at
    g = 0.3443655 x m_cd4, `VCC_CD4_SIZING=carrying`). H1 union and partial are carrying-sized in both."""
    g_h1 = UNION_G * m_h1
    if ctx.rules == "abc":                      # defsvd falls through to the D/E/F (carrying) sizing below
        return compose_arm(whitened, cd4, partial, h1, wcols, pp_pos, m_cd4=m_cd4, g_h1=g_h1, g_part=PARTIAL_G)
    f_cd4, _ = solve_fraction(whitened, cd4, CD4_G_DEF * m_cd4, sizing="carrying")
    merged = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0) * f_cd4
    support = np.asarray(cd4.gene_votes) > 0
    fr = {"cd4": float(f_cd4), "h1_union": 0.0, "cd4_partial": 0.0}
    for tag, res, g in (("h1_union", h1, g_h1), ("cd4_partial", partial, PARTIAL_G)):
        if not g or not (np.asarray(res.gene_votes) > 0).any():
            continue
        f, _ = solve_fraction_carrying(whitened, res, g)
        v = np.asarray(res.gene_votes) > 0
        merged += np.where(v, np.asarray(res.mean_effect), 0.0) * f
        support |= v
        fr[tag] = float(f)
    added = E062.CD4Prepared(mode="merged", mean_effect=merged, gene_votes=support.astype(np.uint8),
                             target_condition_count=np.asarray(cd4.target_condition_count), condition_details=[])
    final, info = E062.add_post_whiten_residual(whitened, added, 1.0, wcols)
    final[:, pp_pos] = 0.0
    radd = np.where(support, merged, 0.0)
    return final, dict(radd=radd, s=float(info["post_add_global_preserve_scale"]), fractions=fr)


# ---------------------------------------------------------------------------------------------- the design grid
H1_GRID = (0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)
CD4_GRID = (0.75, 1.0, 1.25, 1.5, 1.75, 2.0)


def cname(mh, mc):
    return f"h{mh:g}c{mc:g}"


COMPS = {cname(mh, mc): (mc, mh) for mh in H1_GRID for mc in CD4_GRID}     # name -> (m_cd4, m_h1)
ALIAS = {"PP": "h1c1", "h1x2": "h2c1", "cd4x1.5": "h1c1.5", "stack": "h2c1.5"}
CROSS = [cname(mh, 1.5) for mh in H1_GRID] + [cname(2.0, mc) for mc in CD4_GRID] + ["h1c1", "h2c1"]
CROSS = list(dict.fromkeys(CROSS))
KEYS = {"K3": ("K410", 0.3, 3.0), "K5_3": ("K410", 0.5, 3.0), "V1": ("V1", 0.3, 3.0), "V1_5_3": ("V1", 0.5, 3.0)}


def grid():
    """(name, comp, key, q) for every design priced on a fold / the full panel. `name` = comp|key[+shq]."""
    out = [("h1c1|none", "h1c1", None, None), ("h1c1|K3", "h1c1", "K3", None)]
    for c in COMPS:
        for q in (1.25, 1.5):
            out.append((f"{c}|K3+sh{q:g}", c, "K3", q))
    for c in CROSS:
        for q in (1.25, 1.5):
            out.append((f"{c}|V1_5_3+sh{q:g}", c, "V1_5_3", q))
    for c in ("h1c1", "h2c1.5"):
        for key in ("V1", "K5_3"):
            for q in (1.25, 1.5):
                out.append((f"{c}|{key}+sh{q:g}", c, key, q))
        for q in (1.35, 1.75, 2.0):
            out.append((f"{c}|K3+sh{q:g}", c, "K3", q))
    seen, uniq = set(), []
    for r in out:
        if r[0] not in seen:
            seen.add(r[0])
            uniq.append(r)
    return uniq


def parse(name):
    comp, rest = name.split("|")
    key, q = (rest.split("+sh") + [None])[:2]
    return comp, (None if key == "none" else key), (float(q) if q else None)


# ---------------------------------------------------------------------------------------------- read-back
def row_score(dls, scored, truth, betas, nonpanel):
    """comp420.score_rows (identical arithmetic): per-frame endpoints, cos_agg, board-emulated pds lam 0.5/1, render A/B."""
    out, X, Aa, pdsp, pdsl, AB = {}, [], [], [], [], {}
    for ctx in SLOTS:
        dl = dls[ctx]
        P = dl[:, scored]
        e = endpoints(P, truth)
        e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
        for lam in (0.5, 1.0):
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
    summ = {k: float(np.mean([out[c][k] for c in SLOTS])) for k in keys}
    summ["per_frame"] = {c: {k: out[c][k] for k in ("pds", "cos_agg", "pds_lam0.5")} for c in SLOTS}
    summ["AB"] = AB
    arrays = dict(x=np.asarray(X), a=np.asarray(Aa), pds=np.asarray(pdsp), pdsl=np.asarray(pdsl))
    return summ, arrays


def cmp(s, r, fields=("cos_agg", "pds", "pds_lam0.5", "pds_lam1")):
    d = max(abs(s[k] - r[k]) for k in fields)
    return max(d, *[abs(s["AB"][c][k] / r["AB"][c][k] - 1) for c in SLOTS for k in ("A", "B")])


def carrier_profiles(carrier_design):
    """lib420.de_profiles on the bound frames: the carriers' design profile (PP unkeyed, unsharpened) per frame."""
    return {c: [profile_from_effect(carrier_design[i], b, g, "minnorm") for i in range(carrier_design.shape[0])]
            for c, (b, g) in arms414.frames().items()}


def realised_rows(bulk, P_DE, q=None, scale=None, kmaxes=(16,), zero_cols=None, want_design=True):
    """lib420.design_and_realised (identical arithmetic; kmaxes/design optional to save time)."""
    out = {}
    for c, (b, gate) in arms414.frames().items():
        F = bulk if q is None else sharpen(bulk, q, gate, mode="row")
        if scale is not None:
            F = F * float(scale[c])
        n = F.shape[0]
        des = np.empty_like(F) if want_design else None
        rea = {k: np.empty_like(F) for k in kmaxes}
        st = {k: dict(K=np.empty(n, int), r=np.empty(n), clip_mass=np.empty(n), clip_genes=np.empty(n, int)) for k in kmaxes}
        for i in range(n):
            p_bk = profile_from_effect(F[i], b, gate, "minnorm")
            if want_design:
                des[i] = np.log1p(TS * p_bk) - b
            for k in kmaxes:
                p_re, r, K, cm, cg = realise_k(P_DE[c][i], p_bk, k)
                rea[k][i] = np.log1p(TS * p_re) - b
                st[k]["K"][i], st[k]["r"][i], st[k]["clip_mass"][i], st[k]["clip_genes"][i] = K, r, cm, cg
        if zero_cols is not None:
            if want_design:
                des[:, zero_cols] = 0.0
            for k in kmaxes:
                rea[k][:, zero_cols] = 0.0
        out[c] = dict(design=des, realised=rea, stats=st)
    return out


def energy(H, nonpanel):
    H = H[:, nonpanel]
    Hc = H - H.mean(0, keepdims=True)
    return dict(A=float((H * H).sum(1).mean()), B=float((Hc * Hc).sum(1).mean()))


# ---------------------------------------------------------------------------------------------- board-emulated beta
def betas(ctx: Ctx, held: str, nonpanel, A_pp_k3=None):
    """414's artifact beta per frame at lam 0.5. On A/B/C: lib418.betas_for verbatim. Elsewhere: the same A/B/C artifact
    pattern (the emitter + the same control cells), zeroed on THIS panel, at the SAME energy relative to the fold's
    PP|K3 design energy as on the A/B/C fold (so beta keeps its A/B/C share of the row energy)."""
    if ctx.key == "abc":
        return lib418.betas_for(held, nonpanel)
    ref = json.loads((W414 / f"arms414-{held}.json").read_text())["arms"]
    abc = json.loads((W420 / f"comp420-{held}.json").read_text())["arms"]["PP|K3"]["design"]["AB"]
    bz = np.load(W414 / "beta414.npz")
    dec = json.loads((W414 / "decomp414.json").read_text())
    out = {}
    for c in SLOTS:
        bb = np.asarray(bz[f"beta_{c}"], float).copy()
        bb[~nonpanel] = 0.0
        s_abc = np.sqrt(ref["base"]["AB"][c]["A"] / (dec["395"][c]["A"] * 18532 / 300))
        out[c] = bb * s_abc * np.sqrt(A_pp_k3[c] / abc[c]["A"])
    return out


# ---------------------------------------------------------------------------------------------- claims / level rule
def ntc_claims(ctx: Ctx):
    """ntclevel417's per-context claim = f * Cp46 * sum(w) with f = 1 (an unfitted context), from the bundle's controls.
    sum(w) = cell_eval2's `_row_weights` over the panel's 300 rows (identical for every artifact of the panel)."""
    cache = ctx.work / "claims.json"
    if cache.exists():
        return json.loads(cache.read_text())
    sys.path.insert(0, str(ROOT / "experiments/417-bulk-sharpen"))
    import ntclevel417
    from cell_eval2.metrics.delta import _row_weights
    cols = {j: OFF_POS[g] for j, g in enumerate(ctx.panel) if g in OFF_POS}
    wsum = float(_row_weights(len(ctx.panel), len(OFFICIAL), cols).sum())
    out = {}
    for slot, c in zip(SLOTS, ctx.contexts):
        m, sd, n = ntclevel417.cp46(ctx.bundle, c)
        out[slot] = dict(context=c, Cp46=m, Cp46_sd=sd, n_constructs=n, wsum=wsum, claim=m * wsum)
    cache.write_text(json.dumps(out, indent=1))
    return out


def level_rho(sd: float) -> float:
    """ntclevel417's EV-optimal margin: s = s* exp(sd Phi^-1(q)), q = 0.31/(0.31+1.78); B scales as s^2."""
    from scipy.stats import norm
    q = 0.31 / (0.31 + 1.78)
    return float(np.exp(2.0 * sd * norm.ppf(q)))


class Finals:
    """Composition rows on demand (42 compositions x (final + merged residual) would hold ~3.4 GB). Pinned compositions
    stay resident; the rest live in a small LRU and are recomposed (deterministically, ~0.1 s) when evicted."""

    def __init__(self, ctx, whit, cd4, partial, h1, wcols, pp_pos, pinned=("h1c1", "h2c1.5"), size=2):
        self.args = (ctx, whit, cd4, partial, h1, wcols, pp_pos)
        self.pinned, self.size = set(pinned), size
        self.cache, self.order, self.meta = {}, [], {}

    def _make(self, c):
        mc, mh = COMPS[c]
        ctx, whit, cd4, partial, h1, wcols, pp_pos = self.args
        f, p = compose(ctx, whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, m_h1=mh)
        self.meta[c] = dict(m_cd4=mc, m_h1=mh, fractions=p["fractions"], s=p["s"])
        return f, p

    def get(self, c):
        if c in self.cache:
            if c in self.order:
                self.order.remove(c)
                self.order.append(c)
            return self.cache[c]
        fp = self._make(c)
        self.cache[c] = fp
        if c not in self.pinned:
            self.order.append(c)
            while len(self.order) > self.size:
                del self.cache[self.order.pop(0)]
        return fp

    def final(self, c):
        return self.get(c)[0]

    def parts(self, c):
        return self.get(c)[1]


# ---------------------------------------------------------------------------------------------- 534 Part C: the toggles
def set_toggles(NORMREL_MODE=None, NORMREL_EXP=None, PARTIAL_G=None, **_):
    """Set the three base toggles (module globals read at call time). Pass the full dict of a design (designs534c)."""
    g = globals()
    g["NORMREL_MODE"], g["NORMREL_EXP"], g["PARTIAL_G"] = NORMREL_MODE, NORMREL_EXP, PARTIAL_G


def toggles() -> dict:
    return dict(NORMREL_MODE=NORMREL_MODE, NORMREL_EXP=NORMREL_EXP, PARTIAL_G=PARTIAL_G)


def pre_base(dn, sources, fallback_rows=None):
    """The base donor mean under the current normrel toggles: lib418.pre_normrel(dn, sources, mode, exp) -- at production
    (cap, -1.5) exactly comp438/full438's call."""
    return pre_normrel(dn, sources, NORMREL_MODE, NORMREL_EXP, fallback_rows=fallback_rows)


def lodo_bases_base(op, dn, sources, pre, wcols):
    """415's LODO bases under the current normrel toggles (the key follows the base, as 418 priced its normrel arms)."""
    return lodo_bases(op, dn, sources, NORMREL_MODE, NORMREL_EXP, pre, wcols)
