"""218 Phase 0b — CONTRIBUTION: every candidate source measured the way CD4 is measured.

The form is fixed by 207: an **un-whitened post-whitening residual** at a solved dose `f`,
scored on 162's PLOLO instrument, primary endpoint the energy-weighted aggregate cosine
`cos_agg` (204: instrument resolution +-0.20% relative), `pds` secondary, `q_w` third.
Nothing about the donor tables, the operators, the PerturbMatch blend, alpha/gamma/lam, the
renderer or the read-back is touched.  The ONLY thing that varies between arms is WHICH
SOURCE supplies the residual (and, for one arm family, how it is normalised).

WHY THIS EXISTS.  `cd4f034` moved the board +0.0131 by raising one source's dose.  The
mandate's hypothesis is that the property which makes CD4 pay is INDEPENDENCE from
Replogle/Orion.  That is testable on data already on disk: put every other source we hold
through the identical construction and see whether the ranking follows independence.

KNOWN-ANSWER ROWS (AGENTS.md 11.1 rule 10) -- the run aborts if any fails:
  KA-1  `cd4_gn` at f=0.10 must reproduce 207's `f0.1` cos_agg/pds per fold to 1e-9.
  KA-2  `cd4_gn` at f=0 must reproduce 207's `f0`.
  KA-3  `kolfmask` must reproduce 207's `kolf0.05/0.1/0.2/0.4` (207's `alt_line_residual`
        construction: the alt line's effects behind CD4's OWN gene-support mask).
  KA-4  my multi-bank residual builder, run on the three CD4 banks in `gene_norm` mode,
        must equal production `prepare_cd4_mean`'s object bit-for-bit.
  KA-5  NON-VACUITY: every arm reports its median residual row norm and its carrying-row
        count.  A source that covers <=50% of a fold's targets has a zero MEDIAN and the
        sizing rule silently returns fraction 0 (200 found it, 207 verified it).  Those
        arms are reported as VACUOUS, not as "no effect".

    uv run python resid218.py --arms all
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, OFFICIAL, PANEL_SET, PROD_OPS, Operator,  # noqa: E402
                   build_fold, cd4_for, deliver, endpoints, solve_fraction)

sys.path.insert(0, str(HERE))
from sources218 import SOURCES  # noqa: E402

ALPHA = 0.982
GAMMA = 0.65
LAM = 0.0
CD4_NAMES = ("cd4_rest", "cd4_stim8hr", "cd4_stim48hr")
CD4_ROOT = ROOT / "data/external/cd4t/cloud-compact/cd4_compact"

DOSE = [0.05, 0.10, 0.20, 0.34, 0.55, 0.80]

# Sources with enough fold coverage for the median sizing rule, plus two deliberately
# below it so the vacuity assertion has something to fire on.
CANDIDATES = ("kolf2", "feng_ipsc", "k562_gwps", "orion_hct116", "orion_hek293t",
              "perturbmatch", "telohaec", "gse295214_h1")


def sha8(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(8 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ------------------------------------------------------------------ residual builders
_BANK_CACHE: dict[str, object] = {}
_PREP_CACHE: dict[tuple, object] = {}


def bank_for(name: str):
    """`Bank.load`, but tolerant of a table with no `b_ctrl.npy`.

    `perturbmatch` is a reprocessing of the same K562 cells `k562_gwps` holds and ships
    only `delta.npy`.  Nothing in the residual construction reads `b_ctrl` (it centres on
    the bank's own row mean), so the Bank is assembled directly rather than skipping the
    single most USEFUL low-independence anchor the regression has.
    """
    if name in _BANK_CACHE:
        return _BANK_CACHE[name]
    path = (ROOT / SOURCES[name]["path"]).resolve()
    if (path / "b_ctrl.npy").exists():
        _BANK_CACHE[name] = E062.Bank.load(name, path)
        return _BANK_CACHE[name]
    genes = E062.read_lines(path / "genes.txt", require_unique=False)
    perts = E062.read_lines(path / "perts.txt")
    delta = np.load(path / "delta.npy", mmap_mode="r")
    assert delta.shape == (len(perts), len(genes)), f"{name}: shape mismatch"
    _BANK_CACHE[name] = E062.Bank(
        name=name, path=path, genes=genes, perts=perts,
        pert_pos={p: i for i, p in enumerate(perts)},
        delta_path=path / "delta.npy", b_ctrl_path=path / "b_ctrl.npy",
        shape=tuple(map(int, delta.shape)))
    return _BANK_CACHE[name]


def multi_source_residual(banks, targets, mode: str):
    """`prepare_cd4_mean` generalised to any bank list, duplicate-symbol safe.

    Construction-identical to production's CD4 path (centre on the bank's own row mean,
    optional per-gene scale, rescale the median centred row norm to EFFECT_SCALE, average
    across banks per gene by vote count), except that the official-axis mapping goes through
    `map_bank_to_axis`, which keeps the LAST source column for a duplicated symbol instead of
    raising.  On banks with unique symbols the two paths are identical -- KA-4 asserts it.
    """
    n_t, n_g = len(targets), len(OFFICIAL)
    value_sum = np.zeros((n_t, n_g), dtype=np.float64)
    votes = np.zeros((n_t, n_g), dtype=np.uint8)
    counts = np.zeros(n_t, dtype=np.uint8)
    details = []
    panel_pos = {g: i for i, g in enumerate(OFFICIAL)}
    for bank in banks:
        amap = E062.map_bank_to_axis(bank, OFFICIAL)
        keep = np.asarray([OFFICIAL[p] not in PANEL_SET
                           for p in amap.unique_target_positions.tolist()], dtype=bool)
        cols = amap.unique_source_columns[keep]
        pos = amap.unique_target_positions[keep]
        mean = E062.stream_column_mean(bank, cols, chunk_rows=256, max_rss_gib=10.0)
        raw_median, gene_scale, norm_median = E062.centered_statistics(
            bank, cols, mean, need_gene_scale=(mode == "gene_norm"),
            chunk_rows=256, max_rss_gib=10.0)
        mult = (E062.EFFECT_SCALE / norm_median) if mode == "gene_norm" \
            else (E062.EFFECT_SCALE / raw_median)
        delta = np.load(bank.delta_path, mmap_mode="r")
        covered = 0
        for i, t in enumerate(targets):
            row = bank.pert_pos.get(t)
            if row is None:
                continue
            centered = np.asarray(delta[row, cols], dtype=np.float64) - mean
            if mode == "gene_norm":
                centered = centered / gene_scale
            value_sum[i, pos] += centered * mult
            votes[i, pos] += 1
            counts[i] += 1
            covered += 1
        details.append(dict(condition=bank.name, mapped_non_panel_genes=int(len(cols)),
                            covered_targets=covered, raw_median_centered_norm=raw_median,
                            gene_normalized_median_before_rescale=norm_median))
    avail = votes > 0
    mean_effect = np.zeros_like(value_sum)
    mean_effect[avail] = value_sum[avail] / votes[avail]
    complete = counts == len(banks)
    mean_effect[~complete] = 0.0
    votes[~complete] = 0
    counts[~complete] = 0
    assert panel_pos is not None
    return E062.CD4Prepared(mode=f"{mode}:{'+'.join(b.name for b in banks)}",
                            mean_effect=mean_effect, gene_votes=votes,
                            target_condition_count=counts, condition_details=details)


def source_residual(name: str, targets, mode: str):
    key = (name, mode, tuple(targets))
    if key not in _PREP_CACHE:
        _PREP_CACHE[key] = multi_source_residual([bank_for(name)], targets, mode)
    return _PREP_CACHE[key]


def kolf_cd4mask(cd4, targets):
    """207's `alt_line_residual` verbatim: the alt line behind CD4's OWN gene mask."""
    key = ("kolf2", "cd4mask", tuple(targets))
    if key not in _PREP_CACHE:
        bank = bank_for("kolf2")
        official_pos = {g: i for i, g in enumerate(OFFICIAL)}
        panel_positions = np.asarray([official_pos[g] for g in PANEL_SET
                                      if g in official_pos], dtype=np.int64)
        prep = E062.prepare_main_source(bank, targets, OFFICIAL, panel_positions,
                                        chunk_rows=256, max_rss_gib=10.0)
        _PREP_CACHE[key] = (np.asarray(prep.effects, dtype=np.float64),
                            np.asarray(prep.target_available, bool))
    eff, avail = _PREP_CACHE[key]
    votes = np.asarray(cd4.gene_votes).copy()
    votes[~avail] = 0
    counts = np.asarray(cd4.target_condition_count).copy()
    counts[~avail] = 0
    return E062.CD4Prepared(mode="alt_kolf2", mean_effect=np.where(votes > 0, eff, 0.0),
                            gene_votes=votes, target_condition_count=counts,
                            condition_details=[{"condition": "kolf2"}])


def shuffled(res, seed: int = 20726):
    rng = np.random.default_rng(seed)
    n = res.mean_effect.shape[0]
    carry = np.flatnonzero(np.asarray(res.target_condition_count) > 0)
    perm = np.arange(n)
    sh = carry.copy()
    while True:
        rng.shuffle(sh)
        if len(carry) < 2 or not np.any(sh == carry):
            break
    perm[carry] = sh
    return E062.CD4Prepared(mode=res.mode + "+shuffled",
                            mean_effect=np.asarray(res.mean_effect)[perm].copy(),
                            gene_votes=np.asarray(res.gene_votes)[perm].copy(),
                            target_condition_count=np.asarray(
                                res.target_condition_count)[perm].copy(),
                            condition_details=res.condition_details)


def self_residual(cd4, pre):
    eff = np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(pre, dtype=np.float64), 0.0)
    return E062.CD4Prepared(mode="self_unwhitened", mean_effect=eff,
                            gene_votes=np.asarray(cd4.gene_votes).copy(),
                            target_condition_count=np.asarray(
                                cd4.target_condition_count).copy(),
                            condition_details=[{"condition": "base_mixture_unwhitened"}])


def oracle_residual(cd4, truth_full):
    """UPPER BOUND, never shippable: the held line's own truth as the residual.

    This is the ceiling of the un-whitened-residual channel -- what a source that carried
    the recipient's exact response would buy.  Every real source is a fraction of it.
    """
    eff = np.where(np.asarray(cd4.gene_votes) > 0,
                   np.asarray(truth_full, dtype=np.float64), 0.0)
    return E062.CD4Prepared(mode="ORACLE_truth", mean_effect=eff,
                            gene_votes=np.asarray(cd4.gene_votes).copy(),
                            target_condition_count=np.asarray(
                                cd4.target_condition_count).copy(),
                            condition_details=[{"condition": "held_truth"}])


# ------------------------------------------------------------------------------ score
def _per_target(P, T):
    A = P / np.maximum(np.linalg.norm(P, axis=1, keepdims=True), 1e-300)
    B = T / np.maximum(np.linalg.norm(T, axis=1, keepdims=True), 1e-300)
    S = A @ B.T
    n = S.shape[0]
    own = np.diag(S).copy()
    pds = ((S < own[:, None]).sum(1) + 0.5 * ((S == own[:, None]).sum(1) - 1)) / (n - 1)
    return pds, own


def score_effects(final, fold, truth) -> dict:
    out: dict[str, object] = {}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        acc, per, perc = {}, [], []
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            dl = deliver(final, b, gate, mass)
            e = endpoints(dl[:, fold["scored"]], truth)
            e["cos_agg"] = float(np.sqrt(max(0.0, 1.0 - e["sse"])))
            acc[ctx] = e
            p, c = _per_target(dl[:, fold["scored"]], truth)
            per.append(p)
            perc.append(c)
        for k in acc["A"]:
            out[f"{tag}_{k}"] = float(np.mean([acc[c][k] for c in "ABC"]))
        out[f"_pds_per_target_{tag}"] = np.mean(per, axis=0).tolist()
        out[f"_cos_per_target_{tag}"] = np.mean(perc, axis=0).tolist()
    out.update({k: out[f"a23_{k}"] for k in ("pds", "cos", "cos_agg", "sgn", "sgn_w", "sse")})
    out["median_row_norm_final"] = float(np.median(np.linalg.norm(final, axis=1)))
    return out


def build_arms() -> dict[str, dict]:
    arms: dict[str, dict] = {"f0": dict(kind="none", note="baseline: no residual (CD4 OFF)")}
    for f in DOSE:
        arms[f"cd4_gn{f:g}"] = dict(kind="cd4", mode="gene_norm", f=f,
                                    note=f"SHIPPED construction, f={f:g} (KA at 0.10)")
        arms[f"cd4_raw{f:g}"] = dict(kind="cd4_multi", mode="raw", f=f,
                                     note=f"CD4, raw (no per-gene scale), f={f:g}")
        arms[f"self{f:g}"] = dict(kind="self", f=f, note=f"099 self-residual, f={f:g}")
        arms[f"shuf{f:g}"] = dict(kind="shuf", f=f, note=f"target-shuffled CD4, f={f:g}")
        arms[f"oracle{f:g}"] = dict(kind="oracle", f=f,
                                    note=f"UPPER BOUND: held truth as residual, f={f:g}")
        for src in CANDIDATES:
            for mode in ("gene_norm", "raw"):
                tag = "gn" if mode == "gene_norm" else "rw"
                arms[f"{src}_{tag}{f:g}"] = dict(kind="src", src=src, mode=mode, f=f,
                                                 note=f"{src} as un-whitened residual "
                                                      f"({mode}), f={f:g}")
    for f in (0.05, 0.10, 0.20, 0.40):
        arms[f"kolfmask{f:g}"] = dict(kind="kolfmask", f=f,
                                      note=f"207 alt_line_residual KA row, f={f:g}")
    return arms


def residual_for(spec, cd4, cd4_banks, pre, truth_full, targets):
    k = spec["kind"]
    if k == "cd4":
        return cd4
    if k == "cd4_multi":
        return multi_source_residual(cd4_banks, targets, spec["mode"])
    if k == "self":
        return self_residual(cd4, pre)
    if k == "shuf":
        return shuffled(cd4)
    if k == "oracle":
        return oracle_residual(cd4, truth_full)
    if k == "kolfmask":
        return kolf_cd4mask(cd4, targets)
    if k == "src":
        return source_residual(spec["src"], targets, spec["mode"])
    raise ValueError(k)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", nargs="*", default=["all"])
    ap.add_argument("--folds", nargs="*", default=list(FOLDS))
    ap.add_argument("--out", default=str(WORK / "resid218.json"))
    args = ap.parse_args()

    arms = build_arms()
    names = list(arms) if args.arms == ["all"] else args.arms
    out_path = Path(args.out)
    rows = json.loads(out_path.read_text()) if out_path.exists() else []
    done = {(r["fold"], r["arm"]) for r in rows}

    ka207 = {}
    for r in json.loads((ROOT / "experiments/207-cd4-dose/work/sweep207.json").read_text()):
        ka207[(r["fold"], r["arm"])] = r

    for held in args.folds:
        fold = build_fold(held)
        targets = fold["targets"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, fold["scored"]]
        truth_full = np.asarray(fold["truth"], dtype=np.float64)
        pre = np.asarray(fold["pre"], dtype=np.float64)
        cd4 = cd4_for(targets)
        cd4_banks = [E062.Bank.load(n, CD4_ROOT / n, require_cd4_scale=True)
                     for n in CD4_NAMES]

        # KA-4: my builder must reproduce production's CD4 object exactly.
        mine = multi_source_residual(cd4_banks, targets, "gene_norm")
        d_eff = float(np.abs(mine.mean_effect - cd4.mean_effect).max())
        d_votes = int(np.abs(mine.gene_votes.astype(int)
                             - np.asarray(cd4.gene_votes, dtype=int)).max())
        assert d_eff < 1e-12 and d_votes == 0, \
            f"KA-4 FAILED on {held}: max|d effect| {d_eff:.3e}, votes {d_votes}"
        print(f"[{held}] KA-4 PASS (own builder == prepare_cd4_mean, "
              f"max abs diff {d_eff:.2e})", flush=True)

        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        op_meta = dict(path=str(op_path.relative_to(ROOT)), sha8=sha8(op_path)[:8],
                       n_wcols=int(len(fold["wcols"])))
        whitened = op.apply(pre * ALPHA, fold["wcols"], GAMMA, LAM)

        for arm in names:
            if (held, arm) in done:
                continue
            spec = arms[arm]
            if spec.get("src") == held:
                continue                                    # never leak the truth
            t0 = time.time()
            if spec["kind"] == "none":
                res, frac, per_unit = None, 0.0, None
                med_row = 0.0
                carrying = 0
            else:
                res = residual_for(spec, cd4, cd4_banks, pre, truth_full, targets)
                added = np.where(np.asarray(res.gene_votes) > 0,
                                 np.asarray(res.mean_effect), 0.0)
                med_row = float(np.median(np.linalg.norm(added, axis=1)))
                carrying = int(np.sum(np.asarray(res.target_condition_count) > 0))
                frac, per_unit = solve_fraction(whitened, res, spec["f"])
            if frac > 0.0:
                final, preserve = E062.add_post_whiten_residual(
                    whitened, res, frac, fold["wcols"])
                pscale = float(preserve["post_add_global_preserve_scale"])
            else:
                final, pscale = whitened.copy(), 1.0
            final[:, fold["panel_positions"]] = 0.0
            e = score_effects(final, fold, truth)
            vacuous = bool(spec["kind"] != "none" and med_row <= 0.0)
            aud = dict(kind=spec["kind"], src=spec.get("src"), mode=spec.get("mode"),
                       f=spec.get("f"), fraction=float(frac),
                       f_per_unit=(float(per_unit) if per_unit else None),
                       median_residual_row_norm=med_row, carrying_rows=carrying,
                       coverage_frac=(carrying / len(targets)),
                       VACUOUS=vacuous, preserve_scale=pscale, operator=op_meta)
            rows.append(dict(fold=held, arm=arm, n_targets=len(targets),
                             note=spec["note"], **e, audit=aud,
                             seconds=round(time.time() - t0, 1)))
            out_path.write_text(json.dumps(rows, indent=2))
            flag = " [VACUOUS]" if vacuous else ""
            print(f"[{held}/{arm}] frac {frac:.5f} carry {carrying}/{len(targets)} "
                  f"cos_agg {e['cos_agg']:.6f} pds {e['pds']:.6f} "
                  f"q_w {e['a15_sgn_w']:.6f}{flag} ({time.time()-t0:.0f}s)", flush=True)

            # KA-1/2/3 as they come up
            ka_map = {"cd4_gn0.05": "f0.05", "cd4_gn0.1": "f0.1", "cd4_gn0.2": "f0.2",
                      "cd4_gn0.34": "f0.34", "cd4_gn0.55": "f0.55", "cd4_gn0.8": "f0.8",
                      "f0": "f0", "kolfmask0.05": "kolf0.05", "kolfmask0.1": "kolf0.1",
                      "kolfmask0.2": "kolf0.2", "kolfmask0.4": "kolf0.4"}
            if arm in ka_map and (held, ka_map[arm]) in ka207:
                ref = ka207[(held, ka_map[arm])]
                dc = abs(ref["cos_agg"] - e["cos_agg"])
                dp = abs(ref["pds"] - e["pds"])
                assert dc < 1e-9 and dp < 1e-9, (
                    f"KA FAILED {held}/{arm} vs 207 {ka_map[arm]}: "
                    f"dcos {dc:.3e} dpds {dp:.3e}")
                print(f"    KA PASS: {arm} == 207 {ka_map[arm]} "
                      f"(dcos {dc:.1e}, dpds {dp:.1e})", flush=True)
        del op, whitened
        _PREP_CACHE.clear()
        _BANK_CACHE.clear()


if __name__ == "__main__":
    main()
