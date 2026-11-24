"""200 — rebuild the h1-test effect matrix per ablation, so the DE members can be priced.

WHY A REBUILD IS NECESSARY.  Every h1-test path this project holds (093 / 095 / 157 / 164 /
169) re-emits from ONE built prediction's per-gene panel totals -- it holds the PROFILE fixed
and moves only the emission law.  That prices the emission-law components (the tiered spread
policy, `--overdispersion`) and nothing else.  Half the shipped stack is geometry, and its
`reach` / `fid` / `jac` / `nmae` have never been measured at all.

WHY h1-test AND NOT THE ORION PAIR.  157 §1: h1-test reproduces BOTH board raw `reach` points
to within 0.014 (a8 0.1890 board / 0.2029 here; a15 0.1811 / 0.1828) while repaired HCT116
misses a15's by 2.6x.  AGENTS.md §6.4 adds that the Orion harnesses carry a target-independent
blob worth 26.4% / 44.1% of their truth energy against h1-test's 3.6%.

THE TWO TRAPS, BOTH RESPECTED.
  * `data/donors/h1_2025` covers 85/85 of h1-test's targets and is pure oracle leakage.  The
    donor set here is `k562_gwps` + `orion_hek293t` ONLY, and h1_2025 never enters, not as a
    donor and not in any ordering key.
  * `Surface(n_real_target=N)` is the wrong power-matching instrument for `reach` (it deletes
    the small-`n_real` tercile carrying 98% of our reach).  It is not used anywhere here.

THE OPERATOR.  `experiments/095-clean-harness/out-h1test/whitener_union.npz`, sha256
`f80fe8d6ba59...`, fitted on k562_gwps + orion_hek293t with all 85 h1-test targets excluded
from the basis.  Its 17,745 columns are RE-DERIVED here from the donor banks and the harness
gene axis and checked against the operator's own recorded `wcols_sha256`
(`9343b796ff2c40f5...`) -- rule 11's "never infer an operator from a filename", run as an
assertion rather than a claim.

THE AMPLITUDE CONVENTION.  The artifact 157/164/169 calibrated `reach` against was built by
031's path, which applies `alpha` LAST as a renormalisation of the median row norm
(`E *= alpha * effect_scale / median|E|`), not before whitening as production does.  Every arm
here is renormalised to the SAME median row norm as that artifact, so a geometry ablation is a
pure SHAPE change and carries no hidden amplitude confound.  `alpha` itself is therefore not
an h1-test arm; 187 closed it on the two members that can see it.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
P162 = os.path.join(ROOT, "experiments/162-retrodicting-validator")
sys.path.insert(0, P162)
sys.path.insert(0, HERE)

import plolo                                                  # noqa: E402
from plolo import CD4_NAMES, CD4_ROOT, DONOR_PATHS, E062, E064, PM_PATH  # noqa: E402

H1_REAL = f"{ROOT}/data/harness/h1-test.h5ad"
H1_OP = f"{ROOT}/experiments/095-clean-harness/out-h1test/whitener_union.npz"
H1_OP_SHA = "f80fe8d6ba59"          # from work/inputs200.log
H1_WCOLS_SHA = "9343b796ff2c40f5c170976acde1b94547069133289473b11a54158547c2e3f9"
H1_SOURCES = ("k562_gwps", "orion_hek293t")     # h1_2025 is oracle leakage; never here
EFFECT_SCALE = 6.817
os.makedirs(WORK, exist_ok=True)


def harness_axis() -> tuple[list[str], list[str]]:
    """h1-test's gene axis and its 85 target genes, read from the harness itself."""
    import h5py
    cache = os.path.join(WORK, "h1-axis.json")
    if os.path.exists(cache):
        d = json.loads(open(cache).read())
        return d["genes"], d["targets"]
    with h5py.File(H1_REAL, "r") as h:
        node = h["var/_index"]
        raw = node["values"][:] if isinstance(node, h5py.Group) else node[:]
        genes = [g.decode() if isinstance(g, bytes) else str(g) for g in raw]
        o = h["obs/target_gene"]
        if isinstance(o, h5py.Group):
            cats = [c.decode() if isinstance(c, bytes) else str(c) for c in o["categories"][:]]
            labs = np.asarray(cats)[o["codes"][:]]
        else:
            labs = np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in o[:]])
    targets = sorted(set(labs.tolist()) - {"non-targeting"})
    json.dump(dict(genes=genes, targets=targets), open(cache, "w"))
    return genes, targets


def parts(limit_gib: float = 24.0) -> dict:
    """Per-donor effects for the 85 h1-test targets on the harness's own gene axis."""
    genes, targets = harness_axis()
    cache = os.path.join(WORK, "parts-h1test.npz")
    pos = {g: i for i, g in enumerate(genes)}
    panel_positions = np.asarray([pos[g] for g in targets if g in pos], dtype=np.int64)
    if os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        print("[h1-test] parts cache hit", flush=True)
        return dict(genes=genes, targets=targets,
                    src={n: z[f"src_{n}"].astype(np.float64) for n in H1_SOURCES},
                    avail={n: z[f"avail_{n}"] for n in H1_SOURCES},
                    mean_on_axis={n: z[f"mean_{n}"] for n in H1_SOURCES},
                    pm_effects=z["pm_effects"].astype(np.float64), pm_avail=z["pm_avail"],
                    shared=z["shared"], wcols=z["wcols"],
                    panel_positions=panel_positions, audit=json.loads(str(z["audit"])))

    banks = {n: E062.Bank.load(n, DONOR_PATHS[n]) for n in H1_SOURCES}
    prepared = {n: E062.prepare_main_source(b, targets, genes, panel_positions,
                                            chunk_rows=256, max_rss_gib=limit_gib)
                for n, b in banks.items()}
    pm_bank = E064.load_perturbmatch(PM_PATH)
    pm_prep = E062.prepare_main_source(pm_bank, targets, genes, panel_positions,
                                       chunk_rows=256, max_rss_gib=limit_gib)
    shared = E064.shared_k562_pm_columns(banks["k562_gwps"], pm_bank, genes, set(targets))
    wcols = E064.operator_columns([banks[n] for n in H1_SOURCES], genes, set(targets))
    got = hashlib.sha256(np.asarray(wcols, dtype=np.int64).tobytes()).hexdigest()
    if got != H1_WCOLS_SHA:
        raise SystemExit(f"ABORT: h1-test operator columns {got[:16]} != {H1_WCOLS_SHA[:16]}")
    audit = dict(genes=len(genes), targets=len(targets), wcols=int(len(wcols)),
                 wcols_sha256=got, operator=H1_OP, operator_sha256_prefix=H1_OP_SHA,
                 per_donor_target_counts={n: int(prepared[n].target_available.sum())
                                          for n in H1_SOURCES},
                 pm_available=int(pm_prep.target_available.sum()),
                 shared_columns=int(np.asarray(shared).size),
                 raw_median_norm={n: float(prepared[n].raw_median_norm) for n in H1_SOURCES})
    np.savez_compressed(
        cache, audit=json.dumps(audit), wcols=wcols, shared=np.asarray(shared),
        pm_effects=pm_prep.effects.astype(np.float32), pm_avail=pm_prep.target_available,
        **{f"src_{n}": prepared[n].effects.astype(np.float32) for n in H1_SOURCES},
        **{f"avail_{n}": prepared[n].target_available for n in H1_SOURCES},
        **{f"mean_{n}": prepared[n].mean_on_axis for n in H1_SOURCES})
    print(f"[h1-test] parts built: {json.dumps(audit)}", flush=True)
    return dict(genes=genes, targets=targets,
                src={n: prepared[n].effects for n in H1_SOURCES},
                avail={n: prepared[n].target_available for n in H1_SOURCES},
                mean_on_axis={n: prepared[n].mean_on_axis for n in H1_SOURCES},
                pm_effects=pm_prep.effects, pm_avail=pm_prep.target_available,
                shared=np.asarray(shared), wcols=wcols,
                panel_positions=panel_positions, audit=audit)


def cd4(targets) -> object:
    cache = os.path.join(WORK, "cd4-h1test.npz")
    genes, _ = harness_axis()
    if os.path.exists(cache):
        z = np.load(cache)
        return E062.CD4Prepared(mode="gene_norm", mean_effect=z["mean_effect"],
                                gene_votes=z["gene_votes"], target_condition_count=z["tcc"],
                                condition_details=[])
    banks = [E062.Bank.load(n, CD4_ROOT / n, require_cd4_scale=True) for n in CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, tuple(targets), genes, set(targets), "gene_norm",
                                 chunk_rows=256, max_rss_gib=24.0)
    complete = np.asarray(prep.target_condition_count) == len(CD4_NAMES)
    me = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    gv = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    me[~complete] = 0.0
    gv[~complete] = 0
    tcc = np.asarray(prep.target_condition_count, dtype=np.uint8).copy()
    np.savez_compressed(cache, mean_effect=me, gene_votes=gv, tcc=tcc)
    print(f"[h1-test] cd4: {int(complete.sum())}/{len(targets)} complete-mean3 targets",
          flush=True)
    return E062.CD4Prepared(mode="gene_norm", mean_effect=me, gene_votes=gv,
                            target_condition_count=tcc, condition_details=[])


class H1Operator:
    """095's held-safe h1-test operator, loaded once and applied like production's."""

    def __init__(self, path: str = H1_OP):
        with np.load(path) as z:
            self.evals = np.asarray(z["evals"], dtype=np.float64)
            self.evecs = np.asarray(z["evecs"], dtype=np.float64)

    def apply(self, M, wcols, gamma, lam):
        out = np.asarray(M, dtype=np.float64).copy()
        before = float(np.median(np.linalg.norm(out[:, wcols], axis=1))) or 1.0
        w = (self.evals + lam) ** (-gamma / 2.0)
        t = ((out[:, wcols] @ self.evecs) * w) @ self.evecs.T
        after = float(np.median(np.linalg.norm(t, axis=1))) or 1.0
        out[:, wcols] = t * (before / after)
        return out


def effects_for(cfg, p, c4, op, target_norm: float) -> np.ndarray:
    """One arm's h1-test effect matrix, renormalised to a fixed median row norm."""
    from folds200 import consensus
    from plolo import solve_fraction
    pre = consensus(p, pm_blend=cfg["pm"], weights=cfg["weights"])
    if cfg["gamma"] is None:
        whitened = np.asarray(pre, dtype=np.float64).copy()
    else:
        whitened = op.apply(pre, p["wcols"], cfg["gamma"], cfg["lam"])
    if cfg["cd4_f"] is not None:
        frac, _ = solve_fraction(whitened, c4, cfg["cd4_f"])
        final, _ = E062.add_post_whiten_residual(whitened, c4, frac, p["wcols"])
    elif cfg["cd4_frac"] is not None:
        final, _ = E062.add_post_whiten_residual(whitened, c4, float(cfg["cd4_frac"]),
                                                 p["wcols"])
    else:
        final = whitened
    final = np.array(final, dtype=np.float64, copy=True)
    final[:, p["panel_positions"]] = 0.0
    med = float(np.median(np.linalg.norm(final, axis=1))) or 1.0
    return final * (target_norm / med)


if __name__ == "__main__":
    p = parts()
    c4 = cd4(p["targets"])
    print(json.dumps(p["audit"], indent=1))
    print(f"cd4 complete-mean3 targets: "
          f"{int((np.asarray(c4.target_condition_count) == 3).sum())}/{len(p['targets'])}")
