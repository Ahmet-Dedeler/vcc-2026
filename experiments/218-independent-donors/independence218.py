"""218 Phase 0a — INDEPENDENCE, measured in the exact frame the residual is added in.

The mandate's hypothesis is that what makes CD4 pay is INDEPENDENCE from Replogle/Orion.
To regress contribution on independence, independence has to be measured on the same objects
the contribution is measured on -- the residual matrix each source contributes, on the
official gene axis, restricted to a PLOLO fold's targets and scored columns.  A pairwise
donor-table cosine on each pair's own shared axis (166's `independence.py`) is a DIFFERENT
statistic: different genes, different targets, different centring per pair.  Both are
reported; only the first can enter the regression.

Per (source, fold) it computes, on the fold's scored columns:

  cos_truth      per-target cosine to the HELD line's truth   -- transfer alignment
  decoy_truth    off-diagonal cosine to the truth             -- target-independent share
  excess_truth   cos_truth - decoy_truth                      -- target-SPECIFIC transfer
  cos_pre        per-target cosine to the base mixture        -- REDUNDANCY with what we
  excess_pre     cos_pre - decoy_pre                            already have
  sub_r          share of the source's energy inside the top-r right-singular subspace of
                 the production consensus (209's shared-subspace energy), r = 10 and 50
  orth_frac      fraction of the source's Frobenius energy orthogonal to the consensus
                 ROW space -- the geometric form of "a view we do not already hold"

KNOWN-ANSWER ROW: 166's `independence.py` gate is re-run verbatim through 166's own code,
so the three 159 line-pairs must reproduce inside their published draw sd.  The run aborts
otherwise.

    uv run python independence218.py
"""
from __future__ import annotations

import gc
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"

from plolo import FOLDS, build_fold, cd4_for  # noqa: E402
from resid218 import CD4_NAMES, CD4_ROOT, multi_source_residual, bank_for, E062  # noqa: E402
from sources218 import SOURCES  # noqa: E402

CANDIDATES = ("kolf2", "feng_ipsc", "k562_gwps", "orion_hct116", "orion_hek293t",
              "perturbmatch", "telohaec", "gse295214_h1")
PROD_TRIO = ("k562_gwps", "orion_hct116", "orion_hek293t")


def _unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def pair_stats(S, T, rows):
    """same-target / decoy cosine between two effect matrices on a common row set."""
    if rows.sum() < 10:
        return dict(n=int(rows.sum()), same=None, decoy=None, excess=None)
    A, B = _unit(S[rows]), _unit(T[rows])
    M = A @ B.T
    n = M.shape[0]
    same = float(np.diag(M).mean())
    decoy = float(M[~np.eye(n, dtype=bool)].mean())
    return dict(n=int(n), same=same, decoy=decoy, excess=same - decoy)


def known_answer_gate() -> dict:
    """166's gate, through 166's own module, unmodified."""
    spec = importlib.util.spec_from_file_location(
        "e166_ind", ROOT / "experiments/166-donor-acquisition/independence.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules["e166_ind"] = m
    spec.loader.exec_module(m)
    loaded = {}
    for k in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        loaded[k] = m.load_donor(m.DONORS[k][0])
    out: dict[str, dict] = {}
    for (a, b), expect in m.KNOWN.items():
        draws = [m.pair(loaded[a], loaded[b], np.random.default_rng(7000 + i))
                 for i in range(m.KNOWN_DRAWS)]
        got = np.array([[d["full"]["same"], d["full"]["diff"], d["resid"]["same"]]
                        for d in draws])
        mean, sd = got.mean(0), got.std(0, ddof=1)
        tol = np.maximum(m.KNOWN_TOL_SD * sd, m.KNOWN_TOL_FLOOR)
        ok = bool(np.all(np.abs(mean - np.array(expect)) <= tol))
        out[f"{a}|{b}"] = dict(expected=list(expect), mean=mean.tolist(),
                               sd=sd.tolist(), tol=tol.tolist(), pass_=ok)
        print(f"  KA {a} x {b}: expect {expect} got {np.round(mean,4).tolist()} "
              f"+-{np.round(sd,4).tolist()} -> {'PASS' if ok else 'FAIL'}")
        assert ok, f"KA gate FAILED for {a} x {b}"
    # 166's loader materialises each full donor table (two Orion tables are 2.8 GB each),
    # which alone pushes peak RSS past 062's 10 GiB streaming guard in the next stage.
    loaded.clear()
    del loaded, m
    sys.modules.pop("e166_ind", None)
    gc.collect()
    return out


def main() -> None:
    if "--ka-only" in sys.argv:
        print("=== KNOWN-ANSWER GATE (166's independence.py, verbatim) ===", flush=True)
        (WORK / "ka218.json").write_text(json.dumps(known_answer_gate(), indent=2))
        return
    # 166's gate materialises three full donor tables (~6 GB) and 062's streaming guard
    # reads `ru_maxrss`, a HIGH-WATER MARK that gc cannot lower -- so the gate must run in
    # its own process or every later stage inherits its peak and aborts.
    import subprocess
    r = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--ka-only"],
                       cwd=str(HERE))
    assert r.returncode == 0, "known-answer gate FAILED (see above)"
    ka = json.loads((WORK / "ka218.json").read_text())

    rows = []
    for held in FOLDS:
        fold = build_fold(held)
        targets = fold["targets"]
        sc = fold["scored"]
        truth = np.asarray(fold["truth"], dtype=np.float64)[:, sc]
        pre = np.asarray(fold["pre"], dtype=np.float64)[:, sc]
        cd4 = cd4_for(targets)

        # the consensus row space, for the subspace / orthogonality metrics
        U, sv, Vt = np.linalg.svd(pre, full_matrices=False)
        energy = float((pre * pre).sum())
        print(f"[{held}] pre {pre.shape}, top-10 energy share "
              f"{float((sv[:10]**2).sum()/max(energy,1e-300)):.4f}", flush=True)

        sources = {"cd4": (np.where(np.asarray(cd4.gene_votes) > 0,
                                    np.asarray(cd4.mean_effect), 0.0)[:, sc],
                           np.asarray(cd4.target_condition_count) > 0)}
        for name in CANDIDATES:
            if name == held:
                continue
            r = multi_source_residual([bank_for(name)], targets, "gene_norm")
            sources[name] = (np.where(np.asarray(r.gene_votes) > 0,
                                      np.asarray(r.mean_effect), 0.0)[:, sc],
                             np.asarray(r.target_condition_count) > 0)
        sources["self_pre"] = (pre.copy(), np.ones(len(targets), dtype=bool))

        truth_rows = np.ones(len(targets), dtype=bool)
        for name, (S, avail) in sources.items():
            both = avail & truth_rows
            t = pair_stats(S, truth, both)
            p = pair_stats(S, pre, avail)
            live = S[avail]
            e_tot = float((live * live).sum())
            proj = live @ Vt.T
            sub = {f"sub_{r}": float((proj[:, :r] ** 2).sum() / max(e_tot, 1e-300))
                   for r in (10, 50, 100)}
            meta = SOURCES.get(name, {})
            rows.append(dict(
                fold=held, source=name, n_carrying=int(avail.sum()),
                coverage=float(avail.mean()),
                cos_truth=t["same"], decoy_truth=t["decoy"], excess_truth=t["excess"],
                cos_pre=p["same"], decoy_pre=p["decoy"], excess_pre=p["excess"],
                **sub, orth_frac=float(1.0 - sub["sub_100"]),
                platform=meta.get("platform"), lineage=meta.get("lineage"),
                transformation=meta.get("transformation"), lab=meta.get("lab"),
                semantics=meta.get("semantics")))
            print(f"  {name:16s} carry {int(avail.sum()):3d} "
                  f"cos_truth {t['same']:+.5f} excess_truth {t['excess']:+.5f} "
                  f"cos_pre {p['same']:+.5f} excess_pre {p['excess']:+.5f} "
                  f"sub10 {sub['sub_10']:.4f} orth {1-sub['sub_100']:.4f}", flush=True)

    (WORK / "independence218.json").write_text(json.dumps(
        dict(known_answer=ka, rows=rows), indent=2))
    print(f"\nwrote {WORK/'independence218.json'} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
