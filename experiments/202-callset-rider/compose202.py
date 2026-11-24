"""202 falsifier 3 — the CALL-SET COMPOSITION, on the real challenge contexts.

`probe191.py` already measures falsifiers 2 (`n_pred`) and 4 (the delivered order) and is
reused unchanged.  Neither it nor anything else in the repo records **which genes** the call
set contains, and falsifier 3 is exactly that question, so this is the one new measurement.

It is `probe191.py`'s I/O loop verbatim -- the same 30-target draw (`default_rng(191)`), the
same 5-CPM scorer gate, the same `wilcoxon_z` kernel out of `vcc2026.rc_emit`, the same
target-row exclusion -- with the per-block DE call recorded as a MASK instead of a count.
Against that mask it reports:

  * the PRESCRIBED set the arm asks for, recomputed here from the effect matrix and the
    context's own control cells, independently of the emitter that built the artifact:
    `top 400 by |H|  U  top 2,200 by control CPM`, both rankings over the 5-CPM universe;
  * how much of it BH actually delivered (`jaccard`, `missed`, `unsuppressible`);
  * the overlap with `rc_hybrid`'s own delivered call set -- the incumbent this arm is one
    variable from;
  * **how many genes the two halves of the union share**, which is the arm's whole premise.
    192 §5.3's claim is that a set keyed to `|H|` and a set keyed to detectability are
    DIFFERENT objects doing different jobs; if the halves largely coincided the union would
    be a re-spelling of one of them and the +0.011 `reach` half of the prediction would have
    no mechanism.  Reported both ways round: where the `|H|` head sits in the CPM order, and
    where the CPM half sits in the `|H|` order.

No ground truth is read and none is needed: every object here is either ours (the effect
matrix, the artifact) or the scorer's own comparison group (`data/challenge/context_*.h5ad`).

    uv run python compose202.py <label>=<pred.h5ad> [...] [--perts N]
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from scipy.stats import false_discovery_control, norm

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.rc_emit import ctrl_tiesums, wilcoxon_z  # noqa: E402
from vcc2026.donor_transfer import TS  # noqa: E402

BUNDLE = ROOT / "data/challenge"
EFFECTS = ROOT / "experiments/178-composed-riders/work/eff.a15.npy"
ALPHA, SCORER_MIN_CPM, BUILDER_MIN_CPM = 0.05, 5.0, 10.0
# The arm, as preregistered.  Read from the builder so the two can never drift apart.
sys.path.insert(0, str(ROOT / "scripts"))


def _arm_constants():
    """`RC_CALLSET_HEAD` / `RC_CALLSET_CPM` off the builder's source, without importing it."""
    src = (ROOT / "scripts/make_submission_pm_cd4.py").read_text()
    out = {}
    for name in ("RC_CALLSET_HEAD", "RC_CALLSET_CPM"):
        for line in src.splitlines():
            if line.startswith(name + " ="):
                out[name] = int(line.split("=", 1)[1].split("#")[0].strip())
    if len(out) != 2:
        raise SystemExit("could not read the arm's constants out of the builder")
    return out["RC_CALLSET_HEAD"], out["RC_CALLSET_CPM"]


N_HEAD, N_CPM = _arm_constants()


def norm_log(x):
    sums = np.asarray(x.sum(1)).ravel().astype(np.float64)
    sums[sums == 0] = 1.0
    return np.log1p(x * (1e6 / sums[:, None]))


def jac(a, b):
    union = int((a | b).sum())
    return float((a & b).sum() / union) if union else float("nan")


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_perts = 30
    for flag in sys.argv[1:]:
        if flag.startswith("--perts"):
            n_perts = int(flag.split("=")[1])
    packages = {a.split("=", 1)[0]: Path(a.split("=", 1)[1]) for a in args}
    if "rc_callset" not in packages:
        raise SystemExit("one of the packages must be labelled `rc_callset`")

    effects = np.load(EFFECTS)
    effect_sha = hashlib.sha256(
        np.ascontiguousarray(effects, dtype="<f8").tobytes()
    ).hexdigest()
    print(f"effect matrix {EFFECTS}  sha256 {effect_sha[:16]}")
    print(f"prescribed call set: top {N_HEAD} by |H|  U  top {N_CPM} by control CPM "
          f"= {N_HEAD + N_CPM}")
    for label, path in packages.items():
        meta = json.loads(Path(str(path) + ".meta.json").read_text())
        built = meta["effects"]["final_sha256"]
        print(f"  {label:<12} package effects.final_sha256 {built[:16]}  "
              f"{'MATCHES' if built == effect_sha else 'DIFFERS'}")
        if built != effect_sha:
            raise SystemExit(f"{label}: wrong effect matrix; the |H| key would be wrong")

    genes = np.asarray(
        [line.strip() for line in (BUNDLE / "gene_names.csv").read_text().splitlines()[1:]]
    )
    perts = np.asarray(
        [line.split(",")[0] for line in
         (BUNDLE / "pert_counts.csv").read_text().splitlines()[1:]]
    )
    assert genes.size == effects.shape[1] and perts.size == effects.shape[0]

    handles, obs = {}, {}
    for label, path in packages.items():
        adata = ad.read_h5ad(path, backed="r")
        obs[label] = (np.asarray(adata.obs["target_gene"].astype(str).values),
                      np.asarray(adata.obs["context"].astype(str).values))
        adata.file.close()
        handles[label] = h5py.File(path, "r")

    rows: list[dict] = []
    t0 = time.time()
    for context in ("A", "B", "C"):
        controls = ad.read_h5ad(BUNDLE / f"context_{context}.h5ad")
        matrix = controls.X
        dense = np.asarray(matrix.todense() if hasattr(matrix, "todense") else matrix,
                           dtype=np.float64)
        del controls, matrix
        pooled = dense.sum(0)
        ctrl = norm_log(dense)
        del dense
        cpm = np.expm1(ctrl).mean(0)
        keep = cpm > SCORER_MIN_CPM
        cs = np.ascontiguousarray(np.sort(ctrl[:, keep].T, axis=1))
        tie = ctrl_tiesums(cs)
        # `ref_lin` IS the emitter's CPM key: the per-gene mean of per-cell CPM over this
        # context's own control cells, restricted to the scorer's DE universe.
        ref_lin = np.expm1(ctrl[:, keep]).mean(0)
        del ctrl
        b_ctrl = np.log1p(TS * pooled / pooled.sum())
        expression = np.expm1(b_ctrl)
        gate_cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
        builder_gate = (gate_cpm > BUILDER_MIN_CPM) & (b_ctrl > 0)
        key_all = np.where(builder_gate[None, :], effects, 0.0)
        n_keep = int(keep.sum())
        del pooled
        print(f"\nctx {context}: {n_keep} genes past the {SCORER_MIN_CPM}-CPM scorer gate "
              f"[{time.time() - t0:.0f}s]", flush=True)

        chosen = sorted(perts.tolist())
        picks = np.random.default_rng(191).choice(
            len(chosen), min(n_perts, len(chosen)), replace=False
        )
        targets = [chosen[i] for i in sorted(picks)]
        pos = {str(p): i for i, p in enumerate(perts)}

        # ---- the prescribed set, recomputed here, once per (context, target) ------------
        # `np.argsort(-x, kind="stable")` is the emitter's own tie rule: equal keys keep
        # ascending gene index, so the set is a deterministic function of the inputs.
        cpm_order_unrestricted = np.argsort(-ref_lin, kind="stable")
        cpm_top_alone = np.zeros(n_keep, dtype=bool)
        cpm_top_alone[cpm_order_unrestricted[:N_CPM]] = True
        cpm_rank_of_gene = np.empty(n_keep, dtype=np.int64)
        cpm_rank_of_gene[cpm_order_unrestricted] = np.arange(n_keep)

        called: dict[str, dict[str, np.ndarray]] = {label: {} for label in packages}
        for label, path in packages.items():
            labels, contexts = obs[label]
            handle = handles[label]
            data, indices, indptr = (handle["X"]["data"], handle["X"]["indices"],
                                     handle["X"]["indptr"])
            for target in targets:
                sel = np.flatnonzero((contexts == context) & (labels == target))
                assert sel.size == 400, (label, context, target, sel.size)
                lo, hi = int(indptr[sel[0]]), int(indptr[sel[-1] + 1])
                block = np.zeros((sel.size, genes.size), dtype=np.float64)
                cols = indices[lo:hi].astype(np.int64)
                vals = data[lo:hi].astype(np.float64)
                ptr = np.asarray(indptr[sel[0]:sel[-1] + 2]) - lo
                for r in range(sel.size):
                    block[r, cols[ptr[r]:ptr[r + 1]]] = vals[ptr[r]:ptr[r + 1]]
                normed = norm_log(block)
                del block
                z = wilcoxon_z(cs, tie, np.ascontiguousarray(normed[:, keep].T))
                del normed
                p = 2 * norm.sf(np.abs(z))
                q = false_discovery_control(p, method="bh")
                sig = q < ALPHA
                own = np.flatnonzero(genes[keep] == target)
                if own.size:
                    sig[own[0]] = False
                called[label][target] = sig
            print(f"  {label:<12} delivered call sets read", flush=True)

        for target in targets:
            key = np.abs(key_all[pos[target]])[keep]
            h_order = np.argsort(-key, kind="stable")
            head = np.zeros(n_keep, dtype=bool)
            head[h_order[:N_HEAD]] = True
            h_rank_of_gene = np.empty(n_keep, dtype=np.int64)
            h_rank_of_gene[h_order] = np.arange(n_keep)
            cpm_rank = np.where(head, -np.inf, ref_lin)
            rest = np.zeros(n_keep, dtype=bool)
            rest[np.argsort(-cpm_rank, kind="stable")[:N_CPM]] = True
            prescribed = head | rest
            assert int(prescribed.sum()) == N_HEAD + N_CPM

            delivered = called["rc_callset"][target]
            row = dict(
                context=context, target=target,
                n_keep=n_keep,
                prescribed=int(prescribed.sum()),
                delivered=int(delivered.sum()),
                jac_delivered_prescribed=jac(delivered, prescribed),
                missed=int((prescribed & ~delivered).sum()),
                unsuppressible=int((delivered & ~prescribed).sum()),
                head=int(head.sum()), rest=int(rest.sum()),
                halves_shared=int((head & cpm_top_alone).sum()),
                head_median_cpm_rank=float(np.median(cpm_rank_of_gene[head])),
                rest_median_H_rank=float(np.median(h_rank_of_gene[rest])),
                rest_beyond_head_H_rank=int((h_rank_of_gene[rest] >= N_HEAD).sum()),
            )
            if "rc_hybrid" in called:
                inc = called["rc_hybrid"][target]
                row.update(
                    rc_hybrid_n_pred=int(inc.sum()),
                    overlap_delivered_rc_hybrid=int((delivered & inc).sum()),
                    jac_delivered_rc_hybrid=jac(delivered, inc),
                    overlap_prescribed_rc_hybrid=int((prescribed & inc).sum()),
                    jac_prescribed_rc_hybrid=jac(prescribed, inc),
                    head_in_rc_hybrid=int((head & inc).sum()),
                    rest_in_rc_hybrid=int((rest & inc).sum()),
                )
            rows.append(row)
        del cs, tie, ref_lin, key_all
    for handle in handles.values():
        handle.close()

    def med(field):
        return float(np.median([r[field] for r in rows]))

    print("\n=== FALSIFIER 3 — CALL-SET COMPOSITION on the challenge "
          f"({len(rows)} blocks over 3 contexts) ===")
    print(f"  prescribed |top {N_HEAD} by |H|  U  top {N_CPM} by CPM|   "
          f"median {med('prescribed'):.0f}   (must be {N_HEAD + N_CPM})")
    print(f"  delivered by BH                                median "
          f"{med('delivered'):.0f}")
    print(f"    jaccard(delivered, prescribed)               median "
          f"{med('jac_delivered_prescribed'):.4f}")
    print(f"    prescribed but MISSED                        median {med('missed'):.0f}"
          f"   max {max(r['missed'] for r in rows)}")
    print(f"    delivered but UNSUPPRESSIBLE                 median "
          f"{med('unsuppressible'):.0f}   p90 "
          f"{np.percentile([r['unsuppressible'] for r in rows], 90):.0f}   max "
          f"{max(r['unsuppressible'] for r in rows)}")
    print("\n  THE TWO HALVES OF THE UNION")
    print(f"    |top {N_HEAD} by |H||                                {med('head'):.0f}")
    print(f"    |top {N_CPM} by CPM, head excluded|                 {med('rest'):.0f}")
    print(f"    |halves SHARED| (head n unrestricted top-{N_CPM} CPM)  median "
          f"{med('halves_shared'):.0f}")
    print(f"    the |H| head's median rank in the CPM order      "
          f"{med('head_median_cpm_rank'):.0f} of {med('n_keep'):.0f}")
    print(f"    the CPM half's median rank in the |H| order      "
          f"{med('rest_median_H_rank'):.0f} of {med('n_keep'):.0f}")
    print(f"    CPM-half genes the |H| head would never reach    "
          f"median {med('rest_beyond_head_H_rank'):.0f} of {N_CPM}")
    if "rc_hybrid" in packages:
        print("\n  AGAINST `rc_hybrid`'s OWN DELIVERED CALL SET")
        print(f"    `rc_hybrid` n_pred                             median "
              f"{med('rc_hybrid_n_pred'):.0f}")
        print(f"    |delivered n rc_hybrid|                        median "
              f"{med('overlap_delivered_rc_hybrid'):.0f}")
        print(f"    jaccard(delivered, rc_hybrid)                  median "
              f"{med('jac_delivered_rc_hybrid'):.4f}")
        print(f"    |prescribed n rc_hybrid|                       median "
              f"{med('overlap_prescribed_rc_hybrid'):.0f}")
        print(f"    of the {N_HEAD}-gene |H| head, in rc_hybrid's set  median "
              f"{med('head_in_rc_hybrid'):.0f}")
        print(f"    of the {N_CPM}-gene CPM half, in rc_hybrid's set  median "
              f"{med('rest_in_rc_hybrid'):.0f}")

    out = HERE / "work/compose202.json"
    out.write_text(json.dumps(dict(
        n_head=N_HEAD, n_cpm=N_CPM, n_perts=n_perts, rows=rows,
        medians={k: med(k) for k in rows[0] if isinstance(rows[0][k], (int, float))},
    ), indent=2))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
