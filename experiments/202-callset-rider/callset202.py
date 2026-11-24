"""202 falsifier 3 — the CALL-SET COMPOSITION on the real challenge contexts.

The arm's whole logic is that the two halves of the union do different jobs: the `|H|` head
carries the order `reach` walks, the control-CPM half carries the membership `jac` reads.
This measures that, on `data/challenge/context_{A,B,C}.h5ad`, and it measures what the
emitter actually delivered against what was prescribed.

It reuses `probe191.py`'s conventions exactly -- the same effect matrix identified by SHA256
against each package's own `effects.final_sha256` (rule 11, never by filename), the same
builder gate, the same scorer gate, the same target subsample -- so its numbers sit in the
same frame as 191's.  `probe191.py` itself is run unmodified for falsifiers 2 and 4.

    uv run python callset202.py <label>=<pred.h5ad> ... [--perts N]

The first label is treated as the incumbent (`rc_hybrid`) and the second as the arm.
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
N_HEAD, N_CPM = 400, 2200


def norm_log(x):
    sums = np.asarray(x.sum(1)).ravel().astype(np.float64)
    sums[sums == 0] = 1.0
    return np.log1p(x * (1e6 / sums[:, None]))


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_perts = 40
    for flag in sys.argv[1:]:
        if flag.startswith("--perts"):
            n_perts = int(flag.split("=")[1])
    packages = {a.split("=", 1)[0]: Path(a.split("=", 1)[1]) for a in args}
    labels = list(packages)
    inc_label, arm_label = labels[0], labels[1]

    effects = np.load(EFFECTS)
    effect_sha = hashlib.sha256(
        np.ascontiguousarray(effects, dtype="<f8").tobytes()
    ).hexdigest()
    print(f"effect matrix {EFFECTS}  shape {effects.shape}  sha256 {effect_sha[:16]}")
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

    out: dict[str, dict] = {}
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
        ref_lin = np.expm1(ctrl[:, keep]).mean(0)     # the emitter's `control.ref_lin`
        del ctrl
        b_ctrl = np.log1p(TS * pooled / pooled.sum())
        expression = np.expm1(b_ctrl)
        gate_cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
        builder_gate = (gate_cpm > BUILDER_MIN_CPM) & (b_ctrl > 0)
        key_all = np.where(builder_gate[None, :], effects, 0.0)
        n_keep = int(keep.sum())
        print(f"\nctx {context}: {cs.shape[1]} controls, {n_keep} genes past the "
              f"{SCORER_MIN_CPM}-CPM scorer gate  [{time.time() - t0:.0f}s]", flush=True)
        del pooled

        chosen = sorted(perts.tolist())
        picks = np.random.default_rng(191).choice(
            len(chosen), min(n_perts, len(chosen)), replace=False
        )
        targets = [chosen[i] for i in sorted(picks)]
        pos = {str(p): i for i, p in enumerate(perts)}

        # the CPM half is a per-CONTEXT object: `ref_lin` does not depend on the target, so
        # the unrestricted top-N_CPM by control CPM is one set for all 300 perturbations.
        raw_cpm = np.zeros(n_keep, bool)
        raw_cpm[np.argsort(-ref_lin, kind="stable")[:N_CPM]] = True

        sig = {}
        for label, path in packages.items():
            labels_, contexts_ = obs[label]
            handle = handles[label]
            data, indices, indptr = (handle["X"]["data"], handle["X"]["indices"],
                                     handle["X"]["indptr"])
            per_target = {}
            for target in targets:
                sel = np.flatnonzero((contexts_ == context) & (labels_ == target))
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
                s = q < ALPHA
                own = np.flatnonzero(genes[keep] == target)
                if own.size:
                    s[own[0]] = False
                per_target[target] = s
            sig[label] = per_target
            print(f"  {label:<12} read", flush=True)

        rows = []
        for target in targets:
            key = np.abs(key_all[pos[target]])[keep]
            head = np.zeros(n_keep, bool)
            head[np.argsort(-key, kind="stable")[:N_HEAD]] = True
            want = head.copy()
            want[np.argsort(-np.where(head, -np.inf, ref_lin),
                            kind="stable")[:N_CPM]] = True
            inc = sig[inc_label][target]
            got = sig[arm_label][target]
            rows.append(dict(
                target=target,
                n_keep=n_keep,
                head=int(head.sum()),
                cpm_half=int(raw_cpm.sum()),
                halves_shared=int((head & raw_cpm).sum()),
                prescribed=int(want.sum()),
                delivered=int(got.sum()),
                delivered_jaccard=float((got & want).sum() / max((got | want).sum(), 1)),
                unsuppressible=int((got & ~want).sum()),
                missed=int((want & ~got).sum()),
                n_pred_incumbent=int(inc.sum()),
                overlap_incumbent=int((want & inc).sum()),
                jaccard_incumbent=float((want & inc).sum() / max((want | inc).sum(), 1)),
                head_in_incumbent=int((head & inc).sum()),
                cpm_in_incumbent=int((raw_cpm & inc).sum()),
                delivered_overlap_incumbent=int((got & inc).sum()),
            ))
        out[context] = rows
        med = {k: float(np.median([r[k] for r in rows])) for k in rows[0] if k != "target"}
        print(f"  |H| head {med['head']:.0f}   top-{N_CPM} cpm {med['cpm_half']:.0f}   "
              f"halves shared {med['halves_shared']:.0f}   union {med['prescribed']:.0f}")
        print(f"  delivered {med['delivered']:.0f} (jaccard to prescribed "
              f"{med['delivered_jaccard']:.4f}, unsuppressible {med['unsuppressible']:.0f}, "
              f"missed {med['missed']:.0f})")
        print(f"  {inc_label} n_pred {med['n_pred_incumbent']:.0f}   overlap "
              f"{med['overlap_incumbent']:.0f}   jaccard {med['jaccard_incumbent']:.4f}   "
              f"|H| head inside it {med['head_in_incumbent']:.0f}/{N_HEAD}   "
              f"cpm half inside it {med['cpm_in_incumbent']:.0f}/{N_CPM}")
        del cs, tie, ref_lin, key_all

    every = [r for rows in out.values() for r in rows]
    med = {k: float(np.median([r[k] for r in every])) for k in every[0] if k != "target"}
    print("\nPOOLED over the three contexts (median over "
          f"{len(every)} (context, perturbation) pairs)")
    for k, v in med.items():
        print(f"  {k:<28}{v:>12.4f}")
    with open(HERE / "work/callset202.json", "w") as fh:
        json.dump(dict(per_context=out, pooled_median=med), fh, indent=1, default=float)
    print(f"\nwrote {HERE}/work/callset202.json")


if __name__ == "__main__":
    main()
