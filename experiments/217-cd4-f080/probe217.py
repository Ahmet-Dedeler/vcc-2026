"""217 falsifier 5 (and 202's falsifier 4) — copied verbatim from 202's `probe191.py`
except for ONE guard, marked `210 CHANGE` in the body.

Original header follows.

191 falsifiers 3 and 4 — the DE table both packages will actually be scored on.

This is not a proxy.  The scorer recomputes DE from our counts by a Wilcoxon rank-sum against
**the reference's real control cells**, and for the challenge those are exactly the cells that
ship in `data/challenge/context_{A,B,C}.h5ad` -- the same file we were given as model input
(INBOX-3).  So every quantity below is the scorer's own, computed with its own kernel, with no
ground truth needed and none used:

  FALSIFIER 3   `n_pred`, our BH-significant call count per (context, perturbation).
  FALSIFIER 4   Spearman(-log10 p, |H|) INSIDE our call set -- the delivered order.
                185 measured a15 at -0.391 and `rc_hybrid` at +0.984 on h1-test.  If the
                order does not deliver on the REAL contexts the arm is pointless.

Also reported, because the same pass pays for them: the expression correlation (the defect
the arm removes), the fraction of the call set UP, and the up-fraction by confidence decile.

`|H|` is the incumbent's own gated predicted effect -- the same object the emitter keys on --
recomputed here from the effect matrix and the context's `b_ctrl`, with the effect matrix
identified by SHA256 against the package's own `effects.final_sha256`, never by filename
(rule 11).

    uv run python probe191.py <label>=<pred.h5ad> [<label>=<pred.h5ad> ...] [--perts N]
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
from scipy.stats import false_discovery_control, norm, spearmanr

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
from vcc2026.rc_emit import ctrl_tiesums, wilcoxon_z  # noqa: E402

BUNDLE = ROOT / "data/challenge"
EFFECTS = ROOT / "experiments/178-composed-riders/work/eff.a15.npy"
from vcc2026.donor_transfer import TS  # the builder's own pseudobulk scale
ALPHA, SCORER_MIN_CPM, BUILDER_MIN_CPM = 0.05, 5.0, 10.0


def norm_log(x):
    sums = np.asarray(x.sum(1)).ravel().astype(np.float64)
    sums[sums == 0] = 1.0
    return np.log1p(x * (1e6 / sums[:, None]))


def main() -> None:
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    n_perts = 30
    for flag in sys.argv[1:]:
        if flag.startswith("--perts"):
            n_perts = int(flag.split("=")[1])
    packages = {a.split("=", 1)[0]: Path(a.split("=", 1)[1]) for a in args}

    # ---- the ordering key, identified by hash and not by name (rule 11) ----------------
    import hashlib

    effects = np.load(EFFECTS)
    effect_sha = hashlib.sha256(
        np.ascontiguousarray(effects, dtype="<f8").tobytes()
    ).hexdigest()
    print(f"effect matrix {EFFECTS}  shape {effects.shape}  sha256 {effect_sha[:16]}")
    for label, path in packages.items():
        meta = json.loads(Path(str(path) + ".meta.json").read_text())
        built = meta["effects"]["final_sha256"]
        mark = "MATCHES" if built == effect_sha else "DIFFERS"
        print(f"  {label:<12} package effects.final_sha256 {built[:16]}  {mark}")
        if built != effect_sha:
            # 210 CHANGE, and the only one in this file.  191 and 202 were profile-PRESERVING
            # arms, so the cached `eff.a15.npy` WAS every package's own effect matrix and a
            # mismatch could only mean the wrong file.  This arm moves the profile by
            # construction, so `cd4f034`'s `final_sha256` MUST differ and aborting here would
            # abort on the arm working.  The key is therefore a15's throughout: one fixed
            # instrument applied to both packages, which is what makes `n_pred` and rho
            # comparable, and which reproduces 202's own +0.9993 on the incumbent as this
            # probe's known-answer row.  rho for the ARM is a CROSS-key number -- the order it
            # delivers read against the f = 0.10 key -- and is reported as such, never as
            # "the delivered order fell".
            print(f"  {label:<12} ^^ EXPECTED for the dose arm: the |H| key below is a15's "
                  f"(f = 0.10) for BOTH packages, by design")

    genes = np.asarray(
        [line.strip() for line in (BUNDLE / "gene_names.csv").read_text().splitlines()[1:]]
    )
    perts = np.asarray(
        [line.split(",")[0] for line in
         (BUNDLE / "pert_counts.csv").read_text().splitlines()[1:]]
    )
    assert genes.size == effects.shape[1] and perts.size == effects.shape[0]

    handles = {}
    obs = {}
    for label, path in packages.items():
        adata = ad.read_h5ad(path, backed="r")
        frame = adata.obs
        obs[label] = (
            np.asarray(frame["target_gene"].astype(str).values),
            np.asarray(frame["context"].astype(str).values),
        )
        adata.file.close()
        handles[label] = h5py.File(path, "r")

    out: dict[str, dict] = {label: {} for label in packages}
    lfc_store: dict[tuple[str, str, str], np.ndarray] = {}
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
        ref_lin = np.expm1(ctrl[:, keep]).mean(0)
        del ctrl
        # the builder's own gate, and the key it hands the emitter.  `b_ctrl` is the POOLED
        # control profile in the scorer's space, exactly as `write_submission_from_effects`
        # builds it from the same file's column sums.
        b_ctrl = np.log1p(TS * pooled / pooled.sum())
        expression = np.expm1(b_ctrl)
        gate_cpm = 1e6 * expression / max(float(expression.sum()), 1e-12)
        builder_gate = (gate_cpm > BUILDER_MIN_CPM) & (b_ctrl > 0)
        key_all = np.where(builder_gate[None, :], effects, 0.0)
        print(f"\nctx {context}: {cs.shape[1]} controls, {int(keep.sum())} genes past the "
              f"{SCORER_MIN_CPM}-CPM scorer gate, {int(builder_gate.sum())} past the "
              f"builder's {BUILDER_MIN_CPM}-CPM gate  [{time.time() - t0:.0f}s]", flush=True)
        del pooled

        chosen = sorted(perts.tolist())
        picks = np.random.default_rng(191).choice(
            len(chosen), min(n_perts, len(chosen)), replace=False
        )
        targets = [chosen[i] for i in sorted(picks)]
        pos = {str(p): i for i, p in enumerate(perts)}

        for label, path in packages.items():
            labels, contexts = obs[label]
            rows = []
            handle = handles[label]
            data, indices, indptr = handle["X"]["data"], handle["X"]["indices"], \
                handle["X"]["indptr"]
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
                p = 2 * norm.sf(np.abs(z))
                q = false_discovery_control(p, method="bh")
                lfc = np.log2((np.expm1(normed[:, keep]).mean(0) + 1e-9) / (ref_lin + 1e-9))
                del normed
                sig = q < ALPHA
                # the target's own row is excluded from every member (BRIEF §4-7)
                own = np.flatnonzero(genes[keep] == target)
                if own.size:
                    sig[own[0]] = False
                key = np.abs(key_all[pos[target]])[keep]
                expr = np.log1p(cpm[keep])
                lp = -np.log10(np.maximum(p, 1e-300))
                rank = np.argsort(np.argsort(q))
                dec = (rank * 10 // rank.size).clip(0, 9)
                rows.append(dict(
                    target=target,
                    n_pred=int(sig.sum()),
                    rho_H_call=(float(spearmanr(lp[sig], key[sig]).correlation)
                                if sig.sum() > 10 else np.nan),
                    rho_E_call=(float(spearmanr(lp[sig], expr[sig]).correlation)
                                if sig.sum() > 10 else np.nan),
                    rho_H_univ=float(spearmanr(lp, key).correlation),
                    frac_up_call=(float((lfc[sig] > 0).mean()) if sig.any() else np.nan),
                    up_by_dec=[float((lfc[dec == k] > 0).mean()) for k in range(10)],
                ))
                lfc_store[(label, context, target)] = lfc.astype(np.float32)
            npred = np.array([r["n_pred"] for r in rows])
            out[label][context] = dict(
                n=len(rows), targets=targets, n_pred=npred.tolist(),
                n_pred_median=float(np.median(npred)),
                rho_H_call=float(np.nanmedian([r["rho_H_call"] for r in rows])),
                rho_E_call=float(np.nanmedian([r["rho_E_call"] for r in rows])),
                rho_H_univ=float(np.nanmedian([r["rho_H_univ"] for r in rows])),
                frac_up_call=float(np.nanmean([r["frac_up_call"] for r in rows])),
                up_by_dec=np.array([r["up_by_dec"] for r in rows]).mean(0).tolist(),
            )
            row = out[label][context]
            print(f"  {label:<12} n_pred median {row['n_pred_median']:.0f} "
                  f"[{npred.min()}, {npred.max()}]   rho(-log p, |H|) in call set "
                  f"{row['rho_H_call']:+.3f}   rho(expr) {row['rho_E_call']:+.3f}   "
                  f"call set {row['frac_up_call'] * 100:.1f}% UP", flush=True)
        del cs, tie, ref_lin, key_all

    for handle in handles.values():
        handle.close()

    print("\n=== FALSIFIER 3 — challenge `n_pred`, pooled over the probed blocks ===")
    pooled = {}
    for label in packages:
        values = np.concatenate([out[label][c]["n_pred"] for c in ("A", "B", "C")])
        pooled[label] = float(np.median(values))
        print(f"  {label:<12} pooled median {pooled[label]:.0f}   "
              + "  ".join(f"{c}={out[label][c]['n_pred_median']:.0f}" for c in "ABC"))
    labels = list(packages)
    if len(labels) == 2:
        base, arm = ("a15_local", "rc_hybrid") if "a15_local" in labels else labels[:2]
        change = (pooled[arm] - pooled[base]) / pooled[base]
        print(f"  {arm} vs {base}: {change * 100:+.2f}%   "
              f"(falsifier 3 bound +-2%)  {'PASS' if abs(change) <= 0.02 else 'FAIL'}")

    print("\n=== FALSIFIER 4 — delivered order on the REAL contexts ===")
    for label in packages:
        for context in ("A", "B", "C"):
            row = out[label][context]
            print(f"  {label:<12} ctx {context}  rho(-log p, |H|) call set "
                  f"{row['rho_H_call']:+.4f}   universe {row['rho_H_univ']:+.4f}   "
                  f"expression {row['rho_E_call']:+.4f}")
        allrho = [out[label][c]["rho_H_call"] for c in "ABC"]
        print(f"  {label:<12} MEAN over contexts {np.mean(allrho):+.4f}")

    np.savez_compressed(HERE / "work" / "probe217_lfc.npz",
                        **{f"{k[0]}|{k[1]}|{k[2]}": v for k, v in lfc_store.items()})
    with open(HERE / "work" / "probe217.json", "w") as handle:
        json.dump(out, handle, indent=2)
    print(f"\nwrote {HERE / 'work' / 'probe217.json'}")


if __name__ == "__main__":
    main()
