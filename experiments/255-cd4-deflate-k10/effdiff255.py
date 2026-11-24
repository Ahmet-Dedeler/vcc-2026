"""255 falsifiers 2 and 3 at the EFFECT-MATRIX level (249 `effdiff249.py`'s exact form, incumbent = cd4svd5).

`effects-svd10.npy` vs `effects-svd5.npy` (the two staged 300 x 18,533 matrices, SAME builder bytes,
one flag value).  Rows CD4 carries (nonzero support, 276) must move in DIRECTION; the rows it does
not (24) must equal the incumbent's times ONE shared scalar -- the merged preserve ratio -- to
machine precision.  Non-vacuity: `||dE||/||E||` and the number of rows with exactly zero
displacement.  Also re-derives the k = 10 deflation audit offline from the saved CD4 residual
through the builder's function and asserts it equals what the staged build recorded, and (for the
record) the same displacement read against 249's staged `h1union` matrix `effects-off.npy`
(array sha `abeb0749...`, verified before use).

    uv run python effdiff255.py
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
S = WORK / "stage"
S249 = ROOT / "experiments/249-cd4-deflate/work/stage"
OFF_WANT = "abeb0749f8ecc9410fcec0d9d0cebd2569c97ca456d73f9e6b3bfd359164c315"
K = 10


def asha(a):
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def split_stats(on, off, carry, s):
    n_off = np.linalg.norm(off, axis=1)
    n_on = np.linalg.norm(on, axis=1)
    d = np.linalg.norm(on - off, axis=1)
    rel = d / n_off
    cos = (on * off).sum(1) / (n_on * n_off)
    corrected = np.linalg.norm(on - s * off, axis=1) / n_off
    ratio = n_on / n_off
    return n_off, n_on, d, rel, cos, corrected, ratio


def main() -> None:
    inc = np.load(S / "effects-svd5.npy")
    on = np.load(S / "effects-svd10.npy")
    cd4 = np.load(S / "cd4-svd10.npz")
    j_inc = json.loads((S / "stage-svd5.json").read_text())
    j_on = json.loads((S / "stage-svd10.json").read_text())
    assert asha(inc) == j_inc["effects_array_sha256"] and asha(on) == j_on["effects_array_sha256"]
    votes = np.asarray(cd4["gene_votes"])
    counts = np.asarray(cd4["target_condition_count"])
    carry = votes.sum(1) > 0
    targets = np.asarray(cd4["targets"]).astype(str)
    s_inc = j_inc["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]
    s_on = j_on["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]
    s = s_on / s_inc
    lines = []

    def P(t=""):
        print(t)
        lines.append(t)

    P(f"effects: cd4svd5 {j_inc['effects_array_sha256'][:16]} (VCC_CD4_DEFLATE_K=5)  cd4svd10 {j_on['effects_array_sha256'][:16]} "
      f"(=10)  shape {on.shape}  same builder bytes: {j_inc['builder_sha256'] == j_on['builder_sha256']} ({j_on['builder_sha256'][:8]})")
    P(f"CD4 rows with support: {int(carry.sum())}/300 (count==3: {int((counts == 3).sum())}; partial count 1-2: "
      f"{int(((counts > 0) & (counts < 3)).sum())}; count 0: {int((counts == 0).sum())})")
    P(f"merged preserve scalar: cd4svd5 {s_inc!r}  cd4svd10 {s_on!r}  ratio s = {s!r}")

    n_off, n_on, d, rel, cos, corrected, ratio = split_stats(on, inc, carry, s)
    P(f"\n{'row class':16s}{'n':>5s}{'rel displ min':>15s}{'median':>10s}{'max':>10s}{'cos min':>14s}{'cos max':>14s}"
      f"{'||dE - s E||/||E|| max':>24s}{'median':>12s}")
    for tag, sel in (("CD4-carried", carry), ("NOT carried", ~carry)):
        P(f"{tag:16s}{int(sel.sum()):5d}{rel[sel].min():15.6f}{np.median(rel[sel]):10.6f}{rel[sel].max():10.6f}"
          f"{cos[sel].min():14.10f}{cos[sel].max():14.10f}{corrected[sel].max():24.3e}{np.median(corrected[sel]):12.3e}")
    P(f"NOT-carried norm ratio: min {ratio[~carry].min()!r} max {ratio[~carry].max()!r} (s = {s!r})")
    sep = corrected[carry].min() / max(corrected[~carry].max(), 1e-300)
    P(f"carried-vs-non-carried separation on the scale-corrected displacement: {sep:.1e}x "
      f"(min carried {corrected[carry].min():.6f} vs max non-carried {corrected[~carry].max():.2e})")
    f2 = bool(corrected[~carry].max() < 1e-12 and cos[carry].max() < 1 - 1e-9 and cos[~carry].min() > 1 - 1e-12)
    P(f"FALSIFIER 2 (effect level) {'PASS' if f2 else 'FAIL'}: carried rows move in direction "
      f"(cos <= {cos[carry].max():.6f}), non-carried rows are cd4svd5's x one scalar to {corrected[~carry].max():.1e}")
    rel_global = float(np.linalg.norm(on - inc) / np.linalg.norm(inc))
    zero_rows = int(np.sum(d == 0.0))
    P(f"\nNON-VACUITY (effect level, vs cd4svd5): ||dE||/||E|| = {rel_global:.6f}; per-row min {rel.min():.6f}; rows with "
      f"displacement exactly 0: {zero_rows}/300; rows < 1e-6: {int(np.sum(rel < 1e-6))}/300")
    f3 = zero_rows == 0 and rel_global > 0.01
    P(f"FALSIFIER 3 (effect level) {'PASS' if f3 else 'FAIL'}")
    P(f"\nwhere the displacement lives: carried rows carry {100 * float((d[carry] ** 2).sum() / (d ** 2).sum()):.3f}% "
      f"of ||dE||^2; top-5 rows by displacement: "
      + ", ".join(f"{targets[i]} {rel[i]:.3f}" for i in np.argsort(-rel)[:5]))

    # ---- for the record: the same read against h1union (249's staged `off`, verified by array sha)
    rec_off = None
    off_path = S249 / "effects-off.npy"
    if off_path.exists():
        off = np.load(off_path)
        got = asha(off)
        if got == OFF_WANT:
            j_off = json.loads((S249 / "stage-off.json").read_text())
            s_off = j_off["audit"]["cd4"]["post_add_preserve"]["post_add_global_preserve_scale"]
            s2 = s_on / s_off
            _, _, d2, rel2, cos2, corr2, ratio2 = split_stats(on, off, carry, s2)
            rel_g2 = float(np.linalg.norm(on - off) / np.linalg.norm(off))
            P(f"\nvs h1union (249 effects-off.npy, array sha {got[:16]} == abeb0749 VERIFIED): preserve ratio s = {s2!r}; "
              f"||dE||/||E|| {rel_g2:.6f}; carried rel displ min/median/max {rel2[carry].min():.6f}/{np.median(rel2[carry]):.6f}/{rel2[carry].max():.6f}, "
              f"cos min {cos2[carry].min():.7f}; non-carried scale-corrected max {corr2[~carry].max():.2e}, norm ratio {ratio2[~carry].min()!r}..{ratio2[~carry].max()!r}; "
              f"rows at zero {int(np.sum(d2 == 0.0))}/300")
            rec_off = dict(shared_scalar=s2, rel_global=rel_g2, carried_rel=[float(rel2[carry].min()), float(np.median(rel2[carry])), float(rel2[carry].max())],
                           carried_cos_min=float(cos2[carry].min()), not_carried_corrected_max=float(corr2[~carry].max()),
                           zero_rows=int(np.sum(d2 == 0.0)))
        else:
            P(f"\n249's effects-off.npy hashes to {got[:16]}, not abeb0749 -- NOT used")

    # ---- the deflation audit, re-derived offline through the builder from the saved residual
    spec = importlib.util.spec_from_file_location("effdiff255_builder", ROOT / "scripts/make_submission_pm_cd4.py")
    M = importlib.util.module_from_spec(spec)
    sys.modules["effdiff255_builder"] = M
    spec.loader.exec_module(M)
    res = M.EXP062.CD4Prepared(mode="gene_norm", mean_effect=np.asarray(cd4["mean_effect"]), gene_votes=votes,
                               target_condition_count=counts, condition_details=[])
    _, aud = M.deflate_cd4_residual(res, K)
    rec = j_on["audit"]["cd4"]["deflation"]
    same = all(aud[k] == rec[k] for k in aud)
    P(f"\ndeflation audit (k={K}) re-derived offline from the saved CD4 residual == staged build's record: {same}")
    P("  " + json.dumps(rec))
    _, aud5 = M.deflate_cd4_residual(res, 5)
    rec5 = j_inc["audit"]["cd4"]["deflation"]
    same5 = all(aud5[k] == rec5[k] for k in aud5)
    P(f"  and k=5 re-derived == cd4svd5's staged record (249's f06e0262...): {same5}  U_5 sha {rec5['basis_u_k_sha256'][:16]}")
    R = np.asarray(cd4["mean_effect"])[counts > 0]
    P(f"  centred energy share of the top-{K} directions {rec['top_k_energy_frac_of_centred']:.4f}; energy removed "
      f"{rec['energy_removed_frac']:.4f} of the residual's {float((R * R).sum()):.4e} (k=5: {rec5['energy_removed_frac']:.4f})")
    # the top-10 subspace CONTAINS the top-5 one (eigh ordering): the incremental 5 directions' energy
    P(f"  incremental energy of directions 6-10: {rec['energy_removed_frac'] - rec5['energy_removed_frac']:.4f} of the residual")
    out = dict(shared_scalar=s, preserve=[s_inc, s_on], n_carried_rows=int(carry.sum()),
               carried=dict(rel_min=float(rel[carry].min()), rel_median=float(np.median(rel[carry])),
                            rel_max=float(rel[carry].max()), cos_min=float(cos[carry].min()),
                            cos_max=float(cos[carry].max()), corrected_min=float(corrected[carry].min()),
                            corrected_median=float(np.median(corrected[carry]))),
               not_carried=dict(rel_min=float(rel[~carry].min()), rel_max=float(rel[~carry].max()),
                                cos_min=float(cos[~carry].min()), corrected_max=float(corrected[~carry].max()),
                                ratio_min=float(ratio[~carry].min()), ratio_max=float(ratio[~carry].max())),
               separation=float(sep), rel_global=rel_global, zero_rows=zero_rows,
               falsifier2_pass=f2, falsifier3_pass=f3, deflation_audit_reproduced=bool(same),
               deflation5_audit_reproduced=bool(same5), deflation=rec, deflation5=rec5, vs_h1union=rec_off)
    (WORK / "effdiff255.json").write_text(json.dumps(out, indent=1))
    (WORK / "effdiff255.log").write_text("\n".join(lines) + "\n")
    if not (f2 and f3 and same and same5):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
