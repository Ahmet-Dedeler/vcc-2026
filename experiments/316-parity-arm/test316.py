"""316 -- the three known-answer rows for `VCC_NORMREL_PM=blend`, on 313's saved fixture.

  PYTHONDONTWRITEBYTECODE=1 uv run python experiments/316-parity-arm/test316.py

Rows 1 and 2 are the FIXTURE versions of the mandate's byte-for-byte rows: the patched builder is
run against a byte copy of the builder as it stood before the patch (`work/builder.PRE316.py`), on
313's fixture, through the same two real functions a build calls.

  KA-1F  VCC_NORMREL_PM unset, exponent -1.5  ==  the PRE-316 builder's arm, BITWISE   (-> submissions/312-normrel15)
  KA-2F  VCC_NORMREL_PM=blend, exponent 0     ==  the PRE-316 builder's INCUMBENT, BITWISE (-> cd4svd5)
  KA-3   VCC_NORMREL_PM=blend, exponent -1.5  ->  delivered row norm / incumbent == 1 for EVERY target

plus the comparison of the built arm against 313's own `work/parity313.json`.
Read-only outside `experiments/316-parity-arm/`.  No bank, no build, no submission.
"""
from __future__ import annotations

import json
import resource
import sys
import types
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
E313 = ROOT / "experiments/313-normrel-parity"

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(E313))

import parity313 as P313  # noqa: E402  (helpers only; its main() does not run on import)

EXP = -1.5
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def load_pre316():
    """The builder as it stood BEFORE the patch, as an independent module.

    `__file__` is forced to the real script path because the module computes `ROOT` from it.
    """
    name = "builder_pre316"
    mod = types.ModuleType(name)
    mod.__file__ = str(ROOT / "scripts/make_submission_pm_cd4.py")
    src = (WORK / "builder.PRE316.py").read_text()
    sys.modules[name] = mod
    exec(compile(src, mod.__file__, "exec"), mod.__dict__)
    return mod


def main() -> None:
    fx = P313.load_fixture()
    targets = [str(t) for t in fx["targets"]]
    official = [str(g) for g in fx["official"]]
    panel = np.asarray(fx["panel_positions"], dtype=np.int64)
    n, G = len(targets), len(official)

    model = P313.build_model(fx)
    pm = P313.build_pm(fx, official)

    import make_submission_pm_cd4 as MB  # noqa: E402  (the PATCHED builder)

    OLD = load_pre316()
    ALPHA, BETA = float(MB.ALPHA), float(MB.BETA)
    uni = MB.donor_weights(model)
    # production's operator columns = the union of the donors' official coverage MINUS the panel
    op_cols = MB.production_operator_columns(model, [official[i] for i in panel])
    wcols = np.asarray([c for c in range(G) if c not in set(panel.tolist())], dtype=np.int64)

    P(f"[316] fixture {P313.sha(E313 / 'work/fixture313.npz')[:16]}  "
      f"{n} targets x {G} genes x {len(model.names)} donors; panel {panel.size}, "
      f"operator columns {op_cols.size}, non-panel {wcols.size}")
    P(f"[316] builder PATCHED {P313.sha(ROOT / 'scripts/make_submission_pm_cd4.py')[:16]}  "
      f"PRE-316 {P313.sha(WORK / 'builder.PRE316.py')[:16]}")
    P(f"[316] operator columns == non-panel columns on this fixture: "
      f"{bool(np.array_equal(np.sort(op_cols), wcols))}")

    KA: list[dict] = []

    def ka(label, got, tol, want=0.0):
        d = abs(float(got) - float(want))
        ok = bool(d <= tol)
        KA.append(dict(label=label, value=float(got), abs_diff=d, tol=tol, ok=ok))
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.6e} (tol {tol:.0e})")
        return ok

    def full(mod, exponent, pmmode):
        """One complete pre-whiten pass through the real functions: the arm, then PerturbMatch."""
        mod.NORMREL_EXP = exponent
        if hasattr(mod, "NORMREL_PM"):
            mod.NORMREL_PM = pmmode
        handled: set[str] = set()
        kw = {}
        if mod is MB:
            kw = dict(pm=pm, operator_columns=op_cols, pm_handled_out=handled)
        pre = mod.build_v9_pre_whiten_effects(
            model, targets, alpha=ALPHA, beta=BETA, weights=uni, **kw)
        akw = dict(already_blended=handled) if mod is MB else {}
        out, audit = mod.apply_perturbmatch_k562_replacement(
            pre, model, pm, targets, alpha=ALPHA, weights=uni, **akw)
        return out, audit, handled, pre

    # --------------------------------------------------------------- KA-1F: unset == 312-normrel15
    new_off, aud_off, h_off, _ = full(MB, EXP, None)
    old_arm, aud_old, _, _ = full(OLD, EXP, None)
    ka("KA-1F  VCC_NORMREL_PM unset, exp -1.5 == the PRE-316 builder's arm (BITWISE)",
       np.abs(new_off - old_arm).max(), 0.0)
    P(f"        targets blended in-place: {len(h_off)} (must be 0 when the knob is unset)")
    ka("KA-1F' the knob unset blends nothing in the pre-whiten path", float(len(h_off)), 0.0)
    same_audit = {k: v for k, v in aud_off.items()} == {k: v for k, v in aud_old.items()}
    P(f"        pm_audit identical to the PRE-316 builder's: {same_audit}")

    # --------------------------------------------------------------- KA-2F: blend + exp 0 == cd4svd5
    new_exp0, aud0, h0, _ = full(MB, 0.0, "blend")
    OLD.NORMREL_EXP = None
    inc_pre = OLD.build_v9_pre_whiten_effects(model, targets, alpha=ALPHA, beta=BETA, weights=uni)
    incumbent, aud_inc = OLD.apply_perturbmatch_k562_replacement(
        inc_pre, model, pm, targets, alpha=ALPHA, weights=uni)
    ka("KA-2F  VCC_NORMREL_PM=blend, exp 0 == the PRE-316 builder's INCUMBENT (BITWISE)",
       np.abs(new_exp0 - incumbent).max(), 0.0)
    P(f"        targets blended in-place at exponent 0: {len(h0)}")
    same_audit0 = ({k: v for k, v in aud0.items()} == {k: v for k, v in aud_inc.items()})
    P(f"        pm_audit identical to the incumbent's: {same_audit0}")
    ka("KA-2F' the audit dict is unchanged by the in-place blend", 0.0 if same_audit0 else 1.0, 0.0)

    # --------------------------------------------------------------- KA-3: exactly shape-only
    arm, aud_arm, h_arm, arm_pre = full(MB, EXP, "blend")
    na = np.linalg.norm(arm, axis=1)
    ni = np.linalg.norm(incumbent, axis=1)
    ratio = na / np.maximum(ni, 1e-300)
    naw = np.linalg.norm(arm[:, wcols], axis=1)
    niw = np.linalg.norm(incumbent[:, wcols], axis=1)
    ratio_w = naw / np.maximum(niw, 1e-300)
    ka("KA-3   delivered row norm / incumbent == 1 for EVERY target (full axis, relative)",
       float(np.abs(ratio - 1.0).max()), 1e-6)
    P(f"        full axis : median {np.median(ratio):.7f}  range "
      f"[{ratio.min():.7f}, {ratio.max():.7f}]  max|ratio-1| {np.abs(ratio - 1).max():.3e}")
    P(f"        non-panel : median {np.median(ratio_w):.7f}  range "
      f"[{ratio_w.min():.7f}, {ratio_w.max():.7f}]   (the panel is {100 * panel.size / G:.1f} % of "
      f"this fixture's columns vs production's 1.62 %; 313 section 5)")
    P(f"        whole-matrix norm ratio {np.linalg.norm(arm) / np.linalg.norm(incumbent):.7f}")
    P(f"        targets blended in-place: {len(h_arm)}/{n}")

    # --------------------------------------------------------------- the shipped arm, for contrast
    shipped = old_arm
    ns = np.linalg.norm(shipped, axis=1)
    nsw = np.linalg.norm(shipped[:, wcols], axis=1)
    rs, rsw = ns / np.maximum(ni, 1e-300), nsw / np.maximum(niw, 1e-300)
    P("")
    P(f"[316] SHIPPED arm (312, for contrast): full axis median {np.median(rs):.6f} range "
      f"[{rs.min():.6f}, {rs.max():.6f}]; non-panel median {np.median(rsw):.6f} range "
      f"[{rsw.min():.6f}, {rsw.max():.6f}]; whole matrix "
      f"{np.linalg.norm(shipped) / np.linalg.norm(incumbent):.6f}")

    def rowcos(A, B, cols):
        na_, nb_ = np.linalg.norm(A[:, cols], axis=1), np.linalg.norm(B[:, cols], axis=1)
        return np.array([float(A[t, cols] @ B[t, cols] / max(na_[t] * nb_[t], 1e-300))
                         for t in range(n)])

    d_vs_shipped = float(np.abs(arm[:, wcols] - shipped[:, wcols]).max())
    cos_vs_shipped = rowcos(arm, shipped, wcols)
    d_vs_incumbent = float(np.abs(arm[:, wcols] - incumbent[:, wcols]).max())
    P(f"[316] PARITY arm vs SHIPPED arm, non-panel columns: max|diff| {d_vs_shipped:.4e}, "
      f"row cosine median {np.median(cos_vs_shipped):.6f} min {cos_vs_shipped.min():.6f}")
    P(f"[316] PARITY arm vs the INCUMBENT (cd4svd5), non-panel: max|diff| {d_vs_incumbent:.4e} "
      f"-- the arm is a real re-shaping, not a cleanup")

    # 313's own parity variant, for the fixture comparison (it froze D5 the OTHER way: the norm was
    # taken over the FULL axis, and on the PM-blended rows).  Recomputed here from parity313's own
    # clean-room code so the two numbers are directly comparable.
    ref = json.loads((E313 / "work/parity313.json").read_text())
    rows_prod, avail = P313.harness_rows(model, targets, np.asarray([], dtype=np.int64))
    pm_row_prod, pm_avail = P313.harness_pm_row(model, pm, targets, np.asarray([], dtype=np.int64))
    k_index = model.names.index(MB.K562)
    BLEND = float(MB.PM_BLEND)
    pm_cols = np.asarray(pm["target_columns"], dtype=np.int64)
    both = avail[k_index] & pm_avail
    br = np.flatnonzero(both)
    rows_pb = [R.copy() for R in rows_prod]
    rows_pb[k_index][np.ix_(br, pm_cols)] = (
        (1 - BLEND) * rows_prod[k_index][np.ix_(br, pm_cols)] + BLEND * pm_row_prod[np.ix_(br, pm_cols)])

    def uniform_mean(rows):
        acc = np.zeros((n, G))
        cnt = avail.sum(0)
        for j in range(len(rows)):
            acc[avail[j]] += rows[j][avail[j]]
        cov = cnt > 0
        acc[cov] /= cnt[cov, None]
        return acc

    pre_pb = uniform_mean(rows_pb)
    W_313 = P313.normrel_W_production_rule_rows(rows_pb, avail, np.arange(G), EXP)
    parity313_variant = ALPHA * P313.ord_A(rows_pb, avail, W_313, pre_pb)
    d313 = float(np.abs(parity313_variant[:, wcols] - shipped[:, wcols]).max())
    c313 = rowcos(parity313_variant, shipped, wcols)
    P("")
    P(f"[316] 313's own parity variant recomputed here (D5 = FULL axis, weights on the BLENDED "
      f"rows): max|diff| vs shipped {d313:.4e} (313 reported "
      f"{ref['parity_variant']['max_abs_diff_vs_shipped']:.4e}), row cosine median "
      f"{np.median(c313):.6f} (313: {ref['parity_variant']['median_row_cosine_vs_shipped']:.6f}) "
      f"min {c313.min():.6f} (313: {ref['parity_variant']['min_row_cosine_vs_shipped']:.6f})")
    ka("KA-4F  313's parity_variant reproduced from its own clean-room code",
       abs(d313 - float(ref["parity_variant"]["max_abs_diff_vs_shipped"])), 1e-9)
    d_316_vs_313 = float(np.abs(arm[:, wcols] - parity313_variant[:, wcols]).max())
    P(f"[316] 316's arm vs 313's variant (the D5 + D3 freeze difference alone): "
      f"max|diff| {d_316_vs_313:.4e}, row cosine median "
      f"{np.median(rowcos(arm, parity313_variant, wcols)):.6f}")

    out = dict(
        fixture=dict(path="experiments/313-normrel-parity/work/fixture313.npz",
                     sha256=P313.sha(E313 / "work/fixture313.npz"), n_targets=n, n_genes=G),
        builder=dict(patched=P313.sha(ROOT / "scripts/make_submission_pm_cd4.py"),
                     pre316=P313.sha(WORK / "builder.PRE316.py")),
        known_answers=KA,
        delivered_norm_parity=dict(
            full_axis=dict(median=float(np.median(ratio)), min=float(ratio.min()),
                           max=float(ratio.max()), max_abs_dev=float(np.abs(ratio - 1).max())),
            non_panel=dict(median=float(np.median(ratio_w)), min=float(ratio_w.min()),
                           max=float(ratio_w.max())),
            whole_matrix=float(np.linalg.norm(arm) / np.linalg.norm(incumbent)),
            per_target=[float(x) for x in ratio]),
        delivered_norm_shipped=dict(
            full_axis=dict(median=float(np.median(rs)), min=float(rs.min()), max=float(rs.max())),
            non_panel=dict(median=float(np.median(rsw)), min=float(rsw.min()), max=float(rsw.max())),
            whole_matrix=float(np.linalg.norm(shipped) / np.linalg.norm(incumbent))),
        arm_vs_shipped=dict(max_abs_diff=d_vs_shipped,
                            median_row_cosine=float(np.median(cos_vs_shipped)),
                            min_row_cosine=float(cos_vs_shipped.min())),
        arm_vs_incumbent=dict(max_abs_diff=d_vs_incumbent),
        parity313_variant=dict(max_abs_diff_vs_shipped=d313,
                               reported_by_313=float(ref["parity_variant"]["max_abs_diff_vs_shipped"]),
                               median_row_cosine=float(np.median(c313)),
                               reported_median_313=float(ref["parity_variant"]["median_row_cosine_vs_shipped"]),
                               min_row_cosine=float(c313.min()),
                               reported_min_313=float(ref["parity_variant"]["min_row_cosine_vs_shipped"]),
                               max_abs_diff_vs_316=d_316_vs_313),
        blended_targets=sorted(h_arm),
        peak_gib=peak_gib(),
    )
    (WORK / "test316.json").write_text(json.dumps(out, indent=2))
    (WORK / "test316.log").write_text("\n".join(LINES) + "\n")
    P("")
    P(f"[316] KA rows {sum(1 for k in KA if k['ok'])}/{len(KA)} pass; peak {peak_gib():.2f} GiB")
    if not all(k["ok"] for k in KA):
        raise SystemExit("KNOWN-ANSWER FAILURE")


if __name__ == "__main__":
    main()
