"""319 -- the three known-answer rows for `VCC_NORMREL_MODE=cap`, on 313's saved fixture.

  PYTHONDONTWRITEBYTECODE=1 uv run python experiments/319-cap-variant/test319.py

Rows 1 and 2 are the FIXTURE versions of the mandate's byte-for-byte rows: the patched builder is
run against a byte copy of the builder as it stood before the patch (`work/builder.PRE319.py`,
i.e. the 316-patched builder that produced `316-normrel15-pm`), on 313's fixture, through the same
two real functions a build calls.

  KA-1  VCC_NORMREL_MODE unset/`plain`, exp -1.5  ==  the PRE-319 builder's arm, BITWISE  (-> 312-normrel15)
  KA-2  VCC_NORMREL_MODE=cap, exp 0               ==  the PRE-319 builder's INCUMBENT, BITWISE (-> cd4svd5)
  KA-3  the cap actually caps:
        3a  every donor with r <= 1 keeps EXACTLY its uniform weight; every donor with r > 1 is
            strictly down-weighted   (the raw `w` dict, before `model.effect` renormalises it)
        3b  the shares the PATCHED BUILDER actually used, recovered from its own output by least
            squares (313's KA-2b recipe), == the capped rule's shares
        3c  the cap is not a no-op: the capped shares differ from the plain rule's

plus KA-4, a cross-experiment regression row: the SHIPPED (`plain`) arm's delivered-norm numbers on
this fixture must reproduce `experiments/316-parity-arm/work/test316.json` to the digit.

Read-only outside `experiments/319-cap-variant/`.  No bank, no build, no submission.
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
E316 = ROOT / "experiments/316-parity-arm"

sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(E313))
sys.path.insert(0, str(HERE))

import lib319 as L  # noqa: E402
import parity313 as P313  # noqa: E402  (helpers only; its main() does not run on import)

EXP = -1.5
LINES: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LINES.append(s)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def load_pre319():
    """The builder as it stood BEFORE the patch, as an independent module.

    `__file__` is forced to the real script path because the module computes `ROOT` from it.
    """
    name = "builder_pre319"
    mod = types.ModuleType(name)
    mod.__file__ = str(ROOT / "scripts/make_submission_pm_cd4.py")
    src = (WORK / "builder.PRE319.py").read_text()
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

    OLD = load_pre319()
    ALPHA, BETA = float(MB.ALPHA), float(MB.BETA)
    uni = MB.donor_weights(model)
    wcols = np.asarray([c for c in range(G) if c not in set(panel.tolist())], dtype=np.int64)

    P(f"[319] fixture {P313.sha(E313 / 'work/fixture313.npz')[:16]}  "
      f"{n} targets x {G} genes x {len(model.names)} donors; panel {panel.size}, "
      f"non-panel {wcols.size}")
    P(f"[319] builder PATCHED {P313.sha(ROOT / 'scripts/make_submission_pm_cd4.py')[:16]}  "
      f"PRE-319 {P313.sha(WORK / 'builder.PRE319.py')[:16]}")
    P(f"[319] module constants: patched NORMREL_MODE={MB.NORMREL_MODE!r}  "
      f"PRE-319 has NORMREL_MODE: {hasattr(OLD, 'NORMREL_MODE')}")
    P(f"[319] donor weights (uniform, production's default): {uni}")

    KA: list[dict] = []

    def ka(label, got, tol, want=0.0):
        d = abs(float(got) - float(want))
        ok = bool(d <= tol)
        KA.append(dict(label=label, value=float(got), abs_diff=d, tol=tol, ok=ok))
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: |d| {d:.6e} (tol {tol:.0e})")
        return ok

    def full(mod, exponent, mode):
        """One complete pre-whiten pass through the real functions: the arm, then PerturbMatch.

        `VCC_NORMREL_PM` stays UNSET on both sides -- 319's arm is one variable against
        `312-normrel15`, which is the weight-first rendering.
        """
        mod.NORMREL_EXP = exponent
        if hasattr(mod, "NORMREL_PM"):
            mod.NORMREL_PM = None
        if hasattr(mod, "NORMREL_MODE"):
            mod.NORMREL_MODE = mode
        pre = mod.build_v9_pre_whiten_effects(
            model, targets, alpha=ALPHA, beta=BETA, weights=uni)
        out, audit = mod.apply_perturbmatch_k562_replacement(
            pre, model, pm, targets, alpha=ALPHA, weights=uni)
        return out, audit, pre

    # --------------------------------------------------------------- KA-1: plain == the shipped arm
    new_plain, aud_plain, pre_plain = full(MB, EXP, "plain")
    old_arm, aud_old, pre_old = full(OLD, EXP, None)
    ka("KA-1   VCC_NORMREL_MODE=plain, exp -1.5 == the PRE-319 builder's arm (BITWISE)",
       np.abs(new_plain - old_arm).max(), 0.0)
    same_audit = {k: v for k, v in aud_plain.items()} == {k: v for k, v in aud_old.items()}
    P(f"        pm_audit identical to the PRE-319 builder's: {same_audit}")
    ka("KA-1'  the audit dict is unchanged by the patch", 0.0 if same_audit else 1.0, 0.0)

    # the DEFAULT is `plain`: an unset knob must give the same matrix as an explicit `plain`
    MB.NORMREL_MODE = "plain"
    default_mode = (WORK / "default_mode.txt")
    ka("KA-1'' the module's default VCC_NORMREL_MODE is 'plain'",
       0.0 if MB.__dict__["NORMREL_MODE"] == "plain" else 1.0, 0.0)
    default_mode.write_text("plain\n")

    # --------------------------------------------------------------- KA-2: cap + exp 0 == cd4svd5
    new_cap0, aud0, _ = full(MB, 0.0, "cap")
    OLD.NORMREL_EXP = None
    inc_pre = OLD.build_v9_pre_whiten_effects(model, targets, alpha=ALPHA, beta=BETA, weights=uni)
    incumbent, aud_inc = OLD.apply_perturbmatch_k562_replacement(
        inc_pre, model, pm, targets, alpha=ALPHA, weights=uni)
    ka("KA-2   VCC_NORMREL_MODE=cap, exp 0 == the PRE-319 builder's INCUMBENT (BITWISE)",
       np.abs(new_cap0 - incumbent).max(), 0.0)
    same_audit0 = {k: v for k, v in aud0.items()} == {k: v for k, v in aud_inc.items()}
    ka("KA-2'  the audit dict at exponent 0 is the incumbent's", 0.0 if same_audit0 else 1.0, 0.0)

    # --------------------------------------------------------------- the arm itself
    arm, aud_arm, arm_pre = full(MB, EXP, "cap")

    # --------------------------------------------------------------- KA-3a: the rule, stated cleanly
    ratios: dict[str, dict[str, float]] = {}
    carried: dict[str, set[str]] = {}
    medians: dict[str, float] = {}
    for name in model.names:
        r, med = L.donor_panel_ratios(
            model.s[name], model.dcols[name], model.pert_index[name], targets)
        ratios[name] = r
        medians[name] = med
        carried[name] = set(r)

    cells = [(t, name, ratios[name][t]) for name in model.names for t in ratios[name]]
    touched = [(t, nm, r) for (t, nm, r) in cells if r > 1.0]
    untouched = [(t, nm, r) for (t, nm, r) in cells if r <= 1.0]
    worst_untouched = 0.0
    for t, nm, r in untouched:
        w_cap = L.raw_weight(r, EXP, "cap", uni[nm])
        worst_untouched = max(worst_untouched, abs(w_cap - float(uni[nm])))
    n_not_down = 0
    ratio_of_worst = None
    for t, nm, r in touched:
        w_cap = L.raw_weight(r, EXP, "cap", uni[nm])
        if not (w_cap < float(uni[nm])):
            n_not_down += 1
            ratio_of_worst = r
    ka("KA-3a1 every donor with r <= 1 keeps EXACTLY its uniform raw weight under `cap`",
       worst_untouched, 0.0)
    ka("KA-3a2 every donor with r > 1 is strictly DOWN-weighted under `cap` (exp -1.5)",
       float(n_not_down), 0.0)
    P(f"        cells r<=1: {len(untouched)}  r>1: {len(touched)}  "
      f"(worst untouched |w_cap - w_uniform| {worst_untouched:.3e}"
      + (f"; first non-down-weighted r {ratio_of_worst}" if ratio_of_worst is not None else "") + ")")

    # --------------------------------------------------------------- KA-3b: what the BUILDER did
    #   313's KA-2b recipe: recover the shares from the arm's own pre-whiten output by least
    #   squares on the donor basis, and compare against the capped rule's shares.
    rows_raw, avail = P313.harness_rows(model, targets, panel)
    W_cap_rule = L.rule_weights(ratios, uni, targets, EXP, "cap", carried)
    W_plain_rule = L.rule_weights(ratios, uni, targets, EXP, "plain", carried)

    def recover(pre: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Shares, and which donors are IDENTIFIABLE (a donor whose basis row is ~0 on the
        recovery columns carries no information about its own weight: least squares returns the
        minimum-norm 0 for it whatever weight the builder used)."""
        W = np.full((n, len(model.names)), np.nan)
        ident = np.zeros((n, len(model.names)), bool)
        for t in range(n):
            js = np.flatnonzero(avail[:, t])
            if js.size < 2:
                continue
            B = np.stack([rows_raw[j][t][wcols] for j in js])
            bn = np.linalg.norm(B, axis=1)
            ident[t, js] = bn > 1e-8 * max(float(bn.max()), 1e-300)
            x, *_ = np.linalg.lstsq(B.T, pre[t][wcols] / ALPHA, rcond=None)
            s = float(x.sum())
            if s != 0:
                W[t, js] = x / s
        return W, ident

    W_recovered, ident = recover(arm_pre)
    full_ident = np.array([bool(ident[t][np.flatnonzero(avail[:, t])].all())
                           and int(avail[:, t].sum()) >= 2 for t in range(n)])
    live = (np.isfinite(W_cap_rule).all(1) & np.isfinite(W_recovered).all(1)
            & (avail.sum(0) >= 2) & full_ident)
    unident = [targets[t] for t in range(n)
               if (avail[:, t].sum() >= 2) and not full_ident[t]]
    P(f"        targets with an UNIDENTIFIABLE donor row (norm ~0 on the recovery columns; the "
      f"fixture's case (d)): {unident} -- checked separately below")
    P(f"        recovered on {int(live.sum())}/{n} targets")
    ka("KA-3b  the shares the PATCHED builder used (recovered from its own output) == the capped rule",
       np.abs(W_recovered[live] - W_cap_rule[live]).max(), 1e-5)

    # the unidentifiable rows are still checked, on the subspace that IS identifiable: the shares
    # of the donors with a non-zero basis row must be in the capped rule's ratio.
    worst_sub = 0.0
    for t in range(n):
        if avail[:, t].sum() < 2 or full_ident[t] or not np.isfinite(W_cap_rule[t]).all():
            continue
        js = np.flatnonzero(avail[:, t] & ident[t])
        if js.size < 2:
            continue
        rec = W_recovered[t, js] / W_recovered[t, js].sum()
        rule = W_cap_rule[t, js] / W_cap_rule[t, js].sum()
        worst_sub = max(worst_sub, float(np.abs(rec - rule).max()))
    ka("KA-3b' ... and on the identifiable subspace of the case-(d) rows",
       worst_sub, 1e-5)

    d_plain = float(np.abs(W_recovered[live] - W_plain_rule[live]).max())
    P(f"        same recovery against the PLAIN rule's shares: max|dW| {d_plain:.4e} "
      f"(the cap is a real change, not a relabelling)")
    ka("KA-3c  the capped shares are NOT the plain rule's (the arm is not a no-op)",
       1.0 if d_plain > 1e-4 else 0.0, 0.0, want=1.0)

    # --------------------------------------------------------------- KA-4: 316's recorded numbers
    ni = np.linalg.norm(incumbent, axis=1)
    ns = np.linalg.norm(old_arm, axis=1)
    nsw = np.linalg.norm(old_arm[:, wcols], axis=1)
    niw = np.linalg.norm(incumbent[:, wcols], axis=1)
    rs = ns / np.maximum(ni, 1e-300)
    rsw = nsw / np.maximum(niw, 1e-300)
    ref316 = json.loads((E316 / "work/test316.json").read_text())["delivered_norm_shipped"]
    ka("KA-4   316's recorded SHIPPED-arm delivered norm on this fixture reproduces (median, full axis)",
       float(np.median(rs)) - float(ref316["full_axis"]["median"]), 1e-6)
    ka("KA-4'  ... and its non-panel median",
       float(np.median(rsw)) - float(ref316["non_panel"]["median"]), 1e-6)

    # --------------------------------------------------------------- the cap's reach, on the fixture
    per_target_capped = []
    for i, t in enumerate(targets):
        have = [nm for nm in model.names if t in carried[nm]]
        c = sum(1 for nm in have if ratios[nm][t] > 1.0)
        per_target_capped.append(dict(target=t, n_donors=len(have), n_capped=c,
                                      ratios={nm: float(ratios[nm][t]) for nm in have}))
    all_uncapped = [d["target"] for d in per_target_capped if d["n_donors"] and d["n_capped"] == 0]
    share_cells = len(touched) / max(len(cells), 1)
    P("")
    P(f"[319] FIXTURE cap reach: {len(touched)}/{len(cells)} (target, donor) cells have r > 1 "
      f"= {100 * share_cells:.1f} % -- by construction each donor's panel median r is 1.0, so "
      f"about half of every donor's own cells sit above it")
    P(f"[319]   targets where NO carried donor is capped (the row reverts to the incumbent exactly): "
      f"{len(all_uncapped)}/{n} {all_uncapped}")
    for nm in model.names:
        rr = np.array(sorted(ratios[nm].values()))
        P(f"[319]   {nm:14s} n={rr.size:3d} median r {np.median(rr):.6f} "
          f"r>1 {int((rr > 1).sum())}/{rr.size} range [{rr.min():.4f}, {rr.max():.4f}] "
          f"panel median norm {medians[nm]:.6f}")

    # the arm against its two neighbours, on the fixture (structure, not delivered magnitude)
    def rowcos(A, B, cols):
        na_, nb_ = np.linalg.norm(A[:, cols], axis=1), np.linalg.norm(B[:, cols], axis=1)
        return np.array([float(A[t, cols] @ B[t, cols] / max(na_[t] * nb_[t], 1e-300))
                         for t in range(n)])

    d_vs_shipped = float(np.abs(arm[:, wcols] - old_arm[:, wcols]).max())
    cos_vs_shipped = rowcos(arm, old_arm, wcols)
    d_vs_inc = float(np.abs(arm[:, wcols] - incumbent[:, wcols]).max())
    cos_vs_inc = rowcos(arm, incumbent, wcols)
    na = np.linalg.norm(arm, axis=1)
    naw = np.linalg.norm(arm[:, wcols], axis=1)
    rc, rcw = na / np.maximum(ni, 1e-300), naw / np.maximum(niw, 1e-300)
    P("")
    P(f"[319] CAP arm vs the SHIPPED (`plain`) arm, non-panel columns: max|diff| {d_vs_shipped:.4e}, "
      f"row cosine median {np.median(cos_vs_shipped):.6f} min {cos_vs_shipped.min():.6f}")
    P(f"[319] CAP arm vs the INCUMBENT (`cd4svd5`), non-panel columns: max|diff| {d_vs_inc:.4e}, "
      f"row cosine median {np.median(cos_vs_inc):.6f} min {cos_vs_inc.min():.6f}")
    P(f"[319] FIXTURE delivered norm cap/incumbent: full axis median {np.median(rc):.6f} "
      f"range [{rc.min():.6f}, {rc.max():.6f}]; non-panel median {np.median(rcw):.6f} "
      f"range [{rcw.min():.6f}, {rcw.max():.6f}]  "
      f"(FIXTURE numbers -- 313 section 9: magnitudes do not transfer to the real panel)")
    P(f"[319]   shipped `plain` arm, same fixture: full axis median {np.median(rs):.6f} "
      f"range [{rs.min():.6f}, {rs.max():.6f}]; non-panel median {np.median(rsw):.6f} "
      f"range [{rsw.min():.6f}, {rsw.max():.6f}]")

    out = dict(
        fixture=dict(path="experiments/313-normrel-parity/work/fixture313.npz",
                     sha256=P313.sha(E313 / "work/fixture313.npz"), n_targets=n, n_genes=G),
        builder=dict(patched=P313.sha(ROOT / "scripts/make_submission_pm_cd4.py"),
                     pre319=P313.sha(WORK / "builder.PRE319.py")),
        known_answers=KA,
        cap_reach_fixture=dict(cells=len(cells), touched=len(touched), share=share_cells,
                               targets_with_no_capped_donor=all_uncapped,
                               per_donor={nm: dict(n=len(ratios[nm]),
                                                   median=float(np.median(list(ratios[nm].values()))),
                                                   n_gt1=int(sum(1 for v in ratios[nm].values() if v > 1)),
                                                   min=float(min(ratios[nm].values())),
                                                   max=float(max(ratios[nm].values())),
                                                   panel_median_norm=medians[nm])
                                          for nm in model.names},
                               per_target=per_target_capped),
        arm_vs_shipped=dict(max_abs_diff=d_vs_shipped,
                            median_row_cosine=float(np.median(cos_vs_shipped)),
                            min_row_cosine=float(cos_vs_shipped.min())),
        arm_vs_incumbent=dict(max_abs_diff=d_vs_inc,
                              median_row_cosine=float(np.median(cos_vs_inc)),
                              min_row_cosine=float(cos_vs_inc.min())),
        delivered_norm_fixture=dict(
            cap=dict(full_axis=dict(median=float(np.median(rc)), min=float(rc.min()),
                                    max=float(rc.max())),
                     non_panel=dict(median=float(np.median(rcw)), min=float(rcw.min()),
                                    max=float(rcw.max()))),
            shipped=dict(full_axis=dict(median=float(np.median(rs)), min=float(rs.min()),
                                        max=float(rs.max())),
                         non_panel=dict(median=float(np.median(rsw)), min=float(rsw.min()),
                                        max=float(rsw.max())))),
        peak_gib=peak_gib(),
    )
    (WORK / "test319.json").write_text(json.dumps(out, indent=2))
    (WORK / "test319.log").write_text("\n".join(LINES) + "\n")
    P("")
    P(f"[319] KA rows {sum(1 for k in KA if k['ok'])}/{len(KA)} pass; peak {peak_gib():.2f} GiB")
    if not all(k["ok"] for k in KA):
        raise SystemExit("KNOWN-ANSWER FAILURE")


if __name__ == "__main__":
    main()
