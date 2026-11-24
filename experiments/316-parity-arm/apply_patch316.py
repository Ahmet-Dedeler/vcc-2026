#!/usr/bin/env python3
"""316 -- wire `VCC_NORMREL_PM=blend` into the builder: the PARITY-CORRECT rendering of 312's arm.

313 measured that the SHIPPED `VCC_NORMREL_EXP` arm is not the object 308 scored: PerturbMatch is
applied at the UNIFORM donor share and OUTSIDE the norm match, so the delivered rows change
amplitude (per-target delivered norm 0.8062-2.1549 on scored columns, whole-matrix +6.00 %) where
308's arm is exactly shape-only.  `VCC_NORMREL_PM=blend` implements 313 section 8's specification:

    corr(t)  = alpha * effect_scale * PM_BLEND * (pm_standardized(t) - incumbent_k562(t))
    base_pm  = effect(t, uniform)        + k562_share * corr(t)      # the norm reference
    cand_pm  = effect(t, normrel weights)+ w_k562(t)   * corr(t)
    row(t)   = cand_pm * ||base_pm|| / ||cand_pm||                   # norm match LAST

Default UNSET = today's shipped behaviour, bit for bit.  Setting it without `VCC_NORMREL_EXP` is a
hard error.  D5 is frozen to 308's support (the operator columns, i.e. the donor's coverage MINUS
the panel); D7 keeps production's broad PerturbMatch column policy.
"""
import hashlib
import pathlib
import shutil
import sys

ROOT = pathlib.Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/316-parity-arm"
B = ROOT / "scripts/make_submission_pm_cd4.py"
D = ROOT / "scripts/def_pipeline.py"

src = B.read_text()
before_sha = hashlib.sha256(src.encode()).hexdigest()
if "VCC_NORMREL_PM" in src:
    print("already patched; nothing to do")
    sys.exit(0)
shutil.copyfile(B, HERE / "work/builder.PRE316.py")

# ------------------------------------------------------------------ 1) the env knob
anchor_env = '''NORMREL_EXP_ENV = os.environ.get("VCC_NORMREL_EXP", "").strip()
NORMREL_EXP = float(NORMREL_EXP_ENV) if NORMREL_EXP_ENV else None
if NORMREL_EXP is not None and not (-4.0 <= NORMREL_EXP <= 4.0):
    raise SystemExit(f"VCC_NORMREL_EXP out of range [-4,4]: {NORMREL_EXP_ENV!r}")'''
assert anchor_env in src, "NORMREL_EXP env anchor not found"
env_block = anchor_env + '''

# `VCC_NORMREL_PM` -- the PARITY-CORRECT rendering of that arm (316, to 313 section 8).
# `blend` moves the PerturbMatch correction INSIDE the norm match and applies it at the arm's own
# per-target K562 weight instead of the uniform donor share.  That is the object 308 measured
# (`ORD_A`), and it makes the emitted row EXACTLY norm-preserving -- shape-only, AGENTS 8.1d's
# +0.0019 gate rather than the +0.0089 amplitude gate.  Under `blend` the per-target norms are also
# taken over the OPERATOR columns (308's `wcols`: the donor's coverage minus the panel), 313's D5.
# UNSET (the default) is 312's shipped weight-first behaviour, bit for bit.
NORMREL_PM_ENV = os.environ.get("VCC_NORMREL_PM", "").strip()
NORMREL_PM = NORMREL_PM_ENV or None
if NORMREL_PM is not None and NORMREL_PM != "blend":
    raise SystemExit(f"VCC_NORMREL_PM must be 'blend' when set: {NORMREL_PM_ENV!r}")
if NORMREL_PM is not None and NORMREL_EXP is None:
    raise SystemExit(
        "VCC_NORMREL_PM is set without VCC_NORMREL_EXP: it changes nothing when the per-target "
        "weighting is off, so this is refused rather than silently building the incumbent"
    )'''
src = src.replace(anchor_env, env_block, 1)

# ------------------------------------------------------------------ 2) the dispatcher's signature
anchor_sig = '''def build_v9_pre_whiten_effects(
    model: MultiDonorModel,
    targets: Sequence[str],
    *,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
) -> np.ndarray:
    """Construct the exact pre-whitening effect matrix used by v9."""'''
assert anchor_sig in src, "build_v9_pre_whiten_effects signature not found"
new_sig = '''def build_v9_pre_whiten_effects(
    model: MultiDonorModel,
    targets: Sequence[str],
    *,
    alpha: float = ALPHA,
    beta: float = BETA,
    weights: Mapping[str, float] | None = None,
    pm: Mapping[str, object] | None = None,
    operator_columns: np.ndarray | None = None,
    pm_handled_out: set[str] | None = None,
) -> np.ndarray:
    """Construct the exact pre-whitening effect matrix used by v9.

    `pm` / `operator_columns` / `pm_handled_out` are used ONLY by `VCC_NORMREL_PM=blend` (316):
    that path folds the PerturbMatch correction into the row before the norm match, and reports
    back through `pm_handled_out` which targets it handled so that
    `apply_perturbmatch_k562_replacement` does not apply the same correction a second time."""'''
src = src.replace(anchor_sig, new_sig, 1)

anchor_call = '''    return _normrel_pre_whiten_effects(
        model, targets, alpha=alpha, beta=beta, active_weights=active_weights,
        exponent=NORMREL_EXP,
    )'''
assert anchor_call in src, "_normrel dispatch call not found"
src = src.replace(anchor_call, '''    return _normrel_pre_whiten_effects(
        model, targets, alpha=alpha, beta=beta, active_weights=active_weights,
        exponent=NORMREL_EXP, pm=pm, operator_columns=operator_columns,
        handled_out=pm_handled_out,
    )''', 1)

# ------------------------------------------------------------------ 3) the arm itself
old_fn = '''def _normrel_pre_whiten_effects(
    model: MultiDonorModel,
    targets: Sequence[str],
    *,
    alpha: float,
    beta: float,
    active_weights: Mapping[str, float],
    exponent: float,
) -> np.ndarray:
    """097 a6 / 308 `normrel1.5`: per-target donor weights, norm-matched to the uniform row."""
    panel = [str(t) for t in targets]
    median_norm: dict[str, float] = {}
    for name in model.names:
        idx = model.pert_index[name]
        cols = model.dcols[name]
        vals = [
            float(np.linalg.norm(model.s[name][idx[t]][cols])) for t in panel if t in idx
        ]
        median_norm[name] = float(np.median(vals)) if vals else 0.0

    rows: list[np.ndarray] = []
    for target in panel:
        base = np.asarray(
            model.effect(target, alpha, beta, weights=active_weights), dtype=np.float64
        )
        have = [n for n in model.names if target in model.pert_index[n]]
        w: dict[str, float] = {}
        for name in have:
            med = median_norm.get(name, 0.0)
            nv = float(
                np.linalg.norm(
                    model.s[name][model.pert_index[name][target]][model.dcols[name]]
                )
            )
            if med > 0.0 and nv > 0.0:
                w[name] = float((nv / med) ** exponent) * float(
                    active_weights.get(name, 1.0)
                )
            else:
                w[name] = float(active_weights.get(name, 1.0))
        if not w or sum(w.values()) <= 0.0:
            rows.append(base)
            continue
        cand = np.asarray(
            model.effect(target, alpha, beta, weights=w), dtype=np.float64
        )
        nb = float(np.linalg.norm(base))
        nc = float(np.linalg.norm(cand))
        if nb > 0.0 and nc > 0.0:
            cand = cand * (nb / nc)
        rows.append(cand)
    return np.vstack(rows).astype(np.float64)'''
assert old_fn in src, "_normrel_pre_whiten_effects body not found verbatim"

new_fn = '''def _normrel_pre_whiten_effects(
    model: MultiDonorModel,
    targets: Sequence[str],
    *,
    alpha: float,
    beta: float,
    active_weights: Mapping[str, float],
    exponent: float,
    pm: Mapping[str, object] | None = None,
    operator_columns: np.ndarray | None = None,
    handled_out: set[str] | None = None,
) -> np.ndarray:
    """097 a6 / 308 `normrel1.5`: per-target donor weights, norm-matched to the uniform row.

    Two renderings of the same arm:

    * `VCC_NORMREL_PM` UNSET (312, shipped) -- the norm is taken over the donor's whole official
      coverage (the panel columns, i.e. the target's own knockdown, included) and the row is
      norm-matched to the PerturbMatch-FREE uniform row.  `apply_perturbmatch_k562_replacement`
      then adds the PerturbMatch correction on top, at the UNIFORM donor share.  313 measured that
      this is amplitude-changing (delivered norm 0.81-2.15 on the scored columns).
    * `VCC_NORMREL_PM=blend` (316, 313 section 8) -- the correction is folded into BOTH the
      reference and the candidate first, the candidate carrying it at the arm's own K562 weight,
      and the norm match is applied LAST.  `||row|| == ||incumbent row||` identically, so the arm
      is exactly shape-only.  The per-target norms are taken over the operator columns (308's
      `wcols` = the donor's coverage minus the panel; 313's D5, frozen to 308's choice).
    """
    panel = [str(t) for t in targets]
    blend_pm = NORMREL_PM == "blend"
    if blend_pm:
        if pm is None:
            raise SystemExit("VCC_NORMREL_PM=blend requires the PerturbMatch mapping")
        if operator_columns is None:
            raise SystemExit(
                "VCC_NORMREL_PM=blend requires the operator columns (313 D5: the per-target "
                "norms are taken over them)"
            )
        if K562 not in model.names:
            raise SystemExit("VCC_NORMREL_PM=blend requires the K562 slot")

    # The support the per-target norms are taken over.  Shipped: the donor's whole coverage.
    # `blend`: that coverage intersected with the operator columns, which are panel-free by
    # construction (`production_operator_columns`), i.e. exactly 308's panel-zeroed `wcols`.
    if blend_pm:
        op = np.asarray(operator_columns, dtype=np.int64)
        norm_cols: dict[str, np.ndarray] = {
            name: np.asarray(model.dcols[name])[
                np.isin(np.asarray(model.tcols[name], dtype=np.int64), op)
            ]
            for name in model.names
        }
    else:
        norm_cols = {name: model.dcols[name] for name in model.names}

    median_norm: dict[str, float] = {}
    for name in model.names:
        idx = model.pert_index[name]
        cols = norm_cols[name]
        vals = [
            float(np.linalg.norm(model.s[name][idx[t]][cols])) for t in panel if t in idx
        ]
        median_norm[name] = float(np.median(vals)) if vals else 0.0

    if blend_pm:
        k_scale = float(model._scale[K562]) or 1.0
        pm_scale = float(pm["scale"])
        pm_target_columns = np.asarray(pm["target_columns"], dtype=np.int64)
        pm_source_columns = np.asarray(pm["source_columns"], dtype=np.int64)
        pm_mean = np.asarray(pm["mean"], dtype=np.float64)
        pm_delta = pm["delta"]
        pm_pert_index = pm["pert_index"]

    rows: list[np.ndarray] = []
    for target in panel:
        base = np.asarray(
            model.effect(target, alpha, beta, weights=active_weights), dtype=np.float64
        )
        have = [n for n in model.names if target in model.pert_index[n]]
        w: dict[str, float] = {}
        for name in have:
            med = median_norm.get(name, 0.0)
            nv = float(
                np.linalg.norm(
                    model.s[name][model.pert_index[name][target]][norm_cols[name]]
                )
            )
            if med > 0.0 and nv > 0.0:
                w[name] = float((nv / med) ** exponent) * float(
                    active_weights.get(name, 1.0)
                )
            else:
                w[name] = float(active_weights.get(name, 1.0))
        if not w or sum(w.values()) <= 0.0:
            rows.append(base)
            continue
        cand = np.asarray(
            model.effect(target, alpha, beta, weights=w), dtype=np.float64
        )
        if (
            blend_pm
            and target in model.pert_index[K562]
            and target in pm_pert_index
        ):
            # `apply_perturbmatch_k562_replacement`'s own `k562_share`, computed by the same
            # expression on the same objects, so that at exponent 0 (where `w` collapses to
            # `active_weights`) this path is bitwise identical to the incumbent's.
            pm_have = [
                name
                for name in model.names
                if target in model.pert_index[name]
                and active_weights.get(name, 0.0) > 0
            ]
            total_weight = float(sum(active_weights[name] for name in pm_have))
            if total_weight <= 0:
                raise ValueError(f"{target}: no positive incumbent source weight")
            k562_share = float(active_weights[K562]) / total_weight
            k562_weight = float(w[K562]) / float(sum(w.values()))

            incumbent = np.zeros(len(model.target_genes), dtype=np.float64)
            k_row = model.pert_index[K562][target]
            incumbent[model.tcols[K562]] = (
                model.s[K562][k_row][model.dcols[K562]] / k_scale
            )
            pm_row = pm_pert_index[target]
            pm_standardized = (
                np.asarray(pm_delta[pm_row, pm_source_columns], dtype=np.float64)
                - pm_mean
            ) / pm_scale
            drive = pm_standardized - incumbent[pm_target_columns]
            # The scalar chain is production's, left to right, so the exponent-0 row reproduces
            # `cd4svd5` to the last bit rather than to the last ulp.
            base[pm_target_columns] += (
                alpha * k562_share * float(model.effect_scale) * PM_BLEND * drive
            )
            cand[pm_target_columns] += (
                alpha * k562_weight * float(model.effect_scale) * PM_BLEND * drive
            )
            if handled_out is not None:
                handled_out.add(target)
        nb = float(np.linalg.norm(base))
        nc = float(np.linalg.norm(cand))
        if nb > 0.0 and nc > 0.0:
            cand = cand * (nb / nc)
        rows.append(cand)
    return np.vstack(rows).astype(np.float64)'''
src = src.replace(old_fn, new_fn, 1)

# ------------------------------------------------------------------ 4) do not blend twice
anchor_pm_sig = '''    alpha: float = ALPHA,
    weights: Mapping[str, float] | None = None,
    blend: float = PM_BLEND,
) -> tuple[np.ndarray, dict[str, object]]:
    """Replace one K562 slot while preserving target-level source availability.'''
assert anchor_pm_sig in src, "PM function signature not found"
src = src.replace(anchor_pm_sig, '''    alpha: float = ALPHA,
    weights: Mapping[str, float] | None = None,
    blend: float = PM_BLEND,
    already_blended: Iterable[str] | None = None,
) -> tuple[np.ndarray, dict[str, object]]:
    """Replace one K562 slot while preserving target-level source availability.

    `already_blended` (316) lists the targets whose PerturbMatch correction was already folded in
    by `_normrel_pre_whiten_effects` under `VCC_NORMREL_PM=blend`.  Those rows are left untouched
    here -- applying the correction twice would double-count it -- but they are still counted and
    reported in the audit exactly as before, so the postbuild drift check is unchanged.''', 1)

anchor_pm_body = '''    baseline = np.asarray(baseline, dtype=np.float64)
    if baseline.shape != (len(targets), len(model.target_genes)):'''
assert anchor_pm_body in src, "PM function body anchor not found"
src = src.replace(anchor_pm_body, '''    blended_already = set(map(str, already_blended)) if already_blended else set()
    baseline = np.asarray(baseline, dtype=np.float64)
    if baseline.shape != (len(targets), len(model.target_genes)):''', 1)

anchor_hist = '''        source_vote_histogram[str(len(have))] = (
            source_vote_histogram.get(str(len(have)), 0) + 1
        )
'''
assert anchor_hist in src, "PM histogram anchor not found"
src = src.replace(anchor_hist, anchor_hist + '''        if target in blended_already:
            # 316: the row already carries this correction, applied at the arm's own per-target
            # K562 weight and INSIDE the norm match.  The audit is emitted unchanged.
            support["both_supported_modified"].append(target)
            continue
''', 1)

# ------------------------------------------------------------------ 5) thread it through the build
anchor_build = '''    active_weights = dict(weights) if weights is not None else donor_weights(model)
    baseline_pre = build_v9_pre_whiten_effects(
        model,
        targets,
        alpha=alpha,
        beta=beta,
        weights=active_weights,
    )
    if pm is None:
        candidate_pre = baseline_pre.copy()
        pm_audit: dict[str, object] = {"enabled": False}
    else:
        candidate_pre, pm_audit = apply_perturbmatch_k562_replacement(
            baseline_pre,
            model,
            pm,
            targets,
            alpha=alpha,
            weights=active_weights,
        )'''
assert anchor_build in src, "build_production_effects anchor not found"
src = src.replace(anchor_build, '''    active_weights = dict(weights) if weights is not None else donor_weights(model)
    pm_handled: set[str] = set()
    baseline_pre = build_v9_pre_whiten_effects(
        model,
        targets,
        alpha=alpha,
        beta=beta,
        weights=active_weights,
        pm=pm,
        operator_columns=operator_columns,
        pm_handled_out=pm_handled,
    )
    if pm is None:
        candidate_pre = baseline_pre.copy()
        pm_audit: dict[str, object] = {"enabled": False}
    else:
        candidate_pre, pm_audit = apply_perturbmatch_k562_replacement(
            baseline_pre,
            model,
            pm,
            targets,
            alpha=alpha,
            weights=active_weights,
            already_blended=pm_handled,
        )''', 1)

# ------------------------------------------------------------------ 6) the meta config
anchor_cfg = '''            # Written only when the per-target normrel weighting is ON (312).
            **({} if NORMREL_EXP is None else {"normrel_exp": NORMREL_EXP}),'''
assert anchor_cfg in src, "normrel_exp config anchor not found"
src = src.replace(anchor_cfg, anchor_cfg + '''
            # Written only when that weighting uses the parity-correct PerturbMatch
            # rendering (316), the same convention.
            **({} if NORMREL_PM is None else {"normrel_pm": NORMREL_PM}),''', 1)

B.write_text(src)
print("patched builder")

# ------------------------------------------------------------------ 7) register the arm
dsrc = D.read_text()
if "VCC_NORMREL_PM" not in dsrc:
    danchor = '''                   "VCC_CD4_DEFLATE_BASIS": "+0.0042634 (`cd4svd5`) carried to a new panel as the PINNED "'''
    assert danchor in dsrc, "ARM_BOARD_DELTA anchor not found"
    dnew = ('''                   "VCC_NORMREL_PM": "the PARITY-CORRECT rendering of VCC_NORMREL_EXP (316, to 313 "
                                     "section 8): PerturbMatch folded in at the arm's own per-target K562 "
                                     "weight and INSIDE the norm match, so the row is exactly shape-only. "
                                     "The shipped weight-first form scored +0.0016245 (Q519I2pLWRCW7PAEx5ul) "
                                     "with score_nmae -0.007705 -- UNTESTED on the board",
''') + danchor
    shutil.copyfile(D, HERE / "work/def_pipeline.PRE316.py")
    dsrc = dsrc.replace(danchor, dnew, 1)
    D.write_text(dsrc)
    print("def_pipeline.py: VCC_NORMREL_PM registered as a known arm")

print("builder before sha256:", before_sha)
print("builder after  sha256:", hashlib.sha256(src.encode()).hexdigest())
print("def_pipeline after sha256:", hashlib.sha256(D.read_text().encode()).hexdigest())
