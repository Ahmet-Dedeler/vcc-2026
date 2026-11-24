#!/usr/bin/env python3
"""312 -- wire VCC_NORMREL_EXP into the builder: per-target donor weights (097 a6 / 308 normrel1.5).

w_j(t) proportional to (||row_j(t)|| / median_t ||row_j(t)||) ** exp, exp = -1.5, renormalised over
the donors carrying t, then the row is NORM-MATCHED to the uniform-weight row (shape-only, 263 ':nm').

Default OFF: with VCC_NORMREL_EXP unset the code path is the original expression, unchanged, so an
unset build is bit-identical to cd4svd5. Exponent 0 gives uniform weights and must also be identical.
"""
import re, sys, hashlib, shutil, pathlib

B = pathlib.Path("~/Code/virtual-cell-challenge-2026/scripts/make_submission_pm_cd4.py")
src = B.read_text()
before_sha = hashlib.sha256(src.encode()).hexdigest()
if "VCC_NORMREL_EXP" in src:
    print("already patched; nothing to do"); sys.exit(0)

# 1) the env knob, placed right after the DONOR_WEIGHTS block
anchor_env = 'DONOR_WEIGHTS_ENV = os.environ.get("VCC_DONOR_WEIGHTS", "").strip()'
assert anchor_env in src, "donor-weights env anchor not found"
env_block = anchor_env + '''

# `VCC_NORMREL_EXP` -- per-target donor weighting (097 route (a) a6; 308 `normrel1.5`).
# w_j(t) is proportional to (||row_j(t)|| / median over that donor's panel rows)**exp, renormalised
# over the donors carrying t; the resulting row is then rescaled to the uniform-weight row's norm,
# so the arm is SHAPE-ONLY (263's `:nm`; AGENTS 8.1d -- signs pay, magnitudes do not).
# UNSET (the default) leaves the original expression untouched and the build bit-identical.
NORMREL_EXP_ENV = os.environ.get("VCC_NORMREL_EXP", "").strip()
NORMREL_EXP = float(NORMREL_EXP_ENV) if NORMREL_EXP_ENV else None
if NORMREL_EXP is not None and not (-4.0 <= NORMREL_EXP <= 4.0):
    raise SystemExit(f"VCC_NORMREL_EXP out of range [-4,4]: {NORMREL_EXP_ENV!r}")'''
src = src.replace(anchor_env, env_block, 1)

# 2) the per-target path inside build_v9_pre_whiten_effects
old_body = '''    active_weights = dict(weights) if weights is not None else donor_weights(model)
    return np.vstack(
        [
            model.effect(
                str(target),
                alpha,
                beta,
                weights=active_weights,
            )
            for target in targets
        ]
    ).astype(np.float64)'''
assert old_body in src, "build_v9_pre_whiten_effects body not found verbatim"
new_body = '''    active_weights = dict(weights) if weights is not None else donor_weights(model)
    if NORMREL_EXP is None:
        return np.vstack(
            [
                model.effect(
                    str(target),
                    alpha,
                    beta,
                    weights=active_weights,
                )
                for target in targets
            ]
        ).astype(np.float64)
    return _normrel_pre_whiten_effects(
        model, targets, alpha=alpha, beta=beta, active_weights=active_weights,
        exponent=NORMREL_EXP,
    )


def _normrel_pre_whiten_effects(
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
        base = model.effect(target, alpha, beta, weights=active_weights)
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
        cand = model.effect(target, alpha, beta, weights=w)
        nb = float(np.linalg.norm(base))
        nc = float(np.linalg.norm(cand))
        if nb > 0.0 and nc > 0.0:
            cand = cand * (nb / nc)
        rows.append(cand)
    return np.vstack(rows).astype(np.float64)'''
src = src.replace(old_body, new_body, 1)


# 3) record the knob in the meta config, same "written only when set" convention
anchor_cfg = '            **({} if CD4_PARTIAL_G is None else {"cd4_partial_g": CD4_PARTIAL_G}),'
assert anchor_cfg in src, "config anchor not found"
cfg_add = (
    anchor_cfg
    + "\n            # Written only when the per-target normrel weighting is ON (312).\n"
    + '            **({} if NORMREL_EXP is None else {"normrel_exp": NORMREL_EXP}),'
)
src = src.replace(anchor_cfg, cfg_add, 1)

shutil.copyfile(B, "~/Code/virtual-cell-challenge-2026/experiments/312-normrel-arm/work/builder.PRE312.py")
B.write_text(src)
print("patched")

# 4) register the knob as a known arm in the driver, so provenance is announced and the
#    postbuild drift check treats it as a declared arm rather than an unexplained key.
D = pathlib.Path("~/Code/virtual-cell-challenge-2026/scripts/def_pipeline.py")
dsrc = D.read_text()
if "VCC_NORMREL_EXP" not in dsrc:
    danchor = '"VCC_CD4_DEFLATE_BASIS": "+0.0042634 (`cd4svd5`) carried to a new panel as the PINNED "'
    assert danchor in dsrc, "ARM_BOARD_DELTA anchor not found"
    dnew = ('"VCC_NORMREL_EXP": "per-target donor weighting, 097 route (a) a6 / 308 `normrel1.5`; "\n'
            '                                       "local PLOLO +0.00259 [+0.00000, +0.00527] 3/3 folds, "\n'
            '                                       "priced +0.00104 pessimistic -- UNTESTED on the board",\n'
            '                   ') + danchor
    dsrc = dsrc.replace(danchor, dnew, 1)
    shutil.copyfile(D, "~/Code/virtual-cell-challenge-2026/experiments/312-normrel-arm/work/def_pipeline.PRE312.py")
    D.write_text(dsrc)
    print("def_pipeline.py: VCC_NORMREL_EXP registered as a known arm")

print("before sha256:", before_sha)
print("after  sha256:", hashlib.sha256(src.encode()).hexdigest())
