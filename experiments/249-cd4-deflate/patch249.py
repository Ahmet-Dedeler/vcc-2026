"""249 — the ONE builder change, applied as exact-anchor replacements so the diff is auditable.

Adds `VCC_CD4_DEFLATE_K` (unset = byte-inert) to `scripts/make_submission_pm_cd4.py`:
  1. the env flag block (after 247's `VCC_CD4_TARGET_G` validation);
  2. `deflate_cd4_residual()` -- 246 `lib246.deflate(res, "svd", k)` verbatim;
  3. the call in `apply_frozen_cd4_residual`, BEFORE `solve_cd4_dose`, plus the audit key;
  4. `config.cd4_deflate_k`, written only when on.
Every anchor must occur exactly once; the script refuses to run twice.

    uv run python experiments/249-cd4-deflate/patch249.py            # the working tree
    uv run python experiments/249-cd4-deflate/patch249.py <copy>     # any byte copy of HEAD
(the one follow-up edit, the supplied-basis-safe log line, is applied by `intended249.py`)
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
import sys

PATH = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / "scripts/make_submission_pm_cd4.py"
src = PATH.read_text()
assert "VCC_CD4_DEFLATE_K" not in src, "already patched"


def replace_once(text: str, old: str, new: str) -> str:
    n = text.count(old)
    assert n == 1, f"anchor found {n} times:\n{old[:200]}"
    return text.replace(old, new)


# ------------------------------------------------------------------ 1. the env flag
OLD1 = '''if CD4_TARGET_G is not None and not 0.0 < CD4_TARGET_G < 5.0:
    raise ValueError(f"VCC_CD4_TARGET_G must be in (0, 5); got {CD4_TARGET_G}")
'''
NEW1 = OLD1 + '''# --- `VCC_CD4_DEFLATE_K` -- remove the CD4 residual's top-k TARGET-SHARED directions (246/249).
#
# THE MECHANISM (246 section 3).  The CD4 residual enters un-whitened at `f` 0.34, so nothing
# ever removes its target-shared component -- the T-cell-programme hub the whitener strips from
# the donors (209) and nothing strips from CD4 because it bypasses the operator (207 section 5).
# Its top-5 left-singular directions carry ~20% of its energy across the carried rows.  246
# measured removing them on 162's PLOLO instrument, three cancer folds, with the fraction
# re-solved on the deflated object (`B2_cd4_svd5_resolved`): `pds` +0.01226 raw (3/3 folds,
# 8.2 LOO-sigma) for -0.46% aggregate / -0.81% macro cosine.  Energy-matched random subspaces
# (the same ~21% of energy) buy no `pds` (-0.0001 / -0.0050) and cost 3x the cosine, so it is
# the DIRECTION that pays, not the norm cut.
#
# THE CONSTRUCTION (`experiments/246-direction-hypotheses/lib246.py` `deflate(res, "svd", k)`,
# reproduced verbatim in `deflate_cd4_residual` below).  On the carried rows
# (`target_condition_count > 0`): subtract the per-gene mean over those rows, `G = R_c R_c^T`,
# `np.linalg.eigh`, remove the top-k eigen-directions `R_k = (I - U_k U_k^T) R_c`, re-zero
# outside the gene support.  Deterministic, no seed.
#
# THE FLAG.  Unset (the default) is BYTE-INERT: `apply_frozen_cd4_residual` calls nothing new
# and every prior build reproduces bit-for-bit (249's known-answer rows: the staged effects
# array `abeb0749...` and the shipped package `8df4dadd...`).  Set to an integer k >= 1 and the
# CD4 residual is replaced by its deflated form BEFORE `solve_cd4_dose` sizes it, so the
# existing rule sizes the object actually added (246's "resolved" variant; the solved fraction
# rises because the deflated residual's median row norm is lower).  The H1 union residual never
# sees the flag; `solve_cd4_fraction` is not modified.  `effects.cd4.deflation` and
# `config.cd4_deflate_k` are written only when on (the `h1_union_g` convention).
CD4_DEFLATE_K = os.environ.get("VCC_CD4_DEFLATE_K", "").strip()
if CD4_DEFLATE_K:
    try:
        CD4_DEFLATE_K = int(CD4_DEFLATE_K)
    except ValueError as exc:
        raise ValueError(
            f"VCC_CD4_DEFLATE_K must be an integer >= 1; got {CD4_DEFLATE_K!r}"
        ) from exc
    if CD4_DEFLATE_K < 1:
        raise ValueError(f"VCC_CD4_DEFLATE_K must be an integer >= 1; got {CD4_DEFLATE_K}")
else:
    CD4_DEFLATE_K = None
'''
src = replace_once(src, OLD1, NEW1)

# ------------------------------------------------------------------ 2. the function
OLD2 = '''def single_source_residual(
    bank,
'''
NEW2 = '''def deflate_cd4_residual(cd4, k: int, *, basis: np.ndarray | None = None):
    """246's `lib246.deflate(res, "svd", k)`, verbatim: remove the top-k target-shared directions.

    On the CARRIED rows (`target_condition_count > 0`, 246's definition -- production's
    `freeze_complete_mean3` zeroes incomplete rows' effects but keeps their counts, so those
    rows enter the centring as zero rows exactly as they do under 162's `cd4_for`):

        R   = residual on the carried rows, masked to its gene support
        R_c = R - mean over rows(R)          per-gene centring, re-zeroed outside the support
        G   = R_c R_c^T                      (n_carried x n_carried)
        U_k = top-k eigenvectors of G        (`np.linalg.eigh`; deterministic, no seed)
        R_k = (I - U_k U_k^T) R_c            re-zeroed outside the support

    `basis` exists for CONTROLS only: an orthonormal `(n_carried, k')` matrix substituted for
    `U_k`, so an energy-matched random subspace can be pushed through the identical path
    (249's known-answer row 3: a projector, not a norm cut).  Production never passes it.

    Returns a new `CD4Prepared` with the SAME mode, votes, counts and condition details -- so
    `effects.cd4.mode`, `supported_entries` and `fully_supported_targets` cannot move, and the
    rule-10 diff reads only the sizing consequences -- plus the audit the meta records under
    `effects.cd4.deflation`.
    """
    if int(k) < 1:
        raise ValueError(f"CD4 deflation needs k >= 1; got {k}")
    eff = np.asarray(cd4.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(cd4.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(np.asarray(cd4.target_condition_count) > 0)
    if rows.size == 0:
        raise ValueError(
            "VCC_CD4_DEFLATE_K is set but the CD4 residual carries no rows; nothing to deflate"
        )
    R = eff[rows]
    e0 = float((R * R).sum())
    if e0 <= 0.0:
        raise ValueError(
            "VCC_CD4_DEFLATE_K is set but the CD4 residual has zero energy; nothing to deflate"
        )
    audit: dict[str, object] = {"k": int(k), "n_carried": int(rows.size)}
    mean = R.mean(axis=0, keepdims=True)
    Rc = R - mean
    Rc[~sup[rows]] = 0.0
    if basis is None:
        G = Rc @ Rc.T
        w, U = np.linalg.eigh(G)
        Uk = U[:, np.argsort(w)[::-1][:k]]
        audit["top_k_energy_frac_of_centred"] = float(
            np.sort(w)[::-1][:k].sum() / max(w.sum(), 1e-300)
        )
    else:
        Uk = np.asarray(basis, dtype=np.float64)
        if Uk.ndim != 2 or Uk.shape[0] != rows.size:
            raise ValueError(
                f"control basis must be (n_carried={rows.size}, k'); got {Uk.shape}"
            )
        audit["basis"] = f"supplied control, {Uk.shape[1]} dims"
    R = Rc - Uk @ (Uk.T @ Rc)
    R[~sup[rows]] = 0.0
    e1 = float((R * R).sum())
    audit["energy_removed_frac"] = float(1.0 - e1 / e0)
    audit["basis_u_k_sha256"] = array_sha256(np.ascontiguousarray(Uk))
    eff[rows] = R
    return (
        EXP062.CD4Prepared(
            mode=str(cd4.mode),
            mean_effect=eff,
            gene_votes=np.asarray(cd4.gene_votes).copy(),
            target_condition_count=np.asarray(cd4.target_condition_count).copy(),
            condition_details=cd4.condition_details,
        ),
        audit,
    )


def single_source_residual(
    bank,
'''
src = replace_once(src, OLD2, NEW2)

# ------------------------------------------------------------------ 3. the call site
OLD3 = '''    fraction, frame = solve_cd4_dose(whitened, cd4)
    if extra is None:
        added = cd4
'''
NEW3 = '''    # 249: with `VCC_CD4_DEFLATE_K` set, the CD4 residual is replaced by its deflated form
    # BEFORE it is sized, so the sizing rule (247's router, 099's all-rows rule underneath)
    # sizes the object actually added.  Unset: nothing here runs and `cd4` is untouched.
    # The second residual (`extra`, the H1 union) is never deflated.
    deflation_audit: dict[str, object] | None = None
    if CD4_DEFLATE_K is not None:
        cd4, deflation_audit = deflate_cd4_residual(cd4, CD4_DEFLATE_K)
        if np.any(cd4.mean_effect[:, panel_positions] != 0):
            raise AssertionError("deflation put a VCC panel coordinate into the CD4 residual")
        print(
            f"[cd4] deflation k={CD4_DEFLATE_K}: {deflation_audit['n_carried']} carried rows, "
            f"energy removed {deflation_audit['energy_removed_frac']:.6f} (top-k share of the "
            f"centred energy {deflation_audit['top_k_energy_frac_of_centred']:.6f}), "
            f"U_k sha256 {str(deflation_audit['basis_u_k_sha256'])[:16]}",
            flush=True,
        )
    fraction, frame = solve_cd4_dose(whitened, cd4)
    if extra is None:
        added = cd4
'''
src = replace_once(src, OLD3, NEW3)

OLD3b = '''    cd4_audit = {
        **cd4_support_audit(cd4, fraction),
        "residual_frame": frame,
        "post_add_preserve": reported_preserve,
        "panel_residual_nonzero_entries": int(
            np.count_nonzero(cd4.mean_effect[:, panel_positions])
        ),
    }
'''
NEW3b = '''    cd4_audit = {
        **cd4_support_audit(cd4, fraction),
        "residual_frame": frame,
        "post_add_preserve": reported_preserve,
        "panel_residual_nonzero_entries": int(
            np.count_nonzero(cd4.mean_effect[:, panel_positions])
        ),
        # Written only when `VCC_CD4_DEFLATE_K` is set (249), the `h1_union` convention.
        **({} if deflation_audit is None else {"deflation": deflation_audit}),
    }
'''
src = replace_once(src, OLD3b, NEW3b)

# ------------------------------------------------------------------ 4. the config key
OLD4 = '''            **(
                {}
                if CD4_SIZING is None
                else {"cd4_sizing": CD4_SIZING, "cd4_target_g": CD4_TARGET_G}
            ),
'''
NEW4 = OLD4 + '''            # Written only when the CD4 residual is deflated (249), the same convention.
            **({} if CD4_DEFLATE_K is None else {"cd4_deflate_k": CD4_DEFLATE_K}),
'''
src = replace_once(src, OLD4, NEW4)

PATH.write_text(src)
print(f"patched {PATH}")
