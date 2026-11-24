"""417: VCC_BULK_SHARPEN_Q -- 169's magnitude sharpening on the BALLAST bulk render only (414's recommended arm).

sign(h)|h|^q on each context's CPM gate (profile_from_effect's own gate), the row's gated NON-PANEL L2 norm restored,
panel columns (incl. the row's knockdown) untouched. 414 priced it with panel coordinates zeroed, so protecting them keeps
the scored coordinates exactly as priced and stops the knockdown (the largest |h| in most rows) from soaking up energy.
Only the bulk_effects render gets it; the carriers keep the unsharpened design, so the DE channel is untouched.
Apply only when no build has the builder pinned. Apply AFTER patch_gain410d2.py.
usage: uv run python patch_sharpen417.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
"""
import difflib
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
P = Path(os.environ.get("PATCH_TARGET", ROOT / "scripts/make_submission_pm_cd4.py"))
src = P.read_text()
out = src


def sub(old, new, label):
    global out
    n = out.count(old)
    if n != 1:
        raise SystemExit(f"[{label}] expected exactly one match, found {n}")
    out = out.replace(old, new)


sub('''if BULK_GENE_GAIN is not None and PERT_SCALE_BALLAST is None:
    raise SystemExit("VCC_BULK_GENE_GAIN is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
''', '''if BULK_GENE_GAIN is not None and PERT_SCALE_BALLAST is None:
    raise SystemExit("VCC_BULK_GENE_GAIN is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
# `VCC_BULK_SHARPEN_Q` (417): 169's magnitude sharpening, sign(h)|h|^q on each context's CPM gate with the row's own gated
# NON-PANEL L2 norm restored (panel columns, incl. the row's knockdown, untouched), applied to the BULK render only (the
# ballast profile: pds/mse). The carriers keep the unsharpened design, so the DE channel is untouched. 169 closed it on nmae
# only (a weld cost); 414 priced it bulk-only at +0.006..+0.016 overall vs 410c at q 1.25 (4/4 PLOLO folds positive; the
# flatten and permuted-multiplier controls negative). Requires VCC_PERT_SCALE_BALLAST. Unset is a no-op.
BULK_SHARPEN_Q = os.environ.get("VCC_BULK_SHARPEN_Q", "").strip() or None
if BULK_SHARPEN_Q is not None:
    BULK_SHARPEN_Q = float(BULK_SHARPEN_Q)
    if not 1.0 < BULK_SHARPEN_Q <= 2.0:
        raise SystemExit(f"VCC_BULK_SHARPEN_Q must be in (1, 2] when set: {BULK_SHARPEN_Q!r}")
    if PERT_SCALE_BALLAST is None:
        raise SystemExit("VCC_BULK_SHARPEN_Q is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")
''', "knob")

sub('''def profile_from_effect(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    min_cpm: float = MIN_CPM,
    effect_scale: float = 1.0,
    gene_scale: np.ndarray | None = None,
    shrink_center: np.ndarray | None = None,
) -> np.ndarray:''', '''def bulk_sharpen(gated: np.ndarray, q: float, gate: np.ndarray, protect: np.ndarray) -> np.ndarray:
    """417: sign(h)|h|^q on the gated NON-PANEL coordinates, rescaled to their original L2 norm (169's `sharpen(mode="row")`
    on the panel-zeroed row 414 priced). `protect` (the panel columns, incl. the row's own knockdown) and ungated coordinates
    are returned unchanged."""
    h = np.asarray(gated, dtype=np.float64)
    on = np.asarray(gate, dtype=bool) & ~np.asarray(protect, dtype=bool)
    g = np.where(on, h, 0.0)
    s = np.sign(g) * np.abs(g) ** float(q)
    n0 = float(np.linalg.norm(g))
    n1 = float(np.linalg.norm(s))
    if not np.isfinite(n1) or n1 <= 0.0:
        return h
    return np.where(on, s * (n0 / n1), h)


def profile_from_effect(
    effect: np.ndarray,
    b_ctrl: np.ndarray,
    min_cpm: float = MIN_CPM,
    effect_scale: float = 1.0,
    gene_scale: np.ndarray | None = None,
    shrink_center: np.ndarray | None = None,
    sharpen_q: float | None = None,
    protect: np.ndarray | None = None,
) -> np.ndarray:''', "signature")

sub('''    if MASS_FIX is not None:
        gated = mass_consistent_effect(gated, b_ctrl, gate, mode=MASS_FIX)
    profile = np.clip(np.expm1(b_ctrl + gated), 0.0, None) / TS''', '''    # 417: bulk-only sharpening, after every amplitude knob (it commutes with a global scale: the row norm is restored).
    if sharpen_q is not None:
        if protect is None:
            raise ValueError("sharpen_q needs the panel mask (protect)")
        gated = bulk_sharpen(gated, float(sharpen_q), gate, protect)
    if MASS_FIX is not None:
        gated = mass_consistent_effect(gated, b_ctrl, gate, mode=MASS_FIX)
    profile = np.clip(np.expm1(b_ctrl + gated), 0.0, None) / TS''', "apply")

sub('''        if emitter not in RC_EMITTERS or DE_SIGN_PRIOR is None or context_effects is not None:
            raise ValueError("bulk_effects (ballast) needs an RC emitter, VCC_DE_SIGN_PRIOR and one shared effect matrix")
''', '''        if emitter not in RC_EMITTERS or DE_SIGN_PRIOR is None or context_effects is not None:
            raise ValueError("bulk_effects (ballast) needs an RC emitter, VCC_DE_SIGN_PRIOR and one shared effect matrix")
    # 417: the panel columns bulk sharpening leaves untouched (the scorer drops each row's own target; 414 priced the arm
    # with every panel coordinate zeroed).
    _panel_set = set(map(str, spec.perts))
    panel_mask = np.fromiter((str(g) in _panel_set for g in spec.genes), dtype=bool, count=len(spec.genes))
    if BULK_SHARPEN_Q is not None:
        if bulk_effects is None:
            raise ValueError("VCC_BULK_SHARPEN_Q needs the ballast path (bulk_effects)")
        print(f"[bulk-sharpen] q={BULK_SHARPEN_Q} on the ballast bulk render only; {int(panel_mask.sum())} panel columns "
              "protected", flush=True)
''', "panel_mask")

sub('''                                shrink_center=shrink_center,
                            )
                            ballast_k, ballast_ratio = ballast_count(''', '''                                shrink_center=shrink_center,
                                sharpen_q=BULK_SHARPEN_Q,
                                protect=panel_mask,
                            )
                            ballast_k, ballast_ratio = ballast_count(''', "bulk_call")

sub('''            **({} if BULK_GENE_GAIN is None else {"bulk_gene_gain": BULK_GENE_GAIN, "bulk_gene_gain_sha256": BULK_GENE_GAIN_ARTIFACTS[BULK_GENE_GAIN][1]}),
''', '''            **({} if BULK_GENE_GAIN is None else {"bulk_gene_gain": BULK_GENE_GAIN, "bulk_gene_gain_sha256": BULK_GENE_GAIN_ARTIFACTS[BULK_GENE_GAIN][1]}),
            **({} if BULK_SHARPEN_Q is None else {"bulk_sharpen_q": BULK_SHARPEN_Q}),
''', "config")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
