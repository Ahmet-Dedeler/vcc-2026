"""410d: add VCC_BULK_GENE_GAIN (228's per-gene reliability gain, BULK ONLY under ballast) to the builder.

Run ONLY when no build has the builder pinned (after a run*.sh logs "builder after restore"). Exact edits, each asserted:
  1. module-level knob + artifact pin + apply_bulk_gene_gain()
  2. run(): after the per-target scale (and only with VCC_PERT_SCALE_BALLAST), multiply the BULK matrix by m per gene,
     energy-preserving over non-panel columns; the carriers keep the unscaled, ungained matrix
  3. config: bulk_gene_gain + its sha
usage: uv run python patch_gain410d.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
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


KNOB = '''

# `VCC_BULK_GENE_GAIN` -- arm 410d (2026-09-23): 228's per-GENE cross-line reliability gain (224's diagonal, fitted on 1,758
# NON-panel targets of six lines through the production operator; gain_reg, dose t = 0.25: m_g = 1 + t (g/median - 1)),
# applied to the BULK matrix only (pds/mse), never to the carriers. 228 priced it at pds +0.0018 / mse +0.0023 overall and
# killed it on nmae (-0.0104), which was a WELD cost; ballast removes the weld. PLOLO: it stacks multiplicatively with the
# per-target key (+0.52 relative aggregate cosine uncapped + gain, stack410.py). Requires VCC_PERT_SCALE_BALLAST. Unset is a
# no-op.
BULK_GENE_GAIN = os.environ.get("VCC_BULK_GENE_GAIN", "").strip().lower() or None
BULK_GENE_GAIN_ARTIFACTS = {
    "g228reg_t025": (
        ROOT / "experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz",
        "da937cfa69e0730b215b8780c5ca0505d7d4acd6f0ae957891ca74731971fd9f",
    ),
}
if BULK_GENE_GAIN is not None and BULK_GENE_GAIN not in BULK_GENE_GAIN_ARTIFACTS:
    raise SystemExit(f"VCC_BULK_GENE_GAIN must be one of {tuple(BULK_GENE_GAIN_ARTIFACTS)} when set: {BULK_GENE_GAIN!r}")
if BULK_GENE_GAIN is not None and PERT_SCALE_BALLAST is None:
    raise SystemExit("VCC_BULK_GENE_GAIN is a bulk-only arm and needs VCC_PERT_SCALE_BALLAST (welded it reaches the DE channel)")


def apply_bulk_gene_gain(effects: np.ndarray, perts, genes) -> tuple[np.ndarray, dict]:
    """Multiply every bulk row by the pinned per-gene gain m, then restore the non-panel Frobenius norm."""
    import hashlib as _hashlib  # noqa: PLC0415
    path, pinned = BULK_GENE_GAIN_ARTIFACTS[str(BULK_GENE_GAIN)]
    got = _hashlib.sha256(path.read_bytes()).hexdigest()
    if got != pinned:
        raise SystemExit(f"VCC_BULK_GENE_GAIN artifact sha256 {got} != pinned {pinned}")
    with np.load(path, allow_pickle=False) as z:
        if [str(g) for g in z["genes"]] != [str(g) for g in genes]:
            raise SystemExit("VCC_BULK_GENE_GAIN gene axis differs from the official gene order")
        m = np.asarray(z["m"], dtype=np.float64)
    if not np.all(np.isfinite(m)):
        raise SystemExit("VCC_BULK_GENE_GAIN multiplier is not finite")
    E = np.asarray(effects, dtype=np.float64)
    panel = set(map(str, perts))
    nonpanel = np.asarray([str(g) not in panel for g in genes], dtype=bool)
    before = float(np.sqrt((E[:, nonpanel] ** 2).sum()))
    out = E * m[None, :]
    after = float(np.sqrt((out[:, nonpanel] ** 2).sum()))
    out *= before / after
    audit = {"mode": str(BULK_GENE_GAIN), "artifact": str(path.relative_to(ROOT)), "sha256": got,
             "energy_constant": before / after, "m_quantiles": [float(v) for v in np.quantile(m, [0, .01, .5, .99, 1])],
             "scaled_sha256": _hashlib.sha256(np.ascontiguousarray(out).tobytes()).hexdigest()}
    print(f"[bulk-gene-gain] {BULK_GENE_GAIN} sha={got[:16]} energy constant {before / after:.6f}", flush=True)
    return out, audit

'''
sub("\n\nDE_SIGN_PRIOR_BY_TARGET: dict[str, np.ndarray] = {}   # filled by de_sign_prior_full for the per-target knobs",
    KNOB + "\nDE_SIGN_PRIOR_BY_TARGET: dict[str, np.ndarray] = {}   # filled by de_sign_prior_full for the per-target knobs",
    "knob")
sub('''        if PERT_SCALE_BALLAST is not None:
            effect_audit["pert_scale"]["ballast"] = {"lb": BALLAST_LB, "kmax": BALLAST_KMAX}''',
    '''        if PERT_SCALE_BALLAST is not None:
            effect_audit["pert_scale"]["ballast"] = {"lb": BALLAST_LB, "kmax": BALLAST_KMAX}
        if BULK_GENE_GAIN is not None:
            effects, effect_audit["bulk_gene_gain"] = apply_bulk_gene_gain(effects, spec.perts, spec.genes)''',
    "run")
sub('''            **({} if PERT_SCALE_BALLAST is None else {"pert_scale_ballast": {"lb": BALLAST_LB, "kmax": BALLAST_KMAX}}),''',
    '''            **({} if PERT_SCALE_BALLAST is None else {"pert_scale_ballast": {"lb": BALLAST_LB, "kmax": BALLAST_KMAX}}),
            **({} if BULK_GENE_GAIN is None else {"bulk_gene_gain": BULK_GENE_GAIN, "bulk_gene_gain_sha256": BULK_GENE_GAIN_ARTIFACTS[BULK_GENE_GAIN][1]}),''',
    "config")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
