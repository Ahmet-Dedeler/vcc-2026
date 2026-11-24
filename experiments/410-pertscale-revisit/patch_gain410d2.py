"""410d2: render-gated energy matching for VCC_BULK_GENE_GAIN (fixes 414's URGENT-410d finding).

410d normalised the gained bulk to its EFFECT-space Frobenius norm over non-panel columns. The render zeroes every gene
below the context's 10-CPM control gate, and 228's gain favours expressed genes, so the energy the scorer reads rose
~+14..20 % (an amplitude step past context A's mse cliff). This adds knob value `g228reg_t025rg`: same artifact
(da937cfa), normalised so the energy on each context's RENDERED genes (control CPM > min_cpm and b_ctrl > 0, the
profile_from_effect gate, non-panel columns), summed over the bundle's contexts, equals the pre-gain bulk's.
`g228reg_t025` (410d as built) is kept unchanged for provenance. Apply only when no build has the builder pinned.
usage: uv run python patch_gain410d2.py [--write]      (PATCH_TARGET=<path> to dry-run on a copy)
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


sub('''    "g228reg_t025": (
        ROOT / "experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz",
        "da937cfa69e0730b215b8780c5ca0505d7d4acd6f0ae957891ca74731971fd9f",
    ),
}''', '''    "g228reg_t025": (
        ROOT / "experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz",
        "da937cfa69e0730b215b8780c5ca0505d7d4acd6f0ae957891ca74731971fd9f",
    ),
    # 410d2: the same artifact, energy matched on each context's RENDERED genes (414, URGENT-410d.md)
    "g228reg_t025rg": (
        ROOT / "experiments/410-pertscale-revisit/work/bulk_gene_gain_g228reg_t025.npz",
        "da937cfa69e0730b215b8780c5ca0505d7d4acd6f0ae957891ca74731971fd9f",
    ),
}''', "artifact")

sub('''def apply_bulk_gene_gain(effects: np.ndarray, perts, genes) -> tuple[np.ndarray, dict]:
    """Multiply every bulk row by the pinned per-gene gain m, then restore the non-panel Frobenius norm."""''',
    '''def control_render_gate(bundle: Path, context: str, genes, min_cpm: float) -> np.ndarray:
    """profile_from_effect's gate for one context: control CPM > min_cpm and b_ctrl > 0, from the summed control counts."""
    import h5py  # noqa: PLC0415
    with h5py.File(Path(bundle) / f"context_{context}.h5ad", "r") as h:
        idx = h["X/indices"][:]
        dat = h["X/data"][:].astype(np.float64)
    totals = np.bincount(idx, weights=dat, minlength=len(genes))
    if totals.shape[0] != len(genes):
        raise SystemExit(f"context {context}: control gene axis {totals.shape[0]} != {len(genes)}")
    cpm = 1e6 * totals / totals.sum()
    return (cpm > float(min_cpm)) & (totals > 0)


def apply_bulk_gene_gain(effects: np.ndarray, perts, genes, gates=None) -> tuple[np.ndarray, dict]:
    """Multiply every bulk row by the pinned per-gene gain m, then restore the pre-gain energy: over non-panel columns
    (`g228reg_t025`, 410d as built), or -- with `gates` (one boolean gene mask per context) -- over each context's
    RENDERED non-panel genes summed over contexts (`g228reg_t025rg`, 410d2)."""''', "helper")

sub('''    E = np.asarray(effects, dtype=np.float64)
    panel = set(map(str, perts))
    nonpanel = np.asarray([str(g) not in panel for g in genes], dtype=bool)
    before = float(np.sqrt((E[:, nonpanel] ** 2).sum()))
    out = E * m[None, :]
    after = float(np.sqrt((out[:, nonpanel] ** 2).sum()))
    out *= before / after''', '''    E = np.asarray(effects, dtype=np.float64)
    panel = set(map(str, perts))
    nonpanel = np.asarray([str(g) not in panel for g in genes], dtype=bool)
    out = E * m[None, :]
    if gates is None:
        before = float(np.sqrt((E[:, nonpanel] ** 2).sum()))
        after = float(np.sqrt((out[:, nonpanel] ** 2).sum()))
    else:
        cols = [np.asarray(g, dtype=bool) & nonpanel for g in gates]
        before = float(np.sqrt(sum((E[:, c] ** 2).sum() for c in cols)))
        after = float(np.sqrt(sum((out[:, c] ** 2).sum() for c in cols)))
    out *= before / after''', "normalise")

sub('''    audit = {"mode": str(BULK_GENE_GAIN), "artifact": str(path.relative_to(ROOT)), "sha256": got,''',
    '''    audit = {"mode": str(BULK_GENE_GAIN), "artifact": str(path.relative_to(ROOT)), "sha256": got,
             "normalisation": "non-panel effect-space" if gates is None else "rendered non-panel genes, summed over contexts",''',
    "audit")

sub('''            effects, effect_audit["bulk_gene_gain"] = apply_bulk_gene_gain(effects, spec.perts, spec.genes)''',
    '''            gates = None
            if str(BULK_GENE_GAIN).endswith("rg"):
                gates = [control_render_gate(bundle, str(c), spec.genes, min_cpm) for c in spec.contexts]
            effects, effect_audit["bulk_gene_gain"] = apply_bulk_gene_gain(effects, spec.perts, spec.genes, gates=gates)''',
    "run")

sys.stdout.writelines(difflib.unified_diff(src.splitlines(True), out.splitlines(True), "a/" + P.name, "b/" + P.name, n=1))
if "--write" in sys.argv:
    P.write_text(out)
    print("\nWROTE", P)
