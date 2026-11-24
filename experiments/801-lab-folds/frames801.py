"""801 -- bind a WORKING panel once (lib801.bind with dump801's spec) so its frame caches (work/lib/<L>/bctrl-*.npy) exist
before the three fold jobs of that panel start in parallel (they would otherwise race to write the same cache files).
bind asserts the synthetic bundle's frames equal A/B/C's.

    uv run python experiments/801-lab-folds/frames801.py <L>
"""
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib801 as L  # noqa: E402

k = sys.argv[1]
if k == "efp":                              # the KA panel: lib801's registry spec
    ctx = L.bind("efp", tag="efp")
    print(f"efp: frames bound and cached in {ctx.work}")
    raise SystemExit(0)
assert k not in ("L7", "L8", "sealed"), k                 # any working / kit panel; never the sealed reserve
W = HERE / "work" / k
spec = dict(bundle=str((W / "bundle").resolve()), rules="def", cd4_root=str((W / "cd4/cd4_compact").resolve()),
            full_op=str(Path("none").resolve()), folds=["k562_gwps", "orion_hct116", "orion_hek293t"],
            synthetic=not (W / "bundle/frames_src.txt").exists(),          # kit --frames: the bundle's own controls
            meta=None, fold_ops={})
ctx = L.bind(k, spec=spec, tag=k)
print(f"{k}: frames bound and cached in {ctx.work} ({len(ctx.panel)} targets, contexts {ctx.contexts})")
