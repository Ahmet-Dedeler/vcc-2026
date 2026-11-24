"""285 -- the Orion / K562-GWPS REFERENCE rows, read from 267's cached half-truth arrays: the batch (d) and
random-cell (a) halves at the SAME cells per side, paired per target.  267 is the one place in the archive
where a batch split reads clearly below a random split (HEK293T); this row makes that contrast paired and
gives it an interval.  No new pass over the 118 GB of Orion parquet is made -- see RESULTS section 7.
"""
import sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib285 as L
P = L.Logger(L.WORK / "orion_285.log")
E267 = L.ROOT / "experiments/267-replicate-validated-oracle/work/half"
out = {}
for line in ("k562_gwps", "orion_hct116", "orion_hek293t"):
    z = np.load(E267 / f"truth-{line}.npz", allow_pickle=True)
    b = np.asarray(z["rel|batch|cos_resid"], dtype=np.float64)
    r = np.asarray(z["rel|rand|cos_resid"], dtype=np.float64)
    g = np.asarray(z["rel|guide|cos_resid"], dtype=np.float64)
    nb = np.asarray(z["rel|batch|ncells_A"], dtype=np.float64); nr = np.asarray(z["rel|rand|ncells_A"], dtype=np.float64)
    ng = np.asarray(z["rel|guide|ncells_A"], dtype=np.float64)
    okg = np.isfinite(g)
    out[line] = dict(n=int(np.isfinite(b).sum()),
                     cells_batch=float(np.median(nb[np.isfinite(nb)])), cells_rand=float(np.median(nr[np.isfinite(nr)])),
                     batch=dict(**L.stats(b), boot=L.boot_median(b)), rand=dict(**L.stats(r), boot=L.boot_median(r)),
                     rand_minus_batch=L.boot_median_diff(r, b),
                     guide=dict(**L.stats(g[okg]), boot=L.boot_median(g[okg]), n_guide=int(okg.sum()),
                                cells=float(np.median(ng[okg]))),
                     rand_minus_guide_on_guide_targets=L.boot_median_diff(r[okg], g[okg]))
    P(f"{line}: {out[line]['n']} targets; RANDOM cell halves ({out[line]['cells_rand']:.0f} cells/side) "
      f"{L.fmt(out[line]['rand']['boot'])}  BATCH halves ({out[line]['cells_batch']:.0f}) "
      f"{L.fmt(out[line]['batch']['boot'])}  RANDOM − BATCH {L.fmt(out[line]['rand_minus_batch'])}")
    P(f"    GUIDE halves on {int(okg.sum())} two-construct targets ({out[line]['guide']['cells']:.0f} cells/side) "
      f"{L.fmt(out[line]['guide']['boot'])};  RANDOM − GUIDE on those targets "
      f"{L.fmt(out[line]['rand_minus_guide_on_guide_targets'])}")
L.Logger.dump(P, out, L.WORK / "orion_285.json"); P.save()
