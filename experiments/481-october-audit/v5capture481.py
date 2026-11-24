"""481: the pinned CD4 deflation basis V5_abc on a new panel. The builder records `captured_frac_of_centred` (the share of
the day's centred CD4 residual energy inside the A/B/C-fitted 5-dim basis); DEF-RUNBOOK: 0.10-0.14 measured, outside
0.05-0.19 is unmeasured territory. Recomputed with lib438's verbatim `deflate_pinned` on the complete-mean3 CD4 residual
(plolo.cd4_for's arithmetic, CD4 banks from a given root). KA: 453's and 465's metas record the builder's own value.
usage: uv run python v5capture481.py -> work/v5capture481.json"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/438-october-stack-repricing"))
import lib438 as L
SP = Path("/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/ad434c9f-b5f8-46bb-a96a-97fdd8cb3a3b/scratchpad")
E062 = L.E062
def cd4_for(targets, root: Path):
    banks = [E062.Bank.load(n, root / n, require_cd4_scale=True) for n in L.plolo.CD4_NAMES]
    prep = E062.prepare_cd4_mean(banks, targets, L.OFFICIAL, set(targets), "gene_norm", chunk_rows=256, max_rss_gib=3.5)
    complete = np.asarray(prep.target_condition_count) == len(L.plolo.CD4_NAMES)
    me = np.asarray(prep.mean_effect, dtype=np.float64).copy(); votes = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    me[~complete] = 0.0; votes[~complete] = 0
    return E062.CD4Prepared(mode=prep.mode, mean_effect=me, gene_votes=votes,
                            target_condition_count=np.asarray(prep.target_condition_count, dtype=np.uint8).copy(),
                            condition_details=prep.condition_details)
V = L.load_v5()
cases = {"abc (453 banks)": ("abc", ROOT / "experiments/453-abc-dress-rehearsal/work/run/cd4_compact",
                             ROOT / "experiments/453-abc-dress-rehearsal/work/run/453dr_final.h5ad.meta.json"),
         "efr (465 banks)": ("efr", ROOT / "experiments/465-both-switches-efr/work/run/cd4_compact",
                             ROOT / "experiments/465-both-switches-efr/work/run/465efr_final.h5ad.meta.json"),
         "2025 H1 panel (stand-in)": ("h1_2025", SP / "cd4_h1arc", None)}
out = {}
for name, (pn, root, meta) in cases.items():
    targets = [str(t) for t in pd.read_csv(HERE / "work/panels" / f"{pn}.csv")["target_gene"]]
    cd4 = cd4_for(targets, root)
    _, aud = L.deflate_pinned(cd4, V)
    rec = {k: round(float(v), 5) for k, v in aud.items() if isinstance(v, (int, float))}
    if meta is not None:
        rec["builder_meta_captured"] = json.load(open(meta))["effects"]["cd4"]["deflation"]["captured_frac_of_centred"]
    out[name] = rec
    print(name, rec, flush=True)
(HERE / "work/v5capture481.json").write_text(json.dumps(out, indent=1))
