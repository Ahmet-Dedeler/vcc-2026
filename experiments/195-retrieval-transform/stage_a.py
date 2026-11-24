"""195 Stage A — is the hand-chosen family already at the optimum of the direct objective?

Before fitting anything, ask the cheaper question the sweep never asked.  `gamma = 0.65` was
selected by maximising `pds` on the 259 held-LINE panel targets.  The proposal in this
experiment is to fit the operator on *donor* retrieval problems instead.  So: run the same
one-parameter family through the donor objective and see where that objective puts its
optimum.

Two outcomes, both informative.  If the donor-side optimum lands on 0.65, the training signal
is faithful and any fitted operator that fails to transfer has failed for capacity reasons.
If it lands somewhere else, the donor objective is the wrong teacher and no amount of
optimisation on it can be trusted -- which would by itself explain experiment 080.

Reported alongside: macro cosine, which 4.2 requires be allowed to fall.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from foldio import Fold  # noqa: E402
from harness import LODO  # noqa: E402

HERE = Path(__file__).resolve().parent
GAMMAS = [0.0, 0.2, 0.35, 0.45, 0.55, 0.65, 0.75, 0.9, 1.1, 1.3]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--out", default=str(HERE / "work" / "stage_a.json"))
    args = ap.parse_args()

    res = {"gammas": GAMMAS, "batch": args.batch, "folds": {}}
    for held in LODO:
        t0 = time.time()
        f = Fold(held)
        row = {"n_fit": int(sum(len(d["fit"]) for d in f.dirs)),
               "n_dev": int(sum(len(d["dev"]) for d in f.dirs)),
               "panel_pds": [], "panel_cos": [], "train_pds": [], "dev_pds": []}
        for g in GAMMAS:
            p, c = f.panel_scores(g)
            row["panel_pds"].append(p)
            row["panel_cos"].append(c)
            row["train_pds"].append(f.batch_pds(
                f.batches("fit", g, 0, args.batch, seed=1), None, 0))
            row["dev_pds"].append(f.batch_pds(
                f.batches("dev", g, 0, args.batch, seed=2), None, 0))
            print(f"  {held:14} gamma={g:<5g} panel_pds={p:.4f} panel_cos={c:.4f} "
                  f"train_pds={row['train_pds'][-1]:.4f} dev_pds={row['dev_pds'][-1]:.4f}",
                  flush=True)
        res["folds"][held] = row
        print(f"  ({time.time()-t0:.0f}s)", flush=True)

    macro = {kind: [float(np.mean([res["folds"][h][kind][i] for h in LODO]))
                    for i in range(len(GAMMAS))]
             for kind in ["panel_pds", "panel_cos", "train_pds", "dev_pds"]}
    res["macro"] = macro
    print("\n gamma | panel pds | panel cos | donor train pds | donor dev pds")
    for i, g in enumerate(GAMMAS):
        print(f" {g:<5g} |  {macro['panel_pds'][i]:.4f}   |  {macro['panel_cos'][i]:.4f}   |"
              f"     {macro['train_pds'][i]:.4f}      |    {macro['dev_pds'][i]:.4f}")
    for kind in ["panel_pds", "train_pds", "dev_pds"]:
        i = int(np.argmax(macro[kind]))
        print(f"argmax {kind:10} -> gamma {GAMMAS[i]}  ({macro[kind][i]:.4f})")
        res.setdefault("argmax", {})[kind] = GAMMAS[i]
    Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
