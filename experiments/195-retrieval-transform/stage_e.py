"""195 Stage E — where exactly does the fitted correction stop generalising?

Stage B shows a fit that gains on donor perturbations it has never seen and loses on the held
cell line.  Two different explanations produce that:

  (a) the correction is real but line-pair-specific -- it cannot cross a change of cell line;
  (b) the correction is real and line-general, and something about the panel is different.

Only (a) is testable without spending the held line, and it is testable because RPE1 and
Nadig-Jurkat are never scored.  Each fold therefore has four training lines.  Fit on the six
directions that avoid a training line X, then score the three directions whose GALLERY is X --
an unseen cell line, exactly the relationship the panel has to the fit.  Three rungs:

  dev        unseen perturbations, seen lines
  unseen-L   unseen gallery LINE, donor perturbations
  panel      unseen gallery line AND the 259 VCC panel targets  (the real number)
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fitlib import Adapter, BatchStats, batch_forward, make_batches  # noqa: E402
from foldio import Fold, _pds  # noqa: E402
from harness import LODO  # noqa: E402

HERE = Path(__file__).resolve().parent
PREP = HERE / "work" / "prep_multi"


class MultiFold(Fold):
    def __init__(self, held: str, dev_frac: float = 0.2, seed: int = 195):
        super().__init__(held, dev_frac, seed)
        self.mmeta = json.loads((PREP / f"{held}.meta.json").read_text())
        rng = np.random.default_rng(seed + 1)
        self.mdirs = []
        for name, d in self.mmeta["directions"].items():
            Yq = np.load(PREP / f"{d['tag']}.Yq.npy", mmap_mode="r")
            Yt = np.load(PREP / f"{d['tag']}.Yt.npy", mmap_mode="r")
            perm = rng.permutation(Yq.shape[0])
            ndev = int(round(dev_frac * Yq.shape[0]))
            self.mdirs.append(dict(name=name, Yq=Yq, Yt=Yt, tgt=d["target_line"],
                                   pred=d["query_line"], dev=np.sort(perm[:ndev]),
                                   fit=np.sort(perm[ndev:])))

    def mbatches(self, dirs, split: str, gamma: float, k: int, batch: int,
                 seed: int) -> list[BatchStats]:
        w = ((self.evals) ** (-gamma / 2.0)).astype(np.float32)
        rng = np.random.default_rng(seed)
        out = []
        for d in dirs:
            rows = d[split]
            for b in make_batches(len(rows), batch, rng):
                idx = rows[b]
                Z = np.asarray(d["Yq"][idx], np.float64) * w
                T = np.asarray(d["Yt"][idx], np.float64)
                bs = BatchStats.__new__(BatchStats)
                bs.Zp, bs.Tb = Z[:, -k:].copy(), T[:, -k:].copy()
                bs.R = Z @ T.T - bs.Zp @ bs.Tb.T
                bs.c = np.einsum("ij,ij->i", Z, Z) - np.einsum("ij,ij->i", bs.Zp, bs.Zp)
                bs.tn = np.linalg.norm(T, axis=1)
                bs.B = len(idx)
                out.append(bs)
        return out


def fit(bss, k, rank, epochs, lr, reg, tau, seed):
    ad = Adapter(k, rank, True, seed=seed)
    rng = np.random.default_rng(seed)
    order = np.arange(len(bss))
    for _ in range(epochs):
        rng.shuffle(order)
        for i in order:
            _, _, dM = batch_forward(bss[i], ad.M(), tau)
            ad.step(dM, lr, reg)
    return ad.M()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--k", type=int, default=256)
    ap.add_argument("--rank", type=int, default=8)
    ap.add_argument("--epochs", type=int, default=25)
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--out", default=str(HERE / "work" / "stage_e.json"))
    args = ap.parse_args()
    k, r, g = args.k, args.rank, 0.65
    res = {"cfg": vars(args), "folds": {}}

    for held in LODO:
        f = MultiFold(held)
        lines = f.mmeta["lines"]
        rows = {}
        for X in lines:
            fitd = [d for d in f.mdirs if d["tgt"] != X and d["pred"] != X]
            evald = [d for d in f.mdirs if d["tgt"] == X]
            b_fit = f.mbatches(fitd, "fit", g, k, args.batch, seed=21)
            b_dev = f.mbatches(fitd, "dev", g, k, args.batch, seed=22)
            b_uns = f.mbatches(evald, "dev", g, k, args.batch, seed=23)
            M = fit(b_fit, k, r, args.epochs, 0.01, 1e-3, 0.05, 195)
            base_panel, _ = f.panel_scores(g)
            row = {"n_fit_dirs": len(fitd), "n_fit_q": int(sum(b.B for b in b_fit)),
                   "n_unseen_q": int(sum(b.B for b in b_uns)),
                   "d_train": f.batch_pds(b_fit, M, k) - f.batch_pds(b_fit, None, k),
                   "d_dev": f.batch_pds(b_dev, M, k) - f.batch_pds(b_dev, None, k),
                   "d_unseen_line": f.batch_pds(b_uns, M, k) - f.batch_pds(b_uns, None, k),
                   "d_panel": f.panel_scores(g, M, k)[0] - base_panel}
            rows[X] = row
            print(f"  {held:14} hold training line {X:14} dirs={row['n_fit_dirs']:2d} "
                  f"q={row['n_fit_q']:6d}  d_train {row['d_train']:+.4f}  "
                  f"d_dev {row['d_dev']:+.4f}  d_unseenLINE {row['d_unseen_line']:+.4f}  "
                  f"d_panel {row['d_panel']:+.4f}", flush=True)

        # and the all-lines fit, which is what a deployable operator would actually use
        b_fit = f.mbatches(f.mdirs, "fit", g, k, args.batch, seed=21)
        b_dev = f.mbatches(f.mdirs, "dev", g, k, args.batch, seed=22)
        M = fit(b_fit, k, r, args.epochs, 0.01, 1e-3, 0.05, 195)
        base_panel, base_cos = f.panel_scores(g)
        p, c = f.panel_scores(g, M, k)
        allrow = {"n_fit_dirs": len(f.mdirs), "n_fit_q": int(sum(b.B for b in b_fit)),
                  "d_train": f.batch_pds(b_fit, M, k) - f.batch_pds(b_fit, None, k),
                  "d_dev": f.batch_pds(b_dev, M, k) - f.batch_pds(b_dev, None, k),
                  "d_panel": p - base_panel, "panel": p, "base_panel": base_panel,
                  "panel_cos": c, "base_cos": base_cos}
        print(f"  {held:14} ALL 4 lines              dirs={allrow['n_fit_dirs']:2d} "
              f"q={allrow['n_fit_q']:6d}  d_train {allrow['d_train']:+.4f}  "
              f"d_dev {allrow['d_dev']:+.4f}  {'':22} d_panel {allrow['d_panel']:+.4f}",
              flush=True)
        res["folds"][held] = {"lines": lines, "lolo": rows, "all_lines": allrow}
        Path(args.out).write_text(json.dumps(res, indent=2))

    keys = ["d_train", "d_dev", "d_unseen_line", "d_panel"]
    macro = {kk: float(np.mean([res["folds"][h]["lolo"][X][kk]
                                for h in LODO for X in res["folds"][h]["lines"]]))
             for kk in keys}
    res["macro_lolo"] = macro
    res["macro_all_lines"] = {kk: float(np.mean([res["folds"][h]["all_lines"][kk] for h in LODO]))
                              for kk in ["d_train", "d_dev", "d_panel"]}
    print("\n== inner leave-one-training-LINE-out, macro over 12 (fold, held line) cells ==")
    for kk in keys:
        print(f"   {kk:15} {macro[kk]:+.5f}")
    print("== fit on all four training lines ==")
    for kk, v in res["macro_all_lines"].items():
        print(f"   {kk:15} {v:+.5f}")
    Path(args.out).write_text(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
