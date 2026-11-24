"""195 — one prepared LODO fold, with every score computed in the fold's eigenbasis."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from fitlib import BatchStats, make_batches

HERE = Path(__file__).resolve().parent
PREP = HERE / "work" / "prep"


def _pds(cos: np.ndarray, lab: np.ndarray) -> np.ndarray:
    sd = cos[np.arange(cos.shape[0]), lab][:, None]
    return (((cos < sd).sum(1) + 0.5 * ((cos == sd).sum(1) - 1)) / (cos.shape[1] - 1))


class Fold:
    def __init__(self, held: str, dev_frac: float = 0.2, seed: int = 195):
        self.held = held
        self.meta = json.loads((PREP / f"{held}.meta.json").read_text())
        z = np.load(PREP / f"{held}.panel.npz")
        self.evals = z["evals"]
        self.panel_T = z["T"]                                    # (259, G) float64, eigen coords
        others = self.meta["others"]
        self.panel_Y = 0.5 * (z[f"p_{others[0]}"] + z[f"p_{others[1]}"])   # equal weights
        self.G = self.evals.shape[0]

        self.dirs = []
        rng = np.random.default_rng(seed)
        for name, d in self.meta["directions"].items():
            Yq = np.load(PREP / f"{d['tag']}.Yq.npy", mmap_mode="r")
            Yt = np.load(PREP / f"{d['tag']}.Yt.npy", mmap_mode="r")
            n = Yq.shape[0]
            perm = rng.permutation(n)
            ndev = int(round(dev_frac * n))
            self.dirs.append(dict(name=name, Yq=Yq, Yt=Yt,
                                  dev=np.sort(perm[:ndev]), fit=np.sort(perm[ndev:])))

    # ---------- panel (the held LINE: the number that decides everything) ----------
    def panel_scores(self, gamma: float, M: np.ndarray | None = None, k: int = 0,
                     lam: float = 0.0) -> tuple[float, float]:
        w = np.ones(self.G) if gamma == 0.0 else (self.evals + lam) ** (-gamma / 2.0)
        Z = self.panel_Y * w
        if M is not None and k:
            Z = Z.copy()
            Z[:, -k:] = Z[:, -k:] + Z[:, -k:] @ M
        Zn = Z / np.linalg.norm(Z, axis=1, keepdims=True)
        Tn = self.panel_T / np.linalg.norm(self.panel_T, axis=1, keepdims=True)
        cos = Zn @ Tn.T
        lab = np.arange(cos.shape[0])
        return float(_pds(cos, lab).mean()), float(np.diag(cos).mean())

    # ---------- training / dev batches ----------
    def batches(self, split: str, gamma: float, k: int, batch: int, seed: int,
                lam: float = 0.0, limit: int | None = None) -> list[BatchStats]:
        w = (np.ones(self.G) if gamma == 0.0
             else (self.evals + lam) ** (-gamma / 2.0)).astype(np.float32)
        rng = np.random.default_rng(seed)
        out: list[BatchStats] = []
        for d in self.dirs:
            rows = d[split]
            for b in make_batches(len(rows), batch, rng):
                idx = rows[b]
                Z = np.asarray(d["Yq"][idx], np.float64) * w
                T = np.asarray(d["Yt"][idx], np.float64)
                bs = BatchStats.__new__(BatchStats)
                bs.Zp = Z[:, -k:].copy() if k else np.zeros((len(idx), 0))
                bs.Tb = T[:, -k:].copy() if k else np.zeros((len(idx), 0))
                bs.R = Z @ T.T - (bs.Zp @ bs.Tb.T if k else 0.0)
                bs.c = np.einsum("ij,ij->i", Z, Z) - (np.einsum("ij,ij->i", bs.Zp, bs.Zp)
                                                      if k else 0.0)
                bs.tn = np.linalg.norm(T, axis=1)
                bs.B = len(idx)
                out.append(bs)
                if limit and len(out) >= limit:
                    return out
        return out

    def batch_pds(self, bss: list[BatchStats], M: np.ndarray | None, k: int) -> float:
        vals = []
        for bs in bss:
            G_ = bs.Zp if M is None else bs.Zp + bs.Zp @ M
            num = bs.R + (G_ @ bs.Tb.T if k else 0.0)
            q = np.sqrt(bs.c + (np.einsum("ij,ij->i", G_, G_) if k else 0.0))
            cos = num / (q[:, None] * bs.tn[None, :])
            vals.append(_pds(cos, np.arange(bs.B)).mean())
        return float(np.mean(vals))
