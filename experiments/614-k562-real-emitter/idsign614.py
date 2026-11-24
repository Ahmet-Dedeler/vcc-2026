#!/usr/bin/env python3
"""614: the hook-form signs on the K562 harness (ungated forms, as idall3 / the plain hook ship).
  c4h    = CLEAN4 (Nadig Jurkat + HepG2 ONLY: both K562 banks dropped) at tau 0.05: -1 where P_down > 1/2, else +1.
  idall3 = 604 v3's rule: score = z(P_C4 - 1/2) + 0.5 z(I) on CLEAN4-covered genes (zi 0 where identity is undefined),
           sign = -sign(score) for down; uncovered genes 0 (the emitter keeps the delivered sign).
           I = log1p(K562 control mean CPM) - mean log1p(the Flex references A/B/C + H1), z over the genes with a value.
Only K562's own controls, the Flex references and the public CLEAN4 banks are used; no K562 truth is read.
Output: work/idsign614.npz (genes = the reference's full axis; c4h, idall3)."""
import sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/639-identity-axis-lolo"))
sys.path.insert(0, str(ROOT / "experiments/636-lolo-sign-instrument"))
import decomp639 as DC  # noqa: E402
r = np.load(HERE / "work/ref_k562s/reference.npz", allow_pickle=False)
genes = r["genes"].astype(str)
keep = np.asarray(r["keep"], bool)
erl = np.asarray(r["emit_ref_lin"], float)
cpm = np.full(len(genes), np.nan)
cpm[keep] = erl if erl.shape[0] == int(keep.sum()) else erl[keep]
cp = np.load(ROOT / "experiments/371-consensus-block-prior-lolo/work/consensus_prior.npz", allow_pickle=True)
names = list(cp["bank_names"].astype(str))
P = cp["P_banks"].astype(float)[[names.index(b) for b in ("nadig_jurkat", "nadig_hepg2")]].mean(0)
pos = {g: j for j, g in enumerate(cp["genes"].astype(str))}
col = np.array([pos.get(g, -1) for g in genes])
P4 = np.where(col >= 0, P[np.maximum(col, 0)], np.nan)
c4h = np.where(np.isfinite(P4) & (np.abs(P4 - 0.5) >= 0.05), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
prof = DC.flex_profiles()
own = {g: v for g, v in zip(genes, cpm) if np.isfinite(v) and not g.startswith("SUBGATE_")}
Iv = DC.identity(own, [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")], list(genes))
zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
cov = np.isfinite(P4)
zc = np.where(cov, (P4 - 0.5) / np.nanstd((P4 - 0.5)[cov]), np.nan)
sc = np.where(cov, zc + 0.5 * np.nan_to_num(zi), np.nan)
idall3 = np.where(np.isfinite(sc) & (sc != 0), np.where(sc > 0, -1.0, 1.0), 0.0)
np.savez(HERE / "work/idsign614.npz", genes=genes, c4h=c4h, idall3=idall3)
both = (c4h != 0) & (idall3 != 0)
print(f"c4h signs {int((c4h != 0).sum())} (down {int((c4h < 0).sum())}); idall3 signs {int((idall3 != 0).sum())} "
      f"(CLEAN4-covered {int(cov.sum())}); agree where both sign {float((c4h[both] == idall3[both]).mean()):.3f}; "
      f"identity defined on {int(np.isfinite(Iv).sum())} genes")
