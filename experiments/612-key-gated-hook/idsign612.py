#!/usr/bin/env python3
"""612: the identity-augmented low-half sign for the H1 surface in 622's RECOMMENDED form ([7772dc], 622/RESULTS.md):
S = z(P_C4 - 0.5) + lam z(I), positive -> DOWN; lam 0.5; I = log1p(ctrl CPM of h1) - mean over the FLEX references only
(A, B, C challenge controls; never public 3' profiles), z over the genes all of them measure; identity alone where CLEAN4
is undefined. Exactly hyb622.py's s_low for name='h1', refs flex. Written on 323's reference gene axis (18,080) for run612's
`idkg` mode, which passes it on each target's weak-key half (622's key-gated delivery). < 2 GB, seconds."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument")); sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
import os
os.chdir(ROOT / "experiments/480-nmae-instrument")
import fold480 as F          # noqa: E402
import foldsign537 as FS     # noqa: E402
import hyb622 as H           # noqa: E402  (profiles, identity: read-only)
LAM = 0.5
prof = H.profiles()
S = F.load_surface("h1")
genes = [str(g) for g in S["genes"]]
P4 = FS.clean4_P(np.array(genes), ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"])
zc = np.where(np.isfinite(P4), (P4 - 0.5) / np.nanstd(P4 - 0.5), np.nan)
zi = H.identity(prof, "h1", ["ctxA", "ctxB", "ctxC"], genes)
sc_low = np.where(np.isfinite(zc), zc + LAM * np.nan_to_num(zi), np.nan_to_num(zi))
s_low = np.where(np.isfinite(sc_low) & (sc_low != 0), np.where(sc_low > 0, -1.0, 1.0), 0.0)
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
GF = ref["genes"].astype(str)
pos = {g: i for i, g in enumerate(genes)}
sign = np.array([s_low[pos[g]] if g in pos else 0.0 for g in GF])
score = np.array([sc_low[pos[g]] if g in pos else np.nan for g in GF])
np.savez_compressed(ROOT / "experiments/612-key-gated-hook/work/idsign612.npz", genes=GF, sign=sign, score=score,
                    lam=LAM, refs=np.array(["ctxA", "ctxB", "ctxC"]), form=np.array("622 flex lam0.5"))
c4s = np.where(np.isfinite(P4), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
both = (c4s != 0) & (s_low != 0)
print(f"622 form: lam {LAM}, refs A/B/C; h1 axis {len(genes)}: signed {int((s_low != 0).sum())} (CLEAN4 signed "
      f"{int((c4s != 0).sum())}); identity-only genes {int(((c4s == 0) & (s_low != 0)).sum())}; agreement with CLEAN4 "
      f"{float((s_low[both] == c4s[both]).mean()):.3f}; on ref axis signed {int((sign != 0).sum())}")
