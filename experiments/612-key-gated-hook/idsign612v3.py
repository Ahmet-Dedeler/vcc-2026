#!/usr/bin/env python3
"""612 (v3 form, [d4cd5f] 10-02 ~13:55Z): the H1 hook source for KGID in 604 v3's rule. score = z(P_C4 − ½) + 0.5·z(I) on
CLEAN4-COVERED genes only (a covered gene with no identity value gets zc alone), side = sign(score), and NO prior on
CLEAN4-uncovered genes (sign 0 = the emitter keeps the delivered sign). It is idsign612.py with the identity-only
extension removed: 639d measured that extension wrong on uncalled genes of cancer lines (HEK293T 0.27–0.34; HCT116 0.49),
which is where the hook acts. Written on 323's reference gene axis for run612's `idkg3` mode. < 2 GB, seconds."""
import os
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument")); sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
sys.path.insert(0, str(ROOT / "experiments/622-identity-track"))
os.chdir(ROOT / "experiments/480-nmae-instrument")
import fold480 as F          # noqa: E402
import foldsign537 as FS     # noqa: E402
import hyb622 as H           # noqa: E402  (profiles, identity: read-only, as idsign612)
LAM = 0.5
prof = H.profiles()
S = F.load_surface("h1")
genes = [str(g) for g in S["genes"]]
P4 = FS.clean4_P(np.array(genes), ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"])
zc = np.where(np.isfinite(P4), (P4 - 0.5) / np.nanstd(P4 - 0.5), np.nan)
zi = H.identity(prof, "h1", ["ctxA", "ctxB", "ctxC"], genes)
sc = np.where(np.isfinite(zc), zc + LAM * np.nan_to_num(zi), np.nan)          # v3: no identity-only extension
s3 = np.where(np.isfinite(sc) & (sc != 0), np.where(sc > 0, -1.0, 1.0), 0.0)
ref = np.load(ROOT / "experiments/323-current-family-validation/work/reference.npz", allow_pickle=False)
GF = ref["genes"].astype(str)
pos = {g: i for i, g in enumerate(genes)}
sign = np.array([s3[pos[g]] if g in pos else 0.0 for g in GF])
score = np.array([sc[pos[g]] if g in pos else np.nan for g in GF])
np.savez_compressed(ROOT / "experiments/612-key-gated-hook/work/idsign612v3.npz", genes=GF, sign=sign, score=score,
                    lam=LAM, refs=np.array(["ctxA", "ctxB", "ctxC"]), form=np.array("604 v3 rule: flex lam0.5, CLEAN4-covered only"))
old = np.load(ROOT / "experiments/612-key-gated-hook/work/idsign612.npz")["sign"]
c4s = np.where(np.isfinite(P4), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
both = (c4s != 0) & (s3 != 0)
print(f"v3 form: lam {LAM}; h1 axis {len(genes)}: signed {int((s3 != 0).sum())} (CLEAN4 signed {int((c4s != 0).sum())}); "
      f"agreement with CLEAN4 {float((s3[both] == c4s[both]).mean()):.3f}; ref axis signed {int((sign != 0).sum())}; "
      f"vs idsign612 (622 form): equal where v3 signs {bool(np.array_equal(sign[sign != 0], old[sign != 0]))}, "
      f"extension genes dropped {int(((old != 0) & (sign == 0)).sum())}")
