"""612b: the fold check in a DELIVERY-FAITHFUL form. The real emitter delivers base_lfc (key-signed design) + sign_hook x MDE
add on uncalled genes; the hook only sets the ADD's direction (diag612: it flips 2-6 % of signs and mostly shrinks |lfc|).
480's `small_base+mde` construction is exactly that shape: P = BETA * cons (BETA = 480's median |our base| / median |cons|
on H1 gate cells) + s * add. Here s = key sign (no hook) / CLEAN4 sign (plain hook) / key-gated (CLEAN4 on the weak half,
key sign on the strong half). Δ = mean(ref - X) per target bootstrap, positive = X better."""
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "experiments/480-nmae-instrument")); sys.path.insert(0, str(ROOT / "experiments/537-arm-pipeline"))
import fold480 as F
from common480 import boot_delta, load_donor_rows, nmae_rows, load_h1
import foldsign537 as FS
C4ALL = ["k562_gwps", "k562_essential", "nadig_jurkat", "nadig_hepg2"]
C4B = {"h1": C4ALL, "orion_hct116": C4ALL, "orion_hek293t": C4ALL, "k562_gwps": ["nadig_jurkat", "nadig_hepg2"]}
h = load_h1(); G = np.isfinite(h["T"])
L0 = np.nan_to_num(np.asarray(np.load(ROOT / "experiments/426-mde-magnitude/work/h1/stats-base_s0_bps.npz")["lfc"], float))
c3 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, h["targets"], h["genes"], pc_cpm=1.0)[0] for d in F.DON]), 0))
BETA = float(np.median(np.abs(L0[G])) / np.median(np.abs(c3[G])))
print(f"BETA {BETA:.4f} (480's construction)")
for name in C4B:
    S = F.load_surface(name)
    T, cnt, genes, targets = S["T"], S["cnt"], S["genes"], S["targets"]
    cons0 = np.nan_to_num(np.nanmean(np.stack([load_donor_rows(d, targets, genes, pc_cpm=1.0)[0] for d in S["donors"]]), 0))
    P4 = FS.clean4_P(genes, C4B[name])
    s4 = np.where(np.isfinite(P4) & (np.abs(P4 - 0.5) >= 0.05), np.where(P4 > 0.5, -1.0, 1.0), 0.0)
    S4 = np.broadcast_to(s4, T.shape); sc = np.sign(cons0); add = F.mde_add(cnt)[None, :]
    weak = np.zeros(T.shape, bool)
    for i in range(T.shape[0]):
        a = np.abs(cons0[i]); nz = a > 0
        weak[i] = (~nz) | (a < (np.quantile(a[nz], 0.5) if nz.any() else 0.0))
    base = BETA * cons0
    s_plain = np.where(S4 != 0, S4, sc)                    # the hook: CLEAN4 where covered, else the delivered (key) sign
    s_kg = np.where(weak & (S4 != 0), S4, sc)              # key-gated
    P = {"nohook": base + sc * add, "plain": base + s_plain * add, "kg": base + s_kg * add}
    rows = {k: nmae_rows(v, T) for k, v in P.items()}
    d1 = boot_delta(rows["plain"], rows["nohook"]); d2 = boot_delta(rows["kg"], rows["nohook"]); d3 = boot_delta(rows["kg"], rows["plain"])
    print(f"{name:14s} nmae nohook {np.nanmean(rows['nohook']):.5f} plain {np.nanmean(rows['plain']):.5f} kg {np.nanmean(rows['kg']):.5f} | "
          f"plain vs nohook {d1[0]:+.5f} [{d1[1]:+.5f},{d1[2]:+.5f}]  kg vs nohook {d2[0]:+.5f} [{d2[1]:+.5f},{d2[2]:+.5f}]  "
          f"kg vs plain {d3[0]:+.5f} [{d3[1]:+.5f},{d3[2]:+.5f}]")
