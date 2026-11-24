"""517 debug: why does the complement hook move fid on the Arc base? Which called genes change sign?"""
import sys, numpy as np
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import price517 as P
import recon_h1 as RC
zr = np.load(RC.REF, allow_pickle=False); genes_full = zr["genes"]; keep = np.asarray(zr["keep"], bool); T = zr["targets"].astype(str); zr.close()
SRC, info = P.load_sources(genes_full, keep)
arc_block = P.sign_from_P(SRC["arc"]); c4_block = RC.prior_signs_full(genes_full)[keep]
for base, pr in (("c4", lambda i, t: c4_block), ("arc", lambda i, t: arc_block)):
    ctx = P.Ctx(base, 0, pr)
    s = np.where(np.isfinite(SRC["arc"]) & (np.abs(SRC["arc"] - 0.5) >= 0.05), np.where(SRC["arc"] > 0.5, -1.0, 1.0), 0.0)
    lfc_new, rec = ctx.hook(lambda i, t: s)
    called = ctx.st["q"] < 0.05
    ch = called & (np.sign(lfc_new) != np.sign(ctx.lfc_del))
    tot = {"head": 0, "blk": 0, "comp": 0}
    for B in ctx.blocks:
        i = B["i"]
        c = ch[i]
        tot["head"] += int((c & B["o"]["head"]).sum()); tot["blk"] += int((c & B["blk"]).sum()); tot["comp"] += int((c & B["comp"]).sum())
    print(base, "called genes whose sign flips:", int(ch.sum()), tot, "called per target", float(called.sum(1).mean()))
    # among flipped called comp genes: truth agreement before/after
    lr, qr = ctx.lfc_real, ctx.q_real
    adj = np.isfinite(lr) & (lr != 0)
    m = ch & adj
    print("  flipped & adjudicable:", int(m.sum()), "now correct:", int((np.sign(lfc_new[m]) == np.sign(lr[m])).sum()),
          "before correct:", int((np.sign(ctx.lfc_del[m]) == np.sign(lr[m])).sum()))
    # delivered |lfc| of the flipped genes
    print("  median |lfc_del| of flipped:", float(np.median(np.abs(ctx.lfc_del[ch]))) if ch.any() else None)
