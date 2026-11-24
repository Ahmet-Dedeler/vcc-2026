#!/usr/bin/env python3
"""613d: a pre-registered OUT-OF-SAMPLE test of a hook-class instrument against 496a's board read (10-05). Committed
before its first run and before 496a exists on the board.

The hook (VCC_DE_MDE_SIGN, rc_emit_k) changes only the direction of the MDE add on UNCALLED genes (not head, not block).
426's rule: on a gate gene with add < |lfc|, a correct add direction lowers |error| by add and a wrong one raises it.
Against the key-signed add (no hook), flipping the add on a gene where the key had a sign moves |error| by ±2·add, and
signing a gene the key left at 0 moves it by ±add. So the instrument needs no base magnitudes, only sign accuracy and the
add size (fold480.mde_add(cnt): c 0.1 / cap 0.1 / min count 0.3).
  predicted Δnmae(target) = Σ_cells Δ|err| / Σ|T| over the target's pool; then the mean over scored targets
  (cell-eval2's >= 10-gene rule; 521's pools via ds612b.pool_T).
Surfaces: 636's Inst (h1 and the 3 cancer folds). Uncalled = tested & ~head & ~block (decomp639.callset_for). Key sign =
block613c.signed_key (design key on h1, 3-donor consensus on the folds). Hook = CLEAN4 (clean4_lolo) at τ 0.05, i.e.
c4h / 496a's form.

Pre-registered, fixed before any number exists:
  RETRODICTION: the h1 cell must give c4h − base Δnmae < 0 at f 1 AND f 0.35. The real emitter has −0.0015 / −0.0021
    raw, 3 seeds (score612 c4h-base). If not, nothing below is registered.
  PREDICTION for 496a − 446bj2 on the board (10-05): sign(board Δnmae) = the sign of the MEAN of the 3 cancer folds'
    predicted Δnmae at f 0.35. Reported beside it: each fold's value at f 1 / 0.5 / 0.35, and the KGID form (identity
    v3 on weak-key uncalled genes) vs c4h, for the record.
  The read: if 496a's board nmae moves in the predicted direction (|Δ| > 0.0005 raw, outside the DE re-roll noise for a
  prior-class arm), the instrument earns a held-out retrodiction for hook-class levers. If it moves the other way, the
  instrument is closed for hooks. A |Δ| <= 0.0005 read is uninformative.
Output: work/block613d.out, work/block613d.json.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import block613 as B  # noqa: E402
import block613c as BC  # noqa: E402
import decomp639 as DC  # noqa: E402
import ds612b as DB  # noqa: E402
from inst636 import G480, SURF, Inst, clean4_lolo  # noqa: E402

sys.path.insert(0, str(B.ROOT / "experiments/480-nmae-instrument"))
import fold480 as F  # noqa: E402

REG = (1.0, 0.5, 0.35)


def pred_dnmae(hook_sign, key_sign, add, Tpool, uncalled):
    """Mean over scored targets of Σ Δ|err| / Σ|T| on the pool, 426's add rule, hook vs the key-signed add."""
    G = np.isfinite(Tpool)
    den = np.where(G, np.abs(Tpool), 0.0).sum(1)
    ok = (G.sum(1) >= 10) & (den > 0)
    tdn = np.sign(np.nan_to_num(Tpool))
    hs = np.broadcast_to(hook_sign, Tpool.shape)
    act = G & uncalled & (hs != 0)
    flip = act & (key_sign != 0) & (hs != key_sign)
    new = act & (key_sign == 0)
    corr = hs == tdn
    d = np.where(flip, np.where(corr, -2.0, 2.0) * add, 0.0) + np.where(new, np.where(corr, -1.0, 1.0) * add, 0.0)
    rows = d.sum(1) / np.where(den > 0, den, 1.0)
    return float(np.mean(rows[ok])), int(ok.sum()), int(flip.sum()), int(new.sum())


def main():
    I = Inst()
    prof = DC.flex_profiles()
    out, lines = {}, []
    for s in ("h1",) + B.CANCER:
        D = I.S[s]
        T, G, lfc, gate, tested = D["T"], D["G"], D["lfc"], D["gate"], D["tested"]
        z = np.load(G480 / SURF[s], allow_pickle=True)
        p = z["p"].astype(float)
        cnt = z["ctrl_cpm"].astype(float) * float(z["lib_med"]) / 1e6
        add = F.mde_add(cnt)[None, :]
        head, block = DC.callset_for(I, s)
        unc = tested & ~head & ~block
        Ttrue = np.where(gate, lfc, np.nan)
        key = np.sign(BC.signed_key(I, s))
        C4, _ = clean4_lolo(I, s)
        c4h = np.where(np.isfinite(C4) & (np.abs(C4 - 0.5) >= 0.05), np.where(C4 > 0.5, -1.0, 1.0), 0.0)
        # the KGID form for the record: v3 identity sign on CLEAN4-covered genes, only on each target's weak-key half
        if s == "h1":
            own, refs = prof["h1"], [prof[k] for k in ("ctxA", "ctxB", "ctxC")]
        else:
            own, refs = dict(zip(z["genes"].astype(str), z["ctrl_cpm"].astype(float))), [prof[k] for k in ("ctxA", "ctxB", "ctxC", "h1")]
        Iv = DC.identity(own, refs, list(G))
        zi = (Iv - np.nanmean(Iv)) / np.nanstd(Iv)
        zc = (C4 - 0.5) / np.nanstd(C4 - 0.5)
        sc = np.where(np.isfinite(zc), zc + 0.5 * np.nan_to_num(zi), np.nan)
        v3 = np.where(np.isfinite(sc) & (sc != 0), np.where(sc > 0, -1.0, 1.0), 0.0)
        ak = np.abs(BC.signed_key(I, s))
        weak = np.zeros(ak.shape, bool)
        for i in range(ak.shape[0]):
            nz = ak[i] > 0
            weak[i] = (~nz) | (ak[i] < (np.quantile(ak[i][nz], 0.5) if nz.any() else 0.0))
        kg = np.where(weak, np.broadcast_to(v3, weak.shape), 0.0)
        rec = {}
        for f in REG:
            Tp = DB.pool_T(Ttrue, p, f)
            a = pred_dnmae(c4h, key, add, Tp, unc)
            b = pred_dnmae(kg, key, add, Tp, unc)
            rec[str(f)] = {"c4h_minus_base": a[0], "kgid_minus_base": b[0], "kgid_minus_c4h": b[0] - a[0], "n_scored": a[1],
                           "c4h_flips": a[2], "c4h_new": a[3]}
        out[s] = rec
        lines.append(f"== {s}: " + " | ".join(f"f {f}: c4h−base {rec[str(f)]['c4h_minus_base']:+.5f}  KGID−c4h "
                                               f"{rec[str(f)]['kgid_minus_c4h']:+.5f} (n {rec[str(f)]['n_scored']})" for f in REG))
    retro = out["h1"]["1.0"]["c4h_minus_base"] < 0 and out["h1"]["0.35"]["c4h_minus_base"] < 0
    lines.append(f"\nRETRODICTION (h1 c4h − base < 0 at f 1 and f 0.35; real emitter −0.0015 / −0.0021): {'PASS' if retro else 'FAIL'}")
    if retro:
        m = float(np.mean([out[s]["0.35"]["c4h_minus_base"] for s in B.CANCER]))
        lines.append(f"REGISTERED PREDICTION for 496a − 446bj2 board Δnmae (raw): sign {'NEGATIVE (the hook improves nmae)' if m < 0 else 'POSITIVE (the hook worsens nmae)'}; "
                     f"cancer-fold mean at f 0.35 {m:+.5f}; per fold " + ", ".join(f"{s} {out[s]['0.35']['c4h_minus_base']:+.5f}" for s in B.CANCER))
        out["prediction"] = {"sign": int(np.sign(m)), "mean_f035": m}
    else:
        lines.append("Nothing registered: the instrument fails its h1 retrodiction.")
    out["retrodiction_pass"] = retro
    print("\n".join(lines))
    (HERE / "work/block613d.out").write_text("\n".join(lines) + "\n")
    (HERE / "work/block613d.json").write_text(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
