"""213 — the MECHANISM: a transport law is a per-gene, sign-preserving gain.

To first order in the effect, with `dp = p^pert - p^ctrl`,

    g_c   =  f_c(p^pert_don) - f_c(p^ctrl_don)  ~=  dp / (p_don + c)
    dh    =  the delivered log1p effect         ~=  dp_rec / (p_rec + c0),  dp_rec ~= g_c (p_rec + c)

so, composing the two halves of the law,

    dh_c / dh_c0  ~=  [ (p_don + c0) / (p_don + c) ] * [ (p_rec + c) / (p_rec + c0) ]        (*)

**a per-gene gain, determined entirely by where the donor's and the recipient's basal levels sit
relative to `c`, and strictly positive.**  Three consequences, all measured here:

  1. the law CANNOT change a sign, so it cannot touch the ordering channel (8.1) at all -- it
     is the same structural objection 184 raised against the Tahoe operator, arriving from the
     opposite direction;
  2. the whole action on the DIRECTION channel is a re-weighting of which genes dominate the
     cosine, and (*) says the weight is 1 wherever donor and recipient agree in expression --
     so the law can only ever pay for donor/recipient basal MISMATCH;
  3. the delivered AMPLITUDE moves with `c` even though the pre-render amplitude does not,
     because the second bracket of (*) is > 1 for every `c > c0` -- which is why an
     alternative law is an amplitude-changing arm under rule 13.

    uv run python mech213.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
sys.path.insert(0, str(HERE))

import plolo  # noqa: E402
plolo.WORK = WORK / "plolo"
import law213 as L  # noqa: E402
import prep213 as P  # noqa: E402

TS = L.TS
OUT: dict = {}


def donor_half(law_key="logc40", donor="orion_hct116", n_rows=400):
    """Measure `g_c / h` against the closed form `(p_don + c0)/(p_don + c)`."""
    law = L.all_laws()[law_key]
    blk = P.DonorBlock(donor, load_block=False)
    d = np.load(blk.bank.delta_path, mmap_mode="r")
    rows = np.random.default_rng(0).choice(d.shape[0], n_rows, replace=False)
    cols = blk.amap.unique_source_columns
    h = np.asarray(d[np.sort(rows)][:, cols], dtype=np.float64)
    del d
    p = blk.p_uniq
    g = law.to_law(h, p)
    live = np.abs(h) > 1e-6
    ratio = np.where(live, g / np.where(live, h, 1.0), np.nan)
    pred = (p + L.C0) / (p + law.param)
    med = np.nanmedian(ratio, axis=0)
    ok = np.isfinite(med) & (np.sum(live, axis=0) >= 20)
    r = float(np.corrcoef(np.log(med[ok]), np.log(pred[ok]))[0, 1])
    signflip = float(np.mean((np.sign(g[live]) != np.sign(h[live]))))
    OUT["donor_half"] = dict(
        law=law_key, donor=donor, n_rows=int(n_rows), n_genes=int(ok.sum()),
        corr_log_ratio_vs_closed_form=r,
        median_ratio_range=[float(np.nanmin(med[ok])), float(np.nanmax(med[ok]))],
        closed_form_range=[float(pred[ok].min()), float(pred[ok].max())],
        sign_flip_rate=signflip,
        median_abs_rel_error=float(np.nanmedian(np.abs(med[ok] / pred[ok] - 1.0))))
    print(f"[donor half] {law_key} on {donor}: corr(log measured, log closed form) = {r:.5f}; "
          f"median |rel err| {OUT['donor_half']['median_abs_rel_error']:.4f}; "
          f"SIGN FLIPS {100*signflip:.4f}% of {int(live.sum())} live entries")


def recipient_half(law_key="logc40"):
    """Render the SAME effect matrix under the incumbent and under `law_key`.

    This isolates the recipient half of (*): the donor side is held fixed by construction.
    """
    law = L.all_laws()[law_key]
    inc = L.all_laws()["logc20"]
    eff = np.load(WORK / "eff" / "chal_logc20.npz")["effect"]
    res = {}
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate = (cpm > plolo.MIN_CPM) & (b > 0)
        d0 = L.deliver_law(eff, b, gate, inc, None)
        d1 = L.deliver_law(eff, b, gate, law, None)
        p = expr / TS
        pred = (p + law.param) / (p + L.C0)
        live = gate[None, :] & (np.abs(d0) > 1e-8)
        ratio = np.where(live, d1 / np.where(live, d0, 1.0), np.nan)
        med = np.nanmedian(ratio, axis=0)
        ok = np.isfinite(med) & gate
        r = float(np.corrcoef(np.log(np.abs(med[ok])), np.log(pred[ok]))[0, 1])
        flip = float(np.mean(np.sign(d1[live]) != np.sign(d0[live])))
        res[ctx] = dict(corr=r, sign_flip_rate=flip,
                        median_gain=float(np.median(med[ok])),
                        gain_p05=float(np.percentile(med[ok], 5)),
                        gain_p95=float(np.percentile(med[ok], 95)),
                        norm_ratio=float(np.median(np.linalg.norm(d1, axis=1))
                                         / np.median(np.linalg.norm(d0, axis=1))),
                        cos_d0_d1=float(np.mean(np.sum(
                            d0 / np.linalg.norm(d0, axis=1, keepdims=True)
                            * d1 / np.linalg.norm(d1, axis=1, keepdims=True), axis=1))))
        print(f"[recipient half] ctx {ctx}: corr(log gain, log closed form) = {r:.5f}; "
              f"gain median {res[ctx]['median_gain']:.4f} "
              f"[{res[ctx]['gain_p05']:.3f}, {res[ctx]['gain_p95']:.3f}]; "
              f"SIGN FLIPS {100*flip:.4f}%; delivered norm x{res[ctx]['norm_ratio']:.4f}; "
              f"row cosine to the incumbent {res[ctx]['cos_d0_d1']:.5f}")
    OUT["recipient_half"] = dict(law=law_key, per_context=res)


def full_delivered(law_key="logc40"):
    """Both halves at once: the shipped effect vs the law's own, delivered, per context."""
    law = L.all_laws()[law_key]
    inc = L.all_laws()["logc20"]
    e0 = np.load(WORK / "eff" / "chal_logc20.npz")["effect"]
    e1 = np.load(WORK / "eff" / f"chal_{law_key}.npz")["effect"]
    res = {}
    for ctx in ("A", "B", "C"):
        b = plolo.challenge_bctrl(ctx)
        expr = np.expm1(b)
        cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
        gate = (cpm > plolo.MIN_CPM) & (b > 0)
        d0 = L.deliver_law(e0, b, gate, inc, None)
        d1 = L.deliver_law(e1, b, gate, law, None)
        live = gate[None, :] & (np.abs(d0) > 1e-8)
        flip = float(np.mean(np.sign(d1[live]) != np.sign(d0[live])))
        u0 = d0 / np.linalg.norm(d0, axis=1, keepdims=True)
        u1 = d1 / np.linalg.norm(d1, axis=1, keepdims=True)
        res[ctx] = dict(sign_flip_rate=flip,
                        row_cosine=float(np.mean(np.sum(u0 * u1, axis=1))),
                        norm_ratio=float(np.median(np.linalg.norm(d1, axis=1))
                                         / np.median(np.linalg.norm(d0, axis=1))))
        print(f"[full] ctx {ctx}: row cosine to the shipped delivered effect "
              f"{res[ctx]['row_cosine']:.5f}; SIGN FLIPS {100*flip:.4f}%; "
              f"delivered norm x{res[ctx]['norm_ratio']:.4f}")
    OUT["full_delivered"] = dict(law=law_key, per_context=res)


def basal_mismatch():
    """(*) is 1 wherever donor and recipient agree, so the law can only pay for MISMATCH."""
    b = plolo.challenge_bctrl("A")
    e = np.expm1(b)
    p_rec = e / TS
    rows = []
    for donor in ("k562_gwps", "orion_hct116", "orion_hek293t"):
        blk = P.DonorBlock(donor, load_block=False)
        tp = blk.amap.unique_target_positions
        p_don = blk.p_uniq
        cpmr = 1e6 * p_rec[tp]
        cpmd = 1e6 * p_don
        m = (cpmr > 10) & (cpmd > 0)
        lr = np.log10(cpmr[m] / cpmd[m])
        rows.append(dict(donor=donor, n=int(m.sum()),
                         median_log10_ratio=float(np.median(lr)),
                         iqr=[float(np.percentile(lr, 25)), float(np.percentile(lr, 75))],
                         frac_within_2x=float(np.mean(np.abs(lr) < np.log10(2.0))),
                         corr_log_cpm=float(np.corrcoef(np.log(cpmr[m]),
                                                        np.log(cpmd[m]))[0, 1])))
        print(f"[mismatch] {donor} vs context A: corr(log CPM) {rows[-1]['corr_log_cpm']:.4f}, "
              f"{100*rows[-1]['frac_within_2x']:.1f}% of genes within 2x, "
              f"IQR of log10 ratio {rows[-1]['iqr']}")
    OUT["basal_mismatch"] = rows


if __name__ == "__main__":
    donor_half()
    recipient_half()
    full_delivered()
    basal_mismatch()
    (WORK / "mech213.json").write_text(json.dumps(OUT, indent=1, default=float))
    print(f"\nwrote {WORK}/mech213.json")
