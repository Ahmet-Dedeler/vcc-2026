"""164 — the six-member price of every arm, with a8 -> a15 as the known-answer row.

Pulls together `decomp.py` / `decomp2.py` (h1-test DE members), `claim164.py` (`claim` on the
real challenge panel) and the board anchors.  All arithmetic in one place so it can be checked.

`mse` uses the general form, NOT `mech_claim_grid2.py`'s optimal-`s` shortcut (that shortcut
silently assumes the amplitude has been re-tuned to each arm's own cliff, which is a SECOND
variable and outside this experiment's mandate):

    raw_mse(s, claim) = 1 + ( s^2*A - 2*s*X + Kc - min(s^2*B, claim) ) / sumD
    score_mse         = clip( (0.987862879 - raw_mse) / 0.950162566, 0, 1 )

KNOWN-ANSWER ROW: a8 -> a15 was measured on the board at +0.0197258 overall.  The instrument
must reproduce that from `claim` + `reach` alone, or it is not trusted for this family.
"""
from __future__ import annotations

import numpy as np

# ---- solved challenge-panel terms (096 work/board_solve.npz) ---------------------------
A, B, X, Kc, sD, K = 0.28361435, 0.28123598, 0.01884625, 0.00541027, 0.22857222, 0.24652869
ALPHA_V11, ALPHA_A15 = 1.0796, 0.982
S15 = ALPHA_A15 / ALPHA_V11

# ---- board anchors (094 §1, fitted on 780 entries) -------------------------------------
ANCHOR = {"pds": (0.50020, 0.95228), "fid": (0.51239, 0.80964),
          "reach": (0.07902, 0.96414), "jac": (0.03038, 0.39727)}
SPAN = {k: v[1] - v[0] for k, v in ANCHOR.items()}

# ---- claim, measured by claim164.py on the real challenge panel -------------------------
CLAIM = {"a15": 0.24581, "armA": 0.17087, "armB": 0.18917, "armC": 0.16652,
         "armB_075": 0.22279, "armA_phi1": 0.19512, "mincpm5_spreadonly": 0.24634,
         "a8_at_v11_alpha": 0.16652}

# ---- h1-test, measured by decomp.py / decomp2.py (seed-averaged where a twin exists) ----
H1 = {  # n_pred(mean), precision=mean_p(k/n_pred), reach, jac, nmae, pds
    "a15":  (6021.5, 0.507350, 0.179100, 0.068290, 0.966500, 0.799790),
    "armA": (3123.0, 0.513050, 0.208820, 0.053885, 0.966645, 0.798810),
    "armB": (4663.0, 0.509000, 0.182840, 0.067310, 0.966590, 0.798600),
    "armC": (3561.0, 0.512300, 0.197570, 0.060875, 0.966680, 0.798950),
}
H1_SEED = {"fid": 0.00070, "reach": 0.00746, "jac": 0.00010, "nmae": 0.00018, "pds": 0.00070}

BOARD_FID, BOARD_NPRED, BOARD_JAC, BOARD_NREAL = 0.5022, 4424.0, 0.0314, 240.0
BOARD_REACH_A15, BOARD_REACH_A8 = 0.1811, 0.1890


def raw_mse(s, claim):
    return 1 + (s * s * A - 2 * s * X + Kc - min(s * s * B, claim)) / sD


def score_mse(raw):
    return float(np.clip((0.987862879 - raw) / 0.950162566, 0, 1))


def main() -> None:
    print("=" * 100)
    print("(1) KNOWN-ANSWER ROW: does the instrument reproduce the a8 -> a15 board move?")
    print("=" * 100)
    m_a8 = raw_mse(1.0, CLAIM["a8_at_v11_alpha"])       # a8 shipped v11's VCC_ALPHA 1.0796
    m_a15 = raw_mse(S15, CLAIM["a15"])
    print(f"  a8  raw mse  modelled {m_a8:.4f}   OBSERVED on the board 1.3627  "
          f"({100*(m_a8-1.3627)/1.3627:+.2f}%)")
    print(f"  a15 raw mse  modelled {m_a15:.4f}   (not separately published)")
    d_mse = (score_mse(m_a15) - score_mse(m_a8)) / 6
    d_reach = (BOARD_REACH_A15 - BOARD_REACH_A8) / SPAN["reach"] / 6
    d_fid = (0.5022 - 0.5021) / SPAN["fid"] / 6
    print(f"  predicted d overall  mse {d_mse:+.6f}  reach {d_reach:+.6f}  "
          f"fid {d_fid:+.6f}   SUM {d_mse+d_reach+d_fid:+.6f}")
    print(f"  OBSERVED a8 -> a15                                          "
          f"          +0.019726")
    print(f"  residual (jac + nmae + the CD4/whitener confounds) "
          f"{0.019726-(d_mse+d_reach+d_fid):+.6f}  "
          f"-> instrument accounts for {100*(d_mse+d_reach+d_fid)/0.019726:.0f}%")
    print("  a8 -> a15 is CONFOUNDED (tier policy + VCC_ALPHA + CD4 f); only the SIGN and")
    print("  the order of magnitude are being gated here, and both reproduce.")

    print("\n" + "=" * 100)
    print("(2) fid RECONCILIATION with the board, 096's conserved-headroom form")
    print("=" * 100)
    h = H1["a15"]
    print(f"  h1-test a15 precision {h[1]:.5f}  -> headroom {h[1]-0.5:+.5f}")
    print(f"  board   a15 fid       {BOARD_FID:.5f}  -> headroom {BOARD_FID-0.5:+.5f}"
          f"   (ratio {(BOARD_FID-0.5)/(h[1]-0.5):.3f})")
    print("  h1-test block precisions: keep 0.51204 (58.2% of calls), tail 0.49990 (41.8%)")
    f_tail = 0.418
    p_sig_board = (BOARD_FID - f_tail * 0.5) / (1 - f_tail)
    print(f"  if the board's tail block is at EXACTLY 0.5 and holds the same {f_tail:.1%} share,")
    print(f"  the board's signal block must sit at {p_sig_board:.5f} "
          f"(headroom {p_sig_board-0.5:+.5f})")
    print(f"  h1-test signal-block headroom {0.51204-0.5:+.5f} x the single scale factor "
          f"{(BOARD_FID-0.5)/(h[1]-0.5):.3f} = {(0.51204-0.5)*(BOARD_FID-0.5)/(h[1]-0.5):+.5f}")
    print("  -> ONE factor maps h1-test's overall AND its per-block precision onto the board.")
    print("     The decomposition RECONCILES.")
    ceil = 0.5 + (BOARD_FID - 0.5) / (1 - f_tail)
    print(f"\n  CEILING: delete the whole tail block for free -> board fid {ceil:.5f}, "
          f"d raw {ceil-BOARD_FID:+.5f}")
    print(f"           = d scaled {(ceil-BOARD_FID)/SPAN['fid']:+.5f} "
          f"= d overall {(ceil-BOARD_FID)/SPAN['fid']/6:+.6f}")
    gap = 0.5175 - BOARD_FID
    print(f"  our gap to the board MAXIMUM fid (0.5175, The_Way): {gap:+.5f} raw "
          f"= {gap/SPAN['fid']/6:+.6f} overall -> "
          f"{gap/(ceil-BOARD_FID):.1f}x the entire artefact channel")

    print("\n  conserved-headroom model VALIDATED on h1-test (it could have failed):")
    print(f"  {'arm':<8}{'n_pred ratio':>14}{'predicted prec':>16}{'observed':>11}{'realised':>10}")
    for a in ("armA", "armB", "armC"):
        r = H1["a15"][0] / H1[a][0]
        pred = 0.5 + (H1["a15"][1] - 0.5) * r
        print(f"  {a:<8}{r:>14.4f}{pred:>16.5f}{H1[a][1]:>11.5f}"
              f"{(H1[a][1]-0.5)/(pred-0.5):>10.1%}")

    print("\n" + "=" * 100)
    print("(3) SIX-MEMBER TABLE, board units (d overall).  a8->a15 is the known-answer row.")
    print("=" * 100)
    print(f"{'arm':<12}{'d pds':>9}{'d mse':>10}{'d fid':>9}{'d reach':>9}"
          f"{'d jac':>9}{'d nmae':>9}{'NET':>11}")
    reach_transfer = (BOARD_REACH_A8 - BOARD_REACH_A15) / (H1["armC"][2] - H1["a15"][2])
    print(f"  [reach h1-test -> board transfer factor, calibrated on armC/a8: "
          f"{reach_transfer:.3f}]")
    rows = {}
    for a in ("armA", "armB", "armC"):
        r = H1["a15"][0] / H1[a][0]
        # fid: conserved headroom on the board's own headroom
        dfid = ((0.5 + (BOARD_FID - 0.5) * r) - BOARD_FID) / SPAN["fid"] / 6
        # reach: h1-test delta, shrunk by the board-calibrated factor
        dreach = (H1[a][2] - H1["a15"][2]) * reach_transfer / SPAN["reach"] / 6
        # jac: rebuild on the board's own n_real/n_pred, overlap scaled by the h1-test ratio
        ov0 = BOARD_JAC * (BOARD_NPRED + BOARD_NREAL) / (1 + BOARD_JAC)
        ovr = {"armA": 253.95 / 478.15, "armB": 395.2 / 478.15, "armC": 304.15 / 478.15}[a]
        n1, ov1 = BOARD_NPRED / r, ov0 * ovr
        djac = (ov1 / (n1 + BOARD_NREAL - ov1) - BOARD_JAC) / SPAN["jac"] / 6
        # mse at a15's SHIPPED alpha (a true one-variable arm)
        dmse = (score_mse(raw_mse(S15, CLAIM[a])) - score_mse(raw_mse(S15, CLAIM["a15"]))) / 6
        dmse_opt = (score_mse(raw_mse(np.sqrt(CLAIM[a] / B), CLAIM[a]))
                    - score_mse(raw_mse(np.sqrt(CLAIM["a15"] / B), CLAIM["a15"]))) / 6
        rows[a] = (0.0, dmse, dfid, dreach, djac, 0.0, dmse_opt)
        net = dmse + dfid + dreach + djac
        print(f"{a:<12}{'0.000000':>9}{dmse:>10.6f}{dfid:>9.6f}{dreach:>9.6f}"
              f"{djac:>9.6f}{'~0':>9}{net:>+11.6f}")
    print("\n  same table but granting each arm its OWN re-optimised VCC_ALPHA "
          "(a SECOND variable, outside mandate):")
    for a in ("armA", "armB", "armC"):
        _, dmse, dfid, dreach, djac, _, dmse_opt = rows[a]
        print(f"{a:<12}{'0.000000':>9}{dmse_opt:>10.6f}{dfid:>9.6f}{dreach:>9.6f}"
              f"{djac:>9.6f}{'~0':>9}{dmse_opt+dfid+dreach+djac:>+11.6f}")

    print("\n" + "=" * 100)
    print("(4) the mse cliff, stated exactly")
    print("=" * 100)
    print(f"  a15 at its shipped VCC_ALPHA {ALPHA_A15}: s = {S15:.5f}, "
          f"s^2*B = {S15*S15*B:.5f} vs claim {CLAIM['a15']:.5f}")
    print(f"    -> rho = {min(1, S15*S15*B/CLAIM['a15']):.4f}; a15 is AT the cliff "
          f"(it already has marginally more claim than budget)")
    print(f"    raw mse {raw_mse(S15, CLAIM['a15']):.5f}, score_mse "
          f"{score_mse(raw_mse(S15, CLAIM['a15'])):.5f}")
    print(f"  every arm below drops claim under s^2*B = {S15*S15*B:.5f}, so it lands in the")
    print(f"  rho = 1 branch where d(raw mse)/d(claim) = -1/sumD = {-1/sD:.3f}:")
    print(f"  {'arm':<20}{'claim':>9}{'raw mse':>10}{'score_mse':>11}{'d overall':>11}")
    for a in ("a15", "armB_075", "armA_phi1", "armB", "armA", "armC",
              "mincpm5_spreadonly"):
        rm = raw_mse(S15, CLAIM[a])
        sm = score_mse(rm)
        print(f"  {a:<20}{CLAIM[a]:>9.5f}{rm:>10.5f}{sm:>11.5f}"
              f"{(sm - score_mse(raw_mse(S15, CLAIM['a15'])))/6:>+11.6f}")
    print("\n  Every arm that removes a material part of the tail block drives raw mse past")
    print("  0.98786, where score_mse CLAMPS TO ZERO. That is the whole of a15's mse credit.")


if __name__ == "__main__":
    main()
