"""227 — Phase 2: read Aginglab.com's six-tuple off the board and say what produces it.

Three questions, all pure board arithmetic on `data/leaderboard/2026-09-04T0433Z.json`:

  1. Reproduce the mandate's board facts as known-answer rows: our and their per-member
     scaled values, their implied transfer cosine below ours, and
     corr(implied cosine, scaled `nmae`) = -0.394 over the 55 teams above overall 0.15.
  2. Rule amplitude in or out AS AN EXPLANATION OF THEIR ROW, using 187's own board law
     `raw_nmae(alpha) = 0.944330 - 0.0416 ln(alpha/0.982)` and its `mse` counterpart.  If the
     alpha that buys their `nmae` destroys their `mse`, amplitude is not what they did.
  3. Characterise the eight teams above them on `nmae` -- the ones with catastrophic cosine --
     and ask whether Aginglab sits in that cluster or apart from it, on the four members that
     read a call set (`fid`, `jac`, `reach`) rather than a profile.
"""
from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
BOARD = os.path.join(REPO, "data/leaderboard/2026-09-04T0433Z.json")

RAW = dict(pds="pds_cosine", mse="expr_mse_unbiased_capped_norm", nmae="de_wilcoxon_lfc_nmae",
           fid="de_wilcoxon_direction_fidelity_yield_raw",
           reach="de_wilcoxon_direction_reach_raw", jac="de_wilcoxon_sig_jaccard")
SC = {k: f"score_{k}" for k in RAW}
US, THEM = "Ahmet Dedeler", "Aginglab.com"

# 197's inversion of the board's `mse` into a transfer cosine (AGENTS.md 4.1c).
COS_A, COS_B = 1.023611, 2.074495
# 187 sec 5: the board's own measured amplitude response of each member, per unit ln(alpha).
NMAE_SLOPE = -0.0416
NMAE_A23 = 0.944330


def cos_of(e):
    return (COS_A - e[RAW["mse"]]) / COS_B


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    def ok(e):
        for v in list(RAW.values()) + list(SC.values()) + ["score_avg"]:
            x = e.get(v)
            if x is None or not np.isfinite(float(x)):
                return False
        return True

    ALL = json.load(open(BOARD))["entries"]
    E = [e for e in ALL if ok(e)]
    best = {}
    for e in E:
        t = e["team_name"]
        if t not in best or e["score_avg"] > best[t]["score_avg"]:
            best[t] = e
    teams = sorted(best.values(), key=lambda e: -e["score_avg"])
    out = {"n_entries": len(E), "n_teams": len(teams)}

    us, them = best[US], best[THEM]
    out["KA_members"] = {
        m: dict(us_raw=us[RAW[m]], us_sc=us[SC[m]], them_raw=them[RAW[m]],
                them_sc=them[SC[m]], their_edge_sc=them[SC[m]] - us[SC[m]])
        for m in RAW}
    out["KA_members"]["overall"] = dict(us_sc=us["score_avg"], them_sc=them["score_avg"],
                                        their_edge_sc=them["score_avg"] - us["score_avg"])
    out["KA_cosine"] = dict(us=cos_of(us), them=cos_of(them),
                            their_edge=cos_of(them) - cos_of(us))

    # ---- corr(implied cosine, scaled nmae) over the teams above overall 0.15 --------------
    top = [e for e in teams if e["score_avg"] > 0.15]
    c = np.array([cos_of(e) for e in top])
    nm = np.array([e[SC["nmae"]] for e in top])
    out["KA_corr"] = dict(n=len(top), corr_cos_vs_scaled_nmae=float(np.corrcoef(c, nm)[0, 1]))
    for m in RAW:
        v = np.array([e[SC[m]] for e in top])
        out["KA_corr"][f"corr_cos_vs_{m}"] = float(np.corrcoef(c, v)[0, 1])
        out["KA_corr"][f"corr_nmae_vs_{m}"] = float(np.corrcoef(nm, v)[0, 1])

    # ---- can amplitude explain their row? -------------------------------------------------
    d_raw = them[RAW["nmae"]] - us[RAW["nmae"]]
    ln_alpha = d_raw / NMAE_SLOPE
    out["amplitude_test"] = dict(
        d_raw_nmae=d_raw, ln_alpha_required=float(ln_alpha),
        alpha_multiple=float(np.exp(ln_alpha)),
        note=("187 sec 6 measured the h1-test nmae(s) curve turning over at s = 1.29 with a total "
              "gain of 0.0027 raw; the required multiple is far past the turnover, where the "
              "curve is RISING.  187 sec 3's global-alpha table charges mse -0.1216 scaled at "
              "alpha = 1.10 alone."),
        their_mse_edge_sc=them[SC["mse"]] - us[SC["mse"]],
    )

    # ---- the cluster above them on nmae ---------------------------------------------------
    by_nmae = sorted(teams, key=lambda e: -e[SC["nmae"]])[:12]
    out["nmae_leaders"] = [
        dict(team=e["team_name"], overall=e["score_avg"], cos=cos_of(e),
             **{f"sc_{m}": e[SC[m]] for m in RAW},
             **{f"raw_{m}": e[RAW[m]] for m in RAW})
        for e in by_nmae]

    # ---- where Aginglab sits relative to that cluster -------------------------------------
    cluster = [e for e in by_nmae if cos_of(e) < 0.0]
    out["cluster_summary"] = dict(
        n_negative_cosine=len(cluster),
        med={m: float(np.median([e[SC[m]] for e in cluster])) for m in RAW},
        them={m: them[SC[m]] for m in RAW},
        us={m: us[SC[m]] for m in RAW},
    )

    # ---- the (q, k) frontier implied by a raw nmae ----------------------------------------
    def frontier(raw):
        credit = 1.0 - raw
        return dict(raw=raw, credit=credit, q_at_k1=1.0 - raw / 2.0,
                    k_at_q_0p537=credit / (2 * 0.537 - 1.0))
    out["identity_frontier"] = {n: frontier(best[n][RAW["nmae"]])
                                for n in (US, THEM, "cqawesome") if n in best}

    json.dump(out, open(os.path.join(WORK, "board227.json"), "w"), indent=1)

    # ------------------------------------------------------------------ report
    print(f"board 2026-09-04T0433Z: {out['n_entries']} scored entries, "
          f"{out['n_teams']} teams")
    print(f"\n== KA: our row vs Aginglab.com ==")
    print(f"  {'member':8s} {'us raw':>10s} {'us sc':>10s} {'them raw':>10s} {'them sc':>10s} "
          f"{'their edge':>11s}")
    for m in ("pds", "mse", "nmae", "fid", "reach", "jac"):
        d = out["KA_members"][m]
        print(f"  {m:8s} {d['us_raw']:10.5f} {d['us_sc']:+10.5f} {d['them_raw']:10.5f} "
              f"{d['them_sc']:+10.5f} {d['their_edge_sc']:+11.5f}")
    d = out["KA_members"]["overall"]
    print(f"  {'overall':8s} {'':10s} {d['us_sc']:+10.5f} {'':10s} {d['them_sc']:+10.5f} "
          f"{d['their_edge_sc']:+11.5f}")
    kc = out["KA_cosine"]
    print(f"  implied transfer cosine: us {kc['us']:.5f}  them {kc['them']:.5f}  "
          f"edge {kc['their_edge']:+.5f}")
    k = out["KA_corr"]
    print(f"\n== KA: over the {k['n']} teams above overall 0.15 ==")
    print(f"  corr(implied cosine, scaled nmae) = {k['corr_cos_vs_scaled_nmae']:+.4f}")
    print("  corr(implied cosine, .) : " + "  ".join(
        f"{m} {k[f'corr_cos_vs_{m}']:+.3f}" for m in RAW))
    print("  corr(scaled nmae,   .) : " + "  ".join(
        f"{m} {k[f'corr_nmae_vs_{m}']:+.3f}" for m in RAW))

    a = out["amplitude_test"]
    print(f"\n== can a global amplitude change explain their nmae? ==")
    print(f"  d raw nmae {a['d_raw_nmae']:+.5f} at 187's board law -0.0416/ln(alpha) needs "
          f"alpha x {a['alpha_multiple']:.2f}")
    print(f"  their mse edge is only {a['their_mse_edge_sc']:+.5f} scaled; 187 charges "
          f"-0.1216 for alpha = 1.10 alone.  RULED OUT.")

    print(f"\n== the top 12 by scaled nmae ==")
    print(f"  {'team':22s} {'overall':>8s} {'cos':>8s} {'nmae':>8s} {'pds':>8s} {'mse':>8s} "
          f"{'fid':>8s} {'reach':>8s} {'jac':>8s}")
    for e in out["nmae_leaders"]:
        print(f"  {e['team'][:22]:22s} {e['overall']:8.5f} {e['cos']:+8.4f} "
              f"{e['sc_nmae']:+8.5f} {e['sc_pds']:+8.5f} {e['sc_mse']:+8.5f} "
              f"{e['sc_fid']:+8.5f} {e['sc_reach']:+8.5f} {e['sc_jac']:+8.5f}")
    cs = out["cluster_summary"]
    print(f"\n  the {cs['n_negative_cosine']} of those 12 with NEGATIVE implied cosine, median:")
    print("    " + "  ".join(f"{m} {cs['med'][m]:+.4f}" for m in RAW))
    print("    Aginglab: " + "  ".join(f"{m} {cs['them'][m]:+.4f}" for m in RAW))
    print("    us:       " + "  ".join(f"{m} {cs['us'][m]:+.4f}" for m in RAW))

    print(f"\n== the uniform-k identity frontier (187 sec 5: right as identity, wrong as "
          f"measurement) ==")
    for n, f in out["identity_frontier"].items():
        print(f"  {n:16s} raw {f['raw']:.5f}  credit {f['credit']:.5f}  "
              f"min q at k=1 {f['q_at_k1']:.5f}  k needed at q=0.537 {f['k_at_q_0p537']:.3f}")


if __name__ == "__main__":
    main()
