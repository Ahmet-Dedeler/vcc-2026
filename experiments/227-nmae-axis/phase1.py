"""227 Phase 1 — where the `nmae` L1 mass actually sits, and what each channel can pay.

THE STRUCTURAL FACT THIS IS BUILT ON (spec sec 6, and production's own source).  The scorer forms
both sides' lfc against the SAME measured control, so in the numerator that control cancels:

    lfc_hat - lfc  =  log2(m_hat_pert + eps) - log2(m_real_pert + eps)

`nmae` is therefore an L1 error on our predicted arithmetic-mean-CPM profile of the PERTURBED
state, over the reference's own significant gate, divided by the reference's |lfc|.  It reads
neither our call set (the gate is reference-only) nor our emitter (per-gene totals are the design's
count-for-count and every emitted cell carries the same library total, so the arithmetic mean of
per-cell CPM IS the design pseudobulk).  So the only variable is the design profile.

This script asks, on h1-test with the reference's real BH gate:
  1. how the L1 mass is distributed -- by control CPM, by |t|, by sign-correctness, by row;
  2. what each channel's CEILING is, measured with oracles, so no arm is chased below its cap;
  3. what an HONEST (leave-one-perturbation-out) version of each channel actually pays, and
     what it costs in `pds` and aggregate cosine ON THE SAME RUN.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
sys.path.insert(0, os.path.join(REPO, "experiments/153-de-bench"))
from debench import Surface  # noqa: E402

E169 = os.path.join(REPO, "experiments/169-sharpening-six-member")
TABLE = os.path.join(E169, "de169.q1.00_mf.npz")     # a23-family production geometry, h1-test
EFFECT = os.path.join(E169, "work/h1test-effect.npz")


# ------------------------------------------------------------------ gate extraction
def gate_of(S):
    """The scored gate, as a list of (row, gene-index array).  Reference-only, by definition."""
    out = []
    for p in range(len(S.perts)):
        s = np.flatnonzero(S.real_sig[p] & np.isfinite(S.real_lfc[p]))
        if s.size < 10:
            continue
        if not np.isfinite(np.abs(S.real_lfc[p, s]).sum()) or np.abs(S.real_lfc[p, s]).sum() == 0:
            continue
        out.append((p, s))
    return out


def nmae_of(G, P, T):
    """mean over rows of Sum|P-T| / Sum|T| on the gate `G`."""
    return float(np.mean([np.abs(P[p, s] - T[p, s]).sum() / np.abs(T[p, s]).sum()
                          for p, s in G]))


def rows_of(G, P, T):
    return np.array([np.abs(P[p, s] - T[p, s]).sum() / np.abs(T[p, s]).sum() for p, s in G])


def cos_agg(P, T, cols):
    """197's energy-weighted aggregate cosine: <P,T> / (||P|| ||T||) over the stacked matrix."""
    a, b = P[:, cols].ravel(), T[:, cols].ravel()
    return float(a @ b / max(np.linalg.norm(a) * np.linalg.norm(b), 1e-300))


def main() -> None:
    os.makedirs(WORK, exist_ok=True)
    out = {}
    S = Surface(TABLE)
    T = S.real_lfc.astype(np.float64)
    P0 = np.nan_to_num(S.pred_lfc.astype(np.float64))
    G = gate_of(S)
    E = np.load(EFFECT, allow_pickle=True)
    egenes = np.asarray([str(g) for g in E["genes"]])
    ecpm = E["cpm"].astype(np.float64)
    pos = {g: i for i, g in enumerate(egenes)}
    cpm = np.array([ecpm[pos[g]] if g in pos else np.nan for g in S.genes])
    assert np.isfinite(cpm).all()

    base = nmae_of(G, P0, T)
    scored_cols = np.flatnonzero(np.isfinite(T).all(0) & np.isfinite(P0).all(0))
    out["base"] = dict(
        nmae=base, n_rows=len(G), n_genes=int(S.genes.size),
        gate_median=float(np.median([s.size for _, s in G])),
        gate_mean=float(np.mean([s.size for _, s in G])),
        pds=S.pds(), cos_agg=cos_agg(P0, T, scored_cols),
    )

    # ------------------------------------------------------------- 1. row distribution
    r = rows_of(G, P0, T)
    out["rows"] = dict(
        mean=float(r.mean()), median=float(np.median(r)),
        q=[float(np.quantile(r, x)) for x in (0.05, 0.25, 0.5, 0.75, 0.95)],
        max=float(r.max()), frac_above_1=float((r > 1.0).mean()),
        mass_above_1=float(np.clip(r - 1.0, 0, None).sum() / len(r)),
        nmae_if_rows_above_1_zeroed=float(np.minimum(r, 1.0).mean()),
    )

    # ------------------------------------------------------------- 2. strata
    def strata(labels, name, nbin):
        rec = []
        for b in range(nbin):
            num = den = kmass = 0.0
            ok_w = 0.0
            nz = ng = 0
            for p, s in G:
                m = labels[p, s] == b if labels.ndim == 2 else labels[s] == b
                if not m.any():
                    continue
                t, h = T[p, s][m], P0[p, s][m]
                dp = np.abs(T[p, s]).sum()               # the ROW's full denominator
                num += np.abs(h - t).sum() / dp
                den += np.abs(t).sum() / dp
                kmass += np.abs(h).sum() / dp
                w = np.abs(t)
                ok_w += w[np.sign(h) == np.sign(t)].sum() / dp
                nz += int((np.abs(h) < 1e-9).sum())
                ng += int(m.sum())
            n = len(G)
            rec.append(dict(bin=b, n_genes=ng, den_share=den / n, num_share=num / n,
                            K=kmass / n, q_w_share=ok_w / n,
                            q_w=(ok_w / den) if den > 0 else float("nan"),
                            frac_h_zero=nz / max(ng, 1)))
        return {name: rec}

    cpm_edges = np.array([5.0, 10.0, 20.0, 50.0, 150.0, 1e18])
    cpm_bin = np.searchsorted(cpm_edges, cpm, side="left")
    out.update(strata(cpm_bin, "by_control_cpm", len(cpm_edges)))
    out["cpm_edges"] = ["<=5", "5-10", "10-20", "20-50", "50-150", ">150"]

    tb = np.zeros(T.shape, np.int64)
    for p, s in G:
        a = np.abs(T[p, s])
        tb[p, s] = np.searchsorted(np.quantile(a, [0.2, 0.4, 0.6, 0.8]), a, side="left")
    out.update(strata(tb, "by_abs_t_quintile", 5))

    # ------------------------------------------------------------- 3. channel ceilings
    sgnT, sgnP = np.sign(T), np.sign(P0)
    ceil = {}
    ceil["as_shipped"] = base
    ceil["predict_zero"] = nmae_of(G, np.zeros_like(P0), T)
    ceil["oracle_everything"] = nmae_of(G, T, T)
    ceil["oracle_sign_our_magnitude"] = nmae_of(G, np.abs(P0) * sgnT, T)
    ceil["our_sign_oracle_magnitude"] = nmae_of(G, np.abs(T) * np.where(sgnP == 0, 1, sgnP), T)

    ss = np.linspace(0.2, 3.0, 141)
    vv = [nmae_of(G, P0 * float(x), T) for x in ss]
    ceil["best_global_scalar"] = float(np.min(vv))
    out["best_global_s"] = float(ss[int(np.argmin(vv))])

    # per-row oracle scalar (a ceiling for ANY per-target amplitude key, 158's object)
    P = P0.copy()
    for p, s in G:
        cand = [np.abs(P0[p, s] * x - T[p, s]).sum() for x in ss]
        P[p, s] = P0[p, s] * ss[int(np.argmin(cand))]
    ceil["oracle_per_row_scalar"] = nmae_of(G, P, T)

    # per-row oracle ADDITIVE offset on the gate (the compositional/renormalisation channel)
    P = P0.copy()
    offs = []
    for p, s in G:
        d = np.median(T[p, s] - P0[p, s])          # L1 optimum for a constant shift
        offs.append(float(d))
        P[p, s] = P0[p, s] + d
    ceil["oracle_per_row_offset"] = nmae_of(G, P, T)
    out["row_offsets"] = dict(mean=float(np.mean(offs)), median=float(np.median(offs)),
                              sd=float(np.std(offs)),
                              q=[float(np.quantile(offs, x)) for x in (0.05, 0.5, 0.95)])

    # per-CPM-stratum oracle scalar (one scalar per stratum, shared across all rows)
    P = P0.copy()
    st = {}
    for b in range(len(cpm_edges)):
        m = cpm_bin == b
        if not m.any():
            continue
        cand = []
        for x in ss:
            v = 0.0
            for p, s in G:
                mm = m[s]
                if mm.any():
                    v += np.abs(P0[p, s][mm] * x - T[p, s][mm]).sum() / np.abs(T[p, s]).sum()
            cand.append(v)
        x = float(ss[int(np.argmin(cand))])
        st[out["cpm_edges"][b]] = x
        P[:, m] = P0[:, m] * x
    ceil["oracle_cpm_stratum_scalar"] = nmae_of(G, P, T)
    out["cpm_stratum_scalars"] = st
    out["cpm_stratum_arm"] = dict(nmae=ceil["oracle_cpm_stratum_scalar"],
                                  pds=S.pds(pred_lfc=P), cos_agg=cos_agg(P, T, scored_cols))

    # per-gene target-independent additive correction, leave-one-perturbation-out
    acc = np.zeros(S.genes.size)
    cnt = np.zeros(S.genes.size)
    for p, s in G:
        acc[s] += T[p, s] - P0[p, s]
        cnt[s] += 1
    P = P0.copy()
    for p, s in G:
        c = np.maximum(cnt[s] - 1, 0)
        d = np.where(c > 0, (acc[s] - (T[p, s] - P0[p, s])) / np.maximum(c, 1), 0.0)
        P[p, s] = P0[p, s] + d
    ceil["lopo_per_gene_offset"] = nmae_of(G, P, T)
    out["lopo_gene_arm"] = dict(nmae=ceil["lopo_per_gene_offset"], pds=S.pds(pred_lfc=P),
                                cos_agg=cos_agg(P, T, scored_cols))
    out["ceilings"] = ceil

    # ------------------------------------------------------------- 4. zero / saturation audit
    nz = sum(int((np.abs(P0[p, s]) < 1e-9).sum()) for p, s in G)
    ng = sum(int(s.size) for p, s in G)
    huge = sum(int((np.abs(P0[p, s]) > 5.0).sum()) for p, s in G)
    l1_zero = float(np.mean([np.abs(T[p, s][np.abs(P0[p, s]) < 1e-9]).sum()
                             / np.abs(T[p, s]).sum() for p, s in G]))
    out["zero_audit"] = dict(n_gate_pairs=ng, frac_h_exactly_zero=nz / ng,
                             l1_share_of_zero_h_genes=l1_zero,
                             frac_h_above_5_log2=huge / ng)

    json.dump(out, open(os.path.join(WORK, "phase1.json"), "w"), indent=1)

    # ------------------------------------------------------------- report
    b = out["base"]
    print(f"h1-test / de169.q1.00_mf (a23-family production geometry)")
    print(f"  rows scored {b['n_rows']}   gate median {b['gate_median']:.0f} "
          f"mean {b['gate_mean']:.0f} genes")
    print(f"  nmae {b['nmae']:.6f}   pds {b['pds']:.6f}   cos_agg {b['cos_agg']:.6f}")
    rr = out["rows"]
    print(f"\n== per-row nmae distribution ==")
    print(f"  q05 {rr['q'][0]:.4f}  q25 {rr['q'][1]:.4f}  med {rr['q'][2]:.4f}  "
          f"q75 {rr['q'][3]:.4f}  q95 {rr['q'][4]:.4f}  max {rr['max']:.4f}")
    print(f"  rows above 1.0: {100*rr['frac_above_1']:.1f}%   they carry "
          f"{rr['mass_above_1']:.5f} of the mean; clipping every row at 1.0 -> "
          f"{rr['nmae_if_rows_above_1_zeroed']:.6f}")

    for key, edges in (("by_control_cpm", out["cpm_edges"]),
                       ("by_abs_t_quintile", [f"Q{i+1}" for i in range(5)])):
        print(f"\n== L1 mass {key} ==")
        print(f"  {'bin':8s} {'genes':>7s} {'den share':>10s} {'num share':>10s} "
              f"{'K':>8s} {'q_w':>7s} {'h==0':>7s}")
        for rec in out[key]:
            print(f"  {edges[rec['bin']]:8s} {rec['n_genes']:7d} {rec['den_share']:10.4f} "
                  f"{rec['num_share']:10.4f} {rec['K']:8.4f} {rec['q_w']:7.4f} "
                  f"{rec['frac_h_zero']:7.4f}")

    print(f"\n== channel ceilings (raw nmae; lower is better; shipped = {base:.6f}) ==")
    for k in ("predict_zero", "as_shipped", "best_global_scalar", "oracle_per_row_scalar",
              "oracle_per_row_offset", "oracle_cpm_stratum_scalar", "lopo_per_gene_offset",
              "oracle_sign_our_magnitude", "our_sign_oracle_magnitude", "oracle_everything"):
        print(f"  {k:32s} {ceil[k]:9.6f}   gain vs shipped {base - ceil[k]:+.6f}")
    print(f"  best global s = {out['best_global_s']:.3f}   "
          f"cpm-stratum scalars {out['cpm_stratum_scalars']}")
    ca = out["cpm_stratum_arm"]
    lo = out["lopo_gene_arm"]
    print(f"\n== the two realisable-shaped arms, priced with their cosine cost on the SAME run ==")
    print(f"  {'arm':28s} {'nmae':>9s} {'d nmae':>9s} {'pds':>9s} {'d pds':>9s} "
          f"{'cos_agg':>9s} {'d cos %':>9s}")
    for nm, d in (("oracle cpm-stratum scalar", ca), ("LOPO per-gene offset", lo)):
        print(f"  {nm:28s} {d['nmae']:9.6f} {base-d['nmae']:+9.6f} {d['pds']:9.6f} "
              f"{d['pds']-b['pds']:+9.6f} {d['cos_agg']:9.6f} "
              f"{100*(d['cos_agg']/b['cos_agg']-1):+9.3f}")
    za = out["zero_audit"]
    print(f"\n== zero / saturation audit on {za['n_gate_pairs']} gate pairs ==")
    print(f"  h exactly 0: {100*za['frac_h_exactly_zero']:.2f}% of pairs, carrying "
          f"{za['l1_share_of_zero_h_genes']:.4f} of the nmae denominator")
    print(f"  |h| > 5 log2: {100*za['frac_h_above_5_log2']:.3f}% of pairs")
    ro = out["row_offsets"]
    print(f"  per-row optimal offset: mean {ro['mean']:+.4f} median {ro['median']:+.4f} "
          f"sd {ro['sd']:.4f}  q05/q95 {ro['q'][0]:+.4f}/{ro['q'][2]:+.4f}")


if __name__ == "__main__":
    main()
