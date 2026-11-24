"""269 -- a HELD-LINE-FREE target-level lookup, from 264's matched three-donor arrays (K562 / HCT116 / HEK293T + KOLF2 at
matched amplitude on every cancer fold; the KOLF2 fold with the three production donors).  For held X with donors (Y, Z, K):
the candidate rows are mean(Y,Z) (= drop K, the KOLF2-free incumbent-like row), mean(Z,K) (= drop Y), mean(Y,K) (= drop Z).
Each candidate row's ABSOLUTE per-target `pds` is available on a fold whose truth is NOT X and whose row does not contain X:
  mean(Z,K) scored against Y's truth  = arm `drop:X` on fold Y;   mean(Y,K) against Z's truth = arm `drop:X` on fold Z;
  mean(Y,Z) against KOLF2's truth     = arm `drop:X` on the KOLF2 fold.
So the training scores never touch X's response.  Decision: the row with the best X-free score, acting when it beats the
keep row's score by tau; realised utility = 264m's per-target `pds` of that arm on fold X minus the KOLF2-free row's.

    PYTHONDONTWRITEBYTECODE=1 uv run python xfree269.py
"""
import json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
ROOT = Path(__file__).resolve().parents[2]; WORK = Path(__file__).resolve().parent / "work"
E264 = ROOT / "experiments/264-majority-rule-cancer-folds/work"
FOLD = {"K562": "k562_gwps", "HCT116": "orion_hct116", "HEK293T": "orion_hek293t", "KOLF2": "kolf2"}
ARM = {"K562": "k562_gwps_pm50", "HCT116": "orion_hct116", "HEK293T": "orion_hek293t", "KOLF2": "kolf2"}
CANCER = ["K562", "HCT116", "HEK293T"]
B_BOOT, N_PERM, SEED = 20_000, 500, 20260905
rng = np.random.default_rng(SEED)
LINES = []
def P(s=""):
    print(s, flush=True); LINES.append(s)
def load(line):
    tag = "264m" if line != "KOLF2" else "264"
    pt = np.load(E264 / f"pertarget{tag}-{FOLD[line]}.npz")
    tg = np.load(E264 / f"cov{tag}-{FOLD[line]}.npz", allow_pickle=True)["targets"].tolist()
    return pt, tg
PT = {l: load(l) for l in CANCER + ["KOLF2"]}
def boot(r):
    idx = rng.integers(0, r.size, size=(B_BOOT, r.size)); m = r[idx].mean(1)
    return dict(point=float(r.mean()), lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)), p_le0=float((m <= 0).mean()))
def boot3(rs):
    m = np.zeros(B_BOOT)
    for r in rs:
        idx = rng.integers(0, r.size, size=(B_BOOT, r.size)); m += r[idx].mean(1) / len(rs)
    return dict(point=float(np.mean([r.mean() for r in rs])), lo=float(np.percentile(m, 2.5)), hi=float(np.percentile(m, 97.5)), p_le0=float((m <= 0).mean()))
out = {}
# KA: 264's fold totals for the arms used
for X in CANCER:
    pt, tg = PT[X]
    c = pt["cd4svd5|pds"]; dk = pt[f"drop:{ARM['KOLF2']}:all:nm|pds"]
    P(f"[{X}] cd4svd5 {c.mean():.4f} (264: {dict(K562=0.8154, HCT116=0.8017, HEK293T=0.7810)[X]}); drop-KOLF2 row {dk.mean():.4f} (264 s5: {dict(K562=0.8179, HCT116=0.8015, HEK293T=0.7776)[X]}); "
      f"oracle drop whitened {pt['oracle_drop:wht:all:nm|pds'].mean():.4f} (264: {dict(K562=0.8623, HCT116=0.8475, HEK293T=0.8166)[X]})")
def run(tau, *, baseline="dropK", centre=False, keep_score=True):
    rs, per, acted_all = [], {}, []
    for X in CANCER:
        pt, tg = PT[X]
        Y, Z = [l for l in CANCER if l != X]
        ref = pt[f"drop:{ARM['KOLF2']}:all:nm|pds"] if baseline == "dropK" else pt["cd4svd5|pds"]
        real = {"drop Y": pt[f"drop:{ARM[Y]}:all:nm|pds"] - ref, "drop Z": pt[f"drop:{ARM[Z]}:all:nm|pds"] - ref}
        # X-free training scores, aligned by target name
        def score(line_fold, arm_line):
            p, t2 = PT[line_fold]
            d = dict(zip(t2, p[f"drop:{ARM[arm_line]}:all:nm|pds"]))
            return np.array([d.get(t, np.nan) for t in tg])
        sY = score(Y, X)            # mean(Z,K) vs y
        sZ = score(Z, X)            # mean(Y,K) vs z
        sK = score("KOLF2", X)      # mean(Y,Z) vs KOLF2
        if centre:                  # remove the target's difficulty: each score minus that target's mean over the three X-free scores
            m = np.nanmean(np.column_stack([sY, sZ, sK]), 1)
            sY, sZ, sK = sY - m, sZ - m, sK - m
        n = len(tg); r = np.zeros(n); a = np.zeros(n, bool); act = []
        for t in range(n):
            opts = [(sY[t], "drop Y"), (sZ[t], "drop Z")]
            opts = [(s, k) for s, k in opts if np.isfinite(s)]
            if not opts:
                continue
            s, k = max(opts)
            thr = (sK[t] if (keep_score and np.isfinite(sK[t])) else -np.inf)
            if s - thr > tau if keep_score else s > tau:
                r[t] = real[k][t]; a[t] = True; act.append(k)
        rs.append(r); per[X] = dict(**boot(r), acted_frac=float(a.mean()), precision=float((r[a] > 0).mean()) if a.any() else np.nan,
                                    oracle=float(np.maximum(0, np.maximum(real["drop Y"], real["drop Z"])).mean()))
    return boot3(rs), per, rs
P("\n== held-line-free lookup on 264's matched construction (realised vs the KOLF2-free two-donor row; oracle = best of keep / drop Y / drop Z) ==")
for keep_score in (True, False):
    for centre in (False, True):
        for tau in (0.0, 0.02, 0.05, 0.1):
            b, per, rs = run(tau, centre=centre, keep_score=keep_score)
            key = f"{'vs keep score' if keep_score else 'best drop only'} | {'centred' if centre else 'raw'} | tau {tau}"
            out[key] = dict(mean3=b, per_ctx=per)
            P(f"  {key:44s} 3-ctx {b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}] P(<=0) {b['p_le0']:.3f}  per ctx " + " / ".join(f"{per[c]['point']:+.4f}" for c in CANCER) +
              f"  acted {np.mean([per[c]['acted_frac'] for c in CANCER]):.2f}  precision {np.nanmean([per[c]['precision'] for c in CANCER]):.2f}  oracle {np.mean([per[c]['oracle'] for c in CANCER]):+.4f}")
# the per-target correlation between the X-free score difference and the realised utility difference (the transfer statistic itself)
P("\n  transfer statistic: Spearman between (s(drop Y) - s(drop Z)) [X-free] and (u(drop Y) - u(drop Z)) [held X]:")
for X in CANCER:
    pt, tg = PT[X]; Y, Z = [l for l in CANCER if l != X]
    ref = pt[f"drop:{ARM['KOLF2']}:all:nm|pds"]
    du = pt[f"drop:{ARM[Y]}:all:nm|pds"] - pt[f"drop:{ARM[Z]}:all:nm|pds"]
    def score(line_fold, arm_line):
        p, t2 = PT[line_fold]; d = dict(zip(t2, p[f"drop:{ARM[arm_line]}:all:nm|pds"])); return np.array([d.get(t, np.nan) for t in tg])
    ds = score(Y, X) - score(Z, X)
    ok = np.isfinite(ds)
    rho = spearmanr(ds[ok], du[ok]).correlation
    null = [spearmanr(rng.permutation(ds[ok]), du[ok]).correlation for _ in range(N_PERM)]
    p = float((1 + np.sum(np.abs(null) >= abs(rho))) / (N_PERM + 1))
    out[f"transfer rho {X}"] = dict(n=int(ok.sum()), rho=float(rho), p_perm_two_sided=p)
    P(f"    held {X:8s} n {int(ok.sum())} rho {rho:+.3f} (permutation p {p:.3f}); for reference the LEAKY analogue (scores from fold Y's and Z's 'choose' labels involving X) is +0.50 (lookup269 T2b)")
# and the same with the leaky version in this construction: scores = the per-target utility of the SAME action on the other cancer folds (involves X)
P("\n  leaky analogue in the same construction (action 'drop Y' scored by its utility on the fold where Y... is not the truth but X is a donor):")
for X in CANCER:
    pt, tg = PT[X]; Y, Z = [l for l in CANCER if l != X]
    ref = pt[f"drop:{ARM['KOLF2']}:all:nm|pds"]
    du = pt[f"drop:{ARM[Y]}:all:nm|pds"] - pt[f"drop:{ARM[Z]}:all:nm|pds"]
    # on fold Z (donors X, Y, K): 'drop Y' utility = mean(X,K) - mean(X,Y); on fold Y (donors X, Z, K): 'drop Z' = mean(X,K) - mean(X,Z)
    pz, tz = PT[Z]; py, ty = PT[Y]
    dY = dict(zip(tz, pz[f"drop:{ARM[Y]}:all:nm|pds"] - pz[f"drop:{ARM['KOLF2']}:all:nm|pds"]))
    dZ = dict(zip(ty, py[f"drop:{ARM[Z]}:all:nm|pds"] - py[f"drop:{ARM['KOLF2']}:all:nm|pds"]))
    ds = np.array([dY.get(t, np.nan) - dZ.get(t, np.nan) for t in tg]); ok = np.isfinite(ds)
    rho = spearmanr(ds[ok], du[ok]).correlation
    out[f"leaky rho {X}"] = dict(n=int(ok.sum()), rho=float(rho))
    P(f"    held {X:8s} n {int(ok.sum())} rho {rho:+.3f}")
(WORK / "xfree269.json").write_text(json.dumps(out, indent=1)); (WORK / "xfree269.log").write_text("\n".join(LINES) + "\n")
