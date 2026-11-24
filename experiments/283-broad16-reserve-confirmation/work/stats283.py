"""283: the preregistered statistics for claim C2, computed from lolo283_RESERVE.csv (the single access) with
lolo283_KA.csv (the 12 working lines) as the comparison. Paired target bootstrap, 4,000 draws, seed 283."""
import json, numpy as np, pandas as pd
from pathlib import Path
W = Path("~/Code/virtual-cell-challenge-2026/experiments/283-broad16-reserve-confirmation/work")
NB = 4000; SEED = 283
res = pd.read_csv(W / "lolo283_RESERVE.csv"); ka = pd.read_csv(W / "lolo283_KA.csv")
foldsR = json.load(open(W / "lolo283_RESERVE_folds.json")); foldsK = json.load(open(W / "lolo283_KA_folds.json"))
RESERVED = ["JVE127", "MG63", "SLR23", "UMRC3"]
out = {"n_boot": NB, "seed": SEED}

def boot_stat(x, stat=np.median, n=NB, seed=SEED):
    r = np.random.default_rng(seed); x = np.asarray(x, float); m = ~np.isnan(x); x = x[m]
    idx = r.integers(0, len(x), size=(n, len(x))); b = stat(x[idx], axis=1)
    return float(stat(x)), float(np.quantile(b, .025)), float(np.quantile(b, .975))

def boot_paired_mean(d, n=NB, seed=SEED):
    """d: per-target values (already averaged over folds). Returns mean, CI, P(mean <= 0)."""
    r = np.random.default_rng(seed); d = np.asarray(d, float); idx = r.integers(0, len(d), size=(n, len(d)))
    b = d[idx].mean(1)
    return float(d.mean()), float(np.quantile(b, .025)), float(np.quantile(b, .975)), float((b <= 0).mean())

# ---------------- C2a ----------------
print("=" * 110); print("C2a  split-half cosine of the line-specific residual (response - 12-working-line shared response)")
print("     frozen rule: median > 0.028 AND bootstrap 95% lower > 0.028, on >= 3 of 4 reserved lines")
print("-" * 110)
print(f"{'line':9s} {'n':>3s} {'median':>8s} {'95% CI':>18s} {'>0.028?':>8s} | {'in [.0211,.0652]':>17s} {'frac>permnull p95':>18s} {'OR-ctrl med':>12s} {'guide med':>10s} {'UMI/cell':>9s}")
rowsA = []
for l in RESERVED:
    g = res[res.held == l]; m, lo, hi = boot_stat(g.ceiling_cos_AB)
    f = [x for x in foldsR if x["held"] == l][0]
    ok = (m > 0.028) and (lo > 0.028)
    inrange = 0.0211 <= m <= 0.0652
    rowsA.append(dict(line=l, n=len(g), median=m, lo=lo, hi=hi, passes=bool(ok), in_working_range=bool(inrange),
                      frac_above_perm_p95=f["frac_targets_above_null_p95"], perm_null_median=f["null_cos_median"], perm_null_p95=f["null_cos_p95"],
                      or_control_median=f["or_control_cos_median"], guide_median=f["guide_cos_median"], n_guide=f["n_guide_targets"],
                      cells_per_ko_median=f["cells_per_ko_median"]))
    print(f"{l:9s} {len(g):3d} {m:8.4f} [{lo:+.4f},{hi:+.4f}] {'YES' if ok else 'no':>8s} | {'yes' if inrange else 'NO':>17s} {f['frac_targets_above_null_p95']:18.3f} {f['or_control_cos_median']:+12.4f} {f['guide_cos_median']:10.4f}")
nA = sum(r["passes"] for r in rowsA)
print(f"\n  reserved lines passing the frozen 0.028 rule: {nA} of 4  ->  C2a {'PASS' if nA >= 3 else 'FAIL'} (rule needs >= 3)")
kaA = []
for l in sorted(ka.held.unique()):
    g = ka[ka.held == l]; m, lo, hi = boot_stat(g.ceiling_cos_AB); f = [x for x in foldsK if x["held"] == l][0]
    kaA.append(dict(line=l, median=m, lo=lo, hi=hi, passes=bool(m > 0.028 and lo > 0.028), frac_above_perm_p95=f["frac_targets_above_null_p95"],
                    or_control_median=f["or_control_cos_median"], guide_median=f["guide_cos_median"]))
print(f"  same rule on the 12 WORKING lines (disclosed before the access): {sum(r['passes'] for r in kaA)} of 12 pass; medians "
      f"[{min(r['median'] for r in kaA):.4f}, {max(r['median'] for r in kaA):.4f}]")
pooledR = boot_stat(res.ceiling_cos_AB); pooledK = boot_stat(ka.ceiling_cos_AB)
print(f"  pooled median residual cosine: reserved {pooledR[0]:.4f} [{pooledR[1]:.4f}, {pooledR[2]:.4f}] (n={len(res)}) vs working {pooledK[0]:.4f} [{pooledK[1]:.4f}, {pooledK[2]:.4f}] (n={len(ka)})")
dpool = boot_paired_mean(res.groupby("target").ceiling_cos_AB.mean() - res.groupby("target").null_cos_AB.mean())
print(f"  reserved: residual cosine MINUS its own within-line permutation null, paired over targets: {dpool[0]:+.4f} [{dpool[1]:+.4f}, {dpool[2]:+.4f}] P(<=0)={dpool[3]:.4f}")
gR = boot_stat(res.guide_cos_v1v2); gK = boot_stat(ka.guide_cos_v1v2)
print(f"  S4 guide halves (v1 vs v2): reserved pooled median {gR[0]:.4f} [{gR[1]:.4f}, {gR[2]:.4f}], working {gK[0]:.4f}; ratio to cell halves "
      f"{gR[0]/pooledR[0]:.2f}x reserved vs {gK[0]/pooledK[0]:.2f}x working")
orR = [r["or_control_median"] for r in rowsA]; orK = [r["or_control_median"] for r in kaA]
print(f"  S3 OR-control residual cosine (within-line target-independent null): reserved medians {np.round(orR,4).tolist()} (median {np.median(orR):.4f}); working median {np.median(orK):.4f}")
out["C2a"] = {"reserved": rowsA, "working": kaA, "n_pass": int(nA), "verdict": "PASS" if nA >= 3 else "FAIL",
              "pooled_reserved": pooledR, "pooled_working": pooledK, "minus_perm_null_paired": dpool,
              "guide_pooled_reserved": gR, "guide_pooled_working": gK}

# ---------------- C2b ----------------
print("\n" + "=" * 110); print("C2b  cos(12-working-line shared response, reserved line's full response), per-target median")
print("     frozen rule: median inside 268's working-line LOLO per-fold range [0.0400, 0.0937], on >= 3 of 4 lines")
print("-" * 110)
print(f"{'line':9s} {'n':>3s} {'median':>8s} {'95% CI':>18s} {'in range?':>10s} | {'mean':>8s} {'retrieval mean':>14s}")
rowsB = []
for l in RESERVED:
    g = res[res.held == l]; m, lo, hi = boot_stat(g.base_cos_full); f = [x for x in foldsR if x["held"] == l][0]
    inr = 0.0400 <= m <= 0.0937
    rowsB.append(dict(line=l, median=m, lo=lo, hi=hi, in_range=bool(inr), mean=float(g.base_cos_full.mean()), retrieval_mean=f["base_retrieval_mean"]))
    print(f"{l:9s} {len(g):3d} {m:8.4f} [{lo:+.4f},{hi:+.4f}] {'YES' if inr else 'NO':>10s} | {g.base_cos_full.mean():8.4f} {f['base_retrieval_mean']:14.4f}")
nB = sum(r["in_range"] for r in rowsB)
print(f"\n  reserved lines inside the frozen range: {nB} of 4  ->  C2b {'PASS' if nB >= 3 else 'FAIL'} (rule needs >= 3)")
kb = [dict(line=l, median=float(np.median(ka[ka.held == l].base_cos_full)), mean=float(ka[ka.held == l].base_cos_full.mean()),
           retrieval_mean=[x for x in foldsK if x["held"] == l][0]["base_retrieval_mean"]) for l in sorted(ka.held.unique())]
print(f"  working-line comparison: per-fold medians [{min(r['median'] for r in kb):.4f}, {max(r['median'] for r in kb):.4f}], "
      f"retrieval means [{min(r['retrieval_mean'] for r in kb):.3f}, {max(r['retrieval_mean'] for r in kb):.3f}]; "
      f"reserved retrieval means [{min(r['retrieval_mean'] for r in rowsB):.3f}, {max(r['retrieval_mean'] for r in rowsB):.3f}]")
out["C2b"] = {"reserved": rowsB, "working": kb, "n_in_range": int(nB), "verdict": "PASS" if nB >= 3 else "FAIL"}

# ---------------- C2c ----------------
print("\n" + "=" * 110); print("C2c  basal-conditioned model (c), fitted on the 12 WORKING lines, vs (a) = shared response alone")
print("     frozen rule: PASS (268's negative confirmed) iff for BOTH c1 and c2 the 95% CI of dcos_full AND of")
print("     dretrieval covers 0 or lies below 0.  FAIL iff any of those four intervals lies entirely above 0.")
print("-" * 110)
arms = ["c1_target_x_basal_ridge", "c2_target_x_basal_programmes"]
tg = res.groupby("target")
cres = {}
print(f"{'arm':34s} {'metric':12s} {'mean':>9s} {'95% CI':>20s} {'P(<=0)':>8s}   268 working-line value")
ref = {("c1_target_x_basal_ridge", "dcos_full"): "-0.0071 [-0.0093, -0.0050]", ("c1_target_x_basal_ridge", "dretrieval"): "+0.0003 [-0.0097, +0.0110]",
       ("c2_target_x_basal_programmes", "dcos_full"): "+0.0000 [-0.0016, +0.0017]", ("c2_target_x_basal_programmes", "dretrieval"): "-0.0010 [-0.0109, +0.0099]",
       ("c1_target_x_basal_ridge", "cos_res"): "+0.0181 [+0.0134, +0.0229]", ("c2_target_x_basal_programmes", "cos_res"): "+0.0174 [+0.0110, +0.0240]"}
fail = False
for a in arms:
    for met in ("dcos_full", "dretrieval"):
        d = tg[f"{a}|{met}"].mean(); m, lo, hi, P = boot_paired_mean(d)
        cres[f"{a}|{met}"] = [m, lo, hi, P]
        above = lo > 0
        fail |= above
        print(f"{a:34s} {met:12s} {m:+9.4f} [{lo:+.4f},{hi:+.4f}] {P:8.3f}   {ref[(a,met)]}{'   <-- ENTIRELY ABOVE 0' if above else ''}")
print()
for a in arms:
    d = tg[f"{a}|cos_res"].mean(); m, lo, hi, P = boot_paired_mean(d); cres[f"{a}|cos_res"] = [m, lo, hi, P]
    dw = tg[f"{a}__WRONG_recipient|cos_res"].mean(); mw, low, hiw, Pw = boot_paired_mean(dw)
    ds = tg[f"{a}__SHUFFLED_labels|cos_res"].mean(); ms, los, his, Ps = boot_paired_mean(ds)
    dd = boot_paired_mean(d - dw); dd2 = boot_paired_mean(d - ds)
    cres[f"{a}__WRONG|cos_res"] = [mw, low, hiw, Pw]; cres[f"{a}__SHUFFLED|cos_res"] = [ms, los, his, Ps]
    cres[f"{a}|minus_wrong"] = list(dd); cres[f"{a}|minus_shuffled"] = list(dd2)
    print(f"{a:34s} cos_res      {m:+9.4f} [{lo:+.4f},{hi:+.4f}] {P:8.3f}   {ref[(a,'cos_res')]}")
    print(f"{'':34s}   WRONG recip {mw:+9.4f} [{low:+.4f},{hiw:+.4f}]   |  minus WRONG {dd[0]:+.4f} [{dd[1]:+.4f},{dd[2]:+.4f}] P={dd[3]:.3f}")
    print(f"{'':34s}   SHUFFLED    {ms:+9.4f} [{los:+.4f},{his:+.4f}]   |  minus SHUF  {dd2[0]:+.4f} [{dd2[1]:+.4f},{dd2[2]:+.4f}] P={dd2[3]:.3f}")
verdict = "FAIL" if fail else "PASS"
print(f"\n  no interval lies entirely above 0: {not fail}  ->  C2c {verdict} (PASS = 268's negative confirmed on the reserve)")
# per line, for the record
print(f"\n  per reserved line (fold means): {'line':>9s} " + "  ".join(f"{a.split('_')[0]}|dcos  {a.split('_')[0]}|dret" for a in arms))
perline = {}
for l in RESERVED:
    g = res[res.held == l]; perline[l] = {f"{a}|{m}": float(g[f"{a}|{m}"].mean()) for a in arms for m in ("cos_res", "dcos_full", "dretrieval")}
    print(f"{'':33s}{l:>9s} " + "  ".join(f"{g[f'{a}|dcos_full'].mean():+.4f}  {g[f'{a}|dretrieval'].mean():+.4f}" for a in arms))
# the other arms, for the record
other = {}
for a in ("b_recipient_offset", "nn_basal_nearest_line", "oracle_line_offset"):
    other[a] = {m: list(boot_paired_mean(tg[f"{a}|{m}"].mean())) for m in ("cos_res", "dcos_full", "dretrieval")}
    print(f"  {a:26s} cos_res {other[a]['cos_res'][0]:+.4f}  dcos_full {other[a]['dcos_full'][0]:+.4f} [{other[a]['dcos_full'][1]:+.4f},{other[a]['dcos_full'][2]:+.4f}]  dretrieval {other[a]['dretrieval'][0]:+.4f} [{other[a]['dretrieval'][1]:+.4f},{other[a]['dretrieval'][2]:+.4f}]")
out["C2c"] = {"primary": cres, "verdict": verdict, "per_line": perline, "other_arms": other}
out["headline"] = f"C2a {out['C2a']['verdict']}, C2b {out['C2b']['verdict']}, C2c {out['C2c']['verdict']}"
print("\n" + "=" * 110); print(out["headline"]); print("=" * 110)
json.dump(out, open(W / "c2_results.json", "w"), indent=1, default=float)
