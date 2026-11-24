"""275 — re-verify every number quoted in docs/FINDINGS-2026-09.md against the saved JSON/CSV it came from.
Read-only: opens the experiments' saved analysis files and the captured board entries; writes nothing but stdout.
Each row: id | source file : path | stored value | value as quoted in the note | OK (|stored - quoted| <= tol).
Run:  uv run python experiments/275-findings-writeup/work/verify275.py > experiments/275-findings-writeup/work/verify275.log
"""
import csv, json, math, statistics as st, sys
ROOT = "~/Code/virtual-cell-challenge-2026/"
E = ROOT + "experiments/"
A267 = E + "267-replicate-validated-oracle/work/analysis267.json"
A271 = E + "271-donor-reliability-estimation/work/analysis271.json"
D271 = E + "271-donor-reliability-estimation/work/diagwht271.json"
P271 = E + "271-donor-reliability-estimation/work/price271.json"
A272 = E + "272-noise-corrected-operator/work/analysis272.json"
Q272 = E + "272-noise-corrected-operator/work/prod272.json"
A273 = E + "273-flex-reference-instrument/work/analysis273.json"
R273 = E + "273-flex-reference-instrument/work/reliability273.json"
I273 = E + "273-flex-reference-instrument/work/instrument273.json"
A262 = E + "262-where-pds-is-lost/work/analysis262.json"
X262 = E + "262-where-pds-is-lost/work/extra262.json"
B265 = E + "265-majority-permutation-test/work/boot265.json"
D268 = E + "268-broad16-context-residuals/work/decomp268_cpm5.json"
L268 = E + "268-broad16-context-residuals/work/lolo268_K2.json"
A266 = E + "266-k562-de-matched-controls/work/analysis266.json"
R159 = E + "159-transfer-ceiling/work/reliability.json"
T162 = E + "162-retrodicting-validator/work/retrodict.json"
P257 = E + "257-soft-deflation-price/work/price257.json"
BE = ROOT + "data/board-entries/"
B_H1U, B_ARM6, B_SVD5 = BE + "iUziqE1gewfcPKzYUDaB.json", BE + "ydAemVdWmsTx9xhDWcou.json", BE + "Xv0ayvuoHpsKYWv2iu42.json"

_cache = {}
def load(p):
    if p not in _cache:
        with open(p) as f:
            _cache[p] = json.load(f)
    return _cache[p]

def get(p, path):
    o = load(p)
    for tok in path.split("/"):
        if tok.startswith("[") and tok.endswith("]"):
            o = o[int(tok[1:-1])]
        else:
            o = o[tok]
    return o

rows = []
def chk(label, p, path, quoted, tol, cite=""):
    try:
        v = get(p, path)
    except Exception as ex:  # noqa
        rows.append((label, p.replace(ROOT, ""), path, "MISSING", quoted, "??", cite)); return
    ok = "OK" if (isinstance(v, (int, float)) and abs(v - quoted) <= tol) else "MISMATCH"
    rows.append((label, p.replace(ROOT, ""), path, v, quoted, ok, cite))

def val(label, v, quoted, tol, src, cite=""):
    ok = "OK" if abs(v - quoted) <= tol else "MISMATCH"
    rows.append((label, src, "(computed)", v, quoted, ok, cite))

# ---------------- F1 reference reliability (267 §3, 271 §1, 273 §2, 159 §1) ----------------
for line, med, sh, p10, p90, cells, e, m, s, rp, rc in [
    ("k562_gwps", 0.0471, 0.056, 0.0079, 0.162, 81, 0.035, 0.053, 0.121, 0.58, 0.16),
    ("orion_hct116", 0.0276, 0.321, -0.0494, 0.119, 70, 0.014, 0.025, 0.107, 0.23, 0.09),
    ("orion_hek293t", 0.0282, 0.344, -0.122, 0.214, 95, 0.036, 0.032, 0.110, 0.17, 0.12)]:
    b = f"reliability/{line}|batch/"
    chk(f"F1 267 {line} batch median", A267, b + "median", med, 5e-4, "R267 §3")
    chk(f"F1 267 {line} batch share<=0", A267, b + "share_le0", sh, 5e-3, "R267 §3")
    chk(f"F1 267 {line} batch p10", A267, b + "p10", p10, 5e-4, "R267 §3")
    chk(f"F1 267 {line} batch p90", A267, b + "p90", p90, 1e-3, "R267 §3")
    chk(f"F1 267 {line} cells/half", A267, b + "median_cells_half", cells, 0.5, "R267 §3")
    chk(f"F1 267 {line} group (e) reliability", A267, b + "by_group/e_anti", e, 1e-3, "R267 §3 (iii)")
    chk(f"F1 267 {line} marginal reliability", A267, b + "by_group/marginal", m, 1e-3, "R267 §3 (iii)")
    chk(f"F1 267 {line} saturated reliability", A267, b + "by_group/saturated", s, 1e-3, "R267 §3 (iii)")
    chk(f"F1 267 {line} rho(rel, cd4svd5 pds)", A267, b + "rho_pds262", rp, 5e-3, "R267 §3 (iii)")
    chk(f"F1 267 {line} rho(rel, cells)", A267, b + "rho_cells", rc, 5e-3, "R267 §3 (iv)")
chk("F1 267 HEK293T random-cell median", A267, "reliability/orion_hek293t|rand/median", 0.0706, 5e-4, "R267 §3 (v)")
chk("F1 267 HEK293T random-cell share<=0", A267, "reliability/orion_hek293t|rand/share_le0", 0.146, 5e-3, "R267 §3 (v)")
chk("F1 267 K562 guide median (n 8)", A267, "reliability/k562_gwps|guide/median", 0.0103, 5e-4, "R267 §3")
chk("F1 267 HCT116 guide median (n 10)", A267, "reliability/orion_hct116|guide/median", -0.0010, 5e-4, "R267 §3")
chk("F1 267 HEK293T guide median (n 10)", A267, "reliability/orion_hek293t|guide/median", 0.0042, 5e-4, "R267 §3")
chk("F1 267 guide share<=0", A267, "reliability/orion_hct116|guide/share_le0", 0.5, 1e-6, "R267 §3")
for fold, donor, raw, sh, sb, wht in [("k562_gwps", "orion_hct116", 0.027, 0.315, 0.053, -0.141), ("k562_gwps", "orion_hek293t", 0.028, 0.335, 0.054, -0.141),
                                      ("orion_hct116", "k562_gwps", 0.048, 0.058, 0.091, -0.142), ("orion_hct116", "orion_hek293t", 0.028, 0.343, 0.054, -0.243)]:
    b = f"reliability/{fold}|{donor}/"
    chk(f"F1 271 donor {donor} (fold {fold}) raw median", A271, b + "raw/median", raw, 1e-3, "R271 §1.2 D1")
    chk(f"F1 271 donor {donor} share<=0", A271, b + "raw/share_le0", sh, 5e-3, "R271 §1.2 D1")
    chk(f"F1 271 donor {donor} SB full-row median", A271, b + "sb/median", sb, 1e-3, "R271 §1.2 D1")
    chk(f"F3 271 donor {donor} in-sample whitened median", A271, b + "wht/median", wht, 1e-3, "R271 §1.2")
chk("F1 273 Flex batch resid median (200/200)", A273, "reliability/batch|resid|disjoint/median", 0.243, 5e-4, "R273 §2 R")
chk("F1 273 Flex batch resid mean", A273, "reliability/batch|resid|disjoint/mean", 0.331, 5e-4, "R273 §2 R")
chk("F1 273 Flex p10", A273, "reliability/batch|resid|disjoint/p10", 0.067, 5e-4, "R273 §2 R")
chk("F1 273 Flex p90", A273, "reliability/batch|resid|disjoint/p90", 0.748, 5e-4, "R273 §2 R")
chk("F1 273 Flex share<=0", A273, "reliability/batch|resid|disjoint/share_le0", 0.0, 1e-9, "R273 §2 R")
chk("F1 273 Flex share>=0.05", A273, "reliability/batch|resid|disjoint/share_ge005", 0.953, 5e-3, "R273 §2 R")
chk("F1 273 Flex SB of median -> 400", A273, "reliability/batch|resid|disjoint/sb_of_median", 0.391, 1e-3, "R273 §2 (v)")
chk("F1 273 Flex SB of mean -> 400", A273, "reliability/batch|resid|disjoint/sb_of_mean", 0.497, 1e-3, "R273 §2 (v)")
chk("F1 273 Flex full-space median", A273, "reliability/batch|full|disjoint/median", 0.214, 5e-4, "R273 §0.3")
chk("F1 273 Flex full-space share<=0", A273, "reliability/batch|full|disjoint/share_le0", 0.035, 5e-3, "R273 §0.3")
chk("F1 273 Flex random 24 pairings median", A273, "reliability/rand24_m200|resid/median", 0.244, 5e-4, "R273 §2 R")
chk("F1 273 Flex random 24 pairings mean (159 KA-1)", A273, "reliability/rand24_m200|resid/mean", 0.3300, 5e-4, "R273 §2 R, KA-1")
chk("F1 273 Flex matched depth 80/80 median", A273, "reliability/rand24_m80|resid/median", 0.114, 5e-4, "R273 §2 R")
chk("F1 273 Flex matched depth 80/80 mean", A273, "reliability/rand24_m80|resid/mean", 0.206, 5e-4, "R273 §2 R")
chk("F1 273 Flex matched depth share<=0", A273, "reliability/rand24_m80|resid/share_le0", 0.0, 1e-9, "R273 §2 R")
chk("F1 273 guide P1 vs P2 median (n 5)", A273, "reliability/guide|resid/median", 0.0135, 5e-4, "R273 §2 R")
chk("F1 273 guide share<=0", A273, "reliability/guide|resid/share_le0", 0.4, 1e-9, "R273 §2 R")
chk("F1 273 same 5 targets batch median", A273, "reliability/guide_targets_batch|resid/median", 0.243, 5e-4, "R273 §2 R")
chk("F1 273 batch vs random Spearman", A273, "reliability/batch_vs_random24_spearman", 0.994, 1e-3, "R273 §2 (iii)")
chk("F1 273 low tercile cut", A273, "terciles/T1 low/cut", 0.118, 1e-3, "R273 §0.5")
chk("F1 159 h1-test resid disjoint r at m=200 (measured)", R159, "h1test/resid_disjoint/r_mean/[9]", 0.3300, 5e-4, "R159 §1.2")
chk("F1 159 h1-test resid LINEAR-law rho_400 (fitted)", R159, "h1test/resid_disjoint/linear/rho_400", 0.478, 1e-3, "R159 §1.3 quotes 0.481 as 'measured'")
chk("F1 159 h1-test resid POWER-law rho_400 (fitted)", R159, "h1test/resid_disjoint/power/rho_400", 0.458, 1e-3, "R159 §1.3")
chk("F1 159 HCT116 full power rho_400 (the only 0.480)", E + "159-transfer-ceiling/work/crossline.json", "rho150/A_hct116/full/power/rho_400", 0.480, 1e-3, "R273 §2 (v)")
chk("F1 273 KA-1 159 rows reproduced max|d|", R273, "ka/[0]/max_abs_diff", 0.0, 1e-12, "R273 KA-1")
chk("F1 273 KA-2 per-target pds vs cell-eval2 max|d|", R273, "ka/[1]/max_abs_diff_per_target_vs_cloud", 5.6e-17, 1e-17, "R273 KA-2")
chk("F1 273 KA-2 mean pds local", R273, "ka/[1]/mean_local", 0.8655462184873949, 1e-15, "R273 KA-2")
chk("F1 273 KA-2 mean pds scorer", R273, "ka/[1]/mean_scorer", 0.8655462184873949, 1e-15, "R273 KA-2")
chk("F1 273 KA-2a Frobenius ratio", I273, "ka/[0]/value", 1.25e-17, 5e-19, "R273 KA-2a")
chk("F1 273 emission gap a23 mean diff", I273, "emission_gap/cd4frz:f010|a23|vs a15/diff_mean", -0.0025, 1e-4, "R273 §1")
chk("F1 273 emission gap Spearman", I273, "emission_gap/cd4frz:f010|a23|vs a15/spearman", 0.952, 1e-3, "R273 §1")
chk("F1 273 emission gap targets identical", I273, "emission_gap/cd4frz:f010|a23|vs a15/n_equal", 58, 0, "R273 §1")

# ---------------- F2 donor-selection oracle (267 §4, 271 §3, 265) ----------------
O = "oracle/batch|all/"
chk("F2 267 in-sample oracle sym pooled", A267, O + "(ii) in-sample oracle (chosen on the scoring half)/sym/pooled/point", 0.0564, 5e-4, "R267 §4.1")
chk("F2 267 in-sample lo", A267, O + "(ii) in-sample oracle (chosen on the scoring half)/sym/pooled/lo", 0.0499, 5e-4, "R267 §4.1")
chk("F2 267 in-sample hi", A267, O + "(ii) in-sample oracle (chosen on the scoring half)/sym/pooled/hi", 0.0629, 5e-4, "R267 §4.1")
chk("F2 267 cross-validated sym pooled", A267, O + "(iii) CROSS-VALIDATED oracle (chosen on the other half)/sym/pooled/point", -0.0180, 5e-4, "R267 §4.1")
chk("F2 267 cross-validated lo", A267, O + "(iii) CROSS-VALIDATED oracle (chosen on the other half)/sym/pooled/lo", -0.0268, 5e-4, "R267 §4.1")
chk("F2 267 cross-validated hi", A267, O + "(iii) CROSS-VALIDATED oracle (chosen on the other half)/sym/pooled/hi", -0.0094, 5e-4, "R267 §4.1")
chk("F2 267 random pick sym pooled", A267, O + "(iv) random pick, 3 seeds/sym/pooled/point", -0.0225, 5e-4, "R267 §4.1")
chk("F2 267 random lo", A267, O + "(iv) random pick, 3 seeds/sym/pooled/lo", -0.0283, 5e-4, "R267 §4.1")
chk("F2 267 random hi", A267, O + "(iv) random pick, 3 seeds/sym/pooled/hi", -0.0166, 5e-4, "R267 §4.1")
chk("F2 267 cross minus random", A267, O + "(iii)−(iv) cross-validated oracle MINUS random pick (information above a coin flip)/sym/pooled/point", 0.0044, 5e-4, "R267 §4.1")
chk("F2 267 cross minus random lo", A267, O + "(iii)−(iv) cross-validated oracle MINUS random pick (information above a coin flip)/sym/pooled/lo", -0.0045, 5e-4, "R267 §4.1")
chk("F2 267 cross minus random hi", A267, O + "(iii)−(iv) cross-validated oracle MINUS random pick (information above a coin flip)/sym/pooled/hi", 0.0133, 5e-4, "R267 §4.1")
chk("F2 267 cross minus random P(<=0)", A267, O + "(iii)−(iv) cross-validated oracle MINUS random pick (information above a coin flip)/sym/pooled/p_le0", 0.166, 5e-3, "R267 §4.1")
chk("F2 267 retained fraction", A267, O + "retained_fraction/point", -0.32, 5e-3, "R267 §4.1")
chk("F2 267 retained lo", A267, O + "retained_fraction/lo", -0.49, 5e-3, "R267 §4.1")
chk("F2 267 retained hi", A267, O + "retained_fraction/hi", -0.16, 5e-3, "R267 §4.1")
chk("F2 267 abstaining oracle cross-validated", A267, O + "(iii-3) cross-validated abstaining oracle/sym/pooled/point", -0.0102, 5e-4, "R267 §4.1")
chk("F2 267 random-split cross-validated", A267, "oracle/rand|all/(iii) CROSS-VALIDATED oracle (chosen on the other half)/sym/pooled/point", -0.0141, 5e-4, "R267 §4.2")
chk("F2 267 random-split retained", A267, "oracle/rand|all/retained_fraction/point", -0.25, 5e-3, "R267 §4.2")
chk("F2 267 tau0 cross-validated", A267, "oracle/batch|tau0/(iii) CROSS-VALIDATED oracle (chosen on the other half)/sym/pooled/point", -0.0008, 5e-4, "R267 §4.2")
chk("F2 267 T3 cross minus random (batch)", A267, "strata/batch|reliability tercile T3 (high)/cross_minus_random/point", 0.0209, 5e-4, "R267 §4.3")
chk("F2 267 T3 cross minus random lo", A267, "strata/batch|reliability tercile T3 (high)/cross_minus_random/lo", 0.0073, 5e-4, "R267 §4.3")
chk("F2 267 T3 cross minus random hi", A267, "strata/batch|reliability tercile T3 (high)/cross_minus_random/hi", 0.0348, 5e-4, "R267 §4.3")
chk("F2 267 T3 cross-validated vs mean", A267, "strata/batch|reliability tercile T3 (high)/crs/point", -0.0020, 5e-4, "R267 §4.3")
chk("F2 267 T3 cross-validated lo", A267, "strata/batch|reliability tercile T3 (high)/crs/lo", -0.0150, 5e-4, "R267 §4.3")
chk("F2 267 T3 cross-validated hi", A267, "strata/batch|reliability tercile T3 (high)/crs/hi", 0.0114, 5e-4, "R267 §4.3")
chk("F2 267 T3 cross minus random (random split)", A267, "strata/rand|reliability tercile T3 (high)/cross_minus_random/point", 0.0270, 5e-4, "R267 §4.3")
chk("F2 267 T3 random split lo", A267, "strata/rand|reliability tercile T3 (high)/cross_minus_random/lo", 0.0144, 5e-4, "R267 §4.3")
chk("F2 267 T3 random split hi", A267, "strata/rand|reliability tercile T3 (high)/cross_minus_random/hi", 0.0400, 5e-4, "R267 §4.3")
chk("F2 267 group (e) in-sample share", A267, "strata/batch|group (e) base·truth ≤ 0/ins_share", 0.48, 5e-3, "R267 §4.3")
chk("F2 267 group (e) cross-validated", A267, "strata/batch|group (e) base·truth ≤ 0/crs/point", 0.038, 5e-4, "R267 §4.3")
chk("F2 267 group (e) random pick", A267, "strata/batch|group (e) base·truth ≤ 0/rnd_mean", 0.040, 5e-4, "R267 §4.3")
chk("F2 267 pick agreement K562 batch", A267, "agree/k562_gwps/batchA~batchB:all/agree", 0.56, 5e-3, "R267 §4.2 (b)")
chk("F2 267 pick agreement HCT116 batch", A267, "agree/orion_hct116/batchA~batchB:all/agree", 0.58, 5e-3, "R267 §4.2 (b)")
chk("F2 267 pick agreement HEK293T batch", A267, "agree/orion_hek293t/batchA~batchB:all/agree", 0.53, 5e-3, "R267 §4.2 (b)")
chk("F2 271 263's oracle ceiling reproduced", A271, "arms/oracle_wht:all:nm/pooled/point", 0.0459, 5e-4, "R271 §3.1 / R263 §3")
chk("F2 271 oracle ceiling lo", A271, "arms/oracle_wht:all:nm/pooled/lo", 0.0382, 5e-4, "R271 §3.1")
chk("F2 271 oracle ceiling hi", A271, "arms/oracle_wht:all:nm/pooled/hi", 0.0539, 5e-4, "R271 §3.1")
chk("F2 265 kolf2 majority - incumbent", B265, "kolf2/boot/majority - incumbent/point", 0.01209, 5e-5, "R265 §3")
chk("F2 265 kolf2 majority - incumbent lo", B265, "kolf2/boot/majority - incumbent/lo", -0.00346, 5e-5, "R265 §3")
chk("F2 265 kolf2 majority - incumbent hi", B265, "kolf2/boot/majority - incumbent/hi", 0.02783, 5e-5, "R265 §3")
chk("F2 265 kolf2 P(<=0)", B265, "kolf2/boot/majority - incumbent/frac_le0", 0.064, 1e-3, "R265 §3")
chk("F2 265 kolf2 freq-preserving p", B265, "kolf2/perm/majority:raw:all:nm/freq/p_upper", 0.011, 1e-3, "R265 §3")
chk("F2 265 kolf2 freq z", B265, "kolf2/perm/majority:raw:all:nm/freq/z", 2.3, 0.05, "R265 §3")
chk("F2 265 kolf2 null SD", B265, "kolf2/perm/majority:raw:all:nm/freq/null_sd", 0.0068, 1e-4, "R265 §3")
chk("F2 265 kolf2 random-drop mean - incumbent", B265, "kolf2/boot/mean(random drop) - incumbent/point", 0.0068, 1e-4, "R265 §3 / R264 §6")
chk("F2 265 kolf2 three seeds joint p", B265, "kolf2/seeds/all/p_mean_of_three_ge_obs", 0.0047, 1e-4, "R265 §3")
chk("F2 265 cancer K562 freq p", B265, "cancer_264m/folds/k562_gwps/perm/majority:raw:all:nm/freq/p_upper", 0.68, 5e-3, "R265 §4")
chk("F2 265 cancer HCT116 freq p", B265, "cancer_264m/folds/orion_hct116/perm/majority:raw:all:nm/freq/p_upper", 0.12, 5e-3, "R265 §4")
chk("F2 265 cancer HEK293T freq p", B265, "cancer_264m/folds/orion_hek293t/perm/majority:raw:all:nm/freq/p_upper", 0.50, 5e-3, "R265 §4")
chk("F2 265 cancer (b)-(c) pooled", B265, "cancer_264m/pooled/(b) - (c)/point", -0.0107, 5e-4, "R265 §4")
chk("F2 265 cancer (b)-(c) lo", B265, "cancer_264m/pooled/(b) - (c)/lo", -0.0187, 5e-4, "R265 §4")
chk("F2 265 cancer (b)-(c) hi", B265, "cancer_264m/pooled/(b) - (c)/hi", -0.0030, 5e-4, "R265 §4")
chk("F2 265 263 best arbiter union:wht:tau0:nm", B265, "two_donor_263/union:wht:tau0:nm/pooled/point", -0.0004, 5e-5, "R265 §5")

# ---------------- F3 operator noise (271 §1.2, 272, 273 §4.1) ----------------
for op, sz, we, tn in [("k562_gwps", 0.769, 0.875, 0.954), ("orion_hct116", 0.775, 0.877, 0.880), ("orion_hek293t", 0.783, 0.879, 0.902), ("production", 0.787, 0.915, 1.003)]:
    chk(f"F3 272 {op} share s=0 (noise > eigenvalue)", A272, f"spectra/{op}/share_s_zero", sz, 1e-3, "R272 §1.2 S1")
    chk(f"F3 272 {op} whitened energy in noise-dominated", A272, f"spectra/{op}/whitened_energy_share_noise_dominated", we, 1e-3, "R272 §1.2 S1")
    chk(f"F3 272 {op} sum n / sum lambda", A272, f"spectra/{op}/total_noise_over_total_lambda", tn, 1e-3, "R272 §1.2 S1")
chk("F3 272 production share noise-dominated (s<0.5)", A272, "spectra/production/share_noise_dominated", 0.950, 1e-3, "R272 §3")
chk("F3 272 production whitened energy in s=0", A272, "spectra/production/whitened_energy_share_s_zero", 0.698, 1e-3, "R272 §3")
chk("F3 272 in-sample K562 fold / HCT116 donor", A272, "anti/k562_gwps|orion_hct116/wht_median", -0.141, 1e-3, "R272 KA-0 / R271 §1.2")
chk("F3 272 in-sample K562 fold / HEK293T donor", A272, "anti/k562_gwps|orion_hek293t/wht_median", -0.140, 1e-3, "R272 KA-0")
chk("F3 272 in-sample HCT116 fold / HEK293T donor", A272, "anti/orion_hct116|orion_hek293t/wht_median", -0.244, 1e-3, "R272 KA-0")
chk("F3 272 in-sample HEK293T fold / HCT116 donor", A272, "anti/orion_hek293t|orion_hct116/wht_median", -0.241, 1e-3, "R272 KA-0")
chk("F3 272 in-sample share<=0 (HCT116 fold, HEK donor)", A272, "anti/orion_hct116|orion_hek293t/wht_share_le0", 1.0, 1e-9, "R272 §2.2")
chk("F3 272 held K562 through its operator", A272, "anti/k562_gwps|k562_gwps/wht_median", 0.074, 1e-3, "R272 KA-0")
chk("F3 272 held HCT116", A272, "anti/orion_hct116|orion_hct116/wht_median", 0.041, 1e-3, "R272 KA-0")
chk("F3 272 held HEK293T", A272, "anti/orion_hek293t|orion_hek293t/wht_median", 0.029, 1e-3, "R272 KA-0")
chk("F3 272 production operator own donor K562", A272, "anti/production|k562_gwps/wht_median", -0.134, 1e-3, "R272 §3")
chk("F3 272 cross-fit W_A K562 donor (HCT116 fold)", A272, "anti/refit|orion_hct116|k562_gwps/A/median", 0.022, 1e-3, "R272 §2.2 AC2")
chk("F3 272 cross-fit W_A K562 share<=0", A272, "anti/refit|orion_hct116|k562_gwps/A/share_le0", 0.092, 5e-3, "R272 §2.2 AC2")
chk("F3 272 cross-fit W_A HEK293T donor (HCT116 fold)", A272, "anti/refit|orion_hct116|orion_hek293t/A/median", 0.007, 1e-3, "R272 §2.2 AC2")
chk("F3 272 cross-fit W_A HEK293T share<=0", A272, "anti/refit|orion_hct116|orion_hek293t/A/share_le0", 0.221, 5e-3, "R272 §2.2 AC2")
chk("F3 272 cross-fit W_A HCT116 donor (K562 fold)", A272, "anti/refit|k562_gwps|orion_hct116/A/median", 0.011, 1e-3, "R272 §2.3")
chk("F3 272 cross-fit W_A HEK293T donor (K562 fold)", A272, "anti/refit|k562_gwps|orion_hek293t/A/median", 0.009, 1e-3, "R272 §2.3")
chk("F3 272 Wiener k=0.25 pooled", A272, "arms/wiener:k0.25/pooled/point", -0.0285, 5e-4, "R272 §1.3")
chk("F3 272 Wiener k=1 pooled", A272, "arms/wiener:k1/pooled/point", -0.0337, 5e-4, "R272 §1.3")
chk("F3 272 Wiener k=1 lo", A272, "arms/wiener:k1/pooled/lo", -0.0437, 5e-4, "R272 §1.3")
chk("F3 272 Wiener k=1 hi", A272, "arms/wiener:k1/pooled/hi", -0.0238, 5e-4, "R272 §1.3")
chk("F3 272 Wiener k=4 pooled", A272, "arms/wiener:k4/pooled/point", -0.0611, 5e-4, "R272 §1.3")
chk("F3 272 hard cut pooled", A272, "arms/hard:0.5/pooled/point", -0.0693, 5e-4, "R272 §1.3")
chk("F3 272 Wiener k=1 minus shuffled", A272, "vs_controls/wiener:k1/wiener:k1:shuf/point", 0.0318, 5e-4, "R272 §1.3")
chk("F3 272 Wiener k=1 minus shuffled lo", A272, "vs_controls/wiener:k1/wiener:k1:shuf/lo", 0.0151, 5e-4, "R272 §1.3")
chk("F3 272 floor c=0.25 pooled", A272, "arms/floor:c0.25/pooled/point", 0.0001, 5e-5, "R272 §1.3")
chk("F3 272 floor c=1 pooled", A272, "arms/floor:c1/pooled/point", -0.0022, 5e-5, "R272 §1.3")
chk("F3 272 floor c=1 lo", A272, "arms/floor:c1/pooled/lo", -0.0040, 5e-5, "R272 §1.3")
chk("F3 272 floor c=1 hi", A272, "arms/floor:c1/pooled/hi", -0.0004, 5e-5, "R272 §1.3")
chk("F3 272 floor c=1 minus mirrored", A272, "vs_controls/floor:c1/floor:c1:mirror/point", 0.0003, 5e-5, "R272 §1.3")
chk("F3 272 floor c=1 minus shuffled", A272, "vs_controls/floor:c1/floor:c1:shuf/point", -0.0012, 5e-5, "R272 §1.3")
chk("F3 272 floor c=1 whitened energy kept (production)", Q272, "families/floor:c1/whitened_energy_kept", 0.677, 1e-3, "R272 §3 (c = 1: 67.7 %)")
chk("F3 272 floor c=0.5 whitened energy kept (production)", Q272, "families/floor:c0.5/whitened_energy_kept", 0.831, 1e-3, "R272 §3 (c = 0.5: 83.1 %)")
chk("F3 272 floor c=0.25 whitened energy kept (production)", Q272, "families/floor:c0.25/whitened_energy_kept", 0.949, 1e-3, "R272 §3")
chk("F3 272 xfit_avg HCT116 fold", A272, "refit/orion_hct116/arms/xfit_avg|cd4svd5/boot/point", 0.0001, 5e-5, "R272 §2.3 B1")
chk("F3 272 xfit_avg HCT116 lo", A272, "refit/orion_hct116/arms/xfit_avg|cd4svd5/boot/lo", -0.0033, 5e-5, "R272 §2.3 B1")
chk("F3 272 xfit_avg HCT116 hi", A272, "refit/orion_hct116/arms/xfit_avg|cd4svd5/boot/hi", 0.0035, 5e-5, "R272 §2.3 B1")
chk("F3 272 proper cross-fit HCT116 fold", A272, "refit/orion_hct116/arms/half:xfit|half:inc/boot/point", 0.0003, 5e-5, "R272 §2.3 B1")
chk("F3 272 proper cross-fit HCT116 lo", A272, "refit/orion_hct116/arms/half:xfit|half:inc/boot/lo", -0.0050, 5e-5, "R272 §2.3 B1")
chk("F3 272 proper cross-fit HCT116 hi (R272 prints +0.0057; stored rounds to +0.0056)", A272, "refit/orion_hct116/arms/half:xfit|half:inc/boot/hi", 0.0056, 5e-5, "R272 §2.3 B1")
chk("F3 272 proper cross-fit HCT116 rel cos_agg", A272, "refit/orion_hct116/arms/half:xfit|half:inc/rel_agg", -0.0328, 5e-4, "R272 §2.3 B1")
chk("F3 272 cross-fit minus in-fit HCT116", A272, "refit/orion_hct116/arms/half:xfit|half:infit/boot/point", -0.0009, 5e-5, "R272 §2.3 B1")
chk("F3 272 cross-fit minus in-fit HCT116 rel cos_agg", A272, "refit/orion_hct116/arms/half:xfit|half:infit/rel_agg", -0.0645, 5e-4, "R272 §2.3 B1")
chk("F3 272 xfit_avg K562 fold", A272, "refit/k562_gwps/arms/xfit_avg|cd4svd5/boot/point", 0.0003, 5e-5, "R272 §2.3")
chk("F3 272 proper cross-fit K562 fold", A272, "refit/k562_gwps/arms/half:xfit|half:inc/boot/point", 0.0022, 5e-5, "R272 §2.3")
chk("F3 272 proper cross-fit K562 lo", A272, "refit/k562_gwps/arms/half:xfit|half:inc/boot/lo", -0.0034, 5e-5, "R272 §2.3")
chk("F3 272 proper cross-fit K562 hi", A272, "refit/k562_gwps/arms/half:xfit|half:inc/boot/hi", 0.0078, 5e-5, "R272 §2.3")
chk("F3 272 proper cross-fit K562 rel cos_agg", A272, "refit/k562_gwps/arms/half:xfit|half:inc/rel_agg", -0.0373, 5e-4, "R272 §2.3")
chk("F3 272 cross-fit minus in-fit K562", A272, "refit/k562_gwps/arms/half:xfit|half:infit/boot/point", 0.0051, 5e-5, "R272 §2.3")
chk("F3 273 floor c=1 pds", A273, "a/floor:c1 vs INCUMBENT/pds/point", -0.0001, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 pds lo", A273, "a/floor:c1 vs INCUMBENT/pds/lo", -0.0053, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 pds hi", A273, "a/floor:c1 vs INCUMBENT/pds/hi", 0.0043, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 own cosine", A273, "a/floor:c1 vs INCUMBENT/cos/point", 0.0024, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 own cosine lo", A273, "a/floor:c1 vs INCUMBENT/cos/lo", 0.0014, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 own cosine hi", A273, "a/floor:c1 vs INCUMBENT/cos/hi", 0.0034, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 rel cos_agg %", A273, "a/floor:c1 vs INCUMBENT/cos_agg_pct", 4.4, 0.05, "R273 §4.1")
chk("F3 273 floor c=1 without CD4 vs BASE pds", A273, "a/floor:c1:noCD4 vs BASE/pds/point", 0.0043, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 without CD4 lo", A273, "a/floor:c1:noCD4 vs BASE/pds/lo", -0.0010, 5e-5, "R273 §4.1")
chk("F3 273 floor c=1 without CD4 hi", A273, "a/floor:c1:noCD4 vs BASE/pds/hi", 0.0099, 5e-5, "R273 §4.1")
chk("F3 273 Wiener k=1 pds", A273, "a/wiener:k1 vs INCUMBENT/pds/point", -0.0126, 5e-5, "R273 §4.1")
chk("F3 273 Wiener k=1 lo", A273, "a/wiener:k1 vs INCUMBENT/pds/lo", -0.0352, 5e-5, "R273 §4.1")
chk("F3 273 Wiener k=1 hi", A273, "a/wiener:k1 vs INCUMBENT/pds/hi", 0.0119, 5e-5, "R273 §4.1")
chk("F3 273 Wiener k=1 rel cos_agg %", A273, "a/wiener:k1 vs INCUMBENT/cos_agg_pct", -7.7, 0.05, "R273 §4.1")
chk("F3 273 Wiener k=1 group (e) mean", A273, "a/wiener:k1 vs INCUMBENT/groups/e (base.truth <= 0)/[1]", 0.097, 1e-3, "R273 §4.1")
chk("F3 273 Wiener k=1 marginal mean", A273, "a/wiener:k1 vs INCUMBENT/groups/marginal (0, 0.05]/[1]", -0.041, 1e-3, "R273 §4.1")
chk("F3 273 Wiener k=1 minus shuffled", A273, "controls/wiener:k1 - wiener:k1:shuf/point", 0.050, 1e-3, "R273 §4.1")
chk("F3 273 Wiener k=1 minus shuffled lo", A273, "controls/wiener:k1 - wiener:k1:shuf/lo", 0.007, 1e-3, "R273 §4.1")
chk("F3 271 EB row shrinkage pooled", A271, "arms/ebrow/pooled/point", 0.0014, 5e-5, "R271 §3.1")
chk("F3 271 EB row shrinkage lo", A271, "arms/ebrow/pooled/lo", 0.0005, 5e-5, "R271 §3.1")
chk("F3 271 EB row shrinkage hi", A271, "arms/ebrow/pooled/hi", 0.0024, 5e-5, "R271 §3.1")
chk("F3 271 EB minus shuffled", A271, "vs_controls/ebrow/ebrow:shuf/point", 0.0013, 5e-5, "R271 §3.2")
chk("F3 271 EB minus shuffled lo", A271, "vs_controls/ebrow/ebrow:shuf/lo", -0.0001, 5e-5, "R271 §3.2")
chk("F3 271 EB minus shuffled hi", A271, "vs_controls/ebrow/ebrow:shuf/hi", 0.0027, 5e-5, "R271 §3.2")
chk("F3 271 EB leave-top-3-out retention", A271, "retention/ebrow/retained", 0.31, 5e-3, "R271 §3.5")
chk("F3 271 EB pessimistic price", P271, "table/ebrow/scenarios/PESSIMISTIC/total", 0.00050, 5e-6, "R271 §4")
chk("F3 271 reliability-weighted mean pooled", A271, "arms/relw/pooled/point", -0.0037, 5e-5, "R271 §3.1")
chk("F3 271 reliability-weighted lo", A271, "arms/relw/pooled/lo", -0.0103, 5e-5, "R271 §3.1")
chk("F3 271 reliability-weighted hi", A271, "arms/relw/pooled/hi", 0.0028, 5e-5, "R271 §3.1")
chk("F3 271 relw minus reversed", A271, "vs_controls/relw/relw:rev/point", 0.0157, 5e-5, "R271 §3.2")
chk("F3 271 a13 re-implemented pooled", A271, "arms/a13:nm/pooled/point", -0.0069, 5e-5, "R271 §3.1")
chk("F3 273 EB row shrinkage pds (Flex)", A273, "b/ebrow vs mean:rebuilt/pds/point", 0.0000, 5e-5, "R273 §4.2")
chk("F3 273 EB pds lo", A273, "b/ebrow vs mean:rebuilt/pds/lo", -0.0022, 5e-5, "R273 §4.2")
chk("F3 273 EB pds hi", A273, "b/ebrow vs mean:rebuilt/pds/hi", 0.0025, 5e-5, "R273 §4.2")
chk("F3 273 EB own cosine (Flex)", A273, "b/ebrow vs mean:rebuilt/cos/point", 0.0009, 5e-5, "R273 §4.2")
chk("F3 273 EB shuffled control (Flex)", A273, "b/ebrow:shuf vs mean:rebuilt/pds/point", 0.0014, 5e-5, "R273 §4.2")
chk("F3 273 relw reversed (Flex)", A273, "b/relw:rev vs mean:rebuilt/pds/point", -0.0396, 5e-5, "R273 §4.2")
chk("F3 273 deflation on vs off pds (Flex)", A273, "c/INCUMBENT vs cd4h1:off/pds/point", 0.0006, 5e-5, "R273 §4.3")
chk("F3 273 deflation lo", A273, "c/INCUMBENT vs cd4h1:off/pds/lo", -0.0116, 5e-5, "R273 §4.3")
chk("F3 273 deflation hi", A273, "c/INCUMBENT vs cd4h1:off/pds/hi", 0.0125, 5e-5, "R273 §4.3")
chk("F3 273 deflation rel cos_agg %", A273, "c/INCUMBENT vs cd4h1:off/cos_agg_pct", -4.8, 0.05, "R273 §4.3")
chk("F3 273 deflation group (e) mean", A273, "c/INCUMBENT vs cd4h1:off/groups/e (base.truth <= 0)/[1]", -0.060, 1e-3, "R273 §4.3")
chk("F3 273 deflation marginal mean", A273, "c/INCUMBENT vs cd4h1:off/groups/marginal (0, 0.05]/[1]", 0.013, 1e-3, "R273 §4.3")
chk("F3 273 deflation minus shuffled residual", A273, "controls/INCUMBENT - cd4h1:shuf/point", 0.032, 1e-3, "R273 §4.3")
chk("F3 273 deflation minus shuffled lo", A273, "controls/INCUMBENT - cd4h1:shuf/lo", 0.002, 1e-3, "R273 §4.3")
chk("F3 273 deflation minus random subspace", A273, "controls/INCUMBENT - cd4h1:randE5/point", 0.010, 1e-3, "R273 §4.3")
chk("F3 273 CD4 residual vs BASE rel cos_agg %", A273, "c/cd4h1:off vs BASE/cos_agg_pct", 30.3, 0.05, "R273 §4.3")
chk("F3 273 CD4 residual vs BASE own cosine", A273, "c/cd4h1:off vs BASE/cos/point", 0.0087, 5e-5, "R273 §4.3")
chk("F3 273 CD4 residual vs BASE pds", A273, "c/cd4h1:off vs BASE/pds/point", 0.0116, 5e-5, "R273 §4.3")

# ---------------- F4 conversion curve and deficit (262 §4–5, 273 §3) ----------------
for i, (b, m, d) in enumerate([("(-1,-0.02]", 0.079, 0.064), ("(-0.02,-0.01]", 0.147, 0.164), ("(-0.01,0]", 0.374, 0.343), ("(0,0.01]", 0.659, 0.284),
                                ("(0.01,0.02]", 0.860, 0.116), ("(0.02,0.03]", 0.970, 0.023), ("(0.03,0.05]", 0.993, 0.005), ("(0.05,0.08]", 0.998, 0.001)]):
    chk(f"F4 262 PLOLO bin {b} mean pds", X262, f"conversion_curve/[{i}]/mean_pds", m, 1e-3, "R262 §4.2")
    chk(f"F4 262 PLOLO bin {b} deficit share", X262, f"conversion_curve/[{i}]/deficit_share", d, 1e-3, "R262 §4.2")
chk("F4 262 PLOLO slope inside (0, 0.05]", X262, "slope_band", 9.5, 0.05, "R262 §4.2")
for i, (b, m) in enumerate([("(-0.02,-0.01]", 0.162), ("(-0.01,0]", 0.448), ("(0,0.01]", 0.762), ("(0.01,0.02]", 0.923), ("(0.02,0.03]", 0.979), ("(0.03,0.05]", 0.996), ("(0.05,inf]", 0.999)], start=1):
    chk(f"F4 273 Flex bin {b} mean pds", A273, f"conversion/curve/INCUMBENT/[{i}]/[2]", m, 1e-3, "R273 §3 C")
chk("F4 273 Flex slope inside (0, 0.05]", A273, "conversion/slope_band", 6.0, 0.05, "R273 §3")
chk("F4 273 Flex share below chance", A273, "conversion/regimes_delivered/below chance (own cos <= 0)", 0.129, 1e-3, "R273 §3")
chk("F4 273 Flex share marginal", A273, "conversion/regimes_delivered/marginal (0, 0.05]", 0.48, 5e-3, "R273 §3")
chk("F4 273 Flex share saturated", A273, "conversion/regimes_delivered/saturated > 0.05", 0.39, 5e-3, "R273 §3")
chk("F4 273 Flex share at pds = 1", A273, "conversion/share_pds1", 0.55, 5e-3, "R273 §3")
chk("F4 273 Flex mean pds", A273, "conversion/mean_pds", 0.883, 1e-3, "R273 §3")
chk("F4 262 group (e) n", X262, "objects_by_group/base·truth <= 0 (oracle group e)/n", 149, 0, "R262 §5.2 T5d")
chk("F4 262 group (e) cd4svd5 pds", X262, "objects_by_group/base·truth <= 0 (oracle group e)/cd4svd5", 0.389, 1e-3, "R262 §5.2 T5d")
chk("F4 262 marginal n", X262, "objects_by_group/base·truth in (0, 0.05] (the marginal band)/n", 564, 0, "R262 §5.2 T5d")
chk("F4 262 marginal cd4svd5 pds", X262, "objects_by_group/base·truth in (0, 0.05] (the marginal band)/cd4svd5", 0.858, 1e-3, "R262 §5.2 T5d")
chk("F4 262 saturated n", X262, "objects_by_group/base·truth > 0.05/n", 140, 0, "R262 §5.2 T5d")
chk("F4 262 saturated cd4svd5 pds", X262, "objects_by_group/base·truth > 0.05/cd4svd5", 0.999, 1e-3, "R262 §5.2 T5d")
try:
    ne, pe = get(X262, "objects_by_group/base·truth <= 0 (oracle group e)/n"), get(X262, "objects_by_group/base·truth <= 0 (oracle group e)/cd4svd5")
    nm, pm = get(X262, "objects_by_group/base·truth in (0, 0.05] (the marginal band)/n"), get(X262, "objects_by_group/base·truth in (0, 0.05] (the marginal band)/cd4svd5")
    ns, ps = get(X262, "objects_by_group/base·truth > 0.05/n"), get(X262, "objects_by_group/base·truth > 0.05/cd4svd5")
    de, dm, ds = ne * (1 - pe), nm * (1 - pm), ns * (1 - ps)
    tot = de + dm + ds
    val("F4 262 group (e) share of targets", ne / (ne + nm + ns), 0.175, 1e-3, "extra262.json objects_by_group", "R262 §5.2 T5d")
    val("F4 262 group (e) share of deficit", de / tot, 0.531, 2e-3, "extra262.json objects_by_group", "R262 §5.2 T5d")
    val("F4 262 marginal share of deficit", dm / tot, 0.468, 2e-3, "extra262.json objects_by_group", "R262 §5.2 T5d")
    val("F4 262 saturated share of deficit", ds / tot, 0.001, 2e-3, "extra262.json objects_by_group", "R262 §5.2 T5d")
except Exception as ex:
    rows.append(("F4 262 deficit shares", "extra262.json", "(computed)", "ERR " + str(ex), "", "??", ""))
for fold, e, r, d in [("k562_gwps", 0.575, 0.117, 0.087), ("orion_hct116", 0.515, 0.161, 0.106), ("orion_hek293t", 0.490, 0.181, 0.123)]:
    chk(f"F4 262 partition (e) share {fold}", A262, f"decomposition/{fold}/partition/shares/e_base_away", e, 1e-3, "R262 §5.2 T5c")
    chk(f"F4 262 partition (r) share {fold}", A262, f"decomposition/{fold}/partition/shares/r_resid_breaks", r, 1e-3, "R262 §5.2 T5c")
    chk(f"F4 262 partition render share {fold}", A262, f"decomposition/{fold}/partition/shares/d_render", d, 1e-3, "R262 §5.2 T5c")
chk("F4 262 base·truth Spearman pooled", A262, "covariate_ranking/cov|base_cos_truth/rho_pooled", 0.875, 1e-3, "R262 §4.1")
chk("F4 262 base·truth AUC pooled", A262, "covariate_ranking/cov|base_cos_truth/auc_pooled", 0.050, 1e-3, "R262 §4.1")
chk("F4 262 held-vs-donor agreement Spearman", A262, "covariate_ranking/agree_held_vs_donors/rho_pooled", 0.630, 1e-3, "R262 §4.1")
chk("F4 262 donor-donor agreement Spearman", A262, "covariate_ranking/agree_donor_donor/rho_pooled", 0.381, 1e-3, "R262 §4.1")
chk("F4 262 truth norm Spearman", A262, "covariate_ranking/cd4svd5|truth_norm/rho_pooled", -0.029, 1e-3, "R262 §4.1")
chk("F4 262 prediction hubness Spearman", A262, "covariate_ranking/cd4svd5|pred_mean/rho_pooled", -0.034, 1e-3, "R262 §4.1")
chk("F4 262 bottom-decile overlap K562/HCT116", A262, "agreement/k562_gwps|orion_hct116/overlap", 9, 0, "R262 §3.2 T2")
chk("F4 262 bottom-decile overlap K562/HEK293T", A262, "agreement/k562_gwps|orion_hek293t/overlap", 12, 0, "R262 §3.2 T2")
chk("F4 262 bottom-decile overlap HCT116/HEK293T", A262, "agreement/orion_hct116|orion_hek293t/overlap", 12, 0, "R262 §3.2 T2")
chk("F4 262 expected overlap (HCT116/HEK293T)", A262, "agreement/orion_hct116|orion_hek293t/expected", 2.2, 0.05, "R262 §3.2 T2")
chk("F4 262 Spearman pds across folds K562/HCT116", A262, "agreement/k562_gwps|orion_hct116/spearman_pds", 0.756, 1e-3, "R262 §3.2 T2")
chk("F4 262 Spearman pds across folds K562/HEK293T", A262, "agreement/k562_gwps|orion_hek293t/spearman_pds", 0.705, 1e-3, "R262 §3.2 T2")
chk("F4 262 Spearman pds across folds HCT116/HEK293T", A262, "agreement/orion_hct116|orion_hek293t/spearman_pds", 0.741, 1e-3, "R262 §3.2 T2")
chk("F4 262 cd4svd5 pds K562", A262, "distribution/k562_gwps/mean", 0.8154, 5e-5, "R262 §3.1 T1")
chk("F4 262 cd4svd5 pds HCT116", A262, "distribution/orion_hct116/mean", 0.8017, 5e-5, "R262 §3.1 T1")
chk("F4 262 cd4svd5 pds HEK293T", A262, "distribution/orion_hek293t/mean", 0.7810, 5e-5, "R262 §3.1 T1")

# ---------------- F5 Broad-16 (268 §2–4) ----------------
chk("F5 268 shared fraction (pooled, full depth)", D268, "cell_halves/pooled_full_depth/shared", 0.091, 1e-3, "R268 §2")
chk("F5 268 line-specific fraction", D268, "cell_halves/pooled_full_depth/line", 0.115, 1e-3, "R268 §2")
chk("F5 268 noise fraction", D268, "cell_halves/pooled_full_depth/noise", 0.794, 1e-3, "R268 §2")
chk("F5 268 line-specific direction-changing", D268, "cell_halves/pooled_full_depth_line_perp", 0.099, 1e-3, "R268 §2")
chk("F5 268 line-specific amplitude-only", D268, "cell_halves/pooled_full_depth_line_par", 0.019, 1e-3, "R268 §2")
chk("F5 268 per-target median line-specific", D268, "cell_halves/median_full_depth/line", 0.080, 1e-3, "R268 §2")
chk("F5 268 per-target median shared", D268, "cell_halves/median_full_depth/shared", 0.031, 1e-3, "R268 §2")
chk("F5 268 OR-control null line-specific", D268, "or_control_null/pooled_full_depth/line", 0.028, 1e-3, "R268 §2")
chk("F5 268 OR-control null shared", D268, "or_control_null/pooled_full_depth/shared", 0.019, 1e-3, "R268 §2")
chk("F5 268 guide halves line-specific", D268, "guide_halves/pooled_full_depth/line", 0.091, 1e-3, "R268 §2")
chk("F5 268 guide halves shared", D268, "guide_halves/pooled_full_depth/shared", 0.080, 1e-3, "R268 §2")
chk("F5 268 residual split-half cosine median", D268, "cell_halves/cos_res_median", 0.046, 1e-3, "R268 §2")
chk("F5 268 residual (t,l) pairs above null p95", D268, "cell_halves/frac_tl_above_null95", 0.262, 1e-3, "R268 §2")
chk("F5 268 within-batch residual similarity", D268, "batch/within_batch_mean_sim", 0.020, 1e-3, "R268 §2")
chk("F5 268 between-batch residual similarity", D268, "batch/between_batch_mean_sim", -0.106, 1e-3, "R268 §2")
chk("F5 268 batch permutation p", D268, "batch/perm_p", 0.001, 1e-9, "R268 §2")
chk("F5 268 residual components for 50 %", D268, "svd/residual/k_for_50_80_90/[0]", 27, 0, "R268 §3")
chk("F5 268 residual components for 80 %", D268, "svd/residual/k_for_50_80_90/[1]", 167, 0, "R268 §3")
chk("F5 268 shared components for 50 %", D268, "svd/shared/k_for_50_80_90/[0]", 3, 0, "R268 §3")
chk("F5 268 residual components above p99 (of 60)", D268, "svd/residual/n_sig_components_p99", 58, 0, "R268 §3")
chk("F5 268 C2 share (proliferation vs p53)", D268, "components/[1]/cross_energy_share", 0.135, 1e-3, "R268 §3")
chk("F5 268 C1 share", D268, "components/[0]/cross_energy_share", 0.094, 1e-3, "R268 §3")
chk("F5 268 LOLO ceiling median", L268, "ceiling_cos_AB_median", 0.045, 1e-3, "R268 §4")
chk("F5 268 LOLO ceiling mean", L268, "ceiling_cos_AB_mean", 0.061, 1e-3, "R268 §4")
chk("F5 268 LOLO target-only full cosine", L268, "base_cos_full_mean", 0.113, 1e-3, "R268 §4")
chk("F5 268 LOLO target-only retrieval", L268, "base_retrieval_mean", 0.766, 1e-3, "R268 §4")
chk("F5 268 c1 residual cosine", L268, "arms/c1_target_x_basal_ridge/cos_res", 0.0181, 5e-4, "R268 §4")
chk("F5 268 c1 residual cosine lo", L268, "arms/c1_target_x_basal_ridge/cos_res_ci/[0]", 0.0134, 5e-4, "R268 §4")
chk("F5 268 c1 residual cosine hi", L268, "arms/c1_target_x_basal_ridge/cos_res_ci/[1]", 0.0229, 5e-4, "R268 §4")
chk("F5 268 c1 delta full cosine", L268, "arms/c1_target_x_basal_ridge/dcos_full", -0.0071, 5e-4, "R268 §4")
chk("F5 268 c1 delta full lo", L268, "arms/c1_target_x_basal_ridge/dcos_full_ci/[0]", -0.0093, 5e-4, "R268 §4")
chk("F5 268 c1 delta full hi", L268, "arms/c1_target_x_basal_ridge/dcos_full_ci/[1]", -0.0050, 5e-4, "R268 §4")
chk("F5 268 c1 delta retrieval", L268, "arms/c1_target_x_basal_ridge/dretrieval", 0.0003, 5e-4, "R268 §4")
chk("F5 268 c1 delta retrieval lo", L268, "arms/c1_target_x_basal_ridge/dretrieval_ci/[0]", -0.0098, 5e-4, "R268 §4")
chk("F5 268 c1 delta retrieval hi", L268, "arms/c1_target_x_basal_ridge/dretrieval_ci/[1]", 0.0105, 5e-4, "R268 §4")
chk("F5 268 c2 residual cosine", L268, "arms/c2_target_x_basal_programmes/cos_res", 0.0174, 5e-4, "R268 §4")
chk("F5 268 c2 delta full cosine", L268, "arms/c2_target_x_basal_programmes/dcos_full", 0.0000, 5e-4, "R268 §4")
chk("F5 268 c2 delta retrieval", L268, "arms/c2_target_x_basal_programmes/dretrieval", -0.0010, 5e-4, "R268 §4")
chk("F5 268 c1 WRONG recipient residual cosine", L268, "arms/c1_target_x_basal_ridge__WRONG_recipient/cos_res", 0.0012, 5e-4, "R268 §4")
chk("F5 268 c1 SHUFFLED labels residual cosine", L268, "arms/c1_target_x_basal_ridge__SHUFFLED_labels/cos_res", -0.0073, 5e-4, "R268 §4")
chk("F5 268 oracle line offset residual cosine", L268, "arms/oracle_line_offset/cos_res", 0.297, 1e-3, "R268 §4")
chk("F5 268 oracle line offset delta full cosine", L268, "arms/oracle_line_offset/dcos_full", 0.169, 1e-3, "R268 §4")
chk("F5 268 oracle line offset delta retrieval", L268, "arms/oracle_line_offset/dretrieval", -0.019, 1e-3, "R268 §4")
chk("F5 268 nearest basal line delta full cosine", L268, "arms/nn_basal_nearest_line/dcos_full", -0.049, 1e-3, "R268 §4")

# ---------------- F6 DE sparsity (260 §2c, 266) ----------------
R260 = E + "260-raw-de-diagnostic-feasibility/work/"
def csvcol(fn, col):
    return [float(r[col]) for r in csv.DictReader(open(R260 + fn))]
v = csvcol("de_power_panel_targets.csv", "n_real_tc")
val("F6 260 panel targets n", len(v), 272, 0, "de_power_panel_targets.csv", "R260 §2c")
val("F6 260 panel median n_real (tie-corrected)", st.median(v), 5, 0, "de_power_panel_targets.csv", "R260 §2c")
val("F6 260 panel targets >= 10 genes", sum(x >= 10 for x in v), 85, 0, "de_power_panel_targets.csv", "R260 §2c")
val("F6 260 panel targets with zero genes", sum(x == 0 for x in v), 58, 0, "de_power_panel_targets.csv", "R260 §2c")
lfc = [float(r["target_knockdown_lfc"]) for r in csv.DictReader(open(R260 + "de_power_panel_targets.csv")) if r["target_in_axis"] in ("True", "1", "true")]
val("F6 260 targets with own gene on axis", len(lfc), 229, 0, "de_power_panel_targets.csv", "R260 §2c")
val("F6 260 median target-gene LFC", st.median(lfc), -3.05, 0.01, "de_power_panel_targets.csv", "R260 §2c")
val("F6 260 share of target LFC negative", sum(x < 0 for x in lfc) / len(lfc), 1.0, 1e-9, "de_power_panel_targets.csv", "R260 §2c")
c = csvcol("candidates_de_power.csv", "n_real_tc")
val("F6 260 candidates median n_real", st.median(c), 6, 0, "candidates_de_power.csv", "R260 §2c")
val("F6 260 candidates >= 10 genes", sum(x >= 10 for x in c), 17, 0, "candidates_de_power.csv", "R260 §2c")
val("F6 260 candidates >= 89 genes", sum(x >= 89 for x in c), 4, 0, "candidates_de_power.csv", "R260 §2c")
curve = {(r["target"], int(r["cells"])): int(r["n_real_tc"]) for r in csv.DictReader(open(R260 + "de_power_curve.csv"))}
val("F6 260 MYRF at 400 cells", curve[("MYRF", 400)], 7, 0, "de_power_curve.csv", "R260 §2c")
val("F6 260 NABP2 at 400 cells", curve[("NABP2", 400)], 9, 0, "de_power_curve.csv", "R260 §2c")
val("F6 260 PINK1 at 400 cells", curve[("PINK1", 400)], 6, 0, "de_power_curve.csv", "R260 §2c")
ka = {(r["target"], int(r["cells"])): int(r["n_real_tc"]) for r in csv.DictReader(open(R260 + "ka_strong_targets.csv"))}
val("F6 260 KA RPS3 n_real", ka[("RPS3", 163)], 1205, 0, "ka_strong_targets.csv", "R260 §2c KA")
val("F6 260 KA HSPA5 n_real", ka[("HSPA5", 299)], 425, 0, "ka_strong_targets.csv", "R260 §2c KA")
val("F6 260 KA LMO2 n_real", ka[("LMO2", 273)], 342, 0, "ka_strong_targets.csv", "R260 §2c KA")
chk("F6 260 control median UMI", E + "260-raw-de-diagnostic-feasibility/work/de_power_summary.json", "ctrl_median_umi", 11406, 0, "R260 §2c")
for s, med, g in [("pooled_s0", 6.0, 17), ("pooled_s1", 6.0, 18), ("pooled_s2", 6.0, 18), ("gem_s0", 6.0, 18), ("gem_s1", 6.0, 18), ("gem_s2", 5.5, 17), ("gemdepth_s0", 6.5, 18), ("gemdepth_s1", 6.5, 18), ("gemdepth_s2", 6.0, 18)]:
    chk(f"F6 266 {s} median n_real", A266, f"summary/{s}/median", med, 1e-9, "R266 §3")
    chk(f"F6 266 {s} targets >= 10", A266, f"summary/{s}/ge10", g, 0, "R266 §3")
chk("F6 266 pooled gemgroup TV distance", A266, "diagnostics/pooled/tv", 0.39, 5e-3, "R266 §4")
chk("F6 266 gem-matched TV distance", A266, "diagnostics/gem/tv", 0.005, 5e-4, "R266 §4")
chk("F6 266 pooled control UMI CV", A266, "diagnostics/pooled/cv", 0.440, 1e-3, "R266 §4")
chk("F6 266 depth-matched control UMI CV", A266, "diagnostics/gemdepth/cv", 0.378, 1e-3, "R266 §4")
chk("F6 266 KA 260 counts reproduced (50/50)", A266, "ka_exact/n_real_tc", 50, 0, "R266 §2")

# ---------------- F7 board calibrations (BR237, BR249, 162, 257) ----------------
for name, p, pds, mse, nmae, avg, spds, smse, snmae in [
    ("h1union", B_H1U, 0.83572, 0.81726, 0.93389, 0.2118310, 0.74537, 0.17957, 0.11074),
    ("ARM6", B_ARM6, 0.84238, 0.81227, 0.96016, 0.2081907, 0.76027, 0.18482, 0.06778),
    ("cd4svd5", B_SVD5, 0.85064, 0.82181, 0.93397, 0.2160944, 0.77809, 0.17477, 0.11047)]:
    chk(f"F7 board {name} raw pds", p, "entry/pds_cosine", pds, 5e-6, "BR237/BR249 §1")
    chk(f"F7 board {name} raw mse", p, "entry/expr_mse_unbiased_capped_norm", mse, 5e-6, "BR237/BR249 §1")
    chk(f"F7 board {name} raw nmae", p, "entry/de_wilcoxon_lfc_nmae", nmae, 5e-6, "BR237/BR249 §1")
    chk(f"F7 board {name} overall", p, "entry/score_avg", avg, 5e-8, "BR237/BR249")
    chk(f"F7 board {name} scaled pds", p, "entry/score_pds", spds, 5e-6, "BR237/BR249 §1")
    chk(f"F7 board {name} scaled mse", p, "entry/score_mse", smse, 5e-6, "BR237/BR249 §1")
    chk(f"F7 board {name} scaled nmae", p, "entry/score_nmae", snmae, 5e-6, "BR237/BR249 §1")
h, a6, s5 = load(B_H1U)["entry"], load(B_ARM6)["entry"], load(B_SVD5)["entry"]
val("F7 ARM6 overall delta vs h1union", a6["score_avg"] - h["score_avg"], -0.0036403, 5e-8, "board entries", "BR237")
val("F7 ARM6 nmae contribution to overall", (a6["score_nmae"] - h["score_nmae"]) / 6, -0.00716, 5e-6, "board entries", "BR237 §1")
val("F7 ARM6 pds contribution", (a6["score_pds"] - h["score_pds"]) / 6, 0.00248, 5e-6, "board entries", "BR237 §1")
val("F7 ARM6 mse contribution", (a6["score_mse"] - h["score_mse"]) / 6, 0.00088, 5e-6, "board entries", "BR237 §1")
val("F7 ARM6 raw pds delta", a6["pds_cosine"] - h["pds_cosine"], 0.00666, 5e-6, "board entries", "BR237 §1")
val("F7 ARM6 nmae realised / identity -0.004783", ((a6["score_nmae"] - h["score_nmae"]) / 6) / -0.004783, 1.50, 0.005, "board entries + BR237 §2", "BR237 §2")
val("F7 cd4svd5 overall delta vs h1union", s5["score_avg"] - h["score_avg"], 0.0042634, 5e-8, "board entries", "BR249")
val("F7 cd4svd5 raw pds delta", s5["pds_cosine"] - h["pds_cosine"], 0.01492, 5e-6, "board entries", "BR249 §1")
val("F7 cd4svd5 pds contribution", (s5["score_pds"] - h["score_pds"]) / 6, 0.00545, 5e-6, "board entries", "BR249 §1")
val("F7 cd4svd5 raw mse delta", s5["expr_mse_unbiased_capped_norm"] - h["expr_mse_unbiased_capped_norm"], 0.00454, 5e-6, "board entries", "BR249 §1")
val("F7 cd4svd5 mse contribution", (s5["score_mse"] - h["score_mse"]) / 6, -0.00080, 5e-6, "board entries", "BR249 §1")
val("F7 cd4svd5 nmae contribution", (s5["score_nmae"] - h["score_nmae"]) / 6, -0.00005, 5e-6, "board entries", "BR249 §1")
chk("F7 257 PLOLO cd4svd5 - h1union pds", P257, "retrodiction/R1/compare/d_pds", 0.012263, 5e-6, "R257 R2")
val("F7 k_pds re-shape class = board raw / PLOLO", (s5["pds_cosine"] - h["pds_cosine"]) / get(P257, "retrodiction/R1/compare/d_pds"), 1.22, 0.005, "board entries + price257.json", "BR249 §1")
val("F7 k_pds new-source class (ARM6, A/B/C functional +0.002687)", (a6["pds_cosine"] - h["pds_cosine"]) / 0.002687, 2.48, 0.005, "board entries + R258 §3 (PLOLO +0.002687)", "R258 §3")
val("F7 k_pds new-source class (ARM6, frame-A functional +0.001755)", (a6["pds_cosine"] - h["pds_cosine"]) / 0.001755, 3.79, 0.005, "board entries + R258 §3 (PLOLO +0.001755)", "R258 §3")
chk("F7 257 board dcos for cd4svd5 (197 identity)", P257, "retrodiction/R1/price/point/dcos_board", -0.001863, 5e-6, "R257 R1 / BR249 §2")
val("F7 k_cos re-shape class = board dcos / (PLOLO rel cos_agg x cos)", get(P257, "retrodiction/R1/price/point/dcos_board") / (-0.00461748 * 0.107332), 3.76, 0.01, "price257.json + R257 R2 (-0.461748 %) + 243's 0.107332", "BR249 §2")
val("F7 k_cos new-source class (ARM6): +0.003791 / (3.6728 % x 0.107332)", 0.003791 / (0.036728 * 0.107332), 0.962, 0.005, "R258 §3 numbers", "BR237 §2 / R258 §3")
chk("F7 257 fitted-factor retrodiction of BR249 (point)", P257, "retrodiction/R1/price/point/total", 0.00472, 5e-6, "R257 R1")
chk("F7 257 conservative-factor retrodiction", P257, "retrodiction/R1/price/conservative/total", 0.00257, 5e-6, "R257 R1")
chk("F7 271 price recipe retrodiction of BR249", P271, "retrodiction/price/total", 0.00472, 5e-6, "R271 §4")
chk("F7 162 PLOLO LOO RMS on pds (scaled)", T162, "pds/loo_rms", 0.001494, 5e-7, "R162 verdict")
chk("F7 162 predict-zero null RMS", T162, "pds/null_rms", 0.003481, 5e-7, "R162 verdict")
chk("F7 162 fitted scale", T162, "pds/k_in", 0.899, 5e-4, "R162 verdict")
chk("F7 162 sign hits of 6", T162, "pds/sign_hits", 6, 0, "R162 verdict")
chk("F7 162 permutation p", T162, "pds/perm_p", 0.0028, 5e-5, "R162 verdict")
chk("F7 162 mse LOO RMS (rejected)", T162, "mse/loo_rms", 0.1124, 5e-5, "R162 §2")
chk("F7 162 reach sign hits of 6", T162, "reach/sign_hits", 2, 0, "R162 §2")

# ---------------- print ----------------
def fmt(v):
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)
n_ok = sum(r[5] == "OK" for r in rows); n_bad = sum(r[5] != "OK" for r in rows)
print(f"# verify275 — {len(rows)} numbers checked: {n_ok} OK, {n_bad} not OK\n")
print("| id | source : path | stored | quoted | ok | RESULTS section |")
print("|---|---|---:|---:|---|---|")
for r in rows:
    print(f"| {r[0]} | `{r[1]}` : `{r[2]}` | {fmt(r[3])} | {fmt(r[4])} | {r[5]} | {r[6]} |")
if n_bad:
    print("\n## not OK\n")
    for r in rows:
        if r[5] != "OK":
            print(f"- {r[0]}: stored {fmt(r[3])} vs quoted {fmt(r[4])} ({r[1]} : {r[2]})")
