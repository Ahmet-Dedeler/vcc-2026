"""205 gate 3 -- the STRUCTURAL test 182 used to kill pretrained State, applied to X-Cell,
plus the measurable ceiling on what the PROMISED X-Atlas/Pisces upload could be worth.

182's decisive leg was structural: every published Replogle State checkpoint conditions on the
perturbation through `pert_rep: onehot` over a frozen 2,024-gene vocabulary whose intersection
with the VCC 2026 panel is 0/300, and a one-hot has no out-of-vocabulary path.

X-Cell has no checkpoint to read a config off (gate 1 failed), so the vocabulary question is
answered against the artefacts that DO exist: the paper defines its gene/perturbation vocabulary
as "19,082 protein-coding genes, defined by the intersection of GENCODE annotations and the ESM-2
vocabulary", and conditions the perturbation through cross-attention to six prior-knowledge
sources keyed by the perturbed gene. So the measurable question is:

    can those sources address our 300 panel targets at all?

Everything here is a measurement on files already on disk plus one 9 MB HGNC download.
No cell-eval2, no submission, no writes outside this directory.
"""

import json
import pathlib
import sys

import h5py
import numpy as np
import pandas as pd

ROOT = pathlib.Path("~/Code/virtual-cell-challenge-2026")
WORK = ROOT / "experiments/205-xcell-pisces/work"
WORK.mkdir(parents=True, exist_ok=True)
OUT = {}


def say(*a):
    print(*a, flush=True)


# ----------------------------------------------------------------------------- panel
panel = pd.read_csv(ROOT / "data/challenge/pert_counts.csv")["target_gene"].astype(str).tolist()
assert len(panel) == 300 and len(set(panel)) == 300, len(panel)
panel_set = set(panel)
say(f"panel: {len(panel)} targets, first={panel[:3]}")

genes = pd.read_csv(ROOT / "data/challenge/gene_names.csv")["gene_name"].astype(str).tolist()
assert len(genes) == 18533, len(genes)


# ================================================================= 0. KNOWN-ANSWER ROWS
say("\n" + "=" * 78)
say("0. KNOWN-ANSWER ROWS (AGENTS.md rule 10) -- abort if any fails")
say("=" * 78)
ka = {}
fail = []

# KA-1: reproduce 166 sec.1's panel-overlap census off the donor tables themselves.
EXPECT = {"orion_hek293t": 294, "orion_hct116": 290, "kolf2": 273, "k562_gwps": 269,
          "h1_2025": 20, "nadig_jurkat": 0, "rpe1": 0, "k562_essential": 0}
say("\nKA-1  166 sec.1 panel-overlap census (donor tables on disk)")
PERTS = {d: ROOT / "data/donors" / d / "perts.txt" for d in EXPECT}
PERTS["kolf2"] = ROOT / "data/external/kolf2/perts.txt"  # kolf2 lives outside data/donors
for donor, exp in EXPECT.items():
    p = PERTS[donor]
    if not p.exists():
        say(f"  {donor:<16} MISSING {p}")
        fail.append(f"KA-1 {donor} missing")
        continue
    got = len(set(p.read_text().split()) & panel_set)
    ok = got == exp
    say(f"  {donor:<16} expected {exp:>3}  measured {got:>3}   {'OK' if ok else 'FAIL'}")
    ka[f"panel_overlap_{donor}"] = {"expected": exp, "measured": got, "ok": ok}
    if not ok:
        fail.append(f"KA-1 {donor}: {got} != {exp}")

# KA-2: reproduce 197 sec.4.2's two anchor rows from 159's law 1/cos^2(k) = 66.5 + 1920/k.
say("\nKA-2  159 sec.5.5 law reproduces 197 sec.4.2's anchors")


def cos_of_k(k):
    return float(1.0 / np.sqrt(66.5 + 1920.0 / k))


for k, exp in ((3.09, 0.03814), (4.00, 0.04278)):
    got = cos_of_k(k)
    ok = abs(got - exp) < 5e-5
    say(f"  k={k:<5} expected surface cos {exp:.5f}  measured {got:.5f}   {'OK' if ok else 'FAIL'}")
    ka[f"law_k{k}"] = {"expected": exp, "measured": got, "ok": ok}
    if not ok:
        fail.append(f"KA-2 k={k}")

if fail:
    say("\nABORT -- known-answer rows failed: " + "; ".join(fail))
    sys.exit(1)
say("\nall known-answer rows reproduce; numbers below are readable.")


# ============================================ 1. CAN X-CELL'S CONDITIONING ADDRESS THE PANEL?
say("\n" + "=" * 78)
say("1. VOCABULARY / OUT-OF-VOCABULARY TEST  (182's leg 1, asked of X-Cell)")
say("=" * 78)
cov = {}

# 1a. HGNC protein-coding -- the paper's vocabulary is GENCODE protein-coding INTERSECT ESM-2.
hg = pd.read_csv(ROOT / "data/external/xcell/hgnc_protein_coding.txt", sep="\t",
                 low_memory=False, dtype=str)
pc_symbols = set(hg["symbol"].dropna())
alias = set()
for col in ("alias_symbol", "prev_symbol"):
    for v in hg[col].dropna():
        alias.update(x.strip() for x in v.split("|") if x.strip())
hit_direct = panel_set & pc_symbols
hit_alias = (panel_set - hit_direct) & alias
say(f"\n1a. HGNC protein-coding locus group: {len(pc_symbols):,} approved symbols")
say(f"    panel in approved symbols          {len(hit_direct):>3}/300")
say(f"    panel resolvable via alias/prev     {len(hit_alias):>3}/300")
say(f"    panel UNRESOLVABLE                  {len(panel_set - hit_direct - hit_alias):>3}/300"
    f"   {sorted(panel_set - hit_direct - hit_alias)}")
cov["hgnc_protein_coding"] = {"n_vocab": len(pc_symbols), "direct": len(hit_direct),
                              "via_alias": len(hit_alias),
                              "missing": sorted(panel_set - hit_direct - hit_alias)}

# 1b. STRING -- one of the six cross-attention priors, and the paper's own ablation names it
#     (with DepMap) as the two that matter most. On-disk proxy: SLIM's v0.3 symbol-keyed
#     STRING embeddings. NOT the exact v12 h5 X-Cell uses -- labelled as a proxy.
with h5py.File(ROOT / "data/external/slim/gene_string_embeddings.v0.3.h5", "r") as f:
    string_syms = set(f.keys())
say(f"\n1b. STRING (proxy: slim v0.3 symbol-keyed embeddings): {len(string_syms):,} genes")
say(f"    panel covered                      {len(panel_set & string_syms):>3}/300")
say(f"    panel missing                      {sorted(panel_set - string_syms)}")
cov["string_proxy"] = {"n_vocab": len(string_syms), "covered": len(panel_set & string_syms),
                       "missing": sorted(panel_set - string_syms)}

# 1c. DepMap -- the other prior the ablation names. On-disk 26Q1 local build.
dm = pd.read_csv(ROOT / "data/external/depmap26q1/extracted/processed/crispr_table.csv",
                 usecols=["guide"], dtype=str)
dm_syms = set(dm["guide"].dropna().unique())
say(f"\n1c. DepMap 26Q1 (local build): {len(dm_syms):,} genes with a gene_effect row")
say(f"    panel covered                      {len(panel_set & dm_syms):>3}/300")
say("    -> the on-disk build is a 100-gene local subset, NOT DepMap 24Q4 Public"
        " (X-Cell footnote 3). This row is UNINFORMATIVE about X-Cell's DepMap prior."
        " Not downloaded: it cannot change the verdict, because a gene missing from ONE of"
        " six priors is degraded conditioning, not an out-of-vocabulary wall.")
cov["depmap_local"] = {"n_vocab": len(dm_syms), "covered": len(panel_set & dm_syms)}

# 1d. contrast: State's actual released vocabulary, for the side-by-side 182 asks for.
say("\n1d. CONTRAST -- 182's measured figure for published pretrained State:")
say("    pert_rep = onehot over 2,024 perturbations, panel intersection  0/300")

OUT["vocabulary"] = cov


# ============================== 2. THE PANEL-COVERAGE CEILING OF A GENOME-WIDE JURKAT SCREEN
# Pisces' promised upload is: Orion HCT116+HEK293T (held), HepG2+iPSC restricted to a 200-gene
# test set, and the two Jurkat arms in full. Only the Jurkat arms can carry the panel. Their
# libraries are 10,872 / 10,878 targets out of ~18,903 protein-coding genes, and every library
# in 166 sec.1 selects on basal expression in the screened line. Jurkat IS context A (exp 041),
# so the ceiling on that coverage is directly measurable from context_A's own control cells.
say("\n" + "=" * 78)
say("2. CEILING ON A GENOME-WIDE JURKAT SCREEN'S PANEL COVERAGE (measured on context A)")
say("=" * 78)

import anndata as ad  # noqa: E402

A = ad.read_h5ad(ROOT / "data/challenge/context_A.h5ad")
tot = np.asarray(A.X.sum(axis=1)).ravel()
gsum = np.asarray(A.X.sum(axis=0)).ravel()
cpm = gsum / tot.sum() * 1e6
gidx = {g: i for i, g in enumerate(A.var_names.astype(str))}
say(f"context A: {A.shape[0]:,} control cells, median library {np.median(tot):,.0f}")

on_axis = [g for g in panel if g in gidx]
say(f"panel targets on the 18,533-gene readout axis: {len(on_axis)}/300")
pc = np.array([cpm[gidx[g]] for g in on_axis])

# negative control: a random 300-gene draw. Section 2.1 says the panel was selected for
# "sufficiently high basal expression of the target gene", so the panel MUST beat a random draw.
rng = np.random.default_rng(20260903)
rand_cpm = np.stack([cpm[rng.choice(len(cpm), size=300, replace=False)] for _ in range(200)])

say("\n  CPM gate   panel targets expressed   random-300 draw (mean +- sd)")
gate_rows = {}
for gate in (1.0, 5.0, 20.0):
    p_n = int((pc >= gate).sum())
    r_n = (rand_cpm >= gate).sum(axis=1)
    say(f"  >= {gate:>5.1f}   {p_n:>3}/300 ({p_n/3:.1f}%)          "
        f"{r_n.mean():.1f} +- {r_n.std():.1f}")
    gate_rows[str(gate)] = {"panel": p_n, "random_mean": float(r_n.mean()),
                            "random_sd": float(r_n.std())}
ok_neg = gate_rows["5.0"]["panel"] > gate_rows["5.0"]["random_mean"] + 5 * gate_rows["5.0"]["random_sd"]
say(f"\n  negative control (panel must beat a random draw at CPM>=5): {'OK' if ok_neg else 'FAIL'}")
OUT["jurkat_ceiling"] = {"n_on_axis": len(on_axis), "gates": gate_rows,
                         "negative_control_ok": bool(ok_neg)}

# Library-size -> panel-coverage relation, measured on the four panel-covering donors we hold.
say("\n  library size -> panel coverage, measured on the four donors we hold:")
lib_rows = []
for donor in ("orion_hek293t", "orion_hct116", "k562_gwps", "kolf2"):
    perts = set(PERTS[donor].read_text().split())
    lib_rows.append((donor, len(perts), len(perts & panel_set)))
    say(f"    {donor:<16} {len(perts):>6,} targets -> {len(perts & panel_set):>3}/300")
OUT["library_size_relation"] = [{"donor": d, "n_targets": n, "panel": p} for d, n, p in lib_rows]
say("    Pisces Jurkat Resting/Active   10,872 / 10,878 targets -> UNMEASURABLE (not released);"
    " k562_gwps at 9,787 is the nearest held analogue at 269/300")


# ================================= 3. WHAT THE PROMISED UPLOAD WOULD BE WORTH, IF IT LANDED
say("\n" + "=" * 78)
say("3. PROJECTION: value of the promised Pisces upload in cosine / overall units")
say("=" * 78)

# 197 sec.4.2 anchors, all measured, none re-derived here:
OURS_BOARD = 0.079249      # rc_hybrid, measured by inversion, confirmed on the board to 2.5e-7
OURS_SURFACE = 0.03814     # 197 sec.4.2, k_eff 3.09
FOUR_BOARD = 0.0873        # best equal consensus of all four panel-covering lines, k_eff 4.00
LEADER_FLOOR = 0.13117     # best direction on the live board (197 sec.2.6)
RATE = 0.00365             # overall per +0.01 of board cosine (197)
ANCHORS = (1.954, 2.041, 2.078)   # 197 sec.4.1 board<->surface conversion range (abc-knee/a23/a15)

say(f"\n  anchor: ours board {OURS_BOARD:.6f} at k_eff 3.09; four-line consensus {FOUR_BOARD:.4f}"
    f" at k_eff 4.00; leader floor {LEADER_FLOOR:.5f}")
say(f"  board<->surface conversion range {ANCHORS} (197 sec.4.1)\n")
say("  k_eff   surface cos   board cos (a23 anchor)   [range]           vs ours    "
    "overall vs rc_hybrid")
proj = {}
for k in (3.09, 4.00, 5.00, 5.30, 6.00, 10.12):
    cs = cos_of_k(k)
    band = [cs * a for a in ANCHORS]
    cb = cs * 2.041
    say(f"  {k:>5.2f}   {cs:.5f}       {cb:.5f}                 "
        f"[{band[0]:.4f},{band[2]:.4f}]   {cb/OURS_BOARD:>5.2f}x     "
        f"{(cb - OURS_BOARD) * 100 * RATE:+.4f}")
    proj[f"k={k}"] = {"surface_cos": cs, "board_cos_a23": cb,
                      "board_cos_range": band, "ratio_vs_ours": cb / OURS_BOARD,
                      "d_overall_vs_rc_hybrid": (cb - OURS_BOARD) * 100 * RATE}
OUT["projection"] = proj

say(f"\n  AGENTS.md sec.8.2 bar for a single direction arm: +31% cosine")
for label, k in (("+1 line  (Pisces Jurkat, one effective line)", 5.00),
                 ("+1.3 lines (both Jurkat arms, correlated)", 5.30),
                 ("+2 lines (an implausible best case)", 6.00)):
    cb = cos_of_k(k) * 2.041
    say(f"    {label:<44} {cb/OURS_BOARD - 1:+.1%} cosine, "
        f"{'CLEARS' if cb/OURS_BOARD - 1 >= 0.31 else 'MISSES'} the +31% bar")
say(f"    k needed to reach the board's best direction ({LEADER_FLOOR:.5f}): "
    f"{1920 / ((1 / (LEADER_FLOOR / 2.041) ** 2) - 66.5):.1f} effective lines")

(WORK / "gate3.json").write_text(json.dumps(OUT, indent=1))
say(f"\nwrote {WORK / 'gate3.json'}")
