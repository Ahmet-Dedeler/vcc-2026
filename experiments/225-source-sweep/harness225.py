"""225 — shared harness.  This is 223's instrument, imported not rewritten.

Everything that prices an arm here comes from modules that already exist:
  * `experiments/162-retrodicting-validator/plolo.py`  — PLOLO folds, operator, scoring
  * `experiments/218-independent-donors/resid218.py`   — source residual, shuffle, scoring
  * `experiments/218-independent-donors/compose218.py` — carrying-row sizing, union, merge
  * `experiments/223-silent-zero-audit/analyse223.py`  — the six-member rate card

The ONLY new thing here is the source registry: 225 registers every donor-shaped table it
builds (in `work/tables/`) plus the ones 223 built, so `source_residual` can read them
unchanged.  The sizing rule is ALWAYS `solve_fraction_carrying` — the mandate's trap.
`plolo.Surface.pds()` never sizes anything (188 Phase 3).
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(ROOT / "experiments/162-retrodicting-validator"))
import plolo  # noqa: E402

plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)

from plolo import (E062, FOLDS, PROD_OPS, Operator, build_fold, cd4_for,  # noqa: E402
                   solve_fraction)

sys.path.insert(0, str(ROOT / "experiments/218-independent-donors"))
import sources218  # noqa: E402

_T223 = "experiments/223-silent-zero-audit/work/tables"
_T225 = "experiments/225-source-sweep/work/tables"

sources218.SOURCES = dict(sources218.SOURCES)


def _reg(name: str, path: str, **kw) -> None:
    sources218.SOURCES[name] = dict(path=path, **kw)


# ---- tables 223 already built (read-only; 225 writes nothing into 223) ---------------
_reg("dld1_lowmoi", f"{_T223}/dld1_lowmoi",
     platform="DLD-1 low-MOI CRISPRi", chemistry="10x", lineage="epithelial (colorectal)",
     transformation="immortalised cancer line", lab="116's acquisition",
     semantics="author log2FC, split-mean", role="223's live source")
for _n in ("dld1_low1", "dld1_low2"):
    _reg(_n, f"{_T223}/{_n}", platform="DLD-1 low-MOI CRISPRi, ONE library",
         chemistry="10x", lineage="epithelial (colorectal)",
         transformation="immortalised cancer line", lab="116's acquisition",
         semantics="author log2FC, single split", role="replication control")
_reg("altair_a549", f"{_T223}/altair_a549", platform="Altair OPS A549", chemistry="10x",
     lineage="epithelial (lung adeno)", transformation="immortalised cancer line",
     lab="Altair", semantics="official OPS log2FC", role="223 death (F1 fail)")
for _n, _lab, _sem in (
        ("cmap_sh_cgs", "LINCS 2020 trt_sh.cgs, 9 core lines", "MODZ consensus z"),
        ("cmap_sh_cgs_lm", "LINCS 2020 trt_sh.cgs, LANDMARK only", "MODZ z, inferred zeroed"),
        ("cmap_xpr", "LINCS 2020 trt_xpr CRISPR-KO", "MODZ consensus z")):
    _reg(_n, f"{_T223}/{_n}", platform="L1000 bulk, 96 h", chemistry="L1000 Luminex",
         lineage="mixed cancer lines", transformation="immortalised cancer lines",
         lab=_lab, semantics=_sem, role="223 death")

# ---- tables 225 builds itself -------------------------------------------------------
NEW_SOURCES = {
    "igvf_cardiomyocyte": dict(
        platform="IGVF TF-Perturb-seq CRISPRi (IGVFFI0830FXFI / IGVFDS6332VCTO)",
        chemistry="10x, published significant-hit DE table",
        lineage="cardiac muscle cell (WTC11-derived)",
        transformation="NON-transformed, differentiated iPSC", lab="Gary Hon, UT (IGVF)",
        semantics="published log2FC x ln2, SIGNIFICANT HITS ONLY (median 6 genes/panel row)"),
    "gasperini_k562": dict(
        platform="Gasperini 2019 CRISPRi (low MOI)", chemistry="10x", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="Gasperini/Shendure 2019",
        semantics="log-pseudobulk delta"),
    "tian2019_ipsc": dict(
        platform="Tian 2019 CRISPRi (iPSC)", chemistry="10x", lineage="pluripotent (iPSC)",
        transformation="NON-transformed iPSC", lab="Tian/Kampmann 2019",
        semantics="log-pseudobulk delta"),
    "tian_ipsc_neuron": dict(
        platform="Tian 2021 CRISPRi (iPSC-derived neuron)", chemistry="10x",
        lineage="neuron (iPSC-derived)", transformation="NON-transformed, differentiated",
        lab="Kampmann 2021", semantics="log-pseudobulk delta"),
    "tian_ipsc_neuron_crispra": dict(
        platform="Tian 2021 CRISPRa (iPSC-derived neuron)", chemistry="10x",
        lineage="neuron (iPSC-derived)", transformation="NON-transformed, differentiated",
        lab="Kampmann 2021", semantics="log-pseudobulk delta (CRISPRa = OPPOSITE sign)"),
    "xucao_hek293": dict(
        platform="Xu 2023 CRISPRi", chemistry="10x", lineage="epithelial (kidney, HEK293)",
        transformation="immortalised", lab="Xu/Cao 2023", semantics="log-pseudobulk delta"),
    "sunshine_calu3": dict(
        platform="Sunshine 2023 CRISPRi (Calu-3)", chemistry="10x",
        lineage="epithelial (lung, Calu-3)", transformation="immortalised cancer line",
        lab="Sunshine/Hein 2023", semantics="log-pseudobulk delta"),
    "norman_k562": dict(
        platform="Norman 2019 CRISPRa (K562)", chemistry="10x", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="Norman/Weissman 2019",
        semantics="log-pseudobulk delta (CRISPRa)"),
    "wessels_thp1": dict(
        platform="Wessels 2023 CRISPRi", chemistry="10x", lineage="myeloid",
        transformation="immortalised cancer line", lab="Wessels/Satija 2023",
        semantics="log-pseudobulk delta"),
    "papalexi_thp1": dict(
        platform="Papalexi 2021 ECCITE-seq", chemistry="10x", lineage="myeloid (THP-1)",
        transformation="immortalised cancer line", lab="Papalexi/Satija 2021",
        semantics="log-pseudobulk delta"),
    "datlinger_jurkat": dict(
        platform="Datlinger 2017 CROP-seq", chemistry="10x/BD", lineage="lymphoid T (Jurkat)",
        transformation="immortalised cancer line", lab="Datlinger/Bock 2017",
        semantics="log-pseudobulk delta"),
    "adamson_k562": dict(
        platform="Adamson 2016 Perturb-seq", chemistry="10x", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="Adamson/Weissman 2016",
        semantics="log-pseudobulk delta"),
    "dixit_k562": dict(
        platform="Dixit 2016 Perturb-seq", chemistry="drop-seq", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="Dixit/Regev 2016",
        semantics="log-pseudobulk delta"),
    "jiang_A549": dict(
        platform="Jiang 2025, A549 ONLY", chemistry="10x", lineage="epithelial (lung adeno)",
        transformation="immortalised cancer line", lab="Jiang 2025",
        semantics="mean log2FC for A549 over 5 stimuli x ln2 — single-line replication control"),
    "jiang_BXPC3": dict(
        platform="Jiang 2025, BXPC3 ONLY", chemistry="10x", lineage="epithelial (pancreas)",
        transformation="immortalised cancer line", lab="Jiang 2025",
        semantics="mean log2FC for BXPC3 over 5 stimuli x ln2 — single-line replication control"),
    "jiang_HAP1": dict(
        platform="Jiang 2025, HAP1 ONLY", chemistry="10x", lineage="near-haploid CML derivative",
        transformation="immortalised cancer line", lab="Jiang 2025",
        semantics="mean log2FC for HAP1 over 5 stimuli x ln2 — single-line replication control"),
    "jiang_HT29": dict(
        platform="Jiang 2025, HT29 ONLY", chemistry="10x", lineage="epithelial (colorectal)",
        transformation="immortalised cancer line", lab="Jiang 2025",
        semantics="mean log2FC for HT29 over 5 stimuli x ln2 — single-line replication control"),
    "jiang_K562": dict(
        platform="Jiang 2025, K562 ONLY", chemistry="10x", lineage="myeloid (CML)",
        transformation="immortalised cancer line", lab="Jiang 2025",
        semantics="mean log2FC for K562 over 5 stimuli x ln2 — single-line replication control"),
    "jiang_MCF7": dict(
        platform="Jiang 2025, MCF7 ONLY", chemistry="10x", lineage="epithelial (breast)",
        transformation="immortalised cancer line", lab="Jiang 2025",
        semantics="mean log2FC for MCF7 over 5 stimuli x ln2 — single-line replication control"),
    "jiang_5line_noK562": dict(
        platform="Jiang 2025, K562 column REMOVED", chemistry="10x",
        lineage="mixed epithelial (A549/BXPC3/HAP1/HT29/MCF7)",
        transformation="immortalised cancer lines", lab="Jiang 2025",
        semantics="mean log2FC over 5 non-K562 lines x 5 stimuli x ln2 — K562 dropped "
                  "because it is the HELD line on the k562_gwps fold"),
    "jiang_6line": dict(
        platform="Jiang 2025 (6 lines pooled)", chemistry="10x",
        lineage="mixed epithelial/myeloid", transformation="immortalised cancer lines",
        lab="Jiang 2025", semantics="author log2FC"),
    "replogle_rpe1_full": dict(
        platform="Replogle 2022 RPE1 genome-wide", chemistry="10x 3' v3",
        lineage="epithelial (retinal pigment)", transformation="hTERT-immortalised primary",
        lab="Weissman 2022", semantics="log-pseudobulk delta"),
    "hcasmc_igvf": dict(
        platform="IGVF HCASMC CRISPRi", chemistry="10x", lineage="smooth muscle (coronary)",
        transformation="primary", lab="IGVF", semantics="log-pseudobulk delta"),
}
for _n, _m in NEW_SOURCES.items():
    _reg(_n, f"{_T225}/{_n}", role="225 candidate", **_m)
# sign-flipped twins of the two CRISPRa sources — activation is the opposite intervention
# to the CRISPRi knockdown every recipient fold measures, so which sign (if either) helps
# has to be a measurement.  `build225.py:flip` writes them.
for _n in ("tian_ipsc_neuron_crispra", "norman_k562"):
    _m = dict(NEW_SOURCES[_n])
    _m["semantics"] += " — SIGN-FLIPPED twin"
    _reg(f"{_n}_flip", f"{_T225}/{_n}_flip", role="225 candidate (CRISPRa flipped)", **_m)

from resid218 import (ALPHA, GAMMA, LAM, score_effects, shuffled,  # noqa: E402
                      source_residual)
from compose218 import (merged, restrict_rows, solve_fraction_carrying,  # noqa: E402
                        union_rows)

CD4_F = 0.34            # board incumbent cd4f034, PINNED (216/217 own the dose)
UNION_G = 0.28          # 218/222's headline dose for the Pan/h1_2025 union, PINNED
UNION = ("gse295214_h1", "h1_2025")
CANCER = ("k562_gwps", "orion_hct116", "orion_hek293t")
PANEL_SET = plolo.PANEL_SET

# 218 §7.4's rate card == 223's, unchanged.  All six members; never four at zero (§8.4).
RATE = {"mse": 0.003634, "reach": 0.002446, "nmae": 0.000980, "fid": 0.000847,
        "jac": 0.000100, "pds_board": -0.000748}
RATE_EX_PDS = RATE["mse"] + RATE["reach"] + RATE["nmae"] + RATE["fid"] + RATE["jac"]
PDS_SPAN = 0.452
GATE_SHAPE_ONLY = 0.0019        # §8.4: +2.5-2.9% of a 0.0792 cosine -> ~+0.0019 overall
GATE_AMPLITUDE = 0.0089


def price(d_cos_abs: float, d_pds_raw: float, k: float = 1.0) -> dict:
    """All six members, 223's function verbatim."""
    cos_gain = d_cos_abs * k
    five = RATE_EX_PDS * (cos_gain / 0.01)
    pds_local = (d_pds_raw / PDS_SPAN) / 6.0 * 1.011 * 2.0
    return dict(five_members=five, pds=pds_local, total=five + pds_local)


def carrying_mask(res) -> np.ndarray:
    return np.asarray(res.target_condition_count) > 0


def audit(res, frac) -> dict:
    added = np.where(np.asarray(res.gene_votes) > 0, np.asarray(res.mean_effect), 0.0)
    norms = np.linalg.norm(added, axis=1)
    live = norms[norms > 0]
    return dict(fraction=float(frac),
                median_row_norm_all_rows=float(np.median(norms)),
                median_row_norm_carrying=float(np.median(live)) if live.size else 0.0,
                carrying_rows=int(live.size), n_rows=int(norms.size),
                coverage_frac=float(live.size / max(norms.size, 1)))


def row_split(cos_per_target, base_cos_per_target, mask) -> dict:
    a = np.asarray(cos_per_target, dtype=np.float64)
    b = np.asarray(base_cos_per_target, dtype=np.float64)
    m = np.asarray(mask, dtype=bool)
    out = {}
    for tag, sel in (("carrying", m), ("non_carrying", ~m)):
        if not sel.any():
            out[f"{tag}_n"] = 0
            continue
        out[f"{tag}_n"] = int(sel.sum())
        out[f"{tag}_d_cos"] = float(a[sel].mean() - b[sel].mean())
        out[f"{tag}_rel_pct"] = float(
            100.0 * (a[sel].mean() - b[sel].mean()) / max(abs(b[sel].mean()), 1e-300))
    return out


class Fold:
    """One PLOLO fold with its operator applied — built once, reused by every arm."""

    def __init__(self, held: str):
        self.held = held
        self.fold = build_fold(held)
        self.targets = self.fold["targets"]
        self.truth = np.asarray(self.fold["truth"], dtype=np.float64)[:, self.fold["scored"]]
        pre = np.asarray(self.fold["pre"], dtype=np.float64)
        self.cd4 = cd4_for(self.targets)
        op_path = (Path(FOLDS[held]["op"]) if FOLDS[held]["op_kind"] != "production"
                   else PROD_OPS / "whitener_union_pruned.npz")
        op = Operator(op_path)
        self.whitened = op.apply(pre * ALPHA, self.fold["wcols"], GAMMA, LAM)
        del op
        self.frac_cd4, _ = solve_fraction(self.whitened, self.cd4, CD4_F)
        self.h1 = union_rows([source_residual(n, self.targets, "gene_norm") for n in UNION])
        self.frac_union, _ = solve_fraction_carrying(self.whitened, self.h1, UNION_G)
        assert self.frac_union > 0.0, "the Pan/h1_2025 union sized to zero — broken instrument"
        self.base = {
            "f0": (None, 0.0),
            "cd4": (merged(self.cd4, self.frac_cd4, None, 0.0), 1.0),
            "base222": (merged(self.cd4, self.frac_cd4, self.h1, self.frac_union), 1.0),
        }
        self.scored = {n: self.score(*v) for n, v in self.base.items()}

    def finish(self, res_or_none, frac):
        if res_or_none is None or frac == 0.0:
            f = self.whitened.copy()
        else:
            f, _ = E062.add_post_whiten_residual(self.whitened, res_or_none, 1.0,
                                                self.fold["wcols"])
        f[:, self.fold["panel_positions"]] = 0.0
        return f

    def score(self, res, frac) -> dict:
        return score_effects(self.finish(res, frac), self.fold, self.truth)

    def arm(self, res, g: float, base: str = "base222") -> dict:
        """Add `res` at carrying-row dose `g` on top of `base`; return the delta dict."""
        m = merged(self.cd4, self.frac_cd4,
                   self.h1 if base == "base222" else None,
                   self.frac_union if base == "base222" else 0.0)
        fr, per_unit = solve_fraction_carrying(self.whitened, res, g)
        if fr > 0.0:
            m = merged(m, 1.0, res, fr)
        e = self.score(m, 1.0)
        b = self.scored[base]
        aud = audit(res, fr)
        aud["per_unit"] = float(per_unit)
        aud["VACUOUS"] = bool(fr == 0.0)
        if fr > 0.0:
            assert aud["median_row_norm_carrying"] > 0.0, "F5 non-vacuity FAILED"
        return dict(
            base=base, g=g,
            d_cos=e["cos"] - b["cos"],
            rel_cos_pct=100.0 * (e["cos"] - b["cos"]) / max(abs(b["cos"]), 1e-300),
            d_cos_agg=e["cos_agg"] - b["cos_agg"],
            d_pds=e["pds"] - b["pds"],
            d_qw=e["a15_sgn_w"] - b["a15_sgn_w"],
            bit_identical=bool(e["cos_agg"] == b["cos_agg"] and e["pds"] == b["pds"]),
            audit=aud,
            rowsplit=row_split(e["_cos_per_target_a23"], b["_cos_per_target_a23"],
                               carrying_mask(res)),
            **{k: v for k, v in e.items() if not k.startswith("_")})
