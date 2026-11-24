"""273 -- shared construction for the RELIABLE-REFERENCE instrument: h1-test (2025 H1 hESC, 10x Flex, exactly 400
cells per target, 85 targets) scored with the scorer's own `pds_cosine` statistic, the H1 union OFF.

Nothing mathematical is new.  Everything is imported from the archive:
  * `plolo`        162's instrument (Operator with the 252 registry guard, `deliver` = production `profile_from_effect`
                   + the scorer's read-back, `solve_fraction`, `cd4_for`-style complete-mean3 policy).
  * `lib246`       246's `deflate` (svd / rand modes == the builder's `deflate_cd4_residual`), `merged_many`.
  * `resid218`     `_per_target` -- the scorer's per-target rank statistic (fraction of the other truths FARTHER than
                   the own truth, ties at half; == cell-eval2 `discrimination_score` with rank_denominator n-1,
                   tie_policy midrank, exclusion_scope panel, control_source real -- verified in truth273.py KA-2).
  * `op272`        272's `DampedOperator` (Wiener s^kappa and the eigenvalue floor max(lambda, c n); kappa = 0 / c = 0
                   is `plolo.Operator.apply` bit for bit), 271's `Sums` / `write_half_bank` verbatim.
  * `boot265`      265's paired target bootstrap.

The H1 context (AGENTS.md "TWO TRAPS"): the union sources (`gse295214_h1`, `h1_2025`) cover h1-test's targets and are
NEVER used here -- every arm is base (three production donors, PerturbMatch 0.5 on K562's shared columns, alpha 0.982,
`union_pruned` at gamma 0.65 / lam 0) + CD4 only.  The base rows are 170's staged `cand_h1` (the builder's own
`build_v9_pre_whiten_effects` + `apply_perturbmatch_k562_replacement` for the 85 targets, gated on reproducing a15's
shipped `baseline_pre` / `candidate_pre` sha256), i.e. exactly what 221 whitened, emitted and scored on the cloud.
"""
from __future__ import annotations

import json
import resource
import sys
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/273-flex-reference-instrument"
WORK = HERE / "work"
WORK.mkdir(parents=True, exist_ok=True)
E170 = ROOT / "experiments/170-gamma-mse/work"
E221 = ROOT / "experiments/221-cloud-scorer/cloud-output"
E095 = ROOT / "experiments/095-clean-harness/out-h1test/runs/g0.65_q1.0_a1.08_spread0.25"
E159 = ROOT / "experiments/159-transfer-ceiling"
E272 = ROOT / "experiments/272-noise-corrected-operator"
E271 = ROOT / "experiments/271-donor-reliability-estimation"
E267 = ROOT / "experiments/267-replicate-validated-oracle"
E265 = ROOT / "experiments/265-majority-permutation-test"

for p in (ROOT / "experiments/162-retrodicting-validator", ROOT / "experiments/246-direction-hypotheses", E272, E265,
          E159):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import plolo  # noqa: E402
import lib246 as L  # noqa: E402  (pins plolo.WORK to 246's directory at import; re-pinned below)
import resid218  # noqa: E402
import op272 as O  # noqa: E402
from boot265 import boot  # noqa: E402

# every cache path plolo might touch is pinned INSIDE this experiment; nothing here calls build_fold
plolo.WORK = WORK / "plolo"
plolo.WORK.mkdir(parents=True, exist_ok=True)
L.WORK = WORK
E062, E064 = plolo.E062, plolo.E064
OFFICIAL, PANEL, PANEL_SET = plolo.OFFICIAL, list(plolo.PANEL), plolo.PANEL_SET
TS = plolo.TS                       # 50,000 == the scorer's bulk_target_sum (221 run_params.yaml)
MIN_CPM = plolo.MIN_CPM             # 10.0, a15's render gate
ALPHA, GAMMA, LAM = L.ALPHA, L.GAMMA, L.LAM
CD4_F = L.CD4_F                     # 0.34
CTRL = "non-targeting"
SEED = 20260906
RAND_SEED = L.RAND_SEED             # 246's random-subspace seed
DONORS = ("k562_gwps", "orion_hct116", "orion_hek293t")
CD4_NAMES = plolo.CD4_NAMES
# a15's CD4 fraction solved on the CHALLENGE frame (221 gate 1, `cd4_solved_fraction`) and its f = 0.34 multiple
# (221 gate 3: fraction is exactly linear in f): the fractions 221's cloud-scored h1-test candidates carried.
CHAL_FRACTION_F010 = 0.09472706196346517
CHAL_FRACTION_F034 = CHAL_FRACTION_F010 * 0.34 / 0.10
OFFICIAL_POS = {g: i for i, g in enumerate(OFFICIAL)}
PANEL_POS = np.asarray([OFFICIAL_POS[g] for g in PANEL], dtype=np.int64)


def peak_gib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**3


def sha256_file(p: Path) -> str:
    return E062.sha256_file(Path(p))


def unit(A):
    return A / np.maximum(np.linalg.norm(A, axis=1, keepdims=True), 1e-300)


def row_cos(A, B):
    na, nb = np.linalg.norm(A, axis=1), np.linalg.norm(B, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        c = (A * B).sum(1) / (na * nb)
    return np.where((na > 0) & (nb > 0), c, np.nan)


# ------------------------------------------------------------------------------------------ 170's staged objects
def staged170() -> dict:
    z = np.load(E170 / "pre.npz", allow_pickle=True)
    c = np.load(E170 / "cd4.npz", allow_pickle=True)
    genes = z["genes"].astype(str)
    assert tuple(genes) == tuple(OFFICIAL), "170's gene axis is not the official axis"
    h1_targets = z["h1_targets"].astype(str)
    assert np.array_equal(h1_targets, np.sort(h1_targets)), "170's h1 targets are not sorted (the scorer sorts)"
    return dict(cand_h1=np.asarray(z["cand_h1"], dtype=np.float64), h1_targets=h1_targets,
                h1_genes=z["h1_genes"].astype(str), cols=np.asarray(z["operator_columns"], dtype=np.int64),
                panel_positions=np.asarray(z["panel_positions"], dtype=np.int64),
                h1_positions=np.asarray(z["h1_positions"], dtype=np.int64),
                cd4_frozen=E062.CD4Prepared(mode=str(c["mode"][0]), mean_effect=np.asarray(c["h1_mean"], dtype=np.float64),
                                            gene_votes=np.asarray(c["h1_votes"]), target_condition_count=np.asarray(c["h1_counts"]),
                                            condition_details=[]),
                sha=dict(pre=sha256_file(E170 / "pre.npz"), cd4=sha256_file(E170 / "cd4.npz")))


def load_truth() -> dict:
    """truth273.py's output: the scorer-space truth of h1-test and the render inputs."""
    z = np.load(WORK / "truth273.npz", allow_pickle=True)
    return dict(targets=z["targets"].astype(str), h1_genes=z["h1_genes"].astype(str),
                b_ctrl=np.asarray(z["b_ctrl"], dtype=np.float64), truth=np.asarray(z["truth"], dtype=np.float64),
                scored=np.asarray(z["scored"], dtype=np.int64), gate=np.asarray(z["gate"], dtype=bool),
                sha=str(z["h1_sha256"]))


def to_h1_axis(E_official: np.ndarray, h1_genes) -> tuple[np.ndarray, int]:
    """221 / 170 `to_h1_axis`: 18,533 production columns -> h1-test's 18,080-gene axis (zero where absent)."""
    idx = np.asarray([OFFICIAL_POS.get(str(g), -1) for g in h1_genes], dtype=np.int64)
    out = np.zeros((E_official.shape[0], idx.size), dtype=np.float64)
    have = idx >= 0
    out[:, have] = E_official[:, idx[have]]
    return out, int((~have).sum())


# ------------------------------------------------------------------------------------------ the statistic + score
def pds_own(P_scored: np.ndarray, T_scored: np.ndarray):
    """The scorer's per-target rank statistic and the own cosine (resid218._per_target)."""
    return resid218._per_target(np.asarray(P_scored, dtype=np.float64), np.asarray(T_scored, dtype=np.float64))


def cos_agg(P_scored: np.ndarray, T_scored: np.ndarray) -> float:
    """The energy-weighted aggregate cosine plolo.endpoints reports as sqrt(1 - sse)."""
    a = float((P_scored * P_scored).sum()); x = float((P_scored * T_scored).sum()); d = float((T_scored * T_scored).sum())
    return float(np.sqrt(max(0.0, x * x / (a * d)))) if a > 0 and d > 0 else 0.0


def score_h1(final_official: np.ndarray, tr: dict) -> dict:
    """Render the effect rows as the fold instrument renders (profile_from_effect on h1-test's OWN control profile,
    gate cpm > 10 & b > 0; a23 = minnorm mass fix, a15 = none), read back the way the scorer does, and score the
    per-target `pds` / own cosine on the scored columns (every h1 gene minus the 85 target genes).  Also the
    effect-level (un-rendered) read."""
    E, _ = to_h1_axis(final_official, tr["h1_genes"])
    sc, T = tr["scored"], tr["truth"][:, tr["scored"]]
    out: dict = {}
    for tag, mass in (("a23", "minnorm"), ("a15", None)):
        dl = plolo.deliver(E, tr["b_ctrl"], tr["gate"], mass)
        P = dl[:, sc]
        p, c = pds_own(P, T)
        out[f"pds_{tag}"], out[f"cos_{tag}"] = p, c
        out[f"{tag}_cos_agg"] = cos_agg(P, T)
        out[f"{tag}_delivered_norm"] = float(np.median(np.linalg.norm(dl, axis=1)))
    p, c = pds_own(E[:, sc], T)
    out["pds_eff"], out["cos_eff"] = p, c
    out["eff_cos_agg"] = cos_agg(E[:, sc], T)
    out["pds"], out["cos"] = float(out["pds_a23"].mean()), float(out["cos_a23"].mean())
    return out


def summary(e: dict) -> dict:
    return {k: (float(np.mean(v)) if isinstance(v, np.ndarray) else float(v)) for k, v in e.items()}


# ------------------------------------------------------------------------------------------ composition
def compose(base: np.ndarray, parts, cols: np.ndarray) -> np.ndarray:
    """Whitened base + ONE merged residual add + ONE preserve (E062.add_post_whiten_residual; lib246.compose's
    arithmetic without the panel zeroing: on the H1 context the A/B/C panel genes are ordinary scored genes, and 221's
    cloud-scored candidates carried them as the builder left them)."""
    parts = [(r, f) for r, f in parts if r is not None and f != 0.0]
    if not parts:
        return np.asarray(base, dtype=np.float64).copy()
    final, _ = E062.add_post_whiten_residual(base, L.merged_many(parts), 1.0, cols)
    return final


def complete_mean3(prep) -> "E062.CD4Prepared":
    """plolo.cd4_for's policy (production `freeze_complete_mean3`): zero the effects of rows carried by fewer than
    three conditions, keep their counts."""
    complete = np.asarray(prep.target_condition_count) == len(CD4_NAMES)
    me = np.asarray(prep.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(prep.gene_votes, dtype=np.uint8).copy()
    me[~complete] = 0.0
    votes[~complete] = 0
    return E062.CD4Prepared(mode=prep.mode, mean_effect=me, gene_votes=votes,
                            target_condition_count=np.asarray(prep.target_condition_count, dtype=np.uint8).copy(),
                            condition_details=prep.condition_details)


def deflate_pinned(cd4, V: np.ndarray):
    """The builder's `deflate_cd4_residual_pinned` (251 / 250 `deflate_gene_basis`), verbatim arithmetic: on the
    carried rows, per-gene centring on THIS panel's rows, R_k = R_c - (R_c V) V^T, re-zeroed outside the support."""
    V = np.asarray(V, dtype=np.float64)
    eff = np.asarray(cd4.mean_effect, dtype=np.float64).copy()
    votes = np.asarray(cd4.gene_votes)
    sup = votes > 0
    eff[~sup] = 0.0
    rows = np.flatnonzero(np.asarray(cd4.target_condition_count) > 0)
    R = eff[rows]
    e0 = float((R * R).sum())
    mean = R.mean(axis=0, keepdims=True)
    Rc = R - mean
    Rc[~sup[rows]] = 0.0
    ec = float((Rc * Rc).sum())
    proj = Rc @ V
    captured = float((proj * proj).sum() / max(ec, 1e-300))
    Rk = Rc - proj @ V.T
    Rk[~sup[rows]] = 0.0
    e1 = float((Rk * Rk).sum())
    w = np.linalg.eigvalsh(Rc @ Rc.T)
    own_top_k = float(np.sort(w)[::-1][:V.shape[1]].sum() / max(w.sum(), 1e-300))
    eff[rows] = Rk
    audit = dict(mode="pinned_V5_abc", k=int(V.shape[1]), n_carried=int(rows.size), energy_removed_frac=float(1.0 - e1 / e0),
                 captured_frac_of_centred=captured, own_top_k_energy_frac_of_centred=own_top_k,
                 centred_energy_frac_of_raw=float(ec / e0))
    return L._as_res(cd4, eff, "V5_abc", V.shape[1]), audit


def load_v5_abc() -> tuple[np.ndarray, dict]:
    p = ROOT / "data/donors/cd4_deflation_basis/V5_abc.npz"
    meta = json.loads(p.with_suffix(".json").read_text())
    sha = sha256_file(p)
    assert sha == meta["sha256"], f"V5_abc.npz sha {sha} != json {meta['sha256']}"
    z = np.load(p, allow_pickle=True)
    genes = z["genes"].astype(str)
    assert tuple(genes) == tuple(OFFICIAL), "V5_abc was fitted on a different gene axis"
    V = np.asarray(z["V"], dtype=np.float64)
    return V, dict(sha256=sha, k=int(z["k"]), fitted_on=meta["provenance"]["fitted_on"])


def solve_fraction(whitened, cd4, f: float, *, sizing: str | None = None):
    """plolo.solve_fraction (production's all-rows rule by default; raises on a silent zero)."""
    fr, pu = plolo.solve_fraction(whitened, cd4, f, sizing=sizing)
    return float(fr), float(pu)


# ------------------------------------------------------------------------------------------ groups + bins
BINS_262 = [(-np.inf, -0.02), (-0.02, -0.01), (-0.01, 0.0), (0.0, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.05), (0.05, np.inf)]
PLOLO_262 = {"(-inf, -0.02]": 0.08, "(-0.02, -0.01]": 0.15, "(-0.01, 0]": 0.37, "(0, 0.01]": 0.66, "(0.01, 0.02]": 0.86,
             "(0.02, 0.03]": 0.97, "(0.03, 0.05]": 0.99, "(0.05, inf]": 1.0}


def bin_label(lo, hi) -> str:
    def f(x):
        return "-inf" if x == -np.inf else ("inf" if x == np.inf else f"{x:g}")
    return f"({f(lo)}, {f(hi)}]"


def groups_262(base_cos: np.ndarray) -> dict[str, np.ndarray]:
    """262's groups defined ON THIS REFERENCE from the base's own cosine to the truth: (e) <= 0, marginal (0, 0.05],
    saturated > 0.05."""
    return {"e (base.truth <= 0)": base_cos <= 0, "marginal (0, 0.05]": (base_cos > 0) & (base_cos <= 0.05),
            "saturated > 0.05": base_cos > 0.05}


def terciles(x: np.ndarray) -> dict[str, np.ndarray]:
    q1, q2 = np.nanpercentile(x, [100 / 3, 200 / 3])
    return {"T1 low": x <= q1, "T2 mid": (x > q1) & (x <= q2), "T3 high": x > q2}


def fmt(b: dict) -> str:
    return f"{b['point']:+.4f} [{b['lo']:+.4f}, {b['hi']:+.4f}]"
