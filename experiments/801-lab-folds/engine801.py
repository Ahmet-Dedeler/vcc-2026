"""801 -- the lab's bulk fold engine: prices ANY bulk arm on a held-out PLOLO fold from that fold's dump, member by member,
per target, exactly as production's arithmetic would (dump708's / 522's path). Light (~1-2.5 GB per fold), no operator.

VENDORED (arithmetic verbatim; only paths, guards and the extra outputs listed below are new):
  * fold522.py  (experiments/522-october-bulk-key, git blob 070ffac13e49ffae5bfe416ce4c816669fad77fc): `pct_hi`, `wclip`,
    P_PROD / P_546, the Oct loader (here `DumpFold`), and `Abc` (imported from 522 itself: score518.Fold's path).
  * sp2_522.py  (git blob 1c89e5c30c42d37d318d769c8bac7168e38bfd44): `score_q` (per-target q; rows sharpened per q group),
    `terc`, P2. Here `render` is score_q's loop body, returning the realised rows and K as well.
  * arms522c.py (git blob f7a42e63f0415876d88a28c5ff1eddf58c5beeb0): `p2g_oct` (== lib551.q_p2g, the carried-only
    tercilec form the October recipe ships) and `p2g_abc` (arms540c's terc form, used ONLY for 540c's known answers).
The arithmetic underneath (key_scale, sharpen, profile_from_effect, realise_k, row_score) is lib801 (lib438's copy) for dump
folds and lib420 / comp420 for A/B/C folds, exactly as in 522.

A FOLD  is a PLOLO fold: `load(fold_id)` with fold_id "<panel>/<held line>" (L1/k562_gwps ... from work/registry.json;
        efp/<held> = the KA fold built from 238's efp inputs; ABC/<held> = A/B/C's own panel through 518's dumps, held in
        k562_gwps / orion_hct116 / orion_hek293t / kolf2). Attributes: targets, n, kv (the fold's LODO-V1 key), carried,
        strength (S_P2's statistic), frames {A,B,C: (b, gate)}, cpm {A,B,C: control CPM per gene}, scored, nonpanel.
AN ARM  is `Arm(w, q, fin=None)`: w = per-target key weights (n), q = scalar or per-target sharpening exponents (n),
        fin = optional replacement composition rows (n x 18533, the space of the dump's h2c1 `fin`; a candidate that blends
        or swaps external rows per target builds them with `blend_rows`). `named(F, name)` gives the known arms:
        cap3, cap6, pct, pct+P2g (= INCUMBENT), pct+S_P2, cap3+S_P2, cap3+P2g, cap6+P2g, cap6+S_P2 (+ the KA forms).
OUTPUT  `price(F, arm, base=None) -> (summary, arrays)`; arrays (numpy, frames x targets unless noted):
        x = <P, T>, a = <P, P> (scored genes; P = realised log-ratio row, T = the held line's truth), d = <T, T> (targets),
        pds = literal PDS per target (the scorer's rank rule: share of other targets' truths less cosine-similar, ties 1/2),
        pdsl = board-emulated PDS (414's beta at lambda 0.5), K = ballast K (KMAX 16),
        and, when `base` is given, the change energy |P_arm - P_base|^2 per target: dE_scored, dE_unscored (non-panel genes
        off the held line's scored axis), dE_panel (the panel's own columns), dE_cpm0..3 (all non-panel genes, by the
        frame's control CPM bin [0,5), [5,20), [20,100), [100,inf)); dE_total = dE_scored + dE_unscored + dE_panel.
        summary: cos_agg (mean over frames of sum x / sqrt(sum a sum d)), pds, pds_lam0.5, pds_lam1, K_mean, AB (render
        energies), per_frame; with a base: r = cos_agg / cos_agg_base - 1, dpds, dpds_lam0.5, dK_mean, the energy shares.
KNOWN ANSWERS on load (abort unless |d| < 1e-12):
  * every fold: KA-PT0, the per-target path at a constant q 1.75 == fold522's scalar path (cap3 and pct), bit for bit;
  * efp/<held>: 546o's `base q1.75`, `fixed 6`, `p96.6`, `p96.6 + P2g`, `p96.6 + S_P2`; 551's `h2c1|V1_5_3+sh1.75`,
    `V1pct+sh1.75`, `V1pct+ptP2g` (the incumbent), `V1pct+ptS_P2`, `V1_5_3+ptP2g`, `V1_5_3+ptS_P2` (L.cmp: cos_agg, pds,
    pds lam 0.5 / 1, render A / B per frame); 522's arms522 `cap3`, `cap6`, `pctPROD` (cos_agg, pds, K_mean);
  * ABC/<held>: fold522.Abc's own (arms540's V1_0.5_3 / V1_0.5_6), 540c's `V1_0.5_6 + P2g` / `+ S_P2` / `V1_0.5_3 + P2g` /
    `+ S_P2` (cos_agg, pds; L.cmp where 540c stores AB), 522's arms522 `cap3`, `cap6`, `pctPROD`.

    uv run python experiments/801-lab-folds/engine801.py ka <fold_id>...      # load with every KA, print the table
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import guard801  # noqa: E402,F401
import json  # noqa: E402
import time  # noqa: E402
from dataclasses import dataclass  # noqa: E402

import numpy as np  # noqa: E402

import lib801 as L  # noqa: E402
from lib801 import SLOTS  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

WORK = HERE / "work"
REGISTRY = WORK / "registry.json"
W546 = ROOT / "experiments/546-october-cap-pertarget/work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
W551 = ROOT / "experiments/551-r6-c10b/work"
W522 = ROOT / "experiments/522-october-bulk-key/work"
TOL = 1e-12
P_PROD = 96.64697522769474          # 497's unrounded P (== cap 6 on the A/B/C artifact); fold522.P_PROD
P_546 = 96.647                      # the rounded P 543/546 priced (known-answer rows); fold522.P_546
P2 = (2.5, 2.0, 1.75)               # LOW, MID, HIGH (479 / 546 / lib551.P2)
FLOOR = 0.5
CPM_EDGES = (5.0, 20.0, 100.0)      # DESIGN section 7 / 702 rule 3: [0,5), [5,20), [20,100), [100,inf)


# ---------------------------------------------------------------------------------------------- fold522 (verbatim)
def pct_hi(kv, carried, p, mask=None):
    """543's percentile cap: the p-th percentile (numpy linear) of key / median over the carried targets."""
    c = carried if mask is None else (carried & mask)
    kn = kv[c] / float(np.median(kv[c]))
    return float(np.percentile(kn, p))


def wclip(kv, carried, lo, hi, mask=None):
    """key410.weights (clip(key / median_carried, lo, hi) on carried, 1 elsewhere), median over the kept targets."""
    c = carried if mask is None else (carried & mask)
    med = float(np.median(kv[c]))
    w = np.ones(len(kv))
    w[carried] = np.clip(kv[carried] / med, lo, hi)
    return w


# ---------------------------------------------------------------------------------------------- sp2_522 / arms522c (verbatim)
def terc(v, spec):
    t = np.digitize(v, np.quantile(v, [1 / 3, 2 / 3]))
    return np.asarray([spec[k] for k in t], float)


def p2g_oct(W, carried):
    """comp546o's P2g: carried-only cut points, floor -> LOW, non-carried -> MID (== lib551.q_p2g; the October tercilec)."""
    fl = W <= FLOOR + 1e-12
    ec = np.quantile(W[carried], [1 / 3, 2 / 3])
    tk = np.where(fl, 0, np.where(~carried, 1, np.digitize(W, ec)))
    return np.asarray([P2[k] for k in tk], float)


def p2g_abc(W):
    """arms540c's P2g: terc(W, P2, floor 0.5). Used only to reproduce 540c's A/B/C rows."""
    t = np.digitize(W, np.quantile(W, [1 / 3, 2 / 3]))
    t = np.where(W <= FLOOR + 1e-12, 0, t)
    return np.asarray([P2[k] for k in t], float)


# ---------------------------------------------------------------------------------------------- folds
def cpm_of(b):
    expr = np.expm1(b)
    return 1e6 * expr / max(float(expr.sum()), 1e-12)


class DumpFold:
    """fold522.Oct's loader for any dump708/dump522-format file (both write the same fields)."""
    kind = "oct"

    def __init__(self, fold_id: str, path: Path):
        self.fold_id, self.path = fold_id, Path(path)
        self.panel, self.held = fold_id.split("/")
        z = np.load(self.path)
        self.targets = [str(t) for t in z["targets"]]
        self.sources = [str(s) for s in z["sources"]]
        self.fin, self.fPP = z["fin"], z["fPP"]
        self.kv, self.carried = z["kv"], z["carried"]
        self.avail = z["avail"]
        self.V = z["V"]
        self.scored, self.pp_pos = z["scored"], z["panel_positions"]
        self.truth_full = np.asarray(z["truth"], np.float64)
        self.truth = self.truth_full[:, self.scored]
        self.betas = {c: z[f"beta_{c}"] for c in SLOTS}
        self.frames = {c: (z[f"b_{c}"], z[f"gate_{c}"]) for c in SLOTS}
        self.nonpanel = np.ones(self.fin.shape[1], bool)
        self.nonpanel[self.pp_pos] = False
        self.n = len(self.targets)
        self.P_DE = {c: [profile_from_effect(self.fPP[i], b, g, "minnorm") for i in range(self.n)]
                     for c, (b, g) in self.frames.items()}
        self.d = (self.truth * self.truth).sum(1)
        self.strength = np.linalg.norm(self.fin[:, self.nonpanel], axis=1)   # 546o / comp551: the unkeyed h2c1 row's norm
        self.cpm = {c: cpm_of(b) for c, (b, _) in self.frames.items()}
        self.base_fin = self.fin


class AbcFold:
    """A/B/C's own panel through 522's fold522.Abc (518's dumps; score518.Fold), with the same interface as DumpFold."""
    kind = "abc"

    def __init__(self, held: str, variant: str | None = None):
        sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
        import fold522
        if variant == "ABCp":
            # A/B/C's folds rebuilt with DAY-recipe held-out operators (dump518p.py; work/dumps/abcp/dump518-<held>.npz):
            # score518.Fold on that file, exactly fold522.Abc's route with the dump directory swapped. Its known answers
            # live at the dump stage (dump518p with 062's operators == 518's dumps); arms540's rows do not apply here.
            import plolo
            plolo.WORK = ROOT / "experiments/420-bulk-stack/work/plolo"
            import score518
            keep = score518.WORK
            score518.WORK = WORK / "dumps" / "abcp"
            try:
                Fo = score518.Fold(held, check=False)
            finally:
                score518.WORK = keep
            S = type("AbcpS", (), {})()
            S.F, S.targets, S.n, S.kv, S.carried = Fo, Fo.targets, Fo.n, Fo.kv1, Fo.carried
            S.sources = [s for s in fold522.SRC if s != held]
            S.ka = {"dump518p KA (062 ops == 518 dumps)": "see work/dumps/abcp/dump518p-*.json"}
        else:
            S = fold522.Abc(held, check=True)                     # its own KA: arms540's V1_0.5_3 / V1_0.5_6
        self.S, F = S, S.F
        self.variant = variant
        # panel stays "ABC" for ABCp (same panel, same rules; ka804 keys its A/B/C arm set on F.panel); fold_id says ABCp
        self.fold_id, self.panel, self.held = f"{variant or 'ABC'}/{held}", "ABC", held
        self.targets, self.n, self.sources = list(S.targets), S.n, S.sources
        self.kv, self.carried = S.kv, S.carried
        self.fin = F.fST
        self.base_fin = F.fST
        self.scored, self.nonpanel = F.scored, F.nonpanel
        self.pp_pos = np.flatnonzero(~F.nonpanel)
        self.truth, self.betas = F.truth, F.betas
        self.frames = F.fr
        self.P_DE = F.P_DE
        self.d = (F.truth * F.truth).sum(1)
        self.strength = np.linalg.norm(F.fST[:, F.nonpanel], axis=1)        # arms522c / sp2_522's A/B/C strength
        self.cpm = {c: cpm_of(b) for c, (b, _) in self.frames.items()}


def registry() -> dict:
    return json.loads(REGISTRY.read_text()) if REGISTRY.exists() else {"folds": {}}


def load(fold_id: str, check: bool = True):
    panel, held = fold_id.split("/")
    if panel == "ABC":
        F = AbcFold(held)
    elif panel == "ABCp":
        F = AbcFold(held, variant="ABCp")
    else:
        reg = registry()["folds"]
        if fold_id in reg:
            path = ROOT / reg[fold_id]["dump"]
        else:                                                  # the KA folds: dump801 --registry efp writes dump708-efp-<held>
            path = WORK / "dumps" / f"dump708-{panel}-{held}.npz"
        F = DumpFold(fold_id, path)
    F.ka = known_answers(F) if check else {"unchecked": True}
    return F


# ---------------------------------------------------------------------------------------------- arms
@dataclass
class Arm:
    w: np.ndarray
    q: object = 1.75
    fin: np.ndarray | None = None
    name: str = ""


def named(F, name: str) -> Arm:
    kv, car = F.kv, F.carried
    w3, w6 = wclip(kv, car, 0.5, 3.0), wclip(kv, car, 0.5, 6.0)
    wp = wclip(kv, car, 0.5, pct_hi(kv, car, P_PROD))
    w546 = wclip(kv, car, 0.5, pct_hi(kv, car, P_546))
    s_p2 = terc(F.strength, P2)
    arms = {"cap3": (w3, 1.75), "cap6": (w6, 1.75), "pct": (wp, 1.75),
            "pct+P2g": (wp, p2g_oct(wp, car)), "pct+S_P2": (wp, s_p2),
            "cap3+S_P2": (w3, s_p2), "cap3+P2g": (w3, p2g_oct(w3, car)),
            "cap6+S_P2": (w6, s_p2), "cap6+P2g": (w6, p2g_oct(w6, car)),
            "p96.6": (w546, 1.75), "p96.6+P2g": (w546, p2g_oct(w546, car)), "p96.6+S_P2": (w546, s_p2),
            "cap6+P2g(abc)": (w6, p2g_abc(w6)), "cap3+P2g(abc)": (w3, p2g_abc(w3))}
    arms["incumbent"] = arms["pct+P2g"]
    w, q = arms[name]
    return Arm(w=w, q=q, name=name)


def blend_rows(fin, ext, lam, match_norm=True):
    """Per-target blend of external rows into the composition rows: (1 - lam_i) fin_i + lam_i ext_i', where ext_i' is ext_i
    rescaled to |fin_i| when match_norm (panel columns zeroed, as in fin). lam = 1 swaps. A helper for candidate modules."""
    ext = np.asarray(ext, np.float64).copy()
    lam = np.broadcast_to(np.asarray(lam, float), (fin.shape[0],))[:, None]
    if match_norm:
        nf, ne = np.linalg.norm(fin, axis=1), np.linalg.norm(ext, axis=1)
        ext *= (nf / np.maximum(ne, 1e-300))[:, None]
    out = (1 - lam) * fin + lam * ext
    out[:, np.flatnonzero(np.all(fin == 0, axis=0))] = 0.0
    return out


# ---------------------------------------------------------------------------------------------- render + score
def render(F, arm: Arm):
    """sp2_522.score_q's loop: key_scale -> per-q-group row sharpening on each frame's gate -> render (minnorm) -> ballast
    realisation at KMAX 16 against the h1c1 carriers. Returns ({frame: realised rows}, K (frames x n))."""
    scalar = np.ndim(arm.q) == 0                               # fold522's scalar path (whole-bulk sharpen)
    qv = np.broadcast_to(np.asarray(arm.q, float), (F.n,)).copy()
    fin = F.fin if arm.fin is None else np.asarray(arm.fin, np.float64)
    if F.kind == "oct":
        key_scale, sharpen, pfe, realise_k, ts = L.key_scale, L.sharpen, profile_from_effect, L.realise_k, TS
    else:
        from lib420 import TS as TS4, key_scale, profile_from_effect as pfe, realise_k, sharpen
        ts = TS4
    bulk = key_scale(fin, np.asarray(arm.w, float))[0]
    rea, K = {}, np.zeros((len(SLOTS), F.n), int)
    for j, (c, (b, gate)) in enumerate(F.frames.items()):
        if scalar:
            X = sharpen(bulk, float(arm.q), gate, mode="row")
        else:
            X = np.empty_like(bulk)
            for q in np.unique(qv):
                sel = qv == q
                X[sel] = sharpen(bulk[sel], float(q), gate, mode="row")
        R = np.empty_like(X)
        for i in range(F.n):
            p_bk = pfe(X[i], b, gate, "minnorm")
            p_re, _, k, _, _ = realise_k(F.P_DE[c][i], p_bk, 16)
            R[i] = np.log1p(ts * p_re) - b
            K[j, i] = k
        rea[c] = R
    return rea, K


def _row_score(F, rea):
    if F.kind == "oct":
        return L.row_score(rea, F.scored, F.truth, F.betas, F.nonpanel)
    from comp420 import score_rows
    return score_rows(rea, F.scored, F.truth, F.betas, F.nonpanel)


def score_rows_of(F, rea, K):
    s, arr = _row_score(F, rea)
    s["K_mean"] = float(K.mean())
    out = dict(x=np.asarray(arr["x"]), a=np.asarray(arr["a"]), d=F.d.copy(), pds=np.asarray(arr["pds"]),
               pdsl=np.asarray(arr["pdsl"]), K=K)
    return s, out


def change_energy(F, rea, rea0):
    sc = np.zeros(F.nonpanel.size, bool)
    sc[F.scored] = True
    uns = F.nonpanel & ~sc
    keys = ["dE_scored", "dE_unscored", "dE_panel"] + [f"dE_cpm{k}" for k in range(len(CPM_EDGES) + 1)]
    E = {k: np.zeros((len(SLOTS), F.n)) for k in keys}
    for j, c in enumerate(F.frames):
        D2 = (rea[c] - rea0[c]) ** 2
        E["dE_scored"][j] = D2[:, sc].sum(1)
        E["dE_unscored"][j] = D2[:, uns].sum(1)
        E["dE_panel"][j] = D2[:, ~F.nonpanel].sum(1)
        bins = np.digitize(F.cpm[c], CPM_EDGES)
        for k in range(len(CPM_EDGES) + 1):
            E[f"dE_cpm{k}"][j] = D2[:, F.nonpanel & (bins == k)].sum(1)
    E["dE_total"] = E["dE_scored"] + E["dE_unscored"] + E["dE_panel"]
    return E


def cos_agg_arr(x, a, d, idx=None):
    """row_score's identity: mean over frames of sum x / sqrt(sum a * sum d) (boot522.cos_agg)."""
    idx = np.arange(d.size) if idx is None else idx
    return float(np.mean([x[j, idx].sum() / np.sqrt(a[j, idx].sum() * d[idx].sum()) for j in range(x.shape[0])]))


def price(F, arm: Arm, base: Arm | None = None):
    rea, K = render(F, arm)
    s, arr = score_rows_of(F, rea, K)
    if base is not None:
        rea0, K0 = render(F, base)
        s0, arr0 = score_rows_of(F, rea0, K0)
        arr.update(change_energy(F, rea, rea0))
        tot = float(arr["dE_total"].sum())
        s["vs_base"] = dict(
            base=base.name, r=s["cos_agg"] / s0["cos_agg"] - 1, dpds=s["pds"] - s0["pds"],
            dpds_lam05=s["pds_lam0.5"] - s0["pds_lam0.5"], dK_mean=s["K_mean"] - s0["K_mean"],
            dK_rel=s["K_mean"] / max(s0["K_mean"], 1e-300) - 1,
            share_unscored=float((arr["dE_unscored"].sum() + arr["dE_panel"].sum()) / tot) if tot > 0 else 0.0,
            share_cpm={f"cpm{k}": (float(arr[f"dE_cpm{k}"].sum() / tot) if tot > 0 else 0.0) for k in range(len(CPM_EDGES) + 1)})
        for k in ("x", "a", "pds", "pdsl", "K"):
            arr[f"base_{k}"] = arr0[k]
    return s, arr


# ---------------------------------------------------------------------------------------------- known answers
def _cmp(s, r):
    """L.cmp where the reference stores per-frame render energies, else the scalar fields it has."""
    if "AB" in r and isinstance(r["AB"], dict) and all(c in r["AB"] for c in SLOTS):
        return L.cmp(s, r)
    keys = [k for k in ("cos_agg", "pds") if k in r]
    return max(abs(s[k] - r[k]) for k in keys)


def known_answers(F) -> dict:
    t0 = time.time()
    ka = {}

    def run(name):
        arm = named(F, name)
        rea, K = render(F, arm)
        return score_rows_of(F, rea, K)[0]

    def check(label, s, ref, fields=None):
        if fields is None:
            d = _cmp(s, ref)
        else:
            d = max(abs(s[f] - ref[g]) for f, g in fields)
        ka[label] = d
        assert d < TOL, f"KA FAILED {F.fold_id} {label}: {d}"

    # KA-PT0: the per-target path at a constant q == fold522's scalar path (the same render with q a scalar)
    for name in ("cap3", "pct"):
        arm = named(F, name)
        r1, K1 = render(F, arm)
        r2, K2 = render(F, Arm(w=arm.w, q=np.full(F.n, 1.75)))
        d = max(float(np.max(np.abs(r1[c] - r2[c]))) for c in SLOTS) + float(np.max(np.abs(K1 - K2)))
        ka[f"PT0 {name}"] = d
        assert d == 0.0, f"KA FAILED {F.fold_id} PT0 {name}: {d}"
    S = {nm: run(nm) for nm in ("cap3", "cap6", "pct")}
    if F.panel in ("efp", "efr", "real"):
        ref = json.loads((W546 / f"comp546o-{F.panel}-{F.held}.json").read_text())["arms"]
        S.update({nm: run(nm) for nm in ("p96.6", "p96.6+P2g", "p96.6+S_P2")})
        for mine, theirs in (("cap3", "base q1.75"), ("cap6", "fixed 6"), ("p96.6", "p96.6"),
                             ("p96.6+P2g", "p96.6 + P2g"), ("p96.6+S_P2", "p96.6 + S_P2")):
            check(f"546o {theirs}", S[mine], ref[theirs]["realised"]["16"])
        p551 = W551 / F.panel / f"comp551-{F.held}.json"
        if p551.exists():
            r5 = json.loads(p551.read_text())["arms"]
            S.update({nm: run(nm) for nm in ("pct+P2g", "pct+S_P2", "cap3+P2g", "cap3+S_P2")})
            for mine, theirs in (("cap3", "h2c1|V1_5_3+sh1.75"), ("pct", "h2c1|V1pct+sh1.75"), ("pct+P2g", "h2c1|V1pct+ptP2g"),
                                 ("pct+S_P2", "h2c1|V1pct+ptS_P2"), ("cap3+P2g", "h2c1|V1_5_3+ptP2g"),
                                 ("cap3+S_P2", "h2c1|V1_5_3+ptS_P2")):
                check(f"551 {theirs}", S[mine], r5[theirs]["realised"]["16"])
    elif getattr(F, "variant", None) == "ABCp":
        ka["ABCp: stored-row KAs are at the dump stage (dump518p, 062 operators == 518's dumps)"] = 0.0
    elif F.panel == "ABC":
        ka["fold522.Abc (arms540 V1_0.5_3 / V1_0.5_6)"] = 0.0 if F.S.ka else None
        ref = json.loads((W540 / f"arms540c-{F.held}.json").read_text())["arms"]
        S.update({nm: run(nm) for nm in ("cap6+P2g(abc)", "cap3+P2g(abc)", "cap6+S_P2", "cap3+S_P2")})
        for mine, theirs in (("cap6+P2g(abc)", "V1_0.5_6 + P2g"), ("cap6+S_P2", "V1_0.5_6 + S_P2"),
                             ("cap3+P2g(abc)", "V1_0.5_3 + P2g"), ("cap3+S_P2", "V1_0.5_3 + S_P2"),
                             ("cap3", "V1_0.5_3"), ("cap6", "V1_0.5_6")):
            check(f"540c {theirs}", S[mine], ref[theirs]["realised"]["16"])
    p522 = W522 / f"arms522-{'abc' if F.panel == 'ABC' else F.panel}-{F.held}.json"
    if p522.exists() and getattr(F, "variant", None) != "ABCp":
        r22 = json.loads(p522.read_text())["arms"]
        for mine, theirs in (("cap3", "cap3"), ("cap6", "cap6"), ("pct", "pctPROD")):
            check(f"522 arms522 {theirs}", S[mine], r22[theirs],
                  fields=[("cos_agg", "cos_agg"), ("pds", "pds"), ("K_mean", "K_mean")])
    ka["_seconds"] = round(time.time() - t0, 1)
    return ka


def main(argv):
    if argv[0] == "ka":
        out = {}
        for fid in argv[1:]:
            t0 = time.time()
            F = load(fid, check=True)
            out[fid] = F.ka
            worst = max(v for k, v in F.ka.items() if not k.startswith("_") and isinstance(v, float))
            print(f"[{fid}] {len(F.ka) - 1} known answers PASS, worst |d| {worst:.1e} ({time.time() - t0:.0f}s, "
                  f"peak RSS {L.E062.peak_rss_gib():.2f} GiB)", flush=True)
            for k, v in F.ka.items():
                print(f"   {k:45s} {v}", flush=True)
        dst = WORK / "ka" / f"ka801-{'_'.join(f.replace('/', '-') for f in argv[1:])}.json"
        dst.parent.mkdir(parents=True, exist_ok=True)
        dst.write_text(json.dumps(out, indent=1, default=float))
    else:
        raise SystemExit(f"unknown command {argv[0]}")


if __name__ == "__main__":
    main(sys.argv[1:])
