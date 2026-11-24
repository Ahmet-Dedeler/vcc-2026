"""267 -- item 4: THE REPLICATE-VALIDATED ORACLE on one PLOLO fold.  MEASUREMENT ONLY.

Construction = 263's, imported (`lib263.stage`: the per-donor prepared rows, the fold operator resident, 261's
`cd4svd5_parts` re-solved against every arm's own base, 188's nmae identity NOT computed here).  What is new is only
WHICH truth the oracle sees and WHICH truth the arm is scored against:

  arm                      chosen on          scored on (every truth below, in one pass)
  mean (= cd4svd5)         --                 cache (KA vs 263), batchA, batchB, randA, randB, guide-mixes
  oracle:<T>:<gate>        truth T            all of the above  -> in-sample when T == the scoring truth, cross-
                                              validated when T is the other half (A->B and its mirror B->A)
  oracle3:<T>:all          truth T            the abstaining oracle (best of {donor rows, mean row})
  random:<gate>:s<k>       263's picks        control (KA vs 263 on the cache truth)
  union:wht:tau0:nm,       263's picks        263's best observable arbiters (build-time; no truth)
  basal:raw:tau0:nm

Every selection is norm-matched (263's `:nm`, the shape-only reading); the rank statistic is the scorer's
(resid218._per_target on the rendered rows, a23 = mass-fixed render and a15, three challenge frames averaged), with
the SAME half's truths of ALL targets as the distractor set.  Outputs work/sweep267-<fold>.{json,log} and
work/pertarget267-<fold>.npz (per-target pds / own-cosine per arm per truth) and work/picks267-<fold>.npz.

The `orion_hek293t` half-B truths are read from work/reserve/ (sha256 checked and logged) -- their ONE use.
    PYTHONDONTWRITEBYTECODE=1 uv run python sweep267.py orion_hct116
"""
from __future__ import annotations

import hashlib
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("~/Code/virtual-cell-challenge-2026")
HERE = ROOT / "experiments/267-replicate-validated-oracle"
WORK = HERE / "work"
E263 = ROOT / "experiments/263-donor-selection-arbiters"
sys.path.insert(0, str(E263))
import lib263 as Y  # noqa: E402   (imports lib257 -> lib246 -> plolo / resid218 / compose218)

# every cache read from OUR byte copies; nothing is ever written by the staging path (cache hits)
Y.plolo.WORK = WORK / "plolo"
Y.L.WORK = WORK
Y.X.WORK = WORK
Y.WORK = WORK
plolo, L = Y.plolo, Y.L
from resid218 import _per_target, deliver  # noqa: E402

TOL = 1e-9
LOG: list[str] = []


def P(s: str = "") -> None:
    print(s, flush=True)
    LOG.append(s)


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for blk in iter(lambda: f.read(8 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ truths
def load_truths(held: str, n: int, sc: np.ndarray, st_truth: np.ndarray) -> tuple[dict, dict, dict]:
    z = np.load(WORK / f"half/truth-{held}.npz", allow_pickle=True)
    assert np.array_equal(np.asarray(z["scored"]), sc)
    truths = {"cache": st_truth}
    avail = {"cache": np.ones(n, bool)}
    meta = {}
    for k in z.files:
        if k.endswith("|truth"):
            name = k.split("|")[0]
            truths[name] = np.asarray(z[k], dtype=np.float64)
            avail[name] = np.asarray(z[f"{name}|avail"], bool)
    if held == "orion_hek293t":
        rp = WORK / "reserve" / f"truthB-{held}.npz"
        sha = sha256_file(rp)
        want = (WORK / "reserve/SHA256SUMS.txt").read_text().split()[0]
        assert sha == want, "reserve sha mismatch"
        zr = np.load(rp, allow_pickle=True)
        for k in zr.files:
            if k.endswith("|truth"):
                name = k.split("|")[0]
                truths[name] = np.asarray(zr[k], dtype=np.float64)
                avail[name] = np.asarray(zr[f"{name}|avail"], bool)
        meta["reserve_access"] = dict(file=str(rp.relative_to(ROOT)), sha256=sha, keys=[k for k in zr.files if k.endswith("|truth")],
                                      purpose="item 4: the single scoring pass of sweep267.py (in-sample and cross-validated oracles, "
                                              "mean, random and arbiter controls scored on the half-B truths)")
        P(f"[{held}] RESERVE ACCESS {rp.relative_to(ROOT)} sha256 {sha[:16]}... keys {meta['reserve_access']['keys']}")
    # guide mixes: the two-construct targets carry the guide-half truth, every other target its batch-half truth
    gt = [str(t) for t in json.loads((WORK / f"half/truth-{held}.json").read_text())["guide_targets"]]
    for side in ("A", "B"):
        if f"guide{side}" in truths:
            T = truths[f"batch{side}"].copy()
            av = avail[f"batch{side}"].copy()
            g = avail[f"guide{side}"]
            T[g] = truths[f"guide{side}"][g]
            truths[f"guidemix{side}"] = T
            avail[f"guidemix{side}"] = av
            del truths[f"guide{side}"], avail[f"guide{side}"]
    meta["guide_targets"] = gt
    return truths, avail, meta


# ------------------------------------------------------------------------------------------------ scoring
class Scorer:
    """resid218.score_effects's loop with the renders shared across truths: the rendered rows are truth-independent."""

    def __init__(self, fold: dict, truths: dict, avail: dict):
        self.sc = np.asarray(fold["scored"], dtype=np.int64)
        self.frames = []
        for ctx in ("A", "B", "C"):
            b = plolo.challenge_bctrl(ctx)
            expr = np.expm1(b)
            cpm = 1e6 * expr / max(float(expr.sum()), 1e-12)
            gate = (cpm > plolo.MIN_CPM) & (b > 0)
            self.frames.append((ctx, b, gate))
        self.truths, self.avail = truths, avail

    def score(self, final: np.ndarray) -> dict:
        out = {}
        for tag, mass in (("a23", "minnorm"), ("a15", None)):
            acc = {name: ([], []) for name in self.truths}
            for ctx, b, gate in self.frames:
                dl = deliver(final, b, gate, mass)[:, self.sc]
                for name, T in self.truths.items():
                    av = self.avail[name]
                    p = np.full(final.shape[0], np.nan)
                    c = np.full(final.shape[0], np.nan)
                    pp, cc = _per_target(dl[av], T[av])
                    p[av], c[av] = pp, cc
                    acc[name][0].append(p)
                    acc[name][1].append(c)
            for name in self.truths:
                out[f"{tag}|{name}|pds"] = np.mean(acc[name][0], axis=0)
                out[f"{tag}|{name}|cos"] = np.mean(acc[name][1], axis=0)
        return out


def row_cos(A, B):
    return Y.row_cos(A, B)


# ------------------------------------------------------------------------------------------------ main
def main(held: str) -> None:
    t0 = time.time()
    st = Y.stage(held)
    n, sc, avail_d = st["n"], st["sc"], st["avail"]
    J = len(st["rows"])
    P(f"[{held}] staged: {n} targets, donors {st['donors']}, KA-0 pre max|d| {st['ka0']['pre_float32_max_abs_diff']:.1e}, "
      f"peak {Y.peak_gib():.2f} GiB, {time.time() - t0:.0f}s")
    truths, avail, tmeta = load_truths(held, n, sc, st["truth"])
    P(f"[{held}] truths: " + ", ".join(f"{k} ({int(v.sum())} avail)" for k, v in avail.items()))
    scorer = Scorer(st["fold"], truths, avail)
    cov = np.load(E263 / "work" / f"cov263-{held}.npz", allow_pickle=True)
    pk263 = np.load(E263 / "work" / f"picks263-{held}.npz")
    pt263 = np.load(E263 / "work" / f"pertarget263-{held}.npz")
    gates = {"all": np.asarray(cov["gate|all"], bool), "tau0": np.asarray(cov["gate|tau0"], bool)}
    assert np.array_equal(gates["all"], avail_d.sum(0) >= 2)
    ka = []
    failures = []

    def check(label, got, want, tol=TOL):
        d = float(np.nanmax(np.abs(np.asarray(got, float) - np.asarray(want, float))))
        ok = d <= tol
        ka.append(dict(label=label, max_abs_diff=d, ok=bool(ok)))
        if not ok:
            failures.append(label)
        P(f"  KA {'PASS' if ok else 'FAIL'} {label}: max|d| {d:.1e}")

    # oracle scores per truth: cosine of each donor's whitened row (alone) to the truth, scored columns
    S = {}
    Smean = {}
    for name, T in truths.items():
        s = np.stack([row_cos(st["wht_rows"][j][:, sc], T) for j in range(J)], 1)
        for j in range(J):
            s[~avail_d[j], j] = np.nan
        s[~avail[name]] = np.nan
        S[name] = s
        Smean[name] = np.where(avail[name], row_cos(st["whitened"][:, sc], T), np.nan)
    check("oracle scores on the cache truth == 263 cov `cos_truth_wht`", np.nan_to_num(S["cache"], nan=-9), np.nan_to_num(cov["cos_truth_wht"], nan=-9))
    check("base.truth on the cache truth == 263 cov `mean_cos_truth_wht`", Smean["cache"], cov["mean_cos_truth_wht"])

    rows = []
    per = {}
    picks = {}

    def run(arm: str, W, pick=None, extra=None):
        t1 = time.time()
        pre_sel = Y.combine(st, W, norm_match=True) if W is not None else st["pre"].copy()
        base = Y.base_from(st, pre_sel)
        parts, aud = Y.arm_parts(st, base)
        final = L.compose(dict(st, whitened=base), parts)
        e = scorer.score(final)
        for k, v in e.items():
            per[f"{arm}|{k}"] = v
        if pick is not None:
            picks[arm] = pick.astype(np.int8)
        row = dict(fold=held, arm=arm, n_selected=int((pick >= 0).sum()) if pick is not None else 0,
                   cd4_fraction=aud["cd4_fraction"], union_fraction=aud["union_fraction"],
                   median_row_norm_final=float(np.median(np.linalg.norm(final, axis=1))), seconds=time.time() - t1)
        for name in truths:
            row[f"pds|{name}"] = float(np.nanmean(e[f"a23|{name}|pds"]))
            row[f"cos|{name}"] = float(np.nanmean(e[f"a23|{name}|cos"]))
            row[f"pds15|{name}"] = float(np.nanmean(e[f"a15|{name}|pds"]))
        if extra:
            row.update(extra)
        rows.append(row)
        P(f"  {arm:28s} sel {row['n_selected']:3d}  " + "  ".join(f"{name} {row[f'pds|{name}']:.4f}" for name in truths) +
          f"  ({row['seconds']:.1f}s)")
        return e

    # ---- the mean (= cd4svd5): KA vs 263 on the cache truth
    e = run("mean", None, extra=dict(control="the incumbent base, equal-weight mean of the donor rows"))
    check("mean|cache per-target pds == 263 `cd4svd5|pds`", e["a23|cache|pds"], pt263["cd4svd5|pds"])
    check("mean|cache per-target cos == 263 `cd4svd5|cos`", e["a23|cache|cos"], pt263["cd4svd5|cos"])

    # ---- the oracles, chosen on each truth, scored on every truth
    for name in truths:
        for g in ("all", "tau0"):
            p = Y.pick_from_scores(S[name], avail_d, gates[g])
            e = run(f"oracle:{name}:{g}", Y.one_hot(st, p), pick=p, extra=dict(chosen_on=name, gate=g, oracle="whitened row cosine"))
            if name == "cache":
                check(f"oracle:cache:{g} picks == 263 `oracle_wht:{g}:nm`", p, pk263[f"oracle_wht:{g}:nm"])
                check(f"oracle:cache:{g}|cache per-target pds == 263 `oracle_wht:{g}:nm|pds`", e["a23|cache|pds"], pt263[f"oracle_wht:{g}:nm|pds"])
        # abstaining oracle (gate all): the best of {donor rows, the mean row}
        p3 = np.full(n, -1, dtype=np.int64)
        for t in np.flatnonzero(gates["all"] & avail[name]):
            s = np.where(avail_d[:, t], S[name][t], -np.inf)
            j = int(np.argmax(s))
            if np.isfinite(s[j]) and s[j] > Smean[name][t]:
                p3[t] = j
        e = run(f"oracle3:{name}:all", Y.one_hot(st, p3), pick=p3, extra=dict(chosen_on=name, gate="all", oracle="best of {donor rows, mean row}"))
        if name == "cache":
            check("oracle3:cache:all picks == 263 `oracle3_wht:all:nm`", p3, pk263["oracle3_wht:all:nm"])
            check("oracle3:cache:all|cache per-target pds == 263 `oracle3_wht:all:nm|pds`", e["a23|cache|pds"], pt263["oracle3_wht:all:nm|pds"])

    # ---- controls: 263's random picks (three seeds) and its best observable arbiters
    for g in ("all", "tau0"):
        for k in range(3):
            p = np.asarray(pk263[f"random:{g}:s{k}:nm"], dtype=np.int64)
            e = run(f"random:{g}:s{k}", Y.one_hot(st, p), pick=p, extra=dict(control="random pick (263's seeds)", gate=g, seed=k))
            check(f"random:{g}:s{k}|cache per-target pds == 263 `random:{g}:s{k}:nm|pds`", e["a23|cache|pds"], pt263[f"random:{g}:s{k}:nm|pds"])
    for arb in ("union:wht:tau0:nm", "basal:raw:tau0:nm", "cd4d:wht:tau0:nm"):
        p = np.asarray(pk263[arb], dtype=np.int64)
        e = run(f"arb:{arb}", Y.one_hot(st, p), pick=p, extra=dict(control="263 observable arbiter", arbiter=arb))
        check(f"arb:{arb}|cache per-target pds == 263 `{arb}|pds`", e["a23|cache|pds"], pt263[f"{arb}|pds"])

    # ---- oracle agreement between halves (the picks themselves)
    agree = {}
    for g in ("all", "tau0"):
        for a, b in (("batchA", "batchB"), ("randA", "randB"), ("guidemixA", "guidemixB")):
            if f"oracle:{a}:{g}" in picks and f"oracle:{b}:{g}" in picks:
                pa, pb = picks[f"oracle:{a}:{g}"], picks[f"oracle:{b}:{g}"]
                m = (pa >= 0) & (pb >= 0)
                agree[f"{a}~{b}:{g}"] = dict(n=int(m.sum()), agree=float((pa[m] == pb[m]).mean()) if m.any() else float("nan"))
            if f"oracle:{a}:{g}" in picks and f"oracle:cache:{g}" in picks:
                pa, pc = picks[f"oracle:{a}:{g}"], picks[f"oracle:cache:{g}"]
                m = (pa >= 0) & (pc >= 0)
                agree[f"{a}~cache:{g}"] = dict(n=int(m.sum()), agree=float((pa[m] == pc[m]).mean()) if m.any() else float("nan"))
    P(f"[{held}] pick agreement: " + ", ".join(f"{k} {v['agree'] * 100:.0f}% (n {v['n']})" for k, v in agree.items()))

    meta = dict(fold=held, n_targets=n, targets=st["targets"], donors=st["donors"], J=J, truths=list(truths),
                avail={k: int(v.sum()) for k, v in avail.items()}, gates={k: int(v.sum()) for k, v in gates.items()},
                ka=ka, ka_pass=sum(1 for k in ka if k["ok"]), ka_total=len(ka), failures=failures, agree=agree,
                frac_cd4_stack=st["frac_cd4"], frac_x=st["frac_x"], seconds=time.time() - t0, peak_rss_gib=Y.peak_gib(),
                fold_cache=st["fold"]["audit"].get("cache"), **tmeta)
    (WORK / f"sweep267-{held}.json").write_text(json.dumps(dict(meta=meta, rows=rows), indent=1, default=str))
    np.savez_compressed(WORK / f"pertarget267-{held}.npz", **per)
    np.savez_compressed(WORK / f"picks267-{held}.npz", **picks)
    (WORK / f"sweep267-{held}.log").write_text("\n".join(LOG) + "\n")
    P(f"[{held}] done: {len(rows)} arms, KA {meta['ka_pass']}/{meta['ka_total']}, failures {failures}, "
      f"{(time.time() - t0) / 60:.1f} min, peak {Y.peak_gib():.2f} GiB")


if __name__ == "__main__":
    main(sys.argv[1])
