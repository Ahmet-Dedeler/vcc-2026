"""810 -- where the strong tier's bulk misdirection lives: the incumbent's composition rows split into their parts, per target,
on every BULK-W (L1-L6 x 3 cancer folds) and BULK-ABC (A/B/C panel x k562 / hct116 / hek293t / kolf2) fold. EC2 only.

The composition (lib801 = lib438's arithmetic; the same calls dump801 / designs801 make) is
    fin = s * (W + f_cd4 C + f_h1 H + f_p P)        (panel columns zero; s = the post-add global preserve scale)
with W the whitened donor base (op.apply of the normrel'd source mean), C the deflated CD4 residual, H the H1 union residual,
P the CD4 partial rows. So the four PARTS base = sW, cd4 = s f_cd4 C, h1 = s f_h1 H, part = s f_p P sum to fin exactly and
    cos(fin, T) = sum_k <part_k, T> / (|fin| |T|)       (on the held line's scored genes)
splits the prediction's direction into additive contributions.

Per target it also records the source donors' rows (lib801.rows: the processed bank rows the base is built from, panel
columns zeroed), 707's raw unit mean U (data/donors delta rows of the fold's sources, sum of unit rows), the unwhitened
normrel'd donor mean `pre`, each source's and the held line's knockdown of the target itself (raw delta at the target's
own gene), 522's robust strength count and per-target noise for every bank (strength522.robust_strength), z_src and its
tertile (fold522.zsrc's rule over the fold's targets), and the Gram matrix of [base, cd4, h1, part, U, pre, d_src...] and
their inner products with the truth T on the scored genes (enough to evaluate any linear recombination's composition-level
cosine locally).

  compose <fold>   bind + compose + parts; KA (abort unless |d| < 1e-12): sum of parts == fin of the fold's registry dump
                   (L: dump708 `fin` = h2c1; ABC: 518's `fST` = h2c1.5) and the recomposed h1c1 == the dump's `fPP`; ABC only:
                   my robust strengths == strength522.npz's abc arrays (NaN positions included).
                   Writes $OUT810/<panel>__<held>-comp.npz/.json; the full part rows go to $PARTS810 (S3 later, not the Mac).
  realise <fold>   engine801 (its own KAs on load) prices the incumbent (pct+P2g) and stores x, a, d, pds, pdsl, K per frame and
                   target. ABC folds also run 706's top-tier beta 0.25 swap (pct cap, q 1.75, U blended, norm kept) through
                   dirswap706's committed code path and through the engine: they must agree (KA-706, |d| < 1e-12).
                   dirswap706-abc.json is reported beside them (the committed script does not reproduce it; dev/ka706dbg.py).

    OUT810=/work/out810 uv run python experiments/810-strong-direction/decomp810.py compose L3/orion_hct116
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
E801 = ROOT / "experiments/801-lab-folds"
sys.path.insert(0, str(E801))
sys.path.insert(0, str(ROOT / "experiments/522-october-bulk-key"))
import csv  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

OUT = Path(os.environ.get("OUT810", "/work/out810"))
PARTS = Path(os.environ.get("PARTS810", "/work/parts810"))
TOL = 1e-12
FOLDS = ["k562_gwps", "orion_hct116", "orion_hek293t"]
BANK_DIRS = {b: ROOT / "data/donors" / b for b in FOLDS}
GENES = [r["gene_name"] for r in csv.DictReader(open(ROOT / "data/challenge/gene_names.csv"))]
GPOS = {g: i for i, g in enumerate(GENES)}


def log(*a):
    print(f"[810 {time.strftime('%H:%M:%S')}]", *a, flush=True)


def sha_arr(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, np.float64).tobytes()).hexdigest()


def crank(v):
    """fold522.crank verbatim: centred percentile rank in [-0.5, 0.5] over the finite entries; 0 elsewhere."""
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    out = np.zeros(v.size)
    if ok.sum() > 1:
        r = np.argsort(np.argsort(v[ok]))
        out[ok] = r / (ok.sum() - 1) - 0.5
    return out


def tiers(v, n=3):
    """dirswap706.FD_tiers verbatim."""
    r = np.argsort(np.argsort(v, kind="stable"), kind="stable")
    return np.minimum((r * n) // len(v), n - 1)


def raw_rows(bank, targets):
    """dirswap706.donor_rows verbatim (data/donors delta on the official axis; zero rows where the bank lacks the target),
    plus the target's own knockdown (delta at the target's own gene; NaN where either is missing)."""
    D = BANK_DIRS[bank]
    dg = (D / "genes.txt").read_text().split()
    dp = (D / "perts.txt").read_text().split()
    delta = np.load(D / "delta.npy", mmap_mode="r")
    gpos = {g: i for i, g in enumerate(dg)}
    ppos = {p: i for i, p in enumerate(dp)}
    cols = np.array([gpos.get(g, -1) for g in GENES])
    ok = cols >= 0
    R = np.zeros((len(targets), len(GENES)))
    have = np.zeros(len(targets), bool)
    kd = np.full(len(targets), np.nan)
    for i, t in enumerate(targets):
        if t in ppos:
            row = np.asarray(delta[ppos[t]], float)
            R[i, ok] = row[cols[ok]]
            have[i] = bool(np.any(R[i]))
            if t in gpos:
                kd[i] = row[gpos[t]]
    return R, have, kd


def unit_mean(rows, nonpanel):
    U = np.zeros_like(rows[0])
    for R in rows:
        nr = np.linalg.norm(R, axis=1, keepdims=True)
        U += np.where(nr > 0, R / np.maximum(nr, 1e-300), 0.0)
    U[:, ~nonpanel] = 0.0
    return U


def fold_ids(fid):
    panel, held = fid.split("/")
    return panel, held, fid.replace("/", "__")


# ------------------------------------------------------------------------------------------------------------ compose
def cmd_compose(fid):
    import lib801 as L
    from lib801 import ALPHA, GAMMA, LAM
    from strength522 import robust_strength
    t0 = time.time()
    panel, held, stem = fold_ids(fid)
    abc = panel == "ABC"
    if abc:
        ctx = L.bind("abc")
        ref_p = ROOT / f"experiments/518-bulk-mse-gap/work/dump518-{held}.npz"
        ref = np.load(ref_p)
        ref_fin, ref_fpp = ref["fST"], ref["fPP"]
        m_cd4 = 1.5
    else:
        W = E801 / "work" / panel
        spec = dict(bundle=str((W / "bundle").resolve()), rules="def", cd4_root=str((W / "cd4/cd4_compact").resolve()),
                    full_op=str(Path("none").resolve()), folds=FOLDS, synthetic=True, meta=None,
                    fold_ops={h: (W / "ops" / f"held_{h}.npz").resolve() for h in FOLDS})
        ctx = L.bind(panel, spec=spec, tag=panel)
        ref_p = E801 / "work/dumps" / f"dump708-{panel}-{held}.npz"
        ref = np.load(ref_p)
        ref_fin, ref_fpp = ref["fin"], ref["fPP"]
        m_cd4 = 1.0
    m_h1 = 2.0
    fold = L.build_fold438(ctx, held, "held")
    pp_pos, wcols, scored = fold["panel_positions"], fold["wcols"], np.asarray(fold["scored"])
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    assert targets == [str(t) for t in ref["targets"]], "targets differ from the reference dump"
    dn = L.Donors438(targets, sources)
    cd4, partial, h1, aud = L.residuals(ctx, targets)
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fold["pre"])
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    log(f"{fid}: composed inputs ({time.time() - t0:.0f}s)")
    # ---- the composition's parts (compose / compose_arm's algebra, kept apart)
    if ctx.rules == "abc":
        f_cd4, _ = L.solve_fraction(whit, cd4, L.CD4_F * m_cd4)
    else:
        f_cd4, _ = L.solve_fraction(whit, cd4, L.CD4_G_DEF * m_cd4, sizing="carrying")
    comp = {"cd4": np.where(np.asarray(cd4.gene_votes) > 0, np.asarray(cd4.mean_effect), 0.0) * f_cd4}
    support = np.asarray(cd4.gene_votes) > 0
    frac = {"cd4": float(f_cd4)}
    g_part = L.lib418.PARTIAL_G if ctx.rules == "abc" else L.PARTIAL_G      # compose_arm's default vs compose's DEF value
    for tag, res, g in (("h1", h1, L.UNION_G * m_h1), ("part", partial, g_part)):
        v = np.asarray(res.gene_votes) > 0
        if not g or not v.any():
            comp[tag] = np.zeros_like(whit)
            frac[tag] = 0.0
            continue
        f, _ = L.solve_fraction_carrying(whit, res, g)
        comp[tag] = np.where(v, np.asarray(res.mean_effect), 0.0) * f
        support |= v
        frac[tag] = float(f)
    merged = comp["cd4"] + comp["h1"] + comp["part"]
    cand = np.asarray(whit, np.float64).copy()
    cand[support] += merged[support]
    s = (float(np.median(np.linalg.norm(whit[:, wcols], axis=1))) or 1.0) / \
        (float(np.median(np.linalg.norm(cand[:, wcols], axis=1))) or 1.0)
    parts = {"base": s * np.asarray(whit, np.float64), **{k: s * v for k, v in comp.items()}}
    for k in parts:
        parts[k][:, pp_pos] = 0.0
    fin = parts["base"] + parts["cd4"] + parts["h1"] + parts["part"]
    # ---- known answers
    ka = {}
    fC, _ = L.compose(ctx, whit, cd4, partial, h1, wcols, pp_pos, m_cd4=m_cd4, m_h1=m_h1)
    ka["parts sum == compose(final)"] = float(np.max(np.abs(fin - fC)))
    ka["compose(final) == reference fin"] = float(np.max(np.abs(fC - ref_fin)))
    ka["parts sum == reference fin"] = float(np.max(np.abs(fin - ref_fin)))
    fPP, _ = L.compose(ctx, whit, cd4, partial, h1, wcols, pp_pos, m_cd4=1.0, m_h1=1.0)
    ka["h1c1 == reference fPP"] = float(np.max(np.abs(fPP - ref_fpp)))
    ka["compose(final) == reference fin, bitwise"] = bool(np.array_equal(fC, ref_fin))
    for k, d in ka.items():
        if isinstance(d, float):
            log(f"   KA {k}: {d:.2e}")
            assert d < TOL, f"KA FAILED {fid} {k}: {d}"
    # ---- truth, donors, U, pre
    truth_full = np.asarray(fold["truth"], np.float64)
    if abc:
        T = np.asarray(ref["truth"], np.float64)            # 518's truth on 518's scored axis
        scored = np.asarray(ref["scored"])
        if T.shape[1] != scored.size:
            T = T[:, scored]
    else:
        T = truth_full[:, scored]
        assert np.allclose(T, np.asarray(ref["truth"], np.float64)[:, scored], atol=1e-6), "truth != dump truth"
    nonpanel = np.ones(len(GENES), bool)
    nonpanel[pp_pos] = False
    raw, haves, kds = {}, {}, {}
    for b in FOLDS:
        raw[b], haves[b], kds[b] = raw_rows(b, targets)
    U = unit_mean([raw[b] for b in sources], nonpanel)
    proc = {}
    for b in sources:
        r = np.asarray(dn.F[b], np.float64).copy()
        r[:, pp_pos] = 0.0
        proc[b] = r
    pre0 = np.asarray(pre, np.float64).copy()
    pre0[:, pp_pos] = 0.0
    vecs = [("base", parts["base"]), ("cd4", parts["cd4"]), ("h1", parts["h1"]), ("part", parts["part"]), ("U", U),
            ("pre", pre0)] + [(f"src{j}", proc[b]) for j, b in enumerate(sources)] \
        + [(f"raw{j}", raw[b]) for j, b in enumerate(sources)]
    names = [k for k, _ in vecs]
    Vs = np.stack([v[:, scored] for _, v in vecs], 1)            # n x K x scored
    G = np.einsum("nkg,nlg->nkl", Vs, Vs)
    bT = np.einsum("nkg,ng->nk", Vs, T)
    tt = np.einsum("ng,ng->n", T, T)
    Vn = np.stack([v[:, nonpanel] for _, v in vecs], 1)
    Gfull = np.einsum("nkg,nlg->nkl", Vn, Vn)                    # the same Gram on every non-panel gene (energy off-axis)
    del Vs, Vn
    # ---- strengths (522's robust strength on THIS panel), z_src over the fold's sources, tertiles
    pset = set(ctx.panel)
    st = {}
    for b in FOLDS:
        sv, cv, nv = robust_strength(BANK_DIRS[b], list(ctx.panel), pset)
        st[b] = (sv, cv, nv)
    if abc:
        z522 = np.load(ROOT / "experiments/522-october-bulk-key/work/strength522.npz")
        assert [str(t) for t in z522["abc__targets"]] == list(ctx.panel)
        for b in FOLDS:
            for j, suf in enumerate(("", "__count", "__noise")):
                a, r = st[b][j], np.asarray(z522[f"abc__{b}{suf}"], float)
                d = 0.0 if np.array_equal(np.isnan(a), np.isnan(r)) else np.inf
                d = max(d, float(np.nanmax(np.abs(a - r))))
                ka[f"strength522 abc {b}{suf}"] = d
                assert d < TOL, f"KA FAILED strength {b}{suf}: {d}"
        log("   KA strength522 (abc) PASS")
    pidx = {t: i for i, t in enumerate(ctx.panel)}
    ii = np.asarray([pidx[t] for t in targets])
    cols = []
    for b in sources:
        v = st[b][0][ii]
        r = crank(v)
        r[~np.isfinite(v)] = np.nan
        cols.append(r)
    with np.errstate(all="ignore"):
        zm = np.nanmean(np.vstack(cols), 0)
    z = np.where(np.isfinite(zm), zm, 0.0)
    out = dict(targets=np.asarray(targets), sources=np.asarray(sources), names=np.asarray(names), G=G, bT=bT, tt=tt,
               Gfull=Gfull, z=z, tier=tiers(z), kv=np.asarray(ref["kv"] if not abc else ref["kv1"], float),
               carried=np.asarray(ref["carried"], bool),
               **{f"kd_{b}": kds[b] for b in FOLDS}, **{f"have_{b}": haves[b] for b in FOLDS},
               **{f"str_{b}": st[b][0][ii] for b in FOLDS}, **{f"cnt_{b}": st[b][1][ii] for b in FOLDS},
               **{f"noise_{b}": st[b][2][ii] for b in FOLDS})
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"{stem}-comp.npz", **out)
    PARTS.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(PARTS / f"{stem}-parts.npz", targets=np.asarray(targets), scored=scored, pp_pos=pp_pos,
                        **{k: v for k, v in parts.items()}, U=U, pre=pre0)
    rec = dict(fold=fid, n=len(targets), sources=sources, rules=ctx.rules, m_cd4=m_cd4, m_h1=m_h1, fractions=frac, s=s,
               ka=ka, ref_dump=str(ref_p.relative_to(ROOT)), fin_sha256=sha_arr(fin), op=str(fold["op_path"]),
               deflation=aud, n_scored=int(scored.size), peak_rss_gib=L.E062.peak_rss_gib(), seconds=time.time() - t0)
    (OUT / f"{stem}-comp.json").write_text(json.dumps(rec, indent=1, default=lambda o: str(o) if isinstance(o, Path) else float(o)))
    log(f"{fid}: done in {time.time() - t0:.0f}s; peak RSS {rec['peak_rss_gib']:.2f} GiB; fractions {frac} s {s:.6f}")


# ------------------------------------------------------------------------------------------------------------ realise
def cmd_realise(fid):
    import engine801 as E
    t0 = time.time()
    panel, held, stem = fold_ids(fid)
    F = E.load(fid, check=True)
    inc = E.named(F, "incumbent")
    s, arr = E.price(F, inc)
    out = dict(targets=np.asarray(F.targets), w=np.asarray(inc.w, float), q=np.asarray(inc.q, float),
               **{k: np.asarray(arr[k]) for k in ("x", "a", "d", "pds", "pdsl", "K")})
    ka = {k: v for k, v in F.ka.items()}
    if panel == "ABC":                                        # KA-706: 706's top tier, beta 0.25, two implementations
        # (a) dirswap706's committed code path (fold522.Abc / score518.Fold.score with fST swapped; its own U, z, tiers);
        # (b) the engine with Arm(fin=...). They must agree to 1e-12. dirswap706-abc.json (written 1 min before the script's
        # last edit) is reported beside them, not asserted: the committed script does not reproduce it (dev/ka706dbg.py).
        sys.path.insert(0, str(ROOT / "experiments/706-strong-tier-budget"))
        import dirswap706 as DS
        import fold522 as FD
        ref = json.loads((ROOT / "experiments/706-strong-tier-budget/work/dirswap706-abc.json").read_text())[f"abc/{held}"]
        S = F.S
        z = FD.zsrc(S)
        U = np.zeros_like(S.F.fST)
        for R in (DS.donor_rows(src, S.targets)[0] for src in S.sources):
            nr = np.linalg.norm(R, axis=1, keepdims=True)
            U += np.where(nr > 0, R / np.maximum(nr, 1e-300), 0.0)
        U[:, ~S.F.nonpanel] = 0.0
        U810 = unit_mean([raw_rows(b, F.targets)[0] for b in F.sources], F.nonpanel)
        ka["KA-706 U (decomp810 == dirswap706)"] = float(np.abs(U - U810).max())
        have = np.linalg.norm(U, axis=1) > 0
        tz = np.full(F.n, -1)
        tz[have] = DS.FD_tiers(z[have])
        fin0 = S.F.fST.copy()
        nf, nu = np.linalg.norm(fin0, axis=1), np.linalg.norm(U, axis=1)
        Fn = fin0.copy()
        for i in np.flatnonzero((tz == 2) & (nf > 0) & (nu > 0)):
            m = 0.75 * fin0[i] / nf[i] + 0.25 * U[i] / nu[i]
            Fn[i] = nf[i] * m / np.linalg.norm(m)
        wp = FD.wclip(S.kv, S.carried, 0.5, FD.pct_hi(S.kv, S.carried, FD.P_PROD))
        s0a, a0a = S.score(wp)
        S.F.fST = Fn
        s1a, a1a = S.score(wp)
        S.F.fST = fin0
        live = np.flatnonzero(a0a["d"] > 0)
        dXa = 100 * (DS.X_of(np.asarray(a1a["x"]), np.asarray(a1a["a"]), a0a["d"], live)
                     / DS.X_of(np.asarray(a0a["x"]), np.asarray(a0a["a"]), a0a["d"], live) - 1)
        pct = E.named(F, "pct")
        s0, a0 = E.price(F, pct)
        s1, a1 = E.price(F, E.Arm(w=pct.w, q=1.75, fin=Fn))
        dX = 100 * (E.cos_agg_arr(a1["x"], a1["a"], a1["d"], live) / E.cos_agg_arr(a0["x"], a0["a"], a0["d"], live) - 1)
        ka["KA-706 w (engine pct == 706 wp)"] = float(np.abs(pct.w - wp).max())
        ka["KA-706 dX engine vs dirswap706 path"] = abs(dX - dXa)
        ka["KA-706 dpds engine vs dirswap706 path"] = abs((s1["pds"] - s0["pds"]) - (s1a["pds"] - s0a["pds"]))
        ka["706 json (reported): dX_pct, dpds"] = [ref["top_b0.25"]["dX_pct"], ref["top_b0.25"]["dpds"]]
        ka["706 committed path: dX_pct, dpds, n_sel"] = [dXa, s1a["pds"] - s0a["pds"], int((tz == 2).sum())]
        log(f"   KA-706 {fid}: engine dX {dX:+.6f}% == 706 path {dXa:+.6f}% (json {ref['top_b0.25']['dX_pct']:+.6f}); "
            f"dpds engine {s1['pds'] - s0['pds']:+.6f} path {s1a['pds'] - s0a['pds']:+.6f}")
        for k in ("KA-706 U (decomp810 == dirswap706)", "KA-706 w (engine pct == 706 wp)", "KA-706 dX engine vs dirswap706 path",
                  "KA-706 dpds engine vs dirswap706 path"):
            assert ka[k] < 1e-12, f"KA-706 FAILED {fid} {k}: {ka[k]}"
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"{stem}-real.npz", **out)
    (OUT / f"{stem}-real.json").write_text(json.dumps(dict(fold=fid, summary=s, ka=ka, seconds=time.time() - t0), indent=1,
                                                      default=float))
    log(f"{fid}: realised in {time.time() - t0:.0f}s; cos_agg {s['cos_agg']:.6f} pds {s['pds']:.6f}")


if __name__ == "__main__":
    {"compose": cmd_compose, "realise": cmd_realise}[sys.argv[1]](sys.argv[2])
