"""801 -- a COPY of experiments/518-bulk-mse-gap/dump518.py (git blob d50961133c2f7bb7572be8da0253a7617b16ed87): A/B/C's PLOLO fold
under the BOARD's A/B/C rules (lib420's path: CD4 x1.5 all-rows sizing, SVD-5 deflation re-fitted on the rows, H1 x2 stack),
with the held-out OPERATOR as an argument. Changes, nothing else:
  * guard801 first (writes outside 801 raise); plolo.WORK -> 801/work/lib/abcp_plolo/<ops> (fold caches rebuilt there);
  * --ops 062 (062's cloud whiteners, the operators 518's dumps used) | day (DAY-recipe re-fits, build_union_whitener.py
    --basis-prune essential, 12k per-donor cap, on data/challenge: work/ABCp/ops/held_<held>.npz), set by overriding
    plolo.FOLDS[held]["op"] in-process;
  * with --ops 062: dump518's own KAs (KA-PP, V1 == 418's V1__PP, arms540 V1_0.5_3 / V1_0.5_6), PLUS every array == 518's
    stored dump518-<held>.npz (|d| < 1e-12); with --ops day: those reference KAs do not apply (another operator), KA-PP still;
  * output work/dumps/abcp/dump518-<held>.npz (--ops day; the file engine801's ABCp/<held> reads) or
    work/dumps/abc062/dump518-<held>.npz (--ops 062), plus dump518p-<ops>-<held>.json (KA, sha256, the fold's pct-cap hi).
dump518's docstring follows.

518 (PLOLO) — one heavy pass per fold: dump everything a per-target KEY arm needs, so every later arm runs without the
whitening operator (<~2 GB instead of ~5.5 GB).

Exactly arms540.py's setup (lib420, 420's path): the fold's unkeyed PP and stack compositions (fPP, fST; panel zeroed), the
LODO-V1 key vector kv1 and its carried mask (KA: weights(kv1, 0.3, 3) == 418's saved `V1__PP`), the scored columns, the
truth, panel positions. A key arm is then  key_scale(fST, w) -> sharpen q 1.75 per frame (row mode on the frame's gate)
-> render (minnorm) -> realise at KMAX 16 against PP's unkeyed carriers (P_DE from fPP).
KNOWN ANSWER (abort unless |d| < 1e-12, run here before saving): recomputing V1_0.5_3 and V1_0.5_6 from the dumped arrays
through lib420's own functions reproduces arms540-<fold>.json's realised cos_agg and pds.

    uv run python experiments/518-bulk-mse-gap/dump518.py <fold>        (~5.5 GB peak, like arms540; one fold at a time)
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "experiments/420-bulk-stack"))
sys.path.insert(0, str(ROOT / "experiments/518-bulk-mse-gap"))   # 801: rows518
sys.path.insert(0, str(HERE))                                        # 801: guard801
import guard801  # noqa: E402,F401
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib420  # noqa: E402
from lib420 import (COMPS, E062, FOLDS, PROD_OPS, Donors, Operator, betas_for, build_fold,  # noqa: E402
                    cd4_frozen_and_partial, compose_arm, compose_pp, deflate, key_scale, pre_normrel, source_residual,
                    union_rows, weights)
import lodo415  # noqa: E402
from comp420 import score_rows  # noqa: E402
import plolo  # noqa: E402

OPS = sys.argv[2] if len(sys.argv) > 2 else "062"                 # 801: which held-out operators
assert OPS in ("062", "day"), OPS
plolo.WORK = HERE / "work" / "lib" / "abcp_plolo" / OPS            # 801: fold caches rebuilt inside 801
plolo.WORK.mkdir(parents=True, exist_ok=True)
DAY_OPS = HERE / "work" / "ABCp" / "ops"
REF518 = ROOT / "experiments/518-bulk-mse-gap/work"
P_PROD = 96.64697522769474


def sha256(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def pct_hi(kv, carried, p=P_PROD):
    kn = kv[carried] / float(np.median(kv[carried]))
    return float(np.percentile(kn, p))

W418 = ROOT / "experiments/418-bulk-composition/work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
WORK = HERE / "work" / "dumps" / ("abcp" if OPS == "day" else "abc062")   # 801
WORK.mkdir(parents=True, exist_ok=True)
TOL = 1e-12


def run(held: str) -> None:
    t0 = time.time()
    if OPS == "day":                                              # 801: the operator swap, in-process
        # op_kind "production": plolo.build_fold's branch for a self-described build_union_whitener fit (sha256 pinned by its
        # json, axis from the sources' banks, checked again by Operator._check_axis); the "held" branch wants 062's
        # metadata format (axis/gamma/wcols_sha256), which a day-recipe fit's json does not carry
        plolo.FOLDS[held] = dict(plolo.FOLDS[held], op=DAY_OPS / f"held_{held}.npz", op_kind="production")
    log = lambda *a: print(f"[{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    build_fold(held)                                              # 801: builds the cache in 801 if absent
    fold = build_fold(held)                                       # 801: the cache hit (float32-rounded pre/truth), as 518 saw it
    scored = np.asarray(fold["scored"])
    pp_pos = np.asarray(fold["panel_positions"])
    wcols = np.asarray(fold["wcols"])
    truth_full = np.asarray(fold["truth"], dtype=np.float64)
    pre_cached = np.asarray(fold["pre"], dtype=np.float64)
    targets = list(fold["targets"])
    sources = list(FOLDS[held]["sources"])
    op = Operator(FOLDS[held]["op"] if OPS == "day" else  # 801: the day operator itself, never the production file
                  (PROD_OPS / "whitener_union_pruned.npz" if FOLDS[held]["op_kind"] == "production" else FOLDS[held]["op"]))
    nonpanel = np.ones(truth_full.shape[1], bool)
    nonpanel[pp_pos] = False
    betas = betas_for(held, nonpanel)
    h1 = union_rows([source_residual(n, targets, "gene_norm") for n in lodo415.UNION])
    dn = Donors(targets, sources)
    frozen, partial = cd4_frozen_and_partial(targets)
    cd4, _ = deflate(frozen, 5)
    pre228 = dn.pre(sources, False)
    pre228[:, pp_pos] = 0.0
    nocov = sum(dn.A[s].astype(int) for s in sources) == 0
    pre228[nocov] = pre_cached[nocov]
    fb = pre228.copy()
    pre = pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fb)
    f415, _, _, radd415, _ = compose_pp(op, pre, cd4, partial, h1, wcols, pp_pos)
    whit = op.apply(pre * lib420.ALPHA, wcols, lib420.GAMMA, lib420.LAM)
    fPP, pPP = compose_arm(whit, cd4, partial, h1, wcols, pp_pos)
    assert np.array_equal(fPP, f415) and np.array_equal(pPP["radd"], radd415), "KA-PP compose"
    mc, gh = COMPS["stack"]
    fST, pST = compose_arm(whit, cd4, partial, h1, wcols, pp_pos, m_cd4=mc, g_h1=gh)
    import lib418 as _l418
    d_rows = {s_: _l418.donor_rows(s_, targets)[0] for s_ in sources}
    A_ = {s_: dn.A[s_] for s_ in sources}
    carried_all = sum(A_[s_].astype(int) for s_ in sources) >= 2
    bases = _l418.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    kv1 = np.nan_to_num(_l418.lodo_key(fPP, pPP, bases, d_rows, A_, carried_all, sources, pp_pos), nan=0.0)
    if OPS == "062":                                              # 801: a reference KA of 062's operators only; dump518
        _w = weights(kv1, carried_all, 0.3, 3.0)[0]                   # asserted it bit for bit on the Mac; across Mac -> EC2
        _r = np.load(W418 / f"arrays418-{held}.npz")["V1__PP"]        # the key differs in the last bits, so |d| < 1e-12
        V1_d, V1_bit = float(np.max(np.abs(_w - _r))), bool(np.array_equal(_w, _r))
        assert V1_d < TOL, f"KA V1: {V1_d}"
    else:
        V1_d, V1_bit = None, None
    # per-source per-target features the key does not use directly (row norms on the scored axis, whitened-base agreement)
    src_norm = np.stack([np.linalg.norm(d_rows[s_][:, scored], axis=1) for s_ in sources])
    src_avail = np.stack([np.asarray(A_[s_], bool) for s_ in sources])
    whit_norm = np.linalg.norm(whit[:, scored], axis=1)
    radd_norm = np.linalg.norm(pST["radd"][:, scored], axis=1)
    del bases, whit, f415, pre228, fb
    log("KA-PP and KA V1 PASS; composing the known-answer rows from the dump")

    # ---- KA: V1_0.5_3 and V1_0.5_6 from the dumped arrays == arms540's realised rows
    import rows518
    frames = rows518.challenge_frames()
    P_DE = rows518.de_profiles(fPP, frames)
    ref = json.loads((W540 / f"arms540-{held}.json").read_text())["arms"]
    ka = {}
    for lo, hi in ((0.5, 3.0), (0.5, 6.0)):
        w = weights(kv1, carried_all, lo, hi)[0]
        rea = rows518.realised_rows(key_scale(fST, w)[0], frames, P_DE, q=1.75)
        sr, _ = score_rows(rea, scored, truth_full[:, scored], betas, nonpanel)
        tag = f"V1_{lo:g}_{hi:g}"
        want = ref[tag]["realised"]["16"]
        d = max(abs(sr[k] - want[k]) for k in ("cos_agg", "pds", "pds_lam0.5"))
        ka[tag] = d
        log(f"{'KA' if OPS == '062' else 'info (another operator)'} {tag}: cos_agg {sr['cos_agg']:.12f} vs arms540 {want['cos_agg']:.12f}; max |d| {d:.1e}")
        if OPS == "062":
            assert d < TOL, f"KA FAILED {held} {tag}: {d}"
    np.savez_compressed(WORK / f"dump518-{held}.npz", targets=np.asarray(targets), sources=np.asarray(sources),
                        fPP=fPP, fST=fST, kv1=kv1, carried=carried_all, scored=scored, panel_positions=pp_pos,
                        truth=truth_full.astype(np.float32), src_norm=src_norm, src_avail=src_avail, whit_norm=whit_norm,
                        radd_norm=radd_norm, **{f"beta_{c}": betas[c] for c in "ABC"})
    (WORK / f"dump518-{held}.json").write_text(json.dumps(dict(fold=held, n=len(targets), sources=sources, ka=ka,
                                                               peak_rss_gib=E062.peak_rss_gib()), indent=1, default=float))
    # ---- 801: the record (and, for 062's operators, every array == 518's stored dump)
    out = WORK / f"dump518-{held}.npz"
    rec = dict(ka_V1_vs_418=V1_d, ka_V1_bitwise=V1_bit, fold=held, ops=OPS, op=str(op.path.relative_to(ROOT)), op_sha256=op.sha256, dump=str(out.relative_to(ROOT)),
               dump_sha256=sha256(out), ka=ka, pct_hi=pct_hi(kv1, carried_all), n_carried=int(carried_all.sum()),
               peak_rss_gib=E062.peak_rss_gib(), seconds=time.time() - t0)
    z, ref = np.load(out), np.load(REF518 / f"dump518-{held}.npz")
    rec["ref_dump"] = str((REF518 / f"dump518-{held}.npz").relative_to(ROOT))
    rec["ref_pct_hi"] = pct_hi(ref["kv1"], ref["carried"])
    cmpd, bitw = {}, {}
    for k in ref.files:
        a, b = z[k], ref[k]
        if a.shape != b.shape:
            cmpd[k] = float("inf")
            bitw[k] = False
            continue
        cmpd[k] = float(np.max(np.abs(a.astype(np.float64) - b.astype(np.float64)))) if a.dtype.kind in "fc" and a.size else \
            (0.0 if np.array_equal(a, b) else float("inf"))
        bitw[k] = bool(np.array_equal(a, b))
    rec["vs_518_max_abs"], rec["vs_518_bitwise"] = cmpd, bitw
    if OPS == "062":
        bad = {k: v for k, v in cmpd.items() if not v < TOL}
        assert not bad, f"KA FAILED {held}: rebuild with 062's operators != 518's dump on {bad}"
        log(f"KA every array == 518's dump: worst |d| {max(cmpd.values()):.1e}")
    (WORK / f"dump518p-{OPS}-{held}.json").write_text(json.dumps(rec, indent=1, default=float))
    log(f"pct-cap hi {rec['pct_hi']:.4f} (518's dump: {rec['ref_pct_hi']:.4f}); op {op.sha256[:16]}")
    log(f"dumped; peak RSS {E062.peak_rss_gib():.2f} GiB")


if __name__ == "__main__":
    run(sys.argv[1])                                              # 801: dump518p.py <held> [062|day]
