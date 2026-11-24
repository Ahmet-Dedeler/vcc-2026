"""801 -- a COPY of experiments/708-def-panel-repricing/dump708.py (git blob 175666add66f50eac782b2143888a583dfe3b974; itself
dump522.py's per-fold arithmetic on any bundle) through lib801 (lib438's copy, caches under 801/work/lib) instead of lib438.
Changes: this paragraph; lib801 for lib438; OUT = 801/work/dumps; and --ka522 <panel>, which runs dump522.py's own known
answers (git blob 42ff0d2f6049a95346b5a0e1d354ca0d6f3c3f3f, lines `checks = [...]`, verbatim) against the stored rows of a
registry panel (efp / efr / real): 438's `h2c1|V1_5_3+sh1.5`, 546o's `base q1.75`, `fixed 6`, `p96.6` and 540o's `V1_0.5_6`.
One fix: the json record serialises Path values in the spec (dump708 raises TypeError there whenever --fold-ops held
builds a spec with Path operators, i.e. on its own CLI path). Everything else is dump708 verbatim. dump708's docstring follows.

708 Part B: dump522.py's per-fold arithmetic on ANY bundle (the real D/E/F panel on 10-22, or the rehearsal bundle),
bound the way dayof438.py binds a day-of panel (a spec dict: bundle, D/E/F rules, the day's CD4 banks and operators).

Per held-out cancer fold it saves what 708's score708.py needs (work/dump708-<tag>-<fold>.npz): the same fields as
dump522 (targets, sources, fin = h2c1 final, fPP = h1c1 carriers, kv = the fold's LODO-V1 key, carried, per-donor LODO
terms, truth, betas, frames) plus fin_<design> for each extra --designs entry (708 AMENDMENT 1's Part A confirmation).

KNOWN ANSWERS (abort unless |d| < 1e-12):
  - always: the per-donor LODO terms' mean == lodo_key (1e-12 relative);
  - --ka438 <dir>: 438's stored "h2c1|V1_5_3+sh1.5" row (realised and design) for this fold, from <dir>/comp438-<fold>-held.json
    (the rehearsal on 438's abcdef inputs; the day has no reference).

    uv run python experiments/708-def-panel-repricing/dump708.py --tag T --bundle B --op OP --cd4-root CD4 \
        --fold-ops held --held-ops DIR [--synthetic] [--designs h2.5c1.5] [--ka438 DIR] <fold>        # ~6 GB peak
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))                                   # 801: lib801 (lib438 copy) + guard801
import argparse  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib801 as L  # noqa: E402  (imports guard801 first)
from lib801 import ALPHA, GAMMA, LAM, SLOTS, E062  # noqa: E402
import arms414  # noqa: E402
from plolo import TS, profile_from_effect  # noqa: E402

W438 = ROOT / "experiments/438-october-stack-repricing/work"
W540 = ROOT / "experiments/540-lodo-clip-scan/work"
W546 = ROOT / "experiments/546-october-cap-pertarget/work"
WORK = HERE / "work"
OUT = Path(__file__).resolve().parent / "work" / "dumps"
OUT.mkdir(parents=True, exist_ok=True)
TOL = 1e-12


def main(panel: str, held: str, spec: dict | None, designs: list, ka438: Path | None, ka522: bool = False) -> None:
    t0 = time.time()
    log = lambda *a: print(f"[{panel}/{held}]", *a, f"({time.time() - t0:.0f}s)", flush=True)
    ctx = L.bind(panel, spec=spec, tag=panel)
    fold = L.build_fold438(ctx, held, FOLD_OPS)
    scored, pp_pos, wcols = fold["scored"], fold["panel_positions"], fold["wcols"]
    truth_full = fold["truth"]
    truth = truth_full[:, scored]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    nonpanel = ctx.nonpanel
    dn = L.Donors438(targets, sources)
    cd4, partial, h1, _ = L.residuals(ctx, targets)
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fold["pre"])
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    key410, carried410, _ = L.shared_energy_key(targets, sources, op, wcols, cols=None)
    W3 = L.weights(key410, carried410, 0.3, 3.0)[0]
    d_rows = {s: L.donor_rows(s, targets)[0] for s in sources}
    A_ = {s: dn.A[s] for s in sources}
    carried_all = sum(A_[s].astype(int) for s in sources) >= 2
    bases = L.lodo_bases(op, dn, sources, "cap", -1.5, pre, wcols)
    fPP, pPP = FC.final("h1c1"), FC.parts("h1c1")
    kv1 = L.lodo_key(fPP, pPP, bases, d_rows, A_, carried_all, sources, pp_pos)
    kvv = np.nan_to_num(kv1, nan=0.0)
    # ---- per-donor LODO terms (lodo_key's arithmetic, cols=None), and their parts
    pp = np.einsum("ij,ij->i", fPP, fPP)
    ns = len(sources)
    V = np.full((ns, len(targets)), np.nan)
    dnorm, pnorm, cosPd = (np.full((ns, len(targets)), np.nan) for _ in range(3))
    for a, j in enumerate(sources):
        Pj = pPP["s"] * (bases[j] + pPP["radd"])
        Pj[:, pp_pos] = 0.0
        num = np.einsum("ij,ij->i", Pj, d_rows[j])
        V[a] = np.where(A_[j], num / pp, np.nan)
        dj = d_rows[j].copy()
        dj[:, pp_pos] = 0.0
        nd, npj = np.linalg.norm(dj, axis=1), np.linalg.norm(Pj, axis=1)
        dnorm[a] = np.where(A_[j], nd, np.nan)
        pnorm[a] = np.where(A_[j], npj, np.nan)
        cosPd[a] = np.where(A_[j], num / np.maximum(nd * npj, 1e-300), np.nan)
    V[:, ~carried_all] = np.nan
    with np.errstate(all="ignore"):
        kmean = np.nan_to_num(np.nanmean(V, axis=0), nan=0.0)
    dmax = float(np.max(np.abs(kmean[carried_all] - kvv[carried_all]) / np.maximum(np.abs(kvv[carried_all]), 1e-300)))
    log(f"KA per-donor mean == lodo_key: max rel |d| {dmax:.1e}")
    assert dmax < 1e-12, "KA per-donor terms"
    d0, d1 = d_rows[sources[0]].copy(), d_rows[sources[1]].copy()
    d0[:, pp_pos] = 0.0
    d1[:, pp_pos] = 0.0
    both = A_[sources[0]] & A_[sources[1]]
    cos_dd = np.where(both, np.einsum("ij,ij->i", d0, d1)
                      / np.maximum(np.linalg.norm(d0, axis=1) * np.linalg.norm(d1, axis=1), 1e-300), np.nan)
    del bases, d_rows, d0, d1
    P_DE = L.carrier_profiles(fPP)
    bulk_k3 = L.key_scale(fPP, W3)[0]
    dk3 = L.realised_rows(bulk_k3, P_DE, q=None, kmaxes=(), zero_cols=None)
    bet = L.betas(ctx, held, nonpanel, A_pp_k3={c: L.energy(dk3[c]["design"], nonpanel)["A"] for c in SLOTS})
    fin = FC.final("h2c1")
    frames = {c: (np.asarray(b), np.asarray(g)) for c, (b, g) in arms414.frames().items()}
    holder = {}

    def rows_q(q):
        out = {}
        for c, (b, gate) in frames.items():
            F = L.sharpen(holder["bulk"], float(q), gate, mode="row")
            des, rea = np.empty_like(F), np.empty_like(F)
            for i in range(F.shape[0]):
                p_bk = profile_from_effect(F[i], b, gate, "minnorm")
                des[i] = np.log1p(TS * p_bk) - b
                p_re, _, _, _, _ = L.realise_k(P_DE[c][i], p_bk, 16)
                rea[i] = np.log1p(TS * p_re) - b
            out[c] = dict(design=des, realised=rea)
        return out

    def rec(w, q):
        holder["bulk"] = L.key_scale(fin, w)[0]
        dr = rows_q(q)
        sd, _ = L.row_score({c: dr[c]["design"] for c in SLOTS}, scored, truth, bet, nonpanel)
        sr, _ = L.row_score({c: dr[c]["realised"] for c in SLOTS}, scored, truth, bet, nonpanel)
        return sd, sr

    ka = {}
    if ka438 is not None:
        ref438 = json.loads((Path(ka438) / f"comp438-{held}-held.json").read_text())["arms"]["h2c1|V1_5_3+sh1.5"]
        sd, sr = rec(L.weights(kvv, carried_all, 0.5, 3.0)[0], 1.5)
        dR = L.cmp(sr, ref438["realised"]["16"])
        dD = L.cmp(sd, ref438["design"])
        ka["438 h2c1|V1_5_3+sh1.5"] = dict(realised=dR, design=dD, cos_agg=sr["cos_agg"], pds=sr["pds"])
        log(f"KA 438 h2c1|V1_5_3+sh1.5 cos {sr['cos_agg']:.9f} pds {sr['pds']:.6f}: realised |d| {dR:.1e} design |d| {dD:.1e}")
        assert dR < TOL and dD < TOL, f"KA FAILED {panel}/{held}: {dR} {dD}"
    if ka522:                                   # 801: dump522.py's own known answers, verbatim (registry panels only)
        ref438 = json.loads((W438 / panel / f"comp438-{held}-held.json").read_text())["arms"]["h2c1|V1_5_3+sh1.5"]
        ref546 = json.loads((W546 / f"comp546o-{panel}-{held}.json").read_text())["arms"]
        ref540 = json.loads((W540 / f"comp540o-{panel}-{held}.json").read_text())["arms"]
        kn = kvv[carried_all] / float(np.median(kvv[carried_all]))
        hi546 = float(np.percentile(kn, 96.647))
        checks = [("438 h2c1|V1_5_3+sh1.5", L.weights(kvv, carried_all, 0.5, 3.0)[0], 1.5, ref438, True),
                  ("546 base q1.75", L.weights(kvv, carried_all, 0.5, 3.0)[0], 1.75, ref546["base q1.75"], False),
                  ("546 fixed 6", L.weights(kvv, carried_all, 0.5, 6.0)[0], 1.75, ref546["fixed 6"], False),
                  ("540o V1_0.5_6", L.weights(kvv, carried_all, 0.5, 6.0)[0], 1.75, ref540["V1_0.5_6"], False),
                  ("546 p96.6", L.weights(kvv, carried_all, 0.5, hi546)[0], 1.75, ref546["p96.6"], False)]
        for name, w, q, ref, with_design in checks:
            sd, sr = rec(w, q)
            dR = L.cmp(sr, ref["realised"]["16"])
            dD = L.cmp(sd, ref["design"]) if with_design else 0.0
            ka[name] = dict(realised=dR, design=dD, cos_agg=sr["cos_agg"], pds=sr["pds"])
            log(f"KA {name:24s} cos {sr['cos_agg']:.9f} pds {sr['pds']:.6f}: realised |d| {dR:.1e} design |d| {dD:.1e}")
            assert dR < TOL and dD < TOL, f"KA FAILED {panel}/{held} {name}: {dR} {dD}"
    extra = {}
    for dsg in designs:
        extra[f"fin_{dsg}"] = FC.final(dsg)
        log(f"extra design {dsg}: |fin| mean {float(np.linalg.norm(extra[f'fin_{dsg}'], axis=1).mean()):.4f}")
    out = OUT / f"dump708-{panel}-{held}.npz"
    np.savez_compressed(out, targets=np.asarray(targets), sources=np.asarray(sources), fin=fin, fPP=fPP, kv=kvv,
                        carried=carried_all, avail=np.stack([np.asarray(A_[s], bool) for s in sources]), V=V,
                        dnorm=dnorm, pnorm=pnorm, cosPd=cosPd, cos_dd=cos_dd, Fnorm=np.sqrt(pp), scored=np.asarray(scored),
                        panel_positions=np.asarray(pp_pos), truth=np.asarray(truth_full, np.float32),
                        **{f"beta_{c}": bet[c] for c in SLOTS}, **{f"b_{c}": frames[c][0] for c in SLOTS},
                        **{f"gate_{c}": frames[c][1] for c in SLOTS}, **extra)
    rec_ = dict(panel=panel, fold=held, n=len(targets), sources=sources, ka=ka, designs=designs, spec=spec,
                n_carried=int(carried_all.sum()), peak_rss_gib=E062.peak_rss_gib(), seconds=time.time() - t0)
    (OUT / f"dump708-{panel}-{held}.json").write_text(json.dumps(rec_, indent=1, default=lambda o: str(o) if isinstance(o, Path) else float(o)))  # 801: dump708 crashes here on a spec with Path fold_ops
    log(f"dumped {out.name}; peak RSS {rec_['peak_rss_gib']:.2f} GiB")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--bundle", type=Path, default=Path("."))
    ap.add_argument("--op", type=Path, default=Path("none"))
    ap.add_argument("--cd4-root", type=Path, default=Path("."))
    ap.add_argument("--fold-ops", default="held", choices=("held", "prod"))
    ap.add_argument("--held-ops", type=Path, default=None)
    ap.add_argument("--synthetic", action="store_true")
    ap.add_argument("--designs", default="")
    ap.add_argument("--ka438", type=Path, default=None)
    ap.add_argument("--registry", action="store_true", help="801: bind the lib's registry panel <tag> (efp / efr / real) instead of a spec")
    ap.add_argument("--ka522", action="store_true", help="801: run dump522's known answers (registry panels)")
    ap.add_argument("fold")
    a = ap.parse_args()
    FOLDS = ["k562_gwps", "orion_hct116", "orion_hek293t"]
    spec_ = dict(bundle=str(a.bundle.resolve()), rules="def", cd4_root=str(a.cd4_root.resolve()), full_op=str(a.op.resolve()),
                 folds=FOLDS, synthetic=a.synthetic, meta=None, fold_ops={})
    if a.fold_ops == "held" and not a.registry:
        spec_["fold_ops"] = {h: Path((a.held_ops / f"held_{h}.npz").resolve()) for h in FOLDS}
    FOLD_OPS = a.fold_ops
    if a.registry:                              # 801: the registry spec (paths in lib801.PANELS)
        spec_ = None
    main(a.tag, a.fold, spec_, [d for d in a.designs.split(",") if d], a.ka438, a.ka522)
