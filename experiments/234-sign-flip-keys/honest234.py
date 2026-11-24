"""234 — the combination key, fitted AND dosed leave-one-fold-out, with `kolf2` and a bootstrap.

`combo234.py` fitted the model out-of-fold but still chose the DOSE on the cancer-fold mean,
i.e. on the evaluation folds themselves.  That is the same free choice 228 and 232 both had
to strip out before quoting a number, so it is stripped here:

  `fit`   the head sees only the OTHER two cancer folds' (key, correctness) pairs
  `dose`  `f*` is argmax of the OTHER two folds' curves; the held fold is never consulted
  `kolf2` fitted on all three cancer folds and applied to the pluripotent line, reported
          APART -- 232's best rule was negative there (-2.11% cosine) and `kolf2` is the
          only fold whose basal state is not a cancer line
  `boot`  a paired bootstrap over TARGETS (4,000 draws, 232 `retain232.py`'s instrument)

It also writes the chosen fold's key matrix to `combokey-<fold>.npy` so `arm234.py` can
realise it as an actual sign-flip arm and price it through all six members.

    uv run python honest234.py
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from common234 import CANCER, QW_NEEDED, RATE_PER_QW, SEED, WORK  # noqa: E402
from combo234 import FEATURES, fit_head, rank_features  # noqa: E402
from screen234 import FRACS, REGIONS, Fold, calibrated_ceiling, curve  # noqa: E402

ALL_FOLDS = CANCER + ("kolf2",)
SUBSAMPLE = 1_500_000
NBOOT = 4000
CONFIGS = [(k, t) for k in ("lda", "logistic") for t in ("pool", "all")]


def build_keys(folds, feats, rng):
    """LOFO combination keys for every fold.  Cancer folds are fitted on the other two
    cancer folds; `kolf2` on all three (it is never a training fold for a cancer fold)."""
    out = {}
    for kind, train_on in CONFIGS:
        heads = {}
        for h in ALL_FOLDS:
            train = [o for o in CANCER if o != h]
            Xs, ys, ws = [], [], []
            for o in train:
                g = folds[o]
                sel = g.pool.reshape(-1) if train_on == "pool" \
                    else np.ones(g.pool.size, dtype=bool)
                idx = np.flatnonzero(sel)
                if idx.size > SUBSAMPLE:
                    idx = rng.choice(idx, SUBSAMPLE, replace=False)
                Xs.append(feats[o][idx])
                ys.append(g.correct.reshape(-1)[idx].astype(np.float64))
                ww = (g.w.reshape(-1)[idx] if train_on == "pool"
                      else np.abs(g.truth).reshape(-1)[idx])
                ws.append(ww / max(ww.mean(), 1e-30))
            sc = fit_head(np.concatenate(Xs), np.concatenate(ys), np.concatenate(ws),
                          kind, rng)
            heads[h] = sc(feats[h]).reshape(folds[h].truth.shape).astype(np.float32)
        out[(kind, train_on)] = heads
        print(f"  keys built for {kind}/{train_on} ({len(heads)} folds)", flush=True)
    return out


def sel_at(f: Fold, key, region: str, frac: float):
    A = f.regions[region]
    k = np.where(A, key.astype(np.float64), np.inf)
    order = np.argsort(k, axis=1, kind="stable")
    m = A.sum(1)
    kk = np.minimum(np.rint(frac * m).astype(int), m)
    sel = np.zeros(f.truth.shape, dtype=bool)
    rows = np.arange(f.truth.shape[0])
    for i in rows:
        sel[i, order[i, :kk[i]]] = True
    return sel


def boot_dqw(f: Fold, sel, n=NBOOT, seed=SEED):
    """Paired bootstrap over TARGETS of the per-target Dq_w contribution."""
    per = (f.g * sel).sum(1)[f.live]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, per.size, size=(n, per.size))
    b = per[idx].mean(1)
    return dict(mean=float(per.mean()), se=float(b.std(ddof=1)),
                lo=float(np.percentile(b, 2.5)), hi=float(np.percentile(b, 97.5)),
                p_positive=float((b > 0).mean()))


def main() -> None:
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    folds = {h: Fold(h) for h in ALL_FOLDS}
    feats = {h: rank_features(folds[h]) for h in ALL_FOLDS}
    print(f"  features ({time.time()-t0:.0f}s)", flush=True)
    keys = build_keys(folds, feats, rng)
    del feats

    rows = []
    for (kind, train_on), kk in keys.items():
        # (per-region curves are built inside the loop below)

        for region in REGIONS:
            c = {h: curve(folds[h], kk[h], region, +1) for h in ALL_FOLDS}
            # HONEST dose: argmax on the OTHER two cancer folds only
            per_fold, doses = {}, {}
            for h in CANCER:
                other = np.mean([c[o] for o in CANCER if o != h], axis=0)
                j = int(np.argmax(other))
                doses[h] = FRACS[j]
                per_fold[h] = float(c[h][j])
            # kolf2 takes the dose the three cancer folds agree on
            jk = int(np.argmax(np.mean([c[o] for o in CANCER], axis=0)))
            doses["kolf2"] = FRACS[jk]
            per_fold["kolf2"] = float(c["kolf2"][jk])
            dq = float(np.mean([per_fold[h] for h in CANCER]))
            rows.append(dict(head=kind, train_on=train_on, region=region,
                             dqw_honest=dq, per_fold=per_fold, doses=doses,
                             overall={k: dq * v for k, v in RATE_PER_QW.items()},
                             n_pos=sum(1 for h in CANCER if per_fold[h] > 0),
                             curve={h: [float(v) for v in c[h]] for h in ALL_FOLDS}))

    hdr = (f"{'head':<10}{'train':<7}{'region':<11}{'f* (per fold)':<24}{'dqw':>11}"
           f"{'overall':>10}{'x short':>9}{'k+':>5}{'kolf2':>11}")
    print(f"\n{'='*len(hdr)}\nTHE COMBINATION KEY — fit AND dose leave-one-fold-out.  "
          f"Gate d q_w >= +{QW_NEEDED['k1.384']:.5f}\n{'='*len(hdr)}")
    print(hdr)
    print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: -x["dqw_honest"]):
        sh = QW_NEEDED["k1.384"] / r["dqw_honest"] if r["dqw_honest"] > 0 else float("inf")
        ds = "/".join(f"{r['doses'][h]:g}" for h in CANCER)
        print(f"{r['head']:<10}{r['train_on']:<7}{r['region']:<11}{ds:<24}"
              f"{r['dqw_honest']:>+11.5f}{r['overall']['k1.384']:>+10.5f}"
              f"{(f'{sh:.2f}x' if np.isfinite(sh) else '-'):>9}{r['n_pos']:>4d}/3"
              f"{r['per_fold']['kolf2']:>+11.5f}")

    best = max(rows, key=lambda x: x["dqw_honest"])
    print(f"\nBEST honest configuration: {best['head']}/{best['train_on']}/{best['region']}"
          f"  d q_w {best['dqw_honest']:+.5f}  "
          f"= {best['overall']['k1.384']:+.5f} .. {best['overall']['k1.000']:+.5f} overall")
    for h in CANCER:
        print(f"    {h:<15} f* {best['doses'][h]:<7g} d q_w {best['per_fold'][h]:+.5f}")
    print(f"    {'kolf2 (apart)':<15} f* {best['doses']['kolf2']:<7g} "
          f"d q_w {best['per_fold']['kolf2']:+.5f}")

    # ------------------------------------------------------------- bootstrap + save the key
    kk = keys[(best["head"], best["train_on"])]
    boots, mass = {}, {}
    for h in ALL_FOLDS:
        f = folds[h]
        sel = sel_at(f, kk[h], best["region"], best["doses"][h])
        boots[h] = boot_dqw(f, sel)
        mass[h] = float(np.mean((f.w * sel).sum(1)[f.live] / f.den[f.live]))
        np.save(WORK / f"combokey-{h}.npy", kk[h])
    print(f"\nPAIRED TARGET BOOTSTRAP ({NBOOT} draws) on the best arm")
    print(f"{'fold':<16}{'d q_w':>10}{'se':>10}{'95% CI':>22}{'P(>0)':>9}{'pool mass':>11}")
    for h in ALL_FOLDS:
        b = boots[h]
        print(f"{h:<16}{b['mean']:>+10.5f}{b['se']:>10.5f}"
              f"{f'[{b[chr(108)+chr(111)]:+.5f}, {b[chr(104)+chr(105)]:+.5f}]':>22}"
              f"{b['p_positive']:>9.3f}{mass[h]:>11.4f}")

    # ---------------------------------------------- Phase 2 ceiling for the combination key
    ceil = {}
    for nb in (20, 50, 200):
        inf = [calibrated_ceiling(folds[h], kk[h], nb)[0] for h in CANCER]
        ceil[f"infold_{nb}"] = float(np.mean(inf))
    print(f"\nPHASE 2 for the combination key (perfect in-fold calibration): "
          + "  ".join(f"@{n} {ceil[f'infold_{n}']:+.5f}" for n in (20, 50, 200)))

    (WORK / "honest234.json").write_text(json.dumps(
        dict(rows=rows, best=dict(head=best["head"], train_on=best["train_on"],
                                  region=best["region"], doses=best["doses"],
                                  dqw=best["dqw_honest"], per_fold=best["per_fold"]),
             bootstrap=boots, pool_mass=mass, ceiling_combo=ceil,
             qw_needed=QW_NEEDED, seconds=time.time() - t0), indent=1))
    print(f"\n-> {WORK/'honest234.json'}  ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
