"""814 -- the frozen 10-22 reader of PROPOSAL-2's two panel predictions (PROPOSAL-2 sec. 3.2-3.3; VERDICT-2 item 4).
Nothing in this file changes after its commit. The lead runs it once on D/E/F's pert_counts.csv.

    uv run python experiments/814-panel-rule/predict814.py <pert_counts.csv> [--out <json>]
    uv run python experiments/814-panel-rule/predict814.py --ka          # KNOWN ANSWER on A/B/C (must pass first)

P-A  count of the panel's targets in the essential-screen LIBRARY union (Replogle RPE1 u Nadig HepG2 u Nadig Jurkat;
     work/flags/esslib.txt, 2,394 genes).               SUPPORT <= 3 | PARTIAL 4-7 | FAIL >= 8.
P-B  count of the panel's targets in prop10's frozen pool deciles 9-10 (labels p8, p9), denominator held at 300.
                                                         SUPPORT >= 130 | PARTIAL 96-129 | FAIL <= 95.
     Score s of a target: the 10-16 fit's own decision value where it exists (work/flags/propc.tsv, 4,825 genes), else
     s = intercept + sum_f coef_f * (x_f - mean_f) / sd_f on the frozen all-gene feature table
     (work/predict/feat6_allgenes.csv.gz, built by build814p.py --table, bit-equal to features814 on its genes).
     Decile = np.digitize(s, EDGES) (edges below, pinned by number: the 10-16 fit pool's deciles, build814p --constants).
     MISSING-FEATURE RULE (VERDICT-2 item 4): a target whose responsiveness (prod_logcnt_mean: no production bank
     carries it) or H1 control expression (h1_ctl_logcpm: not on the VCC-2025 H1 axis) is missing, or that is absent
     from the feature table, counts as OUTSIDE deciles 9-10. The denominator stays 300.
SAME-PANEL BRANCH (PROPOSAL-2 sec. 3.3): if the panel's target set equals A/B/C's, both predictions are VOID (the
     counts are still reported). Partial overlap with A/B/C is reported and changes nothing.
KNOWN ANSWER (--ka, abort on failure): on data/challenge/pert_counts.csv, P-A = 0 and P-B = 164 (RESULTS.md sec. 2/4).
Light: < 20 MB, about a second. Every input file is sha256-pinned; a different file aborts.
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
import argparse  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ESSLIB = (HERE / "work/flags/esslib.txt", "24fffbfee3d3fa3b038d1a7f3f665772ec171a2341ac35c03dafe6eb0a3d55e7")
PROPC = (HERE / "work/flags/propc.tsv", "b0b90368155f2b5ceb5cd5501d720f2b8a0b18da9ea8778cc23c0d399222b6cd")
TABLE = (HERE / "work/predict/feat6_allgenes.csv.gz", "4d6a3320e0dd17ecf01ea50672e18008dee5f582eec1527bd910de710418d16b")
ABC = (ROOT / "data/challenge/pert_counts.csv", None)              # sha recorded in every output; used for the KA and VOID

FEATS = ("prod_logcnt_mean", "weakkd", "in_vip", "pan_cov", "cd4_guide_anysig", "h1_ctl_logcpm")
MEAN = dict(prod_logcnt_mean=1.063575075843451, weakkd=0.23952738990332975, in_vip=0.9151450053705693,
            pan_cov=0.08915145005370569, cd4_guide_anysig=0.8150375939849624, h1_ctl_logcpm=1.5511795752471593)
SD = dict(prod_logcnt_mean=0.5971180038598142, weakkd=0.4268409055241318, in_vip=0.2786957302901241,
          pan_cov=0.28499283694423144, cd4_guide_anysig=0.3883087768126996, h1_ctl_logcpm=0.6134169873032962)
COEF = dict(prod_logcnt_mean=0.5004223374991579, weakkd=-0.2872195904830349, in_vip=0.12120485317295852,
            pan_cov=0.19766737109433136, cd4_guide_anysig=0.8701121730539231, h1_ctl_logcpm=0.11478870612349917)
INTERCEPT = -3.1332234701744532
EDGES = (-5.166920126721613, -3.9834400302754975, -3.4409317329866456, -3.1339519287497692, -2.9276880588060026,
         -2.74737192562082, -2.551923321221108, -2.3317644190722593, -1.9699499384713561)
BINS_A = dict(support=(None, 3), partial=(4, 7), fail=(8, None))     # inclusive ranges on the count
BINS_B = dict(support=(130, None), partial=(96, 129), fail=(None, 95))
DENOM = 300


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def pinned(spec) -> Path:
    p, want = spec
    got = sha256(p)
    assert want is None or got == want, f"ABORT: {p} sha256 {got} != pinned {want}"
    return p


def targets_of(p: Path) -> list[str]:
    df = pd.read_csv(p)
    col = "target_gene" if "target_gene" in df.columns else df.columns[0]
    t = [str(x) for x in df[col].dropna() if str(x) != "non-targeting"]
    return sorted(set(t))


def binned(n: int, bins: dict) -> str:
    for name, (lo, hi) in bins.items():
        if (lo is None or n >= lo) and (hi is None or n <= hi):
            return name.upper()
    raise AssertionError(n)


def read(panel_path: Path) -> dict:
    T = targets_of(panel_path)
    abc = targets_of(pinned(ABC))
    lib = {ln.strip() for ln in pinned(ESSLIB).read_text().splitlines() if ln.strip()}
    propc = {g: float(v) for g, v in (ln.split("\t") for ln in pinned(PROPC).read_text().splitlines() if ln.strip())}
    tab = pd.read_csv(pinned(TABLE), index_col=0, float_precision="round_trip")
    n_a = sum(t in lib for t in T)
    rows, n_b, missing, absent = {}, 0, [], []
    for t in T:
        if t not in tab.index:
            absent.append(t)
            rows[t] = dict(decile=None, rule="absent from the feature table: outside")
            continue
        x = tab.loc[t]
        if not (np.isfinite(x["prod_logcnt_mean"]) and np.isfinite(x["h1_ctl_logcpm"])):
            missing.append(t)
            rows[t] = dict(decile=None, rule="missing feature: outside")
            continue
        if t in propc:
            s, src = propc[t], "propc (10-16 fit)"
        else:
            s = INTERCEPT + sum(COEF[f] * (float(x[f]) - MEAN[f]) / SD[f] for f in FEATS)
            src = "frozen formula"
        d = int(np.digitize([s], EDGES)[0])                       # 0..9; deciles 9-10 = 8, 9
        n_b += d >= 8
        rows[t] = dict(score=s, decile=d + 1, source=src)
    same = set(T) == set(abc)
    hist = np.bincount([r["decile"] - 1 for r in rows.values() if r["decile"]], minlength=10).tolist()
    return dict(
        written=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), script_sha256=sha256(Path(__file__)),
        panel=str(panel_path), panel_sha256=sha256(panel_path), n_targets=len(T), denominator=DENOM,
        inputs=dict(esslib=ESSLIB[1], propc=PROPC[1], table=TABLE[1], abc_pert_counts=sha256(ABC[0])),
        overlap_with_abc=len(set(T) & set(abc)), same_panel_as_abc=same,
        P_A=dict(count=n_a, bins=BINS_A, verdict="VOID (same panel)" if same else binned(n_a, BINS_A)),
        P_B=dict(count=n_b, of=DENOM, deciles_hist=hist, n_missing_feature=len(missing), n_absent=len(absent),
                 missing=missing, absent=absent, bins=BINS_B, verdict="VOID (same panel)" if same else binned(n_b, BINS_B)),
        per_target=rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("panel", nargs="?")
    ap.add_argument("--ka", action="store_true")
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.ka:
        r = read(ABC[0])
        assert r["n_targets"] == 300 and r["P_A"]["count"] == 0 and r["P_B"]["count"] == 164, \
            f"KA FAIL: A/B/C P-A {r['P_A']['count']} (want 0), P-B {r['P_B']['count']} (want 164)"
        print(f"KA PASS: A/B/C P-A 0, P-B 164/300 (deciles {r['P_B']['deciles_hist']}); same-panel flag "
              f"{r['same_panel_as_abc']} -> verdicts {r['P_A']['verdict']} / {r['P_B']['verdict']}")
        return
    assert a.panel, "give a pert_counts.csv path (or --ka)"
    r = read(Path(a.panel))
    out = Path(a.out) if a.out else HERE / "work/predict" / f"predict814-{r['panel_sha256'][:16]}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(r, indent=1))
    print(f"panel {r['panel']} ({r['n_targets']} targets, overlap with A/B/C {r['overlap_with_abc']})")
    print(f"P-A  essential-library targets: {r['P_A']['count']}  -> {r['P_A']['verdict']}  (support <= 3, fail >= 8)")
    print(f"P-B  prop10 deciles 9-10: {r['P_B']['count']} / {DENOM}  -> {r['P_B']['verdict']}  (support >= 130, fail <= 95); "
          f"missing {r['P_B']['n_missing_feature']}, absent {r['P_B']['n_absent']}; deciles {r['P_B']['deciles_hist']}")
    print(f"-> {out}")


if __name__ == "__main__":
    main()
