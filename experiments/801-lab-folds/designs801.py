"""801 -- composition rows for extra stack designs per fold (prospective row 710: H1 x2.5), written as
work/designs/<fold with / as __>.npz {targets, fin_<design>...} + work/designs/<fold>.json (sha256s, the KA).

  L <panel> <held>   L1-L6 folds, D/E/F rules: dump801.main(..., designs=["h2.5c1"]) -- dump708's own extra-design path,
                     the whole fold re-run with the fold's kept operator. KA (abort unless every array |d| < 1e-12): the re-run
                     dump equals the registry dump on EVERY array (fin = h2c1, fPP, kv, V, truth, betas, frames, ...); its bit
                     equality is recorded. Saves fin_h2c1 (the re-run's fin) and fin_h2.5c1.
  ABC <held>         A/B/C's own folds, A/B/C rules: lib801 (lib438) bound to "abc" (438's abc binding: data/challenge,
                     062's held operators / production operator for kolf2, A/B/C's CD4 banks, SVD-5 deflation), the fold
                     composed as dump522 composes (build_fold438, Donors438, residuals, pre_normrel, op.apply, Finals) and
                     FC.final("h2c1.5") / FC.final("h2.5c1.5"). KA (abort unless |d| < 1e-12): h2c1.5 == 518's dump fST (the
                     fin engine801's AbcFold uses), h1c1 == 518's fPP, and the targets equal.
Env DESIGNS801_STUB=1 swaps in discover801's identity operator (numbers meaningless; the data-closure run only; KAs are
recorded, not asserted). Rows are on the official axis with the panel's columns zero (Finals' contract).

    uv run python experiments/801-lab-folds/designs801.py L L3 orion_hct116     |    ... ABC kolf2
"""
from __future__ import annotations

import sys

sys.dont_write_bytecode = True
from pathlib import Path  # noqa: E402

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import guard801  # noqa: E402,F401
import hashlib  # noqa: E402
import json  # noqa: E402
import os  # noqa: E402
import time  # noqa: E402

import numpy as np  # noqa: E402

import lib801 as L  # noqa: E402

WORK = HERE / "work"
OUTD = WORK / "designs"
TOL = 1e-12
STUB = os.environ.get("DESIGNS801_STUB") == "1"
FOLDS = ["k562_gwps", "orion_hct116", "orion_hek293t"]


def sha256(p) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def stub():
    import discover801
    L.plolo.Operator = discover801.StubOperator


def check(ka, label, d):
    ka[label] = d
    print(f"   KA {label}: {d}", flush=True)
    if not STUB:
        assert d < TOL, f"KA FAILED {label}: {d}"


def run_L(panel: str, held: str) -> None:
    import dump801 as D
    if STUB:
        stub()
        D.TOL = float("inf")
    tmp = WORK / "designs_tmp"
    tmp.mkdir(parents=True, exist_ok=True)
    D.OUT = tmp
    D.FOLD_OPS = "held"
    W = WORK / panel
    spec = dict(bundle=str((W / "bundle").resolve()), rules="def", cd4_root=str((W / "cd4/cd4_compact").resolve()),
                full_op=str(Path("none").resolve()), folds=FOLDS, synthetic=True, meta=None,
                fold_ops={h: (W / "ops" / f"held_{h}.npz").resolve() for h in FOLDS})
    D.main(panel, held, spec, ["h2.5c1"], None)
    z = np.load(tmp / f"dump708-{panel}-{held}.npz")
    ref_p = WORK / "dumps" / f"dump708-{panel}-{held}.npz"
    ref = np.load(ref_p)
    ka, bitwise = {}, {}
    for k in ref.files:
        a, b = z[k], ref[k]
        if a.dtype.kind in "fc":
            d = float(np.nanmax(np.abs(a - b))) if a.size else 0.0
            if np.isnan(a).any() or np.isnan(b).any():
                d = max(d, 0.0 if np.array_equal(np.isnan(a), np.isnan(b)) else np.inf)
        else:
            d = 0.0 if np.array_equal(a, b) else np.inf
        bitwise[k] = bool(np.array_equal(a, b, equal_nan=a.dtype.kind in "fc"))
        check(ka, f"re-run == registry dump: {k}", d)
    fid = f"{panel}/{held}"
    save(fid, list(z["targets"]), {"h2c1": z["fin"], "h2.5c1": z["fin_h2.5c1"]}, ka,
         dict(rules="def", path="dump801 --designs h2.5c1 (dump708's extra-design path)", ref_dump=str(ref_p.relative_to(ROOT)),
              ref_dump_sha256=sha256(ref_p), bitwise_equal_arrays=bitwise, operator_json=json.loads(
                  (W / "ops" / f"held_{held}.json").read_text()) if (W / "ops" / f"held_{held}.json").exists() else None))


def run_ABC(held: str) -> None:
    if STUB:
        stub()
    from lib801 import ALPHA, GAMMA, LAM
    t0 = time.time()
    ctx = L.bind("abc")
    fold = L.build_fold438(ctx, held, "held")
    pp_pos, wcols = fold["panel_positions"], fold["wcols"]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    dn = L.Donors438(targets, sources)
    cd4, partial, h1, aud = L.residuals(ctx, targets)
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fold["pre"])
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    rows = {"h1c1": FC.final("h1c1"), "h2c1.5": FC.final("h2c1.5"), "h2.5c1.5": FC.final("h2.5c1.5")}
    ref_p = ROOT / f"experiments/518-bulk-mse-gap/work/dump518-{held}.npz"
    ref = np.load(ref_p)
    ka = {}
    check(ka, "targets == 518's", 0.0 if [str(t) for t in ref["targets"]] == targets else np.inf)
    check(ka, "h2c1.5 == 518 fST", float(np.max(np.abs(rows["h2c1.5"] - ref["fST"]))))
    check(ka, "h1c1 == 518 fPP", float(np.max(np.abs(rows["h1c1"] - ref["fPP"]))))
    ka["h2c1.5 == 518 fST bitwise"] = bool(np.array_equal(rows["h2c1.5"], ref["fST"]))
    ka["fractions"] = {c: FC.meta[c] for c in FC.meta}
    print(f"[ABC/{held}] composed in {time.time() - t0:.0f}s; peak RSS {L.E062.peak_rss_gib():.2f} GiB", flush=True)
    save(f"ABC/{held}", targets, {"h2c1.5": rows["h2c1.5"], "h2.5c1.5": rows["h2.5c1.5"]}, ka,
         dict(rules="abc", path="lib801.bind('abc') + Finals.final (438's abc binding)", ref_dump=str(ref_p.relative_to(ROOT)),
              ref_dump_sha256=sha256(ref_p), op=str(fold["op_path"]), deflation=aud))


def compose_abc(held: str, fold_ops: dict | None):
    """run_ABC's composition (438's abc binding); fold_ops=None keeps PANELS['abc'] (062's held operators), otherwise the
    given held-out operators replace them (tag 'abcp', its own caches)."""
    from lib801 import ALPHA, GAMMA, LAM
    if fold_ops is None:
        ctx = L.bind("abc")
    else:
        spec = dict(L.PANELS["abc"])
        spec["fold_ops"] = dict(spec["fold_ops"], **fold_ops)
        ctx = L.bind("abc", spec=spec, tag="abcp")
    fold = L.build_fold438(ctx, held, "held")
    pp_pos, wcols = fold["panel_positions"], fold["wcols"]
    targets, sources, op = list(fold["targets"]), list(fold["sources"]), fold["op"]
    dn = L.Donors438(targets, sources)
    cd4, partial, h1, aud = L.residuals(ctx, targets)
    pre = L.pre_normrel(dn, sources, "cap", -1.5, fallback_rows=fold["pre"])
    whit = op.apply(pre * ALPHA, wcols, GAMMA, LAM)
    FC = L.Finals(ctx, whit, cd4, partial, h1, wcols, pp_pos)
    rows = {"h1c1": FC.final("h1c1"), "h2c1.5": FC.final("h2c1.5"), "h2.5c1.5": FC.final("h2.5c1.5")}
    return targets, rows, op, fold, aud, {c: FC.meta[c] for c in FC.meta}


def row_sha(v) -> str:
    return hashlib.sha256(np.ascontiguousarray(v, np.float64).tobytes()).hexdigest()


def run_ABC062check(held: str) -> None:
    """With 062's operator, the same code must reproduce the stored ABC__<held> file's h2c1.5 row sha256 (bit for bit)."""
    targets, rows, op, fold, aud, meta = compose_abc(held, None)
    ref = json.loads((OUTD / f"ABC__{held}.json").read_text())
    rec = dict(held=held, op=str(fold["op_path"]), op_sha256=op.sha256, row_sha256={k: row_sha(v) for k, v in rows.items()},
               ref=f"work/designs/ABC__{held}.json", ref_row_sha256=ref["row_sha256"],
               match_h2c1_5=row_sha(rows["h2c1.5"]) == ref["row_sha256"]["h2c1.5"],
               match_h2_5c1_5=row_sha(rows["h2.5c1.5"]) == ref["row_sha256"]["h2.5c1.5"])
    (OUTD / "check").mkdir(parents=True, exist_ok=True)
    (OUTD / "check" / f"ABC062__{held}.json").write_text(json.dumps(rec, indent=1, default=str))
    print(f"[ABC062/{held}] h2c1.5 row sha == ABC__ file: {rec['match_h2c1_5']}; h2.5c1.5: {rec['match_h2_5c1_5']}", flush=True)
    assert rec["match_h2c1_5"], f"KA FAILED {held}: 062 rebuild's h2c1.5 row sha != ABC__{held}.json"


def run_ABCp(held: str) -> None:
    """A/B/C rules with the DAY-recipe held-out operator (work/ABCp/ops, the one dump518p's ABCp fold used)."""
    t0 = time.time()
    day = (WORK / "ABCp" / "ops" / f"held_{held}.npz").resolve()
    targets, rows, op, fold, aud, meta = compose_abc(held, {held: day})
    ref_p = WORK / "dumps" / "abcp" / f"dump518-{held}.npz"           # = engine801.load("ABCp/<held>").fin (fST), fPP
    ref = np.load(ref_p)
    chk = json.loads((OUTD / "check" / f"ABC062__{held}.json").read_text())
    ka = {}
    check(ka, "targets == ABCp fold's", 0.0 if [str(t) for t in ref["targets"]] == targets else np.inf)
    check(ka, "h2c1.5 == ABCp fold fin (fST)", float(np.max(np.abs(rows["h2c1.5"] - ref["fST"]))))
    check(ka, "h1c1 == ABCp fold fPP", float(np.max(np.abs(rows["h1c1"] - ref["fPP"]))))
    ka["h2c1.5 == ABCp fin bitwise"] = bool(np.array_equal(rows["h2c1.5"], ref["fST"]))
    check(ka, "062 operator, same code: h2c1.5 row sha == ABC__ file (0 = match)", 0.0 if chk["match_h2c1_5"] else np.inf)
    ka["062 check record"] = chk
    ka["fractions"] = meta
    assert Path(fold["op_path"]).resolve() == day, fold["op_path"]
    print(f"[ABCp/{held}] composed in {time.time() - t0:.0f}s; peak RSS {L.E062.peak_rss_gib():.2f} GiB", flush=True)
    save(f"ABCp/{held}", targets, {"h2c1.5": rows["h2c1.5"], "h2.5c1.5": rows["h2.5c1.5"]}, ka,
         dict(rules="abc", path="lib801.bind('abc', fold_ops -> the day operator) + Finals.final", operator=str(day.relative_to(ROOT)),
              operator_sha256=op.sha256, operator_s3=f"s3://<your-s3-bucket>/output/lab/ops801/ABCp/held_{held}.npz",
              ref_fold_dump=str(ref_p.relative_to(ROOT)), ref_fold_dump_sha256=sha256(ref_p), deflation=aud))


def save(fid, targets, rows, ka, extra):
    OUTD.mkdir(parents=True, exist_ok=True)
    stem = fid.replace("/", "__") + ("-STUB" if STUB else "")
    out = OUTD / f"{stem}.npz"
    np.savez_compressed(out, targets=np.asarray(targets), **{f"fin_{k}": np.asarray(v, np.float64) for k, v in rows.items()})
    rec = dict(fold=fid, stub=STUB, n=len(targets), designs=list(rows), npz=str(out.relative_to(ROOT)), npz_sha256=sha256(out),
               row_sha256={k: hashlib.sha256(np.ascontiguousarray(v, np.float64).tobytes()).hexdigest() for k, v in rows.items()},
               panel_cols_zero={k: bool(np.all(np.asarray(v)[:, ~nonpanel_of(fid)] == 0)) for k, v in rows.items()},
               ka=ka, peak_rss_gib=L.E062.peak_rss_gib(), **extra)
    (OUTD / f"{stem}.json").write_text(json.dumps(rec, indent=1, default=lambda o: str(o) if isinstance(o, Path) else float(o)))
    print(f"[{fid}] -> {out.relative_to(ROOT)} sha256 {rec['npz_sha256'][:16]}", flush=True)


def nonpanel_of(fid):
    panel = fid.split("/")[0]
    p = ROOT / "data/challenge/pert_counts.csv" if panel in ("ABC", "ABCp") else WORK / panel / "bundle/pert_counts.csv"
    tg = L.read_csv_col(p, "target_gene")
    m = np.ones(len(L.OFFICIAL), bool)
    m[[L.OFF_POS[g] for g in tg if g in L.OFF_POS]] = False
    return m


if __name__ == "__main__":
    if sys.argv[1] == "L":
        run_L(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "ABC":
        run_ABC(sys.argv[2])
    elif sys.argv[1] == "ABC062check":
        run_ABC062check(sys.argv[2])
    elif sys.argv[1] == "ABCp":
        run_ABCp(sys.argv[2])
    else:
        raise SystemExit(sys.argv[1])
