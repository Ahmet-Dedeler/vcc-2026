"""811 M811 arm module + EC2 fold driver (PROPOSAL-2.md). Frozen by sha256 in PROPOSAL-2.md before any M811 number exists.

M811 = the many-line consensus of base811 on the DAY bank list (no KOLF2 bank: its full rows are not on disk, so the D/E/F
day would not have it), optionally with PER-LINE CENTRING (factor C), blended into the composition at lambda.
  C off : U = base811.consensus(F, DAY groups minus the held group)            (bit-identical to base811's all_lX on L folds)
  C on  : each line row L_g(t) (normalised mean of its banks' unit rows, non-panel genes) has the line's mean over the
          PANEL's targets that the line carries subtracted, panel columns re-zeroed, re-normalised; U = mean over carrying
          lines of the centred rows (equal weight per line). On a fold the mean runs over the fold's full panel (the 300
          targets of its pert_counts.csv), as the D/E/F day runs it over its own 300. Lines carrying < MIN_CENTRE (10) panel
          targets are kept uncentred.
  blend : base811.blend (707's norm-preserving formula on the non-panel genes), all targets carrying U; then the incumbent's
          key, sharpening, render and ballast unchanged (engine801.price vs inc).
Arms: c0_l10, c0_l25, c1_l10, c1_l25 (C off/on x lambda 0.10/0.25), plus the KA arms inc and 707.
Known answers (a failed KA voids the fold; the driver exits non-zero):
  inc, 707      == 804's arrays (work/ref804, |d| <= 1e-12), as base811;
  c0_lXX        == base811's all_lXX arrays (work/run-base811) on every L fold (KOLF2 carries 0 L targets), |d| <= 1e-12;
  on 801's ABCp folds no 804 reference exists: the engine's fold KAs (E.load check=True) are the only KA there.
Day path: day_consensus(panel_targets, nonpanel, hold=identity_hold(identity_line)) with the same banks and rule.

    uv run python experiments/811-lolo-model/m811.py <fold-id> --arms c0_l10,c0_l25,c1_l10,c1_l25   # writes $M811_OUT
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402

import base811 as B  # noqa: E402  (imports engine801 and ka804 exactly as the baseline did)

E, K = B.E, B.K
OUT = Path(os.environ.get("M811_OUT", "/work/m811"))
BASE_RUN = HERE / "work" / "run-base811" / "bulk"
DAY_GROUPS = {g: b for g, b in B.GROUPS.items() if g != "KOLF2"}      # the D/E/F day's bank list (PROPOSAL-2 A)
MIN_CENTRE = 10     # a line carrying fewer panel targets than this is kept uncentred (its panel mean is not estimable)
ARMS = {"c0_l10": (False, 0.10), "c0_l25": (False, 0.25), "c1_l10": (True, 0.10), "c1_l25": (True, 0.25)}
# D/E/F identity rule (PROPOSAL-2 A): a context whose identity call matches a bank line holds that line's whole group out
IDENTITY_TO_GROUP = {"K562": "K562", "HCT116": "HCT116", "HCT-116": "HCT116", "HEK293T": "HEK293T", "HEK293": "HEK293T",
                     "JURKAT": "Jurkat", "A549": "A549", "HT29": "HT29", "HT-29": "HT29", "MCF7": "MCF7", "MCF-7": "MCF7",
                     "BXPC3": "BXPC3", "BXPC-3": "BXPC3", "HAP1": "HAP1", "DLD1": "DLD1", "DLD-1": "DLD1",
                     "RPE1": "RPE1", "HTERT-RPE1": "RPE1", "RPE-1": "RPE1"}


def identity_hold(identity_line: str | None) -> set[str]:
    """groups to hold out for a D/E/F context from its identity call (CCLE match); empty if no bank line matches."""
    if not identity_line:
        return set()
    g = IDENTITY_TO_GROUP.get(identity_line.upper().replace(" ", ""))
    return {g} if g else set()


def centred_consensus(targets: list[str], nonpanel: np.ndarray, groups: dict) -> tuple[np.ndarray, dict]:
    """factor C on: per-line centring over `targets` (the panel), then the equal-weight-per-line mean, rows for `targets`."""
    G = nonpanel.size
    acc = np.zeros((len(targets), G))
    nlines = np.zeros(len(targets), int)
    info = {}
    for line, banks in groups.items():
        S = np.zeros((len(targets), G))
        for b in banks:
            S += B.unit(B.bank_rows(b, list(targets)), nonpanel)
        Lr = B.unit(S, nonpanel)
        has = np.linalg.norm(Lr, axis=1) > 0
        if has.sum() >= MIN_CENTRE:
            mu = Lr[has].mean(0)
            sh = float(has.sum() * (mu @ mu) / (Lr[has] * Lr[has]).sum())
            Lc = np.where(has[:, None], Lr - mu[None, :], 0.0)
            Lr = B.unit(Lc, nonpanel)
            has = np.linalg.norm(Lr, axis=1) > 0
        else:
            sh = None
        info[line] = dict(n=int(has.sum()), shared_before=sh)
        acc += Lr
        nlines += has
    U = np.where(nlines[:, None] > 0, acc / np.maximum(nlines, 1)[:, None], 0.0)
    return U, info


def day_consensus(panel_targets: list[str], nonpanel: np.ndarray, centre: bool, hold: set[str] = frozenset()):
    """the D/E/F day's U for its own panel: DAY_GROUPS minus the identity-matched group(s)."""
    groups = {g: b for g, b in DAY_GROUPS.items() if g not in hold}
    if centre:
        return centred_consensus(panel_targets, nonpanel, groups)

    from types import SimpleNamespace                           # base811.consensus reads F.n, F.targets, F.nonpanel
    return B.consensus(SimpleNamespace(n=len(panel_targets), targets=list(panel_targets), nonpanel=nonpanel), groups)


def shared_share(U: np.ndarray) -> float | None:
    M = U[np.linalg.norm(U, axis=1) > 0]
    if len(M) < 3:
        return None
    mu = M.mean(0)
    return float(len(M) * (mu @ mu) / (M * M).sum())


def panel_targets(F) -> list[str]:
    """the fold's full panel: A/B/C's (also for 801's ABCp refit folds) or an L panel's pert_counts.csv."""
    if F.panel.startswith("ABC"):
        return K.panel_targets(type("P", (), {"panel": "ABC"})())
    return K.panel_targets(F)


def fold_U(F, centre: bool) -> tuple[np.ndarray, dict]:
    groups = {g: b for g, b in DAY_GROUPS.items() if g != B.HELD_GROUP.get(F.held)}
    if not centre:
        return B.consensus(F, groups)
    panel = panel_targets(F)
    Up, info = centred_consensus(panel, F.nonpanel, groups)
    pos = {t: i for i, t in enumerate(panel)}
    return Up[np.asarray([pos[t] for t in F.targets])], info


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("fold")
    ap.add_argument("--arms", required=True)
    a = ap.parse_args(argv)
    names = a.arms.split(",")
    assert all(n in ARMS for n in names), f"unknown arm in {names}"
    t0 = time.time()
    log = lambda *x: print(f"[{a.fold}]", *x, f"({time.time() - t0:.0f}s)", flush=True)
    F = E.load(a.fold, check=True)
    K.OUT = OUT
    # VERDICT-2 G1: engine801.AbcFold keeps F.panel == "ABC" on ABCp folds; ABCp is detected by fold_id / variant only
    abcp = F.fold_id.startswith("ABCp/") or getattr(F, "variant", None) == "ABCp"
    abc = F.panel == "ABC" and not abcp                         # 062-operator A/B/C folds: the board's ALLSW bulk is inc (804)
    known = abc or F.panel.startswith("L")                      # panels with 804 reference arrays; ABCp (801 refit) has none
    inc0 = E.named(F, "cap6+P2g(abc)" if abc else "pct+P2g")    # ABCp (day recipe) and L folds: the October incumbent
    inc = E.Arm(w=inc0.w, q=inc0.q, name="inc")
    fin = np.asarray(F.fin, np.float64)
    npn = F.nonpanel
    arms = {"inc": inc}
    if known:
        z = K.zsrc_fold(F)
        U7, have7 = K.donor_unit_mean(F)
        tier = np.full(F.n, -1)
        r = np.argsort(np.argsort(z[have7], kind="stable"), kind="stable")
        tier[np.flatnonzero(have7)] = np.minimum((r * 3) // int(have7.sum()), 2)
        arms["707"] = E.Arm(w=inc.w, q=inc.q, fin=B.blend(fin, U7, K.BETA, npn, rows=(tier == 2)), name="707")
    rec = dict(fold=a.fold, n=F.n, held=F.held, sources=F.sources, ka_engine=F.ka, arms={}, ka811={}, U={},
               day_groups={g: list(b) for g, b in DAY_GROUPS.items()}, script_sha256=B.sha(Path(__file__)),
               base811_sha256=B.sha(HERE / "base811.py"), banks_manifest_sha256=B.sha(B.BANKS / "manifest.json"))
    Ucache = {}
    for nm in names:
        centre, lam = ARMS[nm]
        if centre not in Ucache:
            Ucache[centre] = fold_U(F, centre)
            U, info = Ucache[centre]
            rec["U"][f"c{int(centre)}"] = dict(info=info, shared_share=shared_share(U),
                                               n_targets=int((np.linalg.norm(U, axis=1) > 0).sum()))
        U = Ucache[centre][0]
        arms[nm] = E.Arm(w=inc.w, q=inc.q, fin=B.blend(fin, U, lam, npn), name=nm)
    ref_dir = B.REF804 / a.fold.replace("/", "__")
    base_dir = BASE_RUN / a.fold.replace("/", "__")
    for nm, arm in arms.items():
        b = None if nm == "inc" else "inc"
        s, arr = E.price(F, arm, base=arms[b] if b else None)
        arr["targets"] = np.asarray(F.targets)
        rec["arms"][nm] = K.save(a.fold, nm, s, arr, b)
        keep = {k: v for k, v in arr.items() if not k.startswith("base_")}
        ref = None
        if nm in ("inc", "707") and known:
            ref = (ref_dir, nm)
        elif nm.startswith("c0_") and F.panel.startswith("L"):
            ref = (base_dir, "all_" + nm.split("_")[1])
        if ref:
            try:
                rec["ka811"][nm] = B.compare(a.fold, ref[1], keep, ref[0])
                log(f"KA {nm} == {ref[0].parent.name}/{ref[1]} (worst |d| {rec['ka811'][nm]['worst']:.1e})")
            except Exception as e:                              # a missing reference is a failed KA (810's G3)
                rec["ka811"][nm] = dict(FAILED=f"{type(e).__name__}: {e}")
                log(f"KA {nm} FAILED: {e}")
        if b:
            v = s["vs_base"]
            log(f"{nm:7s} r {100 * v['r']:+.3f}% dpds {v['dpds']:+.5f} dpdsl {v['dpds_lam05']:+.5f} dK {100 * v['dK_rel']:+.2f}% "
                f"off-axis {v['share_unscored']:.2f}")
    if not known:                                               # G1: the note goes in before the meta file is written
        rec["ka811"]["note"] = dict(info="no 804 reference arrays on this panel (ABCp); the engine's own fold KAs (F.ka) stand")
    rec["inc_arm"] = "cap6+P2g(abc)" if abc else "pct+P2g"
    rec["variant"] = "ABCp" if abcp else None
    (OUT / "meta").mkdir(parents=True, exist_ok=True)
    (OUT / "meta" / f"{a.fold.replace('/', '__')}.json").write_text(json.dumps(rec, indent=1, default=float))
    bad = [k for k, v in rec["ka811"].items() if "FAILED" in v]
    assert not bad, f"KA FAILED on {bad}: this fold's numbers are void"
    log("done")


if __name__ == "__main__":
    main()
