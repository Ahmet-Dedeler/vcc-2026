"""200 — the four DE members per ablation on h1-test, MATCHED-SEED paired, in overall units.

CONVENTIONS, all inherited verbatim from 169 §4 so that experiment's numbers are a
known-answer row for this one:

  pool        the reference's OWN BH<0.05 set -- a SIGNIFICANCE-selected pool, which 188 §2
              measured to be reproduced 0.7219 by `|lfc|*sqrt(cpm)` and only 0.3796 by
              `|lfc|`.  A magnitude pool inverts the sign of `nmae` deltas, so it is not used.
              `Surface(n_real_target=N)` is NOT used either: it flattens `n_real` and deletes
              the small-`n_real` tercile that carries 98% of our reach (AGENTS.md §6.4).
  reach       `debench.Surface.reach()` in the production order, times the board/h1-test
              transfer factor K = (0.1890-0.1811)/(0.2029-0.1828) = 0.3930 (157 §1's own
              board pair), divided by the board reach span and by 6.  157's PROPORTIONAL power
              matching (f = 2, f = 3.5) is reported beside it as the rule-9 sensitivity.
  fid         read as PRECISION `mean_p(k/n_pred)`.  Rule 9: `n_real` is 68-322 on every
              harness against ~240 on the challenge with `n_pred` ~ 4,424, so the yield cap is
              saturated there and precision IS `fid`.
  jac         rebuilt on the BOARD's own counts from the measured recall:
              `J = recall*n_real / (n_pred + n_real - recall*n_real)`, n_real 240, n_pred 4424.
  nmae        `debench.Surface.nmae()`, span 0.600 (AGENTS.md §4 range 0.57-0.63), sign flipped.

MATCHED SEEDS.  157 measured a same-policy re-seed spread of 0.008-0.011 raw `reach` on this
surface -- the same size as every arm delta ever measured on it.  So every arm is paired with
its base at IDENTICAL seeds and the paired difference carries the error bar.  A single-seed
reading of `reach` here is not a measurement and is never quoted as one.

KNOWN-ANSWER ROWS (both abort):
  * 164/157's a15 from `de200.KA_a15_raw.npz`, which is byte-identical (sha256
    79c81c2338d4d489...) to `experiments/164-callset-artefact/de164.a15.npz`:
    `reach` 0.18283 (157 `power_h1test.log` f=1.0, 164 `decomp.log`:6) and median `n_pred`
    **6034** (164 `decomp.log`:29 -- the 6035 in `realize164.log` is the same table counted
    BEFORE the self-gene exclusion, and 6035 - 6034 is that one gene).
  * 169's own five-seed ladder, recomputed from `experiments/169-.../de169.*.npz` on 169's
    OWN zero-filled reach: d`nmae` +0.04640 and d`reach` +0.00818 for the q = 1.15 h-space arm.
"""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
WORK = os.path.join(HERE, "work")
E169 = os.path.join(ROOT, "experiments/169-sharpening-six-member")
sys.path.insert(0, os.path.join(ROOT, "experiments/153-de-bench"))
from debench import ANCHOR, Surface                                   # noqa: E402

ALPHA = 0.05
FACTORS = (1.0, 2.0, 3.5)
K_REACH = (0.1890 - 0.1811) / (0.2029 - 0.1828)      # 157 §1's own board pair -> 0.3930
SPAN_R = ANCHOR["reach"][1] - ANCHOR["reach"][0]
SPAN_F = ANCHOR["fid"][1] - ANCHOR["fid"][0]
SPAN_J = ANCHOR["jac"][1] - ANCHOR["jac"][0]
NMAE_B, NMAE_R = 1.0013, 0.400
SPAN_N = NMAE_B - NMAE_R
BOARD_NR, BOARD_NP = 240.0, 4424.0


class PowerSurface(Surface):
    """157's PROPORTIONAL power matching: pool_p = reference top-round(f * n_real_p)."""

    def __init__(self, path, factor=1.0):
        super().__init__(path)
        z = np.load(path, allow_pickle=True)
        rp = np.where(np.isfinite(z["real_p"]), z["real_p"], np.inf)
        base = (z["real_padj"] < ALPHA) & np.isfinite(z["real_padj"])
        sig = np.zeros(rp.shape, bool)
        for i in range(rp.shape[0]):
            k = int(round(factor * int(base[i].sum())))
            k = max(0, min(k, int(np.isfinite(z["real_p"][i]).sum())))
            if k:
                sig[i, np.argsort(rp[i], kind="stable")[:k]] = True
        self.real_sig = sig & np.isfinite(z["real_p"])
        for p, j in enumerate(self.self_idx):
            if j >= 0:
                self.real_sig[p, j] = False


def precision(s):
    out = []
    for p in range(len(s.perts)):
        sel = np.flatnonzero(s.pred_sig[p] & s.adj[p])
        if sel.size == 0:
            continue
        out.append(float((np.sign(s.pred_lfc[p, sel]) ==
                          np.sign(s.real_lfc[p, sel])).sum()) / sel.size)
    return float(np.mean(out))


def recall_of(s):
    out = []
    for p in range(len(s.perts)):
        nr = int(s.real_sig[p].sum())
        if nr == 0:
            continue
        out.append(int((s.real_sig[p] & s.pred_sig[p]).sum()) / nr)
    return float(np.mean(out))


def per_pert_reach(path, factor=1.0):
    """169's ZERO-FILLED per-perturbation reach.

    Kept verbatim so 169's ladder stays a known-answer row against the number 169 itself
    published.  It differs from `debench.Surface.reach()` in exactly one place: a perturbation
    whose reference pool is EMPTY contributes a 0 to the mean here and is DROPPED from the mean
    there.  On h1-test at f = 1 three of the 85 perturbations have a native BH set consisting
    of their own target gene alone, which the self-exclusion (BRIEF 4-7) then empties, so the
    two conventions differ by the exact constant 82/85 = 0.964706 -- a property of the
    REFERENCE, identical for every arm, so it rescales every delta and reverses no sign.  At
    f >= 2 the pools are no longer empty and the two agree to every printed digit.
    """
    s = PowerSurface(path, factor)
    out = np.zeros(len(s.perts))
    for p in range(len(s.perts)):
        pool = np.flatnonzero(s.real_sig[p])
        if pool.size == 0:
            continue
        order = s._order(p, pool, None, None)
        pool = pool[order]
        pool = pool[s.adj[p, pool]]
        if pool.size == 0:
            continue
        ok = np.sign(s.pred_lfc[p, pool]) == np.sign(s.real_lfc[p, pool])
        purity = np.cumsum(ok) / np.arange(1, pool.size + 1)
        deep = np.flatnonzero(purity >= 0.9)
        out[p] = (int(deep[-1] + 1) if deep.size else 0) / int(s.real_sig[p].sum())
    return out


def empty_pools(path, factor=1.0):
    s = PowerSurface(path, factor)
    return int(sum(int(s.real_sig[p].sum()) == 0 for p in range(len(s.perts)))), len(s.perts)


def row_for(path):
    """Both reach conventions, because two different known-answer rows need them.

    `reach_f*`     `debench.Surface.reach()` -- the production instrument.  157's
                   `power_h1test.log`, 164's `decomp.log` and every board-linked `reach`
                   number in `docs/BOARD-LEDGER.md` are on this one, so the ablation table is.
    `reach_f*_zf`  169's zero-filled variant, so 169's published +0.00818 stays checkable.
    """
    s = Surface(path)
    nr, npd = s.counts()
    row = dict(n_pred=float(np.median(npd)), n_real=float(np.median(nr)),
               n_real_mean=float(nr.mean()), precision=precision(s), fid=s.fid(),
               jac=s.jac(), nmae=s.nmae(), sign=s.sign_accuracy(), recall=recall_of(s))
    for f in FACTORS:
        row[f"reach_f{f}"] = float(PowerSurface(path, f).reach())
        row[f"reach_f{f}_zf"] = float(per_pert_reach(path, f).mean())
    return row


def jac_board(recall, n_pred=BOARD_NP):
    hit = recall * BOARD_NR
    return hit / (n_pred + BOARD_NR - hit)


def overall(row, base):
    """The four DE members of one arm against its base, in overall board units.

    TWO `jac` CONVENTIONS, because the arms that move the call-set SIZE make them disagree in
    sign and nothing in this experiment can adjudicate between them:

      `jac`     BOARD_NP held at 4,424 for both members (169's convention, and the incumbent).
      `jac_np`  BOARD_NP scaled by the arm's own measured h1-test `n_pred` ratio, which is what
                you would write down if the board's call set tracked ours proportionally.

    The one board point that constrains the choice is a21 (`resample_x 0.14`), whose h1-test
    `n_pred` fell 6,034 -> 1,913 at recall 0.669 -> 0.271 and whose board `jac` moved **-0.017
    scaled**.  Fixed-NP predicts -0.0586 scaled (right sign, 3.4x too large); scaled-NP predicts
    **+0.0150** (wrong sign).  So fixed-NP is the headline and scaled-NP is the sensitivity, and
    any verdict that flips between them is reported as unresolved and handed to 192.
    """
    dr = row["reach_f1.0"] - base["reach_f1.0"]
    r_np = BOARD_NP
    b_np = BOARD_NP * (base["n_pred"] / row["n_pred"]) if row.get("n_pred") else BOARD_NP
    return dict(
        d_reach_raw=dr,
        reach=dr * K_REACH / SPAN_R / 6,
        fid=(row["precision"] - base["precision"]) / SPAN_F / 6,
        jac=(jac_board(row["recall"]) - jac_board(base["recall"])) / SPAN_J / 6,
        jac_np=(jac_board(row["recall"], r_np)
                - jac_board(base["recall"], b_np)) / SPAN_J / 6,
        nmae=-(row["nmae"] - base["nmae"]) / SPAN_N / 6,
        d_nmae_raw=row["nmae"] - base["nmae"])


def paired(rows, arm_names, base_names, field):
    """Seed-paired difference: the only honest error bar on this surface (157)."""
    n = min(len(arm_names), len(base_names))
    d = np.asarray([rows[arm_names[i]][field] - rows[base_names[i]][field]
                    for i in range(n)])
    sem = float(d.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan")
    return dict(mean=float(d.mean()), sd=float(d.std(ddof=1)) if n > 1 else float("nan"),
                sem=sem, n=n, t=float(d.mean() / sem) if n > 1 and sem > 0 else float("nan"),
                lo=float(d.mean() - 2.776 * sem) if n > 1 else float("nan"),
                hi=float(d.mean() + 2.776 * sem) if n > 1 else float("nan"))


# ---------------------------------------------------------------------------- arm families
def seed_names(stem, n):
    return [stem if i == 0 else f"{stem}_s{i+1}" for i in range(n)]


COMPARISONS = [
    # (label, arm stem, base stem, n seeds, what the component IS)
    ("a23 mass fix",            "art_mf",        "art_nomf",       5, "LOO: mass fix ON vs OFF"),
    ("tiered spread policy",    "art_mf",        "art_flat_mf",    5, "LOO: tier8000 vs flat 0.25"),
    ("--overdispersion 4.0",    "art_mf",        "art_phi0_mf",    5, "LOO: phi 4 vs phi 0"),
    ("tier + phi together",     "art_mf",        "art_flatphi0_mf", 5, "both off = v11's law"),
    ("tier8000 vs flat 1.00",   "art_mf",        "art_flat100_mf", 3, "both clear the mse cliff"),
    ("tier8000 vs tier4000",    "art_mf",        "art_tier4000_mf", 3, "both clear the mse cliff"),
    # a8's OWN board pair for `--overdispersion`: a8 carried a FLAT 0.25 spread and no tier
    # (BOARD-LEDGER §3), so the board's +0.021398 scaled `reach` for phi 4 was measured in the
    # flat context, not in a15's tiered one.  Retrodicting it needs the flat pair.
    ("phi 4 IN a8's flat context", "art_flat_mf", "art_flatphi0_mf", 5, "a8's own board pair"),
    ("flat 1.00 vs flat 0.25",  "art_flat100_mf", "art_flat_mf",   3, "the spread scalar alone"),
    ("whitening (gamma 0.65)",  "reb_ship",      "reb_noWhiten",   3, "LOO on the rebuild"),
    ("PerturbMatch blend 0.5",  "reb_ship",      "reb_noPM",       3, "LOO on the rebuild"),
    ("CD4 residual f=0.10",     "reb_ship",      "reb_noCD4",      3, "LOO (coverage 4/85!)"),
    ("equal donor weights",     "reb_ship",      "reb_w41",        3, "LOO vs a 4:1 K562 mix"),
    ("equal vs 2:1 K562-heavy", "reb_ship",      "reb_w21",        3, "175's proposed rider"),
    ("gamma 0.65 vs 0.80",      "reb_ship",      "reb_g080",       3, "the a19 direction"),
    ("gamma 0.65 vs 0.50",      "reb_ship",      "reb_g050",       3, "the other side of 170"),
    # The REBUILD baseline is built at the builder's DEFAULT --min-cpm 5.0 (realize200.py:112),
    # while the shipped arm uses 10.0 (pds200.py:60, mse200.py:55).  So the shipped gate's own
    # leave-one-out is the gate10-vs-gate0 pair, not anything anchored on `reb_ship`.
    ("CPM gate 10 (shipped) LOO", "reb_ship_gate10", "reb_ship_gate0", 3, "gate 10 vs no gate"),
    ("CPM gate 10 vs 5 (088)",  "reb_ship_gate10", "reb_ship",       3, "088's own decision"),
    ("CPM gate 10 vs 20",       "reb_ship_gate10", "reb_ship_gate20", 3, ""),
    ("CPM gate 5 vs none",      "reb_ship",      "reb_ship_gate0", 3, "the rebuild's own gate"),
    ("rebuild vs artifact geom", "reb_ship",     "reb_e095",       3, "the full stack vs 095's"),
]


def main() -> None:
    paths = {os.path.basename(p)[len("de200."):-4]: p
             for p in sorted(glob.glob(os.path.join(HERE, "de200.*.npz")))}
    if not paths:
        sys.exit("no de200.*.npz -- run realize200.py first")
    cache_path = os.path.join(WORK, "six200v2.json")
    rows = json.loads(open(cache_path).read()) if os.path.exists(cache_path) else {}
    for a, p in paths.items():
        if a in rows:
            continue
        rows[a] = row_for(p)
        json.dump(rows, open(cache_path, "w"), indent=1, default=float)
        print(f"  scored {a}", flush=True)

    print(f"\n{'=' * 128}\nRAW MEMBERS PER ARM   (h1-test, emitted, reference BH<0.05 pool)")
    print(f"{'arm':<24}{'n_pred':>8}{'n_real':>8}{'reach f1':>10}{'reach f2':>10}"
          f"{'reach f3.5':>11}{'precision':>11}{'fid raw':>9}{'jac raw':>9}{'recall':>8}"
          f"{'nmae':>9}{'sign':>8}")
    for a in sorted(rows):
        r = rows[a]
        print(f"{a:<24}{r['n_pred']:>8.0f}{r['n_real']:>8.0f}{r['reach_f1.0']:>10.5f}"
              f"{r['reach_f2.0']:>10.5f}{r['reach_f3.5']:>11.5f}{r['precision']:>11.5f}"
              f"{r['fid']:>9.5f}{r['jac']:>9.5f}{r['recall']:>8.5f}{r['nmae']:>9.5f}"
              f"{r['sign']:>8.5f}")

    # ---------------------------------------------------------------- known-answer gates
    print(f"\n{'=' * 128}\nKNOWN-ANSWER GATES")
    ok = True
    if "KA_a15_raw" in rows:
        g = rows["KA_a15_raw"]
        e, n = empty_pools(paths["KA_a15_raw"], 1.0)
        print(f"  empty reference pools at f=1: {e}/{n}  -> zf/debench ratio {(n - e) / n:.6f}"
              f"   (zf reach {g['reach_f1.0_zf']:.5f} = {g['reach_f1.0']:.5f} x that)")
        for label, got, want, tol in (("median n_pred", g["n_pred"], 6034, 0.5),
                                      ("reach f=1", g["reach_f1.0"], 0.18283, 5e-4)):
            good = abs(got - want) <= tol
            ok &= good
            print(f"  164/157 a15 {label:<16} got {got:.5f}  expect {want:.5f}  "
                  f"{'PASS' if good else 'FAIL'}")
    # 169's own ladder, recomputed here with THIS code
    d169 = {os.path.basename(p)[len("de169."):-4]: p
            for p in sorted(glob.glob(os.path.join(E169, "de169.*.npz")))}
    have = [n for n in ("q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3", "q1.00_mf_s4",
                        "q1.00_mf_s5", "q1.15_mf", "q1.15_mf_b", "q1.15_mf_s3",
                        "q1.15_mf_s4", "q1.15_mf_s5") if n in d169]
    if len(have) == 10:
        c169 = json.loads(open(os.path.join(WORK, "six200-169v2.json")).read()) \
            if os.path.exists(os.path.join(WORK, "six200-169v2.json")) else {}
        for n in have:
            if n not in c169:
                c169[n] = row_for(d169[n])
                json.dump(c169, open(os.path.join(WORK, "six200-169v2.json"), "w"), indent=1,
                          default=float)
        A = ["q1.15_mf", "q1.15_mf_b", "q1.15_mf_s3", "q1.15_mf_s4", "q1.15_mf_s5"]
        B = ["q1.00_mf", "q1.00_mf_s2", "q1.00_mf_s3", "q1.00_mf_s4", "q1.00_mf_s5"]
        pr = paired(c169, A, B, "reach_f1.0_zf")          # 169's own convention
        prd = paired(c169, A, B, "reach_f1.0")            # this table's convention
        pn = paired(c169, A, B, "nmae")
        print(f"  169's d reach under THIS table's convention: {prd['mean']:+.5f} "
              f"(= {pr['mean']:+.5f} / 0.964706)")
        for label, got, want, tol in (("169 d reach (5 seeds)", pr["mean"], 0.00818, 5e-4),
                                      ("169 d nmae  (5 seeds)", pn["mean"], 0.04640, 5e-4)):
            good = abs(got - want) <= tol
            ok &= good
            print(f"  {label:<28} got {got:+.5f}  expect {want:+.5f}  "
                  f"{'PASS' if good else 'FAIL'}")
    else:
        print(f"  (169's ladder not fully on disk: {len(have)}/10 tables)")
    print(f"  transfer factor K_reach = {K_REACH:.4f}; spans reach {SPAN_R:.5f} "
          f"fid {SPAN_F:.5f} jac {SPAN_J:.5f} nmae {SPAN_N:.5f}")
    if not ok:
        print("  *** A GATE FAILED -- everything below is void ***")

    # ---------------------------------------------------------------- the ablation table
    print(f"\n{'=' * 128}\nLEAVE-ONE-OUT, MATCHED SEEDS, in OVERALL board units")
    print("(sign convention: POSITIVE = the component EARNS its place on that member)")
    print(f"{'component':<28}{'seeds':>6}{'d reach raw':>12}{'sem':>9}{'t':>7}"
          f"{'reach':>10}{'fid':>10}{'jac':>10}{'nmae':>10}{'DE sum':>11}"
          f"{'jac(np)':>10}{'sum(np)':>11}")
    out = {}
    for label, astem, bstem, n, note in COMPARISONS:
        A = [x for x in seed_names(astem, n) if x in rows]
        B = [x for x in seed_names(bstem, n) if x in rows]
        m = min(len(A), len(B))
        if m == 0:
            continue
        A, B = A[:m], B[:m]
        pr = paired(rows, A, B, "reach_f1.0")
        mean_a = {k: float(np.mean([rows[x][k] for x in A]))
                  for k in ("reach_f1.0", "reach_f2.0", "reach_f3.5", "precision", "recall",
                            "nmae", "fid", "jac", "n_pred", "sign")}
        mean_b = {k: float(np.mean([rows[x][k] for x in B]))
                  for k in ("reach_f1.0", "reach_f2.0", "reach_f3.5", "precision", "recall",
                            "nmae", "fid", "jac", "n_pred", "sign")}
        o = overall(mean_a, mean_b)
        s = o["reach"] + o["fid"] + o["jac"] + o["nmae"]
        o["de_sum_np"] = o["reach"] + o["fid"] + o["jac_np"] + o["nmae"]
        out[label] = dict(**o, de_sum=s, paired_reach=pr, note=note,
                          arms=A, bases=B, mean_arm=mean_a, mean_base=mean_b,
                          d_reach_f2=mean_a["reach_f2.0"] - mean_b["reach_f2.0"],
                          d_reach_f35=mean_a["reach_f3.5"] - mean_b["reach_f3.5"],
                          d_npred=mean_a["n_pred"] - mean_b["n_pred"],
                          d_sign=mean_a["sign"] - mean_b["sign"])
        print(f"{label:<28}{m:>6}{pr['mean']:>+12.5f}{pr['sem']:>9.5f}{pr['t']:>7.2f}"
              f"{o['reach']:>+10.6f}{o['fid']:>+10.6f}{o['jac']:>+10.6f}{o['nmae']:>+10.6f}"
              f"{s:>+11.6f}{o['jac_np']:>+10.6f}{o['de_sum_np']:>+11.6f}")

    print(f"\nrule-9 sensitivity: the same d reach under 157's PROPORTIONAL power matching")
    print(f"{'component':<28}{'f=1':>12}{'f=2':>12}{'f=3.5':>12}{'d n_pred':>11}"
          f"{'d sign acc':>12}")
    for label in out:
        r = out[label]
        print(f"{label:<28}{r['d_reach_raw']:>+12.5f}{r['d_reach_f2']:>+12.5f}"
              f"{r['d_reach_f35']:>+12.5f}{r['d_npred']:>+11.0f}{r['d_sign']:>+12.5f}")

    json.dump(out, open(os.path.join(WORK, "six200-tablev2.json"), "w"), indent=1, default=float)
    print(f"\nwrote {os.path.join(WORK, 'six200-tablev2.json')}")


if __name__ == "__main__":
    sys.exit(main())
