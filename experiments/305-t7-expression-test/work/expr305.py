"""305 step 4 -- the frozen test (PREREGISTRATION sec 4-9) on the 15 candidate pairs.

    PYTHONDONTWRITEBYTECODE=1 uv run python expr305.py

Refuses to run unless every known-answer row of `ka305.json` passes.  Emits
`work/expression305.csv` and `work/results305.json`.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import lib305 as X  # noqa: E402

E298 = X.ROOT / "experiments/298-reviewer-checks/work"
MIN_NULL = 30          # PREREGISTRATION sec 6
MIN_CELLS_HALF = 30    # PREREGISTRATION sec 7
MIN_CP50K = 0.5        # PREREGISTRATION sec 7
MAG = -0.10            # PREREGISTRATION sec 7, C2
P_RAW = 0.01           # PREREGISTRATION sec 7, C4


def holm(ps: list[float]) -> list[float]:
    """Holm step-down adjusted p-values, order preserved, monotone, capped at 1."""
    m = len(ps)
    order = sorted(range(m), key=lambda i: ps[i])
    adj = [0.0] * m
    run = 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (m - k) * ps[i]))
        adj[i] = run
    return adj


class Split:
    def __init__(self, s: str):
        z = np.load(X.HERE / f"sums305_{s}.npz", allow_pickle=False)
        self.meta = json.loads((X.HERE / f"sums305_{s}.json").read_text())
        self.name = s
        self.genes = z["genes"].astype(str)
        self.gpos = {g: i for i, g in enumerate(self.genes)}
        names, halves, S = z["names"].astype(str), z["halves"], z["sums"]
        self.cells = {}
        self.umis = {}
        self.half = {}
        self.tot = {}
        for i, (c, h) in enumerate(zip(names, halves)):
            self.half[(c, int(h))] = S[i]
            self.tot[c] = self.tot.get(c, 0) + S[i]
            self.cells[(c, int(h))] = int(z["cells"][i])
            self.umis[(c, int(h))] = int(z["umis"][i])
        self.constructs = self.meta["constructs"]
        self.tgt_of = self.meta["target_of_construct"]
        self.ntc = self.meta["ntc_constructs"]
        self.par = self.meta["ntc_parity"]
        self.targeting = [c for c in self.constructs if self.tgt_of[c] != X.CTRL]
        self.cA = sum(self.tot[c] for c in self.ntc if self.par[c] == 0)
        self.cB = sum(self.tot[c] for c in self.ntc if self.par[c] == 1)
        self.cAB = self.cA + self.cB
        self.cidx = {c: i for i, c in enumerate(self.constructs)}
        self.Pmat = X.L.prof(np.stack([self.tot[c] for c in self.constructs]))   # profiles, control-free
        self.lanes = self.meta["lane_composition"]

    def refs(self, c: str):
        """The (A, B, pooled) control sums for testing construct `c`; leave-one-construct-out when `c`
        is itself a non-targeting construct (PREREGISTRATION sec 5)."""
        if c in self.par:
            p = self.par[c]
            a = self.cA - (self.tot[c] if p == 0 else 0)
            b = self.cB - (self.tot[c] if p == 1 else 0)
            return a, b, a + b, True
        return self.cA, self.cB, self.cAB, False

    def col(self, gene: str, ctrl) -> np.ndarray:
        """Every construct's direct effect on `gene`, against the given control sums."""
        j = self.gpos[gene]
        return self.Pmat[:, j] - float(X.L.prof(ctrl)[j])

    def lane_imbalance(self, c: str) -> int:
        a = self.lanes[f"{c}|0"]
        b = self.lanes[f"{c}|1"]
        return max(abs(a.get(k, 0) - b.get(k, 0)) for k in set(a) | set(b))


def evaluate(sp: Split, c: str, gene: str, all_splits: dict) -> dict:
    j = sp.gpos[gene]
    cA, cB, cAB, loo = sp.refs(c)
    pAB = X.L.prof(cAB)
    e_pool = float(X.L.prof(sp.tot[c])[j] - pAB[j])
    e0 = float(X.L.prof(sp.half[(c, 0)])[j] - X.L.prof(cA)[j])
    e1 = float(X.L.prof(sp.half[(c, 1)])[j] - X.L.prof(cB)[j])

    col = sp.col(gene, cAB)
    keep = [c2 for c2 in sp.targeting if c2 != c and sp.tgt_of[c2] != gene]
    null = np.asarray([col[sp.cidx[c2]] for c2 in keep])
    n_null = int(null.size)
    panel_mean = float(np.mean([col[sp.cidx[c2]] for c2 in sp.targeting]))
    resid = e_pool - panel_mean
    p_raw = float((1 + int((null <= e_pool).sum())) / (1 + n_null)) if n_null else float("nan")
    med, mad = float(np.median(null)), float(np.median(np.abs(null - np.median(null))))
    z_rob = float((e_pool - med) / (1.4826 * mad)) if mad > 0 else float("nan")
    rank = int((null <= e_pool).sum()) + 1

    # secondary, all-splits null (declared in advance; the three splits share the same NTC cells)
    nall = []
    for s2, sp2 in all_splits.items():
        c2A, c2B, c2AB, _ = sp2.refs(c) if s2 == sp.name else (sp2.cA, sp2.cB, sp2.cAB, False)
        col2 = sp2.col(gene, c2AB)
        nall += [col2[sp2.cidx[x]] for x in sp2.targeting if not (s2 == sp.name and x == c)
                 and sp2.tgt_of[x] != gene]
    nall = np.asarray(nall)
    p_all = float((1 + int((nall <= e_pool).sum())) / (1 + nall.size))

    ctrl_cp50k = float(X.L.TS * cAB[j] / cAB.sum())
    con_cp50k = float(X.L.TS * sp.tot[c][j] / sp.tot[c].sum())
    n0, n1 = sp.cells[(c, 0)], sp.cells[(c, 1)]
    untestable = []
    if ctrl_cp50k < MIN_CP50K:
        untestable.append(f"control CP50k {ctrl_cp50k:.3f} < {MIN_CP50K}")
    if min(n0, n1) < MIN_CELLS_HALF:
        untestable.append(f"half cells {min(n0, n1)} < {MIN_CELLS_HALF}")
    if n_null < MIN_NULL:
        untestable.append(f"null {n_null} < {MIN_NULL}")

    return dict(
        split=sp.name, construct=c, candidate_gene=gene, testable=not untestable,
        untestable_reason="; ".join(untestable), loo_control=loo,
        e_pool=e_pool, e_half0=e0, e_half1=e1,
        fold_change=(con_cp50k / ctrl_cp50k) if ctrl_cp50k > 0 else float("nan"),
        ctrl_cp50k=ctrl_cp50k, construct_cp50k=con_cp50k,
        panel_mean=panel_mean, residual=resid,
        n_null=n_null, null_median=med, null_mad=mad, null_min=float(null.min()) if n_null else float("nan"),
        rank_in_null=rank, p_raw=p_raw, z_rob=z_rob,
        n_null_all=int(nall.size), p_raw_all_splits=p_all,
        cells_half0=n0, cells_half1=n1, umi_half0=sp.umis[(c, 0)], umi_half1=sp.umis[(c, 1)],
        lane_imbalance=sp.lane_imbalance(c), n_lanes=len(sp.lanes[f"{c}|0"]),
        C1_negative=e_pool < 0, C2_magnitude=e_pool <= MAG, C3_both_halves=(e0 < 0 and e1 < 0),
        C4_unusual=(p_raw <= P_RAW) if n_null else False)


def main() -> None:
    P = X.L.Logger(X.HERE / "expr305.log")
    ka = json.loads((X.HERE / "ka305.json").read_text())
    if ka["failures"]:
        raise SystemExit(f"known-answer failures {ka['failures']} -- conclusions aborted")
    P(f"KA: {len(ka['ka'])}/{len(ka['ka'])} rows pass; proceeding.")

    SP = {s: Split(s) for s in X.SPLITS}
    pairs = list(csv.DictReader(open(X.HERE / "candidates305.csv")))
    arms = list(csv.DictReader(open(E298 / "seed298.csv")))
    partner = {}
    for a in arms:
        for b in arms:
            if a["target"] == b["target"] and a["construct"] != b["construct"]:
                partner[a["construct"]] = b["construct"]

    rows = []
    for r in pairs:
        c, g, grp = r["construct"], r["candidate_gene"], r["group"]
        splits = [r["split"]] if grp != "ntc" else list(X.SPLITS)
        for s in splits:
            if c not in SP[s].tot:
                continue
            d = evaluate(SP[s], c, g, SP)
            d.update(group=grp, pair_id=r["pair_id"], intended_target=r["intended_target"],
                     kd_lfc=r["kd_lfc"], offsets=r["offsets"], n_hits=r["n_hits"])
            # paired-construct comparator (PREREGISTRATION sec 9)
            pc = partner.get(c)
            if pc and pc in SP[s].tot:
                q = evaluate(SP[s], pc, g, SP)
                d.update(partner_construct=pc, partner_e_pool=q["e_pool"], partner_p_raw=q["p_raw"],
                         partner_e_half0=q["e_half0"], partner_e_half1=q["e_half1"],
                         partner_z_rob=q["z_rob"], partner_rank=q["rank_in_null"])
            else:
                d.update(partner_construct="", partner_e_pool=float("nan"), partner_p_raw=float("nan"),
                         partner_e_half0=float("nan"), partner_e_half1=float("nan"),
                         partner_z_rob=float("nan"), partner_rank=-1)
            # control-vs-control comparator: one NTC half referenced to the other
            d["ctrl_vs_ctrl"] = float(X.effect(SP[s].cA, SP[s].cB)[SP[s].gpos[g]])
            rows.append(d)

    # ---- Holm inside each family, over the TESTABLE pairs (family size recorded both ways)
    fam = {"primary": [r for r in rows if r["group"] in ("failing", "working")],
           "comparator": [r for r in rows if r["group"] == "ntc"]}
    holm_meta = {}
    for name, members in fam.items():
        tv = [r for r in members if r["testable"]]
        adj = holm([r["p_raw"] for r in tv])
        adj_all = holm([r["p_raw_all_splits"] for r in tv])
        for r, a, b in zip(tv, adj, adj_all):
            r["p_holm"] = a
            r["p_holm_all_splits"] = b
        for r in members:
            r.setdefault("p_holm", float("nan"))
            r.setdefault("p_holm_all_splits", float("nan"))
        floor = max((1.0 / (1 + r["n_null"]) for r in tv), default=float("nan"))
        best = min((1.0 / (1 + r["n_null"]) for r in tv), default=float("nan"))
        holm_meta[name] = dict(n_members=len(members), n_testable=len(tv),
                               smallest_attainable_p=best,
                               smallest_attainable_holm=best * len(tv),
                               holm005_attainable=bool(best * len(tv) < 0.05),
                               smallest_attainable_holm_all_splits=min(
                                   (1.0 / (1 + r["n_null_all"]) for r in tv), default=float("nan")) * len(tv))

    # ---- positive control: the identical statistic on each arm's OWN target gene
    pos = []
    for a in arms:
        s, c, t = a["split"], a["construct"], a["target"]
        d = evaluate(SP[s], c, t, SP)
        pos.append(dict(construct=c, split=s, gene=t, failing=a["failing"] == "True",
                        e_pool=d["e_pool"], e_half0=d["e_half0"], e_half1=d["e_half1"],
                        p_raw=d["p_raw"], rank=d["rank_in_null"], n_null=d["n_null"],
                        z_rob=d["z_rob"], fold_change=d["fold_change"],
                        C1=d["C1_negative"], C2=d["C2_magnitude"], C3=d["C3_both_halves"],
                        C4=d["C4_unusual"]))
    npass = sum(1 for p in pos if p["C1"] and p["C2"] and p["C3"] and p["C4"])
    nfail_arm = sum(1 for p in pos if p["failing"])
    P(f"\nPOSITIVE CONTROL (same frozen statistic on each arm's OWN target gene): "
      f"{npass}/34 arms meet all four criteria "
      f"({sum(1 for p in pos if p['failing'] and p['C1'] and p['C2'] and p['C3'] and p['C4'])}/{nfail_arm} "
      f"failing, {npass - sum(1 for p in pos if p['failing'] and p['C1'] and p['C2'] and p['C3'] and p['C4'])}"
      f"/{34 - nfail_arm} working)")
    for p in sorted(pos, key=lambda p: p["e_pool"])[:6]:
        P(f"    {p['construct']:26s} {p['gene']:8s} e {p['e_pool']:+.3f} fc {p['fold_change']:.3f} "
          f"rank {p['rank']}/{p['n_null'] + 1} p_raw {p['p_raw']:.4f} z {p['z_rob']:+.1f}")

    P("\nCANDIDATE PAIRS (frozen test):")
    hdr = (f"  {'group':8s} {'construct':26s} {'gene':9s} {'split':11s} {'e_pool':>8s} {'fc':>6s} "
           f"{'half0':>8s} {'half1':>8s} {'rank':>9s} {'p_raw':>7s} {'p_holm':>7s} {'z':>6s}  flags")
    P(hdr)
    for r in rows:
        flags = "".join(k for k, v in (("1", r["C1_negative"]), ("2", r["C2_magnitude"]),
                                       ("3", r["C3_both_halves"]), ("4", r["C4_unusual"])) if v) or "-"
        P(f"  {r['group']:8s} {r['construct'][:26]:26s} {r['candidate_gene']:9s} {r['split']:11s} "
          f"{r['e_pool']:+8.4f} {r['fold_change']:6.3f} {r['e_half0']:+8.4f} {r['e_half1']:+8.4f} "
          f"{r['rank_in_null']:4d}/{r['n_null'] + 1:<4d} {r['p_raw']:7.4f} {r['p_holm']:7.4f} "
          f"{r['z_rob']:+6.2f}  C{flags}")

    corrob = [r for r in rows if r["group"] in ("failing", "working") and r["testable"]
              and r["C1_negative"] and r["C2_magnitude"] and r["C3_both_halves"] and r["C4_unusual"]]
    comp = [r for r in rows if r["group"] == "ntc" and r["testable"]
            and r["C1_negative"] and r["C2_magnitude"] and r["C3_both_halves"] and r["C4_unusual"]]
    P(f"\nCORROBORATED primary pairs (C1&C2&C3&C4): {len(corrob)} of "
      f"{len([r for r in rows if r['group'] in ('failing', 'working')])} "
      f"({[r['pair_id'] for r in corrob]})")
    P(f"NTC comparator pairs meeting the same conjunction: {len(comp)} of "
      f"{len([r for r in rows if r['group'] == 'ntc'])} split-evaluations")
    P(f"Holm attainability: {json.dumps(holm_meta, indent=1, default=float)}")

    cvc = np.asarray([r["ctrl_vs_ctrl"] for r in rows])
    P(f"control-vs-control effects on the 15 candidate genes: median {np.median(cvc):+.4f}, "
      f"range [{cvc.min():+.4f}, {cvc.max():+.4f}]")

    fields = sorted({k for r in rows for k in r})
    order = ["pair_id", "group", "construct", "intended_target", "candidate_gene", "split", "testable",
             "untestable_reason", "e_pool", "fold_change", "e_half0", "e_half1", "C1_negative",
             "C2_magnitude", "C3_both_halves", "C4_unusual", "rank_in_null", "n_null", "p_raw", "p_holm",
             "p_raw_all_splits", "p_holm_all_splits", "n_null_all", "z_rob", "residual", "panel_mean",
             "null_median", "null_mad", "null_min", "ctrl_cp50k", "construct_cp50k", "ctrl_vs_ctrl",
             "partner_construct", "partner_e_pool", "partner_e_half0", "partner_e_half1",
             "partner_p_raw", "partner_rank", "partner_z_rob", "cells_half0", "cells_half1",
             "umi_half0", "umi_half1", "n_lanes", "lane_imbalance", "loo_control", "kd_lfc",
             "offsets", "n_hits"]
    order += [f for f in fields if f not in order]
    with open(X.HERE / "expression305.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=order)
        w.writeheader()
        w.writerows(rows)

    res = dict(
        preregistration_sha256=(X.HERE.parent / "PREREGISTRATION.sha256").read_text().split()[0],
        ka=[dict(label=k["label"], ok=k["ok"]) for k in ka["ka"]],
        thresholds=dict(MAG=MAG, P_RAW=P_RAW, MIN_NULL=MIN_NULL, MIN_CELLS_HALF=MIN_CELLS_HALF,
                        MIN_CP50K=MIN_CP50K),
        counts=dict(pairs_total=len(rows),
                    primary=len([r for r in rows if r["group"] in ("failing", "working")]),
                    failing=len([r for r in rows if r["group"] == "failing"]),
                    working=len([r for r in rows if r["group"] == "working"]),
                    ntc_evaluations=len([r for r in rows if r["group"] == "ntc"]),
                    untestable=len([r for r in rows if not r["testable"]]),
                    corroborated=len(corrob), ntc_meeting_conjunction=len(comp)),
        holm=holm_meta, rows=rows, positive_control=pos,
        positive_control_summary=dict(
            n_pass=npass, n_pass_failing=sum(1 for p in pos if p["failing"] and p["C1"] and p["C2"]
                                             and p["C3"] and p["C4"]), n_arms=34, n_failing=nfail_arm),
        control_vs_control=dict(median=float(np.median(cvc)), min=float(cvc.min()), max=float(cvc.max())),
        ka6=ka["ka6"], peak_rss_gib=round(X.peak_gib(), 3))
    X.dump(res, X.HERE / "results305.json")
    P(f"\nwrote expression305.csv ({len(rows)} rows) and results305.json; peak {X.peak_gib():.2f} GiB")
    P.save()


if __name__ == "__main__":
    main()
