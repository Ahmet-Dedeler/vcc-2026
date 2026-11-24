#!/usr/bin/env python3
"""7772dc's independent read-only review of 630's rank-neutral UP pass (asked by [d4cd5f], 2026-10-02 ~13:20Z).

Property tests on SYNTHETIC blocks (no repo data). Three module versions are loaded side by side:
  v2   = git 65ee4584:experiments/630-rank-neutral-up/src_rank_up.py (sha 5ae935b9, the KA1 module);
  v3   = git ae39c950 (sha 8e1f05a6, the version named in the review request);
  cur  = the working copy (605-fixed; its sha is printed).
Checks, each on many random blocks, over variants b/s × K × v3 suffixes:
  P1  libraries and every per-gene total exact; counts are non-negative integers (independent recomputation).
  P2  every kept gene's net U change vs the input is within capU/2 = floor(1e-4·2·n·n_ctrl)/2, recomputed with scipy's
      Mann-Whitney U (independent of ZS.ranks_of).
  P3  0 call flips under the module's de_stats (the 323 scorer emulation).
  P4  no kept gene's arithmetic mean falls; no delivered lfc sign changes; DOWN, called and own-target genes untouched.
  P5  sub-gate counts only leave carrier rows and only arrive on ballast rows; kept genes only gain on carriers and only
      lose on ballast.
  P6  determinism (two runs bit-identical) and no RNG consumed (numpy's global state unchanged).
  P7  v2 equivalence on deep-ballast blocks: cur(b_k100) == v2(b_k100) and cur(s_k100) == v2(s_k100) bit for bit.
      cur(suffix) == v3(suffix), i.e. the 605 fix is a no-op there.
  P8  hostile blocks: shallow ballast rows, ballast_k = 0, a top carrier with zero sub-gate counts, heavy ties.
      cur must never raise and must keep P1-P5. v2/v3 behaviour is reported (raise / mean falls / sign flips).
  P9  fail-safe guard: with de_stats forced to report a flip on every re-check, cur returns the input block unchanged,
      with ur_guard_failed = 1 and the SAME key set as a normal call (the builder's summary aggregates by key).
Analysis only, < 1 GB, no repo data.
usage: uv run python experiments/600-breakthrough-council/review630_7772dc.py [--rich]
  --rich: production-like richness (G 2,400 with 1,600 sub-gate genes, 400-cell blocks with 10 deep ballast rows,
          ~90 sub-gate counts per carrier, 1,500 controls), so many genes move and v3's later rounds bind; 10 blocks.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import mannwhitneyu

from vcc2026 import rank_mag as ZS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SRC = "experiments/630-rank-neutral-up/src_rank_up.py"
VARIANTS = ["b_k100", "s_k100", "b_k075", "b_k150", "s_k250", "b_k050", "b_k100_q20r3p10", "b_k100_q1r5p3",
            "s_k100_q2r4p1", "b_k300_q50r2p20", "b_k075_q20r3p10"]


def load_version(name: str, text: str, td: Path):
    p = td / f"rank_up_{name}.py"
    p.write_text(text)
    spec = importlib.util.spec_from_file_location(f"rank_up_{name}", p)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    m._SHA = hashlib.sha256(text.encode()).hexdigest()[:16]
    return m


RICH = "--rich" in sys.argv


def controls(rng, G=420, n_ctrl=700, frac_sub=0.35):
    """Control counts: a log-normal profile; ~frac_sub of the genes at 1-5 CPM (sub-gate, never tested)."""
    n_sub = int(G * frac_sub)
    p = np.empty(G)
    p[:n_sub] = rng.uniform(1.0, 4.5, n_sub) * 1e-6
    hi = rng.lognormal(0, 1.2, G - n_sub)
    p[n_sub:] = hi / hi.sum() * (1 - p[:n_sub].sum())
    perm = rng.permutation(G)
    p = p[perm]
    lib = np.round(rng.lognormal(np.log(20000), 0.45, n_ctrl)).astype(np.int64)
    X = rng.poisson(lib[:, None] * p[None, :])
    return sp.csr_matrix(X.astype(np.int64)), p, lib


def block(rng, p, lib_ctrl, n=100, n_bal=10, own=None, shallow=0, zero_sub_top=False, ties=False, keep=None):
    """A carrier + ballast block: carriers at control-like depths, ballast rows deep (~5e5-9.9e5), with per-gene effects."""
    G = p.size
    eff = np.zeros(G)
    m = rng.random(G) < 0.4
    eff[m] = rng.normal(0, 0.45, m.sum())
    if own is not None:
        eff[own] = -2.5
    n_car = n - n_bal
    lib = np.empty(n, dtype=np.int64)
    lib[:n_car] = rng.choice(lib_ctrl, n_car)
    lib[n_car:] = rng.integers(500_000, 990_000, n_bal)
    if shallow:                                                # hostile: some ballast rows shallower than most carriers
        lib[n_car:n_car + shallow] = rng.integers(3_000, 9_000, shallow)
    mu = lib[:, None] * (p * np.exp(eff))[None, :]
    mu = mu / (p * np.exp(eff)).sum()                          # rows sum ~ lib
    X = rng.poisson(mu).astype(np.int64)
    if ties:                                                   # hostile: many exact zeros / ties on carriers
        X[:n_car][rng.random((n_car, G)) < 0.3] = 0
    if zero_sub_top and keep is not None and n_car > 0:         # hostile: the deepest-CPM carrier holds no sub-gate counts
        X[0, ~keep] = 0
    X[X.sum(1) == 0, 0] = 1
    return X, n_car


def u_change(X0, X1, cc, keep):
    """Independent per-gene net 2U change (scipy Mann-Whitney U1, midrank ties) on the genes whose column changed."""
    s0, s1 = X0.sum(1) / 1e6, X1.sum(1) / 1e6
    gidx = np.flatnonzero(keep)
    out = {}
    for j in np.flatnonzero(np.any(X0[:, gidx] != X1[:, gidx], axis=0)):
        c = cc.cs[j]
        u0 = mannwhitneyu(np.log1p(X0[:, gidx[j]] / s0), c, method="asymptotic").statistic
        u1 = mannwhitneyu(np.log1p(X1[:, gidx[j]] / s1), c, method="asymptotic").statistic
        out[int(j)] = 2.0 * (u1 - u0)
    return out


def check(M, X0, X1, diag, n_car, cc, tie, own_col, tag, fails):
    keep = np.asarray(cc.keep, bool)
    gidx, sidx = np.flatnonzero(keep), np.flatnonzero(~keep)
    n, n_ctrl = X0.shape[0], cc.cs.shape[1]
    bad = []
    # P1
    if X1.dtype.kind not in "iu" or X1.min() < 0:
        bad.append("P1 dtype/negative")
    if not np.array_equal(X1.sum(1), X0.sum(1)):
        bad.append("P1 library moved")
    if not np.array_equal(X1.sum(0), X0.sum(0)):
        bad.append("P1 gene total moved")
    # P2
    capU = int(np.floor(1e-4 * 2.0 * n * n_ctrl))
    du = u_change(X0, X1, cc, keep)
    worst = max([abs(v) for v in du.values()], default=0.0)
    if worst > capU + 1e-9:
        bad.append(f"P2 |d2U| {worst} > capU {capU}")
    # P3 / P4
    p0, q0, l0 = M.de_stats(X0, cc, tie)
    p1, q1, l1 = M.de_stats(X1, cc, tie)
    if np.any((q0 < 0.05) != (q1 < 0.05)):
        bad.append(f"P3 call flips {int(np.sum((q0 < 0.05) != (q1 < 0.05)))}")
    s = X0.sum(1) / 1e6
    A0 = (X0[:, gidx] / s[:, None]).mean(0)
    A1 = (X1[:, gidx] / s[:, None]).mean(0)
    if np.any(A1 < A0 - 1e-12 * np.maximum(A0, 1)):
        bad.append(f"P4 mean fell on {int(np.sum(A1 < A0 - 1e-12 * np.maximum(A0, 1)))} genes")
    fin = np.isfinite(l0) & np.isfinite(l1)
    if np.any(np.sign(l0[fin]) != np.sign(l1[fin])):
        bad.append(f"P4 lfc sign changed on {int(np.sum(np.sign(l0[fin]) != np.sign(l1[fin])))} genes")
    frozen = (q0 < 0.05) | ~np.isfinite(l0) | (l0 <= 0)
    if own_col >= 0:
        frozen[own_col] = True
    ch = np.any(X0[:, gidx] != X1[:, gidx], axis=0)
    if np.any(ch & frozen):
        bad.append(f"P4 a frozen gene moved ({int(np.sum(ch & frozen))})")
    # P5
    d = X1 - X0
    if np.any(d[:n_car][:, sidx] > 0) or np.any(d[n_car:][:, sidx] < 0):
        bad.append("P5 sub-gate direction")
    if np.any(d[:n_car][:, gidx] < 0) or np.any(d[n_car:][:, gidx] > 0):
        bad.append("P5 kept-gene direction")
    if bad:
        fails.append((tag, bad))
    return dict(moved=int(ch.sum()), worst_d2U=float(worst), capU=capU, dmean_sum=float((A1 - A0).sum()))


def main() -> None:
    t0 = time.time()
    with tempfile.TemporaryDirectory(prefix="rev630-") as td:
        td = Path(td)
        git = lambda rev: subprocess.run(["git", "-C", str(ROOT), "show", f"{rev}:{SRC}"], check=True,
                                         capture_output=True, text=True).stdout
        V = {"v2": load_version("v2", git("65ee4584"), td), "v3": load_version("v3", git("ae39c950"), td),
             "cur": load_version("cur", (ROOT / SRC).read_text(), td)}
        print({k: m._SHA for k, m in V.items()}, flush=True)
        assert V["v2"]._SHA.startswith("5ae935b9") and V["v3"]._SHA.startswith("8e1f05a6")
        cur = V["cur"]
        rng = np.random.default_rng(630_7772)
        fails, rec = [], {"normal": [], "hostile": [], "equiv": {}, "hostile_old": []}
        n_runs = 0
        eq = {"v2_b_k100": [0, 0], "v2_s_k100": [0, 0], "v3_all": [0, 0]}
        for bi in range(10 if RICH else 24):
            C, p, libc = controls(rng, **(dict(G=2400, n_ctrl=1500, frac_sub=0.67) if RICH else {}))
            cc = ZS.scorer_control_block(C)
            keep = np.asarray(cc.keep, bool)
            tie = cur.scorer_tie(cc)
            lib_med = float(np.median(libc))
            kpos = np.cumsum(keep) - 1
            own_full = int(rng.choice(np.flatnonzero(keep)))
            own_col = int(kpos[own_full])
            X0, n_car = block(rng, p, libc, own=own_full, **(dict(n=400, n_bal=10) if RICH else {}))
            for v in VARIANTS:
                st0 = np.random.get_state()
                X1, d1 = cur.apply(X0.copy(), n_car, cc, v, tie=tie, lib_med=lib_med, own_col=own_col)
                X1b, d1b = cur.apply(X0.copy(), n_car, cc, v, tie=tie, lib_med=lib_med, own_col=own_col)
                st1 = np.random.get_state()
                if not (np.array_equal(X1, X1b) and d1 == d1b):
                    fails.append((f"b{bi} {v}", ["P6 nondeterministic"]))
                if not (st0[0] == st1[0] and np.array_equal(st0[1], st1[1]) and st0[2:] == st1[2:]):
                    fails.append((f"b{bi} {v}", ["P6 RNG consumed"]))
                r = check(cur, X0, X1, d1, n_car, cc, tie, own_col, f"b{bi} {v}", fails)
                later = 0
                m3 = cur.VARIANT_RE.fullmatch(v)
                if m3 and m3.group(3) and int(m3.group(4)) > 1:
                    v1 = f"{v.split('_q')[0]}_q{m3.group(3)}r1p{m3.group(5)}"
                    X11, _ = cur.apply(X0.copy(), n_car, cc, v1, tie=tie, lib_med=lib_med, own_col=own_col)
                    later = int(np.any(X11 != X1, axis=0).sum())
                rec["normal"].append(dict(block=bi, variant=v, cols_changed_by_later_rounds=later, **r, guard_rounds=d1.get("ur_guard_rounds"),
                                          guard_failed=d1.get("ur_guard_failed"), n_neg_donor_rows=d1.get("ur_n_neg_donor_rows"),
                                          n_mean_revert=d1.get("ur_n_mean_revert")))
                n_runs += 1
                # P7 equivalence
                if "_q" not in v:
                    for old in ("v2",):
                        if v in ("b_k100", "s_k100"):
                            Xo, _ = V[old].apply(X0.copy(), n_car, cc, v, tie=tie, lib_med=lib_med, own_col=own_col)
                            k = f"{old}_{v}"
                            eq[k][0] += int(np.array_equal(Xo, X1)); eq[k][1] += 1
                Xo3, _ = V["v3"].apply(X0.copy(), n_car, cc, v, tie=tie, lib_med=lib_med, own_col=own_col)
                eq["v3_all"][0] += int(np.array_equal(Xo3, X1)); eq["v3_all"][1] += 1
            # P8 hostile blocks
            for kind, kw in (("shallow3", dict(shallow=3)), ("shallow10", dict(shallow=10)), ("nobal", dict(n_bal=0)),
                             ("zero_sub_top", dict(zero_sub_top=True, keep=keep)), ("ties", dict(ties=True))):
                Xh, nch = block(rng, p, libc, own=own_full, **({**dict(n=400, n_bal=10), **kw} if RICH else kw))
                for v in ("b_k100", "s_k100", "b_k100_q20r3p10", "b_k300_q50r2p20"):
                    try:
                        Xh1, dh = cur.apply(Xh.copy(), nch, cc, v, tie=tie, lib_med=lib_med, own_col=own_col)
                        r = check(cur, Xh, Xh1, dh, nch, cc, tie, own_col, f"b{bi} {kind} {v}", fails)
                        rec["hostile"].append(dict(block=bi, kind=kind, variant=v, **r, n_no_pos_donor=dh.get("ur_n_no_pos_donor"),
                                                   n_mean_revert=dh.get("ur_n_mean_revert")))
                    except Exception as e:                          # noqa: BLE001
                        fails.append((f"b{bi} {kind} {v}", [f"P8 cur raised {type(e).__name__}: {e}"]))
                    for old in ("v2", "v3"):
                        if old == "v2" and "_q" in v:
                            continue
                        try:
                            Xo, _ = V[old].apply(Xh.copy(), nch, cc, v, tie=tie, lib_med=lib_med, own_col=own_col)
                            ofails = []
                            check(cur, Xh, Xo, {}, nch, cc, tie, own_col, "x", ofails)
                            rec["hostile_old"].append(dict(block=bi, kind=kind, variant=v, version=old,
                                                           issues=ofails[0][1] if ofails else []))
                        except Exception as e:                      # noqa: BLE001
                            rec["hostile_old"].append(dict(block=bi, kind=kind, variant=v, version=old,
                                                           issues=[f"raised {type(e).__name__}: {str(e)[:80]}"]))
            if bi == 0:
                # P9 fail-safe guard: force a flip on every re-check
                real = cur.de_stats
                calls = {"n": 0}

                def fake(x, c, t):
                    p_, q_, l_ = real(x, c, t)
                    calls["n"] += 1
                    if calls["n"] > 1:                              # every post-pass check reports gene 0's call flipped
                        q_ = q_.copy(); q_[0] = 0.01 if q_[0] >= 0.05 else 0.5
                    return p_, q_, l_
                cur.de_stats = fake
                try:
                    Xf, df = cur.apply(X0.copy(), n_car, cc, "b_k100", tie=tie, lib_med=lib_med, own_col=own_col)
                finally:
                    cur.de_stats = real
                Xn, dn = cur.apply(X0.copy(), n_car, cc, "b_k100", tie=tie, lib_med=lib_med, own_col=own_col)
                p9 = dict(returned_input=bool(np.array_equal(Xf, X0)), guard_failed=df.get("ur_guard_failed"),
                          guard_rounds=df.get("ur_guard_rounds"), passes=calls["n"] - 1,
                          same_keys=sorted(df) == sorted(dn), missing=sorted(set(dn) - set(df)), extra=sorted(set(df) - set(dn)))
                rec["P9"] = p9
                print("P9", p9, flush=True)
                if not (p9["returned_input"] and p9["guard_failed"] == 1 and p9["same_keys"]):
                    fails.append(("P9", [str(p9)]))
            print(f"block {bi}: {n_runs} cur runs, fails so far {len(fails)} ({time.time() - t0:.0f}s)", flush=True)
        rec["equiv"] = {k: f"{a}/{b} bit-identical" for k, (a, b) in eq.items()}
        print("P7 equivalence:", rec["equiv"])
        ho = rec["hostile_old"]
        summ_old = {}
        for o in ho:
            key = f"{o['version']} {o['kind']}"
            summ_old.setdefault(key, {"runs": 0, "with_issue": 0, "examples": []})
            summ_old[key]["runs"] += 1
            if o["issues"]:
                summ_old[key]["with_issue"] += 1
                if len(summ_old[key]["examples"]) < 2:
                    summ_old[key]["examples"].append(o["issues"])
        rec["hostile_old_summary"] = summ_old
        print("P8 old versions on hostile blocks:")
        for k, v in sorted(summ_old.items()):
            print(f"  {k:22s} {v['with_issue']}/{v['runs']} with an issue; e.g. {v['examples'][:1]}")
        hn = rec["hostile"]
        print(f"P8 cur on hostile blocks: {len(hn)} runs; moved genes median {np.median([h['moved'] for h in hn]):.0f}; "
              f"n_no_pos_donor > 0 in {sum((h['n_no_pos_donor'] or 0) > 0 for h in hn)}")
        nr = rec["normal"]
        lr = [r["cols_changed_by_later_rounds"] for r in nr if "_q" in r["variant"]]
        print(f"v3 later rounds changed columns in {sum(x > 0 for x in lr)}/{len(lr)} suffix runs (median {np.median(lr):.0f} columns)")
        print(f"normal: {len(nr)} runs; moved median {np.median([r['moved'] for r in nr]):.0f}; worst |d2U| / capU "
              f"{max(r['worst_d2U'] / r['capU'] for r in nr):.3f}; guard rounds max {max(r['guard_rounds'] for r in nr):.0f}; "
              f"neg-donor rows {sum(r['n_neg_donor_rows'] for r in nr):.0f}; mean reverts {sum(r['n_mean_revert'] for r in nr):.0f}")
        rec["fails"] = [(t, b) for t, b in fails]
        rec["module_shas"] = {k: m._SHA for k, m in V.items()}
        print(f"\nFAILS: {len(fails)}")
        for t, b in fails[:20]:
            print(" ", t, b)
        (HERE / "work").mkdir(exist_ok=True)
        (HERE / f"work/review630_7772dc{'_rich' if RICH else ''}.json").write_text(json.dumps(rec, indent=1, default=float))
        print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
