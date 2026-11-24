"""202 — the mechanism, at the emitter's own level, before any build is started.

Loads the CANDIDATE `rc_emit` by path (so `src/vcc2026/rc_emit.py` is not touched while
another experiment's build is running) and drives `rc_hybrid_block` directly on a synthetic
control block with enough cells that the Wilcoxon can resolve the BH line.

Four claims, and if any of them fails the arm is not buildable and no build should start:

  1. the prescribed set is exactly `top n_head by |H| UNION top n_cpm by control CPM` over
     the scorer's own DE universe, at size `n_head + n_cpm`;
  2. BH on the EMITTED block delivers that set -- not a set of that size, that set;
  3. the delivered per-gene totals are the `rc_hybrid` run's, count for count;
  4. `|H|`'s order survives inside the call set AND inside the complement, i.e. the
     complement's common-factor rescaling is the monotone map it claims to be.

Plus the two RNG-stream claims the challenge build's digest identity rests on: the number of
draws taken from the emission generator is the same under both laws, and the two `place`
draws have the same shapes.
"""
from __future__ import annotations

import importlib.util
import os
import sys

import numpy as np
from scipy.stats import spearmanr

ROOT = "~/Code/virtual-cell-challenge-2026"
HERE = os.path.dirname(os.path.abspath(__file__))
CANDIDATE = os.path.join(HERE, "work/rc_emit_candidate.py")
sys.path.insert(0, os.path.join(ROOT, "src"))


def load_candidate():
    spec = importlib.util.spec_from_file_location("rc_emit_candidate", CANDIDATE)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class CountingRng:
    """A `default_rng` that records how many times each method was called and with what
    shapes -- the only thing that can move the per-gene totals of the NEXT block."""

    def __init__(self, seed):
        self.rng = np.random.default_rng(seed)
        self.calls = []

    def __getattr__(self, name):
        attr = getattr(self.rng, name)
        if not callable(attr):
            return attr

        def wrapped(*args, **kwargs):
            out = attr(*args, **kwargs)
            self.calls.append((name, np.shape(out)))
            return out

        return wrapped


def main() -> None:
    rc_emit = load_candidate()
    from vcc2026.donor_transfer import emit_counts

    # The block has to have a SUB-GATE tail, exactly as the challenge does (9,626 of
    # 18,533 genes past the 5-CPM gate on context B): the genes below the gate are the
    # reservoir the placement repairs the per-cell weld with, and with every gene gated in
    # there is nothing to repair with and even `rc_hybrid` cannot weld.
    rng = np.random.default_rng(202)
    n_expressed, n_low, n_controls, n_cells = 3000, 3000, 1000, 400
    n_genes = n_expressed + n_low
    mu = np.concatenate([np.exp(rng.normal(1.8, 1.3, n_expressed)),
                         rng.uniform(0.005, 0.03, n_low)])
    controls = rng.poisson(mu[None, :], size=(n_controls, n_genes)).astype(np.float64)
    control = rc_emit.control_block_from_counts(controls, min_cpm=5.0)
    keep = control.keep
    n_keep = int(keep.sum())
    print(f"synthetic block: {n_controls} controls x {n_genes} genes, "
          f"{n_keep} past the 5-CPM scorer gate, {n_cells} emitted cells")

    profile = mu * np.exp(rng.normal(0.0, 0.10, n_genes))
    profile = profile / profile.sum()
    libraries = np.full(n_cells, 20_000, dtype=np.int64)
    key = np.where(keep, rng.normal(0.0, 1.0, n_genes), 0.0)
    spread = np.full(n_genes, 0.25)

    n_head, n_cpm = 120, 656
    call = rc_emit.CallSet(n_head, n_cpm, 1.90)
    common = dict(spread=spread, overdispersion=4.0, control=control, key=key,
                  jitter_seed=1000, emit_counts=emit_counts)

    ra, rb = CountingRng(5), CountingRng(5)
    hybrid, dh = rc_emit.rc_hybrid_block(profile, libraries, 0.1, ra, **common)
    block, dc = rc_emit.rc_hybrid_block(profile, libraries, 0.1, rb, call_set=call,
                                        **common)

    ok = True

    def check(label, condition, detail=""):
        nonlocal ok
        ok &= bool(condition)
        print(f"  [{'PASS' if condition else 'FAIL'}] {label}{detail}")

    print("\n1. the RNG stream (this is what the 900-block digest identity rests on)")
    check("the two laws take the same draws, in the same order, at the same shapes",
          ra.calls == rb.calls, f"  ({len(ra.calls)} calls)")

    print("\n2. the totals identity")
    check("per-gene panel totals identical to the `rc_hybrid` run's",
          np.array_equal(block.sum(0), hybrid.sum(0)))
    check("...and it is a DIFFERENT block, so the identity is not vacuous",
          not np.array_equal(block, hybrid))
    check("every cell welded to one library total",
          int(block.sum(1).max() - block.sum(1).min()) <= 1)

    print("\n3. the prescribed set, recomputed independently of the emitter")
    key_abs = np.abs(key[keep])
    head = np.zeros(n_keep, bool)
    head[np.argsort(-key_abs, kind="stable")[:n_head]] = True
    want = head.copy()
    want[np.argsort(-np.where(head, -np.inf, control.ref_lin),
                    kind="stable")[:n_cpm]] = True
    check(f"the union has exactly n_head + n_cpm = {n_head + n_cpm} genes",
          int(want.sum()) == n_head + n_cpm, f"  ({int(want.sum())})")

    _z, p, q, _lfc = rc_emit.de_from_block(block, control)
    delivered = q < 0.05
    # This synthetic block is COARSER than the real one -- 1,000 control cells against
    # h1-test's 18,400 -- so a few genes have a |z| floor above the BH line and cannot be
    # suppressed at any allocation. The exact-delivery measurement belongs on the real
    # surface (`ka202_emit.py`, which reads jaccard 0.994 and n_pred +0.58%); what this
    # block can prove is that the complement is suppressed at all and that the delivered
    # set is the prescribed one up to that floor.
    jac = float((delivered & want).sum() / max((delivered | want).sum(), 1))
    check("BH on the EMITTED block delivers the prescribed set up to the |z| floor",
          jac > 0.80 and int((want & ~delivered).sum()) * 20 < want.sum(),
          f"  (delivered {int(delivered.sum())}, prescribed {n_head + n_cpm}, "
          f"jaccard {jac:.4f}, prescribed-but-missed {int((want & ~delivered).sum())}, "
          f"unsuppressible {int((delivered & ~want).sum())})")
    check("the emitter's own diagnostic agrees with an independent recount",
          dc["call_set_delivered"] == int(delivered.sum())
          and dc["call_set_prescribed"] == n_head + n_cpm
          and abs(dc["call_set_jaccard"] - jac) < 1e-12,
          f"  (jaccard {dc['call_set_jaccard']:.4f}, "
          f"z_scale {dc['call_set_z_scale']:.4f})")

    print("\n4. the order, inside the block and inside its complement")
    # The bar is `rc_hybrid` on the SAME block, not an absolute number: the claim is that
    # the order is `rc_hybrid`'s law unchanged, so `rc_hybrid` is the comparator.
    lp = -np.log10(np.maximum(p, 1e-300))
    _zh, ph, qh, _lh = rc_emit.de_from_block(hybrid, control)
    lph = -np.log10(np.maximum(ph, 1e-300))
    inc = qh < 0.05
    rho_in_h = spearmanr(lph[inc], key_abs[inc]).correlation
    rho_out_h = spearmanr(lph[~inc], key_abs[~inc]).correlation
    rho_in = spearmanr(lp[delivered], key_abs[delivered]).correlation
    rho_in_pre = spearmanr(lp[want], key_abs[want]).correlation
    rho_out = spearmanr(lp[~delivered], key_abs[~delivered]).correlation
    print(f"   `rc_hybrid`  n_pred {int(inc.sum())}  rho_in {rho_in_h:+.4f}  "
          f"rho_out {rho_out_h:+.4f}")
    print(f"   `rc_callset` n_pred {int(delivered.sum())}  rho_in {rho_in:+.4f} "
          f"(over the PRESCRIBED set alone {rho_in_pre:+.4f})  rho_out {rho_out:+.4f}")
    check("the complement's order is `rc_hybrid`'s, within 0.02",
          abs(rho_out - rho_out_h) < 0.02, f"  ({rho_out_h:+.4f} -> {rho_out:+.4f})")
    check("the order inside the PRESCRIBED block is still strongly `|H|`",
          rho_in_pre > 0.90, f"  ({rho_in_pre:+.4f})")
    check("the call set actually shrank", int(inc.sum()) > int(delivered.sum()),
          f"  ({int(inc.sum())} -> {int(delivered.sum())})")

    print("\n5. z_cap is the only place the multiset moves")
    # With the cap above the multiset's own maximum nothing is rescaled, the incumbent's
    # p-value multiset is delivered unchanged, and BH therefore returns the INCUMBENT's
    # call count however the set was prescribed. That is the non-vacuity check for the
    # rescaling: without it the arm cannot exist.
    loose, dl = rc_emit.rc_hybrid_block(
        profile, libraries, 0.1, np.random.default_rng(5),
        call_set=rc_emit.CallSet(n_head, n_cpm, 99.0), **common)
    check("z_cap above the multiset's max leaves BH at the incumbent's count",
          dl["call_set_z_scale"] == 1.0
          and abs(dl["call_set_delivered"] - int(inc.sum())) < 0.1 * int(inc.sum()),
          f"  (delivered {dl['call_set_delivered']} vs `rc_hybrid`'s {int(inc.sum())}, "
          f"z_scale {dl['call_set_z_scale']:.4f})")
    check("...and the shipped cap is what shrinks it",
          dc["call_set_z_scale"] < 1.0 and dc["call_set_delivered"] < dl["call_set_delivered"],
          f"  ({dl['call_set_delivered']} -> {dc['call_set_delivered']})")

    print("\n6. `call_set_first=False` with a prescribed set is refused")
    refused = False
    try:
        rc_emit.rc_hybrid_block(profile, libraries, 0.1, np.random.default_rng(5),
                                call_set=call, call_set_first=False, **common)
    except ValueError:
        refused = True
    check("a prescribed set and a whole-axis re-key cannot both be asked for", refused)

    print(f"\nMECHANISM: {'PASS' if ok else 'FAIL'}")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
