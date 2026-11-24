"""Known-answer tests for the reliability gate. No external data."""

from __future__ import annotations

import numpy as np

from gate import (
    FEATURE_NAMES,
    N_FEATURES,
    aggregate_sources,
    energy_weighted_cosine,
    fit_theta,
    fixed_weight,
    gate_weight,
    prepare_task,
    score_weight,
)


def test_fixed_formula_hand_value():
    # z^2 = 8, n = 4 -> 0.15 + 0.85 * (8/16) * (4/8) = 0.3625
    w = fixed_weight(np.array([8.0]), np.array([4.0]))
    assert np.allclose(w, 0.3625), w
    w0 = fixed_weight(np.array([0.0]), np.array([100.0]))
    assert np.allclose(w0, 0.15), w0
    winf = fixed_weight(np.array([1e12]), np.array([1e12]))
    assert np.allclose(winf, 1.0, atol=1e-6), winf
    grid = fixed_weight(np.array([0.0, 1.0, 8.0, 100.0]), np.array([1.0, 2.0, 4.0, 10.0]))
    assert np.all((grid >= 0.15 - 1e-9) & (grid <= 1.0 + 1e-9))


def test_zero_theta_matches_equal_weight():
    rng = np.random.default_rng(0)
    P, G, S = 25, 40, 3
    delta = rng.normal(size=(S, P, G))
    var = np.full((S, P, G), 0.05)
    n_cells = np.full((S, P), 80.0)
    ctrl = np.zeros((S, G))
    truth = rng.normal(size=(P, G))
    own = np.full(P, -1)
    ok = np.ones(P, dtype=bool)
    task = prepare_task(delta, var, n_cells, ctrl, truth, own, ok)
    w = gate_weight(np.zeros(N_FEATURES + 1), task["Z"])
    assert np.allclose(w, 0.15 + 0.85 * 0.5)
    scored = score_weight(task, w)
    base = score_weight(task, np.ones_like(task["mu"]))
    assert np.allclose(scored["cos"], base["cos"])
    assert np.isclose(scored["ew"], base["ew"])
    assert np.isclose(scored["uw"], base["uw"])


def _synthetic_bank(seed: int = 1):
    """Signal genes agree across lines. A second block is precisely measured but
    line-private (high z^2, high MAD). A third block is sampling noise.
    """
    rng = np.random.default_rng(seed)
    lines = ["A", "B", "C", "D"]
    P, G = 36, 30
    signal = rng.normal(scale=1.0, size=(P, 10))
    bank = {}
    for i, name in enumerate(lines):
        delta = np.zeros((P, G))
        var = np.empty((P, G))
        delta[:, :10] = signal + rng.normal(scale=0.05, size=(P, 10))
        var[:, :10] = 0.01
        # Precise but not shared: large effects, tiny sampling variance, independent per line.
        delta[:, 10:20] = rng.normal(scale=1.5, size=(P, 10))
        var[:, 10:20] = 0.005
        # Sampling noise: small effects, huge variance.
        delta[:, 20:] = rng.normal(scale=0.3, size=(P, 10))
        var[:, 20:] = 4.0
        bank[name] = {
            "delta": delta,
            "var": var,
            "n_cells": np.full(P, 100.0),
            "ctrl": np.linspace(-1, 2, G),
            "truth": None,
        }
    # Truth of a line is its own delta (the response we would observe).
    for name in lines:
        bank[name]["truth"] = bank[name]["delta"].copy()
    return lines, bank


def _task_from_bank(bank, sources, heldout, target_ok=None):
    delta = np.stack([bank[s]["delta"] for s in sources], axis=0)
    var = np.stack([bank[s]["var"] for s in sources], axis=0)
    n_cells = np.stack([bank[s]["n_cells"] for s in sources], axis=0)
    ctrl = np.stack([bank[s]["ctrl"] for s in sources], axis=0)
    truth = bank[heldout]["truth"]
    P = truth.shape[0]
    if target_ok is None:
        target_ok = np.ones(P, dtype=bool)
    own = np.full(P, -1)  # synthetic genes are not named targets
    return prepare_task(delta, var, n_cells, ctrl, truth, own, target_ok)


def test_learned_beats_fixed_and_equal_on_synthetic():
    lines, bank = _synthetic_bank(2)
    # Fit on A,B,C (leave D out of the fit). Inner LOO among A,B,C.
    train = ["A", "B", "C"]
    tasks = []
    for h in train:
        sources = [s for s in train if s != h]
        tasks.append(_task_from_bank(bank, sources, h))
    fit = fit_theta(tasks, seed=486, n_restarts=3, maxiter=400)
    # Outer: predict D from A,B,C.
    outer = _task_from_bank(bank, train, "D")
    w_learn = gate_weight(fit["theta"], outer["Z"])
    w_fixed = fixed_weight(outer["z2"], outer["n_src"])
    w_equal = np.ones_like(outer["mu"])
    ew_l = score_weight(outer, w_learn)["ew"]
    ew_f = score_weight(outer, w_fixed)["ew"]
    ew_e = score_weight(outer, w_equal)["ew"]
    # The discordant block is invisible to z^2 (sampling var is tiny) and visible to MAD.
    assert ew_l > ew_e + 0.02, (ew_l, ew_e, fit["coef"])
    assert ew_l > ew_f + 0.01, (ew_l, ew_f, fit["coef"])
    # Feature permutation across genes removes the alignment the fit used.
    rng = np.random.default_rng(486)
    perm = rng.permutation(outer["Z"].shape[1])
    Zp = outer["Z"][:, perm, :]
    w_perm = gate_weight(fit["theta"], Zp)
    ew_p = score_weight(outer, w_perm)["ew"]
    assert ew_p < ew_l - 0.01, (ew_p, ew_l)
    # Target shuffle collapses cosine.
    pt = rng.permutation(outer["mu"].shape[0])
    shuffled = dict(outer)
    shuffled["mu"] = outer["mu"][pt]
    shuffled["Z"] = outer["Z"][pt]
    shuffled["z2"] = outer["z2"][pt]
    shuffled["n_src"] = outer["n_src"][pt]
    # truth stays put; cos_mask rows must follow the prediction's targets.
    shuffled["cos_mask"] = outer["cos_mask"]  # mask is all-true in this synthetic
    ew_s = score_weight(shuffled, gate_weight(fit["theta"], shuffled["Z"]))["ew"]
    assert abs(ew_s) < 0.15, ew_s


def test_fit_does_not_read_heldout_response():
    lines, bank = _synthetic_bank(3)
    train = ["A", "B", "C"]

    def fit_with(held_delta_scale: float):
        bank2 = {k: dict(v) for k, v in bank.items()}
        bank2["D"] = dict(bank["D"])
        bank2["D"]["truth"] = bank["D"]["truth"] * held_delta_scale + 3.0
        bank2["D"]["delta"] = bank2["D"]["truth"]
        tasks = []
        for h in train:
            sources = [s for s in train if s != h]
            tasks.append(_task_from_bank(bank2, sources, h))
        return fit_theta(tasks, seed=486, n_restarts=1, maxiter=200)["theta"]

    t1 = fit_with(1.0)
    t2 = fit_with(-5.0)
    assert np.allclose(t1, t2, atol=1e-8), (t1, t2)


def test_aggregate_z2_and_mad_directions():
    # Two sources, one gene, one target. Agreeing precise effects -> large z^2, small MAD.
    delta = np.array([[[1.0]], [[1.1]]])
    var = np.array([[[0.01]], [[0.01]]])
    n_cells = np.array([[50.0], [50.0]])
    ctrl = np.array([[0.0]])
    agg = aggregate_sources(delta, var, n_cells, ctrl)
    assert agg["n_src"][0, 0] == 2
    assert agg["z2"][0, 0] > 50  # mu~1.05, v_mu = 0.01/2 = 0.005, z2 ~ 220
    assert agg["mu"][0, 0] == np.mean([1.0, 1.1])
    # Discordant.
    delta2 = np.array([[[1.0]], [[-1.0]]])
    agg2 = aggregate_sources(delta2, var, n_cells, ctrl)
    assert abs(agg2["mu"][0, 0]) < 1e-9
    assert agg2["raw"][0, 0, FEATURE_NAMES.index("log_source_mad")] > agg["raw"][
        0, 0, FEATURE_NAMES.index("log_source_mad")
    ]


def test_energy_weight_prefers_large_responses():
    cos = np.array([1.0, 0.0])
    energy = np.array([10.0, 1.0])
    ok = np.array([True, True])
    assert np.isclose(energy_weighted_cosine(cos, energy, ok), 10.0 / 11.0)


def test_merge_inverse_variance():
    from build_replogle import merge_screens

    a = {
        "genes": np.array(["G1", "G2"]),
        "perts": np.array(["T"]),
        "lfc": np.array([[1.0, 5.0]]),
        "var": np.array([[0.01, 0.02]]),
        "n_cells": np.array([10.0]),
        "n_guides": np.array([1]),
        "ctrl_log": np.array([0.0, 1.0]),
        "n_ctrl": 100.0,
    }
    b = {
        "genes": np.array(["G1", "G3"]),
        "perts": np.array(["T", "U"]),
        "lfc": np.array([[3.0, 7.0], [0.0, 1.0]]),
        "var": np.array([[0.03, 0.01], [0.1, 0.1]]),
        "n_cells": np.array([30.0, 5.0]),
        "n_guides": np.array([2, 1]),
        "ctrl_log": np.array([0.0, 3.0]),
        "n_ctrl": 300.0,
    }
    m = merge_screens(a, b)
    gi = list(m["genes"]).index("G1")
    pi = list(m["perts"]).index("T")
    # w = 1/0.01 and 1/0.03 -> 100 and 33.333...; weighted mean = 1.5
    assert np.isclose(m["lfc"][pi, gi], 1.5), m["lfc"][pi, gi]
    assert np.isclose(m["var"][pi, gi], 1.0 / (100.0 + 100.0 / 3.0)), m["var"][pi, gi]
    # Gene only in a.
    g2 = list(m["genes"]).index("G2")
    assert np.isclose(m["lfc"][pi, g2], 5.0)
    # Pert only in b.
    ui = list(m["perts"]).index("U")
    g3 = list(m["genes"]).index("G3")
    assert np.isclose(m["lfc"][ui, g3], 1.0)
    assert np.isclose(m["n_cells"][pi], 40.0)


if __name__ == "__main__":
    test_fixed_formula_hand_value()
    test_zero_theta_matches_equal_weight()
    test_aggregate_z2_and_mad_directions()
    test_energy_weight_prefers_large_responses()
    test_fit_does_not_read_heldout_response()
    test_learned_beats_fixed_and_equal_on_synthetic()
    test_merge_inverse_variance()
    print("ALL TESTS PASSED")
