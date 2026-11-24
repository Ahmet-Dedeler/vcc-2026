"""162 — the realised board record, transcribed from its primary sources.

Every number below was read out of the file named in `src`, not out of a summary table.
`docs/CURRENT_STATE.md` and `docs/LOG.md` both carry secondary copies of several of these
and they agree; where they did not, the experiment's own BOARD-RESULT file wins.

Scaled members only. The board publishes `score_*` directly and their unweighted mean
reproduces `score_avg` to 3e-16, so `overall = mean(6 members)` is exact arithmetic here,
not an approximation. Deltas are therefore `sum(member deltas)/6` exactly, and that identity
is asserted in `check()` as a transcription check.
"""
from __future__ import annotations

MEMBERS = ("pds", "mse", "nmae", "fid", "reach", "jac")

# arm -> (overall, {member: scaled value})
ENTRIES: dict[str, dict] = {
    "v11": dict(
        overall=0.1501827, entry=None, date="2026-08-27",
        src="docs/CURRENT_STATE.md leaderboard section",
        members=None,  # member decomposition not recorded for v11
    ),
    "a2-prune": dict(
        overall=0.1503084, entry="o8rCQvdgjb4I8rSLHKk1", date="2026-08-27",
        src="experiments/094-reach-ordering/A8-BOARD-RESULT.md",
        members=dict(pds=0.736560993522496, mse=0.0, nmae=0.094339,
                     fid=-0.035087, reach=0.102951, jac=0.003086),
    ),
    "a8": dict(
        overall=0.1543607, entry="lfnx5c9vp1gVhlNWcemT", date="2026-08-27",
        src="experiments/094-reach-ordering/A8-BOARD-RESULT.md",
        members=dict(pds=0.736560993522496, mse=0.0, nmae=0.094301,
                     fid=-0.028694, reach=0.124349, jac=-0.000353),
    ),
    "a10": dict(
        overall=0.1541841, entry="lwHUEPPfuUAwvrfwrkne", date="2026-08-28",
        src="experiments/094-reach-ordering/A10-BOARD-RESULT.md",
        members=dict(pds=0.736560993522496, mse=0.0, nmae=0.094338,
                     fid=-0.027897, reach=0.123595, jac=-0.001493),
    ),
    "a11": dict(
        overall=0.16611503, entry="oIL3tgNUMocufeha4Kys", date="2026-08-28",
        src="experiments/096-mse-floor/A11-BOARD-RESULT.md",
        members=dict(pds=0.739889, mse=0.082041, nmae=0.096982,
                     fid=-0.031535, reach=0.112220, jac=-0.002907),
    ),
    "a12": dict(
        overall=0.16798741, entry=None, date="2026-08-29",
        src="experiments/096-mse-floor/A13-BOARD-RESULT.md (a12 column)",
        members=dict(pds=0.737885, mse=0.091163, nmae=0.093974,
                     fid=-0.031473, reach=0.116790, jac=-0.000413),
    ),
    "a13": dict(
        overall=0.1476855, entry="RnN67G6OWerJ1dTMANzQ", date="2026-08-29",
        src="experiments/096-mse-floor/A13-BOARD-RESULT.md",
        members=dict(pds=0.709106, mse=0.145924, nmae=-0.029332,
                     fid=-0.036758, reach=0.097121, jac=0.000051),
    ),
    "a14": dict(
        overall=0.1704242443, entry="0NvkHnMPzLFaJl8phAxS", date="2026-08-30",
        src="experiments/096-mse-floor/A14-BOARD-RESULT.md",
        members=dict(pds=0.738614, mse=0.110615, nmae=0.091972,
                     fid=-0.032586, reach=0.114542, jac=-0.000612),
    ),
    "a15": dict(
        overall=0.1740865087, entry="KpoZt5Rie0EjgMaZQECE", date="2026-08-30",
        src="experiments/096-mse-floor/A15-BOARD-RESULT.md",
        members=dict(pds=0.745871, mse=0.120057, nmae=0.094239,
                     fid=-0.030397, reach=0.115368, jac=-0.000619),
    ),
    "a19": dict(
        overall=0.160569, entry="QGpNVMwzzdSlhUoQtRdp", date="2026-08-31",
        src="docs/LOG.md 2026-08-31 a19 BOARD RESULT table",
        members=dict(pds=0.744688, mse=0.024149, nmae=0.102690,
                     fid=-0.026398, reach=0.118338, jac=-0.000053),
    ),
    "a21": dict(
        overall=0.1526293, entry="5Zvx3QDyb6nJMKXS0h5C", date="2026-09-01",
        src="docs/LOG.md 2026-09-01 a21 table",
        members=dict(pds=0.7458708290005253, mse=0.1200566151966586,
                     nmae=0.030025, fid=-0.059861, reach=0.096886, jac=-0.017201),
    ),
    "a22": dict(
        overall=0.1715689, entry="EFf5zr3H1PEbJW0imuQo", date="2026-09-01",
        src="docs/LOG.md 2026-09-01 a22 table (deltas vs a15)",
        members=dict(pds=0.7458708290005253, mse=0.1200566151966586,
                     nmae=0.094239 + 0.001520, fid=-0.030397 - 0.004447,
                     reach=0.115368 - 0.010656, jac=-0.000619 - 0.001523),
    ),
    "a23": dict(
        overall=0.1733597079221444, entry="IzQVfDpeDWufkPT8Lmd2", date="2026-09-02",
        src="experiments/156-delivered-geometry/A23-BOARD-RESULT.md + live board API",
        members=dict(pds=0.7476077488415546, mse=0.12158129527129787,
                     nmae=0.09413264797250775, fid=-0.030060474991219818,
                     reach=0.10743545998657263, jac=-0.0005384295478467913),
    ),
}

# arm -> (base arm, one-line config diff, channel)
# `channel` is what the change touches, which determines whether an effect-level
# validator can see it at all:
#   geometry  : changes the log-space effect matrix or the rendered profile
#   emission  : changes only the emitted counts (pds/mse invariant by construction)
#   mixed     : both, in one arm
ARMS: dict[str, dict] = {
    "a2-prune": dict(base="v11", diff="whitener_id union -> union_pruned",
                     channel="geometry"),
    "a8": dict(base="a2-prune", diff="overdispersion 0 -> 4.0", channel="emission"),
    "a10": dict(base="a8", diff="tail_spread .25->.15, overdispersion 4->16",
                channel="emission"),
    "a11": dict(base="a10", diff="alpha 1.0796 -> 0.838", channel="geometry"),
    "a12": dict(base="a11", diff="alpha 0.838 -> 1.0108, tier-8000 spread policy",
                channel="mixed"),
    "a13": dict(base="a12", diff="donor-agreement per-gene + per-pert shrinkage",
                channel="geometry"),
    "a14": dict(base="a12", diff="alpha 1.0108 -> 0.982", channel="geometry"),
    "a15": dict(base="a14", diff="cd4 fraction 0.05 -> 0.09473 (target f 0.10)",
                channel="geometry"),
    "a19": dict(base="a15", diff="whiten gamma .65->.80, lambda 0->0.007",
                channel="geometry"),
    "a21": dict(base="a15", diff="resample_mix emitter x=0.14 inverted",
                channel="emission"),
    "a22": dict(base="a15", diff="resample_mix emitter head-only x=0.10",
                channel="emission"),
    "a23": dict(base="a15", diff="mass_fix minnorm", channel="geometry"),
}


def delta(arm: str) -> dict:
    """Realised board delta of `arm` against its base, per member and overall."""
    spec = ARMS[arm]
    a, b = ENTRIES[arm], ENTRIES[spec["base"]]
    out = {"arm": arm, "base": spec["base"], "channel": spec["channel"],
           "diff": spec["diff"], "overall": a["overall"] - b["overall"]}
    if a["members"] and b["members"]:
        for m in MEMBERS:
            out[m] = a["members"][m] - b["members"][m]
    return out


def check() -> list[str]:
    """mean(member deltas) must reproduce the overall delta. Transcription check."""
    problems = []
    for arm in ARMS:
        d = delta(arm)
        if MEMBERS[0] not in d:
            continue
        implied = sum(d[m] for m in MEMBERS) / 6.0
        if abs(implied - d["overall"]) > 5e-6:
            problems.append(
                f"{arm}: member mean {implied:+.7f} != overall {d['overall']:+.7f} "
                f"(residual {implied - d['overall']:+.2e})"
            )
    return problems


if __name__ == "__main__":
    bad = check()
    hdr = f"{'arm':<10} {'base':<10} {'channel':<9} {'overall':>11} " + " ".join(
        f"{m:>10}" for m in MEMBERS)
    print(hdr)
    print("-" * len(hdr))
    for arm in ARMS:
        d = delta(arm)
        cells = " ".join(f"{d[m]:+10.6f}" if m in d else f"{'--':>10}" for m in MEMBERS)
        print(f"{arm:<10} {d['base']:<10} {d['channel']:<9} {d['overall']:+11.7f} {cells}")
    print()
    print("consistency check:", "PASS" if not bad else "FAIL")
    for p in bad:
        print("  ", p)
