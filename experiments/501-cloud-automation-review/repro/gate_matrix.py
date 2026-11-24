#!/usr/bin/env python3
"""Does each chain's awk gate match another long-lived bash argv?

Patterns are copied from the chain scripts. A MATCH means that process keeps
the gate false (the waiter stays in its until-loop).
"""
import re

GATES = {
    "chain479": r"chain471c\.sh|chain475b\.sh|run_arm_generic\.sh|def_day\.sh|run_ka_4",
    "chain495": r"chain479b?\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4",
    "chain479b": r"chain479\.sh|chain495\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4",
    "chain476": r"chain479b?\.sh|chain495\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4",
    "chain497": r"chain47[69]b?\.sh|chain495\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4",
    "chain498": r"chain47[69]b?\.sh|chain49[57]\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4",
    "chain496p": r"chain4(7[69]b?|9[5-8])\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4",
    "chain497e": r"chain4(7[69]b?|9[5-8])\.sh|chain496p\.sh|run_ka_4|run_arm_generic\.sh|def_day\.sh|run_dump4",
}
PROCS = {
    "chain479": "bash experiments/479-pertarget-sharpen/chain479.sh",
    "chain479b": "bash experiments/479-pertarget-sharpen/chain479b.sh",
    "relaunch479b": "bash experiments/479-pertarget-sharpen/relaunch479b.sh",
    "chain1002": "bash experiments/479-pertarget-sharpen/chain1002.sh",
    "arm1002": "bash experiments/479-pertarget-sharpen/arm1002.sh 1 479p2g log pkg name note",
    "chain495": "bash experiments/495-lodo-cap6/chain495.sh",
    "chain1003": "bash experiments/495-lodo-cap6/chain1003.sh",
    "arm1003": "bash experiments/495-lodo-cap6/arm1003.sh 2 471c log pkg name",
    "chain476": "bash experiments/476-def-day-merge/chain476.sh",
    "chain497": "bash experiments/497-october-switches/chain497.sh",
    "chain497e": "bash experiments/497-october-switches/chain497e.sh",
    "chain498": "bash experiments/498-combined-bulk/chain498.sh",
    "chain1004": "bash experiments/498-combined-bulk/chain1004.sh",
    "arm1004": "bash experiments/498-combined-bulk/arm1004.sh 1 498sp2 log pkg name",
    "chain496p": "bash experiments/496-uncalled-sign-hook/chain496p.sh",
    "run_ka_495": "bash experiments/495-lodo-cap6/run_ka_495.sh ka1",
    "run_dump495": "bash experiments/495-lodo-cap6/run_dump495.sh lodo05hi6",
    "run_arm": "bash experiments/428-zero-gene-audit/run_arm_generic.sh 479p2g",
    "def_day": "bash scripts/def_day.sh bundle work tag 1.75 auto pt 4",
    "commit479": "bash scripts/commit_after_chain.sh experiments/479-pertarget-sharpen/commit479.spec",
    "commit495": "bash scripts/commit_after_chain.sh experiments/495-lodo-cap6/commit495.spec",
    "commit476": "bash scripts/commit_after_chain.sh experiments/476-def-day-merge/commit476.spec",
    "commit497": "bash scripts/commit_after_chain.sh experiments/497-october-switches/commit497.spec",
    "commit496": "bash scripts/commit_after_chain.sh experiments/496-uncalled-sign-hook/commit496.spec",
}

def hit(line, pat):
    # $1 ~ /^bash$/ && /pat/  — line is "bash <args>"
    return line.split()[0] == "bash" and re.search(pat, line) is not None

print(f"{'proc':<16}" + "".join(f"{g:>10}" for g in GATES))
for pname, argv in PROCS.items():
    cells = []
    for g, pat in GATES.items():
        cells.append("MATCH" if hit(argv, pat) else ".")
    print(f"{pname:<16}" + "".join(f"{c:>10}" for c in cells))
