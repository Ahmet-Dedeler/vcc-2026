#!/usr/bin/env bash
# Known-answer table for `scripts/aws_reap.sh`'s verdict classifier, the function extracted
# VERBATIM from the script. A reaper that has only ever seen an empty account is untested on
# the branch that matters, so every branch and both sides of every boundary get a row.
#
# Rows 1-2 are the instance this whole experiment exists because of: `i-0af19b9a537e08d4b`,
# r7i.2xlarge (8 vCPU), us-west-2, left RUNNING for 13 days after experiment 061 finished on
# 2026-08-22 at 0.27% average and 1.9% maximum CPU.
# Rows 7-8 are the failure mode the percentage-only version of this script had: a FULLY BUSY
# single-threaded `cell-eval2` scorer, measured at 13.7% on 8 vCPUs (1.10 busy cores) and
# therefore 3.4% on 32 vCPUs — under a 5% threshold. The cores reading keeps it alive.
verdict() { python3 -c "
import sys
age, cpu, max_age = float(sys.argv[1]), sys.argv[2], float(sys.argv[3])
cpu_max, vcpus, busy_cores = float(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
if age <= max_age:
    print('keep (young)')
elif cpu == 'nan':
    print('keep (no CPU data)')
elif float(cpu) < cpu_max and float(cpu) * vcpus / 100.0 < busy_cores:
    print('REAP')
else:
    print('keep (busy)')" "$1" "$2" "$3" "$4" "$5" "$6"; }

printf '%-8s %-8s %-6s %-8s %-6s %-6s  %-20s %-20s %s\n' \
  AGE_H CPU% VCPU MAX_AGE CPUMAX CORES GOT EXPECTED RESULT
fail=0
while read -r age cpu vcpu max_age cpu_max cores expected; do
  got="$(verdict "$age" "$cpu" "$max_age" "$cpu_max" "$vcpu" "$cores")"
  r=PASS; [[ "$got" == "$expected" ]] || { r=FAIL; fail=1; }
  printf '%-8s %-8s %-6s %-8s %-6s %-6s  %-20s %-20s %s\n' \
    "$age" "$cpu" "$vcpu" "$max_age" "$cpu_max" "$cores" "$got" "$expected" "$r"
done <<'ROWS'
312.0 0.27 8 24 5.0 0.25 REAP
312.0 1.9 8 24 5.0 0.25 REAP
2.0 0.0 8 24 5.0 0.25 keep (young)
0.0 nan 8 24 5.0 0.25 keep (young)
48.0 nan 8 24 5.0 0.25 keep (no CPU data)
48.0 87.0 8 24 5.0 0.25 keep (busy)
48.0 13.7 8 24 5.0 0.25 keep (busy)
48.0 3.4 32 24 5.0 0.25 keep (busy)
48.0 5.0 8 24 5.0 0.25 keep (busy)
48.0 4.999 8 24 5.0 0.25 keep (busy)
48.0 3.2 8 24 5.0 0.25 keep (busy)
48.0 3.13 8 24 5.0 0.25 keep (busy)
48.0 3.12 8 24 5.0 0.25 REAP
48.0 0.77 32 24 5.0 0.25 REAP
48.0 0.79 32 24 5.0 0.25 keep (busy)
24.0 0.1 8 24 5.0 0.25 keep (young)
24.1 0.1 8 24 5.0 0.25 REAP
ROWS
echo
[[ $fail -eq 0 ]] && echo "ALL ROWS PASS" || { echo "SOME ROWS FAIL"; exit 1; }
