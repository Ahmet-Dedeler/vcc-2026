import json, sys
from collections import Counter
d = json.load(open('board/leaderboard_2026-09-03.json'))
E = d['entries'] if isinstance(d, dict) else d
print("n entries", len(E), "keys:", list(d.keys()) if isinstance(d, dict) else None)
E = sorted(E, key=lambda e: -(e.get('score_avg') or -9))
def row(i,e):
    return f"{i:>3} | {e['team_name'][:22]:<22} | {str(e.get('organization'))[:18]:<18} | {e['score_avg']:.5f} | {e['pds_cosine']:.4f} | {e['expr_mse_unbiased_capped_norm']:.3f} | {e['de_wilcoxon_lfc_nmae']:.4f} | {e['de_wilcoxon_direction_fidelity_yield_raw']:.4f} | {e['de_wilcoxon_direction_reach_raw']:.4f} | {e['de_wilcoxon_sig_jaccard']:.4f} | {e['submission_count']:>3} | {e['submission_date'][:16]} | {str(e.get('model_name'))[:40]}"
print("rank | team | org | overall | pds | mse | nmae | fid | reach | jac | n_sub | date | model")
for i,e in enumerate(E[:40],1): print(row(i,e))
print()
for i,e in enumerate(E,1):
    if e['team_name'] in ('Ahmet Dedeler',): print("OURS:", row(i,e))
targets = ['luka','cqawesome','The_Way','UM','Mengran Li','Weijing Zhao','pop','AgingLab_BI','Andrew Liu','Hongjian Zhou','jw_13288','GZ YAI','huhuhu','Leo','Aginglab.com','Slurm Shadi','Garrrbage']
print("\nTARGETS")
for i,e in enumerate(E,1):
    if e['team_name'] in targets: print(row(i,e)); print("    members:", e.get('team_members'), "| desc:", (e.get('description') or '')[:500].replace('\n',' '))
print("\nREACH>=0.23 & pds>=0.80")
for i,e in enumerate(E,1):
    if e['de_wilcoxon_direction_reach_raw']>=0.23 and e['pds_cosine']>=0.80: print(row(i,e))
print("\nPositive scaled fid:")
for i,e in enumerate(E,1):
    if (e.get('score_fid') or -1)>0: print(row(i,e), f"score_fid={e['score_fid']:.4f}")
