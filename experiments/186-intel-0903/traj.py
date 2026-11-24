import json, glob, os
files = sorted(glob.glob('../../data/leaderboard/*.json')) + ['board/leaderboard_2026-09-03.json']
teams = ['cqawesome','Aginglab.com','Andrew Liu','Hongjian Zhou','jw_13288','Garrrbage','Slurm Shadi','luka','The_Way','pop','UM','Mengran Li','Weijing Zhao','AgingLab_BI','GZ YAI','huhuhu','Leo','Haystack','Ahmet Dedeler','OpenPert','Vivai','Happy']
hist = {t:[] for t in teams}
for f in files:
    d = json.load(open(f)); E = d['entries'] if isinstance(d, dict) else d
    for e in E:
        if e['team_name'] in hist:
            hist[e['team_name']].append((os.path.basename(f)[:20], e['submission_date'][:16], e['submission_count'], e['score_avg'], e['pds_cosine'], e['expr_mse_unbiased_capped_norm'], e['de_wilcoxon_lfc_nmae'], e['de_wilcoxon_direction_fidelity_yield_raw'], e['de_wilcoxon_direction_reach_raw'], e['de_wilcoxon_sig_jaccard'], e.get('model_name'), (e.get('description') or '')[:200].replace('\n',' ')))
with open('board/trajectories.md','w') as out:
    for t in teams:
        out.write(f"\n## {t}\n")
        seen=set()
        print(f"\n## {t}")
        for h in hist[t]:
            key=(h[1],h[2])
            if key in seen: continue
            seen.add(key)
            line=f"  {h[1]} n={h[2]:>2} score {h[3]:.4f} pds {h[4]:.4f} mse {h[5]:.3f} nmae {h[6]:.4f} fid {h[7]:.4f} reach {h[8]:.4f} jac {h[9]:.4f} `{h[10]}`"
            print(line); out.write(line+"\n")
            if h[11]: print("      desc:", h[11]); out.write("      desc: "+h[11]+"\n")
