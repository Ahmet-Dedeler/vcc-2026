import json, glob, os
new = json.load(open('board/leaderboard_2026-09-03.json'))['entries']
old = json.load(open('../../data/leaderboard/lb_20260828T1325.json'))
old = old['entries'] if isinstance(old, dict) else old
print("old n", len(old), "new n", len(new))
oldmap = {e['team_id']: e for e in old}
# all descriptions in new
out = []
for e in sorted(new, key=lambda e: -(e['score_avg'] or -9)):
    d = (e.get('description') or '').strip()
    m = (e.get('model_name') or '').strip()
    o = oldmap.get(e['team_id'])
    od = (o.get('description') or '').strip() if o else None
    om = (o.get('model_name') or '').strip() if o else None
    status = 'NEW_TEAM' if o is None else ('CHANGED' if (d!=od or m!=om) else 'same')
    out.append((e['score_avg'], e['team_name'], e.get('organization'), m, d, status, om, od, e['submission_count'], e['pds_cosine'], e['de_wilcoxon_direction_reach_raw'], e['expr_mse_unbiased_capped_norm']))
n_desc = sum(1 for r in out if r[4])
print("entries with description:", n_desc)
with open('board/descriptions_all_2026-09-03.md','w') as f:
    f.write("# All model_name + description fields, 2026-09-03 snapshot, sorted by score\n\n")
    for i,r in enumerate(out,1):
        f.write(f"## {i}. {r[1]} ({r[2]}) score={r[0]:.5f} pds={r[9]:.4f} reach={r[10]:.4f} mse={r[11]:.3f} n_sub={r[8]} [{r[5]}]\n**model:** `{r[3]}`\n")
        if r[4]: f.write(f"**desc:** {r[4]}\n")
        if r[5]=='CHANGED' and r[7] and r[7]!=r[4]: f.write(f"**prev desc (Aug28):** {r[7]}\n")
        if r[5]=='CHANGED' and r[6]!=r[3]: f.write(f"**prev model (Aug28):** `{r[6]}`\n")
        f.write("\n")
# Print top-80 changed/new with descriptions
print("\n=== TOP 120: changed or new, with any description text ===")
for i,r in enumerate(out[:120],1):
    if r[5]!='same' or r[4]:
        print(f"{i:>3} {r[1][:20]:<20} {r[0]:.4f} pds={r[9]:.3f} reach={r[10]:.3f} mse={r[11]:.2f} [{r[5]}] model=`{r[3][:50]}`")
        if r[4]: print("      DESC:", r[4][:700].replace('\n',' '))
        if r[5]=='CHANGED' and r[7] and r[7]!=r[4]: print("      PREV:", r[7][:300].replace('\n',' '))
