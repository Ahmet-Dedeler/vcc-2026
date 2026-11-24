"""301 item 3: resolve JoungZhang2023_atlas obs.perturbation (bare codes 0..3368) to TF ORF ids
and gene symbols. Anchors come from GEO GSE217460 per-cell TFmap files joined to the atlas
barcodes; the ORF id -> symbol table is Joung 2023 Table S1A (mmc2)."""
import gzip, csv, numpy as np, collections, json, re, openpyxl
S = '/private/tmp/claude-501/-Users-ahmet-Code-virtual-cell-challenge-2026/42edd049-4f44-4b3c-b669-9213db727327/scratchpad/'
r3 = np.load(S + 'joung_r3.npy'); p1 = np.load(S + 'joung_p1.npy'); codes = np.load(S + 'joung_codes.npy')
wb = openpyxl.load_workbook('joung2023_mmc2.xlsx', read_only=True)
rows = list(wb['Table S1A'].iter_rows(values_only=True)); ix = {h: i for i, h in enumerate(rows[0])}
bc2orf, orf2sym = {}, {}
for r in rows[1:]:
    nm = str(r[ix['Name']]).strip(); sym = str(r[ix['RefSeq Gene Name']]).strip()
    iso = str(r[ix['Isoform ID']]).strip(); rs = str(r[ix['RefSeq and Gencode ID']]).strip()
    orf2sym[nm] = (sym, iso, rs)
    b = r[ix['Barcode Sequence']]
    if b: bc2orf[str(b).strip()] = nm
files = {'S%02d' % i: 'GSM67199%02d' % (69 + i) for i in range(1, 11)}
p1groups = sorted(set(p1.tolist())); idx_by_p1 = {p: np.where(p1 == p)[0] for p in p1groups}
agg = collections.defaultdict(collections.Counter); assign = {}
for s, gsm in sorted(files.items()):
    tf = {}
    with gzip.open(f'{gsm}_210322_TFmap_{s}.csv.gz', 'rt') as fh:
        for r in csv.reader(fh): tf[r[0]] = r[1]
    scores = []
    for p in p1groups:
        by = collections.defaultdict(collections.Counter); tot = 0
        for i in idx_by_p1[p]:
            b = tf.get(r3[i])
            if b is not None: by[int(codes[i])][b] += 1; tot += 1
        if tot: scores.append((sum(c.most_common(1)[0][1] for c in by.values()) / tot, p, tot))
    scores.sort(reverse=True); pure, p, tot = scores[0]
    assign[s] = {'P1': p, 'purity': round(pure, 4), 'runner_up': [scores[1][1], round(scores[1][0], 4)]}
    print(f'{s}: P1 {p} purity {pure:.4f}  (runner-up {scores[1][1]} {scores[1][0]:.4f})', flush=True)
    for i in idx_by_p1[p]:
        b = tf.get(r3[i])
        if b in bc2orf: agg[int(codes[i])][bc2orf[b]] += 1
anchors = {c: cn.most_common(1)[0][0] for c, cn in agg.items()
           if sum(cn.values()) >= 10 and cn.most_common(1)[0][1] / sum(cn.values()) >= 0.8}
segs = sorted((c, int(re.sub(r'\D', '', o))) for c, o in anchors.items())
mono = all(segs[i][1] < segs[i + 1][1] for i in range(len(segs) - 1))
print('anchors', len(anchors), 'strictly monotone in code order?', mono, flush=True)
assert mono
allnums = list(range(1, 3551)); resolved = {c: n for c, n in segs}
c0, n0 = segs[0]; head = [x for x in allnums if x < n0]
if len(head) == c0:
    for j in range(c0): resolved[j] = head[j]
for i in range(len(segs) - 1):
    c0, n0 = segs[i]; c1, n1 = segs[i + 1]
    bt = [x for x in allnums if n0 < x < n1]
    if len(bt) == c1 - c0 - 1:
        for j, x in enumerate(bt): resolved[c0 + 1 + j] = x
c1, n1 = segs[-1]; tail = [x for x in allnums if x > n1]
if len(tail) == 3368 - c1:
    for j, x in enumerate(tail): resolved[c1 + 1 + j] = x
print('resolved', len(resolved), '/3369  ambiguous', 3369 - len(resolved), flush=True)
cnt = collections.Counter(codes.tolist())
with open('joung_orf_map.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['perturbation_code', 'n_cells', 'tforf_id', 'gene_symbol', 'isoform_id',
                'refseq_gencode_id', 'evidence'])
    for c in range(3369):
        n = resolved.get(c)
        if n is None:
            w.writerow([c, cnt.get(c, 0), '', '', '', '', 'AMBIGUOUS']); continue
        orf = 'TFORF%04d' % n
        sym, iso, rs = orf2sym.get(orf, ('', '', ''))
        w.writerow([c, cnt.get(c, 0), orf, sym, iso, rs,
                    'ANCHOR' if c in anchors else 'INTERPOLATED'])
json.dump({'sublibrary_to_P1': assign, 'n_anchors': len(anchors), 'monotone': mono,
           'n_resolved': len(resolved), 'n_ambiguous': 3369 - len(resolved)},
          open('joung_map_meta.json', 'w'), indent=1)
print('WROTE joung_orf_map.csv', flush=True)
