"""301 item 2 (second, correct attempt): promoter windows on the TSS annotation the hCRISPRi-v2
library was actually designed against -- Horlbeck et al. 2016 Supplementary file 2, "TSS annotations
for hg19", the FANTOM-CAGE-derived P1 / P2 / P1P2 TSSs. Sequence from the Ensembl GRCh37 REST
endpoint. Restricted to the 18,080 gene symbols measured in the VCC 2025 h5ads.

The first attempt (`fetch_promoters.py`, GRCh38, Ensembl gene-level start/end as the TSS) FAILED its
known-answer row: only 22/68 library spacers matched their own target's window. The CAGE TSS is what
the library targets, and the gene-level Ensembl boundary is often far from it."""
import gzip, json, time, urllib.request, openpyxl
from concurrent.futures import ThreadPoolExecutor

FLANK = 500
symbols = {l.strip() for l in open('~/Code/virtual-cell-challenge-2026/data/external/'
                                   'vcc2025/gene_names.csv') if l.strip()}
wb = openpyxl.load_workbook('horlbeck2016_supp2_TSS.xlsx', read_only=True)
rows = list(wb['Human hg19'].iter_rows(values_only=True))
ix = {h: i for i, h in enumerate(rows[0])}
def clean(v):
    v = str(v); return v[2:-1] if v.startswith('="') and v.endswith('"') else v
regions = []
for r in rows[1:]:
    g = clean(r[ix['gene']])
    if g not in symbols: continue
    coords = [int(r[ix[k]]) for k in ("Primary TSS, 5'", "Primary TSS, 3'",
                                      "Secondary TSS, 5'", "Secondary TSS, 3'")
              if r[ix[k]] not in (None, '')]
    if not coords: continue
    chrom = str(r[ix['chromosome']]).replace('chr', '')
    strand = 1 if str(r[ix['strand']]) == '+' else -1
    tss = min(coords) if strand == 1 else max(coords)
    regions.append((g, str(r[ix['transcript']]), chrom, strand, tss,
                    max(1, min(coords) - FLANK), max(coords) + FLANK))
print(f'symbols measured in 2025: {len(symbols)}; TSS rows kept: {len(regions)}; '
      f'distinct genes: {len({r[0] for r in regions})}', flush=True)

URL = 'https://grch37.rest.ensembl.org/sequence/region/human'
def fetch(batch):
    body = json.dumps({'regions': [f'{c}:{lo}..{hi}:1' for _, _, c, _, _, lo, hi in batch]}).encode()
    for a in range(6):
        try:
            req = urllib.request.Request(URL, data=body, headers={
                'Content-Type': 'application/json', 'Accept': 'application/json'})
            with urllib.request.urlopen(req, timeout=120) as resp:
                return json.loads(resp.read())
        except Exception:
            if a == 5: raise
            time.sleep(2 ** a)
B = 50
batches = [regions[i:i + B] for i in range(0, len(regions), B)]
out = gzip.open('promoters_cage_tss_hg19.fa.gz', 'wt'); done = 0
with ThreadPoolExecutor(max_workers=4) as ex:
    for batch, seqs in zip(batches, ex.map(fetch, batches)):
        assert len(seqs) == len(batch)
        for (g, t, c, strand, tss, lo, hi), s in zip(batch, seqs):
            seq = s['seq'].upper()
            out.write(f'>{g}|{t}|{c}|{strand}|{tss}|{lo}|{hi}|{len(seq)}\n{seq}\n')
        done += len(batch)
        if done % 2000 < B: print(f'  {done}/{len(regions)}', flush=True)
out.close()
print('WROTE promoters_cage_tss_hg19.fa.gz', flush=True)
