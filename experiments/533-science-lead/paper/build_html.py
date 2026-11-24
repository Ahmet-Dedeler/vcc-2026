"""Assemble the manuscript artifact: pandoc renders PAPER.md and CLAIMS.md (GFM) to HTML fragments; this script wraps them with
the page chrome, a table of contents, the figure gallery (PNGs embedded as data URIs, captions from figs/MANIFEST.md)."""
import base64, html, re, subprocess, sys
from pathlib import Path
ROOT = Path('~/Code/virtual-cell-challenge-2026'); W = ROOT/'experiments/533-science-lead/paper'  # 533's copy (Draft v15); 294's original points at experiments/294-writeup
OUT = Path(sys.argv[1])
def pandoc(md: Path, shift=0):
    args = ['pandoc', str(md), '-f', 'gfm+tex_math_dollars', '-t', 'html5', '--wrap=none']
    if shift: args += [f'--shift-heading-level-by={shift}']
    return subprocess.run(args, capture_output=True, text=True, check=True).stdout
paper = pandoc(W/'PAPER.md'); claims = pandoc(W/'CLAIMS.md', shift=1)
# strip the H1 title from the paper body (rendered in the masthead instead)
paper = re.sub(r'<h1[^>]*>.*?</h1>\s*', '', paper, count=1, flags=re.S)
# table of contents from h2/h3 ids
toc = []
for m in re.finditer(r'<h([23]) id="([^"]+)">(.*?)</h\1>', paper, flags=re.S):
    lvl, hid, txt = m.group(1), m.group(2), re.sub('<[^>]+>', '', m.group(3))
    toc.append(f'<li class="l{lvl}"><a href="#{hid}">{html.escape(txt)}</a></li>')
# figures: captions from MANIFEST.md, split by "F-n" headings
man = (W/'figs/MANIFEST.md').read_text()
caps = {}
parts = re.split(r'\n(?=#{1,4}\s*\**F-\d+)', '\n'+man)
for part in parts:
    m = re.search(r'F-(\d+)', part)
    if m: caps[f'F-{m.group(1)}'] = pandoc_text = subprocess.run(['pandoc','-f','gfm','-t','html5','--wrap=none'], input=re.sub(r'^#{1,4}.*\n', '', part.strip()+'\n', count=1), capture_output=True, text=True).stdout
figs = []
for i in range(1, 14):
    key = f'F-{i}'; png = W/f'figs/{key}.png'
    b64 = base64.b64encode(png.read_bytes()).decode()
    figs.append(f'<figure id="fig-{key}"><img src="data:image/png;base64,{b64}" alt="Figure {key}"><figcaption><span class="figlabel">Figure {key}</span>{caps.get(key, "")}</figcaption></figure>')
# verdict chips in the claims table
def chip(m):
    word = m.group(1); cls = {'NEW':'new','EXTENDED':'ext','KNOWN':'known','CONTRADICTED':'contra','WITHDRAWN':'withdrawn','RECONCILED':'ext'}[word]
    return f'<span class="chip {cls}">{word}</span>'
claims = re.sub(r'<strong>(NEW|EXTENDED|KNOWN|CONTRADICTED|WITHDRAWN|RECONCILED)</strong>', chip, claims)
paper = re.sub(r'<strong>(NEW|EXTENDED|KNOWN|CONTRADICTED|WITHDRAWN|RECONCILED)</strong>', chip, paper)
# wrap tables for horizontal scroll
paper = paper.replace('<table>', '<div class="tbl"><table>').replace('</table>', '</table></div>')
claims = claims.replace('<table>', '<div class="tbl"><table>').replace('</table>', '</table></div>')
CSS = """
:root{--bg:#F4F6F3;--bg2:#EAEEEA;--ink:#1C2320;--muted:#5A6660;--rule:#C7D0CA;--accent:#0E6B60;--accent-ink:#0A4F47;--rust:#A8552F;--chipbg:#DDEDE9;--code:#EEF1EE}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#151A18;--bg2:#1D2421;--ink:#E4E9E5;--muted:#98A49E;--rule:#2D3833;--accent:#5EC2B3;--accent-ink:#8ADACF;--rust:#E39168;--chipbg:#1F3B36;--code:#202826}}
:root[data-theme="dark"]{--bg:#151A18;--bg2:#1D2421;--ink:#E4E9E5;--muted:#98A49E;--rule:#2D3833;--accent:#5EC2B3;--accent-ink:#8ADACF;--rust:#E39168;--chipbg:#1F3B36;--code:#202826}
body{background:var(--bg);color:var(--ink);font-family:"Newsreader","Source Serif 4",Georgia,serif;font-size:17px;line-height:1.55;margin:0}
a{color:var(--accent-ink);text-decoration-thickness:1px;text-underline-offset:2px}
.wrap{max-width:1180px;margin:0 auto;padding:32px 24px 80px;display:grid;grid-template-columns:230px minmax(0,1fr);gap:48px}
@media (max-width:900px){.wrap{grid-template-columns:minmax(0,1fr);gap:24px} nav.toc{position:static;max-height:none}}
header.mast{grid-column:1/-1;border-bottom:1px solid var(--rule);padding-bottom:20px;margin-bottom:8px}
.eyebrow{font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:12px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted)}
header.mast h1{font-family:"Newsreader",Georgia,serif;font-weight:500;font-size:clamp(28px,3.4vw,40px);line-height:1.15;margin:8px 0 12px;text-wrap:balance;max-width:34ch}
.meta{font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:13.5px;color:var(--muted);display:flex;flex-wrap:wrap;gap:8px 22px}
.tally{display:flex;flex-wrap:wrap;gap:8px;margin-top:14px}
nav.toc{position:sticky;top:20px;align-self:start;max-height:calc(100vh - 40px);overflow:auto;font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:13px;line-height:1.35}
nav.toc ul{list-style:none;padding:0;margin:0;display:flex;flex-direction:column;gap:6px}
nav.toc li.l3{padding-left:14px;color:var(--muted)}
nav.toc a{color:inherit;text-decoration:none}nav.toc a:hover,nav.toc a:focus{color:var(--accent-ink);text-decoration:underline}
main{min-width:0;max-width:74ch}
main h2{font-family:"Newsreader",Georgia,serif;font-weight:500;font-size:26px;margin:44px 0 12px;text-wrap:balance;border-top:1px solid var(--rule);padding-top:18px}
main h3{font-family:"Newsreader",Georgia,serif;font-weight:600;font-size:19px;margin:30px 0 8px;text-wrap:balance}
main h4{font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:14px;letter-spacing:.04em;text-transform:uppercase;color:var(--muted);margin:24px 0 6px}
main p{margin:0 0 14px}main li{margin:4px 0}
blockquote{margin:14px 0;padding:8px 18px;border-left:3px solid var(--accent);background:var(--bg2);font-style:italic}
code{font-family:"IBM Plex Mono",ui-monospace,Menlo,monospace;font-size:.86em;background:var(--code);padding:1px 4px;border-radius:3px}
.tbl{overflow-x:auto;margin:12px 0 18px;border:1px solid var(--rule);border-radius:4px}
table{border-collapse:collapse;font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:13px;line-height:1.35;font-variant-numeric:tabular-nums;min-width:100%}
th,td{padding:7px 10px;border-bottom:1px solid var(--rule);vertical-align:top;text-align:left}th{background:var(--bg2);font-weight:600}
tr:last-child td{border-bottom:none}
.chip{display:inline-block;font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:11px;font-weight:600;letter-spacing:.06em;padding:1px 7px;border-radius:3px;background:var(--chipbg);color:var(--accent-ink);white-space:nowrap}
.chip.known{background:var(--bg2);color:var(--muted);border:1px solid var(--rule)}
.chip.contra,.chip.withdrawn{background:transparent;color:var(--rust);border:1px solid var(--rust)}
figure{margin:28px 0;padding:0}figure img{width:100%;height:auto;display:block;border:1px solid var(--rule);background:#fff}
figcaption{font-family:"IBM Plex Sans",system-ui,sans-serif;font-size:13px;color:var(--muted);margin-top:8px;line-height:1.4}
.figlabel{font-weight:600;color:var(--ink);margin-right:8px}
figcaption p{margin:0 0 6px}
section.claims h2{margin-top:56px}
hr{border:0;border-top:1px solid var(--rule);margin:28px 0}
:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
@media (prefers-reduced-motion: no-preference){html{scroll-behavior:smooth}}
"""
title = "Three Reproducibilities"
page = f"""<title>{title}</title>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:ital,opsz,wght@0,6..72,400;0,6..72,500;0,6..72,600;1,6..72,400&family=IBM+Plex+Sans:wght@400;600&family=IBM+Plex+Mono&display=swap">
<style>{CSS}</style>
<div class="wrap">
<header class="mast">
<div class="eyebrow">Local manuscript draft · Virtual Cell Challenge 2026 research programme · experiments 294 and 533 · Draft v15 · 2026-09-30</div>
<h1>Three reproducibilities: what a single-cell perturbation reference measures, which responses are conserved across cell contexts, and what bounds zero-shot transfer</h1>
<div class="meta"><span>Sources: FINDINGS v3.1 · confirmation ledger · audits 287, 297, 304, 320</span><span>Round 7 correction: 320 grades F16 NEW 0 / EXTENDED 1 / KNOWN 4 — a re-measurement, not a discovery. The partner-construct control is credited to Jackson 2003; the headline p is a gene-sampling p; the wet-lab recommendation is withdrawn.</span><span>Round 8: §3(i) is new — an L1 fold-change error scored only on significant genes has a closed-form optimum under sign uncertainty, confirmed on the benchmark's own reference (C-58). Draft v11 (2026-09-24)</span><span>Round 9: §3(j) is new — the DE-sign gap to an oracle is the recipient line's own shared response programme; a screen of the same line from another lab carries it and most of the target-specific response (+0.093 through the benchmark's scorer at 35 % coverage, H1); provenance pending audit (closest prior work: Systema 2025). Draft v12 (2026-09-26)</span><span>Round 10: the first scored test of the same-line programme on a challenge context (Jurkat, an independent Jurkat screen) is null (sign fidelity +0.0002 vs +0.0015 … +0.004 preregistered); dilution does not explain it; the rows remain H1-only evidence, and a bank's norms do not carry effect size. §3(k) new: the scorer's members read two functionals of the same counts; ballast cells decouple them on the board (F19). Draft v13 (2026-09-26)</span><span>Round 11: a programme estimated from one screen partly measures its control sample; judged with control halves, one lab's cancer lines share about a third of it. Draft v14 (2026-09-26/27, with the 503 audit's corrections)</span><span>Round 14: six scored entries and experiments 513–527 and 540–546. A zero-shot pipeline re-fitted from the panel alone reproduces its tuned self on the board; a second pure fold-change delivery transfers at 1.19×; target-specific public sign sources are null on every cancer fold; a same-lab programme's uncalled-gene hook has a lower bar than its block prior; §3(l) new: reach is sign-limited and per-target call-set sizing is closed by payoff (F20); §4.7 open problems. The 09-30 same-lab test of F18 is PENDING. Draft v15 (2026-09-30)</span></div>
</header>
<nav class="toc" aria-label="Contents"><div class="eyebrow" style="margin-bottom:8px">Contents</div><ul>{''.join(toc)}<li class="l2"><a href="#figures">Figures</a></li><li class="l2"><a href="#claims">Claims ledger</a></li></ul></nav>
<main>
{paper}
<h2 id="figures">Figures</h2>
{''.join(figs)}
<section class="claims" id="claims">
{claims}
</section>
</main>
</div>
"""
OUT.write_text(page); print('wrote', OUT, OUT.stat().st_size//1024, 'KB; toc entries', len(toc), '; captions found', sorted(caps))
