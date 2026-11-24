"""Read a remote zip's central directory over HTTP Range requests (no full download).
Usage: uv run python remotezip270.py <url> <listing_out.tsv> [member_regex] [outdir] [max_bytes]
"""
import sys, zipfile, urllib.request, io, re, os, time
class RangeFile(io.RawIOBase):
    def __init__(self, url, chunk=1<<20):
        self.url=url; self.pos=0; self.chunk=chunk; self.cache={}; self.fetched=0
        self._resolve()
    def _resolve(self):
        # figshare redirects to a presigned S3 URL that is method-specific and expires: resolve with a ranged GET
        req=urllib.request.Request(self.url, headers={'Range': 'bytes=0-0'}); r=urllib.request.urlopen(req, timeout=120)
        cr=r.headers.get('Content-Range'); self.size=int(cr.split('/')[-1]) if cr else int(r.headers['Content-Length']); self.final=r.geturl(); r.read()
    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        if whence==0: self.pos=off
        elif whence==1: self.pos+=off
        else: self.pos=self.size+off
        return self.pos
    def _get(self, a, b):
        req=urllib.request.Request(self.final, headers={'Range': f'bytes={a}-{b-1}'})
        for k in range(5):
            try:
                r=urllib.request.urlopen(req, timeout=300); d=r.read(); self.fetched+=len(d); return d
            except Exception as e:
                err=e; time.sleep(2*(k+1))
                if getattr(e, 'code', None) in (403, 400):
                    self._resolve(); req=urllib.request.Request(self.final, headers={'Range': f'bytes={a}-{b-1}'})
        raise err
    def readinto(self, b):
        d=self.read(len(b)); b[:len(d)]=d; return len(d)
    def read(self, n=-1):
        if n<0: n=self.size-self.pos
        a=self.pos; b=min(self.size, a+n); out=bytearray()
        while a<b:
            ck=a//self.chunk
            if ck not in self.cache:
                s=ck*self.chunk; self.cache[ck]=self._get(s, min(self.size, s+self.chunk))
            blk=self.cache[ck]; s=ck*self.chunk; off=a-s; take=min(len(blk)-off, b-a)
            out+=blk[off:off+take]; a+=take
        self.pos=b; return bytes(out)
url, out = sys.argv[1], sys.argv[2]
rf=RangeFile(url)
print('size', rf.size, 'final', rf.final[:120])
zf=zipfile.ZipFile(rf)
infos=zf.infolist()
with open(out,'w') as f:
    f.write('name\tsize\tcompress_size\tdate\n')
    for i in infos: f.write(f'{i.filename}\t{i.file_size}\t{i.compress_size}\t{i.date_time}\n')
print('members', len(infos), 'fetched_bytes_for_listing', rf.fetched)
if len(sys.argv)>3:
    rx=re.compile(sys.argv[3]); outdir=sys.argv[4]; cap=int(sys.argv[5]) if len(sys.argv)>5 else 1_500_000_000
    os.makedirs(outdir, exist_ok=True); n=0; tot=0
    for i in infos:
        if rx.search(i.filename) and not i.is_dir():
            if tot+i.compress_size>cap: print('CAP reached at', n, tot); break
            data=zf.read(i.filename); tot+=i.compress_size; n+=1
            p=os.path.join(outdir, i.filename.replace('/','__'))
            open(p,'wb').write(data)
    print('extracted', n, 'members; compressed bytes', tot, 'total fetched', rf.fetched)
