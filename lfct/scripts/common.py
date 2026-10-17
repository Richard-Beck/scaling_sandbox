"""Shared configuration, checksums and atomic table/JSON output."""
from pathlib import Path
import csv,gzip,hashlib,json,os
BASE=Path(__file__).resolve().parents[1]
REPO=BASE.parent

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
    return h.hexdigest()

def digest(obj):return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':')).encode()).hexdigest()
def config():return json.loads(Path(os.environ.get('LFCT_CONFIG',str(BASE/'config.json'))).read_text())
def read_json(path):return json.loads(Path(path).read_text())
def write_json(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);tmp=path.with_suffix(path.suffix+'.part');tmp.write_text(json.dumps(obj,indent=2,allow_nan=False));tmp.replace(path)
def catalog():
    with (BASE/'data/sources.csv').open() as f:return list(csv.DictReader(f))
def field_row(key):return next(r for r in catalog() if r['key']==key)
def profile_fields(profile):
    fields=config()['profiles'][profile]['fields']
    return [r['key'] for r in catalog()] if fields=='all' else fields

def frame_selection(row,profile):
    n=int(row['nframes'])
    if profile=='full':return list(range(n))
    if profile=='smoke':return list(range(1,13)) if row['key']=='mcf10a/FLD_102' else [1,n//2,n-2]
    if profile=='pilot':
        return sorted({f for start in [1,n//2-5,n-11] for f in range(start,start+10)})
    raise ValueError(profile)

def csv_write(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=fields or list(rows[0]) if rows else fields or []
    tmp=path.with_suffix(path.suffix+'.part')
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(tmp,'wt',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    tmp.replace(path)
def csv_read(path):
    opener=gzip.open if Path(path).suffix=='.gz' else open
    with opener(path,'rt',newline='') as f:return list(csv.DictReader(f))
def output_dir(profile,key):return BASE/'results'/profile/'fields'/key

def code_hash():
    return digest({p.name:sha(p) for p in sorted((BASE/'scripts').glob('*.py')) if p.name in ['common.py','prepare.py','segment.py','geometry.py','analyze.py','pipeline.py']})
