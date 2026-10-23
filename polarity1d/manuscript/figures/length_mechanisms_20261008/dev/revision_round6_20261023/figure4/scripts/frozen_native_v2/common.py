"""Round-six paths and small I/O helpers. Source campaigns are read-only."""
from pathlib import Path
import os,sys,csv,json,hashlib,math
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
sys.dont_write_bytecode=True
OUT=Path('/home/4473331/projects/scaling_sandbox/polarity1d/manuscript/figures/length_mechanisms_20261008/dev/revision_round6_20261023/figure4')
ROOT=next(p for p in OUT.parents if (p/'polarity1d/manuscript/code/mechanisms').is_dir())
CODE=ROOT/'polarity1d/manuscript/code/mechanisms'
sys.path.insert(0,str(CODE))
OLD=OUT.parents[1]/'revision_round5_20261023/figure4'
SOURCE=ROOT/'polarity1d/manuscript/data/runs/signature_figure4_round3_compact_20261022'
RUN=ROOT/'polarity1d/manuscript/data/runs/figure4_round6_modes_20261023'
PYTHON='/home/4473331/.conda/envs/scaling_sandbox/bin/python'
DATA=OUT/'supporting/data';META=OUT/'supporting/metadata'
MODELS=('wave_pinning','holmes_model3','spring_ring','goryachev')
NAMES=dict(wave_pinning='Wave-pinning',holmes_model3='Holmes',spring_ring='Spring',goryachev='Goryachev')
def read(path):return json.loads(Path(path).read_text())
def json_value(value):
    if isinstance(value,dict):return {key:json_value(v) for key,v in value.items()}
    if isinstance(value,(list,tuple)):return [json_value(v) for v in value]
    if hasattr(value,'item'):return json_value(value.item())
    if isinstance(value,float) and not math.isfinite(value):return None
    return value
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(json_value(value),indent=2,default=str,allow_nan=False)+'\n');temp.replace(path)
def csvread(path):
    if not Path(path).exists():return []
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))
def csvwrite(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def num(r,k,default=float('nan')):
    try:return float(r[k])
    except (ValueError,TypeError,KeyError):return default
