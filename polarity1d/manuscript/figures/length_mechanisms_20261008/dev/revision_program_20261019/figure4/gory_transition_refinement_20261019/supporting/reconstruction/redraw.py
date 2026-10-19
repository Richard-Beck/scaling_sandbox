"""Portable Gory refinement redraw from compact field/CI exports only."""
from pathlib import Path
from datetime import datetime,timezone
import argparse,ast,gzip,hashlib,json,sys
sys.dont_write_bytecode=True
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    source=a.source.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    closure=json.loads((source/'supporting/reconstruction/INPUT_CLOSURE.json').read_text())
    for name,digest in closure['inputs'].items():assert sha(source/name)==digest,name
    def block(event,args):
        if event=='open' and args and isinstance(args[0],(str,bytes)):
            name=args[0].decode() if isinstance(args[0],bytes) else args[0]
            if '/polarity1d/manuscript/data/runs/' in str(Path(name).resolve()):raise RuntimeError('Native path opened by portable redraw: '+name)
    sys.addaudithook(block)
    data=source/'supporting/data';ci=pd.read_csv(data/'condition_window_bootstrap.csv',float_precision='round_trip');cm=pd.read_csv(data/'memory_window_bootstrap.csv',float_precision='round_trip');native=[]
    for f in (.01,.015,.02,.03):
        with gzip.open(data/f'native_medoid_f{f:g}.csv.gz','rt') as h:matrix=np.loadtxt(h,delimiter=',',skiprows=1)
        profile=pd.read_csv(data/f'profiles_f{f:g}.csv',float_precision='round_trip')
        native.append(dict(f=f,t=matrix[1:,0],u=matrix[1:,1:],x=profile.x.to_numpy(),basal=matrix[0,1:],snapshot=profile.noisy.to_numpy(),mean=profile.matched_mean.to_numpy()))
    path=source/'supporting/reconstruction/plot_source.py';tree=ast.parse(path.read_text())
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='main')
    start=next(i for i,n in enumerate(function.body) if isinstance(n,ast.Expr) and ast.unparse(n).startswith('plt.rcParams.update'))
    stop=next(i for i,n in enumerate(function.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='receipt' for t in n.targets))
    nodes=function.body[start:stop];assert len(nodes)==2 and isinstance(nodes[1],ast.With)
    plotting=ast.Module(body=nodes,type_ignores=[]);namespace=dict(plt=plt,PdfPages=PdfPages,np=np,pd=pd,OUT=out,ci=ci,cm=cm,native=native)
    exec(compile(ast.fix_missing_locations(plotting),str(path),'exec'),namespace)
    result=dict(generated_at=datetime.now(timezone.utc).isoformat(),native_archive_opens=0,native_integrations=0,native_path_opens_blocked=True,
        source_closure_sha256=sha(source/'supporting/reconstruction/INPUT_CLOSURE.json'),plot_source_sha256=sha(path),
        ast_scope='Only unchanged figure-construction block; legend placement amended separately from frozen scientific collection',
        outputs={name:sha(out/name) for name in ('gory_transition_review.pdf','gory_transition_readouts.png','gory_transition_native.png')})
    (out/'PORTABLE_REDRAW_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
