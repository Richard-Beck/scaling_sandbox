"""Reproduce corrected preparations, bracket extension and three reversal assays.

Use a fresh --out. Both preparation stages precede exact-checkpoint reversal;
unchanged four-model evidence is independently audited with the new deBelly runs.
"""
from pathlib import Path
import argparse,json
import debelly_corrected_reversal as d
import debelly_preparation_refine as p

def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['launch','prep-refine']);parser.add_argument('--out',type=Path,required=True);parser.add_argument('--slot',type=int);parser.add_argument('--offset',type=int,default=0);parser.add_argument('--round',type=int,default=0)
 args=parser.parse_args();out=args.out.resolve();p.OUT=out
 if args.action=='prep-refine':p.worker(args.slot+args.offset);return
 d.prep_init(out);p.setup();count=len(d.prep_verify(out)['cases'])
 d.r.RUNNER=Path(d.__file__);initial=d.r.submit(out,'prep-worker',0,array=f'0-{count-1}')
 d.r.RUNNER=Path(__file__).resolve();refined=d.r.submit(out,'prep-refine',0,array=f'0-{count-1}',dependency=initial)
 d.r.RUNNER=Path(d.__file__);final=d.r.submit(out,'prep-finish',0,dependency=refined)
 d.r.write(out/'preparation_launch.json',dict(generated_at=d.n.now(),initial_array=initial,refinement_array=refined,controller=final,launcher_sha256=d.n.sha(Path(__file__))))
 print(json.dumps(dict(initial_array=initial,refinement_array=refined,controller=final)))
if __name__=='__main__':main()
