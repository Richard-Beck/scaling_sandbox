"""Redraw final F4 main and comprehensive supplement without native archives.

Run the exported helpers/redraw.py --out /path/to/new/output. The lineage paths
are retained as metadata but never opened. Requires numpy,pandas,matplotlib,
pypdf in an existing Python environment; no solver or native adapter is loaded.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse
import gzip
import hashlib
import importlib.util
import json
import shutil
import sys
sys.dont_write_bytecode=True
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True)
    p.add_argument('--recompute-ci',action='store_true',help='Also replay joint whole-seed CI from compact case sufficient statistics')
    a=p.parse_args()
    bundle=Path(__file__).resolve().parent.parent;out=a.out.resolve();assert out!=bundle
    def block_native(event,args):
        if event=='open' and args and isinstance(args[0],(str,bytes)):
            filename=args[0].decode() if isinstance(args[0],bytes) else args[0]
            if '/polarity1d/manuscript/data/runs/' in str(Path(filename).resolve()):
                raise RuntimeError('Portable redraw attempted to open native lineage: '+filename)
    sys.addaudithook(block_native)
    receipt=json.loads((bundle/'metadata/EXPORT_RECEIPT.json').read_text())
    for path,digest in receipt['outputs'].items():assert sha(bundle/path)==digest,path
    for name in ('supporting/data','supporting/metadata','supporting/vectors'):(out/name).mkdir(parents=True,exist_ok=True)
    for source in (bundle/'tables').glob('*.csv'):shutil.copy2(source,out/'supporting/data'/source.name)
    for source in (bundle/'metadata').glob('*.json'):shutil.copy2(source,out/'supporting/metadata'/source.name)
    index=json.loads((bundle/'metadata/FIELD_INDEX.json').read_text())
    def module(name):
        spec=importlib.util.spec_from_file_location('portable_f4_'+name,bundle/'helpers'/(name+'.py'))
        result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result
    def field(case):
        entry=index[case['portable_case_id']]
        with gzip.open(bundle/entry['path'],'rt') as f:data=np.loadtxt(f,delimiter=',',skiprows=1)
        return dict(time=data[:,0],C=data[:,1],R=data[:,2],A=data[:,3],Q=data[:,4],regions=data[:,5].astype(int),contrast_ratio=data[:,6],
            primary=data[:,7:],coordinate=np.array(entry['coordinate']),scalar_time=data[:,0])
    if a.recompute_ci:
        module('assembly').analyze(out)
        for name in ('condition_window_bootstrap.csv','memory_window_bootstrap.csv','paired_halfwindow_bootstrap.csv'):
            assert sha(out/'supporting/data'/name)==sha(bundle/'tables'/name),name
    local=module('local');local.RUN=bundle/'selection_cache';local.field=field
    sys.argv=['portable_local','--out',str(out)];local.main()
    for name,arguments in [('temporal',['--source',str(out/'supporting/data/memory_all_conditions_windows_lags.csv')]),
                           ('organization',['--conditions',str(out/'supporting/data/condition_window_summaries.csv')]),('quarters',[])]:
        imported=module(name);sys.argv=['portable_'+name,'--out',str(out),*arguments];imported.main()
    compact=module('compact_main');sys.argv=['portable_main','--source',str(out),'--out',str(out)];compact.main()
    from pypdf import PdfReader,PdfWriter
    writer=PdfWriter()
    for name in ('all_models_sizes_strengths_temporal_windows.pdf','all_models_sizes_strengths_organization_windows.pdf','all_conditions_second_half_native_atlas.pdf','late_quarter_sensitivity.pdf'):
        writer.append(PdfReader(out/'supporting/vectors'/name),outline_item=name)
    with (out/'figure4_comprehensive_supplement.pdf').open('wb') as handle:writer.write(handle)
    pages=len(writer.pages);writer.close()
    expected=json.loads((bundle/'metadata/EXPECTED_PNG_HASHES.json').read_text())
    omitted=set(expected)-{name for name in expected if (out/name).exists()}
    assert omitted <= {'figure4_signature_readouts.png'},omitted
    checked={name:dict(expected=digest,actual=sha(out/name),exact_png_match=sha(out/name)==digest) for name,digest in expected.items() if (out/name).exists()}
    assert checked and all(v['exact_png_match'] for v in checked.values())
    result=dict(generated_at=datetime.now(timezone.utc).isoformat(),native_archives_opened=0,native_integrations=0,
        export_receipt_sha256=sha(bundle/'metadata/EXPORT_RECEIPT.json'),supplement_pages=pages,
        exact_pngs=len(checked),png_checks=checked,earlier_matrix_preserved_but_not_redrawn=sorted(omitted),
        joint_CI_recomputed=a.recompute_ci,native_path_opens_blocked_by_audit_hook=True,
        source=str(Path(__file__).resolve()),status='passed')
    (out/'supporting/metadata/PORTABLE_REDRAW_CHECK.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
