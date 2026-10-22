"""Single Figure 4 rebuild, defaulting to already frozen compact plot records.

Examples (repository root):
  conda run -n scaling_sandbox python polarity1d/manuscript/code/mechanisms/rebuild_figure4_round3_20261022.py
  .../rebuild_figure4_round3_20261022.py --campaign /path/to/complete/campaign

Expensive --campaign, --prepare-round2 and --nulls analysis belongs in SLURM.
The default portable redraw performs no integration and no data deletion.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import shutil
from datetime import datetime,timezone

import figure4_round3_render_20261022 as R
import figure4_round3_documents_20261022 as D


def rebuild(output=R.OUT,*,campaign=None,prepare_round2=False,nulls='none',prepare_examples=False):
    output=Path(output)
    for folder in ('data','metadata','vectors','viewers'):(output/'supporting'/folder).mkdir(parents=True,exist_ok=True)
    if campaign is not None or prepare_round2 or nulls!='none':
        import figure4_round3_analysis_20261022 as A
        if campaign is not None:
            A.compact(output,Path(campaign))
            shutil.copy2(Path(campaign)/'manifest.json',output/'supporting/metadata/campaign_manifest.json')
            status=json.loads((Path(campaign)/'status.json').read_text())
            progress_path=output/'supporting/metadata/campaign_progress_receipt.json'
            progress=json.loads(progress_path.read_text()) if progress_path.exists() else {}
            progress.update(generated_utc=datetime.now(timezone.utc).isoformat(),status=status.get('status'),complete=status.get('complete',False),precision_qualified=status.get('precision_qualified',False),source_status=status,source_status_path=str(Path(campaign)/'status.json'),source_status_sha256=R.sha(Path(campaign)/'status.json'))
            progress.update(last_integrity_collection_utc=status.get('generated_utc'),verified_compact_seeds=status.get('native_seeds_complete'),final_allocated_seed_count=status.get('native_seeds_expected'),failed_seed_count=len(status.get('failed_seeds',[])) if isinstance(status.get('failed_seeds'),list) else status.get('failed_seeds'),running_seeds_at_last_collection=status.get('running_seeds'),storage_record=dict(generated_utc=status.get('generated_utc'),allocated_bytes=status.get('allocated_bytes'),measurement='Allocated storage recorded by the final campaign collection; exact source status hash retained'),old_native_deletion_performed=status.get('old_native_deletion_performed',False))
            R.write_json(progress_path,progress)
            for name,relative in [('campaign_final_allocation.json','results/allocation.json'),('campaign_pipeline_status.json','results/pipeline_status.json')]:
                path=Path(campaign)/relative
                if path.exists():shutil.copy2(path,output/'supporting/metadata'/name)
        elif prepare_round2:A.fallback(output)
        if nulls!='none':A.nulls(output,nulls)
    if prepare_examples or prepare_round2 or campaign is not None:
        import figure4_round3_examples_20261022 as X
        X.prepare(output,campaign=campaign)
    mandatory=['condition_lag_metrics.csv','local_alpha.csv','noise_tolerance.csv','noise_tolerance_folds.csv','null_A_acf.csv','null_A_kymographs.csv','null_B_local_alpha.csv','null_B_noise_tolerance_folds.csv','model_D_kymographs.csv','drift_decomposition.csv','valid_fraction.csv']
    missing=[name for name in mandatory if not (output/'supporting/data'/name).exists()]
    if missing:raise FileNotFoundError('Compact plot records are incomplete: '+', '.join(missing)+'. Prepare model summaries and nulls in SLURM before a portable redraw.')
    R.render(output)
    legend=D.create(output)
    receipt=output/'supporting/metadata/render_receipt.json'
    record=json.loads(receipt.read_text());record['documentation']={str(p.relative_to(output)):R.sha(p) for p in [output/'legend.json',output/'index.html',output/'supporting/metadata/main_legend.md',output/'supporting/metadata/scientific_readme.md']}
    R.write_json(receipt,record)
    print(json.dumps(dict(status=legend['status'],main=str(output/'figure4_main.pdf'),supplement=str(output/'figure4_supplement.pdf'),metadata=str(receipt)),indent=2))
    return legend


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=R.OUT)
    parser.add_argument('--campaign',type=Path)
    parser.add_argument('--prepare-round2',action='store_true')
    parser.add_argument('--prepare-examples',action='store_true')
    parser.add_argument('--nulls',choices=['none','A','B','all'],default='none')
    args=parser.parse_args()
    rebuild(args.output,campaign=args.campaign,prepare_round2=args.prepare_round2,nulls=args.nulls,prepare_examples=args.prepare_examples)
