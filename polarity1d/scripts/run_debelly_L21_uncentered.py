#!/usr/bin/env python3
"""Repeat the retained L=21 assay with independent, uncentered endpoints."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from polarity1d.debelly import DeBelly, DeBellyParameters, simulate_debelly
from polarity1d.protocols import LocalizedPulse


def main():
    source=ROOT/'campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000.json'
    tasks=json.loads(source.read_text())['tasks']
    task=next(t for t in tasks if t['coupling_variant']=='full' and t['simulation']['length']==21)
    settings=dict(task['model_parameters'], center_boundary_motion=False)
    protocol={k:v for k,v in task['protocol'].items() if k!='name'}
    trajectory=simulate_debelly(DeBelly(DeBellyParameters(**settings)),
        protocol=LocalizedPulse(**protocol),**task['simulation'])
    trajectory.metadata.update(source_task=task['task_id'],source_manifest=str(source),
        protocol=protocol,changed_parameter={'center_boundary_motion':False},
        substrate_clutch_enabled=False)
    out=ROOT/'results/campaigns/debelly_L21_uncentered'
    out.mkdir(parents=True,exist_ok=True)
    trajectory.save(out/'trajectory.npz')
    f=trajectory.fields
    summary=dict(trajectory.metadata,
        left_displacement_range_um=[float(f['left_boundary_displacement'].min()),float(f['left_boundary_displacement'].max())],
        right_displacement_range_um=[float(f['right_boundary_displacement'].min()),float(f['right_boundary_displacement'].max())],
        length_range_um=[float(21+f['boundary_separation'].min()),float(21+f['boundary_separation'].max())],
        right_rac_range=[float(f['rac'][:,-1].min()),float(f['rac'][:,-1].max())],
        protrusion_threshold=1.3*float(f['rac'][0,0]))
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
