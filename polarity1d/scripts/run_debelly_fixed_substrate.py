#!/usr/bin/env python3
"""Run illustrative fixed-substrate clutch extensions of the retained L=21 assay."""
import argparse
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from polarity1d.debelly import DeBelly, DeBellyParameters, simulate_debelly
from polarity1d.protocols import LocalizedPulse
from polarity1d.substrate_clutch import FixedSubstrateClutch, FixedSubstrateParameters


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--density', type=float, default=10)
    parser.add_argument('--end', type=float, default=10000)
    parser.add_argument('--step', type=float, default=0.125)
    parser.add_argument('--spacing', type=float, default=0.025)
    parser.add_argument('--cells', type=int, default=105)
    parser.add_argument('--local-rac', action='store_true')
    parser.add_argument('--uncentered', action='store_true')
    parser.add_argument('--passive-rear', action='store_true', help='Option 2 at the right end of the left-cue assay')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'campaigns/debelly_2026_standalone_figure2_julia_activation_defaults_panel_a_10000.json').read_text())
    matches = [t for t in manifest['tasks'] if t['coupling_variant']=='full' and t['simulation']['length']==21]
    assert len(matches)==1
    task = matches[0]
    parameters = DeBellyParameters(**dict(task['model_parameters'],
        center_boundary_motion=not (args.uncentered or args.passive_rear or args.local_rac),
        passive_boundary='local_rac' if args.local_rac else ('right' if args.passive_rear else 'none')))
    protocol_settings = {k:v for k,v in task['protocol'].items() if k!='name'}
    protocol = LocalizedPulse(**protocol_settings)
    clutch = FixedSubstrateClutch(FixedSubstrateParameters(density=args.density, spacing=args.spacing), keep_history=True)
    started = time.monotonic()
    termination = None
    try:
        trajectory = simulate_debelly(DeBelly(parameters), substrate_clutch=clutch, protocol=protocol,
            length=21, cells=args.cells, end_time=args.end, output_interval=2, step=args.step)
    except FloatingPointError as error:
        if not args.local_rac or 'folded the coordinate map at t=' not in str(error):
            raise
        failure_time = float(str(error).split('t=')[-1])
        end = 2*((failure_time-1e-9)//2)
        termination = dict(reason=str(error), requested_end=args.end, last_saved_time=end)
        print(termination, flush=True)
        clutch = FixedSubstrateClutch(FixedSubstrateParameters(density=args.density, spacing=args.spacing), keep_history=True)
        trajectory = simulate_debelly(DeBelly(parameters), substrate_clutch=clutch, protocol=protocol,
            length=21, cells=args.cells, end_time=end, output_interval=2, step=args.step)
    if termination:
        trajectory.metadata['early_termination'] = termination
    args.output.mkdir(parents=True, exist_ok=True)
    trajectory.metadata['source_task'] = task['task_id']
    trajectory.metadata['protocol'] = protocol_settings
    trajectory.metadata['interpretation'] = 'Exploratory coupled clutch extension, arbitrary clutch parameters; not a calibrated literature prediction.'
    trajectory.save(args.output/'trajectory.npz')
    clutch.save_history(args.output/'fixed_substrate_history.npz')
    summary = dict(trajectory.metadata)
    f = trajectory.fields
    midpoint = f['cell_translation'][:,0] + .5*(f['left_boundary_displacement'][:,0]+f['right_boundary_displacement'][:,0])
    summary.update(final_translation_um=float(clutch.translation),
                   final_midpoint_displacement_um=float(midpoint[-1]),
                   midpoint_displacement_range_um=[float(midpoint.min()),float(midpoint.max())],
                   min_translation_um=float(trajectory.fields['cell_translation'].min()),
                   max_translation_um=float(trajectory.fields['cell_translation'].max()),
                   max_speed_um_s=float(abs(trajectory.fields['cell_velocity']).max()),
                   length_min_um=float(21+trajectory.fields['boundary_separation'].min()),
                   length_max_um=float(21+trajectory.fields['boundary_separation'].max()),
                   total_wall_seconds=time.monotonic()-started)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__=='__main__':
    main()
