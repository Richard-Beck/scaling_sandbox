#!/usr/bin/env python3
"""Audit production trajectories and report actual refinement differences."""
import json
import argparse
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CAMPAIGN=ROOT/'results/campaigns/debelly_fixed_substrate_L21'


def compare(base_name, other_name):
    result={}
    with np.load(CAMPAIGN/base_name/'trajectory.npz') as base, np.load(CAMPAIGN/other_name/'trajectory.npz') as other:
        np.testing.assert_allclose(base['time'],other['time'])
        for field in ('cell_translation','midpoint_displacement','boundary_separation','rac','velocity','clutch_bound'):
            if field=='midpoint_displacement':
                a=base['field__cell_translation']+.5*(base['field__left_boundary_displacement']+base['field__right_boundary_displacement'])
                b=other['field__cell_translation']+.5*(other['field__left_boundary_displacement']+other['field__right_boundary_displacement'])
            else:
                a=base['field__'+field]
                b=other['field__'+field]
            if b.shape[1]!=a.shape[1]:
                b=np.array([np.interp(base['x'],other['x'],row) for row in b])
            scale=max(float(np.max(np.abs(a))),1e-10)
            result[field]={
                'max_abs_error':float(np.max(np.abs(a-b))),
                'max_error_relative_to_base_max':float(np.max(np.abs(a-b)))/scale,
                'rms_error_relative_to_base_max':float(np.sqrt(np.mean((a-b)**2)))/scale,
                'final_mean_difference':float(np.mean(b[-1]-a[-1])),
            }
    return result


def main():
    global CAMPAIGN
    parser=argparse.ArgumentParser()
    parser.add_argument('--uncentered',action='store_true')
    parser.add_argument('--passive-rear',action='store_true')
    args=parser.parse_args()
    if args.uncentered:
        CAMPAIGN=ROOT/'results/campaigns/debelly_fixed_substrate_L21_uncentered'
    if args.passive_rear:
        args.uncentered=True
        CAMPAIGN=ROOT/'results/campaigns/debelly_fixed_substrate_L21_passive_rear'
    report={}
    for density in (1,10,50):
        name=f'density_{density}'
        summary=json.loads((CAMPAIGN/name/'summary.json').read_text())
        assert summary['max_net_substrate_force_pN']<1e-10
        assert summary['mca_conservation_error']<1e-12
        with np.load(CAMPAIGN/name/'trajectory.npz') as z:
            assert 0<=np.min(z['field__clutch_bound'])<=np.max(z['field__clutch_bound'])<=1
            np.testing.assert_allclose(z['field__cortical_velocity_lab'],z['field__velocity']+z['field__cell_velocity'])
            net=np.sum(z['field__clutch_traction'],axis=1)*summary['dx_um']
            assert np.max(np.abs(net))<1e-10
        report[name]={key:summary[key] for key in (
            'final_translation_um','length_min_um','length_max_um',
            'max_net_substrate_force_pN','mca_conservation_error','fixed_step')}
        if args.uncentered:
            assert summary['parameters']['center_boundary_motion'] is False
            report[name]['final_midpoint_displacement_um']=summary['final_midpoint_displacement_um']
        if args.passive_rear:
            assert summary['parameters']['passive_boundary']=='right'
            residuals=summary['maximum_boundary_residuals']
            assert residuals['passive_stress_residual']<1e-6
            assert residuals['global_force_residual']<1e-6
            assert residuals['passive_kinematic_residual']<1e-9
            report[name]['maximum_boundary_residuals']=residuals
    for name in ('check_timestep','check_substrate','check_cell_grid'):
        report[name]=compare('check_base',name)
    fine_summary=json.loads((CAMPAIGN/'check_long_timestep/summary.json').read_text())
    assert fine_summary['fixed_step']==report['density_10']['fixed_step']/2
    report['check_long_timestep']=compare('density_10','check_long_timestep')
    (CAMPAIGN/'numerical_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    # These are exploratory accuracy gates, not claims of exact convergence.
    assert report['check_long_timestep']['cell_translation']['max_error_relative_to_base_max']<.03
    assert report['check_long_timestep']['midpoint_displacement']['max_error_relative_to_base_max']<.03
    assert report['check_substrate']['cell_translation']['max_error_relative_to_base_max']<.005
    assert report['check_cell_grid']['cell_translation']['max_error_relative_to_base_max']<.005
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
