"""Derived-only narrow-size arc collector; never integrates or submits jobs."""
from pathlib import Path
from datetime import datetime, timezone
from types import SimpleNamespace
import argparse
import fcntl
import hashlib
import json
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import pandas as pd
import revision_standard_angles_duration_20261019 as base
import revision_standard_angles_report_20261019 as prior
import revision_stationary_characterization_20261019 as rich
import revision_standard_angles_sizes_20261019 as spring_api

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'data/runs/revision_arc_narrow_grid_20261019'
OUT = ROOT / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure3/narrow_arc_experiment'
SOURCE = {}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def record(p, expected=None):
    p = Path(p)
    digest = sha(p)
    if expected is not None:
        assert digest == expected, str(p)
    SOURCE[str(p)] = digest
    return p


def read(p):
    return json.loads(Path(p).read_text())


def dest(name, group='data'):
    p = OUT / 'supporting' / group / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


def jsonwrite(name, value):
    p = dest(name, 'metadata')
    q = p.with_suffix('.tmp.json')
    q.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    q.replace(p)


def csvwrite(name, rows):
    p = dest(name)
    q = p.with_suffix('.tmp.csv')
    pd.DataFrame(rows).to_csv(q, index=False)
    q.replace(p)


def level_edges(q, u, origin, level):
    order = np.argsort((q - origin + .5) % 1 - .5)
    x = ((q - origin + .5) % 1 - .5)[order]
    values = u[order]
    threshold = float(values.min() + level * np.ptp(values))
    inside = np.flatnonzero(values >= threshold)
    assert inside.size and np.all(np.diff(inside) == 1)
    a, b = int(inside[0]), int(inside[-1])
    assert a > 0 and b < len(x) - 1
    left = x[a-1] + (threshold - values[a-1]) / (values[a] - values[a-1]) * (x[a] - x[a-1])
    right = x[b] + (threshold - values[b]) / (values[b+1] - values[b]) * (x[b+1] - x[b])
    return float(left), float(right), threshold


def geometry(q, u, case):
    C, d, width = case['circumference'], case['physical_distance_at_onset'], case['full_physical_support']
    result = {}
    for level, suffix in ((.5, '50'), (.1, '10')):
        left, right, threshold = level_edges(q, u, case['preparation_axis_fraction'], level)
        left *= C
        right *= C
        intersection = sum(max(0., min(right, d + width / 2 + k*C) - max(left, d - width / 2 + k*C))
                           for k in (-1, 0, 1))
        overlap = float(np.clip(intersection / width, 0., 1.))
        centre = any(left <= d + k*C <= right for k in (-1, 0, 1))
        values = dict(initial_edge_left=left, initial_edge_right=right,
            initial_core_width=right-left, initial_native_activity_threshold=threshold,
            signed_centre_distance_from_positive_edge=d-right,
            signed_cue_near_edge_gap_from_positive_core_edge=d-width/2-right,
            centre_relative_to_core='inside' if centre else 'outside',
            footprint_relative_to_core='inside' if overlap >= 1-1e-10 else 'outside' if overlap <= 1e-10 else 'straddling',
            fraction_cue_support_inside_core=overlap, cue_support_core_intersection=intersection)
        for key, value in values.items():
            result[f'{key}_{suffix}'] = value
    return result


def phase(path, spring, Cb):
    with np.load(path) as z:
        time = z['time'].copy()
        full = z['states'].copy() if 'states' in z else None
        if spring:
            assert full is not None
            n = (full.shape[1]-1)//2
            u, C, q = full[:, :n], full[:, -1], z['material_positions'].copy()
        else:
            u = full[:, 0] if full is not None else z['primary'].copy()
            C, q = np.full(len(time), Cb), z['x'].copy()/Cb
        snapshot_times = z['snapshot_times'].copy() if 'snapshot_times' in z else None
        snapshots = z['snapshot_states'].copy() if 'snapshot_states' in z else None
        states = full if full is not None else snapshots
        assert states is not None
        return dict(time=time, u=u, C=C, q=q, full=full, snapshots=snapshots,
                    snapshot_times=snapshot_times, first_state=states[0], final_state=states[-1])


def selected_rhs(case, part, initial_M1, signal):
    spring = case['model'] == 'spring'
    n, Cb = case['cells'], case['circumference']
    weights = np.exp(2j*np.pi*part['q'])
    target = abs(float(rich.wrap_error(case['target_angle_in_P_coordinates_deg'], 0.)))
    antipodal = abs(target - 180.) < 1e-6
    rows = []
    model = spring_api.Closure(case['resting_circumference'], n) if spring else base.ring.native_model(case['model'])
    grid = None if spring else base.ring.RingGrid(Cb, n)
    for wanted in sorted({0., min(10., case['cue_duration']), float(case['cue_duration'])}):
        source_time = part['time'] if part['full'] is not None else part['snapshot_times']
        source_states = part['full'] if part['full'] is not None else part['snapshots']
        found = np.flatnonzero(np.isclose(source_time, wanted, rtol=0, atol=1e-12))
        row = dict(condition_id=case['condition_id'], requested_time=wanted,
                   state_coverage='full_saved_states' if part['full'] is not None else 'full_snapshots_only',
                   RHS_available=bool(len(found)), unavailable_reason='' if len(found) else 'requested_full_state_not_retained')
        if not len(found):
            rows.append(row)
            continue
        j = int(found[0])
        state = source_states[j]
        if spring:
            u, C = state[:n], float(state[-1])
            current_signal = spring_api.physical_raised_cosine(model, state, case['cue_centre_fraction'], case['requested_peak'], .2)[0]
            derivative = model.rhs(state, current_signal)
            dM1 = C*np.mean(derivative[:n]*weights) + derivative[-1]*np.mean(u*weights)
            dM0 = C*np.mean(derivative[:n]) + derivative[-1]*np.mean(u)
        else:
            u, C = state[0], Cb
            derivative = base.ring.native_rhs(model, grid, state, signal)
            dM1, dM0 = C*np.mean(derivative[0]*weights), C*np.mean(derivative[0])
        M1 = C*np.mean(u*weights)
        P = M1/initial_M1
        strong = abs(P) >= .1
        available = bool(strong and not antipodal)
        omega = float(np.degrees(np.imag(np.conj(M1)*dM1)/abs(M1)**2)) if available else None
        row.update(actual_time=float(source_time[j]), R=float(abs(P)), M0_rhs=float(dM0),
                   M1_rhs_real=float(dM1.real), M1_rhs_imag=float(dM1.imag),
                   native_RHS_angular_rate_deg_per_native_time=omega,
                   native_RHS_rate_per_actual_target=omega/target if available and target > 1e-12 else None,
                   angular_rate_available=available,
                   angular_rate_unavailable_reason='antipodal_axis_rate_withheld' if antipodal else 'weak_R' if not strong else '')
        rows.append(row)
    return rows


def derive(case, cue, off, initial, q, signal):
    spring, n, Cb, T = case['model'] == 'spring', case['cells'], case['circumference'], case['cue_duration']
    assert cue['time'][0] == off['time'][0] == 0.
    assert cue['time'][-1] == T and off['time'][-1] == 13000.
    np.testing.assert_array_equal(cue['q'], q)
    np.testing.assert_array_equal(off['q'], q)
    start_error = float(np.max(abs(cue['first_state']-initial)))
    restart_error = float(np.max(abs(cue['final_state']-off['first_state'])))
    scale = max(1., float(np.max(abs(initial))))
    assert start_error <= 32*np.finfo(float).eps*scale
    assert restart_error <= 32*np.finfo(float).eps*max(scale, float(np.max(abs(cue['final_state']))))
    time = np.r_[cue['time'], T+off['time'][1:]]
    u, C = np.vstack((cue['u'], off['u'][1:])), np.r_[cue['C'], off['C'][1:]]
    joff = len(cue['time'])-1
    u0 = initial[:n] if spring else initial[0]
    weights = np.exp(2j*np.pi*q)
    M0, M1, M2 = C*u.mean(axis=1), C*(u@weights)/n, C*(u@(weights**2))/n
    P, R = M1/M1[0], abs(M1/M1[0])
    heading = prior.heading(P)
    masked = np.where(R >= .1, heading, np.nan)
    q1, q2 = abs(M1)/M0, abs(M2)/M0
    core = base.core_info(q, u0, case['preparation_axis_fraction'])
    w = core['weights']
    A = C*(u@w)
    bg = Cb*core['background']*w.sum()
    E = (A-bg)/(A[0]-bg)
    target = case['target_angle_in_P_coordinates_deg']
    target_size = abs(float(rich.wrap_error(target, 0.)))
    closure = 1-abs(rich.wrap_error(heading, target))/target_size
    closure[R < .1] = np.nan
    delta = np.degrees(np.angle(P[1:joff+1]*np.conj(P[:joff])))
    dt = np.diff(cue['time'])
    antipodal = abs(target_size-180.) < 1e-6
    available = (R[:joff] >= .1) & (R[1:joff+1] >= .1) & (not antipodal)
    secant = np.where(available, delta/dt, np.nan)
    rawsecant = delta/dt
    scales = np.asarray(case['fixed_scales']) if spring else np.asarray(case['fixed_scales'])[:, None]
    fulltail = prior.tail_error(off, scales)
    final_u = off['u'][-1]
    peaks, regions = prior.periodic_counts(final_u)
    factor = np.exp(-2j*np.pi*np.fft.fftfreq(n)*float(heading[-1])/360*n)
    rotated = np.fft.ifft(np.fft.fft(u0)*factor).real
    if spring:
        rotated_full=np.r_[rotated,np.fft.ifft(np.fft.fft(initial[n:-1])*factor).real,initial[-1]]
        final_full_error=float(abs(off['final_state'][:-1]-rotated_full[:-1]).max())
    else:
        rotated_full=np.fft.ifft(np.fft.fft(initial,axis=1)*factor[None,:],axis=1).real
        final_full_error=float(abs(off['final_state']-rotated_full).max())
    endpoint_angles_available=min(R[joff],R[-1])>=.1
    amount_ratio=M0/M0[0]
    pool_error=0.
    native_model=None if spring else base.ring.native_model(case['model'])
    for part in (cue,off):
        states=part['full'] if part['full'] is not None else part['snapshots']
        if spring:
            amount=states[:,-1]*(states[:,:n]+states[:,n:-1]).mean(axis=1)
            reference=Cb*(initial[:n]+initial[n:-1]).mean()
            pool_error=max(pool_error,float(abs(amount/reference-1).max()))
        else:
            for group in native_model.conserved_groups:
                reference=Cb*initial[list(group)].sum()/n
                amount=Cb*states[:,list(group)].sum(axis=(1,2))/n
                pool_error=max(pool_error,float(abs(amount/reference-1).max()))
    assert pool_error<1e-7
    row = dict(M0_initial=float(M0[0]), M1_initial_real=float(M1[0].real), M1_initial_imag=float(M1[0].imag),
        M1_initial_magnitude=float(abs(M1[0])), q1_initial=float(q1[0]),
        M0_withdrawal=float(M0[joff]), M0_final=float(M0[-1]), M0_ratio_withdrawal=float(M0[joff]/M0[0]), M0_ratio_final=float(M0[-1]/M0[0]),
        P_withdrawal_real=float(P[joff].real), P_withdrawal_imag=float(P[joff].imag), P_final_real=float(P[-1].real), P_final_imag=float(P[-1].imag),
        R_withdrawal=float(R[joff]), R_final=float(R[-1]), R_sampled_minimum=float(R.min()), R_sampled_minimum_time=float(time[R.argmin()]),
        R_cue_on_sampled_minimum=float(R[:joff+1].min()), q1_withdrawal=float(q1[joff]), q1_final=float(q1[-1]), q2_withdrawal=float(q2[joff]), q2_final=float(q2[-1]),
        angle_withdrawal_raw_deg=float(heading[joff]), angle_final_raw_deg=float(heading[-1]),
        angle_withdrawal_deg=float(heading[joff]) if R[joff]>=.1 else None, angle_final_deg=float(heading[-1]) if R[-1]>=.1 else None,
        target_error_withdrawal_deg=float(rich.wrap_error(heading[joff],target)) if R[joff]>=.1 else None,
        target_error_final_deg=float(rich.wrap_error(heading[-1],target)) if R[-1]>=.1 else None,
        withdrawal_error_closed=float(closure[joff]) if R[joff]>=.1 else None, retained_error_closed=float(closure[-1]) if R[-1]>=.1 else None,
        fraction_requested_turn_final=float(heading[-1]/case['requested_angle_deg']) if R[-1]>=.1 else None,
        fixed_region_excess_withdrawal=float(E[joff]), fixed_region_excess_final=float(E[-1]), fixed_region_excess_sampled_minimum=float(E.min()),
        fixed_region_physical_amount_withdrawal=float(A[joff]), fixed_region_physical_amount_final=float(A[-1]),
        physical_C_withdrawal=float(C[joff]), physical_C_final=float(C[-1]),
        final_resolved_peaks=peaks, final_halfcontrast_regions=regions,
        final_primary_rotated_basal_relative_contrast_error=float(np.max(abs(final_u-rotated))/np.ptp(u0)),
        final_full_field_tail_excursion=fulltail, full_field_tail_gate=None if fulltail is None else bool(fulltail<=1e-5),
        full_field_tail_status='unavailable_from_retained_primary_only_history' if fulltail is None else 'passes_finite_tail_gate' if fulltail<=1e-5 else 'fails_finite_tail_gate',
        final_primary_tail_relative_contrast_excursion=float(np.ptp(off['u'][off['time']>=12000.],axis=0).max()/np.ptp(u0)),
        cue_saved_frames=len(cue['time']), cue_minimum_saved_interval=float(dt.min()), cue_maximum_saved_interval=float(dt.max()),
        cue_weak_saved_frames=int((R[:joff+1]<.1).sum()), cue_secant_intervals=len(dt), cue_secant_available_intervals=int(available.sum()),
        cue_maximum_saved_wrapped_increment_deg=float(abs(delta).max()), antipodal_axis_rate_withheld=bool(antipodal),
        cue_archive_full_state_coverage='full_saved_states' if cue['full'] is not None else 'primary_history_with_full_snapshots',
        cue_native_RHS_continuous_saved_state_availability=cue['full'] is not None,
        conserved_pool_max_relative_error=pool_error,
        saved_initial_frame_max_abs_error=start_error, saved_withdrawal_restart_max_abs_error=restart_error)
    row.update(M1_withdrawal_real=float(M1[joff].real),M1_withdrawal_imag=float(M1[joff].imag),
        M1_final_real=float(M1[-1].real),M1_final_imag=float(M1[-1].imag),
        M0_ratio_sampled_minimum=float(amount_ratio.min()),M0_ratio_sampled_maximum=float(amount_ratio.max()),
        M0_ratio_sampled_minimum_time=float(time[amount_ratio.argmin()]),M0_ratio_sampled_maximum_time=float(time[amount_ratio.argmax()]),
        M0_ratio_cue_on_sampled_minimum=float(amount_ratio[:joff+1].min()),M0_ratio_cue_on_sampled_maximum=float(amount_ratio[:joff+1].max()),
        R_sampled_maximum=float(R.max()),R_release_sampled_minimum=float(R[joff:].min()),
        M0_ratio_at_R_sampled_minimum=float(amount_ratio[R.argmin()]),q2_at_R_sampled_minimum=float(q2[R.argmin()]),
        release_signed_turn_deg=float(rich.wrap_error(heading[-1],heading[joff])) if endpoint_angles_available else None,
        release_target_error_change_deg=float(abs(rich.wrap_error(heading[-1],target))-abs(rich.wrap_error(heading[joff],target))) if endpoint_angles_available else None,
        fixed_region_excess_initial=1.,fixed_region_excess_cue_on_sampled_minimum=float(E[:joff+1].min()),
        fixed_region_excess_cue_on_sampled_minimum_time=float(time[:joff+1][E[:joff+1].argmin()]),
        fixed_region_excess_release_sampled_minimum=float(E[joff:].min()),
        fixed_region_excess_release_sampled_minimum_time=float(time[joff:][E[joff:].argmin()]),
        fixed_region_excess_sampled_minimum_time=float(time[E.argmin()]),
        fixed_region_release_net_change=float(E[-1]-E[joff]),fixed_region_rebound_from_sampled_minimum=float(E[-1]-E.min()),
        fixed_region_current_background_excess_final=float(C[-1]*((u[-1]-core['background'])@w)/(A[0]-bg)),
        physical_C_min=float(C.min()),physical_C_max=float(C.max()),final_full_rotated_basal_absolute_error=final_full_error)
    events = []
    for threshold in (.9,.5,.1):
        for phase_name, keep in (('cue_on',time<=T),('release',time>=T),('entire_history',np.ones(len(time),bool))):
            tt, values = time[keep], E[keep]
            condition = values<=threshold
            entries=np.flatnonzero(condition)
            local_context=dict(sampled_minimum=float(values.min()),sampled_minimum_time=float(tt[values.argmin()]),
                rises_above_threshold_after_first_saved_entry=bool(np.any(~condition[int(entries[0]):])) if len(entries) else False,
                fraction_at_phase_end=float(values[-1]))
            for name, bracket in (('first_saved_entry',rich.first_bracket(tt,condition)),
                                  ('entry_retained_through_all_later_saved_phase_frames',rich.retained_bracket(tt,condition))):
                events.append(dict(condition_id=case['condition_id'],phase=phase_name,threshold=threshold,event=name,
                    phase_start=float(tt[0]),phase_horizon=float(tt[-1]),**local_context,**bracket,
                    phase_local_lower=bracket['lower']-tt[0],phase_local_upper=None if bracket['upper'] is None else bracket['upper']-tt[0]))
    trace = dict(time=time,phase=np.where(time<=T,0,1),M0=M0,M0_ratio=M0/M0[0],M1_real=M1.real,M1_imag=M1.imag,
        P_real=P.real,P_imag=P.imag,R=R,angle_raw_deg=heading,angle_masked_deg=masked,q1=q1,q2=q2,
        fixed_region_transported_background_excess=E,fixed_region_physical_active_amount=A,C=C,prescribed_physical_offset=case['fraction']*C,
        angular_target_error_closed=closure,cue_rate_left=cue['time'][:-1],cue_rate_right=cue['time'][1:],cue_rate_dt=dt,
        cue_rate_delta_deg=delta,cue_rate_R_left=R[:joff],cue_rate_R_right=R[1:joff+1],
        cue_rate_available=available,cue_sampled_angular_secant=secant,cue_raw_wrapped_secant=rawsecant,
        final_material_positions=q,final_active=final_u,final_rotated_basal_active=rotated)
    rhs = selected_rhs(case,cue,M1[0],signal)
    return row, events, trace, rhs


def case_status(case):
    if case['execution']=='exact_archive_reuse':
        return dict(status='complete',**{k:v for k,v in case['reuse'].items() if k.startswith('actual_')})
    folder=RUN/'workers'/f'{case["worker_slot"]:03d}'
    path=folder/'status.json'
    result=read(path) if path.exists() else dict(status='not_started')
    row=dict(status=result['status'])
    if result['status']=='complete':
        for part in result['phases']:
            prefix='actual_cue' if part['kind']=='challenge' else 'actual_recovery'
            row[prefix+'_path']=str(RUN/part['path']);row[prefix+'_sha256']=part['sha256']
    return row


def collect(terminal=False):
    manifest=read(record(RUN/'manifest.json'))
    assert sha(RUN/'manifest.json')==read(record(RUN/'freeze.json'))['manifest_sha256']
    assert len(manifest['cases'])==645
    for rel,digest in manifest['source_sha256'].items():record(RUN/'frozen_sources'/rel,digest)
    for rel,digest in manifest['input_sha256'].items():record(RUN/rel,digest)
    # Native RHS calculations must use the exact frozen native functions.
    for module in (base.ring,spring_api,spring_api.Closure.__module__ and __import__(spring_api.Closure.__module__)):
        path=Path(module.__file__)
        rel=str(path.relative_to(ROOT))
        if rel in manifest['source_sha256']:assert sha(path)==manifest['source_sha256'][rel]
        record(path)
    for rel,digest in manifest['source_sha256'].items():
        if rel.startswith(('code/core/models/','code/core/src/polarity1d/')) or rel in ('code/mechanisms/spring_ring_closure_20261017.py',):
            record(ROOT/rel,digest)
    rows=[];crossings=[];rates=[];rhsrows=[];traces={};profiles=[];contexts=[];seen=set();new=0;cached=0
    for c in manifest['cases']:
        cid=c['condition_id'];status=case_status(c)
        row={k:v for k,v in c.items() if k not in ('reuse','fixed_scales')};row.update(status)
        row['source_kind']=c['execution']
        if c['execution']=='exact_archive_reuse':
            row['reuse_source_summary']=c['reuse'].get('source_summary')
            row['reuse_source_summary_condition_id']=c['reuse'].get('source_summary_condition_id')
        with np.load(RUN/c['input']) as z:initial=z['state'].copy()
        with np.load(RUN/c['cue_input']) as z:q=z['material_positions'].copy();signal=z['initial_cue'].copy()
        u0=initial[:c['cells']] if c['model']=='spring' else initial[0]
        g=geometry(q,u0,c)
        for key,value in g.items():
            if key in c:
                if isinstance(value,str):assert value==c[key]
                else:assert abs(value-c[key])<1e-10
        row.update(g)
        if c['cue_input'] not in seen:
            contexts.append({k:v for k,v in row.items() if k not in ('condition_id','cue_duration','recovery','status','worker_slot') and not k.startswith('actual_')})
            seen.add(c['cue_input'])
        if status['status']!='complete':rows.append(row);continue
        fingerprint=dict(manifest_sha256=sha(RUN/'manifest.json'),collector_sha256=sha(__file__),
            cue_sha256=status['actual_cue_sha256'],recovery_sha256=status['actual_recovery_sha256'])
        cache=dest(f'cache/c{cid}.json','metadata');cache_npz=cache.with_suffix('.npz')
        if cache.exists() and cache_npz.exists() and read(cache)['fingerprint']==fingerprint:
            saved=read(cache)
            with np.load(cache_npz) as z:trace={k:z[k].copy() for k in z.files}
            measured,events,rhs=saved['row'],saved['events'],saved['rhs'];cached+=1
            if terminal:
                record(status['actual_cue_path'],status['actual_cue_sha256']);record(status['actual_recovery_path'],status['actual_recovery_sha256'])
            else:
                SOURCE[status['actual_cue_path']]=status['actual_cue_sha256'];SOURCE[status['actual_recovery_path']]=status['actual_recovery_sha256']
        else:
            cue=phase(record(status['actual_cue_path'],status['actual_cue_sha256']),c['model']=='spring',c['circumference'])
            off=phase(record(status['actual_recovery_path'],status['actual_recovery_sha256']),c['model']=='spring',c['circumference'])
            with np.load(status['actual_cue_path']) as z:
                if c['model']=='spring':
                    key='cue_centre_fraction' if 'cue_centre_fraction' in z else 'cue_center_fraction'
                    assert float(z[key])==c['cue_centre_fraction'] and float(z['cue_peak'])==c['requested_peak']
                else:np.testing.assert_array_equal(z['stimulus'],signal)
            measured,events,trace,rhs=derive(c,cue,off,initial,q,signal)
            temp=cache_npz.with_suffix('.tmp.npz');np.savez_compressed(temp,**trace);temp.replace(cache_npz)
            cache.write_text(json.dumps(dict(fingerprint=fingerprint,row=measured,events=events,rhs=rhs),allow_nan=False)+'\n');new+=1
        row.update(measured);rows.append(row);crossings.extend(events);rhsrows.extend(rhs)
        for key,value in trace.items():traces[f'c{cid}_{key}']=value
        for j in range(len(trace['cue_rate_dt'])):
            available=bool(trace['cue_rate_available'][j])
            rates.append(dict(condition_id=cid,t_left=float(trace['cue_rate_left'][j]),t_right=float(trace['cue_rate_right'][j]),
                delta_time=float(trace['cue_rate_dt'][j]),wrapped_phase_increment_deg=float(trace['cue_rate_delta_deg'][j]),
                R_left=float(trace['cue_rate_R_left'][j]),R_right=float(trace['cue_rate_R_right'][j]),
                sampled_angular_secant_deg_per_native_time=float(trace['cue_sampled_angular_secant'][j]),rate_available=available,
                unavailable_reason='antipodal_axis_rate_withheld' if row['antipodal_axis_rate_withheld'] else '' if available else 'weak_saved_endpoint'))
        for j in range(c['cells']):profiles.append(dict(condition_id=cid,material_position=float(q[j]),
            final_active=float(trace['final_active'][j]),rotated_basal_active=float(trace['final_rotated_basal_active'][j])))
    csvwrite('narrow_arc_summary.csv',rows);csvwrite('narrow_arc_initial_context.csv',contexts)
    csvwrite('narrow_arc_crossings.csv',crossings);csvwrite('narrow_arc_cue_on_sampled_secants.csv',rates)
    csvwrite('narrow_arc_cue_on_RHS_rates.csv',rhsrows);csvwrite('narrow_arc_final_profiles.csv',profiles)
    tracepath=dest('narrow_arc_traces.npz');temp=tracepath.with_suffix('.tmp.npz');np.savez_compressed(temp,**traces);temp.replace(tracepath)
    generated=datetime.now(timezone.utc).isoformat();frame=pd.DataFrame(rows);done=int((frame.status=='complete').sum())
    for module in (base,prior,rich,spring_api):record(module.__file__)
    record(__file__)
    jsonwrite('source_sha256.json',SOURCE)
    classes={str(level):frame.groupby([f'centre_relative_to_core_{level}',f'footprint_relative_to_core_{level}']).size().to_dict() for level in (50,10)}
    classes={level:{'/'.join(k):int(v) for k,v in counts.items()} for level,counts in classes.items()}
    state='finished' if done==645 else 'terminal_native_outputs_incomplete' if terminal else 'waiting_for_native_submission' if not (RUN/'launch.json').exists() else 'native_campaign_in_progress'
    validation=dict(generated_at=generated,requested=645,complete=done,newly_derived=new,cached_replays=cached,
        passed_for_completed_subset=True,all_requested_complete=done==645,classification_counts=classes,
        secant_intervals=len(rates),selected_RHS_records=len(rhsrows),regional_event_records=len(crossings),
        full_field_tail_counts=frame['full_field_tail_status'].value_counts().to_dict() if 'full_field_tail_status' in frame else {})
    jsonwrite('validation.json',validation)
    text=f'''# Narrow-size stationary-cue experiment: derived native readouts

Generated {generated}. Completed {done}/645 requested protocols. Producer state: {state}. All initial overlaps are retained.

[Current summary](supporting/data/narrow_arc_summary.csv) · [Independent centre/footprint classifications](supporting/data/narrow_arc_initial_context.csv) · [Physical vector histories](supporting/data/narrow_arc_traces.npz) · [Actual cue-on secant intervals](supporting/data/narrow_arc_cue_on_sampled_secants.csv) · [Selected exact-state RHS rates](supporting/data/narrow_arc_cue_on_RHS_rates.csv) · [Regional entries](supporting/data/narrow_arc_crossings.csv) · [Final native profiles](supporting/data/narrow_arc_final_profiles.csv) · [Latest status](supporting/metadata/latest_status.json).

Chemical sizes are 30/35/40 with 11 prescribed initial offsets; Spring resting sizes are 2.5/2.75/3 with 10 offsets. Five native pulse durations and 13,000 cue-free units follow each exact basal. Initial cue centre and FULL footprint are classified independently against 50% and 10% contrast edges, including periodic wrap. Distances are prescribed initial geometry, not observed site motion. Spring's stationary material cue has a physical arc proportional to its evolving circumference.

Physical moments use instantaneous circumference times native-cell means. $P=M_1/M_1(0)$, $R=|P|$, $M_0$ is physical active amount and $q_1=|M_1|/M_0$. Raw and R≥0.1-masked headings remain separate. Rates are activity-resultant angular rates, not peak translation or cell velocity. All sampled secants retain actual endpoints/dt and weak-R gaps, with no phase unwrapping. Exact antipodal rates are withheld; near-antipodal placement and weak vectors remain visible context. RHS values use only exact retained full states at onset, withdrawal and time 10 when available; legacy primary-only histories are never supplied invented inactive fields. Full-state availability is separate from this selected-sample export.

Fixed original material-region excess subtracts transported initial background; it does not track site ancestry. First and retained entries are saved-sample brackets, not refined minimum pulse durations. Scalar activity heading, finite withdrawal fate and complete-field tail qualification remain distinct. Primary-only retained histories keep full-field tail status unavailable. Cached completed readouts require the exact manifest/archive/collector hashes; source archives are immutable. Compact outputs replace atomically and latest status hashes that generation. No native integration, job submission or figure rendering is performed here.
'''
    p=OUT/'README.md';temp=p.with_suffix('.tmp.md');temp.write_text(text);temp.replace(p)
    artifacts={str(p.relative_to(OUT)):sha(p) for p in (OUT/'supporting/data').glob('*') if p.is_file()}
    jsonwrite('latest_status.json',dict(**validation,producer=state,status_counts=frame.status.value_counts().to_dict(),artifact_sha256=artifacts))
    print(json.dumps(validation,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,default=RUN);p.add_argument('--out',type=Path,default=OUT);p.add_argument('--terminal',action='store_true');a=p.parse_args()
    RUN,OUT=a.run.resolve(),a.out.resolve();OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'.collector.lock').open('w') as lock:
        fcntl.flock(lock.fileno(),fcntl.LOCK_EX)
        collect(a.terminal)
