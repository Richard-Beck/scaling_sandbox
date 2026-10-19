"""Assemble replicated F4 signatures from immutable native-case readouts.

Requires the submitted936-case collection to be complete by default. Joint
whole-seed resampling preserves direction/coverage and full/half-window pairing.
No native integration, stationary-process fit or event classifier.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import collections
import hashlib
import json
import os
import time

for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[key] = '1'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript').is_dir())
MS = ROOT / 'polarity1d/manuscript'
RUN = MS / 'data/runs/revision_noise_strength_replication_scoped_20261019'
OUT = MS / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/signature_assembly_20261019'
MODELS = dict(wave_pinning='Wave-pinning', goryachev='Goryachev', holmes_model3='Holmes 3', spring_ring='Spring ring')
MAIN_LEVELS = (.001, .003, .01, .03, .1)
COLORS = ('#4055a8', '#0d8a93', '#66923d', '#d2882c', '#ba445b')
METRICS = ('multiregion_fraction', 'qualified_multiregion_fraction', 'weak_R_fraction',
           'low_contrast_fraction', 'mean_R', 'mean_A', 'mean_Q', 'mean_contrast_ratio',
           'mean_raw_q1', 'mean_q2', 'mean_paired_AQ')
BAGS = 2000
RNG_SEED = 202610195401


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict):
        return {k: clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write(path, value):
    Path(path).write_text(json.dumps(clean(value), indent=2, allow_nan=False)+'\n')


def snapshot(run, out, draft):
    partial = run / 'results/partial'
    names = ('case_window_descriptors.csv', 'case_window_memory.csv',
             'condition_window_summaries.csv', 'memory_all_conditions_windows_lags.csv', 'worker_inventory.csv')
    for _ in range(12):
        before = (partial / 'status.json').read_bytes()
        data = {name: (partial / name).read_bytes() for name in names}
        after = (partial / 'status.json').read_bytes()
        if before == after and all((partial / name).stat().st_mtime_ns <= (partial / 'status.json').stat().st_mtime_ns for name in names):
            break
        time.sleep(.5)
    else:
        raise RuntimeError('Collector is writing; no consistent snapshot acquired')
    status = json.loads(after)
    if not draft:
        assert status['complete'] and status['producer'] == 'terminal'
        assert status['scientific_cases_available'] == 936 and status['new_completed_outputs'] == 360
    for name, raw in data.items():
        (out / 'supporting/data' / name).write_bytes(raw)
    (out / 'supporting/metadata/native_generation_snapshot.json').write_bytes(after)
    return status, {name: hashlib.sha256(raw).hexdigest() for name, raw in data.items()}


def interval(values, minimum_defined=.95):
    defined = np.isfinite(values)
    fraction = float(defined.mean())
    if fraction < minimum_defined:
        return np.nan, np.nan, fraction
    low, high = np.quantile(values[defined], [.025, .975])
    return float(low), float(high), fraction


def analyze(out):
    data = out / 'supporting/data'
    descriptors = pd.read_csv(data / 'case_window_descriptors.csv')
    memory = pd.read_csv(data / 'case_window_memory.csv')
    conditions = []; correlations = []; paired_windows = []; bags_used = {}
    groupkeys = ['model', 'size', 'fraction', 'kind']
    for ordinal, (key, cases) in enumerate(descriptors.groupby(groupkeys, sort=True)):
        ids = cases.sort_values('seed').case_id.drop_duplicates().tolist()
        n = len(ids); noisy = key[3] == 'noise'
        bag = np.random.default_rng(RNG_SEED + ordinal).integers(0, n, size=(BAGS, n)) if noisy and n >= 3 else None
        bags_used[str(key)] = dict(seeds=n, rng_seed=RNG_SEED+ordinal,
                                  bag_sha256=hashlib.sha256(bag.tobytes()).hexdigest() if bag is not None else None)
        descriptor_boots = {}; memory_boots = {}; memory_points = {}
        for window in ('full', 'first_half', 'last_half'):
            d = cases[cases.window == window].set_index('case_id').loc[ids]
            assert len(d) == n
            row = dict(zip(groupkeys, key), window=window, seeds=n,
                       first_time=float(d.first_time.iloc[0]), last_time=float(d.last_time.iloc[0]))
            for metric in METRICS:
                values = d[metric].to_numpy(float)
                row[metric] = float(values.mean()); row[metric+'_seed_min'] = float(values.min()); row[metric+'_seed_max'] = float(values.max())
                if bag is not None:
                    sample = values[bag].mean(axis=1); low, high, defined = interval(sample)
                    descriptor_boots[(window, metric)] = sample
                else:
                    low = high = np.nan; defined = 0.
                row[metric+'_ci_low'] = low; row[metric+'_ci_high'] = high
                row[metric+'_bootstrap_defined_fraction'] = defined
            conditions.append(row)
            relevant = memory[(memory.model == key[0]) & (memory['size'] == key[1])
                              & (memory.fraction == key[2]) & (memory.kind == key[3]) & (memory.window == window)]
            for lag, lagrows in relevant.groupby('lag', sort=True):
                r = lagrows.set_index('case_id').loc[ids]
                eligible = r.pairs_eligible.to_numpy(float); total = r.pairs_available.to_numpy(float)
                cosine = r.direction_cosine_sum.to_numpy(float)
                count = float(eligible.sum()); available = float(total.sum())
                support = count >= 100 and count/available >= .01
                heading = float(cosine.sum()/count) if count else np.nan
                entry = dict(zip(groupkeys, key), window=window, lag=lag, seeds=n,
                    direction_pairs=int(count), all_pairs=int(available), direction_coverage=count/available,
                    group_support=support, direction_cosine_reported=heading if support else np.nan,
                    seeds_with_pairs=int(np.sum(eligible > 0)), seeds_with_declared_support=int(r.support.sum()))
                energy_a = r.energy_a.to_numpy(float); energy_b = r.energy_b.to_numpy(float)
                product = r.raw_product_sum.to_numpy(float)
                energy = np.sqrt(energy_a.sum()*energy_b.sum())
                entry['raw_vector_correlation'] = float(product.sum()/energy) if energy else np.nan
                boot = {}
                if bag is not None:
                    eb = eligible[bag].sum(axis=1); tb = total[bag].sum(axis=1)
                    cb = cosine[bag].sum(axis=1)
                    valid = (eb >= 100) & (eb/tb >= .01)
                    direction = np.full(BAGS, np.nan); direction[valid] = cb[valid]/eb[valid]
                    boot = dict(direction_cosine_reported=direction, direction_coverage=eb/tb,
                        raw_vector_correlation=product[bag].sum(axis=1)/np.sqrt(energy_a[bag].sum(axis=1)*energy_b[bag].sum(axis=1)))
                for metric in ('direction_cosine_reported', 'direction_coverage', 'raw_vector_correlation'):
                    low, high, defined = interval(boot[metric]) if bag is not None else (np.nan, np.nan, 0.)
                    if metric == 'direction_cosine_reported' and not support:
                        low = high = np.nan
                    entry[metric+'_ci_low'] = low; entry[metric+'_ci_high'] = high
                    entry[metric+'_bootstrap_defined_fraction'] = defined
                    if bag is not None:
                        memory_boots[(window, lag, metric)] = boot[metric]
                    memory_points[(window, lag, metric)] = entry[metric]
                correlations.append(entry)
        if bag is not None:
            for metric in METRICS:
                sample = descriptor_boots[('last_half', metric)]-descriptor_boots[('first_half', metric)]
                low, high, defined = interval(sample)
                paired_windows.append(dict(zip(groupkeys, key), metric=metric, lag=np.nan,
                    contrast='last_half_minus_first_half', estimate=float(cases[cases.window=='last_half'][metric].mean()-cases[cases.window=='first_half'][metric].mean()),
                    ci_low=low, ci_high=high, bootstrap_defined_fraction=defined, seeds=n))
            for lag in sorted(relevant.lag.unique()):
                for metric in ('direction_cosine_reported', 'direction_coverage'):
                    sample = memory_boots[('last_half', lag, metric)]-memory_boots[('first_half', lag, metric)]
                    low, high, defined = interval(sample)
                    paired_windows.append(dict(zip(groupkeys, key), metric=metric, lag=lag,
                        contrast='last_half_minus_first_half', estimate=memory_points[('last_half',lag,metric)]-memory_points[('first_half',lag,metric)],
                        ci_low=low, ci_high=high, bootstrap_defined_fraction=defined, seeds=n))
    df = pd.DataFrame(conditions); mem = pd.DataFrame(correlations)
    df.to_csv(data / 'condition_window_bootstrap.csv', index=False)
    mem.to_csv(data / 'memory_window_bootstrap.csv', index=False)
    pd.DataFrame(paired_windows).to_csv(data / 'paired_halfwindow_bootstrap.csv', index=False)
    write(out / 'supporting/metadata/BOOTSTRAP_PROTOCOL.json',
        dict(whole_seed_resamples=BAGS, base_rng_seed=RNG_SEED, percentile=[.025, .975],
            joint_bags='Same condition bags for all metrics, windows and lags; ordered by seed',
            direction='Ratio of resampled cosine sums to resampled eligible counts; same bags for coverage',
            direction_support='>=100 eligible pairs and>=1% all-pair coverage; >=95% bootstrap support for reported CI',
            descriptors='Equal independent-seed mean of case time means; no per-frame pseudo-replication',
            noisy_n2='Observed seed ranges only; no CI', deterministic_controls='Single paths; no bootstrap replicates',
            groups=bags_used))
    # Exact source reconciliation distinguishes the new uncertainty from a new estimand.
    original = pd.read_csv(data / 'condition_window_summaries.csv')
    joined = df.merge(original, on=groupkeys+['window'], suffixes=('_new', '_source'), validate='one_to_one')
    error = max(float(np.max(abs(joined[k+'_new']-joined[k+'_source']))) for k in METRICS)
    original_mem = pd.read_csv(data / 'memory_all_conditions_windows_lags.csv')
    joined_mem = mem.merge(original_mem, on=groupkeys+['window','lag'], suffixes=('_new','_source'), validate='one_to_one')
    for k in ('direction_cosine_reported','direction_coverage','raw_vector_correlation'):
        delta = abs(joined_mem[k+'_new']-joined_mem[k+'_source'])
        error = max(error, float(delta.max(skipna=True)))
    assert error < 1e-12
    return df, mem, error


def main_figure(out, run, df, mem, draft):
    manifest = json.loads((run / 'manifest.json').read_text())
    Cb = {}
    for basal in manifest['basals']:
        with np.load(run / basal['input']) as z:
            Cb[(basal['model'], float(basal['size']))] = float(z['state'][-1]) if basal['model'] == 'spring_ring' else float(basal['size'])
    fig, axes = plt.subplots(4, 4, figsize=(13.3, 11.5), sharey='row',
                             gridspec_kw={'height_ratios':[1.5,.65,1.2,1.2]})
    fig.subplots_adjust(left=.075, right=.98, bottom=.14, top=.775, hspace=.38, wspace=.25)
    fig.suptitle('Size effects on activity-heading memory depend on input regime', x=.075, y=.974, ha='left', fontsize=17)
    fig.text(.075,.938,'Second-half windows: chemicals16000–32000, Spring1600–3200. First half discarded; stationarity is not assumed.',fontsize=10.3)
    handles=[Line2D([],[],color=c,lw=1.8,marker='o',ms=4,label=f'f={f:g}') for f,c in zip(MAIN_LEVELS,COLORS)]
    fig.legend(handles=handles,loc='upper left',bbox_to_anchor=(.068,.919),ncol=5,frameon=False,fontsize=10,columnspacing=2)
    fig.text(.075,.865,'Draft: actual n in point table; intervals use whole-seed resampling only when n≥3.' if draft else 'Every displayed condition:12 independent seeds. Whiskers:95% whole-seed percentile intervals (2000 paired resamples).',fontsize=9.8)
    points=[]
    for column,(model,name) in enumerate(MODELS.items()):
        for f,color in zip(MAIN_LEVELS,COLORS):
            d=df[(df.model==model)&(df.kind=='noise')&(df.window=='last_half')&(df.fraction==f)].sort_values('size')
            m=mem[(mem.model==model)&(mem.kind=='noise')&(mem.window=='last_half')&(mem.fraction==f)&(mem.lag==(800 if model=='spring_ring' else 8000))].sort_values('size')
            assert len(d)==len(m)==3
            if not draft:assert (d.seeds==12).all() and (m.seeds==12).all()
            x=np.array([Cb[(model,float(v))] for v in d['size']])
            for row,(table,metric) in enumerate(((m,'direction_cosine_reported'),(m,'direction_coverage'),(d,'qualified_multiregion_fraction'),(d,'low_contrast_fraction'))):
                ax=axes[row,column]; y=table[metric].to_numpy(float)
                ax.plot(x,y,color=color,lw=1.55,marker='o',ms=4.5)
                ax.vlines(x,table[metric+'_ci_low'],table[metric+'_ci_high'],color=color,alpha=.65,lw=.95)
                for xi,(_,r) in zip(x,table.iterrows()):
                    points.append(dict(model=model,size=r['size'],basal_circumference=xi,fraction=f,metric=metric,
                        estimate=r[metric],ci_low=r[metric+'_ci_low'],ci_high=r[metric+'_ci_high'],seeds=r.seeds))
        for row,ax in enumerate(axes[:,column]):
            ticks=np.array([Cb[(model,float(v))] for v in ((1.,2.,3.) if model=='spring_ring' else (20.,40.,60.))])
            ax.set_xticks(ticks,[f'{v:.2f}' if model=='spring_ring' else f'{v:g}' for v in ticks])
            ax.set_xlim(ticks[0]-.09*np.ptp(ticks),ticks[-1]+.09*np.ptp(ticks))
            ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',color='#e6eaec',lw=.7);ax.tick_params(labelsize=9)
            ax.set_ylim(-1.05,1.05) if row==0 else ax.set_ylim(-.035,1.05)
            if row==1:ax.set_yticks([0,1])
            elif row>1:ax.set_yticks([0,.5,1])
            if row==3:ax.set_xlabel('Basal physical circumference' if model=='spring_ring' else 'Circumference',fontsize=9.8)
            ax.text(.01,1.04,f'{chr(69+row)}{column+1}',transform=ax.transAxes,fontsize=10.5,ha='left')
        axes[0,column].set_title(f'{name}\nlag={800 if model=="spring_ring" else 8000} native units',loc='left',fontsize=10.8,pad=18)
        axes[0,column].axhline(0,color='#89969b',lw=.65)
    for row,label in enumerate(('Conditional direction cosine','Usable-pair fraction','Contrast-qualified multiple regions','Low-contrast fraction')):
        axes[row,0].set_ylabel(label,fontsize=9.8)
    fig.text(.075,.082,'Direction requires both endpoint R≥0.1; usable coverage retains all available lag pairs. Multiple regions require frame contrast≥10% of basal; low contrast is below that threshold.',fontsize=9.1)
    fig.text(.075,.056,'Positive mean and fluctuation amplitude scale together. Native laws/clocks differ. Activity heading is not cell migration; size effects do not isolate physical site motion from geometry.',fontsize=9.1)
    fig.text(.075,.03,'Read with the native histories and full temporal/window/organization supplements. Spring sampled run-duration inference is excluded from the accepted primary scope.',fontsize=9.1)
    fig.savefig(out/'figure4_signature_readouts.png',dpi=180)
    for ext in ('pdf','svg'):fig.savefig(out/'supporting/vectors'/('figure4_signature_readouts.'+ext))
    plt.close(fig);pd.DataFrame(points).to_csv(out/'supporting/data/main_signature_points.csv',index=False)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run',type=Path,default=RUN)
    p.add_argument('--out',type=Path,default=OUT);p.add_argument('--draft',action='store_true');a=p.parse_args()
    for sub in ('data','metadata','vectors','viewers'):(a.out/'supporting'/sub).mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none','pdf.fonttype':42})
    status,sources=snapshot(a.run,a.out,a.draft);df,mem,error=analyze(a.out)
    main_figure(a.out,a.run,df,mem,a.draft)
    write(a.out/'supporting/metadata/SIGNATURE_ASSEMBLY_RECEIPT.json',dict(
        generated_at=datetime.now(timezone.utc).isoformat(),source=str(Path(__file__).relative_to(ROOT)),
        source_sha256=sha(__file__),native_status=status,source_snapshots=sources,
        maximum_point_difference_from_accepted_collector=error,native_integrations=0,
        draft=a.draft,main_rows=4,main_models=4,main_strengths=list(MAIN_LEVELS)))
    print(json.dumps(dict(out=str(a.out),cases=status['scientific_cases_available'],draft=a.draft,
                         maximum_source_point_difference=error),indent=2))


if __name__=='__main__':main()
