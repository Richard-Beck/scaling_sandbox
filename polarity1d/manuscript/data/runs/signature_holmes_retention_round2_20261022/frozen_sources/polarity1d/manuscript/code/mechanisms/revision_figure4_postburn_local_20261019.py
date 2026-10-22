"""Short second-half native F4 histories plus full condition-medoid atlas.

Uses existing trajectories only. Selection is a declared morphology-descriptor
medoid, with smallest-seed tie breaking; n=2 never establishes typicality.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os

for name in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[4]
MS = ROOT / 'polarity1d/manuscript'
RUN = MS / 'data/runs/revision_noise_strength_broad_20261019'
OUT = MS / 'figures/length_mechanisms_20261008/dev/revision_program_20261019/figure4/postburn_signature_candidate_20261019'
MODELS = dict(wave_pinning='Wave-pinning', goryachev='Goryachev',
              holmes_model3='Holmes 3', spring_ring='Spring ring')
FEATURES = ('mean_R', 'mean_A', 'mean_Q', 'mean_contrast_ratio',
            'qualified_multiregion_fraction', 'weak_R_fraction', 'low_contrast_fraction')
GRID = (.0001, .0003, .001, .003, .01, .03, .1, .3, 1.)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(fig, out, stem):
    fig.savefig(out / (stem + '.png'), dpi=175)
    for ext in ('pdf', 'svg'):
        fig.savefig(out / 'supporting/vectors' / (stem + '.' + ext))


def field(case):
    assert sha(case['profile']) == case['profile_sha256']
    assert sha(case['trace']) == case['trace_sha256']
    with np.load(case['profile']) as z:
        t = z['time']; u = z['primary']; C = np.broadcast_to(z['circumference'], len(t))
        x = z['material_positions'] if 'material_positions' in z else z['x']
        result = dict(time=t.copy(), primary=u.copy(), C=C.copy(), coordinate=x.copy())
    with np.load(case['trace']) as z:
        result['scalar_time'] = z['time'].copy()
        result['R'] = np.hypot(z['P_real'], z['P_imag'])
        result['A'] = z['M0']/z['M0'][0]
        result['Q'] = z['q1']/z['q1'][0]
        result['regions'] = z['current_half_contrast_regions'].copy()
        result['contrast_ratio'] = z['contrast']/z['contrast'][0]
    return result


def selection(cache):
    grouped = {}
    for case_id, case in cache.items():
        descriptor = next(d for d in case['descriptors'] if d['window'] == 'last_half')
        if descriptor['kind'] == 'noise':
            key = (descriptor['model'], descriptor['size'], descriptor['fraction'])
            grouped.setdefault(key, []).append((case_id, descriptor))
    selected = {}; scores = []
    for key, values in sorted(grouped.items()):
        values.sort(key=lambda v: v[1]['seed'])
        X = np.array([[d[k] for k in FEATURES] for _, d in values])
        span = np.ptp(X, axis=0); span[span < 1e-12] = 1.
        distance = np.sum(abs(X[:, None, :]-X[None, :, :])/span, axis=2)
        score = distance.sum(axis=1)
        candidates = np.flatnonzero(np.isclose(score, score.min(), rtol=0, atol=1e-12))
        chosen = int(candidates[0]); selected[key] = values[chosen][0]
        for index, (case_id, d) in enumerate(values):
            scores.append(dict(model=key[0], size=key[1], fraction=key[2], case_id=case_id,
                seed=d['seed'], condition_seeds=len(values), medoid_distance_sum=score[index],
                selected=index == chosen, tie_candidates=len(candidates), **{k:d[k] for k in FEATURES}))
    assert len(selected) == 108
    return selected, scores


def main():
    p = argparse.ArgumentParser(description=__doc__); p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(); out = args.out
    for sub in ('data', 'metadata', 'vectors'):
        (out / 'supporting' / sub).mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'none', 'pdf.fonttype': 42})
    cachepath = RUN / 'results/partial/case_cache.json'
    cache = json.loads(cachepath.read_text()); assert len(cache) == 576
    selected, scores = selection(cache)
    pd.DataFrame(scores).to_csv(out / 'supporting/data/local_medoid_scores.csv', index=False)
    contexts = []; cells = []; plotdata = {}
    fig = plt.figure(figsize=(12.5, 10.4))
    grid = fig.add_gridspec(2, 2, left=.075, right=.94, bottom=.15, top=.84, hspace=.55, wspace=.28)
    fig.suptitle('Local organization during a short second-half window', x=.075, y=.975, ha='left', fontsize=18)
    fig.text(.075, .94, 'Largest original size; f=0.1. Fixed display window begins at half the run; white line marks the native profile shown below.', fontsize=10.6)
    fig.text(.075, .91, 'Morphology-descriptor medoid per condition; each currently has two paths, so both tie and the smaller seed is selected. No typicality claim.', fontsize=10)
    for index, (model, label) in enumerate(MODELS.items()):
        size = 3. if model == 'spring_ring' else 60.
        key = (model, size, .1); case_id = selected[key]; case = cache[case_id]
        descriptor = next(d for d in case['descriptors'] if d['window'] == 'last_half')
        f = field(case); T = float(f['time'][-1]); lo, hi, snapshot = T/2, 5*T/8, 9*T/16
        take = (f['time'] >= lo) & (f['time'] <= hi)
        k, = np.flatnonzero(f['time'] == snapshot)
        j, = np.flatnonzero(f['scalar_time'] == snapshot)
        controls = [v for v in cache.values() if v['descriptors'][0]['model'] == model
            and v['descriptors'][0]['size'] == size and v['descriptors'][0]['fraction'] == .1
            and v['descriptors'][0]['kind'] == 'mean']
        assert len(controls) == 1
        control = controls[0]; mean = field(control)
        mk, = np.flatnonzero(mean['time'] == snapshot)
        np.testing.assert_array_equal(f['primary'][0], mean['primary'][0])
        inner = grid[index//2, index%2].subgridspec(2, 1, height_ratios=(1.65, 1.), hspace=.55)
        imageax = fig.add_subplot(inner[0]); profileax = fig.add_subplot(inner[1])
        extent = 1. if model == 'spring_ring' else size
        maximum = max(float(f['primary'][take].max()), float(f['primary'][0].max()))
        im = imageax.imshow(f['primary'][take].T, origin='lower', aspect='auto',
            interpolation='nearest', extent=[lo, hi, 0, extent], vmin=0, vmax=maximum, cmap='viridis')
        imageax.axvline(snapshot, color='white', lw=1.1, ls='--')
        imageax.set_title(f'{chr(65+index)}  {label} · {"C₀" if model == "spring_ring" else "C"}={size:g}', loc='left', fontsize=11)
        # Native-time tick values and the fixed-window header identify this axis;
        # omit a redundant label to leave the linked profile annotations room.
        imageax.set_xlabel('')
        imageax.set_ylabel('Material q' if model == 'spring_ring' else 'Physical arc', fontsize=9)
        cb = fig.colorbar(im, ax=imageax, fraction=.035, pad=.025)
        cb.ax.tick_params(labelsize=7); cb.set_label('Native activity', fontsize=8)
        profileax.plot(f['coordinate'], f['primary'][0], color='#8d979c', lw=1.25, label='Exact basal')
        profileax.plot(f['coordinate'], mean['primary'][mk], color='#168b84', lw=1.25, ls='--', label='Matched mean control')
        profileax.plot(f['coordinate'], f['primary'][k], color='#b3445d', lw=1.45, label='Noisy snapshot')
        profileax.set(xlim=(0, extent), ylim=(0, maximum*1.04))
        profileax.set_xlabel('Material position q' if model == 'spring_ring' else 'Native physical arc', fontsize=9)
        profileax.set_ylabel('Native activity', fontsize=9)
        profileax.set_title(f't={snapshot:g}: {int(f["regions"][j])} raw regions · contrast/basal={f["contrast_ratio"][j]:.2f}\nR={f["R"][j]:.2f}, A={f["A"][j]:.2f}, Q={f["Q"][j]:.2f}', loc='left', fontsize=8.8)
        for ax in (imageax, profileax):
            ax.tick_params(labelsize=8); ax.spines[['top', 'right']].set_visible(False)
        if index == 0:
            profileax.legend(loc='upper right', fontsize=7, frameon=False)
        if model == 'goryachev':
            zoom = profileax.inset_axes([.37, .32, .59, .55])
            zoom.plot(f['coordinate'], mean['primary'][mk], color='#168b84', lw=1., ls='--')
            zoom.plot(f['coordinate'], f['primary'][k], color='#b3445d', lw=1.)
            zoom.set_xlim(0, extent); zoom.tick_params(labelsize=6)
            zoom.set_title('Native-scale enlargement', fontsize=7)
        context = dict(model=model, size=size, fraction=.1, case_id=case_id, seed=descriptor['seed'],
            display_first_time=lo, display_last_time=hi, actual_snapshot_time=snapshot,
            R=f['R'][j], A=f['A'][j], Q=f['Q'][j], raw_regions=int(f['regions'][j]),
            contrast_ratio=f['contrast_ratio'][j], physical_C=f['C'][k],
            second_half_qualified_multiregion_fraction=descriptor['qualified_multiregion_fraction'],
            trace=case['trace'], trace_sha256=case['trace_sha256'], profile=case['profile'], profile_sha256=case['profile_sha256'],
            mean_profile=control['profile'], mean_profile_sha256=control['profile_sha256'])
        contexts.append(context)
        for q, u, basal, mu in zip(f['coordinate'], f['primary'][k], f['primary'][0], mean['primary'][mk]):
            cells.append(dict(model=model, coordinate=q, primary=u, basal=basal, matched_mean=mu, time=snapshot))
        plotdata[model+'_time'] = f['time'][take]; plotdata[model+'_coordinate'] = f['coordinate']
        plotdata[model+'_primary'] = f['primary'][take]; plotdata[model+'_C'] = f['C'][take]
    fig.text(.075, .09, 'A region is a connected arc above half current contrast; counts need amplitude context. Goryachev’s high-input noisy AND mean-control profiles have low contrast.', fontsize=9.3)
    fig.text(.075, .062, 'Snapshots are fixed in advance within the displayed windows, not selected for extra peaks. Regions do not establish competition or site lineage; no duration classifier.', fontsize=9.3)
    fig.text(.075, .034, 'Native activity; no alignment/peak normalization. Saved fields: chemicals40, Spring4 native units; actual C(t) enters Spring M₀/M₁. Laws/clocks differ.', fontsize=9.3)
    export(fig, out, 'postburn_local_examples'); plt.close(fig)
    pd.DataFrame(contexts).to_csv(out / 'supporting/data/local_example_context.csv', index=False)
    pd.DataFrame(cells).to_csv(out / 'supporting/data/local_example_profile_cells.csv', index=False)
    np.savez_compressed(out / 'supporting/data/local_example_kymograph_data.npz', **plotdata)
    atlas = []
    with PdfPages(out / 'supporting/vectors/all_conditions_second_half_native_atlas.pdf') as pdf:
        for model, label in MODELS.items():
            for size in (1., 2., 3.) if model == 'spring_ring' else (20., 40., 60.):
                fields = {strength: field(cache[selected[(model, size, strength)]]) for strength in GRID}
                maximum = max(float(v['primary'][v['time'] >= v['time'][-1]/2].max()) for v in fields.values())
                fig, axes = plt.subplots(3, 3, figsize=(12.2, 10.5))
                fig.subplots_adjust(left=.075, right=.94, bottom=.14, top=.85, hspace=.65, wspace=.3)
                fig.suptitle(f'{label}, {"C₀" if model == "spring_ring" else "C"}={size:g}: all nine input strengths', x=.075, y=.972, ha='left', fontsize=17)
                fig.text(.075, .935, 'Second-half descriptor medoid per condition, with smallest-seed tie breaking. Common zero-based native activity scale across this page.', fontsize=10)
                fig.text(.075, .905, 'Original0.001/0.01: n=12. All other strengths: n=2 and both paths tie. These are trajectory examples, not frequency estimators.', fontsize=10)
                for ax, strength in zip(axes.flat, GRID):
                    f = fields[strength]; T = float(f['time'][-1]); take = f['time'] >= T/2
                    extent = 1. if model == 'spring_ring' else size
                    image = ax.imshow(f['primary'][take].T, origin='lower', aspect='auto',
                        interpolation='nearest', extent=[T/2, T, 0, extent], vmin=0, vmax=maximum, cmap='viridis')
                    cid = selected[(model, size, strength)]
                    d = next(v for v in cache[cid]['descriptors'] if v['window'] == 'last_half')
                    ax.set_title(f'f={strength:g} · path qualified multi={100*d["qualified_multiregion_fraction"]:.1f}%', loc='left', fontsize=9.5)
                    ax.set_xticks([T/2, 3*T/4, T], [f'{v:g}' for v in (T/2, 3*T/4, T)])
                    ax.set_xlabel('Native time', fontsize=9); ax.set_ylabel('Material q' if model == 'spring_ring' else 'Physical arc', fontsize=9)
                    ax.tick_params(labelsize=8)
                    atlas.append(dict(model=model, size=size, fraction=strength, case_id=cid, seed=d['seed'],
                        display_first_time=T/2, display_last_time=T, color_min=0., color_max=maximum,
                        profile=cache[cid]['profile'], profile_sha256=cache[cid]['profile_sha256']))
                cbar = fig.colorbar(image, ax=list(axes.flat), fraction=.016, pad=.02)
                cbar.ax.tick_params(labelsize=8); cbar.set_label('Native active concentration', fontsize=9)
                fig.text(.075, .085, 'No field alignment, per-history peak normalization, birth/establishment or lifetime classifier. Multi-region counts are qualified at contrast≥10% of own basal.', fontsize=9.2)
                fig.text(.075, .055, 'The first half is discarded as an observed-window choice, without assuming stationarity. Model clocks and input laws remain native and uncalibrated.', fontsize=9.2)
                name = f'local_atlas_{model}_{size:g}'
                fig.savefig(out / 'supporting/vectors' / (name + '.pdf'))
                fig.savefig(out / (name + '.png'), dpi=140)
                pdf.savefig(fig); plt.close(fig)
    pd.DataFrame(atlas).to_csv(out / 'supporting/data/all_condition_native_atlas_selection.csv', index=False)
    result = dict(generated_at=datetime.now(timezone.utc).isoformat(), source=str(Path(__file__).relative_to(ROOT)),
        source_sha256=sha(__file__), cache_sha256=sha(cachepath), cases=576, native_integrations=0,
        main_examples=4, condition_medoid_examples=108, atlas_pages=12,
        medoid_features=list(FEATURES), distance='Sum pairwise L1 distances; each feature divided by within-condition range, scale1 if range<1e-12',
        ties='Scores equal within absolute1e-12; choose smallest seed. All n2 choices tie.',
        fixed_main_window='T/2 to5T/8; snapshot9T/16', profiles_aligned=False,
        amplitude_normalization='None; zero-based native units, shared across nine strengths per model-size atlas page',
        outputs={str(q.relative_to(out)):sha(q) for q in sorted(out.rglob('*')) if q.is_file()
                 and (q.name.startswith('local_') or q.name.startswith('postburn_local') or q.name.startswith('all_condition'))})
    (out / 'supporting/metadata/LOCAL_RENDER_RECEIPT.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
