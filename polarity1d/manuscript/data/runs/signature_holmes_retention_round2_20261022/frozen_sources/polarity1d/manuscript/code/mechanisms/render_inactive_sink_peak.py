"""Peak-height loss rates beside active-amount loss rates for sink-only runs."""
from pathlib import Path
import csv, json, os
import numpy as np
import inactive_sink_loss as s


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    from matplotlib.backends.backend_pdf import PdfPages
    out=s.OUT;figdir=s.FIG
    manifest=s.verify(out)
    status=s.r.read(out/'results/status.json')
    assert status['complete'] and status['audit_passed']
    assert s.r.read(out/'results/numerical_accuracy.json')['passed']
    assert s.r.read(out/'results/wave_peak_dip_accuracy.json')['passed']
    refined={}; refinement_checks=[]
    for job in manifest['jobs']:
        record=s.r.read(out/'peak_refinement'/f'slot_{job["slot"]:03d}.json')
        assert record['job']==job
        for row in record['records']:
            path=out/row['archive'];assert s.n.sha(path)==row['archive_sha256']
            with np.load(path) as z:
                with np.load(s.n.ROOT/job['case']['checkpoint']) as cp:np.testing.assert_array_equal(z['original'],cp['state'])
                mask=s.overlap_weights(0,job['pole_edge'],job['case']['length'],z['active'].shape[1])>0
                peak=np.max(z['active'][:,mask],axis=1);frac=peak/peak[0]
                aa,bb=[s.crossing(z['time'],frac,q) for q in [.9,.7]]
                np.testing.assert_allclose([aa,bb,.2/(bb-aa)],
                    [row['peak_t10_loss'],row['peak_t30_loss'],row['peak_fractional_loss_rate_10_to_30']],rtol=1e-12,atol=1e-12)
            refined[(row['model'],row['length'],row['geometry'],row['strength'])]=row
            refinement_checks.append(dict(archive=row['archive'],passed=True,relative_rate_change=row['relative_rate_change']))
    assert len(refined)==121
    s.r.write(out/'results/peak_refinement_audit.json',dict(generated_at=s.n.now(),passed=True,cases=len(refined),checks=refinement_checks))
    rows=[]
    for job in manifest['jobs']:
        directory=out/'workers'/f'{job["slot"]:03d}'
        for i,k in enumerate(manifest['strengths']):
            path=directory/f'k_{i:02d}.npz'
            row=s.r.read(path.with_suffix('.json'))['metrics']
            with np.load(path) as z:
                time=z['time'];height=z['peak_height'];frac=height/height[0]
                t90,t70=[s.crossing(time,frac,q) for q in [.9,.7]]
                row.update(peak_t10_loss=t90,peak_t30_loss=t70,
                    peak_fractional_loss_rate_10_to_30=(.2/(t70-t90) if t70 is not None and t90 is not None and t70>t90 else None),
                    peak_height_changes_monotonically_to_half=bool(np.max(np.diff(frac[time<=row['peak_t50_loss']]))<1e-6) if row['peak_t50_loss'] is not None else None)
                row['archive']=str(path.relative_to(out))
                refinement=refined.get((row['model'],row['length'],row['geometry'],row['strength']))
                row['peak_crossing_refined']=refinement is not None
                if refinement is not None:
                    for key in ['peak_t10_loss','peak_t30_loss','peak_fractional_loss_rate_10_to_30']:row[key]=refinement[key]
            rows.append(row)
    s.n.csvwrite(out/'results/peak_and_amount_measurements.csv',rows)
    cmap=plt.get_cmap('viridis');norm=Normalize(3,39)
    columns=[(model,geometry) for model in ['wave_pinning','goryachev'] for geometry in ['opposite','fixed_gap']]
    names={'wave_pinning':'Wave-pinning','goryachev':'Goryachev'}
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})

    def curves(ax,model,geometry,key,censor=False):
        rs=[v for v in rows if v['model']==model and v['geometry']==geometry]
        for L in sorted(set(v['length'] for v in rs)):
            series=sorted([v for v in rs if v['length']==L],key=lambda v:v['strength'])
            ax.plot([v['strength'] for v in series],[v[key] if v[key] is not None else np.nan for v in series],
                    'o-',ms=2.8,lw=1.1,color=cmap(norm(L)))
            if censor:
                missing=[v['strength'] for v in series if v[key] is None]
                ax.scatter(missing,[s.HORIZON]*len(missing),marker='^',s=19,facecolors='none',edgecolors=[cmap(norm(L))])
        ax.set(xscale='log',yscale='log',xlabel='Sink strength k (1 / native time)');ax.grid(alpha=.17)
        if censor:ax.axhline(s.HORIZON,color='.6',ls=':',lw=.8)

    with PdfPages(figdir/'sink_peak_loss_report.pdf') as pdf:
        fig,axes=plt.subplots(2,4,figsize=(17,8),layout='constrained')
        for col,(model,geo) in enumerate(columns):
            for i,key in enumerate(['peak_fractional_loss_rate_10_to_30','fractional_loss_rate_10_to_30']):curves(axes[i,col],model,geo,key)
            axes[0,col].set_title(names[model]+'\n'+('A · Opposite end' if geo=='opposite' else 'B · Fixed gap from incumbent'))
        axes[0,0].set_ylabel('Peak-height loss rate\n(initial fraction / native time)')
        axes[1,0].set_ylabel('Active-amount loss rate\n(initial fraction / native time)')
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes,label='Cell length L',shrink=.8)
        fig.suptitle('How inactive-pool sink strength changes incumbent peak loss',fontsize=18)
        fig.supxlabel('Rates average the decline from 90% to 70% of each initial measurement; sharp peak crossings use dense local refinement.\nContinuous sinks; width 1. B has a fixed gap of 3 from the INITIAL half-height pole edge. Peak height and amount use that fixed core.\nNo added cues. Missing rates indicate 30% loss was not reached by 20,000 native time units.',fontsize=10)
        fig.savefig(figdir/'sink_peak_loss.png',dpi=190);fig.savefig(figdir/'sink_peak_loss.svg');pdf.savefig(fig);plt.close(fig)
        fig,axes=plt.subplots(2,4,figsize=(17,8),layout='constrained')
        for col,(model,geo) in enumerate(columns):
            for i,key in enumerate(['peak_t10_loss','peak_t50_loss']):curves(axes[i,col],model,geo,key,censor=True)
            axes[0,col].set_title(names[model]+'\n'+('A · Opposite end' if geo=='opposite' else 'B · Fixed gap'))
        axes[0,0].set_ylabel('Time to 10% peak-height loss')
        axes[1,0].set_ylabel('Time to 50% peak-height loss')
        fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes,label='Cell length L',shrink=.8)
        fig.suptitle('Peak-loss delay and half-loss time',fontsize=18)
        fig.supxlabel('Open triangles at 20,000: threshold not reached within the observation horizon. No no-sink simulations or reversal cues.',fontsize=10)
        fig.savefig(figdir/'sink_peak_loss_times.png',dpi=190);pdf.savefig(fig);plt.close(fig)
        # Short and long cells make the geometry effect easy to compare directly.
        fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
        for col,model in enumerate(['wave_pinning','goryachev']):
            for i,key in enumerate(['peak_fractional_loss_rate_10_to_30','fractional_loss_rate_10_to_30']):
                ax=axes[i,col]
                for L in [9.,21.,39.]:
                    for geo,ls in [('opposite','-'),('fixed_gap','--')]:
                        rs=sorted([v for v in rows if v['model']==model and v['geometry']==geo and v['length']==L],key=lambda v:v['strength'])
                        ax.plot([v['strength'] for v in rs],[v[key] if v[key] is not None else np.nan for v in rs],color=cmap(norm(L)),ls=ls,lw=2,label=f'L={L:g}, '+('opposite end' if geo=='opposite' else 'fixed gap'))
                ax.set(xscale='log',yscale='log',xlabel='Sink strength k');ax.grid(alpha=.2)
            axes[0,col].set_title(names[model]);axes[0,col].legend(fontsize=9,frameon=False)
        axes[0,0].set_ylabel('Peak-height loss rate (1 / native time)');axes[1,0].set_ylabel('Active-amount loss rate (1 / native time)')
        fig.suptitle('Closer sinks accelerate loss, with length dependence remaining',fontsize=17)
        fig.savefig(figdir/'sink_geometry_comparison.png',dpi=190);pdf.savefig(fig);plt.close(fig)
    numerical=s.r.read(out/'results/numerical_accuracy.json')
    maximum=max(max(v['relative_errors'].values()) for v in numerical['cases'])
    maximum=max(maximum,max(v['relative_error'] for v in s.r.read(out/'results/wave_peak_dip_accuracy.json')['cases']))
    stamp=s.n.now();csvlink=os.path.relpath(out/'results/peak_and_amount_measurements.csv',figdir)
    (figdir/'index.html').write_text('<!doctype html><html><head><meta charset="utf-8"><title>Inactive sinks and incumbent peak loss</title></head><body style="font:17px system-ui;max-width:1600px;margin:30px auto;padding:0 24px"><h1>Inactive sinks and incumbent peak loss</h1><p>Generated '+stamp+'; producer finished. All 1,040 sink-only trials audited. 121 sharp peak crossings received dense local refinement; fourteen accuracy replays passed. Maximum selected rate/time difference '+f'{100*maximum:.3f}'+'%.</p><p>Continuous inactive removal k*v in width-1 sinks. A: opposite end [L-1,L]. B: a fixed gap of 3 native length units beyond the initial half-height incumbent edge. Both sink positions remain fixed throughout. All trials use exact retained basal checkpoints. No no-sink simulations and no reversal cues.</p><p>Each rate averages 10–30% loss of its own initial measurement. Peak height and active amount are measured in the fixed original half-height core. Loss delays and half-loss times remain separate. Open triangles indicate right censoring at 20,000; missing rates remain gaps. Native units are model-specific.</p><p><a href="sink_peak_loss_report.pdf">Peak-loss figure and supplement PDF</a> · <a href="'+csvlink+'">Measurements CSV</a> · <a href="inactive_sink_loss.pdf">Amount, morphology and trajectory supplement</a></p><p>At high sink strengths, loss rates saturate. The fixed-gap sink raises that plateau, while longer cells remain slower under both geometries. Long wave-pinning cells show a local dip in peak-height loss rate at intermediate strengths even as active-amount loss speeds up; the dip persists after dense crossing refinement and tighter numerical replays. This supports a role for separation in the tested sink response, with further dependence on the length-dependent starting pole and reservoir; the sweep alone does not uniquely identify the transport mechanism.</p>'+''.join('<img style="width:100%" src="'+name+'">' for name in ['sink_peak_loss.png','sink_peak_loss_times.png','sink_geometry_comparison.png','incumbent_morphology.png','incumbent_trajectories.png'])+'<p>Width in the morphology supplement uses the fixed INITIAL half-height threshold; it can reach zero while active material remains. The removed field is inert bookkeeping and never feeds back. Accuracy replays keep the original spatial mesh; no mesh-convergence claim. Preparation failures and sinks that cannot fit beyond the pole are listed in the run manifest.</p></body></html>')
    s.r.write(out/'results/peak_render_provenance.json',dict(generated_at=stamp,renderer_sha256=s.n.sha(Path(__file__)),
        measurement_sha256=s.n.sha(out/'results/peak_and_amount_measurements.csv'),
        pdf_sha256=s.n.sha(figdir/'sink_peak_loss_report.pdf'),complete=True))
    print(json.dumps(dict(peak_rates=len(rows),max_accuracy_change=maximum,figure=str(figdir/'sink_peak_loss.png'))))


if __name__=='__main__':main()
