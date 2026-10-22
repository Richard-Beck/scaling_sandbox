"""Re-render palette/legend corrections without changing scientific evidence.

The established native-data adapters generate the same plots and measurements.
This wrapper changes categorical colors, places each key beside its own panel
group, and checks exported scientific arrays against the retained revision.
"""
from pathlib import Path
from datetime import datetime,timezone
import copy
import hashlib
import importlib.util
import inspect
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/'figures/length_mechanisms_20261008'
REV=REPORT/'revision_20261014'
OUT=REV/'visual_palette'
ASSETS=OUT/'assets'
BASELINE=REPORT/'archive_20261014_before_spring_dashpot_restore/chapters.json'
CODE=Path(__file__).parent
BLUE,ORANGE,GREEN,PINK,GRAY,INK='#0072B2','#D55E00','#009E73','#CC79A7','#7F7F7F','#202020'
LC={9:'#4477AA',21:'#CC6677',39:'#999933'}
MODEL_COLORS={'goryachev':'#8C564B','wave_pinning':'#008477'}
OVERRIDES={}
AUDITS=[]
SEMANTICS={
 'cue_geometry':{'fixed_physical':BLUE,'fixed_relative':ORANGE,'fixed_slope':'#8B3F82'},
 'bound_MCA_rate':{'profile_advection':BLUE,'flow_gradient':ORANGE,'net_binding':GREEN,
     'uniform_replacement':PINK,'bound_diffusion':GRAY,'net_change':INK},
 'mechanical_velocity_source':{'endpoint_forcing':'#6A3D9A','rho_stress_gradient':'#17BECF',
     'membrane_motion':'#8C564B','observed_total':INK},
 'length':{'9':LC[9],'17':'#AAAAAA','18':'#44AA99','19':'#DDCC77','21':LC[21],
     '23':'#332288','31':'#882255','37':'#88CCEE','39':LC[39]},
 'model':MODEL_COLORS,
 'MCA_face_flux':{'full_bound_advection':'#E69F00','bound_diffusion':GRAY,
     'unbound_advection':'#117733','unbound_diffusion':'#5F9E8D','total_flux':INK},
 'scalar_colormaps':{'Rac':'Oranges','Rho':'Purples','bound_MCA':'Greens',
     'cortical_velocity':'RdBu_r','membrane_tension':'PuOr_r','MCA_drag':'cividis'},
 'legend_policy':'Each key names its panel group and quantity; color names do not define a global meaning across unrelated measurement domains.'
}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):return json.loads(Path(path).read_text())


def write(path,value):Path(path).write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')


def module(name):
    path=CODE/(name+'.py')
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec)
    sys.modules[name]=mod
    spec.loader.exec_module(mod)
    return mod


def line_key(color,label,linestyle='-',marker=None):
    return Line2D([],[],color=color,label=label,ls=linestyle,marker=marker,lw=1.8)


def numeric_lines(fig):
    return [(i,j,np.asarray(line.get_xdata(),dtype=float).copy(),np.asarray(line.get_ydata(),dtype=float).copy())
            for i,ax in enumerate(fig.axes) for j,line in enumerate(ax.lines)]


def verify_numeric(before,fig):
    after=numeric_lines(fig);assert len(before)==len(after)
    for (i,j,x,y),(ii,jj,xx,yy) in zip(before,after):
        assert (i,j)==(ii,jj)
        np.testing.assert_array_equal(x,xx);np.testing.assert_array_equal(y,yy)


def artist_audit(fig):
    curves=[];legends=[]
    for i,ax in enumerate(fig.axes):
        for j,line in enumerate(ax.lines):
            curves.append(dict(axis=i,line=j,label=line.get_label(),color=matplotlib.colors.to_hex(line.get_color()),
                linestyle=line.get_linestyle(),marker=line.get_marker(),
                x_sha256=hashlib.sha256(np.asarray(line.get_xdata(),dtype=float).tobytes()).hexdigest(),
                y_sha256=hashlib.sha256(np.asarray(line.get_ydata(),dtype=float).tobytes()).hexdigest()))
        if ax.get_legend():
            legend=ax.get_legend()
            legends.append(dict(scope='axis '+str(i),entries=[dict(label=t.get_text(),color=matplotlib.colors.to_hex(h.get_color()),
                linestyle=h.get_linestyle()) for t,h in zip(legend.get_texts(),legend.legend_handles) if isinstance(h,Line2D)]))
    for legend in fig.legends:
        legends.append(dict(scope=getattr(legend,'_semantic_scope','figure group'),entries=[dict(label=t.get_text(),
            color=matplotlib.colors.to_hex(h.get_color()),linestyle=h.get_linestyle())
            for t,h in zip(legend.get_texts(),legend.legend_handles) if isinstance(h,Line2D)]))
    return dict(passed=True,numeric_curve_arrays_unchanged=True,curves=curves,legends=legends)


def row_keys(fig,ncols,nrows,keys):
    """Add physical space above a row, preserving every plotted axes' size.

    Twins and reaction-map colorbars follow their original parent row. Keys
    sit above the complete model pair/row, without covering plotted data.
    """
    for legend in list(fig.legends):legend.remove()
    fig.canvas.draw()
    old_width,old_height=fig.get_size_inches()
    base=fig.axes[:ncols*nrows]
    row_of={ax:i//ncols for i,ax in enumerate(base)}
    for ax in fig.axes[len(base):]:
        if hasattr(ax,'_colorbar'):
            row_of[ax]=row_of[ax._colorbar.mappable.axes]
        else:
            siblings=ax.get_shared_x_axes().get_siblings(ax)
            parents=[p for p in siblings if p in row_of]
            assert parents,('Unassigned auxiliary axes',ax)
            row_of[ax]=row_of[parents[0]]
    bounds={ax:ax.get_position().frozen() for ax in fig.axes}
    gap={r:(.70 if len(handles)>5 else .53) for r,(_,handles) in keys.items()}
    for r in gap:
        if any(ax.get_title() for ax in base if row_of[ax]==r):gap[r]+=.28
    new_height=old_height+sum(gap.values())
    fig.set_layout_engine(None)
    fig.set_size_inches(old_width,new_height,forward=True)
    for ax,p in bounds.items():
        row=row_of[ax]
        below=sum(amount for r,amount in gap.items() if r>row)
        ax.set_position([p.x0,(p.y0*old_height+below)/new_height,p.width,p.height*old_height/new_height])
    for r,(scope,handles) in keys.items():
        top=max(bounds[ax].y1*old_height for ax in base if row_of[ax]==r)
        below=sum(amount for row,amount in gap.items() if row>r)
        height=top+below+gap[r]-.04
        columns=3 if len(handles)==6 else min(5,len(handles))
        legend=fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.52,height/new_height),
            ncol=columns,frameon=False,fontsize=10,columnspacing=1.8,handlelength=2.0)
        legend._semantic_scope=scope


def chemistry():
    e=module('length_revision_20261010_chemistry');a=e.a
    e.OUT=OUT;e.ASSETS=ASSETS
    a.FIG=ASSETS;a.DATA=OUT/'chemistry'
    a.DATA.mkdir(parents=True,exist_ok=True)
    a.LC=LC.copy()
    a.MODELS={name:(label,model,MODEL_COLORS[name]) for name,(label,model,_) in a.MODELS.items()}
    baseline=read(BASELINE)
    chapters={c['id']:copy.deepcopy(c) for c in baseline[1:7]}
    chapters['material']['caption']=chapters['material']['caption'].replace('solid blue','solid black',1).replace('dashed orange','dashed brown',1)
    # The second pair refers specifically to C, which now uses blue/coral.
    chapters['material']['caption']=chapters['material']['caption'].replace('peak (dashed orange)','peak (dashed coral)')
    chapters['interface']['caption']=chapters['interface']['caption'].replace(
        'blue mid-pulse in A and 250 after withdrawal in B, orange withdrawal, pink 25 after withdrawal and black the fully relaxed endpoint',
        'blue mid-pulse in A, orange withdrawal, pink 25 after withdrawal, green 250 after withdrawal and dashed black the fully relaxed endpoint')
    e.OUT.mkdir(parents=True,exist_ok=True);ASSETS.mkdir(parents=True,exist_ok=True)
    write(e.OUT/'chapter_overrides.json',chapters)
    length_handles=lambda styles=None:[line_key(color,f'L = {L}',(styles or {}).get(L,'-')) for L,color in LC.items()]
    model_handles=lambda:[line_key(MODEL_COLORS[name],label) for name,(label,_,_) in a.MODELS.items()]

    def save(fig,stem,ignored_caption,methods):
        current=e.current;c=e.chapters[current]
        if current=='growth':
            for ax in fig.axes[:2]:
                for line in list(ax.lines):
                    if len(line.get_xdata())==1 and line.get_marker()=='o':line.remove()
        for legend in list(fig.legends):legend.remove()
        before=numeric_lines(fig)
        if current=='material':
            for i,ax in enumerate(fig.axes[:20]):
                r=i//5
                colors={0:[INK,'#8C564B'],1:[GRAY],2:[BLUE,'#EE7733'],3:[GREEN]}[r]
                for line,color in zip(ax.lines,colors):line.set_color(color)
            keys={0:('A: pole amount',[line_key(INK,'Pole amount with background'),line_key('#8C564B','Pole excess amount','--')]),
                  2:('C: concentration',[line_key(BLUE,'Mean pole concentration'),line_key('#EE7733','Peak pole concentration','--')])}
            row_keys(fig,5,4,keys)
        elif current=='initiation':
            keys={0:('A: spatial profiles',length_handles()+[line_key(INK,'Active u','-'),line_key(INK,'Inactive v','--')]),
                  1:('B: paths in reaction maps',length_handles()+[line_key(INK,'Zero net reaction','--'),line_key(INK,'Incumbent end','none','^'),line_key(INK,'Opposite end','none','o')]),
                  2:('C: spatial reaction rates',length_handles()),
                  3:('D: frozen-substrate reaction',length_handles({9:'-',21:'--',39:':'})+[line_key(INK,'Unstable local threshold','none','o')])}
            row_keys(fig,2,4,keys)
        elif current=='restoration':
            row_keys(fig,2,3,{0:('A: damage histories by length',length_handles()),
                1:('B: reaction histories by length',length_handles()),2:('C–D: models',model_handles())})
        elif current=='sink':
            geometry=[line_key(INK,'Fixed gap','-'),line_key(INK,'Opposite end','--')]
            states=[line_key('#777777','Initial'),line_key(BLUE,'10% core loss'),line_key(ORANGE,'30% core loss')]
            forms=[line_key(INK,'Active u','-'),line_key(INK,'Inactive v','--')]
            row_keys(fig,2,4,{0:('A–B: model and sink position',model_handles()+geometry),
                1:('C–D: model; sink position only in C',model_handles()+
                    [line_key(INK,'C: fixed gap','-'),line_key(INK,'C: opposite end','--')]),
                2:('E: L = 9 state and species',states+forms),3:('F: L = 39 state and species',states+forms)})
        elif current=='interface':
            for ax in fig.axes[2:4]:
                ax.lines[2].set_color(GREEN)
                ax.lines[3].set_color(INK);ax.lines[3].set_linestyle('--')
            row_keys(fig,2,2,{0:('A: cue-on states',[line_key('#777777','Initial'),line_key(BLUE,'Mid-pulse'),line_key(ORANGE,'Withdrawal')]),
                1:('B: recovery states',[line_key(ORANGE,'Withdrawal'),line_key(PINK,'25 after withdrawal'),
                    line_key(GREEN,'250 after withdrawal'),line_key(INK,'10000 after withdrawal','--')])})
        elif current=='growth':
            row_keys(fig,2,2,{0:('A: half-mass displacement',length_handles()),1:('B: half-mass speed',length_handles())})
        else:raise AssertionError(current)
        verify_numeric(before,fig)
        for ax in fig.axes:ax.tick_params(labelsize=10)
        audit=artist_audit(fig)
        ap=ASSETS/(stem+'_curve_legend_audit.json');write(ap,audit)
        files=[ap]
        for ext in ['png','pdf','svg']:
            path=ASSETS/(stem+'.'+ext)
            fig.savefig(path,dpi=220,facecolor='white',bbox_inches='tight',pad_inches=.16)
            files.append(path)
        plt.close(fig)
        cp=ASSETS/(stem+'_caption.txt');cp.write_text(c['caption']+'\n');files.append(cp)
        mp=ASSETS/(stem+'_methods.txt');mp.write_text(methods+'\n\nPalette/legend revision only: every plotted x/y array is checked before and after styling. Local keys apply only to their named panel groups.\n');files.append(mp)
        sources=dict(a.sources);sources[str(Path(__file__).relative_to(ROOT))]=sha(__file__)
        record=dict(generated_at=datetime.now(timezone.utc).isoformat(),producer='finished',passed=True,
            renderer=str(Path(__file__).relative_to(ROOT)),renderer_sha256=sha(__file__),artifacts={p.name:sha(p) for p in files},
            sources=sources,caption=c['caption'],audit=audit,visual_revision_only=True,new_integrations=0,
            reproduction='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B '+str(Path(__file__)))
        manifest=ASSETS/(stem+'.manifest.json');write(manifest,record)
        number=int(c['figure_label'])
        OVERRIDES[str(number)]=dict(image=str((ASSETS/(stem+'.png')).relative_to(REPORT)),
            manifest=str(manifest.relative_to(REPORT)),caption=c['caption'],sources=[[str(mp.relative_to(REPORT)),'Methods'],[str(ap.relative_to(REPORT)),'Curve and local-legend associations']])
        c.update(OVERRIDES[str(number)])
        AUDITS.append(dict(figure=number,curves=len(audit['curves']),legend_groups=len(audit['legends']),numeric_curve_arrays_unchanged=True))
    e.save=save
    e.main()
    baseline_evidence=REPORT/'revision_20261010/editorial/plot_evidence.npz'
    with np.load(e.OUT/'plot_evidence.npz') as current,np.load(baseline_evidence) as old:
        assert set(current.files)==set(old.files)
        for key in current.files:np.testing.assert_array_equal(current[key],old[key])
    AUDITS.append(dict(group='Figures 2–7 evidence arrays',array_identity=True,arrays=len(current.files)))


def figure1():
    m=module('length_revision_20261010_figure1')
    m.OUT=OUT/'figure1';m.OUT.mkdir(parents=True,exist_ok=True)
    render=inspect.getsource(m.render)
    render=render.replace('figsize=(15.3, 6.8)','figsize=(13.2, 7.7)')
    render=render.replace('loc="lower center", bbox_to_anchor=(.52, .018)','loc="upper center", bbox_to_anchor=(.52, .925)')
    render=render.replace('top.plot(x, y, color=color, ls=linestyle, lw=1.9, zorder=3)',
        'top.plot(x, y, color=color, ls=linestyle, lw=1.9, zorder=3, label=cue)')
    exec(compile(render,str(Path(__file__))+':Figure1-legend-layout','exec'),m.__dict__)
    plt.rcParams.update({'font.size':11,'axes.labelsize':11,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10})
    original_save=matplotlib.figure.Figure.savefig
    captured=[]
    def save_hook(fig,*args,**kwargs):
        if not captured:captured.append(artist_audit(fig))
        return original_save(fig,*args,**kwargs)
    matplotlib.figure.Figure.savefig=save_hook
    try:m.main()
    finally:matplotlib.figure.Figure.savefig=original_save
    for path in m.OUT.glob('*.csv'):
        old=REPORT/'revision_20261010/figure1'/path.name
        assert old.exists() and sha(path)==sha(old),('Figure 1 data changed',path)
    ap=m.OUT/'curve_legend_audit.json';write(ap,captured[0])
    meta_path=m.OUT/(m.STEM+'.manifest.json');meta=read(meta_path)
    meta.update(renderer=str(Path(__file__).relative_to(ROOT)),renderer_sha256=sha(__file__),visual_revision_only=True,new_integrations=0)
    meta['sources'][str(Path(__file__).relative_to(ROOT))]=sha(__file__)
    meta['artifacts'][ap.name]=sha(ap)
    meta['reproduction']='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/4473331/.conda/envs/scaling_sandbox/bin/python -B '+str(Path(__file__))
    write(meta_path,meta)
    c=read(BASELINE)[0]
    OVERRIDES['1']=dict(image=str((m.OUT/(m.STEM+'.png')).relative_to(REPORT)),manifest=str(meta_path.relative_to(REPORT)),
        caption=c['caption'],sources=[[str((m.OUT/'chemical_reorganization_length_summary.csv').relative_to(REPORT)),'Chemical-distance measurements'],
        ['../../data/runs/three_test_mechanisms_20261006/results/cell_evidence.json','Matched-peak and Spring feedback controls']])
    AUDITS.append(dict(figure=1,data_exports_byte_identical=True,cue_palette_unchanged=True,status_legend_adjacent_to_A=True))


def collect_sequence_and_conclusion():
    seq=read(OUT/'debelly/chapters.json')
    baseline=read(BASELINE)
    for chapter in seq+[read(REV/'conclusion/chapter.json')]:
        n=int(chapter.get('figure_label',12))
        assert chapter['paragraphs']==baseline[n-1]['paragraphs'],('Scientific prose changed',n)
        OVERRIDES[str(n)]={key:chapter[key] for key in ['image','manifest','caption','sources']}
        meta=read(REPORT/chapter['manifest'])
        AUDITS.append(dict(figure=n,palette_only=True,data_identity_checks=meta.get('audit',{}).get('data_identity_checks',[])))


def main():
    OUT.mkdir(parents=True,exist_ok=True);ASSETS.mkdir(parents=True,exist_ok=True)
    write(OUT/'semantic_palette.json',SEMANTICS)
    figure1();print('Figure 1 palette/layout',flush=True)
    chemistry();print('Figures 2–7 palettes/layout',flush=True)
    collect_sequence_and_conclusion()
    assert set(OVERRIDES)=={str(i) for i in [1,2,3,4,5,6,7,8,10,11,12]}
    for number,override in OVERRIDES.items():
        meta_path=REPORT/override['manifest'];meta=read(meta_path)
        for name,digest in meta.get('artifacts',meta.get('sha256',{})).items():
            path=meta_path.parent/name
            if not path.exists():path=ROOT/name
            assert sha(path)==digest,(number,path)
        for name,digest in meta.get('sources',meta.get('source_sha256',{})).items():
            assert sha(ROOT/name)==digest,(number,name)
    write(OUT/'overrides.json',OVERRIDES)
    write(OUT/'palette_audit.json',dict(passed=True,generated_at=datetime.now(timezone.utc).isoformat(),
        scientific_data_unchanged=True,native_integrations=0,figures=AUDITS,semantic_palette_sha256=sha(OUT/'semantic_palette.json')))
    print(json.dumps(dict(passed=True,overrides=str(OUT/'overrides.json'),figures=sorted(map(int,OVERRIDES))),indent=2))


if __name__=='__main__':main()
