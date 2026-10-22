"""Final S1 forcing and canonical condition/seed census from accepted sources."""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import signature_figure4_interim_render_20261020 as P
import signature_figure4_package_20261020 as Q
import signature_figure4_stage_render_20261020 as S

CONTEXTS=('common','E1_extra','E4_long')


def condition_rows(cases):
    grouped=defaultdict(list)
    for case in cases:
        key=(case['atlas_context'],P.model(case),float(case['C_label']),float(case['f']),case['kind'],float(case['sigma']),float(case['native_horizon']))
        grouped[key].append(case)
    rows=[]
    for key,members in sorted(grouped.items(),key=lambda item:(CONTEXTS.index(item[0][0]),P.MODELS.index(item[0][1]),*item[0][2:])):
        context,model,size,f,kind,sigma,horizon=key
        blocks={row['latent_Gaussian_cohort_id'] for row in members if row['latent_Gaussian_cohort_id'] is not None}
        if kind!='noise' and (blocks or any(row['independent_stochastic_seed_count']!=0 for row in members)):
            raise ValueError('Deterministic controls cannot contribute independent stochastic seeds')
        rows.append(dict(context=context,model=model,C_label=size,f=f,kind=kind,sigma=sigma,native_horizon=horizon,
            native_paths=len(members),within_condition_stochastic_seed_blocks=len(blocks) if kind=='noise' else 0,
            scalar_cadences=json.dumps(sorted({row['scalar_cadence'] for row in members})),
            profile_cadences=json.dumps(sorted({row['profile_cadence'] for row in members})),
            origins=json.dumps(sorted({row['origin'] for row in members})),
            source_cases=json.dumps(members,ensure_ascii=False,sort_keys=True,separators=(',',':')),
            independence_scope='Seed blocks counted within this condition only; sigma branches and retained prefixes remain dependent across conditions'))
    if sum(row['native_paths'] for row in rows)!=len(cases):raise ValueError('A canonical source case was omitted')
    return rows


def load_census(root,digest):
    root=Path(root).resolve();path=root/'coverage_master.json';Q.verified_file(dict(path=str(path),sha256=digest),root)
    master=json.loads(path.read_text())
    if master.get('passed') is not True or master.get('all_required_contexts_complete') is not True:
        raise ValueError('S1 requires the complete independently accepted canonical population')
    if set(master['source_case_inventory_hashes'])!=set(CONTEXTS):raise ValueError('Canonical inventory contexts differ')
    sources=dict(canonical_population_master=dict(path=str(path),sha256=digest));cases=[]
    for context in CONTEXTS:
        path=root/'canonical_inputs'/context/'inventory.json';digest=master['source_case_inventory_hashes'][context]
        Q.verified_file(dict(path=str(path),sha256=digest),root);inventory=json.loads(path.read_text())
        if inventory['atlas_context']!=context or len(inventory['cases'])!=inventory['source_case_count']:
            raise ValueError('A canonical whole-population inventory changed')
        if any(row['atlas_context']!=context for row in inventory['cases']):raise ValueError('A case changed canonical context')
        cases.extend(inventory['cases']);sources[context]=dict(path=str(path),sha256=digest)
    identity=[row['physical_run_id'] for row in cases]
    measured=dict(physical_runs=len(set(identity)),native_observations=len({row['source_observation_id'] for row in cases}),
        noisy_physical_runs=sum(row['kind']=='noise' for row in cases),deterministic_control_runs=sum(row['kind']!='noise' for row in cases),
        unique_latent_Gaussian_cohorts=len({row['latent_Gaussian_cohort_id'] for row in cases if row['latent_Gaussian_cohort_id'] is not None}),
        unique_latent_seed_lineages=len({row['latent_seed_lineage_id'] for row in cases if row['latent_seed_lineage_id'] is not None}))
    if len(identity)!=len(set(identity)) or any(master['independence_and_reuse_census'][key]!=value for key,value in measured.items()):
        raise ValueError('S1 physical paths, controls and dependency census must match the accepted crosswalk')
    return cases,condition_rows(cases),master,sources


def render(root,digest,output,name,forcing_brief,forcing_digest):
    cases,conditions,master,sources=load_census(root,digest)
    brief=Q.verified_file(dict(path=str(forcing_brief),sha256=forcing_digest),Path.cwd())
    sources['forcing_and_design_authority']=dict(path=str(brief),sha256=forcing_digest)
    sources['design_display_source']=dict(path=str(Path(__file__).resolve()),sha256=P.sha(__file__))
    book=S.Book(output,name,'Final canonical design and forcing',sources,section='S1')
    try:
        definitions=[dict(quantity='Input law',definition='I(q,t)=f I_ref exp(σ Z(q,t)−σ²/2)',authority_sha256=forcing_digest),
            dict(quantity='Heading qualification',definition='R(t)=|M1(t)|/|M1,basal|; exact no-input native basal reference, not |M1|/M0',authority_sha256=forcing_digest),
            dict(quantity='Chemical forcing',definition='Gaussian OU knots; native relaxation50; physical kernel1.5',authority_sha256=forcing_digest),
            dict(quantity='Spring forcing',definition='Native relaxation5; material-attached correlation ell_phys(t)=0.15 C(t)/C_basal',authority_sha256=forcing_digest),
            dict(quantity='Controls and uncertainty',definition='σ0 matched-mean and zero-input controls are deterministic, with0 independent stochastic seeds and no seed bootstrap. Noisy uncertainty resamples whole seeds.',authority_sha256=forcing_digest),
            dict(quantity='Native clocks and optics',definition='Actual horizons/cadences come from each case. Physical optical Gaussian is fixed1.5, or0.15 for Spring; primary imaging cadence80, or8 for Spring.',authority_sha256=forcing_digest)]
        fig=plt.figure(figsize=(11.7,8.3))
        fig.text(.5,.83,r'$I(q,t)=f I_{\mathrm{ref}}\exp\!\left[\sigma Z(q,t)-\frac{\sigma^2}{2}\right]$',ha='center',fontsize=19)
        fig.text(.96,.83,'(S1)',ha='right',fontsize=11)
        fig.text(.5,.70,r'$R(t)=\frac{|M_1(t)|}{|M_{1,\mathrm{basal}}|}$',ha='center',fontsize=19)
        fig.text(.96,.70,'(S2)',ha='right',fontsize=11)
        body=['R uses the exact no-input native basal first-mode magnitude; it can exceed1.',
              'Chemical Gaussian forcing: physical kernel1.5; native OU relaxation50.',
              'Spring forcing is material attached: physical correlation length=0.15 C(t)/C_basal; relaxation5.',
              'Fixed physical optical smoothing uses1.5 (Spring0.15), separately from the forcing coordinate.',
              'Matched mean σ=0 and zero-input controls have no independent stochastic-seed bootstrap.',
              'Native clocks, grids, source hashes, σ branches and retained-prefix dependencies remain in the census.',
              'Common and longest observation arms select from these same sources; neither adds a duplicate cohort.',
              'One-dimensional activity-heading persistence is distinct from empirical migration persistence.']
        fig.text(.075,.56,'\n\n'.join(body),va='top',fontsize=10.3,linespacing=1.3)
        book.save(fig,'S1_forcing','Native forcing, reporter qualification and observation coordinates',definitions,
            'f uses each model’s native reference input; native kinetic/time units are not ranked across models. '
            'The original R convention is separate from Figure3 instantaneous |M1|/M0 coherence. '
            'Exact E4 prefix repairs and the unresolved historical Spring scalar-cadence sensitivity remain separate S12 qualifications.')
        summary=[]
        for context in CONTEXTS:
            for model in P.MODELS:
                selected=[row for row in cases if row['atlas_context']==context and P.model(row)==model]
                counts=Counter(row['kind'] for row in selected)
                summary.append(dict(context=context,model=model,native_paths=len(selected),noise_paths=counts['noise'],mean_controls=counts['mean'],zero_controls=counts['zero'],
                    unique_Gaussian_cohorts=len({row['latent_Gaussian_cohort_id'] for row in selected if row['latent_Gaussian_cohort_id'] is not None}),
                    unique_seed_lineages=len({row['latent_seed_lineage_id'] for row in selected if row['latent_seed_lineage_id'] is not None}),
                    accepted_global_census=json.dumps(master['independence_and_reuse_census'],sort_keys=True)))
        def table_page(rows,cells,headers,label,title,note,widths):
            fig,ax=plt.subplots(figsize=(11.7,8.3));ax.axis('off')
            table=ax.table(cellText=cells,colLabels=headers,colWidths=widths,cellLoc='left',bbox=[0,0,1,1])
            table.auto_set_font_size(False);table.set_fontsize(7.5)
            for (i,j),cell in table.get_celld().items():
                cell.PAD=.035
                if i==0:cell.set_facecolor('#edf1f4');cell.set_text_props(weight='bold')
            fig.subplots_adjust(left=.055,right=.97,top=.85,bottom=.24)
            fig.canvas.draw();renderer=fig.canvas.get_renderer()
            for cell in table.get_celld().values():
                box=cell.get_window_extent(renderer);text=cell.get_text().get_window_extent(renderer)
                if not (text.x0>=box.x0 and text.x1<=box.x1 and text.y0>=box.y0 and text.y1<=box.y1):raise ValueError('S1 table clips a source value')
            book.save(fig,label,title,rows,note)
        table_page(summary,[[row['context'],P.LABELS[row['model']],row['noise_paths'],row['mean_controls'],row['zero_controls'],row['native_paths'],row['unique_Gaussian_cohorts'],row['unique_seed_lineages']] for row in summary],
            ['Context','Model','Noisy paths','Mean controls','Zero controls','Native paths','Gaussian\ncohorts','Seed\nlineages'],'S1_population','All canonical paths and their dependency census',
            'Counts are native paths, not a pooled independent-seed count. Sigma branches share Gaussian cohorts; E4 prefixes retain lineage dependence. '
            'The accepted global census is3497 paths,3228 noisy paths,269 deterministic controls,2364 Gaussian cohorts and1990 seed lineages. '
            'E4 includes374 retained-prefix noisy paths and250 independent top-ups. Unaccepted attempted prefixes remain in the engineering archive, separately documented in S12.',[.12,.16,.12,.12,.12,.12,.12,.12])
        def cadence(row):
            return '/'.join(','.join(f'{float(value):g}' for value in json.loads(row[key])) for key in ('scalar_cadences','profile_cadences'))
        for start in range(0,len(conditions),18):
            selected=conditions[start:start+18]
            table_page(selected,[[row['context'],P.LABELS[row['model']],f"{row['C_label']:g}",f"{row['f']:g}",f"{row['sigma']:g}",row['kind'],
                str(row['native_paths'])+' / '+str(row['within_condition_stochastic_seed_blocks']),f"{row['native_horizon']:g}",cadence(row),','.join(json.loads(row['origins']))] for row in selected],
                ['Context','Model','C₀','f','σ','Kind','Paths /\nseed blocks','Native H','Scalar / profile Δt','Origin'],
                f'S1_conditions_{start//18+1:02d}','Complete condition, seed and native observation inventory',
                'Seed blocks are counted within each condition only; they must not be summed across sigma branches, inputs or retained-prefix comparisons. '
                'Controls contribute0 independent stochastic seeds. Cadence lists are exact, including mixed native sources. '
                'Each page CSV retains every member’s complete canonical source record, archive hashes, cohort/lineage identity and prefix/top-up provenance. Native units differ across models.',
                [.085,.135,.055,.075,.05,.075,.10,.105,.17,.15])
        return book.close()
    except Exception:
        book.pdf.close();raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('atlas-root','output','forcing-brief'):parser.add_argument('--'+key,type=Path,required=True)
    for key in ('coverage-sha256','forcing-sha256','name'):parser.add_argument('--'+key,required=True)
    args=parser.parse_args()
    print(json.dumps(render(args.atlas_root,args.coverage_sha256,args.output,args.name,args.forcing_brief,args.forcing_sha256)['pdf'],indent=2))
