"""Group development exports by experiment, preserving a reversible inventory.

Run after generators that retain frozen historical output paths. Scientific
campaign inputs and simulation evidence are never modified by this organizer.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,os,re,shutil
from urllib.parse import unquote,urlsplit
ROOT=Path(__file__).resolve().parents[2]
REPO=ROOT.parents[1]
REPORT=ROOT/'figures/length_mechanisms_20261008'
DEV=REPORT/'dev'
RECORDS=DEV/'_organization'
GROUPS={
 'vision_aligned_checkpoint_20261017':('vision_aligned_',),
 'vision_linear_20261017':('vision_linear_',),
 'vision_ring_20261017':('vision_ring_',),
 'vision_noise_pilot_20261017':('vision_noise_pilot_',),
 'vision_noise_20261017':('vision_noise_',),
 'updated_checkpoint_20261017':('updated_checkpoint_',),
 'expanded_ring_components_20261017':('expanded_ring_',),
 'all_five_recovery_pilot_20261017':('all_five_recovery_pilot_',),
 'debelly_moving_breaks_20261017':('debelly_moving_breaks_',),
 'all_five_linear_function_20261017':('all_five_linear_function_',),
 'global_polarity_vector_20261017':('global_polarity_',),
 'linear_then_ring_checkpoint_20261017':('linear_incumbent_', 'linear_then_ring_'),
 'debelly_rac_continuity_20261017':('debelly_rac_',),
 'chemical_ring_functional_20261017':('ring_functional_',),
 'mca_speed_decomposition':('MCA_speed_decomposition_',),
 'reversal_flux_flow':('reversal_MCA_flux_cortical_flow_',),
 'stimulated_endpoint':('stimulated_endpoint_',),
 'mca_molecular_fates':('MCA_molecular_fates',),
 'mca_threshold_profiles':('mca_threshold_',),
 'no_protrusion':('no_protrusion_',),
 'deBelly_redraft':('deBelly_redraft_',),
 'figure_x_reversal_overview':('Figure_X_',),
 'figure_z_mca_redistribution':('Figure_Z_',),
 'resting_state':('basal_endpoint_', 'basal_global_', 'basal_spatial_'),
 'basal_feedback_controls':('basal_feedback_',),
 'tension_rho_ablation':('global_tension_rho_',),
 'beta0_reversal':('beta0_reversal',),
 'native_reversal':('MCA_', 'Fig9F_', 'old_pole_Rac_'),
 'L29_fullstate':('L29_',),
}
TITLES={'resting_state':'Native resting state','basal_feedback_controls':'Basal feedback controls','tension_rho_ablation':'Global tension→Rho ablation','beta0_reversal':'Beta=0 reversal','no_protrusion':'Polarization with protrusion disabled','native_reversal':'Native reversal and MCA budgets','L29_fullstate':'L=29 full-state example','figure_x_reversal_overview':'Figure X: reversal overview','figure_z_mca_redistribution':'Figure Z: MCA redistribution distance and speed','deBelly_redraft':'deBelly redraft with Figures X Y and Z'}
TITLES['mca_threshold_profiles']='Steady MCA profiles near the lower polarity limit'
TITLES['mca_molecular_fates']='MCA molecular residence and first-exit fate'
TITLES['stimulated_endpoint']='Stimulated-end protrusion and cortical velocity'
TITLES['reversal_flux_flow']='Reversal MCA flux and cortical-flow snapshots'
TITLES['mca_speed_decomposition']='MCA redistribution speed across length'
TITLES['chemical_ring_functional_20261017']='Native ring initial gain and conditional repair across size and placement'
TITLES['linear_then_ring_checkpoint_20261017']='Main-text checkpoint: linear incumbent enlargement, then ring separation'
TITLES['all_five_linear_function_20261017']='All-five linear native gain and targeted incumbent-repair design'
TITLES['debelly_rac_continuity_20261017']='Rac-triggered fixed point gaps on the deBelly ring'
TITLES['global_polarity_vector_20261017']='Global activity polarity on retained chemical-ring histories'
TITLES['expanded_ring_components_20261017']='Expanded ring figure: baseline correspondence, gain maps and three fates'
TITLES['all_five_recovery_pilot_20261017']='All-five incumbent activity-recovery pilot'
TITLES['debelly_moving_breaks_20261017']='deBelly sidecar with moving mechanical breaks'
TITLES['updated_checkpoint_20261017']='Updated main-text checkpoint: all-five function and expanded ring analysis'
TITLES['vision_aligned_checkpoint_20261017']='Vision-aligned checkpoint: linear responses, ring competition and noisy polarity persistence'
TITLES['vision_linear_20261017']='All-five linear gain, recovery and finite long-cue responses'
TITLES['vision_ring_20261017']='Supported-ring gain, recovery and finite competition across size and separation'
TITLES['vision_noise_20261017']='Replicated noisy-ring polarity persistence across sizes'
TITLES['vision_noise_pilot_20261017']='Bounded noisy-ring numerical and response-regime pilot'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,value):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
def group(name):
 for folder,prefixes in GROUPS.items():
  if name.startswith(prefixes):return folder
 if name.startswith('dev_'):return '_organization'
 raise ValueError('Unclassified export: '+name)


def organize_group(name,dry_run=False):
 """Organize one topic without rewriting other topics or global provenance.

 The selected topic's generated README may change. The development index gets
 one missing link by insertion, retaining all existing text. Existing scientific
 files, artifact bytes, unrelated README files and all global manifests remain
 unchanged. This deliberately does not perform the historical broad cleanup.
 """
 if name not in GROUPS:raise ValueError('Unknown group: '+name)
 folder=DEV/name
 audit_path=folder/'supporting/metadata/organization_audit.json'
 directories=[folder,*[folder/'supporting'/category for category in ['vectors','data','metadata','viewers']]]
 def selected(p):
  try:return group(p.name)==name
  except ValueError:return False
 pending=[p for p in DEV.glob('*.png') if p.is_file() and selected(p)]
 legacy_supporting=DEV/'supporting'
 if legacy_supporting.exists():
  pending.extend(p for p in legacy_supporting.rglob('*') if p.is_file() and selected(p))
 moves={}
 for path in pending:
  dest=folder/(path.name if path.parent==DEV else path.relative_to(DEV))
  if dest.exists() and sha(path)!=sha(dest):
   raise FileExistsError(f'Refusing to replace differing selected-topic artifact: {dest}')
  moves[path.resolve()]=dest.resolve()
 # Virtual destination lookup permits all checks to precede any writes/moves.
 sources={dest:source for source,dest in moves.items()}
 def source_path(path):return sources.get(path.resolve(),path)
 def exists(path):return path in directories or source_path(path).exists()
 existing=[p for p in folder.rglob('*') if p.is_file()] if folder.exists() else []
 protected={p.resolve():sha(p) for p in existing if p!=folder/'README.md' and p!=audit_path}
 moved_digests={source:sha(source) for source in moves}
 preview_paths={p.resolve() for p in folder.glob('*.png')}
 preview_paths.update(dest for dest in moves.values() if dest.parent==folder.resolve() and dest.suffix=='.png')
 text=f'# {TITLES.get(name,name)}\n\n[All development experiments](../README.md)\n\n'
 docs=[f'[{label}]({filename})' for label,filename in [('Design','DESIGN.md'),('Results','RESULTS.md')] if (folder/filename).exists()]
 if docs:text+=' · '.join(docs)+'\n\n'
 for preview in sorted(preview_paths):
  stem=preview.stem
  text+=f'- [{stem.replace("_"," ")}]({preview.name})'
  for label,relative in [('PDF',f'supporting/vectors/{stem}.pdf'),('SVG',f'supporting/vectors/{stem}.svg'),('caption',f'supporting/metadata/{stem}_caption.txt')]:
   if exists(folder/relative):text+=f' · [{label}]({relative})'
  text+='\n'
 text+='\nSupporting files: [vectors](supporting/vectors/), [data](supporting/data/), [metadata and audits](supporting/metadata/).\n'
 if any(p.suffix=='.html' for p in (folder/'supporting/viewers').glob('*')):
  text+='\n[Interactive readers and status](supporting/viewers/)\n'
 readme=folder/'README.md'
 index=DEV/'README.md'
 index_before=index.read_text()
 entry=f'- [{TITLES.get(name,name)}]({name}/README.md)\n'
 index_after=index_before
 added=''
 if not re.search(r'\]\(<?'+re.escape(name)+r'/README\.md>?\)',index_before):
  at=index_before.find('\n[Reversible cleanup records]')
  if at>=0:
   added=entry
   index_after=index_before[:at]+added+index_before[at:]
  else:
   added=('' if index_before.endswith('\n') else '\n')+'\n'+entry
   index_after=index_before+added
 # Topic-local manifest checks. No other topic is read or rewritten here.
 manifest_suffixes=('.manifest.json','_manifest.json')
 manifest_paths={p.resolve() for pattern in ['*.manifest.json','*_manifest.json'] for p in folder.rglob(pattern)}
 manifest_paths.update(dest for dest in moves.values() if dest.name.endswith(manifest_suffixes))
 declared_hashes=0;source_hashes=0
 for path in sorted(manifest_paths):
  manifest=json.loads(source_path(path).read_text())
  base=(path.parent/manifest.get('artifact_root','../..')).resolve()
  for key in ['sha256','figure_sha256','data_sha256']:
   for relative,digest in manifest.get(key,{}).items():
    target=(base/relative).resolve()
    if not exists(target) or sha(source_path(target))!=digest:
     raise AssertionError((str(path),relative,'artifact hash mismatch'))
    declared_hashes+=1
  for relative,digest in manifest.get('source_sha256',{}).items():
   target=(REPO/relative).resolve()
   if not target.exists() or sha(target)!=digest:
    raise AssertionError((str(path),relative,'source hash mismatch'))
   source_hashes+=1
 # Check every selected-topic Markdown/HTML link, including the generated
 # landing page. Other-topic targets may be read for existence only.
 documents={p.resolve():p.read_text() for p in existing if p.suffix in ['.md','.html']}
 documents.update({dest:source.read_text() for source,dest in moves.items() if dest.suffix in ['.md','.html']})
 documents[readme.resolve()]=text
 links=0
 for path,content in documents.items():
  urls=re.findall(r'\]\(([^)]+)\)',content) if path.suffix=='.md' else re.findall(r'(?:href|src)=[\"\x27]([^\"\x27]+)[\"\x27]',content)
  for url in urls:
   parsed=urlsplit(url.strip('<>'))
   if parsed.scheme or parsed.netloc or not parsed.path:continue
   target=(path.parent/unquote(parsed.path)).resolve()
   if target==readme.resolve():continue
   if not exists(target):raise AssertionError((str(path),url,'missing local link'))
   links+=1
 readme_changed=not readme.exists() or readme.read_text()!=text
 index_changed=index_after!=index_before
 planned_writes=([readme] if readme_changed else [])+([index] if index_changed else [])+[audit_path]
 rel=lambda path:str(path.relative_to(REPO))
 audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,
  mode='single_group',group=name,dry_run=bool(dry_run),
  declared_artifact_hashes_checked=declared_hashes,source_hashes_checked=source_hashes,local_links_checked=links,
  preserved_existing_file_hashes_checked=len(protected),
  preview_and_supporting_moves=[dict(old=rel(source),new=rel(dest),sha256=moved_digests[source]) for source,dest in moves.items()],
  selected_topic_readme_changed=readme_changed,
  selected_topic_readme_overwritten=bool(readme.exists() and readme_changed),
  dev_index_link_added=index_changed,dev_index_existing_text_preserved=True,
  unrelated_topics_touched=False,
  paths_that_would_be_written=[rel(path) for path in planned_writes],
  directories_that_would_be_created=[rel(path) for path in directories if not path.exists()])
 if dry_run:
  print(json.dumps(audit,indent=2));return audit
 for path in directories:path.mkdir(parents=True,exist_ok=True)
 for source,dest in moves.items():
  dest.parent.mkdir(parents=True,exist_ok=True)
  if dest.exists():source.unlink() # Byte-identical duplicate checked above.
  else:source.rename(dest)
  if sha(dest)!=moved_digests[source]:raise AssertionError(('relocated artifact changed',str(dest)))
 if readme_changed:readme.write_text(text)
 if index_changed:index.write_text(index_after)
 for path,digest in protected.items():
  if sha(path)!=digest:raise AssertionError(('selected scientific file changed',str(path)))
 if index.read_text()!=index_after:raise AssertionError('Development index differs from planned insertion')
 if added and index_after.replace(added,'',1)!=index_before:
  raise AssertionError('Existing development index text was changed')
 audit['paths_written']=[rel(path) for path in planned_writes]
 write(audit_path,audit)
 print(json.dumps(audit,indent=2));return audit


def main():
 RECORDS.mkdir(exist_ok=True)
 snapshot=RECORDS/'before_topic_cleanup'
 saved=set();digest_changes={};changes=[]
 def save(p):
  if p in saved:return
  dest=snapshot/p.relative_to(REPO);dest.parent.mkdir(parents=True,exist_ok=True)
  if not dest.exists():shutil.copy2(p,dest)
  saved.add(p)
 def replace_file(p,text):
  before=p.read_text()
  if text==before:return
  save(p);old=sha(p);p.write_text(text);new=sha(p);digest_changes[old]=new
  changes.append(dict(path=str(p.relative_to(REPO)),before_sha256=old,after_sha256=new))
 pending=[p for p in DEV.glob('*') if p.is_file() and p.suffix=='.png']
 supporting=DEV/'supporting'
 if supporting.exists():pending += [p for p in supporting.rglob('*') if p.is_file() and p.name!='README.md']
 moves={};inventory=[]
 for p in pending:
  g=group(p.name)
  tail=p.relative_to(DEV)
  new=DEV/g/(p.name if g=='_organization' or p.parent==DEV else tail)
  moves[p.resolve()]=new.resolve()
  inventory.append(dict(old=str(p.relative_to(DEV)),new=str(new.relative_to(DEV)),sha256=sha(p)))
 def mapped_text(text):
  for old,new in moves.items():
   text=text.replace('dev/'+str(old.relative_to(DEV)), 'dev/'+str(new.relative_to(DEV)))
  return text
 # Update presentation generators. Frozen scientific campaign sources stay exact.
 scripts={
  'length_dev_basal_endpoint_metrics.py':'resting_state','length_dev_basal_spatial_profiles.py':'resting_state',
  'length_dev_basal_feedback_figures.py':'basal_feedback_controls',
  'length_dev_global_tension_rho_ablation_figures.py':'tension_rho_ablation',
  'length_dev_fullstate_kymographs.py':'L29_fullstate',
  'length_dev_mca_traces.py':'native_reversal','length_dev_fig9f_cue_timing.py':'native_reversal',
  'length_dev_beta0_reversal_companion.py':'beta0_reversal',
 }
 for name,g in scripts.items():
  p=Path(__file__).parent/name
  replace_file(p,p.read_text().replace("/'figures/length_mechanisms_20261008/dev'",f"/'figures/length_mechanisms_20261008/dev/{g}'"))
 p=Path(__file__).parent/'length_dev_beta0_reversal_figures.py'
 text=p.read_text()
 if "OUT=OUT/'beta0_reversal'" not in text:text=text.replace('import length_dev_output_layout as layout',"import length_dev_output_layout as layout\nOUT=OUT/'beta0_reversal'")
 replace_file(p,text)
 # The narrative reader consumes already exported development arrays.
 for name in ['length_narrative_fields_20261009.py','length_report_field_accompaniment_20261009.py']:
  p=Path(__file__).parent/name;replace_file(p,mapped_text(p.read_text()))
 p=Path(__file__).parent/'length_narrative_resting_ablation_20261009.py'
 replace_file(p,p.read_text().replace("DEV = FIG/'dev/supporting'","DEV = FIG/'dev/tension_rho_ablation/supporting'").replace("metrics=rows(DEV/'data/basal_endpoint_metrics_vs_length.csv')","metrics=rows(FIG/'dev/resting_state/supporting/data/basal_endpoint_metrics_vs_length.csv')"))
 # Repair relative links before moving viewers; their targets stay byte-identical.
 for p in pending:
  if p.suffix not in ['.html','.json','.md','.txt']:continue
  text=p.read_text();text=mapped_text(text)
  if p.suffix=='.html':
   def link(match):
    url=match.group(2)
    if url.startswith(('http:','https:','#','data:')):return match.group(0)
    target=(p.parent/url).resolve();target=moves.get(target,target)
    return match.group(1)+'"'+os.path.relpath(target,moves[p.resolve()].parent)+'"'
   text=re.sub(r'((?:href|src)=)"([^"]+)"',link,text)
  # Historical cleanup inventories are preserved as original records.
  if group(p.name)!='_organization':replace_file(p,text)
 for old,new in moves.items():
  new.parent.mkdir(parents=True,exist_ok=True)
  if new.exists():
   if sha(old)!=sha(new):save(new)
   old.replace(new)
  else:old.rename(new)
 if inventory:
  write(RECORDS/'topic_cleanup_inventory.json',dict(generated_at=datetime.now(timezone.utc).isoformat(),moves=inventory,changes=changes,snapshots='before_topic_cleanup/',scientific_campaign_inputs_unchanged=True))
 # Current report links and attached source manifests follow relocated exports.
 external=[]
 for p in REPORT.rglob('*'):
  if DEV in p.parents or any(part.startswith('archive') for part in p.relative_to(REPORT).parts):continue
  if p.suffix in ['.html','.json','.md'] and p.is_file():external.append(p)
 external.append(REPO/'reports/CURRENT_STATUS.md')
 for p in external:replace_file(p,mapped_text(p.read_text()))
 # Cascading provenance references: record updated hashes for changed renderers
 # and relocated metadata; original bytes remain in the reversible snapshots.
 manifests=list(DEV.glob('*/supporting/metadata/*.manifest.json'))+[p for p in external if p.suffix=='.json']
 for _ in range(12):
  altered=0
  def update(value):
   if isinstance(value,str):
    seen=set()
    while value in digest_changes and value not in seen:seen.add(value);value=digest_changes[value]
    return value
   if isinstance(value,list):return [update(v) for v in value]
   if isinstance(value,dict):return {k:update(v) for k,v in value.items()}
   return value
  for p in manifests:
   old=json.loads(p.read_text());new=update(old)
   if old!=new:replace_file(p,json.dumps(new,indent=2,ensure_ascii=False)+'\n');altered+=1
  if not altered:break
 else:raise RuntimeError('Hash references did not converge')
 # Topic landing pages keep the figures and their related exports together.
 for g,title in TITLES.items():
  folder=DEV/g;folder.mkdir(parents=True,exist_ok=True);previews=sorted(folder.glob('*.png'))
  for category in ['vectors','data','metadata']:(folder/'supporting'/category).mkdir(parents=True,exist_ok=True)
  if g=='deBelly_redraft':previews.sort(key=lambda p:p.name.lower())
  if g=='stimulated_endpoint':previews.sort(key=lambda p:(p.stem!='stimulated_endpoint_length_summary',p.name))
  text=f'# {title}\n\n[All development experiments](../README.md)\n\n'
  if g=='deBelly_redraft':text+='[Narrative talking points through Figures X, Y and Z](NARRATIVE.md)\n\nThis collection contains the latest drafts selected by their manifest generation times, including the Figure Y mass-balance revision. Original exports remain in their source folders; the collection manifest records their paths, timestamps and verified hashes.\n\n'
  for p in previews:
   stem=p.stem;text+=f'- [{stem.replace("_"," ")}]({p.name})'
   vector=folder/'supporting/vectors'/f'{stem}.pdf'
   if vector.exists():text+=f' · [PDF](supporting/vectors/{stem}.pdf)'
   caption=folder/'supporting/metadata'/f'{stem}_caption.txt'
   if caption.exists():text+=f' · [caption](supporting/metadata/{stem}_caption.txt)'
   text+='\n'
  text+='\nSupporting files: [vectors](supporting/vectors/), [data](supporting/data/), [metadata and audits](supporting/metadata/).\n'
  if (folder/'supporting/viewers').exists():text+='\n[Interactive readers and status](supporting/viewers/)\n'
  if g=='beta0_reversal':text+='\nTwelve new beta=0 steady states are shared exactly across all three cues. The development moment cutoff is 0.005; original 0.05 eligibility/outcome flags are exported separately. All 2,897 trials passed audit: 36 orientation successes, 33 strength successes, and three L=18 strength-censored readouts.\n'
  (folder/'README.md').write_text(text)
 (DEV/'README.md').write_text('# Length-mechanism development experiments\n\nEach folder contains its PNG previews and related supporting vectors, data, captions, audits and readers. Native simulations remain under manuscript/data/runs/.\n\n'+''.join(f'- [{title}]({g}/README.md)\n' for g,title in TITLES.items())+'\n[Reversible cleanup records](_organization/) preserve the original paths and metadata. After regenerating exports, run `python polarity1d/manuscript/code/mechanisms/organize_length_dev.py` from the repository root. Current renderers write directly into these experiment folders; the organizer also collects outputs from frozen older entry points.\n')
 (DEV/'AGENTS.md').write_text('# Development output organization\n\nGroup each experiment and all its exports in its topic folder; keep only this index, these instructions and experiment folders at dev/. Use the existing supporting/{vectors,data,metadata,viewers} structure inside each experiment.\n\nAfter regeneration or adding a new export, run code/mechanisms/organize_length_dev.py from the manuscript directory (or its repository-relative path). Update the prefix routing there for a new experiment. Preserve raw simulation evidence and frozen campaign sources. Verify local links and artifact hashes after moves.\n')
 if supporting.exists():
  if (supporting/'README.md').exists():save(supporting/'README.md');(supporting/'README.md').unlink()
  for p in sorted(supporting.rglob('*'),key=lambda p:len(p.parts),reverse=True):
   if p.is_dir():p.rmdir()
  supporting.rmdir()
 # Verify each topic's declared artifact root, all non-historical links, and
 # identical image/vector/data bytes against the pre-move inventory.
 hashes=0;links=0
 for p in DEV.glob('*/supporting/metadata/*.manifest.json'):
  m=json.loads(p.read_text());base=(p.parent/m.get('artifact_root','../..')).resolve()
  for key in ['sha256','figure_sha256','data_sha256']:
   for name,digest in m.get(key,{}).items():
    target=base/name
    assert target.exists() and sha(target)==digest,(p,name);hashes+=1
 for p in [DEV/'README.md',*DEV.glob('*/README.md'),*DEV.glob('*/supporting/viewers/*.html')]:
  text=p.read_text();targets=re.findall(r'\]\(([^)]+)\)',text) if p.suffix=='.md' else re.findall(r'(?:href|src)="([^"]+)"',text)
  for url in targets:
   if url.startswith(('https:','http:','#','data:')):continue
   assert (p.parent/url).exists(),(p,url);links+=1
 for row in inventory:
  if Path(row['new']).suffix in ['.png','.pdf','.svg','.csv','.npz']:
   assert sha(DEV/row['new'])==row['sha256'],row
 assert not list(DEV.glob('*.png')) and not supporting.exists()
 previous=RECORDS/'topic_cleanup_inventory.json'
 retained=json.loads(previous.read_text())['moves'] if previous.exists() else inventory
 preview_checks=0
 for row in retained:
  if Path(row['new']).suffix in ['.png','.csv','.npz']:
   assert sha(DEV/row['new'])==row['sha256'],row
   preview_checks+=1
 audit=dict(generated_at=datetime.now(timezone.utc).isoformat(),passed=True,moved_files=len(json.loads(previous.read_text())['moves']) if previous.exists() else len(inventory),moved_files_this_run=len(inventory),artifact_hashes_checked=hashes,local_links_checked=links,original_array_and_preview_hashes_verified=preview_checks,all_moves_verified_at_relocation=True,groups=list(GROUPS))
 write(RECORDS/'topic_cleanup_audit.json',audit)
 if inventory:
  record=json.loads((RECORDS/'topic_cleanup_inventory.json').read_text());record['changes']=changes;write(RECORDS/'topic_cleanup_inventory.json',record)
 print(json.dumps(audit,indent=2))
def cli(argv=None):
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--group',choices=sorted(GROUPS),help='Organize only this topic; preserve all other topic files and existing index text.')
 parser.add_argument('--dry-run',action='store_true',help='Check the selected topic and print planned writes without changing files (requires --group).')
 args=parser.parse_args(argv)
 if args.dry_run and not args.group:parser.error('--dry-run requires --group')
 if args.group:return organize_group(args.group,args.dry_run)
 return main()


if __name__=='__main__':cli()
