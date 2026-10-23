"""Read-only inventory of this dated deBelly attempt; never deletes files."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os,re,subprocess

ROOT=next(p for p in Path(__file__).resolve().parents if (p/'polarity1d/manuscript/AGENTS.md').exists())
MAN=ROOT/'polarity1d/manuscript';DEV=MAN/'figures/length_mechanisms_20261008/dev';OUT=Path(__file__).parent
OLD=DEV/'debelly_ring_development_20261019';STORAGE=OLD/'supporting/metadata/storage_footprint.json';CODE=MAN/'code/mechanisms'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def label(p):
 try:return str(Path(p).relative_to(ROOT))
 except ValueError:return str(p)
def now():return datetime.now(timezone.utc).isoformat()
def inventory(p):
 p=Path(p);seen=set();allocated=logical=count=0;links=[]
 def add(q):
  nonlocal allocated,logical,count
  st=q.lstat();key=(st.st_dev,st.st_ino)
  if key in seen:return
  seen.add(key);allocated+=st.st_blocks*512
  if q.is_symlink():links.append(dict(path=label(q),target=os.readlink(q)));return
  if q.is_file():logical+=st.st_size;count+=1
 if p.is_dir() and not p.is_symlink():
  add(p)
  for base,dirs,files in os.walk(p,followlinks=False):
   for name in [*dirs,*files]:add(Path(base)/name)
 else:add(p)
 return dict(allocated_bytes=allocated,logical_bytes=logical,files=count,symlinks=links)
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 snapshot=json.loads(STORAGE.read_text());tracked=set(subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0'))
 groups=[ROOT/g['path'] for g in snapshot['groups']]
 dated=sorted(p for p in (MAN/'data/runs').glob('*debelly*20261019*') if p.is_dir())
 group_set=set(groups);extra=[p for p in dated if p not in group_set];groups.extend(extra)
 code=sorted(p for p in CODE.iterdir() if p.is_file() and 'debelly' in p.name.lower() and '20261019' in p.name)
 stems={p.stem for p in code};caches=[]
 for cache in CODE.glob('__pycache__/*'):
  if any(cache.name.startswith(stem+'.') for stem in stems):caches.append(cache)
 repotmp=[];historical_tmp=[];tmp=ROOT/'tmp'
 if tmp.exists():
  named=sorted(p for p in tmp.rglob('*') if 'debelly' in p.name.lower())
  for p in named:
   # Name matching alone cannot establish attempt ownership. The extant
   # x50_mechanisms_20261003 material predates this attempt and is protected.
   if not any(parent in historical_tmp for parent in p.parents):historical_tmp.append(p)
 systemtmp=[];ambiguous=[]
 for p in Path('/tmp').iterdir():
  if 'debelly' not in p.name.lower():continue
  st=p.lstat();mtime=datetime.fromtimestamp(st.st_mtime,timezone.utc)
  if st.st_uid==os.getuid() and mtime.date().isoformat()=='2026-10-19':systemtmp.append(p)
  else:ambiguous.append(dict(path=str(p),owner_uid=st.st_uid,mtime_utc=mtime.isoformat(),reason='Not automatically classified as this attempt.'))
 candidates=[]
 status=[MAN/'DEBELLY_RING_STATUS_20261019.md']
 for category,paths in [('dated_run_or_old_topic',groups),('dated_sidecar_code',code),('matching_code_cache',caches),('standalone_status',status),('repository_temp',repotmp),('system_temp',systemtmp)]:
  for p in paths:
   if not p.exists() and not p.is_symlink():continue
   row=dict(path=label(p),category=category,exists=True,**inventory(p))
   st=p.lstat();row.update(owner_uid=st.st_uid,mtime_utc=datetime.fromtimestamp(st.st_mtime,timezone.utc).isoformat())
   if category=='system_temp':row['ownership_evidence']='Current user, UTC mtime 2026-10-19, deBelly task-specific name; retained browser helper is copied inside methods/check_viewer.mjs.'
   pref=label(p).rstrip('/')+'/';desc=[x for x in tracked if x==label(p) or x.startswith(pref)]
   row['tracking']=dict(tracked_files=len(desc),tracked_paths=sorted(desc))
   if p.is_file():row['sha256']=sha(p)
   elif p in groups:
    row['provenance_hashes']={label(q):sha(q) for q in [p/'manifest.json',p/'implementation_manifest.json',p/'validation_result.json',p/'inputs.json',p/'README.md'] if q.is_file()}
   row['deletion_condition']='Only after the retained review is verified self-contained and required incoming references are updated.'
   candidates.append(row)
 # Search incoming textual references outside the candidate trees, retained
 # folder and .git. Internal references within a deleted group are not blockers.
 names=sorted({p.name for p in [*groups,*code,*status]}|{p.stem for p in code if p.suffix=='.py'},key=len,reverse=True)
 pattern='|'.join(re.escape(name) for name in names)
 command=['rg','--files-with-matches','--hidden','--no-ignore','--max-filesize','5M','--glob','*.{md,txt,json,jsonl,html,htm,py,js,mjs,sh,slurm,yml,yaml,toml,cfg,ini,rst,tex,R,jl}','--glob','.gitignore','--glob','!.git/**','--glob','!**/node_modules/**','--glob','!**/__pycache__/**','--glob','!'+label(OUT)+'/**']
 for p in [*groups,*code,*status,*repotmp]:command.extend(['--glob','!'+label(p)+('/**' if p.is_dir() else '')])
 command.extend(['-e',pattern,'.']);search=subprocess.run(command,cwd=ROOT,text=True,capture_output=True)
 if search.returncode not in (0,1):raise RuntimeError(search.stderr)
 expression=re.compile(pattern);references=[]
 for rawpath in sorted(set(search.stdout.splitlines())):
  if rawpath.startswith('Binary file '):continue
  p=ROOT/rawpath.removeprefix('./')
  try:lines=p.read_text(errors='replace').splitlines()
  except OSError:continue
  hits=[]
  for number,line in enumerate(lines,1):
   found=sorted(set(expression.findall(line)))
   if found:hits.append(dict(line=number,targets=found,snippet=line[:240]))
  if hits:
   references.append(dict(path=label(p),sha256=sha(p),tracked=label(p) in tracked,
     runtime_consumer=p.suffix in ('.py','.js','.mjs','.sh','.slurm') and p.name!='organize_length_dev.py',
     hits=hits,action='Review and update narrow references; preserve unrelated content and other-goal changes.'))
 # Enumerate protected established deBelly data/topics and independent current
 # work. These lists do not grant deletion authorization for any omitted path.
 protected=[]
 old_data=[p for p in (MAN/'data/runs').glob('*debelly*') if p.is_dir() and p not in group_set and '20261019' not in p.name]
 old_topics=[p for p in DEV.glob('*debelly*') if p.is_dir() and p not in (OLD,OUT)]
 independent=[p for p in [DEV/'spring_ring_validation_20261019',DEV/'revision_program_20261019',MAN/'data/runs/revision_ring_spring_20261019',MAN/'data/runs/revision_spring_noise_main_20261019',MAN/'data/runs/revision_initial_arc_sizes_20261019',ROOT/'literature_models/debelly_2026_mechanochemical',ROOT/'polarity1d/four_assay',ROOT/'current manuscript figs'] if p.exists()]
 for p in [*old_data,*old_topics,*independent]:protected.append(dict(path=label(p),reason='Established/pre-existing deBelly evidence or independently successful/active other-goal work; outside this attempt cleanup scope.'))
 protected.extend(dict(path=label(p),reason='Repository temporary material predating this attempt or lacking demonstrated ownership; explicitly preserve.') for p in historical_tmp)
 codeprotected=[p for p in CODE.iterdir() if p.is_file() and p not in code and ('debelly' in p.name.lower() or 'spring' in p.name.lower() or 'revision_' in p.name)]
 protected.extend(dict(path=label(p),reason='Not an attempt-owned dated deBelly sidecar source; preserve.') for p in codeprotected)
 for p in [MAN/'code/core/models/debelly.py',MAN/'code/core/src/polarity1d/debelly.py']:
  if p.exists():protected.append(dict(path=label(p),reason='Active established model source, unchanged by this sidecar attempt.'))
 for row in protected:
  p=ROOT/row['path'];row.update(inventory(p))
  if p.is_file():row['sha256']=sha(p)
  else:
   sample=sorted(q for q in p.rglob('*') if q.is_file() and not q.is_symlink() and q.stat().st_size<2*1024*1024)
   favored=[q for q in sample if q.name in ('manifest.json','README.md','inputs.json','summary.json','qualified_state.npz')]
   selected=list(dict.fromkeys([*favored[:3],*sample[:2],*sample[-2:]]))
   row['protected_hash_samples']={label(q):sha(q) for q in selected}
 processes=[]
 for line in subprocess.check_output(['ps','-eo','uid,pid,ppid,etime,args','--no-headers'],text=True).splitlines():
  parts=line.split(None,4)
  if len(parts)!=5 or parts[0]!=str(os.getuid()) or int(parts[1])==os.getpid():continue
  if 'debelly' not in parts[4].lower():continue
  if 'cleanup_inventory.py' in parts[4] or 'ps -eo ' in parts[4]:continue
  processes.append(dict(pid=int(parts[1]),ppid=int(parts[2]),elapsed=parts[3],command=parts[4],retained_package_process=label(OUT) in parts[4]))
 queue=subprocess.run(['squeue','-u',os.environ.get('USER','4473331'),'-h','-o','%i|%j|%T|%R'],text=True,capture_output=True)
 jobs=queue.stdout.splitlines();debelly_jobs=[j for j in jobs if 'debelly' in j.lower()]
 baseline=subprocess.check_output(['git','diff','--name-status'],cwd=ROOT,text=True)
 result=dict(created_utc=now(),mode='READ_ONLY_INVENTORY_NO_DELETIONS',authorization='User replaced development with consolidation and cleanup of this attempt. Deletion is conditional on verifying the retained package; no additional permission is requested.',
  retained_folder=label(OUT),storage_snapshot=dict(path=label(STORAGE),sha256=sha(STORAGE),measured_utc=snapshot['measured_utc'],groups=len(snapshot['groups']),dated_run_roots=sum('/data/runs/' in g['path'] for g in snapshot['groups']),additional_dated_roots=[label(p) for p in extra]),
  candidate_delete=candidates,candidate_allocated_bytes=sum(r['allocated_bytes'] for r in candidates),
  protected=protected,ambiguous_system_temp=ambiguous,refs_to_update=references,
  external_runtime_consumers=[r['path'] for r in references if r['runtime_consumer']],
  reference_search_scope=dict(text_extensions='md txt json jsonl html htm py js mjs sh slurm yml yaml toml cfg ini rst tex R jl and .gitignore',maximum_file_bytes=5*1024*1024,excluded='Candidate trees/files, retained folder, .git, node_modules, __pycache__',limitation='Historical identifiers in binary/large archives are not active text dependencies. Local retained runtime closure is checked separately.'),
  live_state=dict(observed_utc=now(),squeue_exit=queue.returncode,squeue_all_jobs=jobs,debelly_named_jobs=debelly_jobs,local_debelly_processes=processes,no_named_scientific_debelly_job=queue.returncode==0 and not debelly_jobs,note='Other-goal f4/conmon/RStudio jobs are protected. Command-line checks cannot exclude hidden generic processes; candidate source owners also confirm all scientific producers terminal.'),
  working_tree_baseline=baseline.splitlines(),
  requirements_before_delete=['Verify local review runtime/render dependency closure and retained per-file provenance hashes.',
    'Resolve external runtime consumers or retain their required dependencies; update only attempt references in docs/ignore/organizer.',
    'Record pre/post allocated storage and update root DATA_MANIFEST.md with the deletion and retained package measurement.',
    'Preserve unrelated tracked modifications/deletions, corrected linear evidence, manuscript figures, independent Spring and other-goal work.',
    'Stop only this attempt’s remaining watchers/processes before removal; queued or running other-goal jobs are protected.'],
  source_sha256={label(Path(__file__)):sha(__file__)})
 target=OUT/'cleanup_plan.json';target.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
 print(json.dumps(dict(path=label(target),sha256=sha(target),groups=len(groups),dated_run_roots=len(dated),code=len(code),caches=len(caches),repo_tmp=len(repotmp),system_tmp=len(systemtmp),candidate_allocated_bytes=result['candidate_allocated_bytes'],incoming_reference_files=len(references),external_runtime_consumers=result['external_runtime_consumers']),indent=2))
if __name__=='__main__':main()
