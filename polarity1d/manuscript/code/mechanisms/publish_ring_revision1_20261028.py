"""Publish only the authorized ring revision after its existing jobs finish.

Uses a private Git index, preserving unrelated staged and working changes.
No simulation is submitted. The expected branch and code hashes are sealed in
revision1_delivery.json before this publisher is enabled.
"""
from pathlib import Path
from datetime import datetime,timezone
import fcntl,hashlib,json,os,subprocess,tempfile
import audit_ring_revision1_20261028 as A

ROOT=A.R1.ROOT;BASE=A.BASE;META=A.META;BRANCH='codex/ring-r1r2-20261027'
START='<!-- ring-r1r2-20261027 begin -->';STOP='<!-- ring-r1r2-20261027 end -->'
def now():return datetime.now(timezone.utc).isoformat()
def write(p,v):A.R1.write(p,v)
def git(args,env=None,input=None,strip=True):
    value=subprocess.run(['git',*args],cwd=ROOT,env=env,input=input,text=True,capture_output=True,check=True).stdout
    return value.strip() if strip else value
def replace_section(text,block):
    a=text.index(START);b=text.index(STOP,a)+len(STOP);return text[:a]+block+text[b:]
def commit_scoped(paths,message,manifest=False):
    """Commit a private index, then align only owned entries in the main one."""
    parent=git(['rev-parse','HEAD'])
    fd,index=tempfile.mkstemp(prefix='ring-rev1-',dir=ROOT/'.git');os.close(fd);Path(index).unlink()
    env={**os.environ,'GIT_INDEX_FILE':index}
    main_manifest=None
    try:
        git(['read-tree',parent],env)
        git(['add','-f','--',*paths],env)
        if manifest:
            current=(ROOT/'DATA_MANIFEST.md').read_text();a=current.index(START);b=current.index(STOP,a)+len(STOP);block=current[a:b]
            committed=replace_section(git(['show',parent+':DATA_MANIFEST.md'],strip=False),block)
            digest=git(['hash-object','-w','--stdin'],input=committed)
            git(['update-index','--add','--cacheinfo','100644',digest,'DATA_MANIFEST.md'],env)
            # Preserve any unrelated staged manifest edits as staged edits.
            staged=replace_section(git(['show',':DATA_MANIFEST.md'],strip=False),block)
            main_manifest=git(['hash-object','-w','--stdin'],input=staged)
        tree=git(['write-tree'],env)
        if tree==git(['rev-parse',parent+'^{tree}']):return parent
        commit=git(['commit-tree',tree,'-p',parent,'-m',message])
        git(['update-ref','-m','ring revision1 scoped release','refs/heads/'+BRANCH,commit,parent])
        git(['reset','-q',commit,'--',*paths])
        if main_manifest:git(['update-index','--add','--cacheinfo','100644',main_manifest,'DATA_MANIFEST.md'])
        return commit
    finally:
        Path(index).unlink(missing_ok=True);Path(index+'.lock').unlink(missing_ok=True)
def push(commit):
    env={**os.environ,'GIT_SSH_COMMAND':'ssh -F /dev/null'}
    git(['push','origin',commit+':refs/heads/'+BRANCH],env)
    remote=git(['ls-remote','origin','refs/heads/'+BRANCH],env).split()[0]
    if remote!=commit:raise RuntimeError('Remote verification differs from the pushed commit: '+remote)
    return remote

def publish():
    lock=ROOT/'.git/ring-revision1-publish.lock'
    with lock.open('a') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX)
        if git(['branch','--show-current'])!=BRANCH:raise RuntimeError('Requested branch is no longer checked out')
        config=json.loads((META/'revision1_delivery.json').read_text())
        if not config.get('commit_and_push_authorized'):raise RuntimeError('Final-data delivery was not enabled')
        git(['merge-base','--is-ancestor',config['source_commit'],'HEAD'])
        for path,digest in config['code_sha256'].items():
            if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:raise RuntimeError('Code changed after the validated delivery setup: '+path)
        audit=A.audit(require_complete=True)
        for figure in ('R1','R2'):
            p=BASE/f'figure{figure}'/'supporting/metadata/visual_review.json';review=json.loads(p.read_text())
            current={ext:A.sha(BASE/f'figure{figure}'/f'figure{figure}_main.{ext}') for ext in ('png','pdf','svg')}
            review['final_data_validation']=dict(generated_utc=now(),scientific_and_display_audit='../../../supporting/metadata/revision1_audit.json',artifact_sha256=current,manual_visual_review_scope='The revision layout was manually reviewed before campaign completion; final-data artist encodings and algebra are checked automatically. No new manual visual inspection is claimed.')
            write(p,review)
        paths=[str((BASE/name).relative_to(ROOT)) for name in ('figureR1','figureR2','supporting/metadata','README.md')]
        artifacts=commit_scoped(paths,'Rebuild revised R2 after both event campaigns and extensions finish',manifest=True)
        push(artifacts)
        write(META/'integration_commit.json',dict(generated_utc=now(),source_commit=config['source_commit'],figure_data_commit=artifacts,remote_branch='origin/'+BRANCH,remote_verified=True,producers_complete=True))
        complete=json.loads((META/'completion_audit.json').read_text());complete.update(goal_complete=True,checkpoint_only=False,figure_data_commit=artifacts,remote_verified=True);write(META/'completion_audit.json',complete)
        final=f'Revision 1 complete.\nSource commit: {config["source_commit"]}\nFinal figure/data commit: {artifacts}\nRemote branch: origin/{BRANCH}\nValidation: {len(audit["checks"])} scientific/display checks passed; all R1 A–H and R2 A–G panels retained; both event producers and extensions complete; final R2 rebuilt.\nRemaining scientific mismatches and finite-record bounds are reported in each scientific_readme.md (Unexpected). No retuning or new simulation campaign.\n'
        (META/'final.txt').write_text(final)
        write(META/'delivery_status.json',dict(generated_utc=now(),state='PUSHED',source_commit=config['source_commit'],figure_data_commit=artifacts,remote_branch='origin/'+BRANCH,remote_verified=True,receipt_commit='See the following delivery-receipt commit in branch history.'))
        receipts=[str((META/name).relative_to(ROOT)) for name in ('integration_commit.json','completion_audit.json','final.txt','delivery_status.json')]
        receipt_commit=commit_scoped(receipts,'Record verified final ring revision delivery')
        push(receipt_commit)
        print(json.dumps(dict(state='PUSHED',figure_data_commit=artifacts,receipt_commit=receipt_commit,branch=BRANCH)))
if __name__=='__main__':
    try:publish()
    except Exception as error:
        write(META/'delivery_status.json',dict(generated_utc=now(),state='NEEDS_ATTENTION',error=repr(error)))
        raise
