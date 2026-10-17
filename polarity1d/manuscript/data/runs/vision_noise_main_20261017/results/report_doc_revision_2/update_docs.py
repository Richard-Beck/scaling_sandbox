"""Close documentation after the accepted actual-result critique; no data edits."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,shutil

FOLDER=Path(__file__).resolve().parent;RUN=FOLDER.parents[1];ROOT=RUN.parents[2];REPO=ROOT.parents[1]
OUT=ROOT/"figures/length_mechanisms_20261008/dev/vision_noise_20261017"

def sha(path):
 d=hashlib.sha256()
 with Path(path).open('rb') as f:
  for chunk in iter(lambda:f.read(2**20),b''):d.update(chunk)
 return d.hexdigest()

report=OUT/'RESULTS.md';manifest_path=OUT/'supporting/metadata/vision_noise_main_report_manifest.json'
review=OUT/'supporting/metadata/FRESH_RESULT_CRITIQUE.md'
assert sha(review)=='7d4ce5f5916ab71534dfb62b0bfd9482be22193de7acdba951d979944eaf2071'
for source,name in ((report,'previous_RESULTS.md'),(manifest_path,'previous_report_manifest.json')):
 target=FOLDER/name;assert not target.exists();shutil.copy2(source,target)
old_report_sha=sha(report);old_manifest_sha=sha(manifest_path)
text=report.read_text()
rationale=('The paired memory and magnitude curves test retained directional asymmetry; '
           'native histories and physical displacement with coverage distinguish '
           'magnitude loss, multisite cancellation and ring geometry from slower site motion.\n\n')
needle='## Question and declared comparison\n';assert text.count(needle)==1
text=text.replace(needle,rationale+needle)
needle='![Direction eligibility and physical geometry](vision_noise_coverage_geometry.png)\n\n';assert text.count(needle)==1
text=text.replace(needle,needle+"The shaded Wave-pinning and Holmes RMS bands are the monotone square roots\nof the stored pointwise 95% whole-seed bootstrap bounds on mean squared\nendpoint displacement.\n\n")
old=('A fresh-context actual-result critique will assess the curves, local histories,\n'
     'uncertainty and scientific claims before this report is folded into the final\n'
     'three-figure checkpoint. Output generation alone does not establish an accepted\n'
     'manuscript figure or a general enlargement effect.')
assert text.count(old)==1
new=('The [fresh actual-result critique](supporting/metadata/FRESH_RESULT_CRITIQUE.md)\n'
     'finds no essential scientific, statistical or display correction remaining for\n'
     'this noisy-ring component of the review checkpoint. It independently checks\n'
     'the complete statistics, seed uncertainty, all seven plots, representative\n'
     'selection and realized-input audits. The finite-window, native-unit, forced-state\n'
     'and physical-coverage limits stated above remain part of the result.')
text=text.replace(old,new);report.write_text(text)
manifest=json.loads(manifest_path.read_text());manifest['sha256']['RESULTS.md']=sha(report)
sources=[Path(__file__),FOLDER/'previous_RESULTS.md',FOLDER/'previous_report_manifest.json',review,
         OUT/'supporting/metadata/FRESH_RESULT_CHECKS.json',OUT/'supporting/metadata/FRESH_INPUT_AUDIT_CHECKS.json',
         OUT/'supporting/metadata/FRESH_RENDER_CHECKS.json',RUN/'results/statistics_centroid_amendment_20261017/report.py']
for path in sources:manifest['source_sha256'][str(path.relative_to(REPO))]=sha(path)
manifest['documentation_revision']=dict(generated_at=datetime.now(timezone.utc).isoformat(),
  source=str(Path(__file__).relative_to(REPO)),previous_report_sha256=old_report_sha,
  previous_manifest_sha256=old_manifest_sha,review_sha256=sha(review),statistical_or_figure_change=False)
manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
completion=dict(completed_at=datetime.now(timezone.utc).isoformat(),report_sha256=sha(report),manifest_sha256=sha(manifest_path),
                review_sha256=sha(review),unchanged_statistics_audit_sha256=sha(OUT/'supporting/metadata/vision_noise_analysis_audit.json'),
                modified_scientific_csvs=0,modified_figures=0,modified_native_archives=0)
(FOLDER/'completion.json').write_text(json.dumps(completion,indent=2)+'\n')
print(json.dumps(completion,indent=2))
