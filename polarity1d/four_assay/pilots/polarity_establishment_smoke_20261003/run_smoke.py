"""One finite patch/release through each native solver at both campaign widths."""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import polarity_establishment as e

out = Path(__file__).resolve().parent
manifest = e.make_manifest()
cases = [c for c in manifest['cases'] if c['length'] == e.REFERENCES[c['model']]
         and c['amplitude_factor'] == 1.]
case = cases[int(sys.argv[1])]
case['case_index'] = int(sys.argv[1])
duration = .1 if case['model'] == 'spring' else 10.
result = e.trial(out, case, duration, manifest)
print(json.dumps(dict(model=case['model'], width=case['width'],
                      verified=result['verified'], final_polarity=result['cue_aligned_strength'],
                      polarity_arrival=result['polarity_arrival'], release_elapsed=result['convergence']['elapsed'])))
