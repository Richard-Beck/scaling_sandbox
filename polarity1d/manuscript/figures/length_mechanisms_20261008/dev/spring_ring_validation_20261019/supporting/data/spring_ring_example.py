"""Small restart example; run from the repository root inside SLURM.

Output is an API example, not a new frozen scientific campaign.
"""
from pathlib import Path
import argparse
import sys
import numpy as np

root = Path.cwd() / "polarity1d/manuscript"
sys.path.insert(0, str(root / "code/mechanisms"))
from spring_ring_model_20261019 import Closure, CuePhase, run_protocol

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, default=Path("/tmp/spring_ring_example.npz"))
args = parser.parse_args()
if args.output.exists():
    raise FileExistsError("Choose a fresh --output path")
model = Closure(resting_circumference=2.0, cells=200)
with np.load(root / "data/runs/spring_ring_preparation_pilot_20261017/workers/qualified_basal.npz") as saved:
    state = saved["state"].copy()
incumbent = model.material_positions[np.argmax(state[:model.cells])]
result = run_protocol(model, [
    CuePhase(duration=0.1, peak=4.0, center_fraction=incumbent + 0.4,
             total_support=0.2, output_interval=0.1),
    CuePhase(duration=13000.0, output_interval=5.0),
], initial_state=state)
args.output.parent.mkdir(parents=True, exist_ok=True)
np.savez_compressed(args.output, time=result["time"], states=result["states"],
                    physical_amount=result["physical_amount"],
                    polarity_magnitude=result["polarity_magnitude"])
print("Saved full states:", args.output)
