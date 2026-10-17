# Current assay explanation figure

New revision following the visual style of the approved September 29
`assay_explanation_v2` keeper. The keeper remains intact in the historical
archive. This new figure illustrates the active revision-10 suite and is not
itself recorded as user-approved.

Files: `assay_explanation.png`, `.pdf`, `.svg`; editable generator
`make_assay_explanation.py`; plotted fields in `assay_explanation_fields/`;
source/output hashes and numerical checks in `assay_explanation_provenance.json`.
The provenance records both figure generation time and results generation time.

## Caption

Preparation and response measurements for the two active polarity assays.
All examples use wave-pinning at L = 21 μm, from current tasks 24 and 25.
Each row shows applied cue (left), active-species profiles coloured by snapshot
time (centre), and first-moment orientation (right). Grayscale denotes cue
intensity from zero (white) to the model's maximum (black). Profiles share one
minimum and range across all five snapshots within a row; they are not scaled
individually. Snapshot colours agree between profile curves, cue-map vertical
lines and response points. Time is in native model units.

(A) Shared conditioning. Starting directly from the native initial state, the
left linear cue establishes left polarity. Both assays use this same cue and
settle after 1,200 units; conditioning appears once. (B) Release / cue-free
settling, required only before finite countercue. The conditioning cue is
withdrawn, and autonomous left polarity settles after 1,300 units. Gray shading
in A and B marks the final 1,000-unit confirmation windows. These endpoints use
full-state relative range and projected drift at tolerance 0.001; orientation
is displayed for explanation and does not alone determine settling. Eligibility
also requires stable left orientation and sufficient spatial contrast.

(C) Gradient reversal. The prepared left cue cross-fades to the matching right
cue over 20 units, with the right cue subsequently maintained. The response
endpoint is 35 units from cross-fade onset: right orientation and active contrast
reach 0.05 and persist for 20 units (gray shading). The plotted interval shows
the response transient; the actual observation horizon is 20,000 units. The
negative-time segment is the last 100 units of conditioning.

(D) Finite countercue. After cue-free release settling, a right linear cue is
applied temporarily. The shortest passing tested duration is 22 units (solid
response); the neighbouring tested duration of 20 fails (dashed response).
Both are followed by 10,000 units of cue-free recovery. The last 500 units
(gray shading) confirm autonomous right orientation above 0.05, sufficient
contrast, and orientation range below 0.005. The endpoint measures applied cue
duration and excludes preparation, release settling and recovery observation.
The negative-time segment is the last 100 units of cue-free preparation.
Cue and profiles show the passing trial. The time axis is symmetric logarithmic
with a linear central region through absolute time 30.

Within each of the seven models, conditioning and challenge use matching linear
cue shapes and amplitudes; left and right are mirror images. Amplitudes and
native units differ between models. The wave-pinning amplitude is 0.1.
Competition and periodic following are deferred and omitted. These examples
explain protocols and readouts; they do not establish length scaling.

## Reconstruction and verification

The gradient spatial field is read directly from the current accepted diagnostic
trace. Revision 10 reuses the non-Spring revision-9 records; their diagnostic and
checkpoint paths are resolved to the preserved v9 archive after relocation.
No historical runtime code is imported. Conditioning and release fields are
reconstructed in the exact recorded integration chunks with the active runner;
their final full states match checkpoint arrays within absolute 1e-10 / relative
1e-9 tolerance. Checkpoint SHA-256 hashes are verified before use. The gradient
endpoint is rescored and matches 35 exactly. Two finite trials are reconstructed
from the verified current checkpoint using the active pulse/recovery solvers;
both outcomes and final orientations match the current recorded trials (the
orientations match exactly in this run). Saved figure fields are separate from
campaign outputs. The PNG was inspected for labels, cue/profile correspondence,
time keys, response markers and recovery timing.

Rebuild using the existing environment (one core, a few seconds):

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
/home/4473331/.conda/envs/scaling_sandbox/bin/python \
polarity1d/four_assay/results/make_assay_explanation.py
```
