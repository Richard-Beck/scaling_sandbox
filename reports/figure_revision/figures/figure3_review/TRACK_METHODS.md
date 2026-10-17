# Figure 3 track annotations

`comparison.png` stacks three visual alternatives for review. The user selected
Option 3, saved as `figure3_annotated_tracks.png` without the review-option
label. The comparison and each option are retained as PNG files. All options use the same selected pair of
Xuan2018 model-generated tracks from the primary Gaussian OU hierarchy fit;
they are not observed individual-cell tracks. The display windows are 12 h at
10 min cadence.

- Option 1 puts the track measurements and formulas beside each path.
- Option 2 overlays both paths at one origin and uses a numbered measurement key.
- Option 3 separates the two example paths from the fitted population MSD curve.
  Dashed chords show equal 2 h displacement windows. The fitted curve is the
  OU model expectation across cells; a track MSD instead squares and averages
  displacements at a fixed lag across tracks and start times.

Each displayed path is selected from the 32 existing simulated paths in its
group as the candidate closest to group medians in sampled speed, directness,
and endpoint displacement, using robust standardized distance. Path length is
excluded from the selection score because it equals speed times the fixed
duration. The same paths are used in every option.

The displayed path length is the sum of its sampled step lengths. Its mean
sampled speed is path length divided by 12 h. Its directness is its endpoint
displacement divided by its path length. The MSD values and curve use the exact
integrated OU displacement function used by the fitting pipeline; a runtime
check verifies the displayed 11.111 h values against the selected fit rows.
These fits use approximate literature summary extractions, and the Xuan2018
WGD+ label denotes that study's giant polyploid contrast.

Rebuild from the repository root with the existing environment:

```sh
OPENBLAS_NUM_THREADS=1 conda run -n scaling_sandbox python reports/figure_revision/figure3_tracks.py
```

`figure3_tracks_provenance.json` records input hashes, fit identity, parameter values, the
selection rationale, and per-track metrics. The two representative coordinates
are also exported as CSV.
