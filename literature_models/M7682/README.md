# M7682 — default simulation and snapshots

Downloaded the unmodified [Morpheus model M7682](https://morpheus.gitlab.io/model/m7682/), “Mixed-mode Oscillations in a Rac/Rho/Paxillin Subsystem,” and ran it with native Morpheus 2.3.9 on 2026-09-22.

Reference: Plazen et al. (2023), [Polarity and mixed-mode oscillations may underlie different patterns of cellular migration](https://doi.org/10.1038/s41598-023-31042-8).

## Outputs

- [Six default snapshots, PNG](results/snapshots/default_snapshots.png) and [PDF](results/snapshots/default_snapshots.pdf): times 0, 1000, 2000, 3000, 4000, 5000, each showing cell ID, R, B, and kb.
- [Enlarged cell shapes](results/snapshots/cell_shapes.png): the same fixed spatial window and scale at every time, cropped from native Cell ID plots.
- [Individual native snapshots and simulation outputs](results/default/).
- [Provenance and checksums](results/snapshots/manifest.json).

## Default conditions retained

`model/model.xml` is the exact downloaded current model file, without edits. LK = 5.77; random seed = 19; time 0–5000; 250 × 250 periodic lattice; node length 0.1; membrane resolution 200; one initial cell of radius 20 lattice sites centered at (125,125); target area 1200; target perimeter 200; biochemical integration step 0.1; CPM step duration 0.5. All reaction, diffusion, and mechanical parameters and initialization expressions remain as supplied.

The repository's default duration is a short demonstration; the model webpage notes that paper runs typically use 80,000 time steps. This run does not establish the long-time oscillatory regime.

Native plots retain the model's fixed color limits: R 0.05–0.3, B 0–10, kb 0.039–0.043. Values outside these limits saturate; notably the initial kb is 0.05. Native kb tick labels are partly clipped by the model's plot layout. Membrane fields are drawn along the cell boundary; the white interior is not a concentration measurement.

## Execution and validation

SLURM job 20510538 completed successfully (exit 0, 34 seconds, approximately 54 MB maximum resident memory), using one CPU, 2 GB requested memory, partition `red`, QoS `small`. Morpheus reported simulation completion; all six PNGs decode successfully; displacement data are finite and extend through time 5000.

Gnuplot 5.4.2 renders the native plots. It reports missing matrix entries for membrane-only fields outside their defined domain. An initial successful simulation using Gnuplot 6.0.3 produced deprecated `set clabel` syntax errors; its outputs are retained in `results/gnuplot6_trial/`. Repeating with Gnuplot 5.4.2 removed those syntax errors and produced identical displacement data.

Reproduce from the repository root:

```bash
sbatch literature_models/M7682/scripts/run_default.slurm
# After the job finishes:
/home/4473331/.conda/envs/scaling_sandbox/bin/python literature_models/M7682/scripts/make_snapshots.py
```

The simulator is the official static Linux binary in `software/morpheus-2.3.9`; no new Conda environment was needed. The snapshot assembly script uses the existing `scaling_sandbox` environment. Rerunning writes to the same output directories.

## Matched area-scaling experiment

The 80,000-time-unit campaign under `results/area_scaling_80k/` compares
0.5×, 1×, and 2× **area** with 10 matched random seeds (19–28) per group.
The source XML is left unchanged; generated XML files are stored with each run.
For area multiplier `a`, the target area is `1200*a`, the initial radius and
target perimeter scale by `sqrt(a)`, and the area-constraint strength is
`0.0004/a`. This preserves the effective area modulus under the project's CPM
scaling convention. Local biochemical concentrations, reaction rates,
diffusivities, lattice spacing, membrane resolution, surface-constraint
strength, and CPM timing remain fixed. The 80,000 value is Morpheus `StopTime`,
not a count of CPM sweeps; the model uses `MCSDuration=0.5`.

The [final report](results/area_scaling_80k/report/index.html) contains the
30 centroid tracks, seed-level and group-mean MSD curves over lags of 5–20,000
model time units, summary tables, and one native-frame GIF per group. It is
generated only after all 30 runs complete and validate. The run manifest
records the source and binary checksums and the exact scaling rule.

Reproduce from the repository root:

```bash
/home/4473331/.conda/envs/scaling_sandbox/bin/python literature_models/M7682/scripts/prepare_area_campaign.py
sbatch literature_models/M7682/scripts/run_area_campaign.slurm
# After all array tasks complete:
/home/4473331/.conda/envs/scaling_sandbox/bin/python literature_models/M7682/scripts/report_area_campaign.py
```
