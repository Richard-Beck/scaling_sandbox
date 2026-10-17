# Experiment configuration

Configurations are organized first by the number of interacting polarity
species and then by the role of the experiment:

```text
configs/
  single_species/
    parameter_mapping/
    calibration/
    testing/
  three_species/
    parameter_mapping/
    calibration/
    testing/
```

`parameter_mapping` records assumptions and calculations used to transfer a
published parameterization into the moving 1D geometry. `calibration` varies
parameters against a fixed experimental panel. `testing` applies an established
parameterization to a growing library of stimulus protocols.

Names such as `chiou`, `holmes`, `eroume`, and `mori` identify parameterization
provenance only. They are leaf-level metadata, not organizing principles.

Files under `_standards/` define the common panels and metrics. They are
documentation/configuration standards and are not submitted as campaigns.
The single-species campaigns implement the v2 standards directly; the
three-species campaigns have not yet been migrated to them.

## Focused motile length sweep

`single_species/testing/focused_motile_length.json` and
`three_species/testing/focused_motile_length.json` form a matched, motile
length-scaling study.  Both use reference lengths
`5,7.5,10,12.5,15,20,25,30,40,50` micrometers, 500-second trajectories, soft
and stiff endpoint mechanics, and four cue protocols: a single cue, exactly
balanced dual cues, dual cues with a 1% bias, and polarity reversal.  Each
configuration is an unthrottled 80-task SLURM array with a dependent validation
and plotting job.  The single-species chemistry uses `mu=2.1875`; the
three-species chemistry uses the finitely polarizable Holmes point
`kappa_R=0.04`, `n=3`.

The campaign varies axial length at fixed concentrations and `h=1`; it is not
an isometric volume-scaling experiment.  Task outputs include compressed
trajectories so the merge can render laboratory-frame and material-frame
kymographs as well as length--polarity, interval-speed, displacement, and
deformation summaries.
