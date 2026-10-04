# Polarity manuscript workflow

For the minimal Git codeset, start with [CODESET.md](CODESET.md): environment,
file map, verification, figure regeneration and fresh SLURM runs. Git retains
the 291 exact reversal input checkpoints and compact figure inputs; bulk trial
histories and the historical personal explorer remain local.

This is the single working folder for the polarity manuscript. The comparison
uses five models: Wave-pinning, Goryachev, deBelly, Holmes 3 and Spring. Otsuji is
excluded because its initial polarization was not clear; pure LEGI is excluded
because proper basal withdrawal removes its polarity.

The new reversal comparison starts from qualifying **corrected finite-patch
establishment endpoints**, at nominal cue amplitude. Its right challenge mirrors
the left preparation patch, using fixed widths 1.5 and 3 or 20% of instantaneous
length (Spring fixed widths 0.1 and 0.2 in native units). Historical
linear-conditioned reversal results do not supply this figure's length trends.

| Location | Contents |
| --- | --- |
| [code/](code/) | Native solvers, corrected experiment, figure builders, tests and independent audits |
| [figures/index.html](figures/index.html) | Main figure, detailed reversal supplement, corrected establishment heatmaps and personal field explorer |
| [data/establishment_basal/](data/establishment_basal/) | Corrected establishment records and selected native checkpoints/histories |
| [data/runs/](data/runs/) | Fresh five-model matched reversal campaign and compact trial evidence |
| [data/reversal_fields/](data/reversal_fields/) | Retained personal diagnostic only, including its original experiment identity |
| [data/provenance/](data/provenance/) | Source/copy hashes, validation and explicit path mappings |

Read [PROTOCOL.md](PROTOCOL.md) and [ASSUMPTIONS.md](ASSUMPTIONS.md) before
interpreting the results. [REPRODUCE.md](REPRODUCE.md) gives the commands.
[RETENTION.md](RETENTION.md) marks the older seven-model suite and exploratory
campaigns deprecated. No unrelated workstream is included in this cleanup.

The main explanation and length summary, and the detailed width supplement,
were generated from the completed fresh campaign after all 73,975 trials passed
independent audit. Across 660 readouts, 375 succeeded, 207 were censored within
the tested range, and 78 were preparation-ineligible. All 291 qualifying no-cue
controls passed. The landing page links final data-generation time and state.
The corrected establishment heatmaps are retained independently. The interactive
old field explorer is for personal inspection; it is outside the manuscript
comparison.
