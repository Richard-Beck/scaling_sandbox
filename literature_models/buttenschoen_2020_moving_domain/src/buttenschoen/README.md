# Buttenschön spatial mechanochemical model

This package implements Eqs. 12a--e of Buttenschön, Liu and
Edelstein-Keshet (2020): mass-conserved active/inactive GTPase on a moving 1D
cell coupled to Kelvin--Voigt endpoint mechanics.

The `figure9_regimes()` presets cross contractile/protrusive force coupling
with stiff/soft mechanics. Figure 9 specifies the elastic moduli and total
GTPase amounts, but not every numerical input. The reconstruction uses values
within Tables 1--2 for the remaining chemistry and mechanics: gamma=25,
D=0.1, Di=100, Fmax=50000, mu=100000, eta=60, Gc=1.25, and m=n=4.
Therefore these are explicit, editable **reconstruction presets**, not a
pixel-exact reproduction. Every assumption is recorded in `metadata.json`.

Run all four presets from the repository root:

```bash
conda run -n scaling_sandbox python scripts/run_buttenschoen_spatial.py
```

Use `--help` to select a single regime, resolution, duration, or output path.
Each regime produces a compressed trajectory, metadata with conservation and
solver diagnostics, and a physical-coordinate active-GTPase kymograph.

Run the matched 1x/2x-volume companion experiment with:

```bash
conda run -n scaling_sandbox python scripts/run_buttenschoen_wgd.py
```

For linear scale `s = 2^(1/3)`, the 2x condition applies `L0 -> s L0`,
`E -> s E`, `eta -> s eta`, `mu -> s^2 mu`, and `Fmax -> s^2 Fmax`.
The native 1D GTPase amount scales by `s`, preserving initial concentrations;
all biochemical kinetics and physical diffusion coefficients are unchanged.
Both members of every pair start with physical centroid `xc(0) = 0`, and the
paired kymographs use identical physical-position limits. The companion runner
uses a 720-second window so the delayed recoil of the scaled soft protrusive
condition is retained.
