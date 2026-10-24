# Reading queue

A paper leaves the queue when it gets a card. Citations here are as found and not yet checked against publisher pages; confirm before carding.

## Priority 1: could change a novelty claim

- **Buttenschön A, Liu Y, Edelstein-Keshet L (2020). Cell size, mechanical tension, and GTPase signaling in the single cell. Bull Math Biol 82:28.** Size in a GTPase polarity model. The repo already implements it (`literature_models/buttenschoen_2020_moving_domain/`, with `run_buttenschoen_wgd.py`). Read the paper and the folder's results.
- **Walther et al. 2012, full text** (card `walther2012stochastic`, currently a stub). Confirm whether width was varied with molecule number, and whether pinning-position fluctuations are reported against width.
- **Altschuler SJ, Angenent SB, Wang Y, Wu LF (2008). On the spontaneous emergence of cell polarity. Nature 454:886–889.** Noise-driven single-cluster polarity; standard citation for stochastic polarity.

## Priority 2: theory and precedent for Figs 2–4

- **Armero J et al. (1996).** Original effective front-diffusion result under external noise, cited by [@rocco2001front]. Find the exact citation.
- **Meerson B, Sasorov PV, Kaplan Y (2011). Phys Rev E 84:011147.** Front diffusion coefficient from intrinsic shot noise in the weak-noise limit.
- **Hecht I, Kessler DA, Levine H (2010). Transient localized patterns in noise-driven reaction-diffusion systems. Phys Rev Lett 104:158301.** Background for transient secondary domains.
- **Weyer H, Brauns F, Frey E (2023). Coarsening and wavelength selection far from equilibrium (arXiv:2203.03892).** Longer treatment of coarsening in nearly mass-conserving systems; may give competition time against domain size.
- **Kuwamura M, Seirin-Lee S, Ei S (2018). Dynamics of localized unimodal patterns in reaction-diffusion systems for cell polarization by extracellular signaling. SIAM J Appl Math 78:3238–3257.** Peak motion under external signals; relevant to the Fig 3 displacing regime.
- **Hadjitheodorou A et al. (2021). Directional reorientation of migrating neutrophils is limited by suppression of receptor input signaling at the cell rear through myosin II activity. Nat Commun 12.** Cellular limits on reorientation; Figs 1–3.
- **Skoge M et al. (2014). Cellular memory in eukaryotic chemotaxis. PNAS 111:14448–14453.** Memory under changing gradients; Fig 2 memory class.
- **Jilkine A, Edelstein-Keshet L (2011). A comparison of mathematical models for polarization of single eukaryotic cells in response to guided cues. PLoS Comput Biol 7:e1001121.** Source of the Wave-pinning and Goryachev implementations; card to be written from `jilkine_ref.md`.

## Priority 3: background

- **Singh AR, Leadbetter T, Camley BA (2022). Sensing the shape of a cell with reaction diffusion and energy minimization. PNAS 119:e2121302119.** Wave-pinning domains minimize perimeter at roughly fixed area.
- **Hughes JM et al. (2025). A mass conserved reaction-diffusion system reveals switching between coexisting polar and oscillatory cell motility states (arXiv:2504.07446).**
- **Copos C, Mogilner A (2020). A hybrid stochastic-deterministic mechanochemical model of cell polarization. Mol Biol Cell 31:1637–1649.**
- **Gamba A et al. (2005).** Stochastic phase-separation model of chemotactic symmetry breaking (PNAS), cited by [@defranciscis2012bounded].
- **Wu YI et al. (2009)** photoactivatable Rac and **de Seze J et al. (2025)** optogenetic RhoA GEF. Citations are in the archived CPM notes; card as `method`.
