# Cell size and polarity

Serves the framing of Figs 1–4, the Methods justification for fixed kinetics, and the Discussion. Drafted 2026-10-23.

## Established

- **Polarity kinetics need not scale with size, so size changes polarity.** PAR protein kinetics and boundary length scales fail to scale with cell size. In every model examined, including wave-pinning and Goryachev, this gives a size threshold below which polarity cannot be held, and small germ-lineage cells lose polarity at the predicted sizes [@hubatsch2019threshold].
- **More polarity protein gives more polarity sites.** In yeast, raising key polarity proteins gives multiple sites and simultaneous budding, and large cytokinesis-defective cells start with several sites [@chiou2021number]. The theory: competition slows near a saturation point [@chiou2018competition].
- **Size can act through shape and mechanics, in either direction.** In one theory, larger contact area eases polarization, giving faster and more persistent cells [@ron2025size]. In a mechanical cell model, larger cells become rounder, slower and less persistent because adding membrane costs energy [@thuroff2019bridging].
- **Copy number has been varied, but not size at fixed concentration.** Stochastic wave-pinning studies vary molecule number at fixed geometry [@kryvoruchko2026intrinsic]. Walther et al. may have varied width with molecule number [@walther2012stochastic]; not confirmed.

## Where we extend or depart

- **The large-size end.** The known size effect at fixed kinetics is loss of polarity in small cells [@hubatsch2019threshold]. Ours concerns large cells: established polarity becomes harder to replace (Figs 1–3) and heading lasts longer under fluctuating input (Fig 4).
- **Persistence and reorientation as functions of size, across models.** Not found in the papers checked. The closest are reorientation studies at a single size [@buttenschoen2022repolarization] and copy-number studies at a single size [@kryvoruchko2026intrinsic].
- **Concentration held fixed.** Chiou et al.'s large cells carry more protein as well as more area [@chiou2021number]. Our models hold concentration and kinetics fixed while length changes, which isolates size.

## Our interpretation (I)

- Hubatsch et al. give the empirical licence for our "no retuning" design: measured polarity kinetics adapted little, if at all, to size. Cite it in Methods.
- Size effects reported in the literature have both signs, depending on how size enters. Scope our claim to polarity biochemistry on a fixed 1D geometry, where all five models give the same sign.
- WGD changes both size and copy number. Our models capture size at fixed concentration; intrinsic-noise effects of higher copy number would come on top (`noise_driven_domain_motion.md`).

## Implications for the manuscript

- **Novel:** the main claim, that larger size at fixed kinetics makes established polarity harder to replace and more persistent, in five models.
- **Pedagogical:** the threshold result and the saturation rule as the context in which size matters.
- **Discussion:** shape and mechanics can reverse the sign [@thuroff2019bridging] or add to it [@ron2025size].

## To check

- **Top priority:** Buttenschön, Liu & Edelstein-Keshet 2020, "Cell size, mechanical tension, and GTPase signaling in the single cell". The repo already implements it (`literature_models/buttenschoen_2020_moving_domain/`, including `run_buttenschoen_wgd.py`). It may already address size in a polarity model and bears directly on the novelty claim.
- Whether our smallest sizes (L19, C20) are far enough above each model's polarity threshold.
- Walther et al. full text (width design).
