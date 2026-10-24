# Competition, coexistence and coarsening of polarity domains

Serves Fig 2 (routes to resisting replacement), Fig 3 (resolution of coexisting domains) and Fig 4E (multiple domains, births versus takeovers). Drafted 2026-10-23.

## Established

- **Saturation sets the speed of competition.** In mass-conserved activator–substrate models, peaks far below a saturation point compete quickly to a single winner; peaks near it compete slowly and become metastable. One factor explains most of the variation in competition time across parameters [@chiou2018competition].
- **The rule holds in cells.** Raising polarity protein levels in yeast gives multiple coexisting sites, as predicted [@chiou2021number].
- **Substrate depletion limits sites.** With total protein fixed, membrane accumulation depletes cytoplasmic substrate, limiting the size, peak concentration and number of active regions; peaks with more content recruit substrate better [@chiou2021number].
- **Coarsening is driven by mass competition.** In two-component mass-conserving systems, self-amplifying mass transport between neighbouring domains drives coarsening, which is generically uninterrupted; weak source terms can arrest or reverse it [@brauns2021wavelength].
- **Domains have effective interface laws.** Cyclic attachment and detachment at the interface act like an interfacial tension, giving non-equilibrium Gibbs–Thomson and related laws [@weyer2025interface].

## Where we extend or depart

- **Size as the driver of saturation.** Fig 2A: Wave-pinning, Holmes, deBelly and Spring hold the extra active material of a larger cell by broadening at a fixed plateau. Goryachev raises its peak. In the saturation rule's terms, the first group stays at saturation while getting wider; Goryachev moves away from it. Neither the rule's authors nor others checked frame this as a size effect.
- **Two of our Fig 2 routes are size versions of known mechanisms.** Goryachev's slower initiation in larger cells matches substrate depletion [@chiou2021number]. The expansion constraint matches saturation [@chiou2018competition]. What is new is that each strengthens with size at fixed kinetics.
- **Secondary domains and heading.** Fig 4E: larger rings spawn more secondary domains than they lose incumbents, and where multiple domains predominate, drift plus counted relocations under-predict heading loss. No paper checked links coexistence to heading persistence or to its size scaling.

## Our interpretation (I)

- **The saturation rule predicts the direction of the Fig 4E trend in each model.** Saturated, broadening models (Wave-pinning, Holmes, Spring) should keep secondary domains longer in larger rings; Goryachev, whose peak moves away from saturation, should resolve them quickly. Fig 4E shows multiple domains becoming more frequent with size in the first group and rarer in Goryachev. State this as consistency, not proof.
- **Separation should matter little once domains are apart.** Mass exchange in these models runs through a fast-diffusing cytosol, so competition is close to global. That would explain Fig 3's finding that separation acted mainly by setting the turning regime. Check against the coarsening law in [@brauns2021wavelength] before stating it.
- **Coexistence may be where heading is lost in strong input.** If secondary domains persist without winning, heading moves through shifts in the balance between domains, which dominant-domain tracking cannot register (narrative, Section 4). The saturation rule makes such long-lived coexistence expected in the broadening models.

## Implications for the manuscript

- **Pedagogical:** saturation rule and substrate depletion, cited as the basis of the Fig 2 routes.
- **Novel:** size dependence of each route; the per-model direction of the multiple-domain trend; coexistence as a channel of heading loss that erodes the size advantage.
- **Hypothesis-generating:** birth rate, position and fate of secondary protrusions as a cellular readout (`measurement_methods.md`).

## To check

- Whether [@chiou2018competition] or [@brauns2021wavelength] give competition time as a function of domain width or system length. If so, compare directly with Fig 2C and Fig 3.
- Whether the two-component coarsening result covers Holmes model 3 and the mechanochemical models.
- Weyer, Brauns & Frey (arXiv:2203.03892), the longer coarsening and wavelength-selection paper. See `queue.md`.
