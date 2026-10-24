# literature/

Literature for the WGD manuscript: what each relevant paper did, how it bears on our claims, and where our work is established, extends prior work, or is new.

The manuscript aims and claims live in the steering docs (manuscript vision; narrative and claims). This folder serves them. It does not set direction.

## Layout

| Path | Holds | Rule |
|---|---|---|
| `cards/<key>.md` | One card per paper. Facts about that paper only. | Fixed template (below). One screen of body text. |
| `syntheses/<theme>.md` | Interpretation across papers, one theme each. | Every literature claim cites a card key. Our own interpretation is marked. |
| `positioning.md` | Each manuscript claim mapped to the literature: established, extends, novel, contradicted, or untested. | One row per claim in the narrative doc. Updated when the narrative changes. |
| `queue.md` | Papers to read, each with the reason and priority. | A paper leaves the queue when it gets a card. |
| `references.bib` | One BibTeX entry per carded paper. | Keys match card filenames. Metadata copied from the publisher or arXiv page, never from memory. |
| `index.csv` | Generated from card frontmatter. | Never edit by hand. Rebuild with `tools/build_index.py`. |
| `evidence/` | Published migration data we extract and re-analyse: per-study dossiers and extraction packages. | Data and provenance live here, not in cards. |
| `tools/build_index.py` | Builds `index.csv`; checks frontmatter and bib keys. | Run before every commit that touches cards. |
| `archive/` | Notes written for retired work. | History only. Each file gets a header line saying what it was for and why it was archived. |

Native implementations of literature models stay in `../literature_models/`. A model paper's card links to its folder there.

## Keys

`<firstauthor><year><word>`, lower case ASCII, no punctuation: `kryvoruchko2026intrinsic`, `prentice2016memory`. The card filename is the key. The bib entry uses the same key.

## Card template

```markdown
---
key: walther2012stochastic
title: Deterministic versus stochastic cell polarisation through wave-pinning
authors: Walther GR; Marée AFM; Edelstein-Keshet L; Grieneisen VA
year: 2012
venue: Bull Math Biol 74:2570-2599
doi: 10.1007/s11538-012-9766-5
url: https://pmc.ncbi.nlm.nih.gov/articles/PMC3480592
type: model
readout: [simulated_polarity]
size_varied: copy_number
kinetics_fixed: yes
noise: intrinsic
role: [precedent]
themes: [noise_driven_domain_motion, size_and_polarity]
sections: [fig4]
access: abstract
status: stub
checked: 2026-10-23
---

# Walther et al. 2012: stochastic wave-pinning

**Why it matters to us.** One line.

## What they did
System or model, geometry, what was varied, noise source.

## Results
- [R] Reported result.
- [F] Our reading of a plotted value.

## What transfers
What carries over to our models or cells, and what does not.

## Relation to our claims
Supports / precedent for / competes with / contradicts, naming the figure or narrative section.

## Caveats
Limits of the paper or of our access to it.
```

### Frontmatter fields

| Field | Values |
|---|---|
| `type` | `theory`, `model`, `experiment`, `review`, `method`, `data` |
| `readout` | Any of `observed_migration`, `ou_inference`, `simulated_polarity`, `observed_polarity`, `none`. The first three are the manuscript's readouts, which are always kept distinct. |
| `size_varied` | `none`, `persistence`, `reorientation`, `threshold`, `number_of_sites`, `copy_number`, `shape`, `other` |
| `kinetics_fixed` | `yes`, `no`, `na` (whether local kinetics were held fixed when size changed) |
| `noise` | `intrinsic`, `extrinsic`, `both`, `none`, `na` |
| `role` | Any of `background`, `precedent`, `competitor`, `counterpoint`, `evidence`, `method` |
| `themes` | Synthesis file stems |
| `sections` | `fig1` to `fig6`, `methods`, `discussion` |
| `access` | `full_text`, `summary` (machine summary of full text), `abstract`, `secondary` (another paper or an existing note) |
| `status` | `stub`, `read`, `verified` |
| `checked` | YYYY-MM-DD of the last check |

### Status

- `stub`: metadata checked against the publisher or arXiv page; content from the abstract, a summary or a secondary source.
- `read`: the card's author read the full text, or the sections the card relies on.
- `verified`: a second pass checked every R, P and F item against the full text.

Never raise status without doing the reading it requires.

### Provenance tags

Same convention as the existing dossiers:

- **R**: reported in the paper's text, abstract or caption.
- **P**: printed on a figure.
- **F**: our approximate reading of a plot or image.
- **D**: our calculation from a reported value.
- **I**: our interpretation. Not in the paper.

Every bullet under Results carries R, P, F or D. Anything under Relation to our claims is I unless tagged otherwise.

## Syntheses

One file per theme. Each opens with the manuscript figures it serves, then:

1. **Established.** What the literature already shows, with card keys.
2. **Where we extend or depart.** What our results add, with figure references.
3. **Our interpretation.** Marked (I). Hypotheses, candidate nulls, links between papers.
4. **Implications for the manuscript.** How to present the result: pedagogical, novel, or hypothesis-generating (the three roles in the vision).
5. **To check.** Open questions and queue items.

Cite cards inline as `[@key]`. Two pages at most.

## Adding a paper

1. Copy the template to `cards/<key>.md`. Fill metadata from the publisher or arXiv page.
2. Add the BibTeX entry to `references.bib`.
3. Write the body at the level your access supports. Set `access` and `status` honestly.
4. Add the key to the relevant syntheses and, if it changes how a claim is positioned, to `positioning.md`.
5. Remove it from `queue.md`.
6. Run `python tools/build_index.py`. Fix any errors it reports.

## Keeping it useful

- Cards hold facts about one paper. Interpretation goes in syntheses.
- A synthesis that no manuscript figure needs gets archived.
- When the narrative doc changes a claim, update `positioning.md` in the same pass.
- Dossiers under `evidence/` stay as long as they need to be. Their cards stay short and point to them.
