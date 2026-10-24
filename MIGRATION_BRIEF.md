# Brief: restructure the literature folder

Replace `literature_review/` with `literature/`, using the attached package. The format is defined in `literature/README.md`. Follow it.

The attached package already contains:
- `literature/README.md`: the format.
- `literature/cards/`: 22 drafted cards from a literature check on 2026-10-23.
- `literature/syntheses/`: seven drafted syntheses.
- `literature/positioning.md`: claims mapped to the literature, from the round-3 narrative.
- `literature/queue.md`: papers still to read.
- `literature/references.bib`: entries for the drafted cards.
- `literature/tools/build_index.py` and a generated `literature/index.csv`.

Copy these in as they are. Do not rewrite their content beyond what this brief asks.

## 1. Moves

Use `git mv` throughout so history follows the files.

| From | To | Then |
|---|---|---|
| `literature_review/xuan_2018.md` | `literature/evidence/xuan_2018/dossier.md` | Write card `xuan2018pgcc` |
| `literature_review/xuan_2020.md` | `literature/evidence/xuan_2020/dossier.md` | Write card `xuan2020vimentin` (pick the key from the title) |
| `literature_review/mallin_2023.md` | `literature/evidence/mallin_2023/dossier.md` | Write card |
| `literature_review/jemaa_2017.md` | `literature/evidence/jemaa_2017/dossier.md` | Write card |
| `literature_review/jemaa_2023.md` | `literature/evidence/jemaa_2023/dossier.md` | Write card |
| `literature_review/hosny_2026.md` | `literature/evidence/hosny_2026/dossier.md` | Write card |
| `literature_review/migration_data/<study>/` | `literature/evidence/<study>/` | Keep each study's `extraction/` and any `source/`, `figures/` directories as they are |
| `literature_review/migration_data/{README.md, EXTRACTION_PLAN.md, EXTRACTION_STATUS.md, OU_FITTING_DISPOSITION.md, manifest.csv, relocation.csv, studies.csv, audit/, .gitignore}` | `literature/evidence/` | |
| `literature_review/bosgraaf_2009.md` | delete after checking | Already carded as `bosgraaf2009pseudopodia`. Merge anything the card misses |
| `literature_review/mitternwallner_2020.md` | delete after checking | Already carded as `mitterwallner2020gle`. Merge anything the card misses |
| `literature_review/various_polarity_references.md` | `literature/archive/cpm_calibration/various_polarity_references.md` | Write stub cards for Wu 2009 (photoactivatable Rac) and de Seze 2025 (optogenetic RhoA GEF), role `method`, themes `measurement_methods`. Take metadata from the file and check it against the publisher pages |
| `literature_review/deSeze_2025/` | `literature/archive/cpm_calibration/deSeze_2025/` | |
| `jilkine_ref.md` (repo root) | `literature_models/specs/jilkine_2011_spec.md` | Write card `jilkine2011comparison` linking it |
| `tap_ref.md` (repo root) | `literature/archive/tao_2020_spec.md` | Default: archived. If anything live still uses the Tao model, move it to `literature_models/specs/` instead and say so in the commit message |

Give each archived file a first line: `> Archived <date>. Written for <purpose>. <why retired>.` The CPM notes were written to calibrate the 2D Cellular Potts branch, retired 2026-10-17. The Tao note's purpose is unknown, so say that.

## 2. Path references

Before moving `migration_data/`, grep the whole repo, code included, for `literature_review`. Update every reference: `wgd_langevin/`, `DATA_MANIFEST.md`, `RETENTION.md` files, the extraction `METHODS.md` files, `reports/`.

After the move, rerun the Fig 5 and OU rebuild entry points and confirm the outputs are unchanged. If you cannot confirm that, revert the `migration_data/` move only, leave `literature_review/migration_data/` in place with a README pointing to `literature/evidence/`, and record why in the commit message. Everything else still moves.

## 3. New cards from existing material

Write cards for these. Use only the repo material named and the publisher pages. Set `access` and `status` to match what you actually read.

- **The six evidence studies**, from their dossiers. Role `evidence`, readout `observed_migration`, theme `wgd_migration_evidence`. The Results section is a short pointer list: the comparison made, the timescales, and which dossier sections hold the migration numbers. Do not copy tables from the dossiers.
- **The five manuscript models**: Wave-pinning (Mori, Jilkine & Edelstein-Keshet 2008), Goryachev (Goryachev & Pokhilko 2008), Holmes model 3, de Belly 2026, and Spring. Take the source papers from the `literature_models/` provenance READMEs and `MODEL_REGISTRY.md`. Role `background`, type `model`. Link the native implementation folder. If a model's source paper is ambiguous in the repo (Holmes model 3, Spring), write the card with what the repo states and add a `Caveats` line saying the source needs confirming.
- **Buttenschön, Liu & Edelstein-Keshet 2020**, "Cell size, mechanical tension, and GTPase signaling in the single cell" (Bull Math Biol 82:28). The repo already has `literature_models/buttenschoen_2020_moving_domain/` with `run_buttenschoen_wgd.py`. Read the paper and that folder's README and results. Report in the card what the paper says about size, and what the repo's WGD run found. This is top priority in `queue.md` because it may bear on our novelty claims.

Add bib entries for every new card. Update the syntheses' citation lists where a new card belongs, without rewriting their prose.

## 4. Repo wiring

- Delete the empty `literature_review/` once it is empty.
- Add one line to the root `README.md` key-paths list: `literature/`: cards, syntheses and claim positioning for the manuscript; published migration data under `literature/evidence/`.
- Run `python literature/tools/build_index.py`. It must exit cleanly.

## 5. Do not

- Invent bibliographic metadata. If a field can't be confirmed, leave it out of the bib entry and add `% unverified: <field>`.
- Raise a card's status without the reading it requires.
- Edit the syntheses' interpretive content or `positioning.md` beyond adding keys. Those are Richard's to revise.

## Deliverables

One commit (or one per section) on a branch off the current working branch. Report: the moves made, the path references updated, the Fig 5 rebuild check result, the new cards written with their status, and anything left undone with the reason.
