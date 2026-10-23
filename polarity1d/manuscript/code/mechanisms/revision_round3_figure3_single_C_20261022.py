"""Landscape Figure 3: C=30 profiles at left; four-model B/C at right.

Reuse every scientific coordinate and completed native protocol. Preserve the
preceding four-placement artwork and caption in the supplement and old edition.
"""
from pathlib import Path
import shutil

import revision_round3_figure3_20261022 as base


def main():
    previous=base.OUT
    destination=previous.parent/'figure3_single_C'
    data=destination/'supporting/data'
    metadata=destination/'supporting/metadata'
    data.mkdir(parents=True,exist_ok=True);metadata.mkdir(parents=True,exist_ok=True)
    for source in (previous/'supporting/data').iterdir():
        if source.is_file():shutil.copy2(source,data/source.name)
    for name in ('collection.json','source_sha256.json','topology_provenance.json'):
        shutil.copy2(previous/'supporting/metadata'/name,metadata/name)
    base.OUT=destination;base.DATA=data;base.META=metadata
    base.MAIN_EXAMPLES=base.EXAMPLES[:2]
    base.SINGLE_C_LAYOUT=True
    base.EXTRA_SUPPLEMENT_SOURCES=[previous/'figure3_round3.pdf',
        base.ROOT/'manuscript_outline/revision_round3_figure3_20261022/figure3_legend.pdf']
    base.render()
    notes=destination/'README.md'
    notes.write_text(notes.read_text().replace(
        'python -B polarity1d/manuscript/code/mechanisms/revision_round3_figure3_20261022.py collect\n','').replace(
        'revision_round3_figure3_20261022.py render','revision_round3_figure3_single_C_20261022.py'))
    base.write(metadata/'layout_source_preservation.json',dict(
        generated_at=base.native.now(),native_integrations=0,
        previous_main=str((previous/'figure3_round3.pdf').relative_to(base.ROOT)),
        previous_main_sha256=base.sha(previous/'figure3_round3.pdf'),
        source_data_sha256={str(p.relative_to(previous)):base.sha(p)
            for p in (previous/'supporting/data').iterdir() if p.is_file()},
        preserved_main_profiles=['i','ii','iii','iv'],displayed_main_profiles=['i','ii'],
        main_sizes=[30.],curve_sizes_unchanged=True))
    print(destination)


if __name__=='__main__':main()
