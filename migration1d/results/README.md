# Generated migration outputs

Campaign outputs are ignored by Git and can be recreated from `../configs/`.
New runs are grouped as `<chemistry_class>/<stage>/<experiment>/`.

The paired WGD-null screens produce separate reports at
`single_species/testing/wgd_null/report.html` and
`three_species/testing/wgd_null/report.html`. Their `shards/` directories are
validated intermediate products; the merged CSV, JSON summary, figures, and
HTML report are the scientific outputs.

The authoritative pre-cleanup output snapshot remains under
`../../reorg/legacy/results/`. Flat `migration1d_*` directories in this active
directory predate the chemistry/stage organization and are not authoritative.
