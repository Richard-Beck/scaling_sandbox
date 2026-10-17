# Holmes 2012 Rac/Cdc42/Rho model

The reduced six-field implementation is
`../../migration1d/src/migration1d/holmes.py`; its persistence and WGD-motility
campaigns use the generic migration runner. This folder retains the earlier
size-scaling report and avoids a duplicate solver copy.

Native bulk-surface Holmes code is shared with the Giese bundle because it
uses that bundle's mesh and coupling engine.
