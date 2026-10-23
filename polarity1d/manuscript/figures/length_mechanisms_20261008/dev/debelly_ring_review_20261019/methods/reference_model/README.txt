Historical reference model sources
==================================

These are exact frozen sources for the viewed assays, retained for inspection.
They are not required to open or regenerate the viewer. No production campaign
runner, alternative cap model, new operator or full simulation history is kept.

common/ring_core.py contains the shared native staged force, MCA and Rac/Rho
updates. common/debelly.py and common/polarity1d contain its native parameter,
rate and contract dependencies. common/peak_detector.py is the original peak
detector. original/gap_events.py is the original gap-insertion policy.
face/gap_events.py with face/native_gap_events.py is the existing-face policy
used by A and B. All exact source hashes and historical identifiers are in
source_manifest.json. Importing the model for independent work requires NumPy
and SciPy and the common plus selected policy directory on the import path.

basal_state.npz and basal_state.json are the exact shared original qualified
state. Their JSON retains the original preparation clock 18,060 s, next_id and
all four flank controller values. The viewed challenger assays reset protocol
time to zero while keeping that full state. parameters.json records the shared
physical parameters, numerical step and cue geometry; its prose is packaging
metadata, not a changed model law. No new replay was run during consolidation.

The package supports exact regeneration of the retained interactive viewer
from viewer/raw. The archived audits document the earlier scientific replay
checks. Recreating complete original campaigns is outside this retained review.
