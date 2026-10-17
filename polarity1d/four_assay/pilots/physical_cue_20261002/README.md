# Physical-gradient response control: A*x

Launched 2026-10-02 with `../../physical_cue.py`. Both active response assays
reuse the exact complete revision-11 preparation records and checkpoint files
from `../../results/records.json.gz`; no conditioning or release is rerun.
The original campaign and accepted baseline are preserved.

- Right response cue: A*x; left crossfade endpoint: A*(L-x).
- The 20-unit reversal crossfade interpolates those two physical cues.
- A retains each model's existing value. Peak input is now A*L, mean A*L/2;
  this is a fixed physical slope control, not a fixed peak or total input control.
- Spring uses instantaneous physical length and position in native author units.
- Conditioning/release columns retain the original A*x/L preparation clocks.
- Original length grids, finite-duration grids, recovery, eligibility, response
  gates, solvers and physical model parameters remain unchanged.
- Finite endpoint remains the shortest passing TESTED duration. A success at the
  smallest tested pulse is bounded by the existing grid's lower resolution;
  no smaller durations are automatically introduced.

SLURM response array: 20677535 (3677 chunks); automatic finalization: 20677536.
QoS small, one CPU/3 GB per chunk, no array concurrency throttle. Finalization
checks exact preparation reuse, checkpoint hashes, complete duration coverage,
and independently validates the 28 plotted regression series.

Outputs: `results/length_response_and_settling.png/pdf/svg`, the companion slope
heatmap, `comparison_to_normalized.csv`, complete records and validation JSON.
`results/status.json` and `work/launch.json` distinguish running from finished.
The final figure title identifies A*x responses and original A*x/L preparations.

Completed 2026-10-02 17:34 UTC. All 308 cases and 34150 finite trials passed
coverage checks, with 234 successful, 50 censored and 24 ineligible readouts.
All 14 settling fits equal their original values; all 28 regression series were
independently checked (26 eligible fits). No execution errors occurred.

Descriptive countercue changes: Spring's fitted slope becomes negative; LEGI
is nearly flat over successful lengths, with six previously successful long
lengths now censored. Wave-pinning, Goryachev, Holmes 3 and deBelly retain or gain
positive fits. Otsuji still has no successful responses. Goryachev's successful
length range expands substantially; comparisons of whole-range slopes therefore
also reflect changed fit membership. Eight countercue successes are lost overall
(six LEGI, two deBelly), while fourteen are gained (Goryachev).

At nominal L=1, fixed-domain cues reproduce the original readouts. Spring is not
an identical-cue control at resting L0=1: its instantaneous physical length
changes, so A*x and A*x/L still differ during the response.
