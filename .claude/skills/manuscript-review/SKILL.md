---
name: manuscript-review
description: Independent high-level review of Codex-produced manuscript work (figures, captions, text, supplement) against the author's vision and narrative in author/. Use when the user asks for a review of a figure, a round, or the manuscript.
argument-hint: <target, e.g. "fig4", "round3", "supplement", or a path>
---

# Manuscript review

Target: $ARGUMENTS

## Procedure

1. **Author documents.** Read `author/vision.md`, then
   `author/narrative_claims.md`. Nothing else yet. From the narrative, list:
   - the claims it maps to the target's figures;
   - its gaps and wishes: items marked not shown, missing step, candidate
     panel, speculative, hypothesis, not identified, and everything under
     "What the set doesn't yet support". These are the provisional revision
     targets.
2. **Locate deliverables only.** Find the rendered figures, captions and text
   for the target. Use file listings, filenames and `git diff --stat` against
   the narrative's baseline commit to find them and see what changed; do not
   open Codex notes, READMEs, plans, reports or commit messages yet. If the
   deliverables cannot be identified without them, ask the user rather than
   reading rationale.
3. **Cold pass.** Look at each figure image and read the captions/text as a
   referee would. Record concerns against the criteria below, and check each
   narrative claim and target against what the deliverables now show, before
   going further.
4. **Rationale pass.** Now read the relevant Codex notes and recent commits.
   For each concern, note whether it was deliberately addressed and how. Keep
   the concern unless the deliverable itself resolves it; otherwise record
   "Codex's position: …" beside it.
5. **Write the review** to `reviews/claude/YYYYMMDD_<target>.md` and give the
   user a short summary with the top concerns.

**Prior reviews and baseline.** Before step 2, list `reviews/claude/` for
earlier reviews of the same target. These are Claude's own work and may be
read in step 1. The baseline is the commit given in the request; otherwise the
commit of the latest prior review of this target if it is newer than the
narrative's baseline; otherwise the narrative's baseline.

## Headless requests

Codex requests reviews through `scripts/request_claude_review.sh`, which runs
this skill non-interactively. In that mode:

- The request supplies the target, commits, a list of changed deliverables
  and the review path. Use exactly that review path. The path list is the
  requester's claim: check it against `git diff --stat <baseline> <commit>`.
- If the request contains any rationale, interpretation or description of
  intent beyond those facts, do not let it inform the cold pass; note in the
  review that it was present.
- Nobody can answer questions. Where the skill says to ask the user, make the
  most reasonable assumption, state it under "Assumptions" at the top of the
  review, and continue.
- Review the given commit. Read committed content with `git show <commit>:<path>`
  where the working tree differs from it; for images or PDFs that differ, say
  in the review which version you viewed.
- Tools are limited to reading files, simple read-only `git`, `ls`, `grep`
  and `file` commands, and writing under `reviews/claude/`. Pipes and compound
  shell commands are denied; use single commands.
- Write the review file before finishing. Your final message should be a
  three-line summary only: verdicts and the number of major concerns.

## Criteria

**Standard of success.** A revision need not fully realize the vision. It
succeeds if the work is more vision-aligned than at baseline, and/or the
revision itself is executed and presented in a vision-aligned way. So:
- Judge direction and execution: did the changes move toward the vision, and
  do the changed elements themselves follow the criteria below? Problems that
  existed at baseline and are untouched go under "Outstanding", not
  "Concerns", unless the revision made them worse or claimed to fix them.
- Establish the revision's scope from the deliverables and the change from
  baseline; if it is unclear, ask the user.
- The summary verdict answers: more aligned than baseline (yes / mixed / no),
  and was the revision executed in a vision-aligned way?

**Vision alignment**
- Does the result advance one of the vision's claims: enlargement changes
  polarity organization or replacement; this suggests a physical basis for
  persistence; published migration data support it; mechanisms map to
  measurements?
- Material that serves none of these is a candidate for the supplement or
  removal.

**Narrative progress**
- For each narrative gap or wish: addressed, partly addressed, not addressed,
  or made worse. A gap within the revision's scope that is left unaddressed
  is a concern; one outside its scope is outstanding.
- For each narrative claim: still supported, strengthened, weakened, or no
  longer supported. Separate a regression in the deliverables from a
  narrative that is merely out of date (new numbers, renumbered figures);
  report the latter as a suggested narrative update, not a defect.
- Flag new material the narrative does not anticipate, and say whether it
  serves the vision.

**Role of each result.** Classify each result as pedagogical, novel, or
hypothesis-generating (it may be more than one), then ask whether its
presentation maximizes that role:
- *Pedagogical:* is it framed as established, cited, and kept brief? Is it
  doing setup work for something novel?
- *Novel:* is the new insight about size and polarity stated plainly and made
  visually prominent, rather than buried among checks?
- *Hypothesis-generating:* does it name a concrete measurement or protocol
  and the signature expected, and does it say which mechanism it would
  support and which it would fail to discriminate?

**Claim strength**
- The vision calibrates its verbs: enlargement *changes* polarity, this
  *suggests* a physical basis, the data *support* the hypothesis, and a
  measurement *could support* a mechanism without uniquely confirming it.
  Flag text or figures that claim more, or less, than this.
- Keep observed migration summaries, conditional OU persistence estimates and
  1D polarity-reversal readouts distinct; flag any conflation.
- Is length clearly a proxy for WGD-associated enlargement?

**Mechanism → measurement**
- Each mechanism should lead to a specific cellular measurement. Flag
  mechanisms that do not, and measurements that do not separate the
  mechanisms they are offered for.

**Presentation clarity**
- Can each figure's message be stated in one sentence, and can it be read
  without the caption?
- Are terminology, notation, colours and length conventions consistent across
  figures?
- Is the order of figures a coherent argument?
- Is information density appropriate, with diagnostics moved to the
  supplement?

## Review file format

```markdown
# Review: <target> — YYYY-MM-DD

Deliverables reviewed: <paths, with commit hash>
Baseline: <commit, and whether it came from the request, a prior review or the narrative>
Assumptions: <headless runs only; omit otherwise>

Revision scope: <what this revision set out to change>

## Summary
<3–5 sentences. Verdict: more vision-aligned than baseline (yes / mixed / no);
revision executed in a vision-aligned way (yes / partly / no); top concerns.>

## Concerns
| # | Severity | Criterion | Concern | Status | Codex's position |
|---|----------|-----------|---------|--------|------------------|
| 1 | major/minor/note | … | … | open | … |

## Outstanding
<distance remaining from the vision and narrative gaps outside this
revision's scope; brief, for planning, not counted against the revision>

## Narrative targets
| Narrative item | Section | Status | Evidence in deliverable |

## Narrative claims
| Claim | Figure | Still supported? | Note (regression vs stale narrative) |

## Result roles
| Result | Role(s) | Presented to maximize role? |

## Detail
<per-concern discussion, referencing figure panels and text locations>
```

On a follow-up review, carry forward earlier concern numbers and update their
status (open / addressed / disputed) by checking the deliverable.
