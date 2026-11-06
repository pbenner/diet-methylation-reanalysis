# Analysis review, 6 November 2026

The review found **no numerical error affecting the current blocked result**:
158 red sites before, 182 after, and an observed increase of 24. Among 9,999
within-block shuffles, 1,619 changes are at least 24, including 77 ties. The
one-sided Monte Carlo p-value is `(1619 + 1) / (9999 + 1) = 0.1620`.

## Evidence checked

- The author workbook was re-extracted and matches all 48 participants' diets and
  block assignments in the analysis CSV. The source workbook, sample annotations,
  raw inputs and prepared beta matrix match their recorded SHA-256 hashes.
- Every one of the **9,999 shuffled fits was replayed** in an isolated directory.
  All saved before/after count pairs match exactly; the primary results were not
  replaced. See the [full replay record](../results/blocked_randomization/audit/full_replay_summary.json).
- The existing independent factorial-coded R audit covers 11 full-data fits and
  matches every checked count. Its numerical-code hashes still match the current
  implementation. The largest observed statistic discrepancy is below 2.75e-13.
- Before/after samples stay together, and the same assignment is used across all
  CpGs. Diet counts are preserved within each of the six supplied blocks.
  Participant 57 remains fixed in the singleton replacement block.
- Plain CSV arithmetic independently reproduces the upper-tail count and p-value.
  All historical summaries also match their respective saved distributions.
- The observed overlap is **85 shared, 73 before-only and 97 after-only**. The
  paper reports 97 after-only sites alongside the 158/182 totals. This is additional
  aggregate agreement; original author CpG tables are still needed to establish
  site-by-site identity. [Paper, Figure 1B–C](https://onlinelibrary.wiley.com/doi/10.1002/mco2.70899).
- **24 automated tests pass**, including four new regression tests for saved-input
  validation. The notebook was executed again and its displayed figures checked
  against the saved exports.

The independent-seed result remains **0.1617**. A 95% Clopper–Pearson interval for
the underlying permutation-tail probability, based on 1,619 successes in 9,999
independent draws, is approximately **0.1547–0.1693**. This interval describes Monte
Carlo sampling uncertainty, not uncertainty about a biological effect.

## Validation issue fixed

Previously, the count/change summary commands could silently reuse saved shuffle
counts after the allocation CSV or sample annotations changed. This was reproduced
in a temporary copy by changing one block assignment. The current data match the
saved hashes, so this gap did **not** affect the reported result.

The shared [input validator](../scripts/run_inputs.py) now checks scenario, group
field, allocation scheme and recorded sample/block hashes before summaries,
plots, illustrative refits or replay runs proceed. Regression tests confirm that
changed inputs are rejected before a new summary is written. The historical
deposited run predates sample hashing; its separate source audit remains the
evidence for input identity rather than inventing a retrospective fingerprint.

## Statistical scope

The model and red-site thresholds match the author's initial Figure 1B code.
Allocation blocks restrict reassignment; they are not covariates in that model.
Every shuffle refits variance moderation. Raw limma p-values are part of the
counting rule; the single permutation p-value tests the selected count-change
statistic. Keeping raw thresholds inside a jointly shuffled statistic does not
require CpG independence. Validity depends on the permitted joint reassignments.
[Permutation inference reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC4010955/).

Including ties and applying the +1 correction are appropriate for this sampled
permutation procedure. A finite random sample of assignments does not justify
reporting zero p-values. [Phipson and Smyth](https://arxiv.org/abs/1603.05766).

The result is conditional on the analyzed subset and its observed diet counts in
each block. It assumes exchangeability within these blocks after selection into
the subset. The supplied workbook does not establish the full recruitment,
replacement or missingness mechanism. The statistic and one-sided direction were
chosen during exploratory review. These qualifications remain material.

The supported conclusion is: **for this count statistic and conditional block
shuffle scheme, the increase of 24 is not statistically significant**. This does
not establish equivalence of the diets, exclude individual methylation effects,
or refute the paper's other analyses. Individual sex corrections and the original
author M-value and CpG-result tables remain unavailable; sex is not used in this
initial model.
