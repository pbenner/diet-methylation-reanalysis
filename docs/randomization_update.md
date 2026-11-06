# Author-supplied diet coding and allocation blocks, 6 November 2026

Max Storz supplied the diet coding **0 = vegan, 1 = meat-rich** and a workbook,
`ID_BlockAssignment_DrBenner.xlsx`, identifying the diet and randomization block
for every participant in the methylation dataset. He explained that an external
statistician created the randomization list: planned blocks had 13 allocations
through observation 52, followed by blocks of 4 for replacements.

The [minimal allocation table](../data/randomization/participant_blocks.csv)
contains the deposited participant ID, verified diet and explicit block membership.
The [provenance record](../data/randomization/provenance.json) records the workbook
SHA-256, extraction details and output hash. The original workbook stays local and
is ignored by Git; the analysis can be reproduced from the extracted CSV.
The correspondence and workbook document metadata are not copied into the repository.

All **48 participant IDs match exactly**, and all **48 diet assignments agree with
the deposited numeric diet codes**. The short IDs in workbook column B match the
deposited IDs, and the full IDs in column A equal 142000 plus the short ID.
No participant identities or numeric diet assignments were changed. The conflicting
deposited named-diet field is retained as historical source data, but is not used
to assign diets in the current analysis.

| Explicit block | Planned size | Analyzed vegan | Analyzed meat-rich | Analyzed total |
| --- | ---: | ---: | ---: | ---: |
| normal_1 | 13 | 6 | 4 | 10 |
| normal_2 | 13 | 4 | 6 | 10 |
| normal_3 | 13 | 7 | 6 | 13 |
| normal_4 | 13 | 4 | 7 | 11 |
| replacement_1 | 4 | 2 | 1 | 3 |
| replacement_2 | 4 | 1 | 0 | 1 |
| Total | 60 planned slots | 24 | 24 | 48 |

Planned block sizes describe the allocation list; they are not the numbers of
included methylation participants. Block membership is read from the workbook,
not inferred by splitting sorted IDs. The two replacement blocks are separate
from the original blocks with the same numerical suffix.

## Conditional permutation procedure

For each shuffle, independently permute the observed diet labels **within each
explicit block**, preserving that block's diet totals among analyzed participants.
The total remains 24 vegan and 24 meat-rich. Participant 57 is the only included
member of `replacement_2`, so this conditional procedure holds that participant's
vegan assignment fixed. There are 74,918,844,000 possible assignments under these
restrictions, including the observed assignment. We sample 9,999 times with
replacement using seed 70899; the separate stability run uses seed 70901.

Each assignment applies jointly to both visits and all CpGs. Measurements and the
four documented visit corrections stay fixed. Every shuffled fit re-estimates
limma's empirical Bayes variance prior and uses raw two-sided p < 0.05 and
abs(logFC) > 0.5. The signed count-change statistic and one-sided upper tail are
the same as in the preceding analysis. See the [README](../README.md) for results.

The **randomization blocks restrict reassignment**. They are not additional model
covariates or participant correlation blocks. The fitted model remains the initial
unblocked four diet-by-visit means supplied for Figure 1B. The observed fit is
therefore identical to the corrected-visits fit: **158 before, 182 after (+24)**.

## Assumptions and preserved history

The test is conditional on the analyzed subset and its observed diet counts in each
block. It assumes joint exchangeability of whole participant profiles within those
blocks under the null. This would follow from uniform random allocation within
blocks together with selection that permits this conditional exchangeability.
The workbook does not establish the full recruitment, replacement or missingness
mechanism. We do not reconstruct an unconditional test of the entire trial or claim
that blocking resolves possible selection bias. The statistic remains exploratory.

Current results are in `results/blocked_randomization/`. The corrected but
unrestricted run remains in `results/corrected_visits/` (p = 0.1881), and the
original-visit run remains in `results/diet_code/` (p = 0.2630 for +18).
The same seed now draws from a different assignment space; shuffled assignments
are not the same as in either unrestricted run.

This clarification resolves the diet-code meaning, participant-to-diet mapping
and analyzed participants' block membership. Individual sex corrections, the
author's original M-value matrix and site-by-site results remain unavailable.
Sex does not enter the initial Figure 1B model. Matching the published counts
does not establish exact identity of every highlighted CpG or reproduce the
paper's other analyses.
