# Visit-label correction recorded on 31 October 2026

An author clarified that errors occurred while manually annotating the public
submission. The author reported that the metadata used for the paper's analysis
were correct, and that rerunning that analysis reproduced Figure 1B. The private
correspondence is not reproduced here; this note records the information relevant
to the reanalysis.

The author identified reversed time labels for donors 5 and 8 and errors in some
sex annotations. We interpret those donor numbers as the public participant IDs.
The local correction reverses four visit labels:

| Array ID | Participant | Deposited visit | Corrected visit |
| --- | ---: | --- | --- |
| ID_204792200112_R01C01 | 5 | After (1) | Before (0) |
| ID_204792170106_R01C01 | 5 | Before (0) | After (1) |
| ID_204792780062_R01C01 | 8 | Before (0) | After (1) |
| ID_204792170105_R01C01 | 8 | After (1) | Before (0) |

The [machine-readable patch](../data/corrections/visit_labels.csv) is applied by
[prepare_corrections.py](../scripts/prepare_corrections.py). It checks each array's
participant and original visit before applying the change, rejects applying the
patch twice, and verifies complete pairs. Measurements, matrix-column order,
participant identities, and all other annotation fields are preserved.
The corrected table is [samples_corrected_visits.csv](../data/processed/samples_corrected_visits.csv);
[visit_correction.json](../data/processed/visit_correction.json) records input/output
hashes and all four changes.

This is an explicit local correction of the archived public annotations. When
checked on 31 October 2026, the public SDRF and IDF files were still byte-for-byte
identical to the archived copies. Their download location had changed; the stable
BioStudies download endpoint now follows that move. The
[public-file check](../data/corrections/public_annotation_check.json) records URLs
and hashes. A later corrected public release should be reconciled as a new version,
not silently substituted into this historical patch.

## Model clarification and scope

The author supplied the initial unblocked model: four diet-by-time means on
M-values, `lmFit`, `contrasts.fit`, and default `eBayes`. The contrasts included
within-group visit changes, between-group comparisons at both visits, and their
interaction. The red rule was `decideTests(adjust.method="none", p.value=0.05,
lfc=0.5)`. Our implementation now includes these five contrasts and checks that
`decideTests` agrees exactly with the highlighted-site rule. Adding the other
contrasts does not change the two visit-specific raw p-values.

Applying only the four time-label changes to the numeric-code model gives
**158 sites before and 182 after (+24)**, matching the published counts. The prior
annotations gave 158/176. This resolves the six-site count discrepancy for this
reconstruction; it does not independently verify every CpG identity, the original
M-value matrix, or the participant-to-diet mapping.

The author also supplied a later model including diet, time, their interaction,
age, sex, and participant blocking through `duplicateCorrelation`. That later
model is outside the Figure 1B permutation analysis. We continue to construct
M-values as `log2(beta/(1-beta))` from the public processed matrix.

The author confirmed 24 participants per diet and sex totals of 30 women and
18 men. Individual sex corrections and a verified diet mapping have not been
provided. We retain the deposited sex annotations (16 men/32 women) and numeric
24/24 group assignment. Sex is not a covariate in this initial model, so its labels
do not enter the Figure 1B calculations. The conflicting named-diet field remains
an explicitly unverified sensitivity analysis.

Original randomization-block membership is still unavailable. All shuffles are
unrestricted and conditional on group sizes. The author's model clarification
does not resolve this allocation limitation.

## Recalculation and preserved history

All 9,999 shuffled models were refitted with the corrected visits, using seed 70899
and the same participant assignments as the historical run. The new null
distribution is not obtained by applying a +24 threshold to the old distribution.
The new primary statistic is `C_after - C_before = 182 - 158 = +24`.

- Current fits, counts, summaries and figures: `results/corrected_visits/`.
- Original fits, counts, summaries and figures: `results/diet_code/` and `results/diet_named/`.
- Original audit evidence: `results/audit/`.
- Current audit evidence: `results/corrected_visits/audit/`.
- Updated interpretation and numerical results: [README](../README.md).

The historical p-value **0.2630** belongs to the original deposited visit labels
and the +18 statistic. It must not be presented as the result for the corrected
+24 analysis. The original notebook and README are available in
[commit b7cccf8](https://github.com/pbenner/diet-methylation-reanalysis/tree/b7cccf8cf00b811d666af29c5194aeece418684f).

The correction was specified from the author's clarification, rather than selected
by searching for annotations that reproduce the counts. The count-change statistic
and its upper-tail direction remain exploratory choices made before this update.
