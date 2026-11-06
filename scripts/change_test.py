"""Test the after-minus-before change in red-site counts using paired shuffles."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scenarios import paths, add_argument, DEFAULT_SCENARIO
from run_inputs import validate_run_inputs

ROOT = Path(__file__).resolve().parents[1]


def summarize_change(observed_before, observed_after, shuffled_before, shuffled_after):
    """Upper-tail test of an increase, subtracting within each shared assignment."""
    before = np.asarray(shuffled_before, dtype=float)
    after = np.asarray(shuffled_after, dtype=float)
    observed = np.asarray([observed_before, observed_after], dtype=float)
    if before.ndim != 1 or not before.size or before.shape != after.shape:
        raise ValueError("Expected nonempty, equally sized vectors of paired counts")
    for counts in (before, after, observed):
        if (not np.isfinite(counts).all() or (counts < 0).any()
                or (counts != np.floor(counts)).any()):
            raise ValueError("Input counts must be finite nonnegative integers")
    changes = after.astype(np.int64) - before.astype(np.int64)
    observed_change = int(observed_after) - int(observed_before)
    extreme = int(np.count_nonzero(changes >= observed_change))
    p = (extreme + 1) / (len(changes) + 1)
    return {
        "observed_before": int(observed_before),
        "observed_after": int(observed_after),
        "observed_change": observed_change,
        "shuffled_median_change": float(np.median(changes)),
        "shuffled_q025_change": float(np.quantile(changes, 0.025)),
        "shuffled_q975_change": float(np.quantile(changes, 0.975)),
        "shuffles_ge_observed": extreme,
        "permutations": len(changes),
        "change_p_upper": p,
        "change_p_mc_se_approx": float(np.sqrt(p * (1-p) / len(changes))),
    }, changes


def run(field="diet_code", scenario=DEFAULT_SCENARIO):
    _, folder = paths(ROOT, field, scenario)
    null = pd.read_csv(folder / "permutation_null_counts.csv")
    observed = pd.read_csv(folder / "limma_summary.csv").set_index("time")
    metadata = json.loads((folder / "permutation_run.json").read_text())
    validate_run_inputs(ROOT, metadata, scenario, field)
    if metadata.get('annotation_scenario','deposited') != scenario:
        raise ValueError('Annotation scenario differs from the fitted run')
    if metadata["group_field"] != field or len(null) != metadata["repetitions"]:
        raise ValueError("Saved count pairs do not match the permutation run metadata")
    for time in ("before", "after"):
        if int(observed.loc[time, "highlighted"]) != metadata["observed_counts"][time]:
            raise ValueError("Observed counts differ from those validated for this permutation run")
    result, changes = summarize_change(observed.loc["before", "highlighted"],
        observed.loc["after", "highlighted"], null.before, null.after)
    paired = null[["before", "after"]].copy()
    paired.insert(0, "shuffle", np.arange(1, len(paired)+1))
    paired["change"] = changes
    paired.to_csv(folder / "permutation_null_changes.csv", index=False)
    pd.DataFrame([result]).to_csv(folder / "change_test_summary.csv", index=False)
    report = {
        "annotation_scenario": scenario,
        "permutation_scheme": metadata.get('permutation_scheme', 'unrestricted'),
        "block_file_sha256": metadata.get('block_file_sha256'),
        "randomization": metadata['assumption'],
        "samples_sha256": metadata.get("samples_sha256"),
        "statistic": "C_after - C_before, with the same participant assignment at both visits",
        "red_site_rule": "limma raw two-sided p < 0.05 and abs(logFC) > 0.5",
        "tail": "upper: shuffled change >= observed change, including ties",
        "p_formula": "(1 + shuffles_ge_observed) / (B + 1)",
        "question": "Is the observed increase in red-site count unusually large under the shuffle scheme?",
        "group_field": field, "seed": metadata["seed"], "permutations": metadata["repetitions"],
        "permutation_unit": metadata["permutation_unit"],
        "variance_prior": metadata["variance_prior"],
        "source": "existing paired rows of permutation_null_counts.csv; no new shuffles needed",
        "limitations": [
            "Validity requires joint exchangeability under the recorded allocation restrictions, conditional on the analyzed subset.",
            "Recruitment, replacement and selection into the methylation subset are not reconstructed by the conditional shuffles.",
            "Each shuffle uses its own baseline count; the test is not conditional on a baseline count of 158.",
            "This is a change in thresholded-site counts, not a per-CpG diet-by-time interaction test.",
            "The statistic and upper-tail direction were chosen during exploratory review of the results.",
            "The p-value is not the probability that the null hypothesis is true or that there is no diet effect.",
        ],
    }
    (folder / "change_test_run.json").write_text(json.dumps(report, indent=2)+"\n")
    print(pd.DataFrame([result]).to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_argument(parser)
    parser.add_argument("--group-field", choices=["diet_code", "diet_named"], default="diet_code")
    args = parser.parse_args()
    run(args.group_field, args.scenario)
