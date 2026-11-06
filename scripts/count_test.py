"""Test the observed red-site count using saved joint participant shuffles.

The red-site rule remains limma raw p < 0.05 and abs(logFC) > 0.5.
The allocation restrictions are recorded with each saved permutation run.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scenarios import paths, add_argument, DEFAULT_SCENARIO
from run_inputs import validate_run_inputs

ROOT = Path(__file__).resolve().parents[1]


def summarize_counts(observed, shuffled):
    """Upper-tail Monte Carlo test; include ties and the observed assignment."""
    values = np.asarray(shuffled, dtype=float)
    if (values.ndim != 1 or not len(values) or not np.isfinite(values).all()
            or (values < 0).any() or (values != np.floor(values)).any()):
        raise ValueError("Shuffled counts must be a nonempty vector of nonnegative integers")
    if not np.isfinite(observed) or observed < 0 or observed != int(observed):
        raise ValueError("Observed count must be a nonnegative integer")
    extreme = int(np.count_nonzero(values >= observed))
    p = (extreme + 1) / (len(values) + 1)
    return {
        "observed_count": int(observed),
        "shuffled_median": float(np.median(values)),
        "shuffled_q025": float(np.quantile(values, 0.025)),
        "shuffled_q975": float(np.quantile(values, 0.975)),
        "shuffles_ge_observed": extreme,
        "permutations": len(values),
        "count_p_upper": p,
        "count_p_mc_se_approx": float(np.sqrt(p * (1-p) / len(values))),
    }


def run(field="diet_code", scenario=DEFAULT_SCENARIO):
    _, folder = paths(ROOT, field, scenario)
    null = pd.read_csv(folder / "permutation_null_counts.csv")
    observed = pd.read_csv(folder / "limma_summary.csv").set_index("time")
    metadata = json.loads((folder / "permutation_run.json").read_text())
    validate_run_inputs(ROOT, metadata, scenario, field)
    if metadata.get('annotation_scenario','deposited') != scenario:
        raise ValueError('Annotation scenario differs from the fitted run')
    if metadata["group_field"] != field or len(null) != metadata["repetitions"]:
        raise ValueError("Saved count table does not match permutation run metadata")
    for time in ("before", "after"):
        if int(observed.loc[time, "highlighted"]) != metadata["observed_counts"][time]:
            raise ValueError("Observed limma counts differ from those validated for this permutation run")
    rows = [{"time": time, **summarize_counts(observed.loc[time, "highlighted"], null[time])}
            for time in ("before", "after")]
    table = pd.DataFrame(rows)
    table.to_csv(folder / "count_test_summary.csv", index=False)
    report = {
        "annotation_scenario": scenario,
        "permutation_scheme": metadata.get('permutation_scheme', 'unrestricted'),
        "block_file_sha256": metadata.get('block_file_sha256'),
        "samples_sha256": metadata.get("samples_sha256"),
        "statistic": "number of CpGs passing limma raw two-sided p < 0.05 and abs(logFC) > 0.5",
        "tail": "upper: shuffled count >= observed count, including ties",
        "p_formula": "(1 + shuffles_ge_observed) / (B + 1)",
        "interpretation": "test for an excess of red sites, separately at each visit; no across-visit adjustment",
        "group_field": field,
        "permutations": metadata["repetitions"],
        "seed": metadata["seed"],
        "permutation_unit": metadata["permutation_unit"],
        "variance_prior": metadata["variance_prior"],
        "randomization": metadata['assumption'],
        "original_trial_methods": "https://doi.org/10.3390/microorganisms9040727",
        "limitations": [
            "Calibration requires joint exchangeability of whole participant profiles under the recorded allocation restrictions.",
            "The analyzed subset and its observed diet counts are held fixed; selection and replacement mechanisms are not reconstructed.",
            "A large p-value does not establish absence of a diet effect or identify individual false positives.",
        ],
    }
    (folder / "count_test_run.json").write_text(json.dumps(report, indent=2) + "\n")
    print(table.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_argument(parser)
    parser.add_argument("--group-field", choices=["diet_code", "diet_named"], default="diet_code")
    args = parser.parse_args()
    run(args.group_field, args.scenario)
