"""Fit one seeded participant reassignment with the original R/limma script.

The resulting per-site p-values are ordinary limma p-values for one shuffled
grouping. The count test still uses all saved shuffled assignments.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd
from scenarios import paths, add_argument

ROOT = Path(__file__).resolve().parents[1]


def run(seed=None, force=False, scenario='corrected_visits'):
    sample_path, folder = paths(ROOT, 'diet_code', scenario)
    out = folder / "shuffled_example"
    metadata = json.loads((folder / "permutation_run.json").read_text())
    if metadata["group_field"] != "diet_code":
        raise ValueError("The saved run must use numeric diet codes")
    manifest = out / "run.json"
    previous = json.loads(manifest.read_text()) if manifest.exists() else {}
    # Keep the selected example seed when plots are regenerated without a new seed.
    if seed is None:
        seed = previous.get("seed", metadata["seed"])
    if seed < 0:
        raise ValueError("The shuffling seed must be nonnegative")
    sources = [
        "scripts/shuffled_example.py", "scripts/figure1b_limma.R", "scripts/run_r.sh",
        "data/processed/audit.json", str(sample_path.relative_to(ROOT)), "data/processed/cpg_ids.txt",
        str((folder/"permutation_run.json").relative_to(ROOT)),
        str((folder/"permutation_null_counts.csv").relative_to(ROOT)),
    ]
    fingerprints = {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in sources}
    matrix = (ROOT / "data/processed/beta.bin").stat()
    inputs = {"source_sha256": fingerprints, "matrix_bytes": matrix.st_size,
              "matrix_mtime_ns": matrix.st_mtime_ns, "example_seed": seed}
    outputs = [out / name for name in ("assignments.csv", "before_limma.csv.gz",
               "after_limma.csv.gz", "limma_summary.csv", "R_session.txt")]
    if (not force and all(p.exists() for p in outputs) and previous.get("inputs") == inputs):
        print(f"Reusing the R/limma fit for example seed {seed}.")
        return

    samples = pd.read_csv(sample_path,
                          dtype={"participant": str, "diet_code": str})
    if not samples.groupby("participant").time.apply(lambda x: sorted(x) == [0, 1]).all():
        raise ValueError("Expected one before and one after sample per participant")
    if not samples.groupby("participant").diet_code.nunique().eq(1).all():
        raise ValueError("Diet code changes within a participant")
    people = sorted(samples.participant.unique())  # Same ordering as permutation.py.
    baseline = samples[samples.time == 0].set_index("participant").loc[people]
    labels = baseline.diet_code.eq(metadata["numerator"]).to_numpy()
    shuffled = np.random.default_rng(seed).permutation(labels)
    assignments = pd.DataFrame({"participant": people,
        "observed_group": baseline.diet_code.to_numpy(),
        "shuffled_group": np.where(shuffled, metadata["numerator"], metadata["denominator"])})
    assert int(shuffled.sum()) == int(labels.sum())
    out.mkdir(exist_ok=True)
    assignments.to_csv(out / "assignments.csv", index=False)
    subprocess.run(["bash", str(ROOT / "scripts/run_r.sh"),
                    str(ROOT / "scripts/figure1b_limma.R"), "diet_code",
                    str(out / "assignments.csv"), scenario], cwd=ROOT, check=True)

    actual = pd.read_csv(out / "limma_summary.csv").set_index("time")
    for time in ("before", "after"):
        if (int(actual.loc[time, "n_numerator"]) != int(shuffled.sum())
                or int(actual.loc[time, "n_denominator"]) != int((~shuffled).sum())):
            raise ValueError(f"{time}: R group sizes differ from the shuffled assignment")
        fitted = pd.read_csv(out / f"{time}_limma.csv.gz", usecols=["logFC", "P.Value"])
        count = int(((fitted["P.Value"] < 0.05) & (fitted.logFC.abs() > 0.5)).sum())
        if count != int(actual.loc[time, "highlighted"]):
            raise ValueError(f"{time}: saved CpG results do not match the R red-site count")
    verification = "group sizes and red-site counts verified against the R fitted results"
    if seed == metadata["seed"]:
        expected = pd.read_csv(folder / "permutation_null_counts.csv", nrows=1).iloc[0]
        for time in ("before", "after"):
            if int(actual.loc[time, "highlighted"]) != int(expected[time]):
                raise ValueError(f"{time}: independent R count differs from saved shuffle #1")
        verification += "; both counts also match the first saved count-test assignment"
    report = {"inputs": inputs, "shuffle_number": 1, "seed": seed, "annotation_scenario": scenario,
        "selection": "first assignment from the specified seed; seed selected before fitting",
        "method": "same R limma 3.54.0 model and thresholds as the observed fit",
        "p_values": "two-sided limma p-values fitted to this one shuffled grouping",
        "assignment_unit": metadata["permutation_unit"],
        "count_test_seed": metadata["seed"],
        "observed_counts_in_this_shuffle": {t: int(actual.loc[t, "highlighted"]) for t in ("before", "after")},
        "verification": verification}
    manifest.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Verified shuffled example for seed {seed}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_argument(parser)
    parser.add_argument("--seed", type=int, help="Seed for the example only; defaults to the last selected example seed")
    parser.add_argument("--force", action="store_true", help="Refit even when the inputs match the saved example")
    args = parser.parse_args()
    run(seed=args.seed, force=args.force, scenario=args.scenario)
