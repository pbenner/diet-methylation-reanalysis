"""Joint participant-label shuffles for the original-rule red-site count test.

Restricted to complete data, unweighted four-cell designs, trend=False,
robust=False. Variance priors are re-estimated across all sites each shuffle.
Original randomization blocks are unavailable: these shuffles are unrestricted.
"""
import argparse
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import special, stats
from scenarios import paths, add_argument

ROOT = Path(__file__).resolve().parents[1]


def variance_prior(variance, df):
    """limma 3.54.0 fitFDist with a common residual df and no covariate."""
    if variance.ndim != 1 or len(variance) < 2 or df <= 0:
        raise ValueError("Need >=2 site variances and positive residual df")
    if not np.isfinite(variance).all() or (variance < -1e-15).any():
        raise ValueError("Invalid residual variance")
    x = np.maximum(variance, 0)
    median = np.median(x)
    x = np.maximum(x, 1e-5 * (median if median > 0 else 1))
    e = np.log(x) - special.digamma(df/2) + np.log(df/2)
    excess = e.var(ddof=1) - special.polygamma(1, df/2)
    if excess <= 0:
        return float(x.mean()), np.inf
    if excess > 1e7:
        y = 1 / np.sqrt(excess)
    elif excess < 1e-6:
        y = 1 / excess
    else:
        y = 0.5 + 1/excess
        for _ in range(51):
            tri = special.polygamma(1, y)
            step = tri * (1-tri/excess) / special.polygamma(2, y)
            y += step
            if -step/y < 1e-8:
                break
        else:
            raise RuntimeError("Inverse trigamma did not converge")
    return float(np.exp(e.mean() + special.digamma(y) - np.log(y))), float(2*y)


def moderated_t(effect, variance, df, unscaled):
    scale, prior_df = variance_prior(variance, df)
    posterior = scale if np.isinf(prior_df) else (df*variance + prior_df*scale)/(df+prior_df)
    total_df = min(df + prior_df, len(variance)*df)
    return effect / np.sqrt(posterior*unscaled), total_df, scale, prior_df


def highlighted_counts(effects, variance, df, unscaled):
    """Refit moderation and count original-rule red sites at each visit."""
    effects = np.asarray(effects)
    t, total_df, _, _ = moderated_t(effects, variance, df, unscaled)
    # Equivalent to limma's raw two-sided p < .05 and |effect| > .5.
    critical = stats.t.isf(0.025, total_df)
    return ((abs(t) > critical) & (abs(effects) > 0.5)).sum(axis=1)


def prepare_matrices(beta, samples, field):
    """Align two arrays by participant ID; never permute the two visits separately."""
    if beta.ndim != 2 or beta.shape[1] != len(samples):
        raise ValueError("Expected one metadata row for each matrix column")
    matrix_indices = samples.matrix_index.to_numpy()
    if (not np.issubdtype(matrix_indices.dtype, np.integer)
            or not np.array_equal(np.sort(matrix_indices), np.arange(beta.shape[1]))):
        raise ValueError("Matrix indices must identify each column exactly once")
    if samples.participant.isna().any() or not samples.time.isin([0, 1]).all():
        raise ValueError("Expected participant IDs and visits 0/1")
    people = sorted(samples.participant.unique())
    indices, aligned = [], []
    for time in (0, 1):
        sub = samples[samples.time == time]
        if sub.participant.duplicated().any() or set(sub.participant) != set(people):
            raise ValueError("Incomplete/duplicated participant pairs")
        sub = sub.set_index("participant").loc[people]
        aligned.append(sub)
        indices.append(sub.matrix_index.to_numpy())
    # Matrix indices identify beta columns, not positions in a reordered DataFrame.
    first, second = aligned
    if not np.array_equal(first[field].to_numpy(), second[field].to_numpy()):
        raise ValueError("Diet varies within participant")
    groups = ["0", "1"] if field == "diet_code" else ["Vegan Diet", "Meat-rich Diet"]
    if set(first[field]) != set(groups):
        raise ValueError("Unexpected group labels")
    labels = (first[field] == groups[0]).to_numpy()
    matrices = []
    for idx in indices:
        b = np.asarray(beta[:, idx])
        matrix = np.log2(b/(1-b))
        # Row-centering avoids cancellation when calculating within-cell SSE.
        matrix -= matrix.mean(axis=1, keepdims=True)
        matrices.append(np.ascontiguousarray(matrix))
    return matrices, labels, groups


def run(field, repetitions, seed, batch_size, workers=4, scenario='corrected_visits'):
    if repetitions < 1 or batch_size < 1 or workers < 1:
        raise ValueError("Repetitions, batch size, and workers must be positive")
    sample_path, out = paths(ROOT, field, scenario)
    audit = json.loads((ROOT / "data/processed/audit.json").read_text())
    shape = tuple(audit["matrix_shape"])
    beta = np.memmap(ROOT / "data/processed/beta.bin", dtype="<f8", mode="r", shape=shape)
    samples = pd.read_csv(sample_path, dtype={"participant": str, "diet_code": str})
    samples["matrix_index"] = np.arange(len(samples))
    matrices, labels, groups = prepare_matrices(beta, samples, field)
    n = len(labels)
    n1, n0 = int(labels.sum()), int((~labels).sum())
    unscaled = 1/n1 + 1/n0
    df = 2*n-4
    sums = [m.sum(axis=1) for m in matrices]
    sumsq = sum((m*m).sum(axis=1) for m in matrices)

    def from_sums(s1, s2):
        effects, residual = [], sumsq.copy()
        for total, treated in zip(sums, (s1, s2)):
            control = total-treated
            effects.append(treated/n1-control/n0)
            residual -= treated*treated/n1 + control*control/n0
        if (residual < -1e-9).any():
            raise ValueError("Negative residual sum of squares")
        return effects, np.maximum(residual, 0)/df

    effects, variance = from_sums(*(m @ labels.astype(float) for m in matrices))
    observed, dfs, scale, prior_df = moderated_t(np.array(effects), variance, df, unscaled)
    limma = [pd.read_csv(out / f"{time}_limma.csv.gz") for time in ("before", "after")]
    cpg = (ROOT / "data/processed/cpg_ids.txt").read_text().splitlines()
    # Independent validation against the actual archived R package, for every site.
    for j, tab in enumerate(limma):
        if tab.cpg.tolist() != cpg:
            raise ValueError("limma CpG order differs")
        np.testing.assert_allclose(effects[j], tab.logFC, atol=1e-10, rtol=1e-9)
        np.testing.assert_allclose(observed[j], tab.t, atol=1e-8, rtol=1e-8)
        np.testing.assert_allclose(2*stats.t.sf(abs(observed[j]), dfs), tab["P.Value"], atol=1e-9, rtol=1e-7)
    observed_counts = highlighted_counts(effects, variance, df, unscaled)
    expected_counts = [int(((tab["P.Value"] < 0.05) & (tab.logFC.abs() > 0.5)).sum()) for tab in limma]
    np.testing.assert_array_equal(observed_counts, expected_counts)
    print(f"Verified every statistic against limma 3.54.0 ({field})", flush=True)
    rng = np.random.default_rng(seed)
    null_counts = []

    def evaluate(pair):
        perm_effect, perm_var = from_sums(*pair)
        return highlighted_counts(perm_effect, perm_var, df, unscaled)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        for start in range(0, repetitions, batch_size):
            size = min(batch_size, repetitions-start)
            assignments = np.column_stack([rng.permutation(labels) for _ in range(size)]).astype(float)
            # A contiguous row per permutation avoids repeatedly reading strided columns.
            products = [np.ascontiguousarray(assignments.T @ m.T) for m in matrices]
            # map preserves assignment order; workers do not draw random numbers.
            for counts in pool.map(evaluate, zip(*products)):
                null_counts.append(counts)
            if start == 0 or (start+size) % 100 < batch_size or start+size == repetitions:
                print(f"{field}: {start+size}/{repetitions} permutations", flush=True)
    pd.DataFrame(null_counts, columns=["before", "after"]).to_csv(out / "permutation_null_counts.csv", index=False)
    metadata = {"group_field": field, "numerator": groups[0], "denominator": groups[1],
        "annotation_scenario": scenario,
        "samples_file": str(sample_path.relative_to(ROOT)),
        "samples_sha256": hashlib.sha256(sample_path.read_bytes()).hexdigest(),
        "repetitions": repetitions, "seed": seed, "batch_size": batch_size, "workers": workers,
        "tail": "upper: shuffled red-site count >= observed count, including ties",
        "correction": "(extreme+1)/(B+1)",
        "observed_counts": dict(zip(("before", "after"), expected_counts)),
        "permutation_unit": "participant; one assignment shared by both visits and every CpG",
        "variance_prior": "re-estimated across all sites for each permutation",
        "saved_count_rule": "limma raw two-sided p < 0.05 and abs(logFC) > 0.5",
        "count_test": "run scripts/count_test.py to test for an excess of red sites",
        "model": "unblocked four group-by-time means; structure and moderation confirmed by author clarification",
        "assumption": "unrestricted participant exchangeability, conditional on group counts; allocation blocks unavailable",
        "status": "provisional: numeric grouping retained; verified diet mapping and allocation blocks unavailable",
        "minimum_count_p": 1/(repetitions+1), "observed_df_prior": prior_df, "observed_s2_prior": scale}
    (out / "permutation_run.json").write_text(json.dumps(metadata, indent=2)+"\n")
    print(f"Saved {repetitions:,} paired red-site counts to {out / 'permutation_null_counts.csv'}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    add_argument(parser)
    parser.add_argument("--group-field", choices=["diet_code", "diet_named"], default="diet_code")
    parser.add_argument("--permutations", type=int, default=9999)
    parser.add_argument("--seed", type=int, default=70899)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    run(args.group_field, args.permutations, args.seed, args.batch_size, args.workers, args.scenario)
