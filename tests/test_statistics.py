"""Independent numerical and pairing checks for the statistical implementation."""
import itertools
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from permutation import highlighted_counts, moderated_t, prepare_matrices
from count_test import summarize_counts
from change_test import summarize_change


class StatisticsTests(unittest.TestCase):
    def test_change_uses_paired_counts_and_includes_ties(self):
        result, changes = summarize_change(10, 12, [10, 100, 20, 30], [12, 98, 23, 29])
        np.testing.assert_array_equal(changes, [2, -2, 3, -1])
        self.assertEqual(result["observed_change"], 2)
        self.assertEqual(result["shuffles_ge_observed"], 2)
        self.assertEqual(result["change_p_upper"], 3/5)
        result, _ = summarize_change(10, 8, [10, 100, 20, 30], [12, 98, 23, 29])
        self.assertEqual(result["observed_change"], -2)
        self.assertEqual(result["change_p_upper"], 1)

    def test_change_rejects_unpaired_count_vectors(self):
        with self.assertRaisesRegex(ValueError, "paired counts"):
            summarize_change(10, 12, [10, 20], [12])

    def test_count_tail_includes_ties_and_observed_assignment(self):
        result = summarize_counts(3, [1, 3, 3, 7])
        self.assertEqual(result["shuffles_ge_observed"], 3)
        self.assertEqual(result["count_p_upper"], 4/5)
        self.assertEqual(summarize_counts(8, [1, 3, 3, 7])["count_p_upper"], 1/5)
        self.assertEqual(summarize_counts(0, [1, 3, 3, 7])["count_p_upper"], 1)

    def test_threaded_moderation_is_deterministic(self):
        rng = np.random.default_rng(4)
        cases = [(rng.normal(size=(2, 100)), rng.lognormal(size=100)) for _ in range(12)]
        def calculate(case):
            return moderated_t(case[0], case[1], 20, 1/6)
        expected = list(map(calculate, cases))
        with ThreadPoolExecutor(max_workers=4) as pool:
            actual = list(pool.map(calculate, cases))
        for a, e in zip(actual, expected):
            np.testing.assert_array_equal(a[0], e[0])
            self.assertEqual(a[1:], e[1:])

    def test_pair_alignment_ignores_row_order(self):
        samples = pd.DataFrame({"participant": ["b", "a", "b", "a"],
            "time": [1, 0, 0, 1], "diet_code": ["1", "0", "1", "0"],
            "matrix_index": range(4)})
        beta = np.array([[0.8, 0.2, 0.4, 0.6], [0.3, 0.4, 0.5, 0.6]])
        m, labels, _ = prepare_matrices(beta, samples, "diet_code")
        np.testing.assert_array_equal(labels, [True, False])
        for actual, idx in zip(m, ([1, 2], [3, 0])):
            expected = np.log2(beta[:, idx]/(1-beta[:, idx]))
            expected -= expected.mean(axis=1, keepdims=True)
            np.testing.assert_allclose(actual, expected)
        samples.loc[0, "diet_code"] = "0"
        with self.assertRaisesRegex(ValueError, "Diet varies"):
            prepare_matrices(beta, samples, "diet_code")

    def test_every_small_assignment_against_actual_limma(self):
        # All 20 assignments of 6 participants to groups of 3, each with 2 visits.
        rng = np.random.default_rng(17)
        matrix = rng.normal(size=(60, 12))*np.exp(rng.normal(size=(60, 1)))
        assignments = np.zeros((20, 6))
        for i, choice in enumerate(itertools.combinations(range(6), 3)):
            assignments[i, list(choice)] = 1
        with tempfile.TemporaryDirectory() as temp:
            temp = Path(temp)
            np.savetxt(temp/"m.csv", matrix, delimiter=",")
            np.savetxt(temp/"assignments.csv", assignments, delimiter=",")
            script = temp/"oracle.R"
            script.write_text('''
.libPaths(c(file.path(getwd(), ".tools/r-library"), .libPaths()))
library(limma)
stopifnot(as.character(packageVersion("limma"))=="3.54.0")
args <- commandArgs(trailingOnly=TRUE)
m <- as.matrix(read.csv(file.path(args[1], "m.csv"), header=FALSE))
a <- as.matrix(read.csv(file.path(args[1], "assignments.csv"), header=FALSE))
ans <- list()
for (i in 1:nrow(a)) {
 group <- rep(a[i,], 2)
 time <- rep(0:1, each=6)
 cell <- factor(paste0("g", group, "t", time), levels=c("g0t0", "g1t0", "g0t1", "g1t1"))
 d <- model.matrix(~0+cell)
 colnames(d) <- levels(cell)
 c <- makeContrasts(before=g1t0-g0t0, after=g1t1-g0t1, levels=d)
 f <- eBayes(contrasts.fit(lmFit(m,d), c))
 ans[[i]] <- cbind(f$t[,1], f$t[,2], f$p.value[,1], f$p.value[,2], f$coefficients[,1], f$coefficients[,2])
}
write.table(do.call(rbind,ans), file.path(args[1], "oracle.csv"), sep=",", row.names=FALSE, col.names=FALSE)
''')
            subprocess.run(["bash", str(ROOT/"scripts/run_r.sh"), str(script), str(temp)], cwd=ROOT, check=True, capture_output=True)
            oracle = np.loadtxt(temp/"oracle.csv", delimiter=",").reshape(20, 60, 6)
        for i, assignment in enumerate(assignments):
            effects = []
            sse = np.zeros(60)
            for m in (matrix[:, :6], matrix[:, 6:]):
                one, zero = m[:, assignment == 1], m[:, assignment == 0]
                effects.append(one.mean(axis=1)-zero.mean(axis=1))
                sse += ((one-one.mean(axis=1, keepdims=True))**2).sum(axis=1)
                sse += ((zero-zero.mean(axis=1, keepdims=True))**2).sum(axis=1)
            t, df, _, _ = moderated_t(np.array(effects), sse/8, 8, 2/3)
            np.testing.assert_allclose(t.T, oracle[i, :, :2], rtol=1e-9, atol=1e-10)
            np.testing.assert_allclose((2*stats.t.sf(abs(t), df)).T, oracle[i, :, 2:4], rtol=1e-8, atol=1e-10)
            # Independently validate the exact red-site count for every assignment.
            counts = highlighted_counts(effects, sse/8, 8, 2/3)
            expected_counts = ((oracle[i, :, 2:4] < 0.05) & (abs(oracle[i, :, 4:6]) > 0.5)).sum(axis=0)
            np.testing.assert_array_equal(counts, expected_counts)

if __name__ == "__main__":
    unittest.main()
