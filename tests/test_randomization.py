import sys
import unittest
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from randomization import align_blocks, allocation, draw_assignment, BLOCK_SCENARIO
from scenarios import paths


class RandomizationTests(unittest.TestCase):
    def setUp(self):
        self.people = ['1', '2', '3', '4', '5']
        self.labels = np.array([True, False, True, False, True])
        self.table = pd.DataFrame({'participant': self.people,
            'diet_code': ['0', '1', '0', '1', '0'],
            'block': ['normal_1', 'normal_1', 'replacement_1', 'replacement_1', 'replacement_2']})

    def test_within_block_counts_and_singleton_are_preserved(self):
        blocks = align_blocks(self.table, self.people, self.labels)
        rng = np.random.default_rng(17)
        seen = Counter()
        for _ in range(4096):
            assignment = draw_assignment(rng, self.labels, blocks)
            self.assertEqual(int(assignment[:2].sum()), 1)
            self.assertEqual(int(assignment[2:4].sum()), 1)
            self.assertTrue(assignment[4])
            seen[tuple(assignment)] += 1
        self.assertEqual(len(seen), 4)  # All and only the four legal assignments.
        self.assertTrue(all(abs(count/4096 - .25) < .04 for count in seen.values()))

    def test_table_order_does_not_control_participant_alignment(self):
        expected = align_blocks(self.table, self.people, self.labels)
        actual = align_blocks(self.table.iloc[::-1], self.people, self.labels)
        for a, b in zip(actual, expected):
            np.testing.assert_array_equal(a, b)

    def test_bad_participant_or_diet_mapping_is_rejected(self):
        for table in (self.table.iloc[:-1], pd.concat([self.table, self.table.iloc[:1]]),
                      self.table.assign(participant=['1', '2', '3', '4', '99'])):
            with self.assertRaises(ValueError):
                align_blocks(table, self.people, self.labels)
        wrong = self.table.copy()
        wrong.loc[0, 'diet_code'] = '1'
        with self.assertRaisesRegex(ValueError, 'disagree'):
            align_blocks(wrong, self.people, self.labels)

    def test_unrestricted_historical_rng_is_unchanged(self):
        a, b = np.random.default_rng(70899), np.random.default_rng(70899)
        for _ in range(50):
            np.testing.assert_array_equal(draw_assignment(a, self.labels), b.permutation(self.labels))

    def test_real_blocks_match_all_diets_and_constrain_every_draw(self):
        samples = pd.read_csv(paths(ROOT)[0], dtype={'participant': str, 'diet_code': str})
        baseline = samples[samples.time == 0].set_index('participant').sort_index()
        labels = baseline.diet_code.eq('0').to_numpy()
        blocks, info = allocation(ROOT, baseline.index, labels, BLOCK_SCENARIO)
        self.assertEqual(len(blocks), 6)
        self.assertEqual(sorted(map(len, blocks)), [1, 3, 10, 10, 11, 13])
        self.assertEqual(info['block_counts']['replacement_2'], {'vegan': 1, 'meat_rich': 0})
        rng = np.random.default_rng(70899)
        for _ in range(100):
            assignment = draw_assignment(rng, labels, blocks)
            for indices in blocks:
                self.assertEqual(int(assignment[indices].sum()), int(labels[indices].sum()))
            self.assertTrue(assignment[list(baseline.index).index('57')])

    def test_blocked_results_are_separate_but_use_same_corrected_samples(self):
        current_samples, current_results = paths(ROOT)
        previous_samples, previous_results = paths(ROOT, scenario='corrected_visits')
        self.assertEqual(current_samples, previous_samples)
        self.assertNotEqual(current_results, previous_results)
        with self.assertRaisesRegex(ValueError, 'numeric'):
            paths(ROOT, field='diet_named')


if __name__ == '__main__':
    unittest.main()
