import contextlib
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
import change_test
import count_test
from run_inputs import validate_run_inputs


class SavedInputTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.folder = self.root/'results/blocked_randomization/diet_code'
        names = [f'results/blocked_randomization/diet_code/{name}' for name in
                 ('permutation_run.json', 'limma_summary.csv', 'permutation_null_counts.csv')]
        names += ['data/processed/samples_corrected_visits.csv',
                  'data/randomization/participant_blocks.csv']
        for name in names:
            destination = self.root/name
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(ROOT/name, destination)
        self.metadata = json.loads((self.folder/'permutation_run.json').read_text())

    def test_current_inputs_pass_and_scenarios_cannot_be_mixed(self):
        validate_run_inputs(self.root, self.metadata, 'blocked_randomization')
        with self.assertRaisesRegex(ValueError, 'scenario'):
            validate_run_inputs(self.root, self.metadata, 'corrected_visits')
        with self.assertRaisesRegex(ValueError, 'Group field'):
            validate_run_inputs(self.root, self.metadata, 'blocked_randomization', 'diet_named')

    def test_summary_commands_refuse_changed_blocks_before_writing(self):
        path = self.root/'data/randomization/participant_blocks.csv'
        path.write_text(path.read_text().replace('normal_1', 'normal_2', 1))
        for module in (change_test, count_test):
            with patch.object(module, 'ROOT', self.root), contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, 'Block-file SHA-256'):
                    module.run(scenario='blocked_randomization')
        self.assertFalse((self.folder/'change_test_summary.csv').exists())
        self.assertFalse((self.folder/'count_test_summary.csv').exists())

    def test_summary_commands_refuse_changed_sample_annotations(self):
        path = self.root/'data/processed/samples_corrected_visits.csv'
        path.write_bytes(path.read_bytes()+b'\n')
        for module in (change_test, count_test):
            with patch.object(module, 'ROOT', self.root):
                with self.assertRaisesRegex(ValueError, 'Sample-file SHA-256'):
                    module.run(scenario='blocked_randomization')

    def test_blocked_run_requires_hashes_and_correct_scheme(self):
        for key in ('samples_sha256', 'block_file_sha256', 'block_file'):
            metadata = dict(self.metadata)
            del metadata[key]
            with self.assertRaises(ValueError):
                validate_run_inputs(self.root, metadata, 'blocked_randomization')
        metadata = dict(self.metadata, permutation_scheme='unrestricted')
        with self.assertRaisesRegex(ValueError, 'Permutation scheme'):
            validate_run_inputs(self.root, metadata, 'blocked_randomization')


if __name__ == '__main__':
    unittest.main()
