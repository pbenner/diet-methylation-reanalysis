import sys
import unittest
import tempfile
import hashlib
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
from prepare_corrections import apply_visit_corrections
from scenarios import paths
from fetch_data import download


class CorrectionTests(unittest.TestCase):
    def setUp(self):
        self.samples = pd.read_csv(ROOT/'data/processed/samples.csv',
                                   dtype={'participant':str,'diet_code':str})
        self.patch = pd.read_csv(ROOT/'data/corrections/visit_labels.csv',dtype={'participant':str})

    def test_correction_changes_only_four_explicit_visits(self):
        corrected = apply_visit_corrections(self.samples,self.patch)
        changed = corrected.time!=self.samples.time
        self.assertEqual(int(changed.sum()),4)
        self.assertEqual(set(corrected.loc[changed,'participant']),{'5','8'})
        pd.testing.assert_frame_equal(corrected.drop(columns='time'),self.samples.drop(columns='time'))
        self.assertEqual(corrected[corrected.time==0].diet_code.value_counts().to_dict(),{'0':24,'1':24})
        # Verify array identity, not DataFrame row position, controls the patch.
        reordered = self.samples.sample(frac=1,random_state=8)
        actual = apply_visit_corrections(reordered,self.patch)
        pd.testing.assert_frame_equal(actual.sort_index(),corrected)

    def test_correction_refuses_wrong_identity_or_double_application(self):
        corrected = apply_visit_corrections(self.samples,self.patch)
        with self.assertRaisesRegex(ValueError,'Original participant/time'):
            apply_visit_corrections(corrected,self.patch)
        wrong = self.patch.copy()
        wrong.loc[0,'participant']='8'
        with self.assertRaisesRegex(ValueError,'Original participant/time'):
            apply_visit_corrections(self.samples,wrong)

    def test_scenarios_keep_original_files_separate(self):
        old_samples, old_results = paths(ROOT,scenario='deposited')
        new_samples, new_results = paths(ROOT,scenario='corrected_visits')
        self.assertNotEqual(old_samples,new_samples)
        self.assertNotEqual(old_results,new_results)
        self.assertEqual(old_results,ROOT/'results/diet_code')

    def test_download_pins_contents_even_when_sizes_match(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            source, destination = directory/'source.txt', directory/'download.txt'
            source.write_bytes(b'original')
            original_hash = hashlib.sha256(b'original').hexdigest()
            download(source.as_uri(),destination,8,original_hash)
            self.assertEqual(destination.read_bytes(),b'original')
            destination.write_bytes(b'modified')
            with self.assertRaisesRegex(ValueError,'SHA-256 differs'):
                download(source.as_uri(),destination,8,original_hash)
            destination.unlink()
            source.write_bytes(b'modified')
            with self.assertRaisesRegex(ValueError,'SHA-256 differs'):
                download(source.as_uri(),destination,8,original_hash)
            self.assertFalse(destination.exists())


if __name__=='__main__':
    unittest.main()
