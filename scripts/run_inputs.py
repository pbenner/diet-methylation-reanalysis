"""Refuse to combine saved permutation results with changed input annotations."""
import hashlib

from randomization import BLOCK_FILE, BLOCK_SCENARIO
from scenarios import paths


def validate_run_inputs(root, metadata, scenario, field='diet_code'):
    if metadata.get('annotation_scenario', 'deposited') != scenario:
        raise ValueError('Analysis scenario differs from the saved permutation run')
    if metadata.get('group_field') != field:
        raise ValueError('Group field differs from the saved permutation run')
    sample_path, _ = paths(root, field, scenario)
    sample_name = str(sample_path.relative_to(root))
    if metadata.get('samples_file', sample_name) != sample_name:
        raise ValueError('Sample file differs from the selected analysis scenario')
    sample_hash = metadata.get('samples_sha256')
    if scenario != 'deposited' and not sample_hash:
        raise ValueError('Saved permutation run has no sample-file hash')
    # The original deposited run predates sample hashing; its separate source audit
    # remains the evidence for input identity. Do not invent a historical fingerprint.
    if sample_hash and hashlib.sha256(sample_path.read_bytes()).hexdigest() != sample_hash:
        raise ValueError('Sample-file SHA-256 differs from the saved permutation run; refit the analysis')
    expected_scheme = 'within_randomization_blocks' if scenario == BLOCK_SCENARIO else 'unrestricted'
    if metadata.get('permutation_scheme', 'unrestricted') != expected_scheme:
        raise ValueError('Permutation scheme differs from the selected analysis scenario')
    if scenario == BLOCK_SCENARIO:
        if metadata.get('block_file') != BLOCK_FILE or not metadata.get('block_file_sha256'):
            raise ValueError('Saved blocked run has no valid allocation-file fingerprint')
        if hashlib.sha256((root/BLOCK_FILE).read_bytes()).hexdigest() != metadata['block_file_sha256']:
            raise ValueError('Block-file SHA-256 differs from the saved permutation run; refit the analysis')
