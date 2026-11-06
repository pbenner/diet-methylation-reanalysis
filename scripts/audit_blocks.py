"""Verify supplied allocations, replay restrictions and independently fit selected shuffles in R."""
import json
import subprocess
from math import comb, prod

import numpy as np
import pandas as pd

from audit_results import ROOT, sha256, validate_saved_summaries
from prepare_corrections import apply_visit_corrections
from scenarios import paths, result_base
from randomization import BLOCK_SCENARIO, BLOCK_FILE


def main():
    sample_path, folder = paths(ROOT, scenario=BLOCK_SCENARIO)
    out = result_base(ROOT, BLOCK_SCENARIO)/'audit'
    out.mkdir(parents=True, exist_ok=True)
    original = pd.read_csv(ROOT/'data/processed/samples.csv', dtype={'participant': str, 'diet_code': str})
    samples = pd.read_csv(sample_path, dtype={'participant': str, 'diet_code': str})
    patch = pd.read_csv(ROOT/'data/corrections/visit_labels.csv', dtype={'participant': str})
    pd.testing.assert_frame_equal(samples, apply_visit_corrections(original, patch))
    old_audit = json.loads((ROOT/'results/audit/audit_summary.json').read_text())
    assert sha256(ROOT/'data/processed/beta.bin') == old_audit['binary_sha256']
    for name, digest in old_audit['source_sha256'].items():
        assert sha256(ROOT/'data/raw'/name) == digest
    metadata = json.loads((folder/'permutation_run.json').read_text())
    provenance = json.loads((ROOT/'data/randomization/provenance.json').read_text())
    assert metadata['samples_sha256'] == sha256(sample_path)
    assert metadata['block_file_sha256'] == provenance['output_sha256'] == sha256(ROOT/BLOCK_FILE)
    assert metadata['permutation_scheme'] == 'within_randomization_blocks'
    workbook = ROOT/provenance['source']
    if workbook.exists():
        from prepare_randomization import read_workbook
        assert sha256(workbook) == provenance['source_sha256']
        pd.testing.assert_frame_equal(read_workbook(workbook),
            pd.read_csv(ROOT/BLOCK_FILE, dtype={'participant': str, 'diet_code': str}))
    baseline = samples[samples.time == 0].set_index('participant').sort_index()
    blocks = pd.read_csv(ROOT/BLOCK_FILE, dtype=str).set_index('participant')
    assert blocks.index.is_unique and set(blocks.index) == set(baseline.index)
    blocks = blocks.loc[baseline.index]
    pd.testing.assert_series_equal(blocks.diet_code, baseline.diet_code)
    # Reconstruct block membership directly from IDs, independently of allocation().
    groups = [list(blocks.index[blocks.block == block]) for block in sorted(blocks.block.unique())]
    observed_labels = baseline.diet_code.eq('0')
    consistency = validate_saved_summaries(folder)
    counts = pd.read_csv(folder/'permutation_null_counts.csv')
    changes = counts.after-counts.before
    d = consistency['observed_change']
    indices = {0, len(counts)-1, int(changes.idxmin()), int(changes.idxmax())}
    indices.update(i for i in (99, 999, 4999) if i < len(counts))
    for value in (d-1, d, d+1):
        matches = np.flatnonzero(changes == value)
        if len(matches):
            indices.add(int(matches[0]))
    rng = np.random.default_rng(metadata['seed'])
    assignments = [{'assignment': 'observed', **baseline.diet_code.to_dict()}]
    for i in range(len(counts)):
        shuffled = observed_labels.copy()
        for people in groups:
            shuffled.loc[people] = rng.permutation(observed_labels.loc[people].to_numpy())
            assert int(shuffled.loc[people].sum()) == int(observed_labels.loc[people].sum())
        assert shuffled['57'] and int(shuffled.sum()) == 24
        if i in indices:
            assignments.append({'assignment': str(i+1),
                **dict(zip(baseline.index, np.where(shuffled, '0', '1')))})
    assignment_file = out/'reference_assignments.csv'
    pd.DataFrame(assignments).to_csv(assignment_file, index=False)
    print('All allocation restrictions verified; checking independent R fits.', flush=True)
    subprocess.run(['bash', str(ROOT/'scripts/run_r.sh'), str(ROOT/'scripts/audit_limma.R'),
                    str(assignment_file), str(out/'independent_limma_counts.csv'),
                    BLOCK_SCENARIO], cwd=ROOT, check=True)
    reference = pd.read_csv(out/'independent_limma_counts.csv', dtype={'assignment': str})
    observed = pd.read_csv(folder/'limma_summary.csv').set_index('time')
    for row in reference.itertuples():
        expected = (observed.loc[['before', 'after'], 'highlighted'].to_numpy()
                    if row.assignment == 'observed' else counts.iloc[int(row.assignment)-1].to_numpy())
        np.testing.assert_array_equal([row.before, row.after], expected)
    errors = {k: float(v) for k, v in reference[reference.assignment == 'observed'].iloc[0].items()
              if k.startswith('max_abs_')}
    assert max(errors.values()) < 1e-8
    prior = ROOT/'results/corrected_visits/diet_code'
    for name in ('limma_summary.csv', 'separate_time_sensitivity.csv'):
        pd.testing.assert_frame_equal(pd.read_csv(folder/name), pd.read_csv(prior/name))
    for visit in ('before', 'after'):
        pd.testing.assert_frame_equal(pd.read_csv(folder/f'{visit}_limma.csv.gz'),
                                      pd.read_csv(prior/f'{visit}_limma.csv.gz'))
    comparison = []
    for scenario in ('deposited', 'corrected_visits', BLOCK_SCENARIO):
        _, location = paths(ROOT, scenario=scenario)
        row = pd.read_csv(location/'change_test_summary.csv').iloc[0].to_dict()
        comparison.append({'scenario': scenario, **row})
    pd.DataFrame(comparison).to_csv(out/'scenario_comparison.csv', index=False)
    space = prod(comb(len(people), int(observed_labels.loc[people].sum())) for people in groups)
    assert space == provenance['conditional_assignment_space']
    report = {
        'annotation_scenario': BLOCK_SCENARIO,
        'samples_sha256': sha256(sample_path), 'block_file_sha256': sha256(ROOT/BLOCK_FILE),
        'source_workbook_present_and_verified': workbook.exists(),
        'all_48_author_diets_match_numeric_codes': True,
        'all_shuffles_preserve_each_block_count_and_singleton': True,
        'independent_assignment_replay': 'Explicit participant-ID groups; does not call production allocation or assignment functions.',
        'conditional_assignment_space': space,
        'observed_site_results_identical_to_corrected_unrestricted_scenario': True,
        'raw_inputs_and_beta_matrix_match_original_audit': True,
        'four_visit_corrections_verified': True,
        'summary_consistency': consistency,
        'independent_real_data_r_fits': len(reference), 'all_independent_r_counts_match': True,
        'observed_max_absolute_errors': errors,
        'permutation_count_sha256': sha256(folder/'permutation_null_counts.csv'),
        'source_sha256': {name: sha256(ROOT/'scripts'/name) for name in
                         ('permutation.py', 'randomization.py', 'prepare_randomization.py', 'audit_limma.R')},
    }
    (out/'audit_summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
