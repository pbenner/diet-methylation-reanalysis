"""Validate the visit correction, new summaries and selected full-data R fits."""
import json
import subprocess

import numpy as np
import pandas as pd

from audit_results import ROOT, sha256, validate_saved_summaries
from prepare_corrections import apply_visit_corrections
from scenarios import paths, result_base


def main():
    sample_path, folder = paths(ROOT, scenario='corrected_visits')
    out = result_base(ROOT, 'corrected_visits')/'audit'
    out.mkdir(parents=True, exist_ok=True)
    original = pd.read_csv(ROOT/'data/processed/samples.csv',
                           dtype={'participant':str,'diet_code':str})
    samples = pd.read_csv(sample_path, dtype={'participant':str,'diet_code':str})
    patch = pd.read_csv(ROOT/'data/corrections/visit_labels.csv',dtype={'participant':str})
    expected = apply_visit_corrections(original, patch)
    pd.testing.assert_frame_equal(samples, expected)
    assert int((samples.time!=original.time).sum())==4
    archived_audit = json.loads((ROOT/'results/audit/audit_summary.json').read_text())
    assert sha256(ROOT/'data/processed/beta.bin')==archived_audit['binary_sha256']
    for name, expected_hash in archived_audit['source_sha256'].items():
        assert sha256(ROOT/'data/raw'/name)==expected_hash
    consistency = validate_saved_summaries(folder)
    metadata = json.loads((folder/'permutation_run.json').read_text())
    assert metadata['samples_sha256']==sha256(sample_path)
    counts = pd.read_csv(folder/'permutation_null_counts.csv')
    changes = counts.after-counts.before
    d = consistency['observed_change']
    indices = {0,len(counts)-1,int(changes.idxmin()),int(changes.idxmax())}
    indices.update(i for i in (99,999,4999) if i<len(counts))
    for value in (d-1,d,d+1):
        matches = np.flatnonzero(changes==value)
        if len(matches): indices.add(int(matches[0]))
    baseline = samples[samples.time==0].set_index('participant').sort_index()
    previous = original[original.time==0].set_index('participant').sort_index()
    pd.testing.assert_series_equal(baseline.diet_code,previous.diet_code)
    labels = baseline.diet_code.eq('0').to_numpy()
    rng = np.random.default_rng(metadata['seed'])
    assignments = [{'assignment':'observed',**baseline.diet_code.to_dict()}]
    for i in range(len(counts)):
        shuffled = rng.permutation(labels)
        if i in indices:
            assignments.append({'assignment':str(i+1),
                **dict(zip(baseline.index,np.where(shuffled,'0','1')))})
    assignment_file = out/'reference_assignments.csv'
    pd.DataFrame(assignments).to_csv(assignment_file,index=False)
    subprocess.run(['bash',str(ROOT/'scripts/run_r.sh'),str(ROOT/'scripts/audit_limma.R'),
                    str(assignment_file),str(out/'independent_limma_counts.csv'),
                    'corrected_visits'],cwd=ROOT,check=True)
    reference = pd.read_csv(out/'independent_limma_counts.csv',dtype={'assignment':str})
    observed = pd.read_csv(folder/'limma_summary.csv').set_index('time')
    for row in reference.itertuples():
        expected_counts = (observed.loc[['before','after'],'highlighted'].to_numpy()
                           if row.assignment=='observed' else counts.iloc[int(row.assignment)-1].to_numpy())
        np.testing.assert_array_equal([row.before,row.after],expected_counts)
    errors = {k:float(v) for k,v in reference[reference.assignment=='observed'].iloc[0].items()
              if k.startswith('max_abs_')}
    assert max(errors.values())<1e-8
    old_folder = ROOT/'results/diet_code'
    old_run = json.loads((old_folder/'permutation_run.json').read_text())
    old_counts = pd.read_csv(old_folder/'permutation_null_counts.csv')
    same_assignments = old_run['seed']==metadata['seed'] and len(old_counts)==len(counts)
    comparison = []
    for scenario, location in [('deposited',old_folder),('corrected_visits',folder)]:
        summary = pd.read_csv(location/'change_test_summary.csv').iloc[0].to_dict()
        comparison.append({'scenario':scenario,**summary})
    pd.DataFrame(comparison).to_csv(out/'scenario_comparison.csv',index=False)
    report = {
        'annotation_scenario':'corrected_visits', 'samples_sha256':sha256(sample_path),
        'exactly_four_visit_labels_changed':True,'other_annotations_and_array_order_unchanged':True,
        'raw_inputs_and_beta_matrix_match_original_audit':True,
        'same_participant_grouping_and_shuffle_assignments_as_deposited_run':same_assignments,
        'shuffle_rows_with_changed_before_or_after_counts':int((counts!=old_counts).any(axis=1).sum()) if same_assignments else None,
        'summary_consistency':consistency,
        'independent_real_data_r_fits':len(reference),'all_independent_r_counts_match':True,
        'observed_max_absolute_errors':errors,
        'permutation_count_sha256':sha256(folder/'permutation_null_counts.csv'),
        'permutation_script_sha256':sha256(ROOT/'scripts/permutation.py'),
    }
    (out/'audit_summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
