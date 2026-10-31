"""Check real-data symmetries, saved summaries and sample-pair distances.

Distance checks are descriptive QC, not a verification of biological identity.
The original labels and all primary results remain unchanged.
"""
import json
import argparse

import numpy as np
import pandas as pd
from scipy import stats

from audit_results import ROOT, sha256, validate_saved_summaries
from permutation import moderated_t, prepare_matrices
from scenarios import paths, result_base, add_argument


def direct_fit(matrices, labels):
    """Four-cell fit using explicitly demeaned residuals, without sum subtraction."""
    effects = []
    residual = np.zeros(matrices[0].shape[0])
    for matrix in matrices:
        means = []
        for group in (labels, ~labels):
            cell = matrix[:, group]
            means.append(cell.mean(axis=1))
            cell -= means[-1][:, None]
            residual += np.einsum('ij,ij->i', cell, cell)
        effects.append(means[0]-means[1])
    df = 2*len(labels)-4
    effect = np.array(effects)
    t, total_df, _, _ = moderated_t(effect, residual/df, df,
                                   1/labels.sum()+1/(~labels).sum())
    p = 2*stats.t.sf(abs(t), total_df)
    counts = ((abs(effect)>.5) & (p<.05)).sum(axis=1)
    return effect, t, p, counts


def pair_distances(beta, samples):
    """All-CpG beta RMSE between every two arrays, with both directional ranks."""
    gram = np.zeros((beta.shape[1], beta.shape[1]))
    for start in range(0, beta.shape[0], 20000):
        chunk = np.asarray(beta[start:start+20000])
        gram += chunk.T @ chunk
    diagonal = np.diag(gram)
    distances = np.sqrt(np.maximum(
        (diagonal[:, None]+diagonal[None, :]-2*gram)/beta.shape[0], 0))
    people = sorted(samples.participant.unique())
    indices = [samples.loc[samples.time==visit].set_index('participant')
               .loc[people, 'matrix_index'].to_numpy() for visit in (0, 1)]
    cross = distances[np.ix_(*indices)]
    rows = []
    for i, person in enumerate(people):
        rows.append({
            'participant':person, 'paired_beta_rmse':cross[i, i],
            'nearest_after_participant':people[int(cross[i].argmin())],
            'nearest_before_participant':people[int(cross[:, i].argmin())],
            'paired_after_distance_rank':1+int((cross[i]<cross[i, i]).sum()),
            'paired_before_distance_rank':1+int((cross[:, i]<cross[i, i]).sum()),
        })
    table = pd.DataFrame(rows)
    info = {
        'metric':'root mean squared beta difference across every deposited CpG',
        'before_arrays_with_annotated_after_as_nearest':int((table.paired_after_distance_rank==1).sum()),
        'after_arrays_with_annotated_before_as_nearest':int((table.paired_before_distance_rank==1).sum()),
        'participants':len(people),
        'minimum_distinct_array_rmse':float(distances[np.triu_indices(len(samples), 1)].min()),
        'interpretation':'descriptive distance QC only; cannot establish identity, swaps or exclusions; labels retained',
    }
    return table, info


def main(scenario='corrected_visits'):
    sample_path, folder = paths(ROOT, 'diet_code', scenario)
    out = result_base(ROOT, scenario)/'audit'
    out.mkdir(exist_ok=True)
    consistency = validate_saved_summaries(folder)
    shape = tuple(json.loads((ROOT/'data/processed/audit.json').read_text())['matrix_shape'])
    beta = np.memmap(ROOT/'data/processed/beta.bin', dtype='<f8', mode='r', shape=shape)
    samples = pd.read_csv(sample_path,
                          dtype={'participant':str, 'diet_code':str})
    samples['matrix_index'] = np.arange(len(samples))
    matrices, labels, _ = prepare_matrices(beta, samples, 'diet_code')
    fitted = direct_fit(matrices, labels)
    saved = [pd.read_csv(folder/f'{visit}_limma.csv.gz') for visit in ('before','after')]
    for actual, column in zip(fitted[:3], ('logFC','t','P.Value')):
        np.testing.assert_allclose(actual, np.array([table[column] for table in saved]),
                                   rtol=1e-9, atol=1e-10)
    observed = pd.read_csv(folder/'limma_summary.csv').set_index('time')
    np.testing.assert_array_equal(fitted[3], observed.loc[['before','after'],'highlighted'])
    invariances = {}
    order = np.random.default_rng(4824).permutation(len(labels))
    cases = [
        ('swap_group_names', matrices, ~labels,
         (-fitted[0], -fitted[1], fitted[2], fitted[3])),
        ('reverse_visits', matrices[::-1], labels, tuple(x[::-1] for x in fitted)),
        ('reorder_participants', [m[:, order] for m in matrices], labels[order], fitted),
    ]
    for name, values, groups, expected in cases:
        actual = direct_fit(values, groups)
        errors = {}
        for key, value, reference in zip(('effect','t','p'), actual[:3], expected[:3]):
            np.testing.assert_allclose(value, reference, rtol=1e-9, atol=1e-10)
            errors['max_abs_'+key+'_error'] = float(np.max(abs(value-reference)))
        np.testing.assert_array_equal(actual[3], expected[3])
        invariances[name] = {**errors, 'red_counts':actual[3].tolist(), 'passed':True}
    table, pair_info = pair_distances(beta, samples)
    table.to_csv(out/'pair_distance_check.csv', index=False)
    report = {
        'annotation_scenario':scenario,
        'samples_sha256':sha256(sample_path),
        'summary_consistency':consistency,
        'direct_residual_fit_matches_saved_limma':True,
        'invariances':invariances,
        'closest_distance_to_observed_cutoffs':{
            visit:{'absolute_effect_0_5':float((tab.logFC.abs()-.5).abs().min()),
                   'raw_p_0_05':float((tab['P.Value']-.05).abs().min())}
            for visit, tab in zip(('before','after'), saved)},
        'pair_distance_qc':pair_info,
        'count_file_sha256':sha256(folder/'permutation_null_counts.csv'),
        'permutation_script_sha256':sha256(ROOT/'scripts/permutation.py'),
    }
    (out/'invariance_summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    add_argument(parser)
    main(parser.parse_args().scenario)
